# LASSO-based covariate and interaction selection for IV analysis (2SRI / control function)

# ==============================================================================
# Internal helpers
# ==============================================================================

#' Extract non-zero variable names from a cv.glmnet fit
#' @keywords internal
.selected_vars <- function(cvfit) {
  coefs <- coef(cvfit, s = "lambda.min")
  setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
}

#' Drop constant columns from a matrix
#' @keywords internal
.drop_constant_cols <- function(X) {
  keep   <- apply(X, 2, function(col) sd(col) > 0)
  dropped <- colnames(X)[!keep]
  if (length(dropped) > 0) {
    warning(
      length(dropped), " constant column(s) removed before LASSO: ",
      paste(dropped, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    X       = X[, keep, drop = FALSE],
    dropped = dropped
  )
}

#' Build a numeric matrix from a tibble and a vector of column names
#' @keywords internal
.build_matrix <- function(data, predictors) {
  predictors <- unique(trimws(predictors))
  predictors <- predictors[!is.na(predictors)]
  predictors <- predictors[predictors != ""]

  missing <- setdiff(predictors, colnames(data))
  if (length(missing) > 0) {
    stop(
      length(missing), " predictor(s) not found in dataset:\n",
      paste0("  - ", missing, collapse = "\n"),
      call. = FALSE
    )
  }
  as.matrix(data[, predictors, drop = FALSE])
}

#' Uniform (0/1) penalty factor — used for Part 1 (main effects)
#'
#' Returns 0 for unpenalised columns, 1 for all others.
#' @keywords internal
.penalty_factor <- function(col_names, unpenalized) {
  unpenalized <- unique(trimws(unpenalized))
  as.numeric(!col_names %in% unpenalized)
}

#' SD-based penalty factor — used for interaction stages
#'
#' Returns 1/SD for penalised columns, 0 for unpenalised.
#' Columns with SD == 0 (constant) are set to 0 (never dropped).
#' This corrects for the tendency of LASSO to drop rare binary interactions
#' due to their small coefficients rather than true unimportance.
#' @keywords internal
.penalty_factor_sd <- function(X, unpenalized) {
  sds           <- apply(X, 2, sd) # calculate sd for each column
  pf            <- 1 / sds
  pf[colnames(X) %in% unpenalized] <- 0
  pf
}

#' Run cv.glmnet with a fixed seed
#' @keywords internal
.cv_lasso <- function(X, y, pf, family, seed, standardize = TRUE) {
  set.seed(seed)
  glmnet::cv.glmnet(
    X, y,
    family         = family,
    alpha          = 1,
    penalty.factor = pf,
    standardize    = standardize, # we use glmnet's internal standardization for part 1; for part 2 (interaction term identification) we use the SD-based penalty factor 
    type.measure   = "deviance"
  )
}

#' Add pairwise X x X interaction columns to a dataset
#'
#' Returns a list with the augmented dataset and the new column names.
#' Warns if any generated name exceeds 32 characters (Stata hard limit).
#' @keywords internal
.add_pairwise_xx_interactions <- function(data, vars) {
  vars <- unique(trimws(vars))
  vars <- vars[!is.na(vars)]
  vars <- vars[vars != ""]

  if (length(vars) < 2) {
    return(list(data = data, col_names = character(0)))
  }

  pairs         <- utils::combn(vars, 2, simplify = FALSE)
  new_col_names <- vapply(pairs, function(p) paste0(p[1], "_x_", p[2]), character(1))

  too_long <- new_col_names[nchar(new_col_names) > 32]
  if (length(too_long) > 0) {
    warning(
      length(too_long), " pairwise X\u00d7X column name(s) exceed 32 characters ",
      "(Stata hard limit). Shorten the base variable names to avoid this.\n",
      "Offending names:\n",
      paste0("  ", too_long, collapse = "\n"),
      call. = FALSE
    )
  }

  for (i in seq_along(pairs)) {
    data[[new_col_names[[i]]]] <- data[[pairs[[i]][1]]] * data[[pairs[[i]][2]]]
  }
  list(data = data, col_names = new_col_names)
}

#' Add interaction columns between a set of base variables and a reference column
#'
#' @param data          Tibble to augment
#' @param base_vars     Character vector of base variable names
#' @param reference_col Name of the column to interact with (e.g. treatment, instrument)
#' @param suffix        Suffix appended to each base name (e.g. "_d", "_iv")
#' @return List with augmented `data` and `col_names` of new columns
#' @keywords internal
.add_interactions <- function(data, base_vars, reference_col, suffix) {
  base_vars <- trimws(base_vars)
  bad_idx   <- which(is.na(base_vars) | base_vars == "")
  if (length(bad_idx) > 0) {
    stop(
      "Empty/NA base var name(s) passed to .add_interactions() for suffix '",
      suffix, "'.\nIndices: ", paste(bad_idx, collapse = ", "),
      call. = FALSE
    )
  }
  if (!reference_col %in% colnames(data)) {
    stop("Reference column not found in data: ", reference_col, call. = FALSE)
  }

  new_col_names <- paste0(base_vars, suffix)
  for (i in seq_along(base_vars)) {
    data[[new_col_names[[i]]]] <- data[[base_vars[[i]]]] * data[[reference_col]]
  }
  list(data = data, col_names = new_col_names)
}

#' Resolve a glmnet family string to a stats family object
#' @keywords internal
.glm_family <- function(family_str) {
  switch(
    family_str,
    "binomial" = stats::binomial(),
    "gaussian"  = stats::gaussian(),
    stop(
      "Unsupported family for glm() re-estimation: '", family_str, "'.\n",
      "Supported values: 'binomial', 'gaussian'.",
      call. = FALSE
    )
  )
}

#' Compute goodness-of-fit statistics for a glm.fit object
#'
#' For gaussian outcomes: R-squared (1 - SS_res / SS_tot).
#' For binomial outcomes: McFadden pseudo-R-squared and Hosmer-Lemeshow test
#' (requires the \pkg{ResourceSelection} package).
#'
#' @param fit        A \code{glm.fit} object.
#' @param y          Numeric response vector used to fit the model.
#' @param family_str One of \code{"gaussian"} or \code{"binomial"}.
#' @return A named list of goodness-of-fit statistics (invisibly printed).
#' @keywords internal
.gof_stats <- function(fit, y, family_str) {
  fitted_vals <- fit$fitted.values

  if (family_str == "gaussian") {
    ss_res  <- sum((y - fitted_vals)^2)
    ss_tot  <- sum((y - mean(y))^2)
    r2      <- 1 - ss_res / ss_tot
    message(sprintf("  R-squared: %.4f", r2))
    return(invisible(list(r_squared = r2)))
  }

  if (family_str == "binomial") {
    # McFadden pseudo-R-squared: 1 - (log-lik model) / (log-lik null)
    ll_model <- sum(stats::dbinom(y, size = 1, prob = fitted_vals,     log = TRUE))
    ll_null  <- sum(stats::dbinom(y, size = 1, prob = mean(y), log = TRUE))
    mcfadden <- 1 - ll_model / ll_null

    message(sprintf("  McFadden pseudo-R-squared: %.4f", mcfadden))

    # Hosmer-Lemeshow test (10 groups by default)
    hl_result <- ResourceSelection::hoslem.test(y, fitted_vals, g = 10)
    message(sprintf(
      "  Hosmer-Lemeshow test: chi-sq=%.3f  df=%d  p=%.4f",
      hl_result$statistic, hl_result$parameter, hl_result$p.value
    ))
    hl <- hl_result

    return(invisible(list(mcfadden_r2 = mcfadden, hoslem_test = hl)))
  }
}

# ==============================================================================
# Exported function
# ==============================================================================

#' Run two-stage LASSO variable selection for IV analysis (2SRI / control function)
#'
#' Implements a three-part LASSO selection strategy for constructing the
#' exposure and outcome models in a two-stage residual inclusion (2SRI)
#' instrumental variable analysis.
#'
#' **Part 1 — Main effect selection (uniform penalty):**
#' Runs LASSO separately for the exposure model (forcing Z and prespecified
#' subgroups) and the outcome model (forcing D and prespecified subgroups).
#' Takes the union of selected main effects; prespecified subgroups are always
#' retained. Squared continuous variables should be included in 
#' \code{penalized_main_effects} as a pre-computed column upstream.
#'
#' **Outcome interaction selection (SD-based penalty):**
#' Models the outcome with unpenalised terms (union main effects + D +
#' prespecified subgroups + D×prespecified) and penalised interaction
#' candidates (D×X for non-prespecified X, plus all X×X pairs from the union).
#' Retains surviving penalised interactions as \code{retained_dx_outcome} and
#' \code{retained_xx_outcome}.
#'
#' **Exposure interaction selection (SD-based penalty):**
#' Models the exposure with unpenalised terms (union main effects + Z +
#' prespecified subgroups + Z×prespecified + Z×X for each X where D×X was
#' retained above) and penalised interaction candidates (remaining Z×X plus
#' all X×X pairs from the union). Surviving penalised terms are returned as
#' \code{retained_zx_exposure} and \code{retained_xx_exposure}.
#' The selected exposure model is refit via \code{glm.fit()} and generalised
#' residuals (D - fitted P(D=1|Z,X)) are computed.
#'
#' **Final outcome model (2SRI control function correction):**
#' Reruns the outcome interaction model augmented with the generalised residual.
#' The D×residual interaction term is included only when
#' \code{include_resid_d = TRUE}.
#'
#' @param dataset A data frame (or tibble) containing all required variables.
#'   All covariates must be pre-coded as numeric (dummy variables); factor
#'   columns are not supported. Age-squared should be pre-computed upstream
#'   and passed in via \code{prespecified_subgroups} or
#'   \code{penalized_main_effects}.
#' @param outcome String. Name of the outcome variable.
#' @param treatment String. Name of the binary treatment variable.
#' @param instrument String. Name of the instrumental variable.
#' @param prespecified_subgroups Character vector. Covariates forced into both
#'   models at every stage and whose interactions with D and Z are always
#'   unpenalised.
#' @param penalized_main_effects Character vector. Additional covariates
#'   eligible for selection in Part 1 (penalised at the main-effect stage).
#' @param family_stage1 String. \code{glmnet} family for the exposure model.
#'   Default \code{"binomial"}.
#' @param family_stage2 String. \code{glmnet} family for the outcome model.
#'   Default \code{"gaussian"}. Note: tobit is not supported by \code{glmnet};
#'   use \code{"gaussian"} as an approximation for censored outcomes.
#' @param include_resid_d Logical. Whether to include the D×(generalised
#'   residual) interaction term in the final outcome model (2SRI control
#'   function correction for heterogeneous treatment effects). Default
#'   \code{TRUE}. Set to \code{FALSE} to include only the residual main effect
#'   and compare model fit using the printed goodness-of-fit statistics.
#' @param seed Single integer. Random seed for reproducibility. Default
#'   \code{1276}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{union_main_effects}}{Character vector. Union of main effects
#'       selected in Part 1, always including prespecified subgroups.}
#'     \item{\code{retained_dx_outcome}}{Character vector. D×X interaction
#'       terms (non-prespecified X only) retained in the outcome interaction
#'       stage.}
#'     \item{\code{retained_xx_outcome}}{Character vector. X×X interaction
#'       terms retained in the outcome interaction stage.}
#'     \item{\code{retained_dx_xx_outcome}}{Character vector. Union of
#'       \code{retained_dx_outcome} and \code{retained_xx_outcome}; these
#'       terms enter the final outcome model.}
#'     \item{\code{retained_zx_exposure}}{Character vector. Z×X interaction
#'       terms (beyond forced terms) retained in the exposure interaction
#'       stage.}
#'     \item{\code{retained_xx_exposure}}{Character vector. X×X interaction
#'       terms retained in the exposure interaction stage.}
#'     \item{\code{final_exposure_terms}}{Character vector. All terms in the
#'       refitted exposure model.}
#'     \item{\code{generalised_residuals}}{Numeric vector. Residuals
#'       D - P̂(D=1|Z,X) from the refitted exposure model
#'       (length = \code{nrow(dataset)}).}
#'     \item{\code{final_outcome_terms}}{Character vector. All terms in the
#'       final 2SRI outcome model.}
#'     \item{\code{final_outcome_fit}}{A \code{glm.fit} object for the final
#'       outcome model.}
#'     \item{\code{gof}}{Named list of goodness-of-fit statistics for the
#'       final outcome model. For gaussian: \code{r_squared}. For binomial:
#'       \code{mcfadden_r2} and \code{hoslem_test}.}
#'     \item{\code{cvfit_p1_exposure}}{cv.glmnet fit — Part 1 exposure model.}
#'     \item{\code{cvfit_p1_outcome}}{cv.glmnet fit — Part 1 outcome model.}
#'     \item{\code{cvfit_outcome_interactions}}{cv.glmnet fit — outcome
#'       interaction selection stage.}
#'     \item{\code{cvfit_exposure_interactions}}{cv.glmnet fit — exposure
#'       interaction selection stage.}
#'   }
#'
#' @importFrom magrittr %>%
#' @importFrom glmnet cv.glmnet
#' @export
run_lasso_iv_selection <- function(
    dataset,
    outcome,
    treatment,
    instrument,
    prespecified_subgroups,
    penalized_main_effects,
    family_stage1  = "binomial",
    family_stage2  = "gaussian",
    include_resid_d = TRUE,
    seed           = 1276
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  outcome                <- trimws(outcome)
  treatment              <- trimws(treatment)
  instrument             <- trimws(instrument)
  prespecified_subgroups <- trimws(prespecified_subgroups)
  penalized_main_effects <- trimws(penalized_main_effects)

  bad_subgroups <- prespecified_subgroups[
    is.na(prespecified_subgroups) | prespecified_subgroups == ""
  ]
  if (length(bad_subgroups) > 0) {
    stop(
      "`prespecified_subgroups` contains ", length(bad_subgroups),
      " empty/NA entry(ies).",
      call. = FALSE
    )
  }

  bad_main <- penalized_main_effects[
    is.na(penalized_main_effects) | penalized_main_effects == ""
  ]
  if (length(bad_main) > 0) {
    stop(
      "`penalized_main_effects` contains ", length(bad_main),
      " empty/NA entry(ies).",
      call. = FALSE
    )
  }

  required_cols <- unique(c(
    outcome, treatment, instrument,
    prespecified_subgroups, penalized_main_effects
  ))
  missing_cols <- setdiff(required_cols, colnames(dataset))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from `dataset`:\n",
      paste0("  - ", missing_cols, collapse = "\n"),
      call. = FALSE
    )
  }

  supported_families <- c("binomial", "gaussian")
  if (!family_stage1 %in% supported_families) {
    stop(
      "`family_stage1` must be one of: ",
      paste(supported_families, collapse = ", "),
      call. = FALSE
    )
  }
  if (!family_stage2 %in% supported_families) {
    stop(
      "`family_stage2` must be one of: ",
      paste(supported_families, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.logical(include_resid_d) || length(include_resid_d) != 1) {
    stop("`include_resid_d` must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("`seed` must be a single integer.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Part 1 — Main effect selection (uniform penalty)
  # ---------------------------------------------------------------------------

  message("\n--- Part 1: Main effect selection (uniform penalty) ---")

  p1_exposure_preds <- unique(c(instrument, prespecified_subgroups, penalized_main_effects))
  p1_outcome_preds  <- unique(c(treatment,  prespecified_subgroups, penalized_main_effects))

  X_p1_exp <- .build_matrix(dataset, p1_exposure_preds)
  X_p1_out <- .build_matrix(dataset, p1_outcome_preds)

  pf_p1_exp <- .penalty_factor(
    colnames(X_p1_exp),
    unpenalized = c(instrument, prespecified_subgroups)
  )
  pf_p1_out <- .penalty_factor(
    colnames(X_p1_out),
    unpenalized = c(treatment, prespecified_subgroups)
  )

  cvfit_p1_exposure <- .cv_lasso(
    X_p1_exp, dataset[[treatment]], pf_p1_exp, family_stage1, seed
  )
  cvfit_p1_outcome <- .cv_lasso(
    X_p1_out, dataset[[outcome]],   pf_p1_out, family_stage2, seed
  )

  sel_p1_exp <- .selected_vars(cvfit_p1_exposure)
  sel_p1_out <- .selected_vars(cvfit_p1_outcome)

  # Union of selected main effects; strip Z and D (they enter separately),
  # always retain prespecified subgroups
  union_main_effects <- unique(c(
    setdiff(union(sel_p1_exp, sel_p1_out), c(instrument, treatment)),
    prespecified_subgroups
  ))
  union_main_effects <- union_main_effects[
    !is.na(trimws(union_main_effects)) & trimws(union_main_effects) != ""
  ]

  if (length(union_main_effects) == 0) {
    stop(
      "union_main_effects is empty after Part 1; cannot proceed.",
      call. = FALSE
    )
  }

  message(sprintf(
    "  Exposure model selected: %d  |  Outcome model selected: %d  |  Union: %d",
    length(sel_p1_exp), length(sel_p1_out), length(union_main_effects)
  ))
  message("  union_main_effects: ", paste(union_main_effects, collapse = ", "))

  # ---------------------------------------------------------------------------
  # Pairwise X x X interactions — created once, shared across both interaction
  # selection stages (each stage selects independently from this pool)
  # ---------------------------------------------------------------------------

  xx_result   <- .add_pairwise_xx_interactions(dataset, union_main_effects)
  dataset     <- xx_result$data
  all_xx_cols <- trimws(xx_result$col_names)

  message(sprintf("\n  Pairwise X\u00d7X terms created: %d", length(all_xx_cols)))

  # ---------------------------------------------------------------------------
  # Outcome interaction selection (SD-based penalty)
  #
  # Unpenalised: union main effects + D + prespecified subgroups + D*prespec
  # Penalised:   D*X (non-prespec X) + all X*X pairs
  # ---------------------------------------------------------------------------

  message("\n--- Outcome interaction selection (SD-based penalty) ---")

  dx_result   <- .add_interactions(dataset, union_main_effects, treatment, "_d")
  dataset     <- dx_result$data
  all_dx_cols <- trimws(dx_result$col_names)

  # D x prespec interactions are always unpenalised
  prespec_in_union <- trimws(intersect(prespecified_subgroups, union_main_effects))
  prespec_in_union <- prespec_in_union[!is.na(prespec_in_union) & prespec_in_union != ""]
  dx_cols_unpen    <- paste0(prespec_in_union, "_d")

  # Penalised pool: D x X (non-prespec) + all X x X
  dx_cols_penalized <- c(setdiff(all_dx_cols, dx_cols_unpen), all_xx_cols)

  out_int_unpenalized <- unique(c(
    union_main_effects, treatment, prespecified_subgroups, dx_cols_unpen
  ))
  out_int_predictors <- unique(c(out_int_unpenalized, dx_cols_penalized))
  out_int_predictors <- out_int_predictors[
    trimws(out_int_predictors) != "" & !is.na(out_int_predictors)
  ]

  X_out_int  <- .build_matrix(dataset, out_int_predictors)
  X_out_int  <- .drop_constant_cols(X_out_int)  
  pf_out_int <- .penalty_factor_sd(X_out_int, unpenalized = out_int_unpenalized)

  cvfit_outcome_interactions <- .cv_lasso(
    X_out_int, dataset[[outcome]], pf_out_int, family_stage2, seed, standardize = FALSE  
  )
  sel_out_int <- .selected_vars(cvfit_outcome_interactions)

  # Separate retained D*X and X*X terms for clarity
  retained_dx_outcome <- intersect(sel_out_int, setdiff(all_dx_cols, dx_cols_unpen))
  retained_xx_outcome <- intersect(sel_out_int, all_xx_cols)

  retained_dx_outcome <- retained_dx_outcome[nchar(trimws(retained_dx_outcome)) > 0]
  retained_xx_outcome <- retained_xx_outcome[nchar(trimws(retained_xx_outcome)) > 0]

  # Union feeds the final outcome model
  retained_dx_xx_outcome <- unique(c(retained_dx_outcome, retained_xx_outcome))

  message(sprintf(
    "  Penalised pool: D\u00d7X (non-prespec)=%d  X\u00d7X=%d  total=%d",
    length(setdiff(all_dx_cols, dx_cols_unpen)),
    length(all_xx_cols),
    length(dx_cols_penalized)
  ))
  message(sprintf(
    "  Retained: D\u00d7X=%d  X\u00d7X=%d  total=%d",
    length(retained_dx_outcome),
    length(retained_xx_outcome),
    length(retained_dx_xx_outcome)
  ))

  # ---------------------------------------------------------------------------
  # Exposure interaction selection (SD-based penalty)
  #
  # Unpenalised: union main effects + Z + prespecified subgroups
  #              + Z*prespec + Z*X for each X where D*X was retained above
  # Penalised:   remaining Z*X + all X*X pairs
  # ---------------------------------------------------------------------------

  message("\n--- Exposure interaction selection (SD-based penalty) ---")

  zx_result   <- .add_interactions(dataset, union_main_effects, instrument, "_iv")
  dataset     <- zx_result$data
  all_zx_cols <- trimws(zx_result$col_names)
  all_zx_cols <- all_zx_cols[!is.na(all_zx_cols) & all_zx_cols != "" & all_zx_cols != "_iv"]

  # Z x prespec: always unpenalised
  zx_prespec_unpen <- paste0(prespec_in_union, "_iv")

  # Z x X also unpenalised for each X where D x X was retained
  retained_dx_base  <- trimws(gsub("_d$", "", retained_dx_outcome))
  retained_dx_base  <- retained_dx_base[!is.na(retained_dx_base) & retained_dx_base != ""]
  retained_zx_unpen <- paste0(retained_dx_base, "_iv")
  retained_zx_unpen <- retained_zx_unpen[nchar(trimws(retained_zx_unpen)) > 0]

  zx_cols_unpen <- unique(trimws(c(zx_prespec_unpen, retained_zx_unpen)))
  zx_cols_unpen <- zx_cols_unpen[
    !is.na(zx_cols_unpen) & zx_cols_unpen != "" & zx_cols_unpen != "_iv"
  ]

  # Penalised pool: remaining Z x X (not forced) + all X x X
  zx_cols_penalized <- c(setdiff(all_zx_cols, zx_cols_unpen), all_xx_cols)

  exp_int_unpenalized <- unique(c(
    union_main_effects, instrument, prespecified_subgroups, zx_cols_unpen
  ))
  exp_int_predictors <- unique(c(exp_int_unpenalized, zx_cols_penalized))
  exp_int_predictors <- exp_int_predictors[
    trimws(exp_int_predictors) != "" & !is.na(exp_int_predictors)
  ]

  X_exp_int  <- .build_matrix(dataset, exp_int_predictors)
  X_exp_int  <- .drop_constant_cols(X_exp_int)       
  pf_exp_int <- .penalty_factor_sd(X_exp_int, unpenalized = exp_int_unpenalized)

  cvfit_exposure_interactions <- .cv_lasso(
    X_exp_int, dataset[[treatment]], pf_exp_int, family_stage1, seed, standardize = FALSE
  )
  sel_exp_int <- .selected_vars(cvfit_exposure_interactions)

  # Surface retained Z*X and X*X terms separately
  retained_zx_exposure <- intersect(sel_exp_int, setdiff(all_zx_cols, zx_cols_unpen))
  retained_xx_exposure <- intersect(sel_exp_int, all_xx_cols)

  retained_zx_exposure <- retained_zx_exposure[nchar(trimws(retained_zx_exposure)) > 0]
  retained_xx_exposure <- retained_xx_exposure[nchar(trimws(retained_xx_exposure)) > 0]

  final_exposure_terms <- unique(c(exp_int_unpenalized, sel_exp_int))

  n_prespec_forced   <- length(zx_prespec_unpen)
  n_from_retained_dx <- length(setdiff(retained_zx_unpen, zx_prespec_unpen))
  n_total_forced     <- length(zx_cols_unpen)
  n_penalized_zx     <- length(setdiff(all_zx_cols, zx_cols_unpen))
  n_penalized_xx     <- length(all_xx_cols)
  n_penalized_total  <- length(zx_cols_penalized)

  message(sprintf(
    "  Forced Z\u00d7X: prespec=%d + from retained D\u00d7X=%d => %d (after dedup)",
    n_prespec_forced, n_from_retained_dx, n_total_forced
  ))
  message(sprintf(
    "  Penalised pool: Z\u00d7X=%d  X\u00d7X=%d  total=%d",
    n_penalized_zx, n_penalized_xx, n_penalized_total
  ))
  message(sprintf(
    "  Retained: Z\u00d7X=%d  X\u00d7X=%d",
    length(retained_zx_exposure), length(retained_xx_exposure)
  ))

  # ---------------------------------------------------------------------------
  # Refit exposure model via glm.fit() and compute generalised residuals
  #
  # Generalised residuals = D - P̂(D=1|Z,X), i.e. response-scale residuals.
  # These are valid for binary treatment regardless of the outcome type.
  # ---------------------------------------------------------------------------

  X_exposure_refit <- cbind(intercept = 1, .build_matrix(dataset, final_exposure_terms))
  X_exposure_refit <- .drop_constant_cols(X_exposure_refit)
  final_exposure_terms <- setdiff(colnames(X_exposure_refit), "intercept")

  exposure_refit <- stats::glm.fit(
    x      = X_exposure_refit,
    y      = dataset[[treatment]],
    family = .glm_family(family_stage1)
  )

  generalised_residuals <- dataset[[treatment]] - exposure_refit$fitted.values

  message(sprintf(
    "\n  Generalised residuals  mean: %.4f  |  SD: %.4f",
    mean(generalised_residuals),
    stats::sd(generalised_residuals)
  ))

  # ---------------------------------------------------------------------------
  # Final outcome model — 2SRI control function correction
  #
  # Always includes the generalised residual main effect.
  # D x residual interaction included only when include_resid_d = TRUE;
  # set to FALSE to compare fit with/without the interaction term.
  # ---------------------------------------------------------------------------

  message("\n--- Final outcome model: 2SRI correction ---")
  message(sprintf("  include_resid_d = %s", include_resid_d))

  resid_col   <- ".ctrl_resid"
  resid_d_col <- ".ctrl_resid_x_d"

  dataset[[resid_col]]   <- generalised_residuals
  dataset[[resid_d_col]] <- generalised_residuals * dataset[[treatment]]

  final_outcome_terms <- unique(c(
    out_int_unpenalized,
    retained_dx_xx_outcome,
    resid_col,
    if (include_resid_d) resid_d_col
  ))

  final_outcome_X <- cbind(intercept = 1, .build_matrix(dataset, final_outcome_terms))
  final_outcome_X <- .drop_constant_cols(final_outcome_X)

  final_outcome_fit <- stats::glm.fit(
    x      = final_outcome_X,
    y      = dataset[[outcome]],
    family = .glm_family(family_stage2)
  )

  message(sprintf(
    "  Final outcome model terms: %d", length(final_outcome_terms)
  ))

  # ---------------------------------------------------------------------------
  # Goodness-of-fit statistics for the final outcome model (always printed)
  # ---------------------------------------------------------------------------

  message("\n--- Goodness-of-fit: final outcome model ---")
  gof <- .gof_stats(final_outcome_fit, dataset[[outcome]], family_stage2)

  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------

  message("\n========== LASSO IV Selection Summary ==========")
  message(sprintf("  Part 1  union main effects:              %d", length(union_main_effects)))
  message(sprintf("  Outcome interactions  D\u00d7X retained:    %d", length(retained_dx_outcome)))
  message(sprintf("  Outcome interactions  X\u00d7X retained:    %d", length(retained_xx_outcome)))
  message(sprintf("  Exposure interactions Z\u00d7X retained:    %d", length(retained_zx_exposure)))
  message(sprintf("  Exposure interactions X\u00d7X retained:    %d", length(retained_xx_exposure)))
  message(sprintf("  Final exposure model terms:              %d", length(final_exposure_terms)))
  message(sprintf("  Final outcome model terms:               %d", length(final_outcome_terms)))
  message(sprintf("  D\u00d7residual interaction included:        %s", include_resid_d))
  message("=================================================\n")

  list(
    union_main_effects          = union_main_effects,
    retained_dx_outcome         = retained_dx_outcome,
    retained_xx_outcome         = retained_xx_outcome,
    retained_dx_xx_outcome      = retained_dx_xx_outcome,
    retained_zx_exposure        = retained_zx_exposure,
    retained_xx_exposure        = retained_xx_exposure,
    final_exposure_terms        = final_exposure_terms,
    generalised_residuals       = generalised_residuals,
    final_outcome_terms         = final_outcome_terms,
    final_outcome_fit           = final_outcome_fit,
    gof                         = gof,
    cvfit_p1_exposure           = cvfit_p1_exposure,
    cvfit_p1_outcome            = cvfit_p1_outcome,
    cvfit_outcome_interactions  = cvfit_outcome_interactions,
    cvfit_exposure_interactions = cvfit_exposure_interactions,

    prespecified_subgroups = prespecified_subgroups,
    instrument             = instrument,
    treatment              = treatment,
    outcome                = outcome,

    dropped_constant_cols  = dropped_constant_cols
  )
}