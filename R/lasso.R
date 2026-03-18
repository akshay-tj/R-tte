# lasso.R
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

#' SD-based penalty factor — used for Parts 2a/2b (interaction terms)
#'
#' Returns 1/SD for penalised columns, 0 for unpenalised.
#' Columns with SD == 0 (constant) are set to 0 (never dropped).
#' This corrects for the tendency of LASSO to drop rare binary interactions
#' due to their small coefficients rather than true unimportance.
#' @keywords internal
.penalty_factor_sd <- function(X, unpenalized) {
  sds        <- apply(X, 2, sd)
  sds[sds == 0] <- NA_real_
  pf         <- 1 / sds
  pf[is.na(pf)] <- 0   # sd == 0: do not penalise (constant column)
  pf[colnames(X) %in% unpenalized] <- 0
  pf
}

#' Run cv.glmnet with a fixed seed
#' @keywords internal
.cv_lasso <- function(X, y, pf, family, seed) {
  set.seed(seed)
  glmnet::cv.glmnet(
    X, y,
    family         = family,
    alpha          = 1,
    penalty.factor = pf,
    standardize    = TRUE,
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

  pairs        <- utils::combn(vars, 2, simplify = FALSE)
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
#' @param data         Tibble to augment
#' @param base_vars    Character vector of base variable names
#' @param reference_col  Name of the column to interact with (e.g. treatment, instrument)
#' @param suffix       Suffix appended to each base name (e.g. "_d", "_iv")
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
    "gaussian" = stats::gaussian(),
    stop(
      "Unsupported family for glm() re-estimation: '", family_str, "'.\n",
      "Supported values: 'binomial', 'gaussian'.",
      call. = FALSE
    )
  )
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
#' retained.
#'
#' **Part 2a — Outcome interaction selection (SD-based penalty):**
#' Models the outcome with unpenalised terms (union main effects + D +
#' prespecified subgroups + D×prespecified) and penalised interaction
#' candidates (D×X for non-prespecified X, plus all X×X pairs). Retains
#' surviving penalised interactions.
#'
#' **Part 2b — Exposure interaction selection (SD-based penalty):**
#' Models the exposure with unpenalised terms (union main effects + Z +
#' prespecified subgroups + Z×prespecified + Z×X for each X where D×X was
#' retained in Part 2a) and penalised interaction candidates (remaining Z×X
#' plus all X×X pairs). Refits the selected exposure model via `glm.fit()`
#' and computes generalised residuals.
#'
#' **Final outcome model:**
#' Reruns the Part 2a outcome model augmented with the generalised residual
#' and its interaction with D (2SRI control function correction).
#'
#' @param dataset A data frame (or tibble) containing all required variables.
#'   All covariates must be pre-coded as numeric (dummy variables); factor
#'   columns are not supported.
#' @param outcome String. Name of the outcome variable.
#' @param treatment String. Name of the binary treatment variable.
#' @param instrument String. Name of the instrumental variable.
#' @param prespecified_subgroups Character vector. Covariates forced into both
#'   models at every stage and whose interactions with D and Z are always
#'   unpenalised.
#' @param penalized_main_effects Character vector. Additional covariates
#'   eligible for selection in Part 1 (penalised at the main-effect stage).
#' @param family_stage1 String. `glmnet` family for the exposure model.
#'   Default `"binomial"`.
#' @param family_stage2 String. `glmnet` family for the outcome model.
#'   Default `"gaussian"`. Note: tobit is not supported by `glmnet`; use
#'   `"gaussian"` as an approximation for censored outcomes.
#' @param seed Single integer. Random seed for reproducibility. Default `1276`.
#'
#' @return A named list:
#'   \describe{
#'     \item{`union_main_effects`}{Character vector. Union of main effects
#'       selected in Part 1, always including prespecified subgroups.}
#'     \item{`retained_dx_interactions`}{Character vector. D×X and X×X
#'       interaction terms retained in Part 2a.}
#'     \item{`final_exposure_terms`}{Character vector. All terms in the
#'       refitted Part 2b exposure model.}
#'     \item{`generalised_residuals`}{Numeric vector. Residuals from the
#'       refitted exposure model (length = `nrow(dataset)`).}
#'     \item{`final_outcome_terms`}{Character vector. All terms in the final
#'       2SRI outcome model (Part 2a terms + residual + D×residual).}
#'     \item{`final_outcome_fit`}{A `glm.fit` object for the final outcome
#'       model.}
#'     \item{`cvfit_p1_exposure`}{`cv.glmnet` fit — Part 1 exposure model.}
#'     \item{`cvfit_p1_outcome`}{`cv.glmnet` fit — Part 1 outcome model.}
#'     \item{`cvfit_p2a_outcome`}{`cv.glmnet` fit — Part 2a outcome model.}
#'     \item{`cvfit_p2b_exposure`}{`cv.glmnet` fit — Part 2b exposure model.}
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
    family_stage1 = "binomial",
    family_stage2 = "gaussian",
    seed          = 1276
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
  # Pairwise X x X interactions (created once, used in both 2a and 2b)
  # ---------------------------------------------------------------------------

  xx_result   <- .add_pairwise_xx_interactions(dataset, union_main_effects)
  dataset     <- xx_result$data
  all_xx_cols <- trimws(xx_result$col_names)

  message(sprintf("\n  Pairwise X\u00d7X terms created: %d", length(all_xx_cols)))

  # ---------------------------------------------------------------------------
  # Part 2a — Outcome interaction selection (SD-based penalty)
  # ---------------------------------------------------------------------------

  message("\n--- Part 2a: Outcome interaction selection (SD-based penalty) ---")

  dx_result   <- .add_interactions(dataset, union_main_effects, treatment, "_d")
  dataset     <- dx_result$data
  all_dx_cols <- trimws(dx_result$col_names)

  # D x prespec interactions: always unpenalised
  prespec_in_union <- trimws(intersect(prespecified_subgroups, union_main_effects))
  prespec_in_union <- prespec_in_union[!is.na(prespec_in_union) & prespec_in_union != ""]
  dx_cols_unpen    <- paste0(prespec_in_union, "_d")

  # Penalised pool: D x X (non-prespec) + all X x X
  dx_cols_penalized <- c(setdiff(all_dx_cols, dx_cols_unpen), all_xx_cols)

  p2a_unpenalized <- unique(c(
    union_main_effects, treatment, prespecified_subgroups, dx_cols_unpen
  ))
  p2a_predictors <- unique(c(p2a_unpenalized, dx_cols_penalized))
  p2a_predictors <- p2a_predictors[trimws(p2a_predictors) != "" & !is.na(p2a_predictors)]

  X_p2a  <- .build_matrix(dataset, p2a_predictors)
  pf_p2a <- .penalty_factor_sd(X_p2a, unpenalized = p2a_unpenalized)

  cvfit_p2a_outcome <- .cv_lasso(
    X_p2a, dataset[[outcome]], pf_p2a, family_stage2, seed
  )
  sel_p2a <- .selected_vars(cvfit_p2a_outcome)

  # Retain penalised D x X (non-prespec) and X x X that survived
  retained_dx_only         <- intersect(sel_p2a, setdiff(all_dx_cols, dx_cols_unpen))
  retained_xx_interactions <- intersect(sel_p2a, all_xx_cols)

  retained_dx_only         <- retained_dx_only[nchar(trimws(retained_dx_only)) > 0]
  retained_xx_interactions <- retained_xx_interactions[nchar(trimws(retained_xx_interactions)) > 0]

  retained_dx_interactions <- unique(c(retained_dx_only, retained_xx_interactions))

  message(sprintf(
    "  Penalised pool: D\u00d7X (non-prespec)=%d  X\u00d7X=%d  total=%d",
    length(setdiff(all_dx_cols, dx_cols_unpen)),
    length(all_xx_cols),
    length(dx_cols_penalized)
  ))
  message(sprintf(
    "  Retained: D\u00d7X=%d  X\u00d7X=%d  total=%d",
    length(retained_dx_only),
    length(retained_xx_interactions),
    length(retained_dx_interactions)
  ))

  # ---------------------------------------------------------------------------
  # Part 2b — Exposure interaction selection (SD-based penalty)
  # ---------------------------------------------------------------------------

  message("\n--- Part 2b: Exposure interaction selection (SD-based penalty) ---")

  zx_result   <- .add_interactions(dataset, union_main_effects, instrument, "_iv")
  dataset     <- zx_result$data
  all_zx_cols <- trimws(zx_result$col_names)
  all_zx_cols <- all_zx_cols[!is.na(all_zx_cols) & all_zx_cols != "" & all_zx_cols != "_iv"]

  # Z x prespec: always unpenalised
  zx_prespec_unpen <- paste0(prespec_in_union, "_iv")

  # Z x X also unpenalised for each X where D x X was retained (non-prespec D x X only)
  retained_dx_base  <- trimws(gsub("_d$", "", retained_dx_only))
  retained_dx_base  <- retained_dx_base[!is.na(retained_dx_base) & retained_dx_base != ""]
  retained_zx_unpen <- paste0(retained_dx_base, "_iv")
  retained_zx_unpen <- retained_zx_unpen[nchar(trimws(retained_zx_unpen)) > 0]

  zx_cols_unpen <- unique(trimws(c(zx_prespec_unpen, retained_zx_unpen)))
  zx_cols_unpen <- zx_cols_unpen[!is.na(zx_cols_unpen) & zx_cols_unpen != "" & zx_cols_unpen != "_iv"]

  # Penalised pool: remaining Z x X (not forced) + all X x X
  zx_cols_penalized <- c(setdiff(all_zx_cols, zx_cols_unpen), all_xx_cols)

  p2b_unpenalized <- unique(c(
    union_main_effects, instrument, prespecified_subgroups, zx_cols_unpen
  ))
  p2b_predictors <- unique(c(p2b_unpenalized, zx_cols_penalized))
  p2b_predictors <- p2b_predictors[trimws(p2b_predictors) != "" & !is.na(p2b_predictors)]

  X_p2b  <- .build_matrix(dataset, p2b_predictors)
  pf_p2b <- .penalty_factor_sd(X_p2b, unpenalized = p2b_unpenalized)

  cvfit_p2b_exposure <- .cv_lasso(
    X_p2b, dataset[[treatment]], pf_p2b, family_stage1, seed
  )
  sel_p2b <- .selected_vars(cvfit_p2b_exposure)

  final_exposure_terms <- unique(c(p2b_unpenalized, sel_p2b))

  # -- Corrected message: compute forced counts without double-counting overlap --
  n_prespec_forced   <- length(zx_prespec_unpen)
  n_from_retained_dx <- length(setdiff(retained_zx_unpen, zx_prespec_unpen))
  n_total_forced     <- length(zx_cols_unpen)   # after deduplication
  n_penalized_zx     <- length(setdiff(all_zx_cols, zx_cols_unpen))
  n_penalized_xx     <- length(all_xx_cols)
  n_penalized_total  <- length(zx_cols_penalized)
  n_additional_sel   <- length(intersect(sel_p2b, zx_cols_penalized))

  message(sprintf(
    "  Forced Z\u00d7X: prespec=%d + from retained D\u00d7X=%d => %d (after dedup)",
    n_prespec_forced, n_from_retained_dx, n_total_forced
  ))
  message(sprintf(
    "  Penalised pool: Z\u00d7X=%d  X\u00d7X=%d  total=%d  |  Additional selected: %d",
    n_penalized_zx, n_penalized_xx, n_penalized_total, n_additional_sel
  ))

  # ---------------------------------------------------------------------------
  # Refit exposure model via glm.fit() and compute generalised residuals
  # ---------------------------------------------------------------------------

  X_exposure_refit <- cbind(1, .build_matrix(dataset, final_exposure_terms))

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
  # ---------------------------------------------------------------------------

  message("\n--- Final outcome model: 2SRI correction ---")

  resid_col   <- ".ctrl_resid"
  resid_d_col <- ".ctrl_resid_x_d"

  dataset[[resid_col]]   <- generalised_residuals
  dataset[[resid_d_col]] <- generalised_residuals * dataset[[treatment]]

  final_outcome_terms <- unique(c(
    p2a_unpenalized,
    retained_dx_interactions,
    resid_col,
    resid_d_col
  ))

  final_outcome_X <- cbind(1, .build_matrix(dataset, final_outcome_terms))

  final_outcome_fit <- stats::glm.fit(
    x      = final_outcome_X,
    y      = dataset[[outcome]],
    family = .glm_family(family_stage2)
  )

  message(sprintf(
    "  Final outcome model terms: %d (incl. generalised residual + D\u00d7residual)",
    length(final_outcome_terms)
  ))

  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------

  message("\n========== LASSO IV Selection Summary ==========")
  message(sprintf("  Part 1  union main effects:           %d", length(union_main_effects)))
  message(sprintf("  Part 2a retained D\u00d7X interactions:   %d", length(retained_dx_only)))
  message(sprintf("  Part 2a retained X\u00d7X interactions:   %d", length(retained_xx_interactions)))
  message(sprintf("  Part 2b final exposure terms:         %d", length(final_exposure_terms)))
  message(sprintf("  Final   outcome model terms:          %d", length(final_outcome_terms)))
  message("=================================================\n")

  list(
    union_main_effects        = union_main_effects,
    retained_dx_interactions  = retained_dx_interactions,
    final_exposure_terms      = final_exposure_terms,
    generalised_residuals     = generalised_residuals,
    final_outcome_terms       = final_outcome_terms,
    final_outcome_fit         = final_outcome_fit,
    cvfit_p1_exposure         = cvfit_p1_exposure,
    cvfit_p1_outcome          = cvfit_p1_outcome,
    cvfit_p2a_outcome         = cvfit_p2a_outcome,
    cvfit_p2b_exposure        = cvfit_p2b_exposure
  )
}