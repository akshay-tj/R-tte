# descriptives.R
# Module: unadjusted outcome summaries (Table 2)
#
# Exported functions:
#   table2_outcomes() — stratified TableOne + Wald t-tests for one horizon


#' Table 2: Unadjusted outcomes by treatment arm (one horizon)
#'
#' Produces a stratified descriptive table (via \pkg{tableone}) and Wald-type
#' mean / risk-difference t-tests for a user-specified set of outcomes.
#' Call once per time horizon with explicit column name vectors.
#'
#' @param data A tibble containing outcome columns and the strata column.
#' @param strata_col Character. Name of the binary treatment indicator column
#'   where \code{1} = intervention arm and \code{0} = control arm. Defaults to
#'   \code{"early_surgery"}.
#' @param cont_vars Character vector of continuous outcome column names to
#'   include in the table and t-tests.
#' @param cat_vars Character vector of binary (0/1) outcome column names to
#'   include. Defaults to \code{character(0)} (no categorical variables).
#' @param horizon_label Character. Label printed in section headers to identify
#'   the time horizon, e.g. \code{"90 days"}. Optional but recommended.
#'
#' @return A tibble with one row per outcome variable and columns:
#'   \code{variable}, \code{mean_control}, \code{mean_interv},
#'   \code{difference} (intervention minus control), \code{ci_lower},
#'   \code{ci_upper} (rounded to 2 d.p.), \code{p_value} (rounded to 4 d.p.).
#'   The TableOne is printed as a side-effect.
#'
#' @details
#' ## Usage
#' Call once per time horizon with fully explicit column names:
#'
#' \preformatted{
#'   results_90 <- table2_outcomes(
#'     data          = df,
#'     strata_col    = "early_surgery",
#'     cont_vars     = c("daoh_bypass_surg_90d",
#'                       "daoh_myles_bypass_surg_90d",
#'                       "bypass_surg_los_no",
#'                       "post_bypass_surg_los_no_90d",
#'                       "total_los_no_90d"),
#'     cat_vars      = c("died_post_bypass_surg_90d"),
#'     horizon_label = "90 days"
#'   )
#' }
#'
#' ## Strata column
#' Must already be present in \code{data} — no join is performed here.
#'
#' ## Continuous outcomes
#' Compared with a Welch two-sample t-test (unequal variances). The reported
#' difference is \strong{intervention minus control} (\code{strata = 1} minus
#' \code{strata = 0}).
#'
#' ## Binary outcomes (mortality, readmission)
#' Also compared via Welch t-test on the 0/1 indicator, yielding a
#' \strong{Wald-type risk difference} (intervention risk minus control risk)
#' with asymptotic 95\% CI. Interpret with caution if events < ~10 per arm.
#'
#' @examples
#' \dontrun{
#'   results_90 <- table2_outcomes(
#'     data          = analysis_df,
#'     strata_col    = "early_surgery",
#'     cont_vars     = c("daoh_bypass_surg_90d",
#'                       "daoh_myles_bypass_surg_90d",
#'                       "bypass_surg_los_no",
#'                       "post_bypass_surg_los_no_90d",
#'                       "total_los_no_90d"),
#'     cat_vars      = c("died_post_bypass_surg_90d"),
#'     horizon_label = "90 days"
#'   )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom tableone CreateTableOne
#' @export
table2_outcomes <- function(data,
                            strata_col    = "early_surgery",
                            cont_vars,
                            cat_vars      = character(0),
                            horizon_label = NULL) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    "`data` must be a data frame" =
      is.data.frame(data),
    "`strata_col` must be a single character string" =
      is.character(strata_col) && length(strata_col) == 1,
    "`cont_vars` must be a non-empty character vector" =
      is.character(cont_vars) && length(cont_vars) >= 1,
    "`cat_vars` must be a character vector" =
      is.character(cat_vars)
  )

  if (!strata_col %in% names(data)) {
    stop(sprintf(
      "Strata column '%s' not found in data.", strata_col
    ), call. = FALSE)
  }

  strata_vals <- sort(unique(data[[strata_col]]))
  if (!identical(strata_vals, c(0L, 1L)) &&
      !identical(strata_vals, c(0, 1))) {
    stop(sprintf(
      "'%s' must be a binary 0/1 column; found values: %s",
      strata_col, paste(strata_vals, collapse = ", ")
    ), call. = FALSE)
  }

  all_vars     <- c(cont_vars, cat_vars)
  missing_cols <- setdiff(all_vars, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "The following columns were not found in data:\n  ",
      paste(missing_cols, collapse = "\n  "),
      call. = FALSE
    )
  }

  # ── Header ──────────────────────────────────────────────────────────────────
  label <- if (!is.null(horizon_label)) horizon_label else "outcomes"
  sep   <- strrep("=", 60)
  cat(sprintf(
    "\n%s\n  Table 2: unadjusted %s\n  %s = 1 (intervention)  |  %s = 0 (control)\n%s\n",
    sep, label, strata_col, strata_col, sep
  ))

  # ── TableOne ─────────────────────────────────────────────────────────────────
  tab <- tableone::CreateTableOne(
    vars       = all_vars,
    data       = data,
    factorVars = cat_vars,
    strata     = strata_col,
    test       = FALSE,
    includeNA  = TRUE
  )

  if (length(cat_vars) > 0) {
    cat("\n--- Categorical outcomes ---\n")
    print(tab$CatTable, digits = 2, smd = TRUE)
  }

  cat("\n--- Continuous outcomes ---\n")
  print(tab$ContTable, digits = 2, smd = TRUE)

  # ── t-tests ──────────────────────────────────────────────────────────────────
  # t.test(y ~ x) with x in {0, 1}:
  #   estimate[1] = mean(group 0 = control)
  #   estimate[2] = mean(group 1 = intervention)
  #   diff(estimate) = estimate[2] - estimate[1] = intervention - control  ✓
  #   conf.int is for (group 0 - group 1) = (control - intervention)
  #   Negate and flip conf.int to report CI for (intervention - control).
  cat(sprintf(
    paste0(
      "\n--- Differences: intervention (%s = 1) minus control (%s = 0) ---\n",
      "  Welch two-sample t-test (unequal variances).\n",
      "  Binary outcomes: Wald risk difference (not OR/RR)",
      " — caution if n_events < 10 per arm.\n\n"
    ),
    strata_col, strata_col
  ))

  results <- purrr::map_dfr(all_vars, function(var) {
    test <- t.test(data[[var]] ~ data[[strata_col]])
    tibble::tibble(
      variable     = var,
      mean_control = test$estimate[[1]],
      mean_interv  = test$estimate[[2]],
      difference   = diff(test$estimate),        # intervention - control
      ci_lower     = round(-test$conf.int[[2]], 3),
      ci_upper     = round(-test$conf.int[[1]], 3),
      p_value      = test$p.value
    )
  })

  results
}
