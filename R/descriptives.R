# descriptives.R
# Module: descriptive tables (Table 1 baseline characteristics, Table 2 outcomes)
#
# Exported functions:
#   descriptive_table() — stratified TableOne + optional Wald t-tests


#' Stratified descriptive table with optional t-tests
#'
#' Produces a stratified descriptive table (via \pkg{tableone}) and optionally
#' Wald-type mean / risk-difference t-tests. Covers both Table 1 (baseline
#' characteristics) and Table 2 (unadjusted outcomes) use cases — controlled
#' by the \code{ttest} argument.
#'
#' @param data A tibble containing all variable and strata columns.
#' @param strata_col Character. Name of the binary treatment indicator column
#'   where \code{1} = intervention arm and \code{0} = control arm. Defaults to
#'   \code{"early_surgery"}.
#' @param cont_vars Character vector of continuous variable column names.
#' @param cat_vars Character vector of categorical variable column names —
#'   covers both binary (0/1) flags \strong{and} multi-level factors.
#'   Defaults to \code{character(0)}.
#'   \strong{Important:} multi-level variables must be pre-converted to
#'   \code{factor} in \code{data} before calling this function, with levels
#'   set in the desired display order. Binary 0/1 columns do not need
#'   pre-conversion. Example:
#'   \preformatted{
#'     data <- data %>%
#'       mutate(
#'         smoking = factor(smoking,
#'                          levels = c("never", "ex", "current")),
#'         asa     = factor(asa, levels = 1:4)
#'       )
#'   }
#' @param ttest Logical. Should Wald-type t-tests be run and returned?
#'   Defaults to \code{FALSE}. Set to \code{TRUE} for Table 2 outcome
#'   comparisons. When \code{FALSE}, returns \code{invisible(NULL)}.
#' @param label Character. Label printed in the section header, e.g.
#'   \code{"90-day outcomes"} or \code{"baseline characteristics"}.
#'   Optional but recommended.
#' @param smd Logical. Should standardised mean differences be printed?
#'   Defaults to \code{TRUE}.
#'
#' @return When \code{ttest = TRUE}: a tibble with one row per variable in
#'   \code{cont_vars} and \code{cat_vars}, with columns \code{variable},
#'   \code{mean_control}, \code{mean_interv}, \code{difference}
#'   (intervention minus control), \code{ci_lower}, \code{ci_upper},
#'   \code{p_value}. The \pkg{tableone} table is printed as a side-effect.
#'   When \code{ttest = FALSE}: \code{invisible(NULL)} — the
#'   \pkg{tableone} table is printed as a side-effect only.
#'
#' @details
#' ## Continuous outcomes (when \code{ttest = TRUE})
#' Compared with a Welch two-sample t-test (unequal variances). The reported
#' difference is \strong{intervention minus control} (\code{strata = 1} minus
#' \code{strata = 0}).
#'
#' ## Binary outcomes (when \code{ttest = TRUE})
#' Also compared via Welch t-test on the 0/1 indicator, yielding a
#' \strong{Wald-type risk difference} (intervention risk minus control risk)
#' with asymptotic 95\% CI. Interpret with caution if events < ~10 per arm.
#'
#' ## Multi-level categorical variables
#' \pkg{tableone} renders multi-level factors as grouped count/\% rows.
#' Pre-convert to \code{factor} with explicit levels before calling — see
#' \code{cat_vars} parameter documentation above.
#'
#' ## t-tests and multi-level factors
#' When \code{ttest = TRUE}, t-tests run over \code{cont_vars} and
#' \code{cat_vars}. Multi-level factor columns will produce nonsensical
#' t-test results — only pass binary 0/1 columns in \code{cat_vars} when
#' \code{ttest = TRUE}.
#'
#' @examples
#' \dontrun{
#'   # Table 1 — no t-tests
#'   descriptive_table(
#'     data       = cohort_df,
#'     strata_col = "early_surgery",
#'     cont_vars  = c("age_at_surgery"),
#'     cat_vars   = c("sex_female", "smoking", "asa"),
#'     ttest      = FALSE,
#'     label      = "baseline characteristics"
#'   )
#'
#'   # Table 2 — with t-tests (binary cat_vars only)
#'   results_90 <- descriptive_table(
#'     data       = outcomes_df,
#'     strata_col = "early_surgery",
#'     cont_vars  = c("daoh_bypass_surg_90d", "total_los_no_90d"),
#'     cat_vars   = c("died_post_bypass_surg_90d"),
#'     ttest      = TRUE,
#'     label      = "90-day outcomes"
#'   )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom tableone CreateTableOne
#' @export
descriptive_table <- function(data,
                              strata_col,
                              cont_vars,
                              cat_vars   = character(0),
                              ttest      = FALSE,
                              label      = NULL,
                              smd        = TRUE) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    "`data` must be a data frame"                          =
      is.data.frame(data),
    "`strata_col` must be a single character string"       =
      is.character(strata_col) && length(strata_col) == 1L,
    "`cont_vars` must be a non-empty character vector"     =
      is.character(cont_vars) && length(cont_vars) >= 1L,
    "`cat_vars` must be a character vector"                =
      is.character(cat_vars),
    "`ttest` must be a single logical"                     =
      is.logical(ttest) && length(ttest) == 1L,
    "`smd` must be a single logical"                       =
      is.logical(smd) && length(smd) == 1L
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

  # ── Warn if ttest = TRUE and multi-level factors present ────────────────────
  if (ttest && length(cat_vars) > 0) {
    multilevel <- cat_vars[sapply(cat_vars, function(v) {
      is.factor(data[[v]]) && nlevels(data[[v]]) > 2L
    })]
    if (length(multilevel) > 0) {
      warning(
        "ttest = TRUE but the following cat_vars are multi-level factors — ",
        "t-test results for these will be meaningless:\n  ",
        paste(multilevel, collapse = "\n  "),
        call. = FALSE
      )
    }
  }

  # ── Header ──────────────────────────────────────────────────────────────────
  label_str <- if (!is.null(label)) label else "descriptive table"
  sep       <- strrep("=", 60)
  cat(sprintf(
    "\n%s\n  %s\n  %s = 1 (intervention)  |  %s = 0 (control)\n%s\n",
    sep, label_str, strata_col, strata_col, sep
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
    cat("\n--- Categorical variables ---\n")
    print(tab$CatTable, digits = 2, smd = smd)
  }

  cat("\n--- Continuous variables ---\n")
  print(tab$ContTable, digits = 2, smd = smd)

  # ── t-tests (Table 2 only) ───────────────────────────────────────────────────
  if (!ttest) return(invisible(NULL))

  cat(sprintf(
    paste0(
      "\n--- Differences: intervention (%s = 1) minus control (%s = 0) ---\n",
      "  Welch two-sample t-test (unequal variances).\n",
      "  Binary outcomes: Wald risk difference (not OR/RR)",
      " — caution if n_events < 10 per arm.\n\n"
    ),
    strata_col, strata_col
  ))

  purrr::map_dfr(all_vars, function(var) {
    test <- t.test(data[[var]] ~ data[[strata_col]])
    tibble::tibble(
      variable     = var,
      mean_control = test$estimate[[1]],
      mean_interv  = test$estimate[[2]],
      difference   = diff(test$estimate),
      ci_lower     = round(-test$conf.int[[2]], 3),
      ci_upper     = round(-test$conf.int[[1]], 3),
      p_value      = test$p.value
    )
  })
}