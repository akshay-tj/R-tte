# krt_interaction_wald.R
#
# Wald test for treatment effect heterogeneity by KRT status.
#
# For each outcome x horizon combination, reads the bootstrap results CSV,
# extracts point estimates and SEs for krt_yn and krt_no, then tests whether
# the difference is statistically significant.
#
# Method: Z = (ATE_krt_yes - ATE_krt_no) / sqrt(SE_yes^2 + SE_no^2)
# Assumes independence between subgroups (valid: krt_yn and krt_no are
# mutually exclusive complements derived from the same sample).
#
# Outputs:
#   - Printed results table
#   - krt_interaction_wald_results.csv

library(dplyr)
library(stringr)
library(purrr)

# =============================================================================
# CONFIGURATION — edit here
# =============================================================================

RESULTS_DIR <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/clinical_effectiveness_results/"
VERSION     <- "model2"
TIME_HORIZONS <- c(90, 180, 365)

# Outcome definitions: col_pattern (with {H} placeholder) and pred_type
OUTCOMES <- list(
  list(
    label       = "DAOH",
    col_pattern = "daoh_bypass_surg_{H}d",
    pred_type   = "ystar"
  ),
  list(
    label       = "Total LOS (day-stays excl.)",
    col_pattern = "total_los_no_{H}d",
    pred_type   = "ystar"
  ),
  list(
    label       = "Post-bypass LOS (day-stays excl.)",
    col_pattern = "post_bypass_surg_los_no_{H}d",
    pred_type   = "ystar"
  ),
  list(
    label       = "Readmission",
    col_pattern = "readmit_post_bypass_surg_{H}d",
    pred_type   = "pr"
  ),
  list(
    label       = "Mortality",
    col_pattern = "died_post_bypass_surg_{H}d",
    pred_type   = "pr"
  ),
  list(
    label       = "ILR",
    col_pattern = "ilr_{H}d",
    pred_type   = "pr"
  ),
  list(
    label       = "ILMA",
    col_pattern = "ilma_{H}d",
    pred_type   = "pr"
  )
)

# =============================================================================
# HELPERS
# =============================================================================

#' Read bootstrap CSV and extract est + se for a named subgroup variable
#' matching a given pred_type.
#'
#' Stat name pattern in CSV: m{pred_type}_{subgroup}Out{N}
#' We match on prefix m{pred_type}_{subgroup} and exclude m0/m1 prefixes
#' (counterfactual means).
.extract_subgroup <- function(csv_path, pred_type, subgroup_var) {

  if (!file.exists(csv_path)) {
    warning("File not found: ", csv_path)
    return(NULL)
  }

  results <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Treatment effect rows: starts with m{pred_type}_{subgroup}Out
  # Exclude counterfactual mean rows starting with m0 or m1
  pattern_keep    <- paste0("^m", pred_type, "_", subgroup_var, "Out")
  pattern_exclude <- paste0("^m[01]", pred_type, "_")

  row <- results %>%
    filter(
      str_detect(stat, pattern_keep),
      !str_detect(stat, pattern_exclude)
    )

  if (nrow(row) == 0) {
    warning(sprintf(
      "No matching row for pred_type='%s', subgroup='%s' in %s",
      pred_type, subgroup_var, basename(csv_path)
    ))
    return(NULL)
  }

  if (nrow(row) > 1) {
    warning(sprintf(
      "%d rows matched for pred_type='%s', subgroup='%s' in %s — using first",
      nrow(row), pred_type, subgroup_var, basename(csv_path)
    ))
    row <- row[1, ]
  }

  list(est = row$est, se = row$se)
}

#' Wald test for difference between two independent estimates.
#' Returns a one-row tibble.
.wald_test <- function(est_yes, se_yes, est_no, se_no,
                       outcome_label, horizon) {
  diff    <- est_yes - est_no
  se_diff <- sqrt(se_yes^2 + se_no^2)
  z       <- diff / se_diff
  p       <- 2 * pnorm(-abs(z))

  tibble::tibble(
    outcome      = outcome_label,
    horizon_days = horizon,
    est_krt_yes  = est_yes,
    se_krt_yes   = se_yes,
    est_krt_no   = est_no,
    se_krt_no    = se_no,
    difference   = diff,
    se_diff      = se_diff,
    z_stat       = z,
    p_value      = p
  )
}

# =============================================================================
# MAIN
# =============================================================================

results <- purrr::map_dfr(OUTCOMES, function(outcome_def) {
  purrr::map_dfr(TIME_HORIZONS, function(h) {

    outcome_col <- str_replace(outcome_def$col_pattern, "\\{H\\}", as.character(h))
    pred_type   <- outcome_def$pred_type

    csv_path <- file.path(
      RESULTS_DIR,
      sprintf("bootstrap_results_%s_%s_%dd.csv", outcome_col, VERSION, h)
    )

    krt_yes <- .extract_subgroup(csv_path, pred_type, "krt_yn")
    krt_no  <- .extract_subgroup(csv_path, pred_type, "krt_no")

    if (is.null(krt_yes) || is.null(krt_no)) {
      warning(sprintf("Skipping %s | %dd — could not extract both subgroups", outcome_col, h))
      return(NULL)
    }

    .wald_test(
      est_yes       = krt_yes$est,
      se_yes        = krt_yes$se,
      est_no        = krt_no$est,
      se_no         = krt_no$se,
      outcome_label = outcome_def$label,
      horizon       = h
    )
  })
})

# =============================================================================
# OUTPUT
# =============================================================================

results_fmt <- results %>%
  mutate(
    across(c(est_krt_yes, se_krt_yes, est_krt_no, se_krt_no,
             difference, se_diff, z_stat), ~ round(.x, 4)),
    p_value = round(p_value, 4),
    p_label = case_when(
      p_value < 0.001 ~ "<0.001",
      p_value < 0.01  ~ "<0.01",
      p_value < 0.05  ~ "<0.05",
      TRUE            ~ as.character(p_value)
    )
  )

print(results_fmt, n = Inf)

output_path <- file.path(RESULTS_DIR, "krt_interaction_wald_results.csv")
write.csv(results_fmt, output_path, row.names = FALSE)
message("Results written to: ", output_path)