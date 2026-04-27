# non_elective_forest_plots.R
#
# Produces forest plots for model2 (heterogeneous 2SRI) subgroup results.
# One plot per outcome x time horizon.
#
# Inputs:
#   - Bootstrap results CSVs: bootstrap_results_{outcome}_{version}_{H}d.csv
#   - Non-elective DTAs:      non_elective_{H}d.dta [from lasso outputs]
#
# Outputs:
#   - One PNG per outcome x horizon: forest_{outcome}_{H}d_{version}.png
#
# Note:
#   - No existence check on DTA files before calling get_subgroup_ns() — the
#     loop will error if a file is missing rather than skip gracefully.
#   - surgyr_2015 is derived as "not any of 2016-2023" — will be incorrect if
#     the data contain surgery years outside 2015-2023.

library(dplyr)
library(ggplot2)
library(haven)
library(stringr)
library(purrr)
library(tidyr)

source("R/plots.R")  # provides make_forest_plot()

# =============================================================================
# CONFIGURATION
# =============================================================================

DATA_DIR    <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/lasso_outputs/"
RESULTS_DIR <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/clinical_effectiveness_results/"

TIME_HORIZONS <- c(90, 180, 365)
VERSION       <- "model2"

OUTCOMES <- list(
  list(
    col_pattern     = "daoh_bypass_surg_{H}d",
    pred_type       = "ystar",
    label           = "DAOH",
    favors_positive = TRUE
  ),
  list(
    col_pattern     = "total_los_no_{H}d",
    pred_type       = "ystar",
    label           = "Total LOS",
    favors_positive = FALSE
  ),
  list(
    col_pattern     = "died_post_bypass_surg_{H}d",
    pred_type       = "pr",
    label           = "Mortality",
    favors_positive = FALSE,
    as_percentage   = TRUE
  ),
  list(
    col_pattern     = "readmit_post_bypass_surg_{H}d",
    pred_type       = "pr",
    label           = "Post-surgery readmission",
    favors_positive = FALSE,
    as_percentage   = TRUE
  ),
  list(
    col_pattern     = "ilr_{H}d",
    pred_type       = "pr",
    label           = "Index limb revascularisation",
    favors_positive = FALSE,
    as_percentage   = TRUE
  ),
  list(
    col_pattern     = "ilma_{H}d",
    pred_type       = "pr",
    label           = "Index limb major amputation",
    favors_positive = FALSE,
    as_percentage   = TRUE
  )
)

SUBGROUP_DEFS <- list(
  list(var = "All",              label = "Overall", group = "Base case"),
  list(var = "gender_F",         label = "Female",  group = "Sex"),
  list(var = "gender_M",         label = "Male",    group = "Sex"),
  list(var = "fontaine_4",       label = "IV",      group = "Fontaine"),
  list(var = "fontaine_not4",    label = "III",     group = "Fontaine"),
  list(var = "comorbidity_1",    label = "Yes",     group = "Diabetes"),
  list(var = "comorbidity_not1", label = "No",      group = "Diabetes"),
  list(var = "krt_yn",           label = "Yes",     group = "KRT"),
  list(var = "krt_no",           label = "No",      group = "KRT"),
  list(var = "surgyr_2015",      label = "2015",    group = "Year of surgery"),
  list(var = "surgyr_2016",      label = "2016",    group = "Year of surgery"),
  list(var = "surgyr_2017",      label = "2017",    group = "Year of surgery"),
  list(var = "surgyr_2018",      label = "2018",    group = "Year of surgery"),
  list(var = "surgyr_2019",      label = "2019",    group = "Year of surgery"),
  list(var = "surgyr_2020",      label = "2020",    group = "Year of surgery"),
  list(var = "surgyr_2021",      label = "2021",    group = "Year of surgery"),
  list(var = "surgyr_2022",      label = "2022",    group = "Year of surgery"),
  list(var = "surgyr_2023",      label = "2023",    group = "Year of surgery")
)

# =============================================================================
# HELPERS
# Note: these are experiment-specific helpers, not package-level functions.
# =============================================================================

# Compute per-subgroup Ns from a DTA file.
# Complement variables (gender_M, fontaine_not4, etc.) are derived manually
# here — this is intentional for this experiment and not worth generalising.
get_subgroup_ns <- function(dta_path, subgroup_defs) {
  dta <- haven::read_dta(dta_path) %>%
    mutate(
      All              = 1L,
      gender_M         = 1L - gender_F,
      fontaine_not4    = 1L - fontaine_4,
      comorbidity_not1 = 1L - comorbidity_1,
      krt_no           = 1L - krt_yn,
      # surgyr_2015 derived as "not any of 2016-2023" — see file-level assumption
      surgyr_2015      = as.integer(
        surgyr_2016 == 0 & surgyr_2017 == 0 & surgyr_2018 == 0 &
        surgyr_2019 == 0 & surgyr_2020 == 0 & surgyr_2021 == 0 &
        surgyr_2022 == 0 & surgyr_2023 == 0
      )
    )

  purrr::map_dfr(subgroup_defs, function(s) {
    tibble::tibble(
      var = s$var,
      n   = sum(dta[[s$var]] == 1L, na.rm = TRUE)
    )
  })
}

# Parse a bootstrap results CSV and return point estimates + CIs for one
# outcome x horizon combination, filtered to the subgroups in subgroup_defs.
get_estimates <- function(results_dir, outcome_col, horizon, version,
                          pred_type, subgroup_defs) {
  csv_path <- file.path(
    results_dir,
    sprintf("bootstrap_results_%s_%s_%dd.csv", outcome_col, version, horizon)
  )

  if (!file.exists(csv_path)) {
    warning("File not found: ", csv_path)
    return(NULL)
  }

  results <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Keep treatment effect rows only — pattern: m{pred_type}_{subgroup}Out{N}
  # Excludes counterfactual mean rows prefixed m0/m1
  effect_rows <- results %>%
    filter(
      str_detect(stat, paste0("^m", pred_type, "_")),
      !str_detect(stat, paste0("^m[01]", pred_type, "_"))
    )

  if (nrow(effect_rows) == 0) {
    warning("No rows found for pred_type '", pred_type, "' in ", csv_path)
    return(NULL)
  }

  subgroup_vars <- purrr::map_chr(subgroup_defs, "var")

  effect_rows %>%
    mutate(
      var = stat %>%
        str_remove(paste0("^m", pred_type, "_")) %>%
        str_remove("Out\\d+$")
    ) %>%
    filter(var %in% subgroup_vars) %>%
    select(var, est, ci_lo, ci_hi)
}

# =============================================================================
# MAIN
# =============================================================================

for (h in TIME_HORIZONS) {

  dta_path <- file.path(DATA_DIR, sprintf("non_elective_%dd.dta", h))
  ns       <- get_subgroup_ns(dta_path, SUBGROUP_DEFS)

  for (outcome_def in OUTCOMES) {

    outcome_col   <- gsub("\\{H\\}", h, outcome_def$col_pattern)
    outcome_label <- outcome_def$label
    pred_type     <- outcome_def$pred_type

    message(sprintf("Plotting: %s | %dd", outcome_label, h))

    estimates <- get_estimates(
      RESULTS_DIR, outcome_col, h, VERSION, pred_type, SUBGROUP_DEFS
    )

    if (is.null(estimates)) next

    subgroup_meta <- purrr::map_dfr(SUBGROUP_DEFS, function(s) {
      tibble::tibble(var = s$var, label = s$label, group = s$group)
    })

    plot_data <- estimates %>%
      left_join(ns,            by = "var") %>%
      left_join(subgroup_meta, by = "var")

    # x_limits derived from observed CI range — see file-level assumption
    x_range  <- range(c(plot_data$ci_lo, plot_data$ci_hi), na.rm = TRUE)
    if (isTRUE(outcome_def$as_percentage)) x_range <- x_range * 100
    x_limits <- c(
      x_range[1] - diff(x_range) * 0.15,
      x_range[2] + diff(x_range) * 0.15
    )

    p <- make_forest_plot(
      plot_data       = plot_data,
      outcome_label   = outcome_label,
      horizon         = h,
      x_limits        = x_limits,
      favors_positive = outcome_def$favors_positive,
      subgroup_defs   = SUBGROUP_DEFS, 
      as_percentage   = isTRUE(outcome_def$as_percentage)
    )

    out_path <- file.path(
      RESULTS_DIR,
      sprintf("forest_%s_%dd_%s.png", outcome_col, h, VERSION)
    )

    n_rows <- length(SUBGROUP_DEFS) +
      length(unique(purrr::map_chr(SUBGROUP_DEFS, "group")))

    ggsave(out_path, p, width = 12, height = n_rows * 0.35 + 2, dpi = 300)
    message("Saved: ", out_path)
  }
}