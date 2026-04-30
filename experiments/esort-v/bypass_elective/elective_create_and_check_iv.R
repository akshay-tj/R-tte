# =============================================================================
# create_and_check_iv.R
#
# Purpose: Build the tendency-to-expedite-surgery (TTES) instrumental variable
#          for the elective bypass cohort, then produce IV diagnostic plots
#          (covariate balance, IV stability across centres).
#
# Assumes in environment: nvr_df, elective_cohort, elective_cohort_lookback_only,
#                         INCLUSION_CRITERIA, HIGH_VOL_HOSPS
#
# Outputs:
#   - elective_IV.csv
#   - covariate_balance.png
#   - iv_stability_overall.png
#   - iv_stability_facet.png
# =============================================================================

library(tidyverse)

source("R/create_iv.R")
source("R/plots.R")

# =============================================================================
# Parameters
# =============================================================================

LOOKBACK_DAYS    <- 365L
IV_OUTPUT_PATH   <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_with_ttes.csv"
RESULTS_DIR      <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/"

# =============================================================================
# SECTION 1: Build lookback cohort (pre-study peers for TTES denominator)
# =============================================================================

names(elective_cohort_lookback_only)
elective_cohort_lookback_only <- elective_cohort_lookback_only %>%
  select(STUDY_ID, early_surgery, valid_op_date, NvrHospitalName) %>%
  mutate(
    include_flag = 0L,
    STUDY_ID     = as.character(STUDY_ID)
  ) 

# =============================================================================
# SECTION 2: Prepare study cohort for combining with lookback cohort
# =============================================================================

study_cohort_slim <- elective_cohort %>%
  select(STUDY_ID, early_surgery, valid_op_date, NvrHospitalName) %>%
  mutate(
    include_flag = 1L,
    STUDY_ID     = as.character(STUDY_ID)
  ) %>%
  rename(`Patient:PatientId` = STUDY_ID)

# =============================================================================
# SECTION 3: Calculate TTES instrumental variable
# =============================================================================

combined_df <- bind_rows(elective_cohort_lookback_only, study_cohort_slim)

n_study <- sum(combined_df$include_flag == 1)
message(sprintf("Calculating TTES for %d study patients...", n_study))

combined_df <- combined_df %>%
  rowwise() %>%
  mutate(
    instrumental_variable = if_else(
      include_flag == 1,
      calc_ttes(
        combined_df,
        `Patient:PatientId`,
        NvrHospitalName,
        valid_op_date,
        lookback_days = LOOKBACK_DAYS, 
        date_col      = "valid_op_date"
      ),
      NA_real_
    )
  ) %>%
  ungroup()

# Impute missing IV values (patients with no peers in lookback window)
n_missing <- sum(is.na(combined_df$instrumental_variable) & combined_df$include_flag == 1)
message(sprintf(
  "IV missing for %d of %d study patients (%.1f%%) — imputing with median.",
  n_missing, n_study, 100 * n_missing / n_study
))

median_iv <- median(combined_df$instrumental_variable, na.rm = TRUE)
combined_df <- combined_df %>%
  mutate(
    instrumental_variable = if_else(
      is.na(instrumental_variable),
      median_iv,
      instrumental_variable
    )
  )

iv_df <- combined_df %>%
  filter(include_flag == 1) %>%
  select(`Patient:PatientId`, instrumental_variable) %>%
  rename(STUDY_ID = `Patient:PatientId`) %>%
  mutate(STUDY_ID = as.character(STUDY_ID))

write_csv(iv_df, IV_OUTPUT_PATH)
message(sprintf("IV written to: %s", IV_OUTPUT_PATH))

# =============================================================================
# SECTION 4: Join IV to study cohort
# =============================================================================

analysis_df <- elective_cohort %>%
  left_join(iv_df, by = "STUDY_ID") %>%
  rename(
    study_id           = STUDY_ID,
    age_at_surgery     = `Patient:AgeAtSurgery`,
    gender_code        = `Patient:GenderCode`,
    nvr_procedure_date = `NvrEpisode:ProcedureStartDate`,
    nvr_hospital_name  = NvrHospitalName
  )

# =============================================================================
# SECTION 5: Covariate balance plot
# =============================================================================

balance_labels <- c(
  age_at_surgery                     = "Age",
  gender_code_1                      = "Male",
  gender_code_2                      = "Female",
  rcs_ch_cat_None                    = "Charlson - No comorbidities",
  rcs_ch_cat_One                     = "Charlson - One comorbidity",
  rcs_ch_cat_Two                     = "Charlson - Two comorbidities",
  `rcs_ch_cat_Three +`               = "Charlson - Three or more comorbidities",
  `scarf_cat_Fit (0/1)`              = "SCARF - Fit",
  `scarf_cat_Mild Frailty (2/3)`     = "SCARF - Mild Frailty",
  `scarf_cat_Moderate Frailty (4/5)` = "SCARF - Moderate Frailty",
  `scarf_cat_Severe Frailty (6+)`    = "SCARF - Severe Frailty"
)

balance_plot <- plot_covariate_balance(
  df               = analysis_df,
  iv_col           = "instrumental_variable",
  covariates       = c("age_at_surgery", "gender_code", "rcs_ch_cat", "scarf_cat"),
  categorical_cols = c("gender_code", "rcs_ch_cat", "scarf_cat"),
  labels           = balance_labels,
  x_label          = "Tendency to Expedite Surgery (Decile)",
  y_label          = "Mean baseline covariate rescaled by its SD",
  title            = "Covariate Balance Plot"
)

balance_plot

# =============================================================================
# SECTION 6: IV stability plots
# =============================================================================

# Restrict to hospitals with > 2 years of data
stability_df <- analysis_df %>%
  mutate(year = format(nvr_procedure_date, "%Y")) %>%
  group_by(nvr_hospital_name) %>%
  filter(n_distinct(year) > 2) %>%
  ungroup()

stability_plots <- plot_iv_stability(
  df              = stability_df,
  iv_col          = "instrumental_variable",
  centre_col      = "nvr_hospital_name",
  date_col        = "nvr_procedure_date",
  y_label         = "Mean TTE Value",
  x_label_overall = "Hospitals"
)

stability_plots$overall
stability_plots$facet

# =============================================================================
# SECTION 7: Save plots
# =============================================================================

ggsave(
  file.path(RESULTS_DIR, "covariate_balance.png"),
  plot   = balance_plot,
  width  = 12,
  height = 7,
  dpi    = 300
)
ggsave(
  file.path(RESULTS_DIR, "iv_stability_overall.png"),
  plot   = stability_plots$overall,
  width  = 8,
  height = 5,
  dpi    = 300
)
ggsave(
  file.path(RESULTS_DIR, "iv_stability_facet.png"),
  plot   = stability_plots$facet,
  width  = 14,
  height = 10,
  dpi    = 300
)