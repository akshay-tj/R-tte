################################################################################
# elective_sens_trimmed_window_setup.R
#
# Sensitivity analysis - Trimmed window - for elective bypass cohort.
#
# Sensitivity treatment definition:
#   early_surg_sens = 1 if daystosurgery is 1–14 days inclusive
#   early_surg_sens = 0 if daystosurgery is 15–28 days inclusive
#   exclude patients outside 1–28 days
#
# This script:
#   1. Loads the existing elective baseline, lookback-only, and outcomes datasets
#   2. Re-defines early/later surgery for the sensitivity analysis
#   3. Recalculates the TTES instrumental variable within the sensitivity cohort
#   4. Uses previously generated LASSO IV covariate selection
#   5. Exports one Stata-ready .dta and one globals .csv per time horizon
################################################################################

library(dplyr)
library(readr)
library(haven)
library(purrr)
library(stringr)

source("R/create_iv.R")
source("R/lasso.R")

# =============================================================================
# PATHS — edit here only
# =============================================================================

ELECTIVE_COHORT_BASELINE_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_with_confounders.csv"
ELECTIVE_COHORT_LOOKBACK_ONLY_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_lookback_only.csv"

# This should be the wide outcome dataset used by the existing elective_run_lasso.R
ELECTIVE_OUTCOMES_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_with_outcomes.csv"

OUTPUT_DIR <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/sensitivity_trimmed_window/"

IV_OUTPUT_PATH <- file.path(
  OUTPUT_DIR,
  "elective_sens_bypass_study_participants_with_ttes.csv"
)

PREVIOUS_LASSO_OUTPUT_DIR <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/lasso_outputs/"
SENS_LASSO_OUTPUT_DIR <- file.path(OUTPUT_DIR, "lasso_outputs")

LOOKBACK_DAYS <- 365L

# =============================================================================
# PARAMETERS
# =============================================================================

TIME_HORIZONS <- c(90, 180, 365)

OUTCOMES <- list(
  daoh = "daoh_bypass_surg_{H}d",
  daoh_myles = "daoh_myles_bypass_surg_{H}d",
  total_los = "total_los_no_{H}d",
  readmission = "readmit_post_bypass_surg_{H}d",
  mortality = "died_post_bypass_surg_{H}d",
  post_bypass_surg_los_no = "post_bypass_surg_los_no_{H}d",
  ilr = "ilr_{H}d",
  ilma = "ilma_{H}d"
)

OUTCOME_FAMILIES <- c(
  daoh = "gaussian",
  daoh_myles = "gaussian",
  total_los = "gaussian",
  readmission = "binomial",
  mortality = "binomial",
  post_bypass_surg_los_no = "gaussian",
  ilr = "binomial",
  ilma = "binomial"
)

# Prespecified subgroups — always forced in, unpenalised
not_penalized <- c(
  "ageatsurgery",
  "gender_F",
  "fontaine_4",
  "comorbidity_1",
  "krt_yn"
)

# Penalised main effects
penalized_vars <- c(
  "age_sq",
  "smoking_ex",
  "smoking_current",
  "asa_2",
  "asa_3",
  "asa_4",
  "scarf_mild",
  "scarf_moderate",
  "scarf_severe",
  "rcs_two",
  "rcs_threeplus",
  "imd_q2",
  "imd_q3",
  "imd_q4",
  "imd_q5",
  "amp_tissueloss",
  "comorbidity_2",
  "comorbidity_3",
  "comorbidity_4",
  "comorbidity_5",
  "comorbidity_6",
  "comorbidity_7",
  "comorbidity_8",
  "medication_1",
  "medication_2",
  "medication_3",
  "medication_4",
  "surgyr_2017",
  "surgyr_2018",
  "surgyr_2019",
  "surgyr_2020",
  "surgyr_2021",
  "surgyr_2022",
  "surgyr_2023",
  "covid_period",
  "postcovid_period"
)

# =============================================================================
# HELPERS
# =============================================================================

define_sensitivity_treatment <- function(df) {

  df %>%
    mutate(
      STUDY_ID = as.character(STUDY_ID),
      early_surg_sens = case_when(
        daystosurgery >= 1 & daystosurgery <= 14 ~ 1L,
        daystosurgery >= 15 & daystosurgery <= 28 ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(!is.na(early_surg_sens))
}

# =============================================================================
# LOAD DATA
# =============================================================================

elective_cohort_baseline <- read_csv(ELECTIVE_COHORT_BASELINE_DF_PATH)
elective_cohort_lookback_only_raw <- read_csv(ELECTIVE_COHORT_LOOKBACK_ONLY_DF_PATH)
elective_outcomes <- read_csv(ELECTIVE_OUTCOMES_DF_PATH)

# =============================================================================
# DEFINE SENSITIVITY COHORT
# =============================================================================

elective_cohort_sens <- define_sensitivity_treatment(
  elective_cohort_baseline
)

elective_cohort_sens_lookback_only <- define_sensitivity_treatment(
  elective_cohort_lookback_only_raw
)

group_counts <- elective_cohort_sens %>%
  count(early_surg_sens, name = "n") %>%
  mutate(
    group = if_else(
      early_surg_sens == 1L,
      "Early surgery, 1–14 days",
      "Later surgery, 15–28 days"
    ),
    pct = 100 * n / sum(n)
  ) %>%
  select(group, early_surg_sens, n, pct)

print(group_counts)

# =============================================================================
# RECALCULATE TTES IV
# =============================================================================

lookback_slim <- elective_cohort_sens_lookback_only %>%
  select(STUDY_ID, early_surg_sens, valid_op_date, NvrHospitalName) %>%
  mutate(
    include_flag = 0L,
    `Patient:PatientId` = STUDY_ID,
    early_surgery = early_surg_sens
  ) %>%
  select(`Patient:PatientId`, early_surgery, valid_op_date, NvrHospitalName, include_flag)

study_slim <- elective_cohort_sens %>%
  select(STUDY_ID, early_surg_sens, valid_op_date, NvrHospitalName) %>%
  mutate(
    include_flag = 1L,
    `Patient:PatientId` = STUDY_ID,
    early_surgery = early_surg_sens
  ) %>%
  select(`Patient:PatientId`, early_surgery, valid_op_date, NvrHospitalName, include_flag)

combined_df <- bind_rows(lookback_slim, study_slim)

n_study <- sum(combined_df$include_flag == 1L)

message(sprintf("Calculating sensitivity TTES for %d study patients...", n_study))

combined_df <- combined_df %>%
  rowwise() %>%
  mutate(
    instrumental_variable = if_else(
      include_flag == 1L,
      calc_ttes(
        combined_df,
        `Patient:PatientId`,
        NvrHospitalName,
        valid_op_date,
        lookback_days = LOOKBACK_DAYS,
        date_col = "valid_op_date"
      ),
      NA_real_
    )
  ) %>%
  ungroup()

n_missing <- sum(
  is.na(combined_df$instrumental_variable) &
    combined_df$include_flag == 1L
)

message(sprintf(
  "IV missing for %d of %d study patients (%.1f%%) — imputing with median.",
  n_missing,
  n_study,
  100 * n_missing / n_study
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
  filter(include_flag == 1L) %>%
  transmute(
    STUDY_ID = as.character(`Patient:PatientId`),
    instrumental_variable
  )

write_csv(iv_df, IV_OUTPUT_PATH)
message(sprintf("Sensitivity IV written to: %s", IV_OUTPUT_PATH))

# =============================================================================
# REUSE EXISTING LASSO OUTPUTS
# =============================================================================

for (h in TIME_HORIZONS) {

  message(sprintf("Preparing sensitivity Stata files for %dd horizon...", h))

  old_dta_path <- file.path(
    PREVIOUS_LASSO_OUTPUT_DIR,
    sprintf("elective_%dd.dta", h)
  )

  old_globals_path <- file.path(
    PREVIOUS_LASSO_OUTPUT_DIR,
    sprintf("elective_%dd_globals.csv", h)
  )

  new_dta_path <- file.path(
    SENS_LASSO_OUTPUT_DIR,
    sprintf("elective_sens_%dd.dta", h)
  )

  new_globals_path <- file.path(
    SENS_LASSO_OUTPUT_DIR,
    sprintf("elective_sens_%dd_globals.csv", h)
  )

  # ---------------------------------------------------------------------------
  # Load original Stata-ready LASSO dataset
  # ---------------------------------------------------------------------------

  dta_df <- read_dta(old_dta_path) %>%
    mutate(study_id = as.character(study_id))

  # ---------------------------------------------------------------------------
  # Sensitivity ID/treatment/IV lookup
  # ---------------------------------------------------------------------------

  sens_lookup <- elective_cohort_sens %>%
    transmute(
      study_id = as.character(STUDY_ID),
      early_surg_sens = as.integer(early_surg_sens)
    ) %>%
    left_join(
      iv_df %>%
        transmute(
          study_id = as.character(STUDY_ID),
          instrumental_variable_sens = instrumental_variable
        ),
      by = "study_id"
    )

  # ---------------------------------------------------------------------------
  # Subset original .dta to sensitivity participants and replace treatment + IV
  # ---------------------------------------------------------------------------

  dta_sens <- dta_df %>%
    inner_join(sens_lookup, by = "study_id") %>%
    mutate(
      early_surgery = early_surg_sens,
      instrumental_variable = instrumental_variable_sens
    ) %>%
    select(-instrumental_variable_sens)

  # Keep an explicit copy as well; useful for checking / clarity
  if (!"early_surg_sens" %in% names(dta_sens)) {
    stop("early_surg_sens not found after joining sensitivity lookup.")
  }

  # ---------------------------------------------------------------------------
  # Important: regenerate any instrument interaction columns
  #
  # The previous .dta may contain z_x_stage1_* columns based on the old IV.
  # These MUST be updated to use the recalculated sensitivity IV.
  # ---------------------------------------------------------------------------

  zx_cols <- grep("^z_x_stage1_", names(dta_sens), value = TRUE)

  for (zx_col in zx_cols) {
    base_var <- sub("^z_x_stage1_", "", zx_col)

    if (!base_var %in% names(dta_sens)) {
      warning(sprintf(
        "Could not regenerate %s because base variable %s is missing.",
        zx_col,
        base_var
      ))
      next
    }

    dta_sens[[zx_col]] <- dta_sens[["instrumental_variable"]] * dta_sens[[base_var]]
  }

  # ---------------------------------------------------------------------------
  # Optional checks
  # ---------------------------------------------------------------------------

  message(sprintf(
    "Horizon %dd: original N = %d; sensitivity N = %d",
    h,
    nrow(dta_df),
    nrow(dta_sens)
  ))

  print(
    dta_sens %>%
      count(early_surgery, name = "n") %>%
      mutate(pct = 100 * n / sum(n))
  )

  if (any(is.na(dta_sens$instrumental_variable))) {
    warning(sprintf(
      "Horizon %dd has %d missing IV values.",
      h,
      sum(is.na(dta_sens$instrumental_variable))
    ))
  }

  # ---------------------------------------------------------------------------
  # Write sensitivity .dta
  # ---------------------------------------------------------------------------

  write_dta(dta_sens, new_dta_path)
  message(sprintf("Written sensitivity .dta: %s", new_dta_path))

  # ---------------------------------------------------------------------------
  # Copy original LASSO globals CSV unchanged
  #
  # This keeps the same selected covariate sets by outcome/time horizon.
  # ---------------------------------------------------------------------------

  file.copy(
    from = old_globals_path,
    to = new_globals_path,
    overwrite = TRUE
  )

  message(sprintf("Copied globals CSV: %s", new_globals_path))
}