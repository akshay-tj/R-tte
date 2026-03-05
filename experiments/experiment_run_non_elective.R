#library(ttemulate)
library(dplyr)
library(readr)
library(lubridate)
library(stringr)

# KNOWN DEVIATIONS FROM ORIGINAL (intentional):
#   1. Mortality flag: original used abs(died_date - index_admidate) <= H, which
#      incorrectly flags deaths BEFORE time zero if within H days.
#      calculate_outcomes() uses directional logic: death_date >= starting_point.
#      Expect a small number of patients to differ on died_within_{H}_days.
#
#   2. Output format: original produced 3 separate per-horizon CSVs.
#      This script produces one wide tibble (all horizons combined).
#
#   3. Column names follow new package convention (see CLAUDE.md §7).

# =============================================================================
# PATHS
# =============================================================================

HES_APC_PATH             <- "Z:/PHP/HSR/ESORT-V/ESORT-V/HES Data - May 2025/HES_data_concatenated_across_years/Cleaned_up_HES_data/HES_APC_2015_to_2023_variables_of_interest_only.qs"
HES_MORT_PATH            <- "Z:/PHP/HSR/ESORT-V/ESORT-V/HES Data - May 2025/HES_data_concatenated_across_years/HES_CIVREG_MORT.txt"
NON_ELECTIVE_COHORT_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_bypass_study_participants_with_covariates.csv"
HES_APC_WANTED_COLS      <- c("STUDY_ID", "ADMIDATE", "EPISTART", "DISDATE", "EPIEND")
HES_APC_ID_COL           <- "STUDY_ID"

# =============================================================================
# NVR COHORT — loaded first; required for HES semi_join below
# =============================================================================
#
# GOTCHA: this CSV is a pre-pipeline output — standardise_nvr_names() has NOT
# been applied. Column names are raw NVR format, e.g.:
#   Patient:PatientId       (not patient_patientid)
#   NvrEpisode.AdmissionDate (not nvr_admission_date)
# Check column names with names(non_elective_cohort) before running joins.

non_elective_cohort <- read.csv(NON_ELECTIVE_COHORT_PATH)
names(non_elective_cohort)
# =============================================================================
# HES APC
# =============================================================================

HES_APC_df <- read_qs_df(HES_APC_PATH, HES_APC_WANTED_COLS) %>%
  clean_HES_admidate(id_col = HES_APC_ID_COL) %>%
  clean_HES_df_id_for_matching(id_col = HES_APC_ID_COL)

# GOTCHA: RHS of join uses raw Patient:PatientId — confirm exact column name
# in non_elective_cohort with names() before running.
HES_APC_non_elective_cohort <- HES_APC_df %>%
  semi_join(
    non_elective_cohort %>% mutate(STUDY_ID = as.character(STUDY_ID)),
    by = c("STUDY_ID_clean" = "STUDY_ID")
  )

# DISDATE cleaning
sentinel_dates <- as.Date(c("1800-01-01", "1801-01-01"))

HES_APC_non_elective_cohort <- HES_APC_non_elective_cohort %>%
  group_by(STUDY_ID_clean, ADMIDATE) %>%
  mutate(
    DISDATE = if (all(is.na(DISDATE) | DISDATE %in% sentinel_dates)) {
      max(EPIEND, na.rm = TRUE)
    } else {
      max(DISDATE[!is.na(DISDATE) & !DISDATE %in% sentinel_dates], na.rm = TRUE)
    }
  ) %>%
  ungroup()

# GOTCHA: NvrEpisode.AdmissionDate is the raw column name — readr/readxl may
# have mangled the colon to a dot on import. Confirm with names() before joining.
HES_APC_non_elective_cohort <- HES_APC_non_elective_cohort %>%
  left_join(
    non_elective_cohort %>%
      select(STUDY_ID, `NvrEpisode.AdmissionDate`) %>%
      rename(
        STUDY_ID_clean     = STUDY_ID,
        nvr_admission_date = `NvrEpisode.AdmissionDate`
      ) %>%
      mutate(
        STUDY_ID_clean     = as.character(STUDY_ID_clean),
        nvr_admission_date = as.Date(nvr_admission_date)
      ),
    by = "STUDY_ID_clean"
  ) %>%
  distinct()

# =============================================================================
# MORTALITY
# =============================================================================

mortality_clean <- read.table(HES_MORT_PATH, header = TRUE, sep = "|") %>%
  select(STUDY_ID, REG_DATE_OF_DEATH) %>%
  mutate(
    study_id   = str_remove(STUDY_ID, "^AB"),
    death_date = as.Date(as.character(REG_DATE_OF_DEATH), format = "%Y%m%d")
  ) %>%
  select(study_id, death_date)

# =============================================================================
# TTE-1: Non-Elective Bypass Outcomes
# =============================================================================
non_elective_outcomes <- calculate_outcomes(
  cohort = non_elective_cohort %>%
    mutate(STUDY_ID = as.character(STUDY_ID)) %>%
    rename(
      study_id           = STUDY_ID,
      nvr_admission_date = NvrEpisode.AdmissionDate
    ),
  hes_admissions = HES_APC_non_elective_cohort %>%
    rename(
      study_id           = STUDY_ID_clean,
      hes_admission_date = ADMIDATE,
      hes_discharge_date = DISDATE
    ),
  mortality                = mortality_clean,
  intervention_name        = "bypass_surg",
  intervention_date_col    = "nvr_admission_date",
  starting_point_col       = "nvr_admission_date",
  time_horizons            = c(90L),
  include_pre_intervention = FALSE,
  mortality_id_prefix      = ""
)

cat("n =", nrow(non_elective_outcomes), "\n")
cat("cols =", ncol(non_elective_outcomes), "\n")
glimpse(non_elective_outcomes)

explore_df <- non_elective_outcomes %>%
  select(study_id, 
        daoh_bypass_surg_90d,
        total_los_no_90d,
        post_bypass_surg_los_no_90d, 
        died_post_bypass_surg_90d) %>% left_join(
          non_elective_cohort %>%
            select(STUDY_ID, early_surgery) %>%
            mutate(study_id = as.character(STUDY_ID)),
          by = "study_id"
        )

explore_df %>% group_by(early_surgery) %>%
  summarise(
    n = n(),
    mean_daoh_90d = mean(daoh_bypass_surg_90d),
    mean_los_90d = mean(total_los_no_90d),
    mortality_rate_90d = mean(died_post_bypass_surg_90d)
  )

#write_csv(
#  non_elective_outcomes,
#  "analysable_subsets/non_elective_bypass_outcomes_all_horizons.csv"
#)