#library(ttemulate)
library(dplyr)
library(readr)
library(lubridate)
library(stringr)

# KNOWN DEVIATIONS FROM ORIGINAL (intentional):
#
#   1. Output format: original produced 3 separate per-horizon CSVs.
#      This script produces one wide tibble (all horizons combined) -- needs modification
#
#   2. Column names follow new package convention (see CLAUDE.md §7) -- needs modification

# temp - source (we can then replace with package installation once finalised)
source("R/hes_utils.R")
source("R/nvr_utils.R")
source("R/utils.R")
source("R/outcomes.R")
source("R/descriptives.R")
source("R/lasso.R")

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
      nvr_admission_date = NvrEpisode.AdmissionDate, 
      nvr_procedure_start_date = NvrEpisode.ProcedureStartDate
    ),
  hes_admissions = HES_APC_non_elective_cohort %>%
    rename(
      study_id           = STUDY_ID_clean,
      hes_admission_date = ADMIDATE,
      hes_discharge_date = DISDATE
    ),
  mortality                = mortality_clean,
  intervention_name        = "bypass_surg",
  intervention_admission_date_col = "nvr_admission_date",
  intervention_date_col    = "nvr_procedure_start_date", # nvr_admission_date
  starting_point_col       = "nvr_procedure_start_date", # nvr_admission_date --- gives the OG results
  time_horizons            = c(90, 180, 365),
  include_pre_intervention = FALSE,
  mortality_id_prefix      = ""
)

non_elective_outcomes_90d <- non_elective_outcomes %>%
  select(study_id, 
         daoh_bypass_surg_90d,
         total_los_no_90d,
         readmit_post_bypass_surg_90d,
         died_post_bypass_surg_90d) %>% 
  left_join(non_elective_cohort %>% select(-c(X)) %>% mutate(STUDY_ID = as.character(STUDY_ID)), 
  by = c("study_id" = "STUDY_ID"))

View(non_elective_outcomes_90d)

cat("n =", nrow(non_elective_outcomes), "\n")
cat("cols =", ncol(non_elective_outcomes), "\n")
# glimpse(non_elective_outcomes)

non_elective_outcomes <- non_elective_outcomes %>%
   left_join(
          non_elective_cohort %>%
            select(STUDY_ID, early_surgery) %>%
            mutate(study_id = as.character(STUDY_ID)),
          by = "study_id"
        )

# 90 days
cont_90 <- c(
  "daoh_bypass_surg_90d",
  "total_los_no_90d",
  "bypass_surg_proc_los_no",
  "post_bypass_surg_los_no_90d",
  "bypass_surg_los_no"
)
cat_90 <- c(
  "readmit_post_bypass_surg_90d",
  "died_post_bypass_surg_90d"
)

outcomes_90_days <- table2_outcomes(non_elective_outcomes, cont_vars = cont_90, cat_vars = cat_90,
                    horizon_label = "90 days")
outcomes_90_days

# 180 days
cont_180 <- c(
  "daoh_bypass_surg_180d",
  "total_los_no_180d",
  "bypass_surg_proc_los_no",
  "post_bypass_surg_los_no_180d",
  "bypass_surg_los_no"
)
cat_180 <- c(
  "readmit_post_bypass_surg_180d",
  "died_post_bypass_surg_180d"
)

outcomes_180_days <- table2_outcomes(non_elective_outcomes, cont_vars = cont_180, cat_vars = cat_180,
                    horizon_label = "180 days")
# 365 days
cont_365 <- c(
  "daoh_bypass_surg_365d",
  "total_los_no_365d",
  "bypass_surg_proc_los_no",
  "post_bypass_surg_los_no_365d",
  "bypass_surg_los_no"
)
cat_365 <- c(
  "readmit_post_bypass_surg_365d",
  "died_post_bypass_surg_365d"
)

outcomes_365_days <- table2_outcomes(non_elective_outcomes, cont_vars = cont_365, cat_vars = cat_365,
                    horizon_label = "365 days")

# LASSO 

non_elective_outcomes_90d <- non_elective_outcomes_90d %>%
  mutate(
    # Sex
    gender_F = if_else(Patient.GenderCode == 2, 1, 0),
    
    # Smoking
    smoking_ex      = if_else(RiskScores.SmokingStatus == 2, 1, 0),
    smoking_current = if_else(RiskScores.SmokingStatus == 3, 1, 0),
    
    # ASA
    asa_2 = if_else(RiskScores.ASA == 2, 1, 0),
    asa_3 = if_else(RiskScores.ASA == 3, 1, 0),
    asa_4 = if_else(RiskScores.ASA == 4, 1, 0),
    
    # SCARF
    scarf_mild     = if_else(scarf_cat == "Mild",     1, 0),
    scarf_moderate = if_else(scarf_cat == "Moderate", 1, 0),
    scarf_severe   = if_else(scarf_cat == "Severe",   1, 0),
    
    # RCS
    rcs_two       = if_else(rcs_ch_cat == 2, 1, 0),
    rcs_threeplus = if_else(rcs_ch_cat == 3, 1, 0),
    
    # IMD
    imd_q2 = if_else(IMD_quintile == 2, 1, 0),
    imd_q3 = if_else(IMD_quintile == 3, 1, 0),
    imd_q4 = if_else(IMD_quintile == 4, 1, 0),
    imd_q5 = if_else(IMD_quintile == 5, 1, 0),
    
    # Fontaine
    fontaine_4 = if_else(Indications.PadFontaineCode == 4, 1, 0),
    
    # Amputation indication
    amp_tissueloss = if_else(Indications.AmpIndicationCode == 4, 1, 0),

    # Surgery year
    surgyr_2016 = if_else(year_of_surgery == 2016, 1, 0),
    surgyr_2017 = if_else(year_of_surgery == 2017, 1, 0),
    surgyr_2018 = if_else(year_of_surgery == 2018, 1, 0),
    surgyr_2019 = if_else(year_of_surgery == 2019, 1, 0),
    surgyr_2020 = if_else(year_of_surgery == 2020, 1, 0),
    surgyr_2021 = if_else(year_of_surgery == 2021, 1, 0),
    surgyr_2022 = if_else(year_of_surgery == 2022, 1, 0),
    surgyr_2023 = if_else(year_of_surgery == 2023, 1, 0),
    
    # Covid period
    covid_period     = if_else(covid_time_period == "covid",      1, 0),
    postcovid_period = if_else(covid_time_period == "post_covid", 1, 0)
  )

  temp_90day <- read.csv("Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_clinical_effectiveness_df_270226_90dayonly.csv")
  temp_90day <- temp_90day %>% select(STUDY_ID, medication_2, medication_1, medication_3, medication_4, instrumental_variable)

  non_elective_outcomes_90d <- non_elective_outcomes_90d %>%
  left_join(
    temp_90day %>% mutate(STUDY_ID = as.character(STUDY_ID)),
    by = c("study_id" = "STUDY_ID")
  )
# Variable lists
################################################################################

not_penalized <- c("Patient.AgeAtSurgery", "gender_F", "fontaine_4", "comorbidity_1")

penalized_vars <- c(
  "smoking_ex", "smoking_current",
  "asa_2", "asa_3", "asa_4",
  "scarf_mild", "scarf_moderate", "scarf_severe",
  "rcs_two", "rcs_threeplus",
  "imd_q2", "imd_q3", "imd_q4", "imd_q5",
  "amp_tissueloss",
  "comorbidity_2", "comorbidity_3",
  "comorbidity_4", "comorbidity_5", "comorbidity_6",
  "comorbidity_7", "comorbidity_8",
  "medication_1", "medication_2", "medication_3", "medication_4",
  "surgyr_2016", "surgyr_2017", "surgyr_2018", "surgyr_2019",
  "surgyr_2020", "surgyr_2021", "surgyr_2022", "surgyr_2023",
  "covid_period", "postcovid_period"
)

# 3. Run
result <- run_lasso_iv_selection(
  dataset                = non_elective_outcomes_90d,
  outcome                = "daoh_bypass_surg_90d",
  treatment              = "early_surgery",
  instrument             = "instrumental_variable",
  prespecified_subgroups = not_penalized,
  penalized_main_effects = penalized_vars,
  family_stage1          = "binomial",
  family_stage2          = "gaussian",
  seed                   = 1276
)

write.csv(non_elective_outcomes_90d, "non_elective_outcomes_90d_temp_06March26.csv")
