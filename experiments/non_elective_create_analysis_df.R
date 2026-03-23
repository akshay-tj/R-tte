# experiment_run_non_elective.R
#
# Builds the non-elective bypass study cohort from raw NVR + HES data,
# then calculates TTE outcomes (90 / 180 / 365 day) and runs LASSO IV
# variable selection.
#
# KNOWN DEVIATIONS FROM ORIGINAL (intentional):
#   1. Output format: original produced 3 separate per-horizon CSVs.
#      This script produces one wide tibble (all horizons combined).
#   2. Column names follow new package convention (see CLAUDE.md §7).

library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(rlang)

source("R/hes_utils.R")
source("R/nvr_utils.R")
source("R/utils.R")
source("R/outcomes.R")
source("R/descriptives.R")
source("R/lasso.R")
source("R/cohort.R")
source("R/Charlson_SCARF_scoring.R")

# =============================================================================
# PATHS
# =============================================================================

NVR_DF_PATH              <- "Z:/PHP/HSR/ESORT-V/ESORT-V/NVR Data - May 2025/NVR data for ESORT.xlsx"
HES_APC_PATH             <- "Z:/PHP/HSR/ESORT-V/ESORT-V/HES Data - May 2025/HES_data_concatenated_across_years/Cleaned_up_HES_data/HES_APC_2015_to_2023_variables_of_interest_only.qs"
HES_MORT_PATH            <- "Z:/PHP/HSR/ESORT-V/ESORT-V/HES Data - May 2025/HES_data_concatenated_across_years/HES_CIVREG_MORT.txt"
AVG_BYPASS_SURGERIES_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/NVR Data - May 2025/Bypass_subsets/Avg_bypass_surgeries_for_clti_per_hospital.csv"
NON_ELECTIVE_COHORT_BASELINE_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_bypass_study_participants_with_covariates_230326.csv"

NVR_WANTED_COLS <- c(
  "ProcedureType", "Patient:PatientId", "Patient:AgeAtSurgery",
  "Patient:GenderCode", "NvrEpisode:AdmissionDate",
  "NvrEpisode:ProcedureStartDate", "NvrHospitalName", "LSOA",
  "RiskScores:SmokingStatus", "Indications:AmpIndicationCode",
  "Indications:PadFontaineCode", "RiskScores:ASA",
  "RiskScores:Comorbidities", "NvrEpisode:AdmissionModeCode"
)

HES_APC_WANTED_COLS <- c(
  "STUDY_ID", "ADMIDATE", "EPISTART", "DISDATE", "EPIEND",
  "DIAG_4_CONCAT", "ADMISORC", "DISDEST", "IMD04_DECILE"
)

NVR_ID_COL        <- "Patient:PatientId"
NVR_ADMISSION_COL <- "NvrEpisode:AdmissionDate"
HES_APC_ID_COL    <- "STUDY_ID"
HES_ID_CLEAN_COL  <- "STUDY_ID_clean"

# =============================================================================
# SECTION 1: HES — load and clean once
# =============================================================================
# All downstream sections reuse HES_APC_df_clean

HES_APC_df_clean <- read_qs_df(HES_APC_PATH, HES_APC_WANTED_COLS) %>%
  clean_HES_admidate(id_col = HES_APC_ID_COL) %>%
  clean_HES_df_id_for_matching(id_col = HES_APC_ID_COL)

# =============================================================================
# SECTION 2: NVR cohort — inclusion, linkage, deduplication
# =============================================================================

HIGH_VOL_HOSPS <- read_csv(AVG_BYPASS_SURGERIES_PATH) %>%
  filter(Avg_2015_2023 >= 20) %>%
  pull(NvrHospitalName)

INCLUSION_CRITERIA <- list(
  procedure_type       = quo(ProcedureType %in% c("Lower-limb Bypass",
                                                   "Lower-limb Bypass (Linked)")),
  clti_only            = quo(`Indications:AmpIndicationCode` %in% c(2, 4)),
  fontaine_score       = quo(`Indications:PadFontaineCode` %in% c(3, 4)),
  asa_score            = quo(`RiskScores:ASA` %in% 1:4),
  english_lsoa         = quo(grepl("^E", LSOA)),
  high_volume_hospital = quo(NvrHospitalName %in% HIGH_VOL_HOSPS)
)

nvr_df <- read_xlsx_df(NVR_DF_PATH, NVR_WANTED_COLS)

nvr_included <- do.call(
  perform_inclusion,
  c(list(df    = nvr_df %>%
           filter(`NvrEpisode:AdmissionDate` >= "2015-01-01") %>%
           distinct(),
         label = "NVR_include"),
    INCLUSION_CRITERIA)
)

check_duplicates(nvr_included, NVR_ID_COL)

nvr_deduped <- handle_duplicates(nvr_included, id_col = NVR_ID_COL)

# HES linkage: step 1 (any HES record) → step 2 (matched admission date)
nvr_linked <- nvr_deduped %>%
  filter_NVR_patients_with_HES_record(
    HES_APC_df_clean, NVR_ID_COL, HES_ID_CLEAN_COL
  ) %>%
  filter_NVR_patients_with_matched_admission(
    HES_APC_df_clean, NVR_ID_COL, HES_ID_CLEAN_COL,
    NVR_ADMISSION_COL
  )

# =============================================================================
# SECTION 3: Non-elective cohort — arm assignment and derived variables
# =============================================================================

non_elective_cohort <- nvr_linked %>%
  filter(`NvrEpisode:AdmissionModeCode` == 2) %>%
  mutate(
    daystosurgery = as.numeric(
      as.Date(`NvrEpisode:ProcedureStartDate`) -
      as.Date(`NvrEpisode:AdmissionDate`)
    )
  ) %>%
  filter(daystosurgery > 0 & daystosurgery <= 21) %>%
  mutate(
    early_surgery = if_else(daystosurgery <= 5, 1L, 0L),
    year_of_surgery = year(as.Date(`NvrEpisode:ProcedureStartDate`)),
    covid_time_period = case_when(
      `NvrEpisode:AdmissionDate` <= as.Date("2020-03-15") ~ "pre_covid",
      `NvrEpisode:AdmissionDate` <= as.Date("2022-03-11") ~ "covid",
      TRUE                                                 ~ "post_covid"
    )
  ) %>%
  separate_out_nvr_comorbidities() %>%
  select(-comorbidity_NA, -daystosurgery)

print(sprintf("Cohort size after all selection criteria/linkage, and complete case filtering: %d", nrow(non_elective_cohort)))

# =============================================================================
# SECTION 4: HES — filter to cohort and compute index-admission covariates
# =============================================================================

# All HES rows for cohort patients
HES_cohort_all <- subset_HES_to_NVR_cohort(
  HES_APC_df_clean, non_elective_cohort, NVR_ID_COL, HES_ID_CLEAN_COL
)

# Join NVR admission date so we can identify the index admission row
HES_cohort_all <- HES_cohort_all %>%
  left_join(
    non_elective_cohort %>%
      select(`Patient:PatientId`, `NvrEpisode:AdmissionDate`) %>%
      mutate(STUDY_ID_clean = as.character(`Patient:PatientId`)) %>%
      rename(index_admidate = `NvrEpisode:AdmissionDate`) %>%
      select(STUDY_ID_clean, index_admidate),
    by = "STUDY_ID_clean"
  )

# Sanity check: every patient should have exactly one index admission date
HES_cohort_all %>%
  group_by(STUDY_ID_clean) %>%
  summarise(n_index = n_distinct(index_admidate)) %>%
  filter(n_index != 1)

HES_index <- HES_cohort_all %>%
  filter(ADMIDATE == index_admidate)

# ── Charlson / RCS ────────────────────────────────────────────────────────────
charlson_flags <- HES_index %>%
  add_flags_from_concat(charlson_all_list, icd10_groups_charlson, "DIAG_4_CONCAT") %>%
  group_by(STUDY_ID_clean) %>%
  summarise(
    across(
      all_of(paste0(charlson_all_list, "_yn")),
      ~ as.integer(any(. %in% c(1L, TRUE), na.rm = TRUE)),
      .names = "pat_level_{.col}"
    ),
    .groups = "drop"
  ) %>%
  rename_with(~ sub("^pat_level_(.*)_yn$", "pat_level_\\1", .x),
            starts_with("pat_level_")) %>%
  add_rcs_score(charlson_all_list) %>%
  select(STUDY_ID_clean, rcs_ch_cat)

# ── SCARF ─────────────────────────────────────────────────────────────────────
scarf_flags <- HES_index %>%
  add_flags_from_concat(scarf_deficits, icd10_groups_frailty, "DIAG_4_CONCAT") %>%
  add_scarf_admin_flags() %>%
  group_by(STUDY_ID_clean) %>%
  summarise(
    across(
      all_of(paste0(c(scarf_deficits, "admi_scarf", "discharge_scarf"), "_yn")),
      ~ as.integer(any(. %in% c(1L, TRUE), na.rm = TRUE)),
      .names = "pat_level_{.col}"
    ),
    .groups = "drop"
  ) %>%
  rename_with(~ sub("^pat_level_(.*)_yn$", "pat_level_\\1", .x),
              starts_with("pat_level_")) %>%
  add_scarf_score(c(scarf_deficits, "admi_scarf", "discharge_scarf")) %>%
  select(STUDY_ID_clean, scarf_cat)

# ── IMD ───────────────────────────────────────────────────────────────────────
imd_flags <- HES_index %>%
  select(STUDY_ID_clean, ADMIDATE, IMD04_DECILE) %>%
  distinct() %>%
  HES_imd_decile_to_quintile() %>%
  select(STUDY_ID_clean, IMD_quintile)

# =============================================================================
# SECTION 5: Assemble final cohort and write CSV
# =============================================================================

non_elective_cohort <- non_elective_cohort %>%
  rename(STUDY_ID = `Patient:PatientId`) %>%
  left_join(charlson_flags %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID") %>%
  left_join(scarf_flags   %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID") %>%
  left_join(imd_flags     %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID") %>%
  filter(complete.cases(.))

write.csv(non_elective_cohort, NON_ELECTIVE_COHORT_BASELINE_DF_PATH, row.names = FALSE)
message(sprintf("  Cohort written: %d patients → %s", nrow(non_elective_cohort), NON_ELECTIVE_COHORT_BASELINE_DF_PATH))

# =============================================================================
# SECTION 6: HES — filter to final cohort + DISDATE cleaning + prep for outcomes
# =============================================================================
# Now that non_elective_cohort is finalised (complete cases only), build the
# cohort-specific HES object used by calculate_outcomes().

sentinel_dates <- as.Date(c("1800-01-01", "1801-01-01"))

HES_APC_non_elective_cohort <- HES_APC_df_clean %>%
  semi_join(
    non_elective_cohort %>% mutate(STUDY_ID = as.character(STUDY_ID)),
    by = c("STUDY_ID_clean" = "STUDY_ID")
  ) %>%
  group_by(STUDY_ID_clean, ADMIDATE) %>%
  mutate(
    DISDATE = if (all(is.na(DISDATE) | DISDATE %in% sentinel_dates)) {
      max(EPIEND, na.rm = TRUE)
    } else {
      max(DISDATE[!is.na(DISDATE) & !DISDATE %in% sentinel_dates], na.rm = TRUE)
    }
  ) %>%
  ungroup() %>%
  left_join(
    non_elective_cohort %>%
      select(STUDY_ID, `NvrEpisode:AdmissionDate`) %>%
      rename(
        STUDY_ID_clean     = STUDY_ID,
        nvr_admission_date = `NvrEpisode:AdmissionDate`
      ) %>%
      mutate(STUDY_ID_clean = as.character(STUDY_ID_clean)),
    by = "STUDY_ID_clean"
  ) %>%
  distinct()

# =============================================================================
# SECTION 7: MORTALITY
# =============================================================================

mortality_clean <- read.table(HES_MORT_PATH, header = TRUE, sep = "|") %>%
  select(STUDY_ID, REG_DATE_OF_DEATH) %>%
  mutate(
    study_id   = str_remove(STUDY_ID, "^AB"),
    death_date = as.Date(as.character(REG_DATE_OF_DEATH), format = "%Y%m%d")
  ) %>%
  select(study_id, death_date)

# =============================================================================
# SECTION 8: TTE outcomes
# =============================================================================

non_elective_outcomes <- calculate_outcomes(
  cohort = non_elective_cohort %>%
    mutate(STUDY_ID = as.character(STUDY_ID)) %>%
    rename(
      study_id                 = STUDY_ID,
      nvr_admission_date       = `NvrEpisode:AdmissionDate`,
      nvr_procedure_start_date = `NvrEpisode:ProcedureStartDate`
    ),
  hes_admissions = HES_APC_non_elective_cohort %>%
    rename(
      study_id           = STUDY_ID_clean,
      hes_admission_date = ADMIDATE,
      hes_discharge_date = DISDATE
    ),
  mortality                       = mortality_clean,
  intervention_name               = "bypass_surg",
  intervention_admission_date_col = "nvr_admission_date",
  intervention_date_col           = "nvr_procedure_start_date",
  starting_point_col              = "nvr_procedure_start_date",
  time_horizons                   = c(90, 180, 365),
  include_pre_intervention        = FALSE,
  mortality_id_prefix             = ""
)

# =============================================================================
# SECTION 9: Descriptive outcomes tables (Table 2)
# =============================================================================

non_elective_outcomes <- non_elective_outcomes %>%
  left_join(
    non_elective_cohort %>%
      select(STUDY_ID, early_surgery) %>%
      mutate(study_id = as.character(STUDY_ID)),
    by = "study_id"
  )

cont_90 <- c("daoh_bypass_surg_90d", "total_los_no_90d",
             "bypass_surg_proc_los_no", "post_bypass_surg_los_no_90d",
             "bypass_surg_los_no")
cat_90  <- c("readmit_post_bypass_surg_90d", "died_post_bypass_surg_90d")
outcomes_90_days <- table2_outcomes(non_elective_outcomes,
                                    cont_vars = cont_90, cat_vars = cat_90,
                                    horizon_label = "90 days")

cont_180 <- c("daoh_bypass_surg_180d", "total_los_no_180d",
              "bypass_surg_proc_los_no", "post_bypass_surg_los_no_180d",
              "bypass_surg_los_no")
cat_180  <- c("readmit_post_bypass_surg_180d", "died_post_bypass_surg_180d")
outcomes_180_days <- table2_outcomes(non_elective_outcomes,
                                     cont_vars = cont_180, cat_vars = cat_180,
                                     horizon_label = "180 days")

cont_365 <- c("daoh_bypass_surg_365d", "total_los_no_365d",
              "bypass_surg_proc_los_no", "post_bypass_surg_los_no_365d",
              "bypass_surg_los_no")
cat_365  <- c("readmit_post_bypass_surg_365d", "died_post_bypass_surg_365d")
outcomes_365_days <- table2_outcomes(non_elective_outcomes,
                                     cont_vars = cont_365, cat_vars = cat_365,
                                     horizon_label = "365 days")
