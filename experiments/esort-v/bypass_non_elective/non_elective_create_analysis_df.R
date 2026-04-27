# experiment_run_non_elective.R
#
# Builds the non-elective bypass study cohort from raw NVR + HES data,
# then calculates TTE outcomes (90 / 180 / 365 day) and writes two CSVs: 
# 1) baseline characteristics + confounders, 2) outcomes.

library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(rlang)

source("R/hes_utils.R")
source("R/nvr_utils.R")
source("R/utils.R")
source("R/outcomes.R")
source("R/lasso.R")
source("R/cohort.R")
source("R/Charlson_SCARF_scoring.R")
source("R/vascular_outcomes.R")

# =============================================================================
# PATHS
# =============================================================================

# Input DF Paths: 
NVR_DF_PATH              <- "Z:/PHP/HSR/ESORT-V/ESORT-V/NVR Data - May 2025/NVR data for ESORT.xlsx"
HES_APC_PATH             <- "Z:/PHP/HSR/ESORT-V/ESORT-V/HES Data - May 2025/HES_data_concatenated_across_years/HES_APC_2015_to_2023.qs"
HES_MORT_PATH            <- "Z:/PHP/HSR/ESORT-V/ESORT-V/HES Data - May 2025/HES_data_concatenated_across_years/HES_CIVREG_MORT.txt"
AVG_BYPASS_SURGERIES_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/NVR Data - May 2025/Bypass_subsets/Avg_bypass_surgeries_for_clti_per_hospital.csv"

# Input codelist paths: 
KRT_ICD_1_PATH  <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/krt_codelists/hes_codelist_additionalkrt.dta"
KRT_OPCS_1_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/krt_codelists/hes_codelist_additionalkrt_opcs.dta"
KRT_ICD_2_PATH  <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/krt_codelists/hes_codelist_dialysis.dta"
KRT_OPCS_2_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/krt_codelists/hes_codelist_dialysis_opcs.dta"
KRT_ICD_3_PATH  <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/krt_codelists/hes_codelist_kidneytransplant.dta"
KRT_OPCS_3_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/krt_codelists/hes_codelist_kidneytransplant_opcs.dta"

# Output paths: 
NON_ELECTIVE_COHORT_BASELINE_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/analysable_subsets/non_elective_bypass_study_participants_with_confounders.csv"
NON_ELECTIVE_COHORT_OUTCOMES_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/analysable_subsets/non_elective_bypass_study_participants_with_outcomes.csv"

# =============================================================================
# Parameters to define 
# =============================================================================

NVR_WANTED_COLS <- c(
  "ProcedureType", "Patient:PatientId", "Patient:AgeAtSurgery",
  "Patient:GenderCode", "NvrEpisode:AdmissionDate",
  "NvrEpisode:ProcedureStartDate", "NvrHospitalName", "LSOA",
  "RiskScores:SmokingStatus", "Indications:AmpIndicationCode",
  "Indications:PadFontaineCode", "RiskScores:ASA", "RiskScores:Medication",
  "RiskScores:Comorbidities", "NvrEpisode:AdmissionModeCode", 
  "PostOp:CompFurtherSurgeryCode", "Indications:IndicationSideCode"
)

HES_APC_WANTED_COLS <- c(
  "STUDY_ID", "ADMIDATE", "EPISTART", "DISDATE", "EPIEND",
  "DIAG_4_CONCAT", "ADMISORC", "DISDEST", "IMD04_DECILE", 
  "OPERTN_4_CONCAT",
  paste0("OPERTN_", sprintf("%02d", 1:24)),
  paste0("OPDATE_", sprintf("%02d", 1:24))
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
  { 
    message("Admission mode counts (before filter):")
    print(count(., `NvrEpisode:AdmissionModeCode`, sort = TRUE))
    . 
  } %>%
  filter(`NvrEpisode:AdmissionModeCode` == 2) %>%
  { message("After admission mode filter — rows: ", nrow(.),
            " | unique patients: ", n_distinct(.[["Patient:PatientId"]])); . } %>%
  mutate(
    daystosurgery = as.numeric(
      as.Date(`NvrEpisode:ProcedureStartDate`) -
      as.Date(`NvrEpisode:AdmissionDate`)
    )
  ) %>%
  {
    message("daystosurgery breakdown (before filter):")
    message("  daystosurgery <= 0  : ", sum(.$daystosurgery <= 0,  na.rm = TRUE), " rows")
    message("  daystosurgery 1-21  : ", sum(.$daystosurgery >= 1 & .$daystosurgery <= 21, na.rm = TRUE), " rows")
    message("  daystosurgery > 21  : ", sum(.$daystosurgery > 21,  na.rm = TRUE), " rows")
    message("  daystosurgery NA    : ", sum(is.na(.$daystosurgery)), " rows")
    .
  } %>%
  filter(daystosurgery > 0 & daystosurgery <= 21) %>%
  { message("After daystosurgery filter — rows: ", nrow(.),
            " | unique patients: ", n_distinct(.[["Patient:PatientId"]])); . } %>%
  mutate(
    early_surgery = if_else(daystosurgery <= 5, 1L, 0L),
    year_of_surgery = year(as.Date(`NvrEpisode:ProcedureStartDate`)),
    covid_time_period = case_when(
      `NvrEpisode:AdmissionDate` <= as.Date("2020-03-15") ~ "pre_covid",
      `NvrEpisode:AdmissionDate` <= as.Date("2022-03-11") ~ "covid",
      TRUE                                                 ~ "post_covid"
    )
  ) %>%
  separate_out_nvr_pipe_col(
    col    = "RiskScores:Comorbidities",
    prefix = "comorbidity_"
  ) %>%
  separate_out_nvr_pipe_col(
    col    = "RiskScores:Medication",
    prefix = "medication_"
  ) %>%
  mutate(
    # Codes 8 and 9 are subtypes of medication group 1
    medication_1 = if_else(medication_8 == 1L | medication_9 == 1L, 1L, medication_1)
  ) %>%
  select(-medication_0, -medication_5, -medication_6,
         -medication_7, -medication_8, -medication_9)

names(non_elective_cohort)
non_elective_cohort <- non_elective_cohort %>% select(-c(`comorbidity_8 `, `comorbidity_7 `))

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

# ── KRT codelist ──────────────────────────────────────────────────────────────
krt_icd_codes <- c(
  gsub("\\.", "", haven::read_dta(KRT_ICD_1_PATH)$icd10),
  gsub("\\.", "", haven::read_dta(KRT_ICD_2_PATH)$icd10),
  gsub("\\.", "", haven::read_dta(KRT_ICD_3_PATH)$icd10)
)

krt_opcs_codes <- c(
  haven::read_dta(KRT_OPCS_1_PATH)$opcs,
  haven::read_dta(KRT_OPCS_2_PATH)$opcs,
  haven::read_dta(KRT_OPCS_3_PATH)$opcs
)

krt_flags <- HES_index %>%
  add_flags_from_concat("krt_icd",  list(krt_icd  = krt_icd_codes),  "DIAG_4_CONCAT") %>%
  add_flags_from_concat("krt_opcs", list(krt_opcs = krt_opcs_codes), "OPERTN_4_CONCAT") %>%
  group_by(STUDY_ID_clean) %>%
  summarise(
    krt_yn = as.integer(any(krt_icd_yn == 1L | krt_opcs_yn == 1L, na.rm = TRUE)),
    .groups = "drop"
  )

# ── IMD ───────────────────────────────────────────────────────────────────────
imd_flags <- HES_index %>%
  select(STUDY_ID_clean, ADMIDATE, IMD04_DECILE) %>%
  distinct() %>%
  HES_imd_decile_to_quintile() %>%
  select(STUDY_ID_clean, IMD_quintile)

# =============================================================================
# SECTION 5: Assemble final cohort and write CSV
# =============================================================================

# Assemble the final cohort-level DF by left-joining the HES-derived covariate flags to the NVR cohort.  
non_elective_cohort <- non_elective_cohort %>%
  rename(STUDY_ID = `Patient:PatientId`) %>%
  left_join(charlson_flags %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID") %>%
  left_join(scarf_flags   %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID") %>%
  left_join(imd_flags     %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID") %>%
  left_join(krt_flags     %>% rename(STUDY_ID = STUDY_ID_clean), by = "STUDY_ID")

# check for missing data in the final cohort DF
cat("Rows:", nrow(non_elective_cohort), "\n\nMissing data summary:\n")
print(colSums(is.na(non_elective_cohort)))

# keep only complete cases (ignoring missingness in "PostOp:CompFurtherSurgeryCode", "Indications:IndicationSideCode" since these are not confounders)
non_elective_cohort <- non_elective_cohort %>%  filter(complete.cases(select(., -`PostOp:CompFurtherSurgeryCode`, -`Indications:IndicationSideCode`)))

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
# SECTION 8: TTE outcomes - Primary 
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
# SECTION 9: TTE outcomes - Secondary  
# =============================================================================

# set OPSC revas code prefixes based on https://doi.org/10.1016/j.ejvs.2023.05.003 
revasc_prefixes   <- c("L161", "L162", "L163", "L168", "L169",  
"L206", "L216",  "L501", "L502", "L503", "L504", "L505", "L506",  
"L511", "L512", "L513", "L514", "L515", "L516", "L518", "L519",  
"L521", "L522", "L528", "L529",  "L581", "L582", "L583", "L584", "L585", "L586", "L587", "L588", "L589", 
"L591", "L592", "L593", "L594", "L595", "L596", "L597", "L598", "L599", 
"L601", "L602", "L603", "L604", "L608", "L609",  
"L651", "L652", "L653", "L658", "L659", 
"L681", "L682", "L683", "L684", "L688", "L689",
"L541", "L544", "L548", "L549",  
"L631", "L635", "L638", "L639", 
"L662", "L665", "L667", "L668", "L669", 
"L711", "L718", "L719")

amp_prefixes <- "X09"

censoring_date    <- as.Date("2024-03-31")  # study end date

# --- Pre-filter HES to post-index episodes only ----------------------------
HES_post_index_df_for_secondary_outcomes <- subset_HES_to_NVR_cohort(HES_cohort_all, non_elective_cohort, "STUDY_ID", HES_ID_CLEAN_COL) %>% 
                          select(-c("EPISTART", "EPIEND", "DISDATE", "DIAG_4_CONCAT", "ADMISORC", "DISDEST", "IMD04_DECILE", "OPERTN_4_CONCAT")) %>% 
                          filter(ADMIDATE >  index_admidate) %>%  # keep only index and post-index admissions 
                          rename(study_id = STUDY_ID_clean) # rename for consistency with limb event flagging function

# --- Extract index limb laterality from NVR --------------------------------
nvr_index_laterality <- non_elective_cohort %>%
                        select(study_id = STUDY_ID, index_side = `Indications:IndicationSideCode`) %>%
                        distinct()

# --- Flag earliest same-side limb events in HES ----------------------------
limb_events <- flag_limb_events(
  hes_post_index  = HES_post_index_df_for_secondary_outcomes,
  amp_prefixes    = amp_prefixes,
  revasc_prefixes = revasc_prefixes,
  index_side_df   = nvr_index_laterality,
  use_laterality  = TRUE
  # opertn_cols / opdate_cols: package defaults (OPERTN_01..24) apply
)

# --- Compute ILR, ILMA, AFS outcomes ---------------------------------------
limb_outcomes <- compute_limb_outcomes(
  cohort         = non_elective_cohort %>% rename(study_id = STUDY_ID),
  limb_events    = limb_events,
  mortality      = mortality_clean,
  start_date_col = "NvrEpisode:ProcedureStartDate",
  censoring_date = censoring_date,
  time_horizons  = c(90L, 180L, 365L)
)

non_elective_outcomes <- non_elective_outcomes %>%
  left_join(limb_outcomes, by = "study_id") %>%
  left_join(
    non_elective_cohort %>%
      select(STUDY_ID, early_surgery) %>%
      mutate(study_id = as.character(STUDY_ID)),
    by = "study_id"
  )

write.csv(non_elective_outcomes,
          NON_ELECTIVE_COHORT_OUTCOMES_DF_PATH,
          row.names = FALSE)

message(sprintf("Outcomes written: %d patients → %s", nrow(non_elective_outcomes), NON_ELECTIVE_COHORT_OUTCOMES_DF_PATH))
