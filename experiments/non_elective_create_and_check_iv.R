#### new iv pipeline df 
library(tidyverse)

source("R/create_iv.R")
source("R/utils.R")
source("R/cohort.R")
# =============================================================================
# Parameters to define 

NVR_WANTED_COLS <- c(
  "ProcedureType", "Patient:PatientId", "Patient:AgeAtSurgery",
  "Patient:GenderCode", "NvrEpisode:AdmissionDate",
  "NvrEpisode:ProcedureStartDate", "NvrHospitalName", "LSOA",
  "RiskScores:SmokingStatus", "Indications:AmpIndicationCode",
  "Indications:PadFontaineCode", "RiskScores:ASA", "RiskScores:Medication",
  "RiskScores:Comorbidities", "NvrEpisode:AdmissionModeCode"
)

NVR_ID_COL        <- "Patient:PatientId"
NVR_ADMISSION_COL <- "NvrEpisode:AdmissionDate"
AVG_BYPASS_SURGERIES_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/NVR Data - May 2025/Bypass_subsets/Avg_bypass_surgeries_for_clti_per_hospital.csv"
NVR_DF_PATH              <- "Z:/PHP/HSR/ESORT-V/ESORT-V/NVR Data - May 2025/NVR data for ESORT.xlsx"
NON_ELECTIVE_COHORT_BASELINE_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_bypass_study_participants_with_covariates_080426.csv"

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

nvr_lookback_only <- do.call(
  perform_inclusion,
  c(list(df    = nvr_df %>%
           filter(`NvrEpisode:AdmissionDate` < "2015-01-01") %>%
           distinct(),
         label = "NVR_include"),
    INCLUSION_CRITERIA)
)

check_duplicates(nvr_lookback_only, NVR_ID_COL)

nvr_deduped <- handle_duplicates(nvr_lookback_only, id_col = NVR_ID_COL)

# =============================================================================
# Non-elective cohort — arm assignment and derived variables
# =============================================================================

nvr_lookback_only_deduped <- nvr_deduped %>%
  filter(`NvrEpisode:AdmissionModeCode` == 2) %>%
  mutate(
    daystosurgery = as.numeric(
      as.Date(`NvrEpisode:ProcedureStartDate`) -
      as.Date(`NvrEpisode:AdmissionDate`)
    )
  ) %>%
  filter(daystosurgery > 0 & daystosurgery <= 21) %>%
  mutate(
    early_surgery = if_else(daystosurgery <= 5, 1L, 0L))

nvr_lookback_only_deduped <- nvr_lookback_only_deduped %>% 
  rename(NvrEpisode.AdmissionDate = `NvrEpisode:AdmissionDate`, 
         Patient.PatientId = `Patient:PatientId`) %>% 
  select(early_surgery, Patient.PatientId, NvrEpisode.AdmissionDate, NvrHospitalName) %>% 
  mutate(include_flag = 0L)

non_elective_cohort_baseline_df <- read_csv(NON_ELECTIVE_COHORT_BASELINE_DF_PATH) %>% 
                                   select(STUDY_ID, early_surgery, `NvrEpisode:AdmissionDate`, NvrHospitalName) %>% 
                                      mutate(include_flag = 1L) %>% 
                                      rename(NvrEpisode.AdmissionDate = `NvrEpisode:AdmissionDate`, 
                                             Patient.PatientId = STUDY_ID)

combined_df <- bind_rows(nvr_lookback_only_deduped, non_elective_cohort_baseline_df)

combined_df <- combined_df %>%
  rowwise() %>%
  mutate(
    instrumental_variable = if_else(
      include_flag == 1,
      calc_ttes(
        combined_df,
        Patient.PatientId,
        NvrHospitalName,
        NvrEpisode.AdmissionDate,
        lookback_days = 365L
      ),
      NA_real_
    )
  ) %>%
  ungroup()

# for missing values fill with median of non missing values
median_iv <- median(combined_df$instrumental_variable, na.rm = TRUE)

combined_df <- combined_df %>%
  mutate(instrumental_variable = if_else(is.na(instrumental_variable), median_iv, instrumental_variable))

# create a df with study id and instrumental variable 
iv_df <- combined_df %>% 
  filter(include_flag == 1) %>%
  select(Patient.PatientId, instrumental_variable) %>% 
  # STUDY_ID should be integer format to match with the main analysis df
  rename(STUDY_ID = Patient.PatientId) %>%
  mutate(STUDY_ID = as.character(STUDY_ID))

write_csv(iv_df, "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_IV_080426.csv")
