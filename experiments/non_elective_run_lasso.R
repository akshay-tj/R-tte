# =============================================================================
# SECTION 10: LASSO — covariate dummy coding + variable selection
# =============================================================================

non_elective_outcomes_90d <- non_elective_outcomes %>%
  select(study_id,
         daoh_bypass_surg_90d,
         total_los_no_90d,
         readmit_post_bypass_surg_90d,
         died_post_bypass_surg_90d) %>%
  left_join(
    non_elective_cohort %>%
      mutate(STUDY_ID = as.character(STUDY_ID)),
    by = c("study_id" = "STUDY_ID")
  ) %>%
  mutate(
    gender_F         = if_else(`Patient:GenderCode` == 2, 1L, 0L),
    smoking_ex       = if_else(`RiskScores:SmokingStatus` == 2, 1L, 0L),
    smoking_current  = if_else(`RiskScores:SmokingStatus` == 3, 1L, 0L),
    asa_2            = if_else(`RiskScores:ASA` == 2, 1L, 0L),
    asa_3            = if_else(`RiskScores:ASA` == 3, 1L, 0L),
    asa_4            = if_else(`RiskScores:ASA` == 4, 1L, 0L),
    scarf_mild       = if_else(scarf_cat == "Mild Frailty (2/3)",     1L, 0L),
    scarf_moderate   = if_else(scarf_cat == "Moderate Frailty (4/5)", 1L, 0L),
    scarf_severe     = if_else(scarf_cat == "Severe Frailty (6+)",    1L, 0L),
    rcs_two          = if_else(rcs_ch_cat == "Two",      1L, 0L),
    rcs_threeplus    = if_else(rcs_ch_cat == "Three +",  1L, 0L),
    imd_q2           = if_else(IMD_quintile == "Q2",               1L, 0L),
    imd_q3           = if_else(IMD_quintile == "Q3",               1L, 0L),
    imd_q4           = if_else(IMD_quintile == "Q4",               1L, 0L),
    imd_q5           = if_else(IMD_quintile == "Q5 (most deprived)", 1L, 0L),
    fontaine_4       = if_else(`Indications:PadFontaineCode` == 4, 1L, 0L),
    amp_tissueloss   = if_else(`Indications:AmpIndicationCode` == 4, 1L, 0L),
    surgyr_2016      = if_else(year_of_surgery == 2016, 1L, 0L),
    surgyr_2017      = if_else(year_of_surgery == 2017, 1L, 0L),
    surgyr_2018      = if_else(year_of_surgery == 2018, 1L, 0L),
    surgyr_2019      = if_else(year_of_surgery == 2019, 1L, 0L),
    surgyr_2020      = if_else(year_of_surgery == 2020, 1L, 0L),
    surgyr_2021      = if_else(year_of_surgery == 2021, 1L, 0L),
    surgyr_2022      = if_else(year_of_surgery == 2022, 1L, 0L),
    surgyr_2023      = if_else(year_of_surgery == 2023, 1L, 0L),
    covid_period     = if_else(covid_time_period == "covid",      1L, 0L),
    postcovid_period = if_else(covid_time_period == "post_covid", 1L, 0L)
  )

temp_90day <- read.csv("Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_clinical_effectiveness_df_270226_90dayonly.csv") %>%
  select(STUDY_ID, instrumental_variable)

non_elective_outcomes_90d <- non_elective_outcomes_90d %>%
  left_join(
    temp_90day %>% mutate(STUDY_ID = as.character(STUDY_ID)),
    by = c("study_id" = "STUDY_ID")
  )

not_penalized <- c("Patient:AgeAtSurgery", "gender_F", "fontaine_4", "comorbidity_1")

penalized_vars <- c(
  "smoking_ex", "smoking_current",
  "asa_2", "asa_3", "asa_4",
  "scarf_mild", "scarf_moderate", "scarf_severe",
  "rcs_two", "rcs_threeplus",
  "imd_q2", "imd_q3", "imd_q4", "imd_q5",
  "amp_tissueloss",
  "comorbidity_2", "comorbidity_3", "comorbidity_4", "comorbidity_5",
  "comorbidity_6", "comorbidity_7", "comorbidity_8",
  "medication_1", "medication_2", "medication_3", "medication_4",
  "surgyr_2016", "surgyr_2017", "surgyr_2018", "surgyr_2019",
  "surgyr_2020", "surgyr_2021", "surgyr_2022", "surgyr_2023",
  "covid_period", "postcovid_period"
)

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
