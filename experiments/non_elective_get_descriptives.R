source("R/descriptives.R")

# =============================================================================
# SECTION X: Table 1 — baseline characteristics
# =============================================================================
# NOTE: comorbidity_7 and comorbidity_8 have duplicate columns with trailing
# spaces ("comorbidity_7 ", "comorbidity_8 ") — data quality issue to fix in
# non_elective_create_analysis_df.R at ingestion. Using clean versions here.

non_elective_cohort_table1 <- non_elective_cohort %>%
  mutate(
    `Patient:GenderCode` = factor(
      `Patient:GenderCode`,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    `RiskScores:SmokingStatus` = factor(
      `RiskScores:SmokingStatus`,
      levels = c(1, 2, 3),
      labels = c("Current", "Ex-smoker", "Never smoked")
    ),
    `RiskScores:ASA` = factor(
      `RiskScores:ASA`,
      levels = c(1, 2, 3, 4),
      labels = c("Normal", "Mild", "Severe (not life-threatening)", "Severe (life-threatening)")
    ),
    `Indications:PadFontaineCode` = factor(
      `Indications:PadFontaineCode`,
      levels = c(3, 4),
      labels = c("Resting pain", "Necrosis/gangrene")
    ),
    `Indications:AmpIndicationCode` = factor(
      `Indications:AmpIndicationCode`,
      levels = c("2", "4"),
      labels = c("Chronic limb ischaemia", "Tissue loss")
    ),
    year_of_surgery = factor(
      year_of_surgery,
      levels = 2015:2024
    )
  )

cont_vars_t1 <- c(
  "Patient:AgeAtSurgery"
)

cat_vars_t1 <- c(
  "Patient:GenderCode",
  "RiskScores:SmokingStatus",
  "comorbidity_1",
  "comorbidity_2",
  "comorbidity_3",
  "comorbidity_4",
  "comorbidity_5",
  "comorbidity_6",
  "comorbidity_7",
  "comorbidity_8",
  "comorbidity_0",
  "RiskScores:ASA",
  "Indications:PadFontaineCode",
  "Indications:AmpIndicationCode",
  "year_of_surgery",
  "krt_yn", 
  "rcs_ch_cat",      # already a factor — levels: None, One, Two, Three+
  "scarf_cat"
)

descriptive_table(
  data       = non_elective_cohort_table1,
  strata_col = "early_surgery",
  cont_vars  = cont_vars_t1,
  cat_vars   = cat_vars_t1,
  ttest      = FALSE,
  label      = "baseline characteristics — LL bypass non-elective"
)

# =============================================================================
# SECTION X: Descriptive outcomes tables (Table 2)
# =============================================================================

cont_90 <- c("daoh_bypass_surg_90d", "total_los_no_90d",
             "bypass_surg_proc_los_no", "post_bypass_surg_los_no_90d",
             "bypass_surg_los_no")
cat_90  <- c("readmit_post_bypass_surg_90d", "died_post_bypass_surg_90d")
outcomes_90_days <- descriptive_table(non_elective_outcomes, strata_col = "early_surgery",
                                    cont_vars = cont_90, cat_vars = cat_90, ttest = TRUE,
                                    label = "90 days")

cont_180 <- c("daoh_bypass_surg_180d", "total_los_no_180d",
              "bypass_surg_proc_los_no", "post_bypass_surg_los_no_180d",
              "bypass_surg_los_no")
cat_180  <- c("readmit_post_bypass_surg_180d", "died_post_bypass_surg_180d")
outcomes_180_days <- descriptive_table(non_elective_outcomes, strata_col = "early_surgery",
                                      cont_vars = cont_180, cat_vars = cat_180, ttest = TRUE,
                                      label = "180 days")

cont_365 <- c("daoh_bypass_surg_365d", "total_los_no_365d",
              "bypass_surg_proc_los_no", "post_bypass_surg_los_no_365d",
              "bypass_surg_los_no")
cat_365  <- c("readmit_post_bypass_surg_365d", "died_post_bypass_surg_365d")
outcomes_365_days <- descriptive_table(non_elective_outcomes, strata_col = "early_surgery",
                                     cont_vars = cont_365, cat_vars = cat_365, ttest = TRUE,
                                     label = "365 days")
