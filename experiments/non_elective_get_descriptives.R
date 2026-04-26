source("R/descriptives.R")

# =============================================================================
# Table 1 — baseline characteristics
# =============================================================================

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
# Descriptive outcomes tables (Table 2)
# =============================================================================

# get stats for incomplete follow-up due to study end date (31st March 2024)
# 180 day outcomes; participants with procedure after 30th September 2023 will not have full 180 day follow-up
non_elective_cohort %>%
  mutate(after_date = `NvrEpisode:ProcedureStartDate` >= as.Date("2023-09-30")) %>%
  group_by(early_surgery) %>%
  summarise(
    n_total = n(),
    n_after = sum(after_date, na.rm = TRUE),
    pct = round(n_after / n_total * 100, 2)
  )

# 365 day outcomes; participants with procedure after 30th March 2023 will not have full 365 day follow-up 
non_elective_cohort %>%
  mutate(after_date = `NvrEpisode:ProcedureStartDate` >= as.Date("2023-03-30")) %>%
  group_by(early_surgery) %>%
  summarise(
    n_total = n(),
    n_after = sum(after_date, na.rm = TRUE),
    pct = round(n_after / n_total * 100, 2)
  )

cont_90 <- c("daoh_bypass_surg_90d", "total_los_no_90d",
             "bypass_surg_proc_los_no", "post_bypass_surg_los_no_90d",
             "bypass_surg_los_no")
cat_90  <- c("readmit_post_bypass_surg_90d", "died_post_bypass_surg_90d", "ilr_90d", "ilma_90d")
outcomes_90_days <- descriptive_table(non_elective_outcomes, strata_col = "early_surgery",
                                    cont_vars = cont_90, cat_vars = cat_90, ttest = TRUE,
                                    label = "90 days")

cont_180 <- c("daoh_bypass_surg_180d", "total_los_no_180d", "post_bypass_surg_los_no_180d")
cat_180  <- c("readmit_post_bypass_surg_180d", "died_post_bypass_surg_180d", "ilr_180d", "ilma_180d")
outcomes_180_days <- descriptive_table(non_elective_outcomes, strata_col = "early_surgery",
                                      cont_vars = cont_180, cat_vars = cat_180, ttest = TRUE,
                                      label = "180 days")

cont_365 <- c("daoh_bypass_surg_365d", "total_los_no_365d", "post_bypass_surg_los_no_365d")
cat_365  <- c("readmit_post_bypass_surg_365d", "died_post_bypass_surg_365d", "ilr_365d", "ilma_365d")
outcomes_365_days <- descriptive_table(non_elective_outcomes, strata_col = "early_surgery",
                                     cont_vars = cont_365, cat_vars = cat_365, ttest = TRUE,
                                     label = "365 days")
