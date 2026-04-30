source("R/descriptives.R")

# =============================================================================
# Descriptive tables — Elective LL bypass cohort
#
# Assumes the following objects are already in the environment,
# produced by elective_create_analysis_df.R:
#   - elective_cohort   : cohort tibble with baseline characteristics
#   - elective_outcomes : tibble with all outcome columns
#
# Produces:
#   - Table 1 : baseline characteristics stratified by early_surgery
#   - Table 2 : unadjusted outcomes at 90, 180, and 365 days
#   - afs_crude_rates  : AFS incidence rates per 1000 person-years by arm
#   - afs_rate_diff    : AFS rate difference (early minus late) with 95% CI
# =============================================================================

# =============================================================================
# Table 1 — baseline characteristics
# =============================================================================

elective_cohort_table1 <- elective_cohort %>%
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
  data       = elective_cohort_table1,
  strata_col = "early_surgery",
  cont_vars  = cont_vars_t1,
  cat_vars   = cat_vars_t1,
  ttest      = FALSE,
  label      = "baseline characteristics — LL bypass elective"
)

# =============================================================================
# Descriptive outcomes tables (Table 2)
# =============================================================================

# get stats for incomplete follow-up due to study end date (31st March 2024)
# 180 day outcomes; participants with procedure after 30th September 2023 will not have full 180 day follow-up
elective_cohort %>%
  mutate(after_date = `NvrEpisode:ProcedureStartDate` >= as.Date("2023-09-30")) %>%
  group_by(early_surgery) %>%
  summarise(
    n_total = n(),
    n_after = sum(after_date, na.rm = TRUE),
    pct = round(n_after / n_total * 100, 2)
  )

# 365 day outcomes; participants with procedure after 30th March 2023 will not have full 365 day follow-up 
elective_cohort %>%
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
outcomes_90_days <- descriptive_table(elective_outcomes, strata_col = "early_surgery",
                                    cont_vars = cont_90, cat_vars = cat_90, ttest = TRUE,
                                    label = "90 days")

cont_180 <- c("daoh_bypass_surg_180d", "total_los_no_180d", "post_bypass_surg_los_no_180d")
cat_180  <- c("readmit_post_bypass_surg_180d", "died_post_bypass_surg_180d", "ilr_180d", "ilma_180d")
outcomes_180_days <- descriptive_table(elective_outcomes, strata_col = "early_surgery",
                                      cont_vars = cont_180, cat_vars = cat_180, ttest = TRUE,
                                      label = "180 days")

cont_365 <- c("daoh_bypass_surg_365d", "total_los_no_365d", "post_bypass_surg_los_no_365d")
cat_365  <- c("readmit_post_bypass_surg_365d", "died_post_bypass_surg_365d", "ilr_365d", "ilma_365d")
outcomes_365_days <- descriptive_table(elective_outcomes, strata_col = "early_surgery",
                                     cont_vars = cont_365, cat_vars = cat_365, ttest = TRUE,
                                     label = "365 days")

## AFS crude rates and rate differences at 90, 180, and 365 days

afs_horizons <- c(90, 180, 365)

afs_crude_rates <- purrr::map_dfr(afs_horizons, function(h) {
  elective_outcomes %>%
    mutate(
      afs_days_h  = pmin(afs_days, h),
      afs_event_h = as.integer(afs_event == 1L & afs_days <= h)
    ) %>%
    group_by(early_surgery) %>%
    summarise(
      horizon           = h,
      n_total           = n(),
      n_events          = sum(afs_event_h, na.rm = TRUE),
      person_years      = sum(afs_days_h,  na.rm = TRUE) / 365.25,
      crude_rate_1000py = round(n_events / person_years * 1000, 2),
      .groups           = "drop"
    )
}) %>%
  dplyr::relocate(horizon, early_surgery)

afs_rate_diff <- afs_crude_rates %>%
  dplyr::select(horizon, early_surgery, n_events, person_years) %>%
  tidyr::pivot_wider(
    names_from  = early_surgery,
    values_from = c(n_events, person_years),
    names_sep   = "_"
  ) %>%
  mutate(
    rate_0 = n_events_0 / person_years_0 * 1000,
    rate_1 = n_events_1 / person_years_1 * 1000,
    rd     = rate_1 - rate_0,
    var_0  = n_events_0 / person_years_0^2 * 1000^2,
    var_1  = n_events_1 / person_years_1^2 * 1000^2,
    se_rd  = sqrt(var_0 + var_1),
    ci_lo  = rd - 1.96 * se_rd,
    ci_hi  = rd + 1.96 * se_rd
  ) %>%
  dplyr::select(horizon, rate_0, rate_1, rd, ci_lo, ci_hi) %>%
  mutate(across(c(rate_0, rate_1, rd, ci_lo, ci_hi), ~ round(.x, 2)))

print(afs_rate_diff)
