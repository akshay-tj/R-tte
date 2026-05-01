library(survival)
library(tidyverse)
library(haven)

source("R/lasso.R")

ELECTIVE_COHORT_BASELINE_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_with_confounders.csv"
ELECTIVE_COHORT_OUTCOMES_DF_PATH <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_with_outcomes.csv"
ELECTIVE_IV_OUTPUT_PATH   <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/analysable_subsets/elective_bypass_study_participants_with_ttes.csv"

elective_cohort <- read.csv(ELECTIVE_COHORT_BASELINE_DF_PATH)
elective_outcomes <- read.csv(ELECTIVE_COHORT_OUTCOMES_DF_PATH)
elective_iv <- read.csv(ELECTIVE_IV_OUTPUT_PATH)

# ── 1. Build analysis dataset ─────────────────────────────────────────────────

afs_df <- elective_cohort %>%
  mutate(study_id = as.character(STUDY_ID)) %>%
  left_join(
    elective_outcomes %>%
      mutate(study_id = as.character(study_id)) %>%
      select(study_id, afs_days, afs_event),
    by = "study_id"
  ) %>%
  left_join(
    elective_iv %>%
      mutate(study_id = as.character(STUDY_ID)) %>%
      select(study_id, instrumental_variable),
    by = "study_id"
  ) %>%
  filter(!is.na(afs_days), afs_days >= 0) %>%
  mutate(
    ageatsurgery     = Patient.AgeAtSurgery,
    age_sq           = Patient.AgeAtSurgery^2,
    gender_F         = if_else(Patient.GenderCode == 2,                         1L, 0L),
    smoking_ex       = if_else(RiskScores.SmokingStatus == 2,                   1L, 0L),
    smoking_current  = if_else(RiskScores.SmokingStatus == 3,                   1L, 0L),
    asa_2            = if_else(RiskScores.ASA == 2,                             1L, 0L),
    asa_3            = if_else(RiskScores.ASA == 3,                             1L, 0L),
    asa_4            = if_else(RiskScores.ASA == 4,                             1L, 0L),
    scarf_mild       = if_else(scarf_cat == "Mild Frailty (2/3)",               1L, 0L),
    scarf_moderate   = if_else(scarf_cat == "Moderate Frailty (4/5)",           1L, 0L),
    scarf_severe     = if_else(scarf_cat == "Severe Frailty (6+)",              1L, 0L),
    rcs_two          = if_else(rcs_ch_cat == "Two",                             1L, 0L),
    rcs_threeplus    = if_else(rcs_ch_cat == "Three +",                         1L, 0L),
    imd_q2           = if_else(IMD_quintile == "Q2",                            1L, 0L),
    imd_q3           = if_else(IMD_quintile == "Q3",                            1L, 0L),
    imd_q4           = if_else(IMD_quintile == "Q4",                            1L, 0L),
    imd_q5           = if_else(IMD_quintile == "Q5 (most deprived)",            1L, 0L),
    fontaine_4       = if_else(Indications.PadFontaineCode == 4,                1L, 0L),
    amp_tissueloss   = if_else(Indications.AmpIndicationCode == 4,              1L, 0L),
    surgyr_2016      = if_else(year_of_surgery == 2016, 1L, 0L),
    surgyr_2017      = if_else(year_of_surgery == 2017, 1L, 0L),
    surgyr_2018      = if_else(year_of_surgery == 2018, 1L, 0L),
    surgyr_2019      = if_else(year_of_surgery == 2019, 1L, 0L),
    surgyr_2020      = if_else(year_of_surgery == 2020, 1L, 0L),
    surgyr_2021      = if_else(year_of_surgery == 2021, 1L, 0L),
    surgyr_2022      = if_else(year_of_surgery == 2022, 1L, 0L),
    surgyr_2023      = if_else(year_of_surgery == 2023, 1L, 0L),
    covid_period     = if_else(covid_time_period == "covid",      1L, 0L),
    postcovid_period = if_else(covid_time_period == "post_covid", 1L, 0L),
    medication_1     = as.integer(medication_1),
    medication_2     = as.integer(medication_2),
    medication_3     = as.integer(medication_3),
    medication_4     = as.integer(medication_4),
    krt_yn           = as.integer(krt_yn)
  )

# ── 2. Covariate vectors ──────────────────────────────────────────────────────
# surgyr_2015 is reference category — not included

not_penalized <- c("ageatsurgery", "gender_F", "fontaine_4", "comorbidity_1", "krt_yn")

penalized_vars <- c(
  "age_sq",
  "smoking_ex", "smoking_current",
  "asa_2", "asa_3", "asa_4",
  "scarf_mild", "scarf_moderate", "scarf_severe",
  "rcs_two", "rcs_threeplus",
  "imd_q2", "imd_q3", "imd_q4", "imd_q5",
  "amp_tissueloss",
  "comorbidity_2", "comorbidity_3", "comorbidity_4", "comorbidity_5",
  "comorbidity_6", "comorbidity_7", "comorbidity_8",
  "medication_1", "medication_2", "medication_3", "medication_4", "surgyr_2016",
  "surgyr_2017", "surgyr_2018", "surgyr_2019",
  "surgyr_2020", "surgyr_2021", "surgyr_2022", "surgyr_2023",
  "covid_period", "postcovid_period"
)

# ── 3. LASSO selection ────────────────────────────────────────────────────────
# family_stage2 = "gaussian" is an approximation used for covariate selection
# only — it determines which terms enter the Cox model, not the estimator
# itself. The actual treatment effect is estimated by the Cox frailty model
# in Step 8.

lasso_afs <- run_lasso_iv_selection(
  dataset                = afs_df,
  outcome                = "afs_days",
  treatment              = "early_surgery",
  instrument             = "instrumental_variable",
  prespecified_subgroups = not_penalized,
  penalized_main_effects = penalized_vars,
  family_stage1          = "binomial",
  family_stage2          = "gaussian",
  seed                   = 1276
)

# ── 4. Pre-generate interaction columns in afs_df ─────────────────────────────

# D*X interactions — retained (penalised pool) + forced (prespecified subgroups)
# Both sets appear in final_outcome_terms but only retained ones are in
# retained_dx_outcome; prespec D*X are forced in unpenalised and must be
# derived separately
all_dx_terms <- grep("_d$", lasso_afs$final_outcome_terms, value = TRUE)

for (term in all_dx_terms) {
  base_var <- sub("_d$", "", term)
  if (!term %in% names(afs_df)) {
    if (!base_var %in% names(afs_df)) {
      warning(sprintf("Base variable '%s' not found in afs_df — skipping '%s'", base_var, term))
      next
    }
    afs_df[[term]] <- afs_df[["early_surgery"]] * afs_df[[base_var]]
  }
}

# X*X interactions
all_xx_terms <- grep("_x_", lasso_afs$final_outcome_terms, value = TRUE)

for (term in all_xx_terms) {
  if (!term %in% names(afs_df)) {
    parts <- strsplit(term, "_x_")[[1]]
    if (length(parts) != 2 || !all(parts %in% names(afs_df))) {
      warning(sprintf("Could not create interaction column: %s", term))
      next
    }
    afs_df[[term]] <- afs_df[[parts[1]]] * afs_df[[parts[2]]]
  }
}

# ── 5. Attach generalised residuals ───────────────────────────────────────────

afs_df <- afs_df %>%
  mutate(
    generalised_residual     = lasso_afs$generalised_residuals,
    generalised_residual_x_d = generalised_residual * early_surgery
  )

# ── 6. Build Cox covariate formula terms ──────────────────────────────────────
# final_outcome_terms from LASSO already contains:
#   union main effects + retained D*X + retained X*X + residual cols
# Strip residual cols — added explicitly below with meaningful names

cox_terms <- setdiff(
  lasso_afs$final_outcome_terms,
  c(".ctrl_resid", ".ctrl_resid_x_d")
)

cox_formula_str <- paste0(
  "Surv(afs_days, afs_event) ~ ",
  paste(cox_terms, collapse = " + "),
  " + generalised_residual + generalised_residual_x_d"
)

# ── 7. Proportional hazards check (no frailty — cox.zph does not support it) ──

cox_ph_check <- coxph(as.formula(cox_formula_str), data = afs_df)

ph_test <- cox.zph(cox_ph_check)
print(ph_test)
plot(ph_test)

# ── 8. Cox 2SRI model with frailty ───────────────────────────────────────────
# Patient-level Gaussian frailty captures unmeasured individual heterogeneity.
# set.seed() is required here because the frailty EM algorithm uses random
# initialisation; fixing it ensures the model is reproducible.

set.seed(1)

cox_afs <- coxph(
  as.formula(paste0(cox_formula_str, " + frailty(study_id, dist = 'gauss')")),
  data = afs_df
)

summary(cox_afs)

# ── 9. Extract HR for early_surgery ──────────────────────────────────────────

c_idx <- which(names(coef(cox_afs)) == "early_surgery")
beta  <- coef(cox_afs)[c_idx]
se    <- sqrt(vcov(cox_afs)[c_idx, c_idx])
HR    <- exp(beta)
CI    <- exp(beta + c(-1, 1) * 1.96 * se)

cat(sprintf("HR: %.3f  95%% CI: %.3f - %.3f\n", HR, CI[1], CI[2]))