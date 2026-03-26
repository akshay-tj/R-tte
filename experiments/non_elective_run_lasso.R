# non_elective_run_lasso.R
#
# Runs LASSO IV variable selection for all outcomes × timepoints for the
# non-elective bypass cohort. Exports one analysis-ready DTA + globals CSV
# per timepoint for Stata modelling.
#
# NOTE: instrumental_variable is currently loaded from an external CSV.
#       This will be computed within the pipeline in a future update.
#
# Inputs:
#   - non_elective_outcomes (wide tibble from non_elective_create_analysis_df.R)
#   - non_elective_cohort   (baseline covariates)
#   - temp IV file          (temporary — see NOTE above)
#
# Outputs (one per timepoint):
#   - non_elective_{H}d.dta          — analysis-ready dataset, union of all
#                                       LASSO-selected columns across outcomes,
#                                       plus pre-generated D*X and Z*X columns
#   - non_elective_{H}d_globals.csv  — outcome -> Xlist/Mlist mapping for Stata

library(dplyr)
library(haven)

source("R/lasso.R")

# =============================================================================
# PATHS
# =============================================================================

IV_PATH    <- "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/non_elective_clinical_effectiveness_df_270226_90dayonly.csv"
OUTPUT_DIR <- "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/march23_lasso_outputs/"

# =============================================================================
# PARAMETERS
# =============================================================================

TIME_HORIZONS <- c(90, 180, 365)

# Outcome column names per timepoint — {H} replaced at runtime
OUTCOMES <- list(
  daoh        = "daoh_bypass_surg_{H}d",
  total_los   = "total_los_no_{H}d",
  readmission = "readmit_post_bypass_surg_{H}d",
  mortality   = "died_post_bypass_surg_{H}d"
)

# glmnet family per outcome
OUTCOME_FAMILIES <- c(
  daoh        = "gaussian",
  total_los   = "gaussian",
  readmission = "binomial",
  mortality   = "binomial"
)

# Prespecified subgroups — always forced in, unpenalised, D*X and Z*X forced
not_penalized <- c("ageatsurgery", "gender_F", "fontaine_4", "comorbidity_1")

# Penalised main effects — age_sq included here per protocol (not forced)
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
  "medication_1", "medication_2", "medication_3", "medication_4",
  "surgyr_2016", "surgyr_2017", "surgyr_2018", "surgyr_2019",
  "surgyr_2020", "surgyr_2021", "surgyr_2022", "surgyr_2023",
  "covid_period", "postcovid_period"
)

# =============================================================================
# PREPARE BASE DATASET (shared across all timepoints)
# =============================================================================

# NOTE: instrumental_variable loaded externally — will move into pipeline later
iv_data <- read.csv(IV_PATH) %>%
  select(STUDY_ID, instrumental_variable) %>%
  mutate(STUDY_ID = as.integer(STUDY_ID))

base_df <- non_elective_outcomes %>%
  left_join(
    non_elective_cohort %>%
      mutate(STUDY_ID = as.integer(STUDY_ID)) %>%
      select(-early_surgery),
    by = c("study_id" = "STUDY_ID")
  ) %>%
  left_join(iv_data, by = c("study_id" = "STUDY_ID")) %>%
  rename(
    ageatsurgery    = `Patient.AgeAtSurgery`,
    nvrhospitalname = NvrHospitalName
  ) %>%
  mutate(
    age_sq           = ageatsurgery^2,  # NOTE: penalised per protocol, not forced
    gender_F         = if_else(`Patient.GenderCode` == 2, 1L, 0L),
    smoking_ex       = if_else(`RiskScores.SmokingStatus` == 2, 1L, 0L),
    smoking_current  = if_else(`RiskScores.SmokingStatus` == 3, 1L, 0L),
    asa_2            = if_else(`RiskScores.ASA` == 2, 1L, 0L),
    asa_3            = if_else(`RiskScores.ASA` == 3, 1L, 0L),
    asa_4            = if_else(`RiskScores.ASA` == 4, 1L, 0L),
    scarf_mild       = if_else(scarf_cat == "Mild Frailty (2/3)",     1L, 0L),
    scarf_moderate   = if_else(scarf_cat == "Moderate Frailty (4/5)", 1L, 0L),
    scarf_severe     = if_else(scarf_cat == "Severe Frailty (6+)",    1L, 0L),
    rcs_two          = if_else(rcs_ch_cat == "Two",      1L, 0L),
    rcs_threeplus    = if_else(rcs_ch_cat == "Three +",  1L, 0L),
    imd_q2           = if_else(IMD_quintile == "Q2",                  1L, 0L),
    imd_q3           = if_else(IMD_quintile == "Q3",                  1L, 0L),
    imd_q4           = if_else(IMD_quintile == "Q4",                  1L, 0L),
    imd_q5           = if_else(IMD_quintile == "Q5 (most deprived)",  1L, 0L),
    fontaine_4       = if_else(`Indications.PadFontaineCode` == 4,    1L, 0L),
    amp_tissueloss   = if_else(`Indications.AmpIndicationCode` == 4,  1L, 0L),
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

# =============================================================================
# RUN LASSO — all outcomes x timepoints
# =============================================================================

# Results stored as lasso_results[[horizon]][[outcome_name]]
lasso_results <- list()

for (h in TIME_HORIZONS) {

  lasso_results[[as.character(h)]] <- list()

  for (outcome_name in names(OUTCOMES)) {

    outcome_col <- gsub("\\{H\\}", h, OUTCOMES[[outcome_name]])
    message(sprintf("\n===== LASSO: %s | %dd =====", outcome_name, h))

    lasso_results[[as.character(h)]][[outcome_name]] <- run_lasso_iv_selection(
      dataset                = base_df,
      outcome                = outcome_col,
      treatment              = "early_surgery",
      instrument             = "instrumental_variable",
      prespecified_subgroups = not_penalized,
      penalized_main_effects = penalized_vars,
      family_stage1          = "binomial",
      family_stage2          = OUTCOME_FAMILIES[[outcome_name]],
      seed                   = 1276
    )
  }
}

# =============================================================================
# EXPORT — one DTA + globals CSV per timepoint
# =============================================================================

for (h in TIME_HORIZONS) {

  h_results <- lasso_results[[as.character(h)]]

  # ── Union of X*X terms across all outcomes at this horizon
  all_xx_outcome_terms  <- unique(unlist(lapply(h_results, `[[`, "retained_xx_outcome")))
  all_xx_exposure_terms <- unique(unlist(lapply(h_results, `[[`, "retained_xx_exposure")))
  all_xx_terms          <- unique(c(all_xx_outcome_terms, all_xx_exposure_terms))

  # ── Union of D*X and Z*X base variables across all outcomes
  all_dx_base <- unique(gsub("_d$",  "", unlist(lapply(h_results, `[[`, "retained_dx_outcome"))))
  all_zx_base <- unique(gsub("_iv$", "", unlist(lapply(h_results, `[[`, "retained_zx_exposure"))))

  # All base variables that need D*X or Z*X pre-generated
  # Always include prespecified subgroups (forced per protocol)
  dx_base_vars <- unique(c(not_penalized, all_dx_base))
  zx_base_vars <- unique(c(not_penalized, all_zx_base))

  # ── Build export dataset
  outcome_cols <- sapply(OUTCOMES, function(x) gsub("\\{H\\}", h, x))

  export_df <- base_df

  # Pre-generate X*X interaction columns
  for (term in all_xx_terms) {
    if (!term %in% names(export_df)) {
      parts <- strsplit(term, "_x_")[[1]]
      export_df[[term]] <- export_df[[parts[1]]] * export_df[[parts[2]]]
    }
  }

  # Pre-generate D*X interaction columns (avoids factor variable issues in Stata)
  # Named d_{base_var} to avoid collision with LASSO _d suffix convention
  dx_col_names <- paste0("d_x_", dx_base_vars)
  for (i in seq_along(dx_base_vars)) {
    export_df[[dx_col_names[i]]] <- export_df[["early_surgery"]] * export_df[[dx_base_vars[i]]]
  }

  # Pre-generate Z*X interaction columns for each version:
  #   z_x_stage2_{var} — Z x D*X base vars (for z_x_stage2_treatment version)
  #   z_x_stage1_{var} — Z x Z*X base vars (for z_x_stage1_instrument version)
  #   z_x_full_{var}   — Z x all main effects (for z_x_stage1_full version)
  all_main_base <- unique(c(not_penalized, penalized_vars))

  # z_x_stage2_: Z x vars selected for Stage 2 (Mlist / D*X base vars)
  #   used in z_x_stage2_treatment version of Stage 1
  # z_x_stage1_: Z x vars selected for Stage 1 (Zlist / Z*X base vars, LASSO Part 2b)
  #   used in z_x_stage1_instrument version of Stage 1
  # z_x_full_:   Z x all main effects — used in z_x_stage1_full version
  z_x_stage2_col_names <- paste0("z_x_stage2_", dx_base_vars)
  z_x_stage1_col_names <- paste0("z_x_stage1_", zx_base_vars)
  zx_full_col_names    <- paste0("z_x_full_",   all_main_base)

  for (i in seq_along(dx_base_vars)) {
    export_df[[z_x_stage2_col_names[i]]] <- export_df[["instrumental_variable"]] * export_df[[dx_base_vars[i]]]
  }
  for (i in seq_along(zx_base_vars)) {
    export_df[[z_x_stage1_col_names[i]]] <- export_df[["instrumental_variable"]] * export_df[[zx_base_vars[i]]]
  }
  for (i in seq_along(all_main_base)) {
    if (all_main_base[i] %in% names(export_df)) {
      export_df[[zx_full_col_names[i]]] <- export_df[["instrumental_variable"]] * export_df[[all_main_base[i]]]
    }
  }

  export_df <- export_df %>%
    select(
      study_id, early_surgery, instrumental_variable, nvrhospitalname,
      all_of(not_penalized), all_of(penalized_vars),
      all_of(unname(outcome_cols)),
      any_of(all_xx_terms),
      any_of(dx_col_names),
      any_of(z_x_stage2_col_names),
      any_of(z_x_stage1_col_names),
      any_of(zx_full_col_names)
    )

  # ── Write DTA
  dta_path <- file.path(OUTPUT_DIR, sprintf("non_elective_%dd.dta", h))
  haven::write_dta(export_df, dta_path)
  message(sprintf("Written: %s", dta_path))

  # ── Build globals CSV
  globals_rows <- lapply(names(h_results), function(outcome_name) {
    res         <- h_results[[outcome_name]]
    outcome_col <- gsub("\\{H\\}", h, OUTCOMES[[outcome_name]])

    # D*X base vars for this outcome (prespec + LASSO-selected)
    dx_vars_this_outcome <- unique(c(
      not_penalized,
      gsub("_d$", "", res$retained_dx_outcome)
    ))
    dx_cols_this_outcome <- paste0("d_x_", dx_vars_this_outcome)

    # Z*X base vars for this outcome (prespec + LASSO-selected, Part 2b)
    zx_vars_this_outcome <- unique(c(
      not_penalized,
      gsub("_iv$", "", res$retained_zx_exposure)
    ))

    # Stage 2 Xlist: union main effects + X*X outcome + D*X cols
    xlist_stage2 <- unique(c(
      res$union_main_effects,
      res$retained_xx_outcome,
      dx_cols_this_outcome
    ))

    # Stage 1 Xlists: controls only (union main effects + X*X)
    # Z*X terms moved to zlists — they are excluded instruments, not controls
    xlist_stage1_z_only <- unique(c(
      res$union_main_effects,
      res$retained_xx_outcome
    ))
    zlist_stage1_z_only <- character(0)

    xlist_stage1_z_x_stage2_treatment <- unique(c(
      res$union_main_effects,
      res$retained_xx_outcome
    ))
    zlist_stage1_z_x_stage2_treatment <- paste0("z_x_stage2_", dx_vars_this_outcome)

    xlist_stage1_z_x_stage1_full <- unique(c(
      res$union_main_effects,
      res$retained_xx_outcome
    ))
    zlist_stage1_z_x_stage1_full <- paste0("z_x_full_", res$union_main_effects)

    # v4: controls use X*X exposure (not outcome) per LASSO Part 2b selection
    xlist_stage1_z_x_stage1_instrument <- unique(c(
      res$union_main_effects,
      res$retained_xx_exposure
    ))
    zlist_stage1_z_x_stage1_instrument <- paste0("z_x_stage1_", zx_vars_this_outcome)

    data.frame(
      outcome      = outcome_col,
      family       = OUTCOME_FAMILIES[[outcome_name]],
      xlist_s2     = paste(xlist_stage2,                            collapse = " "),
      xlist_s1_v1  = paste(xlist_stage1_z_only,                    collapse = " "),
      xlist_s1_v2  = paste(xlist_stage1_z_x_stage2_treatment,      collapse = " "),
      xlist_s1_v3  = paste(xlist_stage1_z_x_stage1_full,           collapse = " "),
      xlist_s1_v4  = paste(xlist_stage1_z_x_stage1_instrument,     collapse = " "),
      zlist_s1_v1  = paste(zlist_stage1_z_only,                    collapse = " "),
      zlist_s1_v2  = paste(zlist_stage1_z_x_stage2_treatment,      collapse = " "),
      zlist_s1_v3  = paste(zlist_stage1_z_x_stage1_full,           collapse = " "),
      zlist_s1_v4  = paste(zlist_stage1_z_x_stage1_instrument,     collapse = " "),
      stringsAsFactors = FALSE
    )
  })

  globals_df   <- do.call(rbind, globals_rows)
  globals_path <- file.path(OUTPUT_DIR, sprintf("non_elective_%dd_globals.csv", h))
  write.csv(globals_df, globals_path, row.names = FALSE)
  message(sprintf("Written: %s", globals_path))
}