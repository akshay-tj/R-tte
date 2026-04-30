# elective_run_lasso.R
#
# Runs LASSO IV variable selection for all outcomes × timepoints for the
# elective bypass cohort. Exports one analysis-ready DTA + globals CSV
# per timepoint for Stata modelling.
#
# Assumes the following objects are already in the environment:
#   - elective_outcomes (wide tibble from elective_create_analysis_df.R)
#   - elective_cohort   (baseline confounders)
#   - iv_df
#
# Outputs (one per timepoint):
#   - elective_{H}d.dta          — analysis-ready dataset, union of all
#                                       LASSO-selected columns across outcomes,
#                                       plus pre-generated D*X and Z*X columns
#   - elective_{H}d_globals.csv  — outcome -> Xlist/Mlist mapping for Stata

library(dplyr)
library(haven)
library(purrr)

source("R/lasso.R")

# =============================================================================
# PATHS
# =============================================================================

OUTPUT_DIR <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/lasso_outputs/"

# =============================================================================
# PARAMETERS
# =============================================================================

TIME_HORIZONS <- c(90, 180, 365)

# Outcome column names per timepoint — {H} replaced at runtime
OUTCOMES <- list(
  daoh        = "daoh_bypass_surg_{H}d",
  daoh_myles              = "daoh_myles_bypass_surg_{H}d",
  total_los   = "total_los_no_{H}d",
  readmission = "readmit_post_bypass_surg_{H}d",
  mortality   = "died_post_bypass_surg_{H}d",
  post_bypass_surg_los_no = "post_bypass_surg_los_no_{H}d",
  ilr                     = "ilr_{H}d",
  ilma                    = "ilma_{H}d"
)

# glmnet family per outcome
OUTCOME_FAMILIES <- c(
  daoh                    = "gaussian",
  daoh_myles              = "gaussian",
  total_los               = "gaussian",
  readmission             = "binomial",
  mortality               = "binomial",
  post_bypass_surg_los_no = "gaussian",
  ilr                     = "binomial",
  ilma                    = "binomial"
)

# Two 2SRI models are run (m1/m2 — see globals CSV building block below):
#   model1: homogeneous effects — Z only in Stage 1, no Z*X or D*X terms
#   model2: heterogeneous effects — Z + Z*X in Stage 1, D*X in Stage 2
# Prespecified subgroups — always forced in, unpenalised, D*X and Z*X forced
not_penalized <- c("ageatsurgery", "gender_F", "fontaine_4", "comorbidity_1", "krt_yn")

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
  "surgyr_2017", "surgyr_2018", "surgyr_2019",
  "surgyr_2020", "surgyr_2021", "surgyr_2022", "surgyr_2023",
  "covid_period", "postcovid_period"
)

# =============================================================================
# PREPARE BASE DATASET (shared across all timepoints)
# =============================================================================

# NOTE: instrumental_variable loaded externally — will move into pipeline later
base_df <- elective_outcomes %>%
  left_join(
    elective_cohort %>%
        select(-early_surgery),
    by = c("study_id" = "STUDY_ID")
  ) %>%
  left_join(iv_df, by = c("study_id" = "STUDY_ID")) %>%
  rename(
    ageatsurgery    = `Patient:AgeAtSurgery`,
    nvrhospitalname = NvrHospitalName
  ) %>%
  mutate(
    age_sq           = ageatsurgery^2,  # NOTE: penalised per protocol, not forced
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
    imd_q2           = if_else(IMD_quintile == "Q2",                  1L, 0L),
    imd_q3           = if_else(IMD_quintile == "Q3",                  1L, 0L),
    imd_q4           = if_else(IMD_quintile == "Q4",                  1L, 0L),
    imd_q5           = if_else(IMD_quintile == "Q5 (most deprived)",  1L, 0L),
    fontaine_4       = if_else(`Indications:PadFontaineCode` == 4,    1L, 0L),
    amp_tissueloss   = if_else(`Indications:AmpIndicationCode` == 4,  1L, 0L),
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
  
  # adding tryCatch as lasso fails to converge for binary outcome mortality 
  lasso_results[[as.character(h)]][[outcome_name]] <- tryCatch(
  run_lasso_iv_selection(
    dataset                = base_df,
    outcome                = outcome_col,
    treatment              = "early_surgery",
    instrument             = "instrumental_variable",
    prespecified_subgroups = not_penalized,
    penalized_main_effects = penalized_vars,
    family_stage1          = "binomial",
    family_stage2          = OUTCOME_FAMILIES[[outcome_name]],
    seed                   = 1276
  ),
  error = function(e) {
    message(sprintf("  SKIPPED: %s | %dd — %s", outcome_name, h, conditionMessage(e)))
    NULL
  }
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

  # ── Union of D*X and Z*X base variables across all outcomes
  all_dx_base <- unique(gsub("_d$",  "", unlist(lapply(h_results, `[[`, "retained_dx_outcome"))))
  all_zx_base <- unique(gsub("_iv$", "", unlist(lapply(h_results, `[[`, "retained_zx_exposure"))))

  # All base variables that need D*X or Z*X pre-generated
  # Always include prespecified subgroups (forced per protocol)
  dx_base_vars <- unique(c(not_penalized, all_dx_base))
  zx_base_vars <- unique(c(not_penalized, all_zx_base, all_dx_base))
  # include all_dx_base: per protocol, any X where D*X was retained in the
  # outcome model is also forced unpenalised as Z*X in the Stage 1 exposure
  # model (model2). The LASSO applies this correctly; this line ensures the
  # corresponding z_x_stage1_ columns are also pre-generated in the DTA for Stata.

  # ── Build export dataset
  outcome_cols <- sapply(OUTCOMES, function(x) gsub("\\{H\\}", h, x))

  export_df <- base_df

  # Pre-generate X*X interaction columns
  for (term in unique(c(all_xx_exposure_terms, all_xx_outcome_terms))) {
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
  
  # Z x Z*X base vars for Stage 1 (model2 only)
  z_x_stage1_col_names <- paste0("z_x_stage1_", zx_base_vars)

  for (i in seq_along(zx_base_vars)) {
    export_df[[z_x_stage1_col_names[i]]] <- export_df[["instrumental_variable"]] * export_df[[zx_base_vars[i]]]
  }

  export_df <- export_df %>%
    select(
      study_id, early_surgery, instrumental_variable, nvrhospitalname,
      all_of(not_penalized), all_of(penalized_vars),
      all_of(unname(outcome_cols)),
      any_of(all_xx_exposure_terms),
      any_of(all_xx_outcome_terms),
      any_of(dx_col_names),
      any_of(z_x_stage1_col_names)
    )

  # ── Write DTA
  dta_path <- file.path(OUTPUT_DIR, sprintf("elective_%dd.dta", h))
  haven::write_dta(export_df, dta_path)
  message(sprintf("Written: %s", dta_path))

  # ── Build globals CSV
  globals_rows <- lapply(names(h_results), function(outcome_name) {
    res         <- h_results[[outcome_name]]
    if (is.null(res)) return(NULL)
    outcome_col <- gsub("\\{H\\}", h, OUTCOMES[[outcome_name]])

    # D*X base vars for this outcome (prespec + LASSO-selected)
    dx_vars_this_outcome <- unique(c(
      not_penalized,
      gsub("_d$", "", res$retained_dx_outcome)
    ))
    dx_cols_this_outcome <- paste0("d_x_", dx_vars_this_outcome)

    # Z*X base vars for this outcome: prespec (always forced) +
    # D*X-mirrored (forced in Part 2 Step 1) + LASSO-selected (penalised pool)
    zx_vars_this_outcome <- unique(c(
      not_penalized,
      gsub("_d$",  "", res$retained_dx_outcome),
      gsub("_iv$", "", res$retained_zx_exposure)
    ))

    # Model definitions (version argument passed to Stata subgroupboot):
    #   model1: homogeneous treatment effects — Z only in Stage 1, no Z*X;
    #           same Xlist in both stages
    #   model2: heterogeneous treatment effects — Z + Z*X in Stage 1,
    #           D*X added to Stage 2 to allow effect modification
    # Both models are exported to globals CSV and run in Stata.

    # Model 1: homogeneous — no D*X in Stage 2, no Z*X in Stage 1
    xlist_s2_m1 <- unique(c(
      res$union_main_effects,
      res$retained_xx_exposure
    ))
    xlist_s1_m1 <- xlist_s2_m1

    # Model 2: heterogeneous — D*X in Stage 2, Z*X in Stage 1
    xlist_s2_m2 <- unique(c(
      res$union_main_effects,
      res$retained_xx_exposure
    ))

    xlist_s1_m2 <- unique(c(
      res$union_main_effects,
      res$retained_xx_exposure
    ))
    zlist_s1_m2 <- paste0("z_x_stage1_", zx_vars_this_outcome)
    dxlist_s2_m2 <- dx_cols_this_outcome

    # no_iv: plain regression — X*X from outcome model (Part 2 Step 1),
    # no D*X, no retained_xx_exposure (that is exposure-model specific)
    xlist_s2_noiv <- unique(c(
      res$union_main_effects,
      res$retained_xx_outcome
    ))
    
    data.frame(
      outcome      = outcome_col,
      family       = OUTCOME_FAMILIES[[outcome_name]],
      xlist_s2_m1  = paste(xlist_s2_m1,                            collapse = " "),
      xlist_s1_m1  = paste(xlist_s1_m1,                            collapse = " "),
      xlist_s2_m2  = paste(xlist_s2_m2,                            collapse = " "),
      xlist_s1_m2  = paste(xlist_s1_m2,                            collapse = " "),
      zlist_s1_m2  = paste(zlist_s1_m2,                            collapse = " "),
      dxlist_s2_m2 = paste(dxlist_s2_m2, collapse = " "),
      xlist_s2_noiv  = paste(xlist_s2_noiv, collapse = " "), 
      stringsAsFactors = FALSE
    )
  })

  globals_df <- do.call(rbind, Filter(Negate(is.null), globals_rows))
  globals_path <- file.path(OUTPUT_DIR, sprintf("elective_%dd_globals.csv", h))
  write.csv(globals_df, globals_path, row.names = FALSE)
  message(sprintf("Written: %s", globals_path))
}

# for failed mortality outcome, take union of the selected covariates across all other outcomes at the same time horizon

build_mortality_union_row <- function(globals_path, horizon) {
  
  globals <- read.csv(globals_path, stringsAsFactors = FALSE)
  
  # Check if mortality row already exists — skip if so
  mortality_col <- sprintf("died_post_bypass_surg_%dd", horizon)
  if (any(globals$outcome == mortality_col)) {
    message(sprintf("Mortality row already exists at %dd — skipping", horizon))
    return(globals)
  }
  
  # Helper to take union of space-separated variable lists across all outcomes
  union_vars <- function(col) {
    paste(
      unique(unlist(strsplit(paste(globals[[col]], collapse = " "), " "))),
      collapse = " "
    )
  }
  
  mortality_row <- data.frame(
    outcome       = mortality_col,
    family        = "binomial",
    xlist_s2_m1   = union_vars("xlist_s2_m1"),
    xlist_s1_m1   = union_vars("xlist_s1_m1"),
    xlist_s2_m2   = union_vars("xlist_s2_m2"),
    xlist_s1_m2   = union_vars("xlist_s1_m2"),
    zlist_s1_m2   = union_vars("zlist_s1_m2"),
    dxlist_s2_m2  = union_vars("dxlist_s2_m2"),
    xlist_s2_noiv = union_vars("xlist_s2_noiv"),
    stringsAsFactors = FALSE
  )
  
  globals_updated <- rbind(globals, mortality_row)
  write.csv(globals_updated, globals_path, row.names = FALSE)
  message(sprintf("Mortality union row added and written: %s", globals_path))
  
  globals_updated
}

# Run for all three horizons
for (h in c(90, 180, 365)) {
  globals_path <- file.path(
    "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/lasso_outputs/",
    sprintf("elective_%dd_globals.csv", h)
  )
  build_mortality_union_row(globals_path, h)
}
