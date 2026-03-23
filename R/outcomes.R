# outcomes.R
#
# Calculates post-operative outcomes for a TTE cohort:
#   - Readmission flags (binary, per time horizon)
#   - Length of stay (LOS): day-stay-exclusive and day-stay-inclusive,
#     at intervention level and per time horizon
#   - Mortality flags (directional; deaths before starting_point excluded)
#   - Days alive and out of hospital (DAOH): primary and Myles formulations
#
# Returns a SINGLE wide tibble — one row per patient — so analysts can inspect
# any column directly with View() or glimpse().
#
# Only calculate_outcomes() is exported. All helpers are internal (.prefix).

# =============================================================================
# 1. Validation
# =============================================================================

#' Validate intervention_name length and format
#'
#' Enforces the Stata 32-char column name constraint. The binding pattern is
#' `readmit_post_{name}_365d` (13 + len(name) + 5 chars), which must be <= 32,
#' giving a maximum name length of 14 characters.
#'
#' @param intervention_name Character scalar.
#' @return Invisibly returns `intervention_name` if valid; stops with a
#'   descriptive error if not.
#' @keywords internal
.validate_intervention_name <- function(intervention_name) {
  # Stata hard limit: 32 chars per column name.
  # Binding pattern: readmit_post_{name}_365d (13 + name + 5 = 18 + name)
  # Max name length: 32 - 18 = 14 chars
  max_len <- 14L

  if (!is.character(intervention_name) || length(intervention_name) != 1L) {
    stop(
      "`intervention_name` must be a single character string.",
      call. = FALSE
    )
  }
  if (!grepl("^[a-z][a-z0-9_]*$", intervention_name)) {
    stop(
      "`intervention_name` must be snake_case: lowercase letters, digits, ",
      "and underscores only, starting with a letter.\n",
      "  Got: \"", intervention_name, "\"",
      call. = FALSE
    )
  }
  n <- nchar(intervention_name)
  if (n > max_len) {
    longest <- paste0("readmit_post_", intervention_name, "_365d")
    stop(
      "`intervention_name` \"", intervention_name, "\" is ", n, " chars ",
      "(max is ", max_len, "). The longest generated column would be:\n",
      "  \"", longest, "\" (", nchar(longest), " chars), ",
      "exceeding Stata's 32-char limit.\n",
      "Use a shorter abbreviation.",
      call. = FALSE
    )
  }
  invisible(intervention_name)
}

#' Check that required columns are present in a data frame
#' @keywords internal
.check_cols <- function(df, required, df_name) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(
      "`", df_name, "` is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(df)
}

# =============================================================================
# 2. HES preparation: overlap merging
# =============================================================================

#' Merge overlapping or contiguous HES admission intervals within each patient
#'
#' HES transfers and multi-episode admissions can produce contiguous intervals
#' where `hes_discharge_date[n] == hes_admission_date[n+1]`. If left unmerged,
#' these shared boundary days are counted twice in LOS, producing inflated LOS
#' and potentially negative DAOH. This function collapses contiguous or
#' overlapping intervals into single intervals before any LOS computation.
#'
#' Contiguous is defined as `hes_admission_date[n+1] == hes_discharge_date[n]`
#' — these share a boundary day and are merged. Strictly adjacent
#' (`hes_admission_date[n+1] > hes_discharge_date[n]`) are kept separate.
#'
#' @param hes_long A tibble with at minimum `study_id`,
#'   `hes_admission_date` (Date), `hes_discharge_date` (Date). One row per
#'   admission after DISDATE cleaning in `hes_utils.R`.
#' @return A tibble in the same three-column format with contiguous and
#'   overlapping intervals merged. Additional columns from the input are
#'   dropped; join dates back from `cohort` downstream.
#' @keywords internal
.merge_overlapping_admissions <- function(hes_long) {
  hes_long %>%
    dplyr::select(study_id, hes_admission_date, hes_discharge_date) %>%
    dplyr::arrange(study_id, hes_admission_date) %>%
    dplyr::group_by(study_id) %>%
    dplyr::mutate(
      prev_dis     = dplyr::lag(hes_discharge_date),
      new_interval = is.na(prev_dis) | hes_admission_date > prev_dis,
      interval_id  = cumsum(new_interval)
    ) %>%
    dplyr::group_by(study_id, interval_id) %>%
    dplyr::summarise(
      hes_admission_date = min(hes_admission_date),
      hes_discharge_date = max(hes_discharge_date),
      .groups = "drop"
    ) %>%
    dplyr::select(-interval_id)
}

# =============================================================================
# 3. Wide pivot of readmission date columns (for inspection)
# =============================================================================

#' Pivot pre or post-intervention admissions from long to wide
#'
#' Produces interleaved admidate/disdate column pairs named
#' `{role}_{intervention_name}_{N}_admidate` /
#' `{role}_{intervention_name}_{N}_disdate`.
#'
#' Warns if the maximum readmissions per patient exceeds 20 — this usually
#' indicates unresolved overlapping admissions or a data quality issue. All
#' readmissions are returned regardless; no truncation occurs.
#'
#' @param hes_long Long tibble already filtered to the correct role and time
#'   window. Must contain `study_id`, `hes_admission_date`,
#'   `hes_discharge_date`.
#' @param intervention_name Character scalar.
#' @param role Character: `"pre"` or `"post"`.
#' @return Wide tibble, one row per patient. Patients with no admissions in
#'   this role have all-NA columns.
#' @keywords internal
.pivot_wide_admissions <- function(hes_long, intervention_name, role) {

  if (nrow(hes_long) == 0L) {
    return(tibble::tibble(study_id = character()))
  }

  df_numbered <- hes_long %>%
    dplyr::select(study_id, hes_admission_date, hes_discharge_date) %>%
    dplyr::rename(
      admidate = hes_admission_date,
      disdate  = hes_discharge_date
    ) %>%
    dplyr::arrange(study_id, admidate) %>%
    dplyr::group_by(study_id) %>%
    dplyr::mutate(read_no = dplyr::row_number()) %>%
    dplyr::ungroup()

  max_n <- max(df_numbered$read_no, na.rm = TRUE)

  col_prefix <- paste0(role, "_", intervention_name, "_")

  wide <- df_numbered %>%
    tidyr::pivot_wider(
      id_cols     = study_id,
      names_from  = read_no,
      values_from = c(admidate, disdate),
      names_glue  = paste0(col_prefix, "{read_no}_{.value}")
    )

  # Enforce interleaved admidate/disdate column ordering
  adm_cols <- paste0(col_prefix, seq_len(max_n), "_admidate")
  dis_cols <- paste0(col_prefix, seq_len(max_n), "_disdate")
  ordered  <- c("study_id", as.vector(rbind(adm_cols, dis_cols)))

  wide %>%
    dplyr::select(dplyr::all_of(ordered))
}

# =============================================================================
# 4. Intervention-level LOS (not horizon-specific)
# =============================================================================

#' @param proc_start_col Character scalar or NULL. Column holding procedure
#'   start date (i.e. `starting_point`). When supplied, two additional columns
#'   are computed: `{name}_proc_los_no` and `{name}_proc_los_w`, measuring LOS
#'   from procedure start → discharge (Option A). When NULL (default), only
#'   full admission LOS columns are produced.
.calc_intervention_los <- function(wide_df, intervention_name,
                                   proc_start_col = NULL) {
  adm_col      <- paste0(intervention_name, "_admidate")
  dis_col      <- paste0(intervention_name, "_disdate")
  no_ds_col    <- paste0(intervention_name, "_los_no")
  w_ds_col     <- paste0(intervention_name, "_los_w")
  proc_no_col  <- paste0(intervention_name, "_proc_los_no")
  proc_w_col   <- paste0(intervention_name, "_proc_los_w")

  result <- wide_df %>%
    dplyr::mutate(
      # Full admission LOS (admission date → discharge)
      !!no_ds_col := as.integer(.data[[dis_col]] - .data[[adm_col]]),
      !!w_ds_col  := dplyr::if_else(
        .data[[adm_col]] == .data[[dis_col]],
        1L,
        as.integer(.data[[dis_col]] - .data[[adm_col]])
      )
    )

  if (!is.null(proc_start_col)) {
    result <- result %>%
      dplyr::mutate(
        # Procedure LOS (procedure start date → discharge, Option A)
        # Identical to full admission LOS when admission date == procedure date.
        !!proc_no_col := as.integer(.data[[dis_col]] - .data[[proc_start_col]]),
        !!proc_w_col  := dplyr::if_else(
          .data[[proc_start_col]] == .data[[dis_col]],
          1L,
          as.integer(.data[[dis_col]] - .data[[proc_start_col]])
        )
      )
  }
  result
}

# =============================================================================
# 5. Horizon-specific outcomes
# =============================================================================

#' Compute all outcome columns for a single time horizon
#'
#' Operates on the merged long HES data to compute, for one horizon H:
#' - Binary readmission and pre-admission flags
#' - LOS components (post, pre if applicable, total) in both day-stay variants
#' - Mortality flag (directional — deaths before starting_point excluded)
#' - Primary DAOH and Myles DAOH
#'
#' All LOS columns use `_no` (same-day admissions excluded) and `_w`
#' (same-day admissions included as 1 day) suffixes consistently throughout.
#'
#' DAOH always uses `total_los_w` (day-stay-inclusive LOS). References:
#' Alexander et al. 2022 (doi:10.1111/ans.18099),
#' Donnelly et al. 2024 (doi:10.1161/JAHA.123.032321),
#' Harrington et al. 2024 (doi:10.1161/JAHA.122.028951),
#' Fanaroff et al. 2019 (doi:10.1161/CIRCOUTCOMES.118.004755),
#' Vaena et al. 2025 (doi:10.1038/s41598-025-14526-7).
#'
#' @param hes_merged Long tibble after deduplication. One row per
#'   `study_id` / `hes_admission_date` / `hes_discharge_date`.
#' @param cohort_dates Tibble with `study_id`, `starting_point`,
#'   `interv_admi_date`, `interv_date`.
#' @param mortality_clean Tibble with `study_id`, `death_date`.
#' @param horizon Integer. Follow-up window in days.
#' @param intervention_name Character scalar.
#' @param include_pre_intervention Logical.
#' @return Tibble with `study_id` and all horizon-specific columns.
#' @keywords internal
.compute_horizon_outcomes <- function(
  hes_merged,
  cohort_dates,
  mortality_clean,
  horizon,
  intervention_name,
  include_pre_intervention, 
  exclude_day_stay_readmissions = TRUE # TO DO: NEED TO CLARIFY WHAT THIS ARGUMENT IS FOR CLEARLY IN THE DOCS
) {
  h_suf <- paste0("_", horizon, "d")
  name  <- intervention_name

  # Filter HES to window [starting_point, starting_point + horizon], classify
   hes_in_window <- hes_merged %>%
    dplyr::inner_join(cohort_dates, by = "study_id") %>%
    dplyr::filter(
      # Include index admission even if hes_admission_date < starting_point
      # (occurs when procedure start date > admission date). All other
      # admissions must start within [starting_point, starting_point + horizon].
      hes_admission_date == interv_admi_date |
        (hes_admission_date >= starting_point &
           hes_admission_date <= starting_point + horizon)
    ) %>%
    dplyr::mutate(
      # Role: index admission identified by interv_admi_date (NVR admission
      # date), not interv_date (procedure start date). Pre/post boundary
      # remains interv_date.
      role = dplyr::case_when(
        hes_admission_date == interv_admi_date ~ "intervention",
        hes_admission_date <  interv_date      ~ "pre",
        hes_admission_date >  interv_admi_date      ~ "post"
      ),
      # LOS for intervention: starts from starting_point (procedure start date)
      # not admission date — Option A: only days in hospital after procedure
      # count. For non-elective where admission == procedure date, this is
      # identical to admission → discharge.
      .effective_start = dplyr::if_else(
        role == "intervention",
        starting_point,
        hes_admission_date
      ),
      los_no = as.integer(hes_discharge_date - .effective_start),
      los_w  = dplyr::if_else(
        .effective_start == hes_discharge_date,
        1L,
        as.integer(hes_discharge_date - .effective_start)
      )
    ) %>%
    dplyr::select(-.effective_start)
  
  # Warn if post-intervention admissions per patient unusually high for
  # this horizon — indicates overlapping admissions or data quality issues.
  max_post_n <- hes_in_window %>%
    dplyr::filter(role == "post") %>%
    dplyr::count(study_id) %>%
    dplyr::pull(n) %>%
    max(na.rm = TRUE)

  if (length(max_post_n) > 0L && max_post_n > 20L) {
    warning(
      "At horizon ", horizon, "d: maximum post-intervention admissions ",
      "per patient is ", max_post_n, " (threshold: 20). ",
      "Check for unresolved overlapping admissions in the HES data.",
      call. = FALSE
    )
  }

  # Summarise by patient and role
  role_summary <- hes_in_window %>%
    dplyr::filter(!is.na(role)) %>%
    dplyr::group_by(study_id, role) %>%
    dplyr::summarise(
      los_no    = sum(los_no, na.rm = TRUE),
      los_w     = sum(los_w,  na.rm = TRUE),
      n_admissions = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols     = study_id,
      names_from  = role,
      values_from = c(los_no, los_w, n_admissions),
      names_glue  = "{role}_{.value}"
    )

  # Ensure all expected columns exist (absent roles become 0)
  expected_cols <- c(
    "intervention_los_no", "intervention_los_w",
    "post_los_no",         "post_los_w",
    "post_n_admissions"
  )
  if (include_pre_intervention) {
    expected_cols <- c(
      expected_cols,
      "pre_los_no", "pre_los_w", "pre_n_admissions"
    )
  }
  for (col in expected_cols) {
    if (!col %in% names(role_summary)) role_summary[[col]] <- 0L
  }

  # Left-join onto full cohort; patients with no HES in window get 0 LOS
  result <- cohort_dates %>%
    dplyr::select(study_id, starting_point) %>%
    dplyr::left_join(role_summary, by = "study_id") %>%
    dplyr::mutate(
      dplyr::across(dplyr::any_of(expected_cols), ~ dplyr::coalesce(.x, 0L))
    )

  # Compute total LOS (intervention + pre + post), cap at horizon
  result <- result %>%
    dplyr::mutate(
      .pre_no = if (include_pre_intervention) pre_los_no else 0L,
      .pre_w  = if (include_pre_intervention) pre_los_w  else 0L,
      !!paste0("total_los_no", h_suf) := pmin(
        intervention_los_no + .pre_no + post_los_no, horizon
      ),
      !!paste0("total_los_w", h_suf) := pmin(
        intervention_los_w + .pre_w + post_los_w, horizon
      )
    ) %>%
    dplyr::select(-.pre_no, -.pre_w)

  # Binary readmission flags
  # Binary readmission flag — excludes day-stay admissions per
  # exclude_day_stay_readmissions. LOS is unaffected: _no naturally
  # contributes 0 for day-stays, _w counts them as 1 day.
  post_n_for_flag <- hes_in_window %>%
    dplyr::filter(role == "post") %>%
    { if (exclude_day_stay_readmissions)
        dplyr::filter(., hes_admission_date != hes_discharge_date)
      else . } %>%
    dplyr::count(study_id, name = "post_n_flag")

  result <- result %>%
    dplyr::left_join(post_n_for_flag, by = "study_id") %>%
    dplyr::mutate(
      !!paste0("readmit_post_", name, h_suf) := dplyr::if_else(
        dplyr::coalesce(post_n_flag, 0L) > 0L, 1L, 0L
      )
    ) %>%
    dplyr::select(-post_n_flag)

  if (include_pre_intervention) {
    result <- result %>%
      dplyr::mutate(
        !!paste0("admit_pre_", name, h_suf) := dplyr::if_else(
          dplyr::coalesce(pre_n_admissions, 0L) > 0L, 1L, 0L
        )
      )
  }

  # Rename internal role-summary LOS columns to output names with
  # intervention name and horizon suffix.
  rename_map <- stats::setNames(
    c("post_los_no", "post_los_w"),
    c(
      paste0("post_", name, "_los_no", h_suf),
      paste0("post_", name, "_los_w",  h_suf)
    )
  )
  if (include_pre_intervention) {
    rename_map <- c(rename_map, stats::setNames(
      c("pre_los_no", "pre_los_w"),
      c(
        paste0("pre_", name, "_los_no", h_suf),
        paste0("pre_", name, "_los_w",  h_suf)
      )
    ))
  }
  result <- result %>%
    dplyr::rename(dplyr::all_of(rename_map))

  # Mortality flag — directional (deaths before starting_point excluded).
  # BUG FIX vs. original scripts: original used abs(died_date - time_zero)
  # which incorrectly counts deaths occurring *before* time zero if the gap
  # fell within the horizon. Fix: require death_date >= starting_point.
  result <- result %>%
    dplyr::left_join(mortality_clean, by = "study_id") %>%
    dplyr::mutate(
      !!paste0("died_post_", name, h_suf) := dplyr::if_else(
        !is.na(death_date) &
          death_date >= starting_point &
          death_date <= starting_point + horizon,
        1L, 0L
      )
    )

  # DAOH (primary formulation — all 5 references):
  # Survivors: horizon - total_los_w
  # Deceased:  horizon - total_los_w - days_dead_before_end_of_followup
  #
  # DAOH (Myles formulation, secondary output):
  # Deceased:  0
  # Survivors: horizon - total_los_w
  #
  # Both formulations use total_los_w (day-stay-inclusive) throughout.
  total_w_col    <- paste0("total_los_w", h_suf)
  died_col       <- paste0("died_post_", name, h_suf)
  daoh_col       <- paste0("daoh_", name, h_suf)
  daoh_myles_col <- paste0("daoh_myles_", name, h_suf)

  result <- result %>%
    dplyr::mutate(
      .follow_up_end = starting_point + horizon,
      .days_dead     = dplyr::if_else(
        .data[[died_col]] == 1L & !is.na(death_date),
        as.integer(.follow_up_end - death_date),
        0L
      ),
      !!daoh_col := dplyr::case_when(
        .data[[died_col]] == 1L ~
          as.integer(horizon - .data[[total_w_col]] - .days_dead),
        .data[[died_col]] == 0L ~
          as.integer(horizon - .data[[total_w_col]]),
        TRUE ~ NA_integer_
      ),
      !!daoh_myles_col := dplyr::case_when(
        .data[[died_col]] == 1L ~ 0L,
        .data[[died_col]] == 0L ~
          as.integer(horizon - .data[[total_w_col]]),
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::select(-.follow_up_end, -.days_dead, -death_date)

  neg_n <- result %>%
    dplyr::filter(
      !is.na(.data[[daoh_col]]),
      .data[[daoh_col]] < 0L
    ) %>%
    nrow()
  if (neg_n > 0L) {
    warning(
      neg_n, " patient(s) have negative ", daoh_col,
      " at horizon ", horizon, " days ",
      "(day-stay-inclusive LOS exceeded horizon). ",
      "Likely cause: overlapping admissions not fully resolved by ",
      ".merge_overlapping_admissions(). ",
      "Flooring to 0 — matches original script behaviour (pmax(DAOH, 0)).",
      call. = FALSE
    )
    result <- result %>%
      dplyr::mutate(
        !!daoh_col       := pmax(.data[[daoh_col]],       0L, na.rm = FALSE),
        !!daoh_myles_col := pmax(.data[[daoh_myles_col]], 0L, na.rm = FALSE)
      )
  }

  # Return only horizon-specific columns
  result %>%
    dplyr::select(
      study_id,
      dplyr::ends_with(h_suf),
      -dplyr::starts_with("intervention_los_"),
      -starting_point
    )
}

# =============================================================================
# 6. Main exported function
# =============================================================================

#' Calculate post-operative outcomes for a TTE cohort
#'
#' Calculates readmission, length of stay (LOS), mortality, and days alive
#' and out of hospital (DAOH) over one or more time horizons. Returns a
#' single wide tibble — one row per patient — suitable for direct inspection
#' with `View()` or `glimpse()`.
#'
#' @param cohort A tibble produced by `build_cohort()`, one row per patient.
#'   Must include a `starting_point` Date column (or the column named by
#'   `starting_point_col`) and a column named by `intervention_date_col`.
#' @param hes_admissions A long tibble of HES APC admissions cleaned by
#'   `clean_hes()`. Must contain `study_id`, `hes_admission_date` (Date),
#'   `hes_discharge_date` (Date).
#' @param mortality A tibble of mortality records cleaned by `clean_hes()`.
#'   Must contain `study_id` and `death_date` (Date).
#' @param intervention_name Character scalar naming the intervention, e.g.
#'   `"bypass_surg"`. Used as a prefix/infix in all output column names.
#'   Maximum 14 characters (Stata 32-char constraint — see Details).
#' @param intervention_admission_date_col Character scalar. Column in `cohort`
#'   holding the NVR admission date for the index intervention. Used solely to
#'   identify the corresponding HES admission row (matched on
#'   `hes_admission_date`). Defaults to `"nvr_admission_date"`.
#' @param intervention_date_col Character scalar. Column in `cohort` holding
#'   the NVR procedure start date. Used as the pre/post boundary — admissions
#'   before this date are classified as pre-intervention, admissions on or
#'   after as post-intervention. Defaults to `"nvr_procedure_start_date"`.
#' @param starting_point_col Character scalar. Column in `cohort` representing
#'   time zero. Defaults to `"starting_point"`, which `build_cohort()` always
#'   creates. The elective bypass study used `"valid_op_date"` as time zero;
#'   the non-elective study used the intervention admission date itself.
#' @param time_horizons Integer vector of follow-up windows in days. Defaults
#'   to `c(90, 180, 365)`. Any vector of positive integers is accepted.
#' @param include_pre_intervention Logical. Should admissions occurring between
#'   `starting_point` and `intervention_date_col` be tracked? Defaults to
#'   `TRUE`. Set to `FALSE` when starting_point equals intervention date (e.g.
#'   non-elective studies with no pre-intervention window).
#' @param mortality_id_prefix Character scalar. Prefix to strip from `study_id`
#'   in the mortality table before joining to `cohort`. Defaults to `"AB"`.
#'   Set to `""` if no prefix stripping is needed.
#'
#' @return A single wide tibble, one row per patient:
#'
#'   **Fixed columns (not horizon-specific):**
#'   - `study_id`, `starting_point`
#'   - `{name}_admidate`, `{name}_disdate`
#'   - `{name}_los_no` (intervention LOS, admission → discharge, day-stays = 0)
#'   - `{name}_los_w`  (intervention LOS, admission → discharge, day-stays = 1)
#'   - `{name}_proc_los_no` (intervention LOS, procedure start → discharge,
#'     day-stays = 0; identical to `_los_no` when admission == procedure date)
#'   - `{name}_proc_los_w`  (intervention LOS, procedure start → discharge,
#'     day-stays = 1)
#'   - `post_{name}_N_admidate`, `post_{name}_N_disdate` (all post-intervention
#'     readmissions within `max(time_horizons)` days — for inspection)
#'   - `pre_{name}_N_admidate`, `pre_{name}_N_disdate` (if
#'     `include_pre_intervention = TRUE`)
#'   - `death_date`
#'
#'   **Per-horizon columns (repeated for each H in `time_horizons`):**
#'   - `readmit_post_{name}_{H}d`  — any post-intervention readmission (binary)
#'   - `admit_pre_{name}_{H}d`     — any pre-intervention admission (binary)
#'   - `died_post_{name}_{H}d`     — death on or after starting_point (binary)
#'   - `post_{name}_los_no_{H}d`   — post-intervention LOS, day-stays = 0
#'   - `post_{name}_los_w_{H}d`    — post-intervention LOS, day-stays = 1
#'   - `pre_{name}_los_no_{H}d`    — pre-intervention LOS, day-stays = 0
#'   - `pre_{name}_los_w_{H}d`     — pre-intervention LOS, day-stays = 1
#'   - `total_los_no_{H}d`         — total LOS, day-stays = 0, capped at H
#'   - `total_los_w_{H}d`          — total LOS, day-stays = 1, capped at H
#'   - `daoh_{name}_{H}d`          — DAOH, primary formulation
#'   - `daoh_myles_{name}_{H}d`    — DAOH, Myles formulation
#'
#' @details
#' **DAOH (primary):** For survivors: `horizon - total_los_w`. For deaths
#' within the horizon: `horizon - total_los_w - days_dead_before_followup_end`.
#' DAOH uses day-stay-inclusive LOS (`total_los_w`) throughout. References:
#' Alexander et al. (2022, \doi{10.1111/ans.18099}),
#' Donnelly et al. (2024, \doi{10.1161/JAHA.123.032321}),
#' Harrington et al. (2024, \doi{10.1161/JAHA.122.028951}),
#' Fanaroff et al. (2019, \doi{10.1161/CIRCOUTCOMES.118.004755}),
#' Vaena et al. (2025, \doi{10.1038/s41598-025-14526-7}).
#'
#' **DAOH (Myles):** Deaths = 0; survivors = `horizon - total_los_w`.
#' (Myles et al., BMJ Open 2017, \doi{10.1136/bmjopen-2017-015828})
#'
#' **Overlap merging:** Currently disabled — overlapping admissions are not
#' merged. DAOH is floored at 0 to match original script behaviour. See TODO
#' in source for planned implementation.
#' 
#' **Index admission LOS:** When `starting_point_col` differs from
#' `intervention_admission_date_col` (e.g. TTEs where procedure
#' start date is after admission date), the index admission LOS is counted from
#' `starting_point` → discharge, not admission → discharge. Only days in
#' hospital after the procedure start contribute to LOS and DAOH outcome calc. For
#' studies where admission date == procedure start date, this is
#' identical to the full admission LOS.
#'
#' **Mortality:** Deaths before `starting_point` are never counted regardless
#' of how close they fall to the horizon. The `abs()` in earlier analysis
#' scripts was a bug that could count pre-baseline deaths.
#'
#' **Stata 32-char constraint:** `intervention_name` is limited to 14 chars.
#' All LOS columns use `_no` / `_w` suffixes consistently.
#'
#' @examples
#' \dontrun{
#'   results <- calculate_outcomes(
#'     cohort                   = bypass_cohort,
#'     hes_admissions           = hes_apc_clean,
#'     mortality                = hes_mortality_clean,
#'     intervention_name        = "bypass_surg",
#'     intervention_date_col    = "nvr_procedure_start_date",
#'     intervention_admission_date_col = "nvr_admission_date",
#'     starting_point_col       = "starting_point", # usually time zero, but can be any date column in `cohort`
#'     time_horizons            = c(90, 180, 365),
#'     include_pre_intervention = TRUE,
#'     mortality_id_prefix      = "AB"
#'   )
#'   glimpse(results)
#' }
#'
#'   summarise ungroup arrange lag coalesce if_else case_when rename
#'   any_of all_of starts_with ends_with across row_number n distinct
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom purrr reduce map
#' @importFrom tibble tibble
#' @export
calculate_outcomes <- function(
  cohort,
  hes_admissions,
  mortality,
  intervention_name,
  intervention_admission_date_col = "nvr_admission_date",
  intervention_date_col    = "nvr_procedure_start_date",
  starting_point_col       = "starting_point",
  time_horizons            = c(90, 180, 365),
  include_pre_intervention = TRUE,
  exclude_day_stay_readmissions = TRUE,
  mortality_id_prefix      = "AB"
) {

  # ---- Input validation -------------------------------------------------------

  .validate_intervention_name(intervention_name)

  stopifnot(
    is.data.frame(cohort),
    is.data.frame(hes_admissions),
    is.data.frame(mortality),
    is.character(intervention_admission_date_col),
    length(intervention_admission_date_col) == 1L,
    is.character(intervention_date_col), length(intervention_date_col) == 1L,
    is.character(starting_point_col),    length(starting_point_col) == 1L,
    is.numeric(time_horizons),           all(time_horizons > 0),
    is.logical(include_pre_intervention),
    length(include_pre_intervention) == 1L,
    is.character(mortality_id_prefix),   length(mortality_id_prefix) == 1L
  )

  time_horizons <- as.integer(time_horizons)

  .check_cols(
    cohort, c("study_id", starting_point_col, intervention_admission_date_col, intervention_date_col),
    "cohort"
  )
  .check_cols(
    hes_admissions,
    c("study_id", "hes_admission_date", "hes_discharge_date"),
    "hes_admissions"
  )
  .check_cols(mortality, c("study_id", "death_date"), "mortality")

  # ---- Prepare cohort date reference -----------------------------------------

  cohort_dates <- cohort %>%
    dplyr::select(
      study_id,
      starting_point = dplyr::all_of(starting_point_col),
      interv_admi_date = dplyr::all_of(intervention_admission_date_col),
      interv_date    = dplyr::all_of(intervention_date_col)
    ) %>%
    dplyr::mutate(
      starting_point = as.Date(starting_point),
      interv_admi_date = as.Date(interv_admi_date),
      interv_date    = as.Date(interv_date)
    )

  # ---- Prepare HES -----------------------------------------------------------

  # Filter to cohort patients and merge overlapping intervals.
  # DISDATE sentinel cleaning (1800-01-01, 1801-01-01) is done upstream in
  # hes_utils.R::clean_hes() — not repeated here.
  hes_cohort <- hes_admissions %>%
    dplyr::filter(study_id %in% cohort$study_id) %>%
    dplyr::mutate(
      hes_admission_date = as.Date(hes_admission_date),
      hes_discharge_date = as.Date(hes_discharge_date)
    )

  # TODO: overlap merging disabled — matches original script behaviour which
  # did not merge overlapping admissions and instead floored DAOH at 0 post hoc.
  # Re-enable once Option B is implemented: merge within the index admission
  # only, not across the index/post boundary (which currently inflates index
  # LOS by absorbing contiguous post-index episodes into the index span).
  hes_merged <- hes_cohort %>%
    dplyr::distinct(study_id, hes_admission_date, hes_discharge_date)

  # ---- Prepare mortality -----------------------------------------------------

  # Strip mortality_id_prefix from study_id before joining.
  # e.g. "AB12345" -> "12345" to match cohort study_id format.
  mortality_clean <- mortality %>%
    dplyr::mutate(
      study_id   = .strip_id_prefix(study_id, paste0("^", mortality_id_prefix)),
      death_date = as.Date(death_date)
    ) %>%
    dplyr::select(study_id, death_date)

  # ---- Build base wide tibble ------------------------------------------------

  max_h <- max(time_horizons)

  # Intervention admission (one row per patient, renamed to intervention_name)
  interv_wide <- hes_merged %>%
    dplyr::inner_join(
      cohort_dates %>%
        dplyr::select(study_id, interv_admi_date, starting_point),
      by = "study_id"
    ) %>%
    dplyr::filter(hes_admission_date == interv_admi_date) %>%
    dplyr::distinct(study_id, .keep_all = TRUE) %>%
    dplyr::select(
      study_id, hes_admission_date, hes_discharge_date, starting_point
    ) %>%
    dplyr::rename(
      !!paste0(intervention_name, "_admidate") := hes_admission_date,
      !!paste0(intervention_name, "_disdate")  := hes_discharge_date
    ) %>%
    .calc_intervention_los(intervention_name, proc_start_col = "starting_point") %>%
    dplyr::select(-starting_point)

  # Add surgical discharge date to cohort_dates — needed to correctly classify
  # post-intervention admissions as strictly after surgical discharge.
  cohort_dates <- cohort_dates %>%
    dplyr::left_join(
      interv_wide %>%
        dplyr::select(
          study_id,
          interv_disdate = !!paste0(intervention_name, "_disdate")
        ),
      by = "study_id"
    )

  # Post-intervention readmission date columns (all within max horizon)
  post_hes <- hes_merged %>%
    dplyr::inner_join(cohort_dates, by = "study_id") %>%
    dplyr::filter(
      hes_admission_date >  interv_admi_date,
      hes_admission_date <= starting_point + max_h
    )
  post_wide <- .pivot_wide_admissions(post_hes, intervention_name, "post")

  # Pre-intervention readmission date columns (if applicable)
  if (include_pre_intervention) {
    pre_hes <- hes_merged %>%
      dplyr::inner_join(cohort_dates, by = "study_id") %>%
      dplyr::filter(
        hes_admission_date >= starting_point,
        hes_admission_date <  interv_date
      )
    pre_wide <- .pivot_wide_admissions(pre_hes, intervention_name, "pre")
  }

  # Assemble base: cohort + intervention + pre (opt) + post + death_date
  base_wide <- cohort %>%
    dplyr::select(study_id, dplyr::all_of(starting_point_col)) %>%
    dplyr::rename(starting_point = dplyr::all_of(starting_point_col)) %>%
    dplyr::left_join(interv_wide, by = "study_id")

  if (include_pre_intervention) {
    base_wide <- base_wide %>%
      dplyr::left_join(pre_wide, by = "study_id")
  }

  base_wide <- base_wide %>%
    dplyr::left_join(post_wide, by = "study_id") %>%
    dplyr::left_join(mortality_clean, by = "study_id")

  # ---- Compute horizon-specific outcomes -------------------------------------

  horizon_tibbles <- purrr::map(time_horizons, function(h) {
    .compute_horizon_outcomes(
      hes_merged               = hes_merged,
      cohort_dates             = cohort_dates,
      mortality_clean          = mortality_clean,
      horizon                  = h,
      intervention_name        = intervention_name,
      include_pre_intervention = include_pre_intervention, 
      exclude_day_stay_readmissions = TRUE
    )
  })

  # Join all horizon-specific columns onto the base tibble
  purrr::reduce(
    horizon_tibbles,
    ~ dplyr::left_join(.x, .y, by = "study_id"),
    .init = base_wide
  )
}
