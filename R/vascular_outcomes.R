# vascular_outcomes.R
#
# Limb-specific vascular outcomes for TTE studies using HES OPCS-4 coded data:
#   - Index limb revascularisation (ILR)
#   - Index limb major amputation (ILMA)
#   - Amputation-free survival (AFS)
#
# Only flag_limb_events() and compute_limb_outcomes() are exported.
# All helpers are internal (.prefix).
#
# Caller responsibilities (NOT handled here):
#   - Pre-filter HES to post-index episodes before passing to flag_limb_events()
#   - Extract index limb laterality from NVR via extract_index_limb_events_from_nvr()
#     (script-level function, not in this package — NVR data dict may change)
#   - Supply study-specific OPCS-4 code prefix vectors (amp_prefixes, revasc_prefixes)

# =============================================================================
# 1. Package constant: OPCS-4 laterality codes
# =============================================================================

#' OPCS-4 laterality qualifier codes
#'
#' Standard OPCS-4 codes used to qualify the operated side when coded
#' adjacently to a procedure code:
#' - `Z941` = Bilateral
#' - `Z942` = Right side
#' - `Z943` = Left side
#'
#' Used by `.extract_opcs_event_and_laterality()` and
#' `.match_opcs_event_laterality()`.
#'
#' @keywords internal
.opcs_laterality_codes <- c(
  bilateral = "Z941",
  right     = "Z942",
  left      = "Z943"
)

# =============================================================================
# 2. Row-level OPCS-4 event extraction
# =============================================================================

#' Extract first matching OPCS-4 procedure and adjacent laterality code
#'
#' Scans a vector of OPCS-4 procedure codes for the first code matching any
#' of `target_prefixes`. When found, searches up to `laterality_window`
#' positions ahead in the code vector for an OPCS-4 laterality qualifier
#' (`Z941`, `Z942`, `Z943`).
#'
#' Designed for row-wise use inside `rowwise() %>% mutate()` on a HES FCE
#' tibble where each `OPERTN_01`...`OPERTN_24` column holds one OPCS-4 code.
#'
#' @param codes Character vector of OPCS-4 codes for one episode (e.g.
#'   `c_across(all_of(opertn_cols))`). `NA` values are tolerated.
#' @param dates Date vector of corresponding operation dates, same length as
#'   `codes` (e.g. `c_across(all_of(opdate_cols))`).
#' @param target_prefixes Character vector of OPCS-4 code prefixes to match
#'   (e.g. `c("X09", "L161")`). Matched via `startsWith()`.
#' @param use_laterality Logical. If `TRUE` (default), searches for an
#'   adjacent laterality qualifier after the matched procedure code. If
#'   `FALSE`, returns `NA_character_` for laterality — use for procedure
#'   types where laterality is not coded (e.g. some endovascular procedures).
#' @param laterality_window Integer. Number of positions ahead of the matched
#'   code to search for a laterality qualifier. Default `4L`. Values beyond
#'   this threshold can be noisy due to unrelated codes intervening.
#'
#' @return A named list with two elements:
#'   - `date`: the `Date` of the matched procedure, or `NA` if no match
#'   - `laterality`: the first laterality code found within the window, or
#'     `NA_character_` if none found (or `use_laterality = FALSE`)
#'
#' @keywords internal
.extract_opcs_event_and_laterality <- function(
  codes,
  dates,
  target_prefixes,
  use_laterality    = TRUE,
  laterality_window = 4L
) {
  for (i in seq_along(codes)) {
    if (!is.na(codes[i]) && any(startsWith(codes[i], target_prefixes))) {
      event_date <- dates[i]

      if (!use_laterality) {
        return(list(date = event_date, laterality = NA_character_))
      }

      window_idx <- seq(i + 1L, min(i + laterality_window, length(codes)))
      window     <- codes[window_idx]
      lat_match  <- window[!is.na(window) & window %in% .opcs_laterality_codes]
      lat        <- if (length(lat_match) > 0L) lat_match[1L] else NA_character_

      return(list(date = event_date, laterality = lat))
    }
  }
  list(date = NA_Date_, laterality = NA_character_)
}

# =============================================================================
# 3. Laterality matching
# =============================================================================

#' Check whether an OPCS-4 event laterality matches the index limb side
#'
#' Bilateral qualifies as a match regardless of index side. A missing index
#' side or event laterality always returns `FALSE`.
#'
#' Index side coding (NVR `IndicationSideCode`):
#' - `1` = Right
#' - `2` = Left
#' - `3` = Bilateral
#'
#' OPCS-4 laterality codes (see [.opcs_laterality_codes]):
#' - `Z941` = Bilateral
#' - `Z942` = Right
#' - `Z943` = Left
#'
#' @param index_side Integer scalar. NVR index limb side code (1, 2, or 3).
#' @param event_lat Character scalar. OPCS-4 laterality code from the event
#'   episode (one of `Z941`, `Z942`, `Z943`, or `NA`).
#'
#' @return Logical scalar. `TRUE` if the event is on the same side as the
#'   index limb (or either side is bilateral); `FALSE` otherwise.
#'
#' @keywords internal
.match_opcs_event_laterality <- function(index_side, event_lat) {
  if (is.na(index_side) || is.na(event_lat)) return(FALSE)
  # Bilateral on either side = match
  if (index_side == 3L || event_lat == .opcs_laterality_codes["bilateral"]) {
    return(TRUE)
  }
  (index_side == 1L && event_lat == .opcs_laterality_codes["right"]) ||
  (index_side == 2L && event_lat == .opcs_laterality_codes["left"])
}

# =============================================================================
# 4. Main event-flagging function
# =============================================================================

#' Identify earliest post-index amputation and revascularisation events
#'
#' Searches pre-filtered post-index HES FCEs for amputation and
#' revascularisation events using OPCS-4 code prefix matching. When
#' `use_laterality = TRUE`, only same-side events (matched against
#' `index_side_df`) are retained; when `FALSE`, the earliest event regardless
#' of laterality is returned.
#'
#' Returns one row per patient with the earliest qualifying event date for
#' each event type.
#'
#' @param hes_post_index A long tibble of HES FCEs, **pre-filtered by the
#'   caller to episodes occurring after the index procedure date**. Must
#'   contain `study_id` plus the columns named in `opertn_cols` and
#'   `opdate_cols`.
#' @param amp_prefixes Character vector of OPCS-4 code prefixes for
#'   amputation events (e.g. `"X09"`).
#' @param revasc_prefixes Character vector of OPCS-4 code prefixes for
#'   revascularisation events (e.g. `c("L161", "L501", ...)`).
#' @param index_side_df A tibble with columns `study_id` (character) and
#'   `index_side` (integer: 1 = Right, 2 = Left, 3 = Bilateral). Required
#'   when `use_laterality = TRUE`; ignored (and may be `NULL`) when
#'   `use_laterality = FALSE`.
#' @param use_laterality Logical. If `TRUE` (default), restricts to same-side
#'   events using `index_side_df`. If `FALSE`, returns the earliest matching
#'   event regardless of laterality — use for procedure types where laterality
#'   is not coded.
#' @param opertn_cols Character vector of OPCS-4 operation code column names.
#'   Defaults to `paste0("OPERTN_", sprintf("%02d", 1:24))`.
#' @param opdate_cols Character vector of operation date column names,
#'   parallel to `opertn_cols`.
#'   Defaults to `paste0("OPDATE_", sprintf("%02d", 1:24))`.
#' @param laterality_window Integer. Positions ahead of matched code to search
#'   for an OPCS-4 laterality qualifier. Default `4L`.
#'
#' @return A tibble with one row per `study_id` and two Date columns:
#'   - `amp_date`: earliest qualifying amputation event date, or `NA`
#'   - `revasc_date`: earliest qualifying revascularisation event date, or `NA`
#'
#' @details
#' The caller is responsible for pre-filtering `hes_post_index` to episodes
#' occurring strictly after the index procedure date. This function does not
#' enforce a lower date bound.
#'
#' Sentinel dates (`1800-01-01`, `1801-01-01`) in operation date columns
#' should be cleaned upstream (see `hes_utils.R`). Any remaining sentinel
#' dates are silently treated as `NA` via `as.Date()` coercion failures.
#'
#' @examples
#' \dontrun{
#'   amp_prefixes    <- "X09"
#'   revasc_prefixes <- c("L161", "L501", "L511")
#'
#'   index_side_df <- extract_index_limb_events_from_nvr(nvr_cohort)
#'
#'   hes_post <- hes_all %>% filter(hes_admission_date > index_procedure_date)
#'
#'   limb_events <- flag_limb_events(
#'     hes_post_index  = hes_post,
#'     amp_prefixes    = amp_prefixes,
#'     revasc_prefixes = revasc_prefixes,
#'     index_side_df   = index_side_df
#'   )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rowwise mutate select filter group_by summarise ungroup
#'   left_join c_across all_of if_else
#' @importFrom tidyr replace_na
#' @export
flag_limb_events <- function(
  hes_post_index,
  amp_prefixes,
  revasc_prefixes,
  index_side_df     = NULL,
  use_laterality    = TRUE,
  opertn_cols       = paste0("OPERTN_", sprintf("%02d", 1:24)),
  opdate_cols       = paste0("OPDATE_", sprintf("%02d", 1:24)),
  laterality_window = 4L
) {

  # ---- Input validation ------------------------------------------------------

  stopifnot(
    is.data.frame(hes_post_index),
    is.character(amp_prefixes),    length(amp_prefixes)    > 0L,
    is.character(revasc_prefixes), length(revasc_prefixes) > 0L,
    is.logical(use_laterality),    length(use_laterality)  == 1L,
    is.integer(laterality_window) || is.numeric(laterality_window),
    laterality_window >= 1L
  )
  laterality_window <- as.integer(laterality_window)

  .check_cols(hes_post_index, c("study_id", opertn_cols, opdate_cols),
              "hes_post_index")

  if (use_laterality) {
    if (is.null(index_side_df)) {
      stop(
        "`index_side_df` must be supplied when `use_laterality = TRUE`.\n",
        "Pass `NULL` and set `use_laterality = FALSE` for procedures ",
        "where laterality is not coded.",
        call. = FALSE
      )
    }
    .check_cols(index_side_df, c("study_id", "index_side"), "index_side_df")
  }

  sentinel_dates <- as.Date(c("1800-01-01", "1801-01-01"))

  # ---- Step 1: row-wise event extraction -------------------------------------

  events <- hes_post_index %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(opdate_cols), ~ as.Date(.x))
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      .amp    = list(.extract_opcs_event_and_laterality(
        codes             = c_across(dplyr::all_of(opertn_cols)),
        dates             = c_across(dplyr::all_of(opdate_cols)),
        target_prefixes   = amp_prefixes,
        use_laterality    = use_laterality,
        laterality_window = laterality_window
      )),
      .revasc = list(.extract_opcs_event_and_laterality(
        codes             = c_across(dplyr::all_of(opertn_cols)),
        dates             = c_across(dplyr::all_of(opdate_cols)),
        target_prefixes   = revasc_prefixes,
        use_laterality    = use_laterality,
        laterality_window = laterality_window
      )),
      amp_date    = as.Date(.amp$date),
      amp_lat     = .amp$laterality,
      revasc_date = as.Date(.revasc$date),
      revasc_lat  = .revasc$laterality
    ) %>%
    dplyr::select(-dplyr::starts_with(".amp"), -dplyr::starts_with(".revasc")) %>%
    dplyr::ungroup() %>%
    # Clean sentinel dates that survived upstream cleaning
    dplyr::mutate(
      amp_date    = dplyr::if_else(amp_date    %in% sentinel_dates,
                                   NA_Date_, amp_date),
      revasc_date = dplyr::if_else(revasc_date %in% sentinel_dates,
                                   NA_Date_, revasc_date)
    ) %>%
    # Keep only FCEs where at least one event was found
    dplyr::filter(!is.na(amp_date) | !is.na(revasc_date))

  # ---- Step 2: laterality matching and earliest same-side event --------------

  if (use_laterality) {
    events <- events %>%
      dplyr::left_join(index_side_df, by = "study_id") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        amp_same_side    = .match_opcs_event_laterality(index_side, amp_lat),
        revasc_same_side = .match_opcs_event_laterality(index_side, revasc_lat)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(study_id) %>%
      dplyr::summarise(
        amp_date    = min(amp_date[amp_same_side],       na.rm = TRUE),
        revasc_date = min(revasc_date[revasc_same_side], na.rm = TRUE),
        .groups     = "drop"
      )
  } else {
    events <- events %>%
      dplyr::group_by(study_id) %>%
      dplyr::summarise(
        amp_date    = min(amp_date,    na.rm = TRUE),
        revasc_date = min(revasc_date, na.rm = TRUE),
        .groups     = "drop"
      )
  }

  # min() over an all-NA group returns Inf — coerce back to NA
  events %>%
    dplyr::mutate(
      amp_date    = dplyr::if_else(is.infinite(amp_date),    NA_Date_, amp_date),
      revasc_date = dplyr::if_else(is.infinite(revasc_date), NA_Date_, revasc_date)
    )
}

# =============================================================================
# 5. Horizon-specific limb outcomes
# =============================================================================

#' Compute index limb revascularisation, amputation, and AFS outcomes
#'
#' Takes a cohort tibble, the output of [flag_limb_events()], and a mortality
#' tibble, and computes:
#' - Binary ILR and ILMA flags at each time horizon
#' - Uncapped AFS time-to-event and event indicator for survival modelling
#'
#' @param cohort A tibble with one row per patient. Must contain `study_id`
#'   and the column named by `start_date_col`.
#' @param limb_events Output of [flag_limb_events()]. Tibble with `study_id`,
#'   `amp_date`, `revasc_date` (all Date).
#' @param mortality Tibble with `study_id` and `death_date` (Date).
#' @param start_date_col Character scalar. Column in `cohort` representing
#'   time zero — the date from which days-to-event are measured. Typically the
#'   procedure start date.
#' @param censoring_date Scalar Date. Study end date; patients event-free
#'   beyond this date are censored here.
#' @param time_horizons Integer vector of follow-up windows in days for the
#'   binary ILR/ILMA outcomes. Default `c(90, 180, 365)`.
#'
#' @return A wide tibble, one row per patient, with columns:
#'   - `study_id`
#'   - `ilr_{H}d` — binary (0/1): revascularisation within H days of
#'     `start_date_col` (one column per horizon)
#'   - `ilma_{H}d` — binary (0/1): amputation within H days (one per horizon)
#'   - `afs_days` — integer: days from `start_date_col` to first of
#'     (amputation, death, censoring); uncapped, for Cox/KM models
#'   - `afs_event` — integer (0/1): 1 if amputation or death occurred before
#'     `censoring_date`; 0 if censored
#'
#' @details
#' `afs_days` and `afs_event` are uncapped — they are not restricted to any
#' horizon. Use `time_horizons` only to control which binary columns are
#' produced.
#'
#' AFS event = amputation OR death, whichever occurs first before
#' `censoring_date`. Patients with neither are censored at `censoring_date`.
#'
#' @examples
#' \dontrun{
#'   limb_outcomes <- compute_limb_outcomes(
#'     cohort          = non_elective_cohort,
#'     limb_events     = limb_events,
#'     mortality       = mortality_clean,
#'     start_date_col  = "nvr_procedure_start_date",
#'     censoring_date  = as.Date("2024-03-31"),
#'     time_horizons   = c(90, 180, 365)
#'   )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate select across all_of if_else
#' @importFrom purrr reduce
#' @export
compute_limb_outcomes <- function(
  cohort,
  limb_events,
  mortality,
  start_date_col,
  censoring_date,
  time_horizons  = c(90L, 180L, 365L)
) {

  # ---- Input validation ------------------------------------------------------

  stopifnot(
    is.data.frame(cohort),
    is.data.frame(limb_events),
    is.data.frame(mortality),
    is.character(start_date_col), length(start_date_col) == 1L,
    inherits(censoring_date, "Date"), length(censoring_date) == 1L,
    is.numeric(time_horizons), all(time_horizons > 0)
  )
  time_horizons <- as.integer(time_horizons)

  .check_cols(cohort,      c("study_id", start_date_col),            "cohort")
  .check_cols(limb_events, c("study_id", "amp_date", "revasc_date"), "limb_events")
  .check_cols(mortality,   c("study_id", "death_date"),              "mortality")

  # ---- Assemble base table ---------------------------------------------------

  base <- cohort %>%
    dplyr::select(study_id, start_date = dplyr::all_of(start_date_col)) %>%
    dplyr::mutate(start_date = as.Date(start_date)) %>%
    dplyr::left_join(limb_events, by = "study_id") %>%
    dplyr::left_join(mortality,   by = "study_id") %>%
    dplyr::mutate(
      days_to_amp    = as.integer(amp_date    - start_date),
      days_to_revasc = as.integer(revasc_date - start_date),
      days_to_death  = as.integer(death_date  - start_date),
      days_to_censor = as.integer(censoring_date - start_date)
    )

  # ---- AFS: uncapped time-to-event (for Cox / KM) ----------------------------
  # Event = amputation OR death, whichever first
  # Censored at censoring_date if neither occurs before study end
  base <- base %>%
    dplyr::mutate(
      afs_event = as.integer(
        (!is.na(days_to_amp)   & days_to_amp   <= days_to_censor) |
        (!is.na(days_to_death) & days_to_death <= days_to_censor)
      ),
      afs_days = as.integer(
        pmin(days_to_amp, days_to_death, days_to_censor, na.rm = TRUE)
      )
    )

  # ---- Binary horizon outcomes -----------------------------------------------

  for (h in time_horizons) {
    base <- base %>%
      dplyr::mutate(
        # ILR: revascularisation within horizon
        !!paste0("ilr_", h, "d")  := as.integer(
          !is.na(days_to_revasc) & days_to_revasc <= h
        ),
        # ILMA: amputation within horizon
        !!paste0("ilma_", h, "d") := as.integer(
          !is.na(days_to_amp) & days_to_amp <= h
        )
      )
  }

  # ---- Return clean output ---------------------------------------------------

  horizon_cols <- c(
    paste0("ilr_",  time_horizons, "d"),
    paste0("ilma_", time_horizons, "d")
  )

  base %>%
    dplyr::select(
      study_id,
      dplyr::all_of(horizon_cols),
      afs_days,
      afs_event
    )
}
