# cohort.R
#
# Cohort construction helpers for NVR / HES-linked studies.
# Covers inclusion/exclusion filtering, deduplication, and HES linkage
# validation. All functions are internal except where noted.
#
# Functions:
#   perform_inclusion()                        — apply named inclusion criteria
#   perform_exclusion()                        — apply named exclusion criteria
#   check_duplicates()                         — report duplicates on a key col
#   handle_duplicates()                        — resolve duplicates via tiebreaker
#   filter_NVR_patients_with_HES_record()      — linkage step 1: any HES record
#   filter_NVR_patients_with_matched_admission() — linkage step 2: matched date
#   subset_HES_to_NVR_cohort()                 — subset HES to final cohort

# =============================================================================
# 1. INCLUSION / EXCLUSION
# =============================================================================

#' Apply inclusion criteria to a data frame
#'
#' Pass named `quo()` conditions via `do.call()`. Each condition identifies
#' rows to **keep**. NAs in referenced columns are treated as failing the
#' condition and are dropped. Prints before/after counts and a breakdown of
#' removed rows per referenced column.
#'
#' @param df    A tibble / data frame.
#' @param ...   Named `quo()` conditions identifying rows to keep.
#' @param label Label prefix for console messages (e.g. `"NVR_include"`).
#' @return Filtered tibble.
#'
#' @examples
#' \dontrun{
#'   INCLUSION_CRITERIA <- list(
#'     procedure_type = quo(ProcedureType %in% c("Lower-limb Bypass")),
#'     clti_only      = quo(`Indications:AmpIndicationCode` %in% c(2, 4))
#'   )
#'   do.call(perform_inclusion,
#'           c(list(df = nvr_df, label = "NVR_include"), INCLUSION_CRITERIA))
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang eval_tidy as_label
#' @importFrom dplyr filter count mutate
#' @importFrom tidyr replace_na
#' @importFrom purrr walk pmap
#' @export
perform_inclusion <- function(df, ..., label = "inclusion") {
  conditions <- list(...)

  for (i in seq_along(conditions)) {
    cond      <- conditions[[i]]
    cond_name <- if (nzchar(names(conditions)[i])) names(conditions)[i] else rlang::as_label(cond)

    keep       <- rlang::eval_tidy(cond, data = df)
    keep       <- tidyr::replace_na(keep, FALSE)
    removed_df <- dplyr::filter(df, !keep)
    df         <- dplyr::filter(df, keep)

    message(sprintf("  [INCLUDE | %s | %s] Kept %d  |  Removed %d",
                    label, cond_name, nrow(df), nrow(removed_df)))

    breakdown_cols <- intersect(all.vars(cond), names(df))
    purrr::walk(breakdown_cols, function(col) {
      message(sprintf("    Breakdown by '%s':", col))
      removed_df %>%
        dplyr::count(.data[[col]], sort = TRUE, name = "n") %>%
        dplyr::mutate(val = ifelse(is.na(.data[[col]]), "<missing>",
                                   as.character(.data[[col]]))) %>%
        purrr::pmap(function(val, n, ...) message(sprintf("      %-30s : %d", val, n)))
    })
  }
  df
}


#' Apply exclusion criteria to a data frame
#'
#' Pass named `quo()` conditions via `do.call()`. Each condition identifies
#' rows to **keep**. Unlike [perform_inclusion()], rows where referenced
#' columns are `NA` are always retained. Prints before/after counts and a
#' breakdown of removed rows per referenced column.
#'
#' @param df    A tibble / data frame.
#' @param ...   Named `quo()` conditions identifying rows to keep.
#' @param label Label prefix for console messages (e.g. `"NVR_exclude"`).
#' @return Filtered tibble.
#'
#' @examples
#' \dontrun{
#'   EXCLUSION_CRITERIA <- list(
#'     english_lsoa = quo(grepl("^E", LSOA))
#'   )
#'   do.call(perform_exclusion,
#'           c(list(df = nvr_df, label = "NVR_exclude"), EXCLUSION_CRITERIA))
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang eval_tidy as_label
#' @importFrom dplyr filter count mutate
#' @importFrom tidyr replace_na
#' @importFrom purrr walk pmap
#' @export
perform_exclusion <- function(df, ..., label = "exclusion") {
  conditions <- list(...)

  for (i in seq_along(conditions)) {
    cond      <- conditions[[i]]
    cond_name <- if (nzchar(names(conditions)[i])) names(conditions)[i] else rlang::as_label(cond)

    breakdown_cols <- intersect(all.vars(cond), names(df))

    result  <- rlang::eval_tidy(cond, data = df)
    na_mask <- if (length(breakdown_cols) > 0)
      rowSums(sapply(breakdown_cols, function(col) is.na(df[[col]]))) > 0
    else
      rep(FALSE, nrow(df))

    keep       <- tidyr::replace_na(result, FALSE) | na_mask
    removed_df <- dplyr::filter(df, !keep)
    df         <- dplyr::filter(df, keep)

    message(sprintf("  [EXCLUDE | %s | %s] Removed %d  |  Retained %d",
                    label, cond_name, nrow(removed_df), nrow(df)))

    purrr::walk(breakdown_cols, function(col) {
      message(sprintf("    Breakdown by '%s':", col))
      removed_df %>%
        dplyr::count(.data[[col]], sort = TRUE, name = "n") %>%
        dplyr::mutate(val = ifelse(is.na(.data[[col]]), "<missing>",
                                   as.character(.data[[col]]))) %>%
        purrr::pmap(function(val, n, ...) message(sprintf("      %-30s : %d", val, n)))
    })
  }
  df
}

# =============================================================================
# 2. DUPLICATES
# =============================================================================

#' Check for duplicates on a key column and print a summary
#'
#' Returns all duplicate rows (both original and duplicate) invisibly for
#' inspection. Does not modify the data frame.
#'
#' @param df     A tibble / data frame.
#' @param id_col Name of the ID column to check for duplicates.
#' @return Invisibly returns a tibble of all rows involved in duplicates.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select group_by filter ungroup arrange n_distinct
#' @export
check_duplicates <- function(df, id_col) {
  n_exact   <- sum(duplicated(df))
  n_id_dups <- sum(duplicated(df %>% dplyr::select(dplyr::all_of(id_col))))
  n_unique  <- dplyr::n_distinct(df[[id_col]])

  message("  ── Duplicate Report ──────────────────────────────────────")
  message(sprintf("  Total rows                  : %d", nrow(df)))
  message(sprintf("  Exact duplicate rows        : %d", n_exact))
  message(sprintf("  Duplicate %s values   : %d", id_col, n_id_dups))
  message(sprintf("  Unique %s values      : %d", id_col, n_unique))
  message("  ──────────────────────────────────────────────────────────")

  dup_rows <- df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data[[id_col]])

  invisible(dup_rows)
}


#' Resolve duplicate NVR patient records via sequential tiebreaker logic
#'
#' Applies four sequential steps to reduce to one row per patient. Any records
#' unresolvable after all steps are dropped and logged.
#'
#' Steps:
#' 1. Remove exact duplicate rows.
#' 2. Keep row(s) with earliest procedure AND admission date.
#' 3. Keep row(s) with most comorbidities (highest pipe count).
#' 4. Keep row(s) with highest indication + Fontaine + ASA score.
#'    True ties after step 4 are dropped entirely.
#'
#' @param df              A tibble / data frame.
#' @param id_col          Name of the patient ID column.
#' @param procedure_date  Name of the procedure date column.
#' @param admission_date  Name of the admission date column.
#' @param comorbidity_col Name of the pipe-separated comorbidity string column.
#' @param indication_col  Name of the amputation indication code column.
#' @param fontaine_col    Name of the Fontaine score column.
#' @param asa_col         Name of the ASA score column.
#' @return De-duplicated tibble with one row per patient.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct group_by mutate filter select ungroup n
#' @importFrom stringr str_count fixed
#' @importFrom purrr walk
#' @export
handle_duplicates <- function(df,
                               id_col,
                               procedure_date  = "NvrEpisode:ProcedureStartDate",
                               admission_date  = "NvrEpisode:AdmissionDate",
                               comorbidity_col = "RiskScores:Comorbidities",
                               indication_col  = "Indications:AmpIndicationCode",
                               fontaine_col    = "Indications:PadFontaineCode",
                               asa_col         = "RiskScores:ASA") {
  before <- nrow(df)

  # ── Step 1: exact duplicates ─────────────────────────────────────────────
  df <- dplyr::distinct(df)
  message(sprintf(
    "  [Dedup Step 1 | exact duplicates] Rows after: %d (removed %d)",
    nrow(df), before - nrow(df)
  ))

  # ── Step 2: earliest procedure + admission date ──────────────────────────
  before_step <- nrow(df)
  df <- df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      .wins      = (.data[[procedure_date]] == min(.data[[procedure_date]], na.rm = TRUE)) +
                   (.data[[admission_date]] == min(.data[[admission_date]], na.rm = TRUE)),
      .best_wins = max(.wins, na.rm = TRUE),
      .best_proc = min(.data[[procedure_date]][.wins == .best_wins], na.rm = TRUE),
      .best_adm  = min(.data[[admission_date]][
                         .wins == .best_wins &
                         .data[[procedure_date]] == .best_proc
                       ], na.rm = TRUE)
    ) %>%
    dplyr::filter(
      .wins == .best_wins,
      .data[[procedure_date]] == .best_proc,
      .data[[admission_date]] == .best_adm
    ) %>%
    dplyr::select(-.wins, -.best_wins, -.best_proc, -.best_adm) %>%
    dplyr::ungroup()
  message(sprintf(
    "  [Dedup Step 2 | earliest procedure/admission] Rows after: %d (removed %d)",
    nrow(df), before_step - nrow(df)
  ))

  # ── Step 3: most comorbidities ───────────────────────────────────────────
  before_step <- nrow(df)
  df <- df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      .pipe_n = dplyr::if_else(
        is.na(.data[[comorbidity_col]]),
        0L,
        stringr::str_count(.data[[comorbidity_col]], stringr::fixed("|"))
      )
    ) %>%
    dplyr::filter(.pipe_n == max(.pipe_n)) %>%
    dplyr::select(-.pipe_n) %>%
    dplyr::ungroup()
  message(sprintf(
    "  [Dedup Step 3 | most comorbidities] Rows after: %d (removed %d)",
    nrow(df), before_step - nrow(df)
  ))

  # ── Step 4: highest indication + Fontaine + ASA ──────────────────────────
  before_step <- nrow(df)
  df <- df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      .wins      = (.data[[indication_col]] == max(.data[[indication_col]], na.rm = TRUE)) +
                   (.data[[fontaine_col]]   == max(.data[[fontaine_col]],   na.rm = TRUE)) +
                   (.data[[asa_col]]        == max(.data[[asa_col]],        na.rm = TRUE)),
      .best_wins = max(.wins, na.rm = TRUE)
    ) %>%
    dplyr::filter(.wins == .best_wins) %>%
    dplyr::select(-.wins, -.best_wins) %>%
    dplyr::ungroup()
  message(sprintf(
    "  [Dedup Step 4 | highest indication/fontaine/ASA] Rows after: %d (removed %d)",
    nrow(df), before_step - nrow(df)
  ))

  # ── Drop true ties ────────────────────────────────────────────────────────
  true_ties      <- df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  n_tie_patients <- dplyr::n_distinct(true_ties[[id_col]])

  if (n_tie_patients > 0) {
    message(sprintf(
      "  [Dedup | true ties dropped] %d rows across %d patients unresolvable",
      nrow(true_ties), n_tie_patients
    ))
    purrr::walk(unique(true_ties[[id_col]]), function(pid) {
      message(sprintf("    %s : %d rows", pid, sum(true_ties[[id_col]] == pid)))
    })
    df <- df %>% dplyr::filter(!(.data[[id_col]] %in% true_ties[[id_col]]))
  }

  message(sprintf(
    "  [Dedup | final] Rows after: %d (removed %d total from input)",
    nrow(df), before - nrow(df)
  ))
  df
}

# =============================================================================
# 3. NVR / HES LINKAGE
# =============================================================================

#' Step 1: Keep NVR patients who have at least one HES record
#'
#' Semi-joins NVR onto HES on patient ID. Patients with no HES record at all
#' are removed and counted.
#'
#' @param nvr_df     NVR tibble (post-deduplication).
#' @param hes_df     HES tibble with cleaned ID column (post
#'   [clean_HES_df_id_for_matching()]).
#' @param nvr_id_col Name of the NVR patient ID column.
#' @param hes_id_col Name of the cleaned HES ID column (e.g. `"STUDY_ID_clean"`).
#' @return Filtered NVR tibble.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate semi_join
#' @export
filter_NVR_patients_with_HES_record <- function(nvr_df,
                                                 hes_df,
                                                 nvr_id_col,
                                                 hes_id_col) {
  before <- nrow(nvr_df)

  nvr_df <- nvr_df %>%
    dplyr::mutate(!!nvr_id_col := as.character(.data[[nvr_id_col]]))
  hes_df <- hes_df %>%
    dplyr::mutate(!!hes_id_col := as.character(.data[[hes_id_col]]))

  nvr_df <- dplyr::semi_join(
    nvr_df, hes_df,
    by = stats::setNames(hes_id_col, nvr_id_col)
  )

  message(sprintf(
    "  [NVR/HES linkage | Step 1] Patients with HES record: %d  |  removed: %d",
    nrow(nvr_df), before - nrow(nvr_df)
  ))
  nvr_df
}


#' Step 2: Keep NVR patients where the NVR admission date matches a HES ADMIDATE
#'
#' For each patient, checks whether any HES admission date equals the NVR
#' admission date. Patients with no matching date across all their HES rows
#' are removed.
#'
#' @param nvr_df            NVR tibble (output of
#'   [filter_NVR_patients_with_HES_record()]).
#' @param hes_df            HES tibble with cleaned ID column.
#' @param nvr_id_col        Name of the NVR patient ID column.
#' @param hes_id_col        Name of the cleaned HES ID column.
#' @param nvr_admission_col Name of the NVR admission date column.
#' @param hes_admission_col Name of the HES admission date column.
#'   Defaults to `"ADMIDATE"`.
#' @return Filtered NVR tibble.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate inner_join group_by summarise filter pull
#' @export
filter_NVR_patients_with_matched_admission <- function(nvr_df,
                                                        hes_df,
                                                        nvr_id_col,
                                                        hes_id_col,
                                                        nvr_admission_col,
                                                        hes_admission_col = "ADMIDATE") {
  before <- nrow(nvr_df)

  matched_ids <- hes_df %>%
    dplyr::select(dplyr::all_of(c(hes_id_col, hes_admission_col))) %>%
    dplyr::mutate(!!hes_id_col := as.character(.data[[hes_id_col]])) %>%
    dplyr::inner_join(
      nvr_df %>%
        dplyr::select(dplyr::all_of(c(nvr_id_col, nvr_admission_col))) %>%
        dplyr::mutate(!!nvr_id_col := as.character(.data[[nvr_id_col]])),
      by = stats::setNames(nvr_id_col, hes_id_col)
    ) %>%
    dplyr::group_by(.data[[hes_id_col]]) %>%
    dplyr::summarise(
      any_date_match = any(
        .data[[hes_admission_col]] == .data[[nvr_admission_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(any_date_match) %>%
    dplyr::pull(.data[[hes_id_col]])

  nvr_df <- nvr_df %>%
    dplyr::filter(as.character(.data[[nvr_id_col]]) %in% matched_ids)

  message(sprintf(
    "  [NVR/HES linkage | Step 2] Patients with matched admission date: %d  |  removed: %d",
    nrow(nvr_df), before - nrow(nvr_df)
  ))
  nvr_df
}


#' Step 3: Subset HES to final NVR cohort
#'
#' Returns **all** HES rows for patients in the final NVR cohort — not just
#' the index admission. Downstream scripts filter further as needed.
#'
#' @param hes_df     HES tibble.
#' @param nvr_df     Final NVR tibble (post-linkage validation).
#' @param nvr_id_col Name of the NVR patient ID column.
#' @param hes_id_col Name of the cleaned HES ID column.
#' @return Filtered HES tibble.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter n_distinct
#' @export
subset_HES_to_NVR_cohort <- function(hes_df,
                                      nvr_df,
                                      nvr_id_col,
                                      hes_id_col) {
  nvr_ids <- as.character(nvr_df[[nvr_id_col]])

  hes_df <- hes_df %>%
    dplyr::filter(as.character(.data[[hes_id_col]]) %in% nvr_ids)

  message(sprintf(
    "  [HES subset] Rows retained: %d  (unique patients: %d)",
    nrow(hes_df), dplyr::n_distinct(hes_df[[hes_id_col]])
  ))
  hes_df
}
