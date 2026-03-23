# hes_utils.R
#
# HES-specific cleaning helpers (internal).
# All functions here are called at data ingestion — never mid-pipeline.
#
# Functions:
#   clean_HES_admidate()           — fixes sentinel / invalid ADMIDATE values
#   clean_HES_df_id_for_matching() — strips "AB" prefix from STUDY_ID to
#                                    produce a matchable _clean column
#
# Higher-level loading + renaming (DISDATE cleaning, column renaming to
# internal names) is done inline in experiment_run.R — see comment scaffolds
# there for the expected output schema.

# =============================================================================
# 1. ADMIDATE cleaning
# =============================================================================

#' Clean sentinel and invalid ADMIDATE values in HES data
#'
#' For rows where ADMIDATE is a sentinel value, computes the minimum valid
#' EPISTART within the same patient-admission group and uses it as a
#' replacement. Rows where both ADMIDATE and all EPISTARTs within the group
#' are sentinel are dropped entirely. Prints full diagnostics before and after.
#'
#' @param hes_df         A tibble of raw HES data.
#' @param id_col         Name of the patient ID column used for grouping.
#'   Defaults to `"STUDY_ID_FOR_NVR_MERGE"`.
#' @param admidate_col   Name of the admission date column to clean.
#'   Defaults to `"ADMIDATE"`.
#' @param epistart_col   Name of the episode start date column used as
#'   fallback. Defaults to `"EPISTART"`.
#' @param sentinel_dates Character vector of sentinel date strings to treat
#'   as invalid. Defaults to `c("1800-01-01", "1801-01-01")`.
#' @return Cleaned tibble with sentinel ADMIDATE values replaced or rows
#'   dropped where no valid EPISTART exists. All other columns preserved.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate ungroup filter if_else select
#' @importFrom purrr walk
#' @keywords internal
clean_HES_admidate <- function(hes_df,
                               id_col,
                               admidate_col   = "ADMIDATE",
                               epistart_col   = "EPISTART",
                               sentinel_dates = c("1800-01-01", "1801-01-01")) {

  hes_df[[admidate_col]] <- as.Date(hes_df[[admidate_col]])
  hes_df[[epistart_col]] <- as.Date(hes_df[[epistart_col]])
  sentinel_vec           <- as.Date(sentinel_dates)
  before                 <- nrow(hes_df)

  # ---- Diagnostics (before) --------------------------------------------------

  message(sprintf("  [HES | ADMIDATE] Total rows                  : %d", nrow(hes_df)))
  message(sprintf("  [HES | ADMIDATE] Missing ADMIDATE            : %d", sum(is.na(hes_df[[admidate_col]]))))
  message(sprintf("  [HES | ADMIDATE] Min ADMIDATE                : %s", min(hes_df[[admidate_col]], na.rm = TRUE)))
  message(sprintf("  [HES | ADMIDATE] Max ADMIDATE                : %s", max(hes_df[[admidate_col]], na.rm = TRUE)))
  purrr::walk(sentinel_dates, ~ {
    n <- sum(hes_df[[admidate_col]] == as.Date(.x), na.rm = TRUE)
    message(sprintf("  [HES | ADMIDATE] Sentinel %s           : %d rows", .x, n))
  })

  # ---- Compute min valid EPISTART per patient per admission ------------------

  hes_df <- hes_df %>%
    dplyr::group_by(.data[[id_col]], .data[[admidate_col]]) %>%
    dplyr::mutate(
      .min_epistart = suppressWarnings(
        min(
          .data[[epistart_col]][!(.data[[epistart_col]] %in% sentinel_vec)],
          na.rm = TRUE
        )
      )
    ) %>%
    dplyr::ungroup()

  # ---- Drop rows where ADMIDATE is sentinel AND no valid EPISTART exists -----

  hes_df <- hes_df %>%
    dplyr::filter(
      !(.data[[admidate_col]] %in% sentinel_vec &
        (is.na(.min_epistart) | is.infinite(.min_epistart)))
    )

  message(sprintf("  [HES | ADMIDATE] Dropped (no valid EPISTART) : %d rows", before - nrow(hes_df)))

  # ---- Replace sentinel ADMIDATE with min valid EPISTART --------------------

  n_to_replace <- sum(hes_df[[admidate_col]] %in% sentinel_vec)
  hes_df <- hes_df %>%
    dplyr::mutate(
      !!admidate_col := dplyr::if_else(
        .data[[admidate_col]] %in% sentinel_vec,
        .min_epistart,
        .data[[admidate_col]]
      )
    ) %>%
    dplyr::select(-.min_epistart)

  # ---- Diagnostics (after) --------------------------------------------------

  message(sprintf("  [HES | ADMIDATE] Replaced with min EPISTART  : %d rows", n_to_replace))
  message(sprintf("  [HES | ADMIDATE] Min ADMIDATE after clean     : %s", min(hes_df[[admidate_col]], na.rm = TRUE)))
  message(sprintf("  [HES | ADMIDATE] Max ADMIDATE after clean     : %s", max(hes_df[[admidate_col]], na.rm = TRUE)))
  message(sprintf("  [HES | ADMIDATE] Remaining sentinel ADMIDATE  : %d", sum(hes_df[[admidate_col]] %in% sentinel_vec)))
  message(sprintf("  [HES | ADMIDATE] Complete. Rows remaining     : %d", nrow(hes_df)))

  hes_df
}

# =============================================================================
# 2. STUDY_ID prefix cleaning
# =============================================================================

#' Strip a string prefix from HES STUDY_ID to allow matching with NVR
#'
#' Creates a new `<id_col>_clean` column alongside the original so the raw
#' ID is preserved. Calls [.strip_id_prefix()] from `utils.R`.
#'
#' @param hes_df           A tibble of HES data.
#' @param id_col           Name of the raw HES ID column (e.g. `"STUDY_ID"`).
#' @param prefix_to_remove Regex prefix to strip. Defaults to `"^AB"`.
#' @return Input tibble with a new `<id_col>_clean` column appended.
#'   The original `id_col` column is preserved unchanged.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @keywords internal
clean_HES_df_id_for_matching <- function(hes_df,
                                         id_col,
                                         prefix_to_remove = "^AB") {
  new_col <- paste0(id_col, "_clean")
  hes_df %>%
    dplyr::mutate(
      !!new_col := .strip_id_prefix(.data[[id_col]], prefix_to_remove)
    )
}

# =============================================================================
# 3. ICD-10 flag generation
# =============================================================================

#' Flag comorbidities from a concatenated ICD-10 diagnosis string
#'
#' For each group in `groups`, searches a comma-separated ICD-10 string column
#' for any code matching the group's prefixes and adds a binary `0`/`1`
#' indicator column named `<group>_yn`.
#'
#' @param data   A tibble / data frame.
#' @param groups Character vector of group names to flag (must be keys in
#'   `dict`).
#' @param dict   Named list mapping group names to ICD-10 code prefixes
#'   (e.g. `icd10_groups_charlson`).
#' @param src    Name of the column containing comma-separated ICD-10 codes.
#'   Defaults to `"DIAG_4_CONCAT"`.
#' @return Tibble with one additional binary column per group: `<group>_yn`.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfc
#' @importFrom stringr str_replace_all str_detect
#' @importFrom tibble tibble
#' @keywords internal
add_flags_from_concat <- function(data, groups, dict, src = "DIAG_4_CONCAT") {
  stopifnot(src %in% names(data))
  s <- toupper(as.character(data[[src]]))
  s <- stringr::str_replace_all(s, "\\s+", "")

  flags <- purrr::map_dfc(groups, function(g) {
    prefs <- dict[[g]]
    if (is.null(prefs)) stop(sprintf("Group '%s' not in dict.", g), call. = FALSE)
    pat <- regex(paste0("(?:^|,)(?:", paste(prefs, collapse = "|"), ")[A-Z0-9]*(?:,|$)"))
    tibble::tibble(!!paste0(g, "_yn") := as.integer(stringr::str_detect(s, pat)))
  })

  dplyr::bind_cols(data, flags)
}

# =============================================================================
# 4. IMD decile to quintile conversion
# =============================================================================

#' Convert HES IMD04 decile labels to quintiles and resolve conflicts at
#' patient level
#'
#' Normalises the IMD04_DECILE string column, maps deciles to quintiles
#' (Q1 least deprived to Q5 most deprived), then resolves conflicts across
#' episodes by taking the highest deprivation value per patient.
#'
#' @param hes_df       A tibble containing HES episode data.
#' @param id_col       Name of the patient ID column. Defaults to
#'   `"STUDY_ID_clean"`.
#' @param decile_col   Name of the IMD decile column. Defaults to
#'   `"IMD04_DECILE"`.
#' @param quintile_col Name of the new quintile column to create. Defaults to
#'   `"IMD_quintile"`.
#' @return Tibble with new factor column `IMD_quintile` (Q1 least deprived to
#'   Q5 most deprived), one row per patient, `decile_col` dropped.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate group_by ungroup select distinct case_when
#' @keywords internal
HES_imd_decile_to_quintile <- function(hes_df,
                                        id_col       = "STUDY_ID_clean",
                                        decile_col   = "IMD04_DECILE",
                                        quintile_col = "IMD_quintile") {
  hes_df %>%
    dplyr::mutate(
      !!decile_col := tolower(as.character(.data[[decile_col]])),
      !!decile_col := gsub("less deprived 10%", "least deprived 10%",
                            .data[[decile_col]]),
      !!quintile_col := dplyr::case_when(
        .data[[decile_col]] %in% c("least deprived 10%",
                                    "less deprived 10-20%")   ~ "Q1 (least deprived)",
        .data[[decile_col]] %in% c("less deprived 20-30%",
                                    "less deprived 30-40%")   ~ "Q2",
        .data[[decile_col]] %in% c("less deprived 40-50%",
                                    "more deprived 40-50%")   ~ "Q3",
        .data[[decile_col]] %in% c("more deprived 20-30%",
                                    "more deprived 30-40%")   ~ "Q4",
        .data[[decile_col]] %in% c("more deprived 10-20%",
                                    "most deprived 10%")      ~ "Q5 (most deprived)"
      ),
      !!quintile_col := factor(
        .data[[quintile_col]],
        levels = c("Q1 (least deprived)", "Q2", "Q3", "Q4",
                    "Q5 (most deprived)")
      )
    ) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      !!quintile_col := {
        mx <- max(as.integer(.data[[quintile_col]]), na.rm = TRUE)
        factor(
          ifelse(is.infinite(mx), NA_integer_, mx),
          levels = 1:5,
          labels = c("Q1 (least deprived)", "Q2", "Q3", "Q4",
                      "Q5 (most deprived)")
        )
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::all_of(decile_col)) %>%
    dplyr::distinct()
}