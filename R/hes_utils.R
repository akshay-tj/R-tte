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
