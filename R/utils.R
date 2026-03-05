# utils.R
#
# Shared internal utilities used across multiple modules.
# All functions here are internal (not exported).
#
# Functions:
#   read_qs_df()           — reads a .qs file, retaining only wanted columns
#   .standardise_names()   — generic snake_case column renamer (used at ingestion
#                            in hes_utils.R and nvr_utils.R)
#   .strip_id_prefix()     — strips a regex prefix from an ID vector; used by
#                            clean_HES_df_id_for_matching() and calculate_outcomes()

# =============================================================================
# 0. .qs file reader
# =============================================================================

#' Read a .qs file, selecting only wanted columns
#'
#' Warns if any requested columns are absent from the file — returns the
#' intersection rather than erroring, so callers see clearly what was dropped.
#'
#' @param file_path   Path to the `.qs` file.
#' @param wanted_cols Character vector of column names to retain.
#' @return A tibble containing only the columns in `wanted_cols` that exist
#'   in the file.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @importFrom qs qread
#' @keywords internal
read_qs_df <- function(file_path, wanted_cols) {
  df      <- qs::qread(file_path)
  missing <- setdiff(wanted_cols, names(df))

  if (length(missing) > 0L) {
    warning(
      "read_qs_df: the following wanted columns were not found: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
    wanted_cols <- intersect(wanted_cols, names(df))
  }

  df %>% dplyr::select(dplyr::all_of(wanted_cols))
}

# =============================================================================
# 1. Generic column name standardiser
# =============================================================================

#' Standardise data frame column names to snake_case
#'
#' Lowercases all column names, replaces any non-alphanumeric character
#' (including `.`, `:`, spaces) with `_`, strips leading/trailing underscores,
#' and truncates to 32 characters to satisfy Stata's hard column name limit.
#'
#' Called immediately after loading raw data in `hes_utils.R` and
#' `nvr_utils.R`. Never called mid-pipeline.
#'
#' @param data A data frame or tibble.
#' @return The same object with standardised column names.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_with
#' @importFrom stringr str_replace_all str_remove str_trunc
#' @keywords internal
.standardise_names <- function(data) {
  data %>%
    dplyr::rename_with(
      ~ .x %>%
        tolower() %>%
        stringr::str_replace_all("[^a-z0-9]+", "_") %>%
        stringr::str_remove("^_|_$") %>%
        stringr::str_trunc(32, ellipsis = "")
    )
}

# =============================================================================
# 2. ID prefix stripper
# =============================================================================

#' Strip a regex prefix from a character ID vector
#'
#' Used wherever a raw ID column needs its prefix removed before joining.
#' Centralised here so the stripping logic is defined once and consistent
#' across `clean_HES_df_id_for_matching()` and `calculate_outcomes()`.
#'
#' @param x      Character vector of IDs (e.g. `c("AB12345", "AB99999")`).
#' @param prefix Regex prefix to remove. Defaults to `"^AB"`.
#' @return Character vector with prefix stripped.
#'
#' @examples
#' \dontrun{
#'   .strip_id_prefix(c("AB001", "AB002"))  # returns c("001", "002")
#'   .strip_id_prefix(c("PT001"), "^PT")    # returns c("001")
#' }
#'
#' @importFrom stringr str_remove
#' @keywords internal
.strip_id_prefix <- function(x, prefix = "^AB") {
  stringr::str_remove(x, prefix)
}