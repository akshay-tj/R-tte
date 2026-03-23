# nvr_utils.R
#
# NVR-specific column cleaning helpers (internal).
# All functions here are called at data ingestion — never mid-pipeline.
#
# Functions:
#   standardise_nvr_names()          — two-step snake_case renamer for NVR data
#   separate_out_nvr_comorbidities() — pivots pipe-delimited comorbidity string
#                                      into binary indicator columns

# =============================================================================
# Internal: generic snake_case renamer (shared with hes_utils.R)
# Defined in utils.R — documented here for reference:
#   .standardise_names(data) — lowercases, replaces non-alphanumeric with _,
#   strips leading/trailing _, truncates to 32 chars.
# =============================================================================

# Known raw NVR column → internal name mappings.
# Applied AFTER .standardise_names() has normalised to snake_case, so keys
# here are the post-normalisation forms (e.g. "nvrepisode_admissiondate"),
# not the raw originals (e.g. "NvrEpisode.AdmissionDate").
#
# Add new mappings here as additional NVR columns enter the pipeline.
.NVR_COLUMN_MAP <- c(
  # Identifiers
  nvrepisode_admissiondate = "nvr_admission_date",

  # Procedure / episode dates — add as needed
  nvrepisode_proceduredate = "nvr_procedure_date",

  # Operator / surgeon — add as needed
  nvrepisode_surgeoncode   = "nvr_surgeon_code",
  nvrepisode_hospitalcode  = "nvr_hospital_code"
)

#' Standardise NVR column names to internal snake_case convention
#'
#' Renames columns in a raw NVR tibble to the internal naming convention used
#' throughout the package. Applied immediately after loading raw NVR data —
#' never called mid-pipeline.
#'
#' Two-step process:
#' 1. Generic normalisation via [.standardise_names()]: lowercase, replace
#'    non-alphanumeric characters with `_`, strip leading/trailing `_`,
#'    truncate to 32 characters. This handles the `.` and `:` separators
#'    common in NVR exports (e.g. `NvrEpisode.AdmissionDate`).
#' 2. Explicit remapping of known NVR column names to unambiguous internal
#'    names (e.g. `nvrepisode_admissiondate` → `nvr_admission_date`).
#'    Columns not in the map are kept as-is after step 1.
#'
#' @param nvr_df A tibble of raw NVR data, as loaded from CSV or database.
#' @return The same tibble with standardised column names. All columns are
#'   retained; only names are changed.
#'
#' @details
#' **Columns always renamed** (step 2 remappings currently defined):
#' | Raw (post step 1)               | Internal name           |
#' |----------------------------------|-------------------------|
#' | `nvrepisode_admissiondate`       | `nvr_admission_date`    |
#' | `nvrepisode_proceduredate`       | `nvr_procedure_date`    |
#' | `nvrepisode_surgeoncode`         | `nvr_surgeon_code`      |
#' | `nvrepisode_hospitalcode`        | `nvr_hospital_code`     |
#'
#' `study_id` and `valid_op_date` require no step 2 remapping — they are
#' already correct after generic normalisation.
#'
#' To add a new mapping, append to `.NVR_COLUMN_MAP` in `nvr_utils.R`.
#'
#' @examples
#' \dontrun{
#'   elective_cohort <- read.csv("elective_bypass_participants.csv") %>%
#'     standardise_nvr_names()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_with any_of
#' @keywords internal
standardise_nvr_names <- function(nvr_df) {
  nvr_df %>%
    .standardise_names() %>%
    dplyr::rename(dplyr::any_of(.NVR_COLUMN_MAP))
}

# =============================================================================
# Comorbidity parsing
# =============================================================================

#' Pivot NVR pipe-delimited comorbidity string into binary indicator columns
#'
#' The NVR `RiskScores:Comorbidities` column stores comorbidities as a
#' pipe-separated string of integer codes (e.g. `"1|3|5"`). This function
#' splits each row on `|`, pivots wide, and produces one binary column per
#' unique code value, prefixed with `"comorbidity_"`.
#'
#' Rows with no comorbidities recorded (empty or `NA`) receive `0` across all
#' indicator columns.
#'
#' @param df              A tibble / data frame.
#' @param comorbidity_col Name of the pipe-separated comorbidity column.
#'   Defaults to `"RiskScores:Comorbidities"`.
#' @return Tibble with the comorbidity string column replaced by binary
#'   indicator columns named `comorbidity_<code>` (e.g. `comorbidity_1`,
#'   `comorbidity_2`, ..., `comorbidity_8`).
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom dplyr mutate all_of
#' @keywords internal
separate_out_nvr_comorbidities <- function(df,
                                            comorbidity_col = "RiskScores:Comorbidities") {
  df %>%
    tidyr::separate_rows(dplyr::all_of(comorbidity_col), sep = "\\|") %>%
    dplyr::mutate(value = 1L) %>%
    tidyr::pivot_wider(
      names_from   = dplyr::all_of(comorbidity_col),
      values_from  = value,
      names_prefix = "comorbidity_",
      values_fill  = list(value = 0L)
    )
}