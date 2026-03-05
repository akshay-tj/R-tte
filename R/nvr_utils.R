# nvr_utils.R
#
# NVR-specific column cleaning helpers (internal).
# All functions here are called at data ingestion — never mid-pipeline.
#
# Only standardise_nvr_names() is currently defined. Additional NVR-specific
# helpers (e.g. procedure code lookups, validity flags) go here as needed.

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
