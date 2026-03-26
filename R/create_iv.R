# TODO: column names here use raw NVR naming convention (pre-standardisation).
#       When this function moves into iv.R, update to standardised names:
#       NvrHospitalName          -> nvr_hospital_name
#       Patient.PatientId        -> study_id
#       NvrEpisode.AdmissionDate -> nvr_admission_date
#       early_surgery remains the same

#' Calculate Tendency to Expedite Surgery (TTES) for a single patient
#'
#' Computes the proportion of other patients at the same hospital who underwent
#' early surgery in the \code{lookback_days} prior to the index patient's time
#' zero date. Intended to be called row-wise via \code{dplyr::rowwise()} in an
#' experiment script — not vectorised.
#'
#' @param nvr_df A tibble of NVR data containing all candidate peer patients.
#'   Must include \code{NvrHospitalName}, \code{Patient.PatientId},
#'   \code{NvrEpisode.AdmissionDate}, and \code{early_surgery}.
#' @param patient_id The ID of the index patient to exclude from calculation.
#' @param hospital_name The hospital to filter on.
#' @param time_zero_date The index date for the patient. Only procedures
#'   strictly before this date and within \code{lookback_days} prior are
#'   included.
#' @param lookback_days Integer. Width of the lookback window in days.
#'   Defaults to \code{365L}.
#'
#' @return A single numeric value (proportion of peers with early surgery),
#'   or \code{NA_real_} if no eligible peers exist in the window.
#'
#' @examples
#' \dontrun{
#'   calc_ttes(nvr_df, patient_id = 123, hospital_name = "Royal Infirmary",
#'             time_zero_date = as.Date("2020-06-01"))
#'
#'   # Shorter lookback window
#'   calc_ttes(nvr_df, patient_id = 123, hospital_name = "Royal Infirmary",
#'             time_zero_date = as.Date("2020-06-01"), lookback_days = 180L)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter summarise pull
calc_ttes <- function(nvr_df,
                      patient_id,
                      hospital_name,
                      time_zero_date,
                      lookback_days = 365L) {

  # Step 1: Get a list of all patient ids for the same hospital before time_zero_date but within the lookback window for that patient

  df_filtered <- nvr_df %>%
    dplyr::filter(
      NvrHospitalName          == hospital_name,
      Patient.PatientId        != patient_id, # Exclude the given patient_id
      NvrEpisode.AdmissionDate <  time_zero_date,
      NvrEpisode.AdmissionDate >= time_zero_date - lookback_days
    )
  # Step 2: Get the sum of urgent_surgery for the ids obtained, divided by the number of IDs
  
  # Return NA if no patients in window
  if (nrow(df_filtered) == 0L) return(NA_real_)

  # If non missing, calculate the proportion of patients with early surgery
  df_filtered %>%
    dplyr::summarise(ttes = mean(early_surgery, na.rm = TRUE)) %>%
    dplyr::pull(ttes)
}
