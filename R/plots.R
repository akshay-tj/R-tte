# plots.R
# Visualisation functions for IV diagnostics and covariate balance.
#
# Functions:
#   plot_covariate_balance()  — standardised covariate balance across IV deciles
#   plot_iv_stability()       — IV stability across centres (overall + facet)

# ==============================================================================
# plot_covariate_balance
# ==============================================================================

#' Plot covariate balance across instrumental variable deciles
#'
#' For each covariate, computes the mean of the SD-rescaled value within each
#' IV decile and plots all covariates on a single panel. Categorical covariates
#' are one-hot encoded internally before rescaling.
#'
#' @param df A tibble containing the IV and all covariates of interest.
#' @param iv_col String. Name of the instrumental variable column.
#' @param covariates Character vector. Names of covariate columns to include.
#'   All columns listed here will appear in the plot.
#' @param categorical_cols Character vector. Subset of \code{covariates} that
#'   are categorical and should be one-hot encoded before rescaling. All
#'   dummy levels are retained (\code{remove_first_dummy = FALSE}).
#'   Pass \code{NULL} (default) if all covariates are already numeric.
#' @param labels Named character vector mapping internal column names (after
#'   one-hot encoding) to display labels. Names are the column names as they
#'   appear post-encoding (e.g. \code{"Gender_Male"}); values are the labels
#'   shown in the legend. Columns not present in \code{labels} are displayed
#'   using their raw column name.
#' @param x_label String. X-axis label. Defaults to
#'   \code{"Instrumental Variable (Decile)"}.
#' @param y_label String. Y-axis label. Defaults to a standard description of
#'   the rescaling applied.
#' @param title String. Plot title. Defaults to
#'   \code{"Covariate Balance Plot"}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Rescaling divides each covariate by its overall SD (no centring), following
#' the approach in Brookhart et al. (2006) for covariate balance assessment.
#' One-hot encoding uses \code{fastDummies::dummy_cols()} with
#' \code{remove_first_dummy = FALSE} so all levels are visible in the plot.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate across group_by summarise all_of ntile
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line labs scale_x_continuous
#'   scale_y_continuous scale_color_discrete theme_minimal theme
#'   element_line element_blank
#' @export
plot_covariate_balance <- function(df,
                                   iv_col,
                                   covariates,
                                   categorical_cols = NULL,
                                   labels           = NULL,
                                   x_label          = "Instrumental Variable (Decile)",
                                   y_label          = "Mean baseline covariate rescaled by its SD",
                                   title            = "Covariate Balance Plot") {

  # --- Input checks -----------------------------------------------------------
  stopifnot(
    is.data.frame(df),
    is.character(iv_col),        length(iv_col)   == 1, iv_col   %in% names(df),
    is.character(covariates),    length(covariates) > 0,
    all(covariates %in% names(df))
  )
  if (!is.null(categorical_cols)) {
    stopifnot(
      is.character(categorical_cols),
      all(categorical_cols %in% covariates)
    )
  }

  # --- Step 1: Subset to IV + covariates --------------------------------------
  working <- df %>%
    dplyr::select(dplyr::all_of(c(iv_col, covariates)))

  # --- Step 2: One-hot encode categorical columns (if any) --------------------
  if (!is.null(categorical_cols) && length(categorical_cols) > 0) {
    if (!requireNamespace("fastDummies", quietly = TRUE)) {
      stop(
        "Package 'fastDummies' is required for one-hot encoding.\n",
        "Install with: install.packages('fastDummies')",
        call. = FALSE
      )
    }
    working <- fastDummies::dummy_cols(
      working,
      select_columns          = categorical_cols,
      remove_first_dummy      = FALSE,
      remove_selected_columns = TRUE
    )
  }

  # --- Step 3: Identify numeric covariate columns after encoding --------------
  covariate_cols <- setdiff(names(working), iv_col)

  # --- Step 4: Rescale each covariate by its overall SD (no centring) ---------
  rescaled_cols <- paste0("sd_", covariate_cols)

  working <- working %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(covariate_cols),
        ~ .x / sd(.x, na.rm = TRUE),
        .names = "sd_{.col}"
      )
    )

  # --- Step 5: Bin IV into deciles --------------------------------------------
  working <- working %>%
    dplyr::mutate(iv_decile = dplyr::ntile(.data[[iv_col]], 10))

  # --- Step 6: Mean rescaled covariate per decile -----------------------------
  decile_means <- working %>%
    dplyr::group_by(iv_decile) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(rescaled_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )

  # --- Step 7: Reshape to long format -----------------------------------------
  plot_data <- decile_means %>%
    tidyr::pivot_longer(
      cols      = dplyr::starts_with("mean_sd_"),
      names_to  = "covariate",
      values_to = "mean_value"
    ) %>%
    dplyr::mutate(covariate = gsub("mean_sd_", "", covariate))

  # --- Step 8: Build plot -----------------------------------------------------
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x     = iv_decile,
      y     = mean_value,
      color = covariate,
      group = covariate
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9, linetype = "dashed") +
    ggplot2::labs(
      title = title,
      x     = x_label,
      y     = y_label
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 10, by = 2)) +
    ggplot2::scale_y_continuous(limits = c(NA, NA), expand = c(0, 0.05)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line    = ggplot2::element_line(color = "black"),
      panel.border = ggplot2::element_blank()
    )

  # --- Step 9: Apply labels if supplied ---------------------------------------
  if (!is.null(labels)) {
    # For unlabelled covariates, fall back to the raw column name
    all_covariates <- unique(plot_data$covariate)
    full_labels    <- stats::setNames(all_covariates, all_covariates)  # identity fallback
    full_labels[names(labels)] <- labels                               # override with supplied

    p <- p +
      ggplot2::scale_color_discrete(
        name   = "Baseline Covariates",
        labels = full_labels,
        breaks = names(full_labels)
      )
  } else {
    p <- p +
      ggplot2::scale_color_discrete(name = "Baseline Covariates")
  }

  p
}


# ==============================================================================
# plot_iv_stability
# ==============================================================================

#' Plot instrumental variable stability across centres
#'
#' Returns two plots as a named list:
#' \describe{
#'   \item{\code{overall}}{A single line showing the mean IV value per centre,
#'     ordered by row index (anonymised). Centre names are not displayed.}
#'   \item{\code{facet}}{One panel per centre showing the mean IV over time,
#'     faceted by centre. Centre names and x-axis tick labels are suppressed
#'     for anonymity.}
#' }
#'
#' @param df A tibble. Must already be filtered to the desired subset
#'   (e.g. include flags, year range, minimum data requirements) before
#'   passing to this function.
#' @param iv_col String. Name of the instrumental variable column.
#' @param centre_col String. Name of the centre/hospital identifier column.
#' @param date_col String. Name of the date column used to extract the year
#'   for the facet plot. \strong{Must be a \code{Date} column} — character
#'   date columns should be converted upstream with \code{as.Date()} before
#'   calling this function.
#' @param x_label_overall String. X-axis label for the overall plot.
#'   Defaults to \code{"Hospital (anonymised index)"}.
#' @param y_label String. Y-axis label shared by both plots.
#'   Defaults to \code{"Mean IV value"}.
#'
#' @return A named list with elements \code{overall} and \code{facet},
#'   each a \code{ggplot} object.
#'
#' @details
#' The overall plot assigns each centre an anonymous integer index based on
#' its row order after grouping — no centre names are shown.
#'
#' The facet plot suppresses both x-axis tick labels and facet strip text,
#' so individual centres cannot be identified from the output.
#'
#' Year extraction uses \code{format(date, "\%Y")} which requires a
#' \code{Date} input. Passing a character column will produce incorrect
#' results without an error — convert upstream.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate row_number n
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#'   theme facet_wrap scale_y_continuous element_blank element_text
#' @importFrom scales number_format
#' @export
plot_iv_stability <- function(df,
                              iv_col,
                              centre_col,
                              date_col,
                              x_label_overall = "Hospital (anonymised index)",
                              y_label         = "Mean IV value") {

  # --- Input checks -----------------------------------------------------------
  stopifnot(
    is.data.frame(df),
    is.character(iv_col),     length(iv_col)     == 1, iv_col     %in% names(df),
    is.character(centre_col), length(centre_col) == 1, centre_col %in% names(df),
    is.character(date_col),   length(date_col)   == 1, date_col   %in% names(df)
  )
  if (!inherits(df[[date_col]], "Date")) {
    warning(
      sprintf(
        "`%s` is not a Date column (class: %s). Year extraction may be incorrect. ",
        date_col, paste(class(df[[date_col]]), collapse = "/")
      ),
      "Convert upstream with as.Date() before calling plot_iv_stability().",
      call. = FALSE
    )
  }

  # --- Overall plot -----------------------------------------------------------
  overall_summary <- df %>%
    dplyr::group_by(.data[[centre_col]]) %>%
    dplyr::summarise(
      mean_iv = mean(.data[[iv_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(hosp_idx = dplyr::row_number())

  overall_plot <- ggplot2::ggplot(
    overall_summary,
    ggplot2::aes(x = hosp_idx, y = mean_iv)
  ) +
    ggplot2::geom_line(color = "darkblue", linewidth = 0.75) +
    ggplot2::labs(
      x = x_label_overall,
      y = y_label
    ) +
    ggplot2::theme_minimal()

  # --- Facet plot -------------------------------------------------------------
  facet_summary <- df %>%
    dplyr::mutate(year = format(.data[[date_col]], "%Y")) %>%
    dplyr::group_by(.data[[centre_col]], year) %>%
    dplyr::summarise(
      mean_iv = mean(.data[[iv_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # Y-axis breaks: min, mean, max of the summarised values
  y_breaks <- c(
    min(facet_summary$mean_iv,  na.rm = TRUE),
    mean(facet_summary$mean_iv, na.rm = TRUE),
    max(facet_summary$mean_iv,  na.rm = TRUE)
  )

  facet_plot <- ggplot2::ggplot(
    facet_summary,
    ggplot2::aes(
      x     = year,
      y     = mean_iv,
      group = .data[[centre_col]]
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = NULL,
      y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", centre_col))
    ) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12),
      strip.text   = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = scales::number_format()
    )

  list(
    overall = overall_plot,
    facet   = facet_plot
  )
}