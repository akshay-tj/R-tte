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

#' Draw a subgroup forest plot
#'
#' Produces a publication-style forest plot for a single outcome and time
#' horizon. Subgroup rows are organised into labelled group sections defined
#' by \code{subgroup_defs}; group headers are rendered in bold. Point
#' estimates are shown as diamonds with horizontal CI lines and end-caps.
#' Sample sizes and formatted CI strings are printed as text columns outside
#' the plot region.
#'
#' @param plot_data A tibble with one row per subgroup. Required columns:
#'   \code{var} (string matching \code{var} fields in \code{subgroup_defs}),
#'   \code{n} (integer sample size), \code{est} (point estimate),
#'   \code{ci_lo} (lower 95\% CI), \code{ci_hi} (upper 95\% CI).
#' @param outcome_label String. Human-readable outcome name used in the plot
#'   title (e.g. \code{"Days alive and out of hospital"}).
#' @param horizon Integer. Follow-up window in days, appended to the title
#'   (e.g. \code{90} produces \code{"... (90 days)"}).
#' @param x_limits Numeric vector of length 2. Data-range limits for the CI
#'   axis, e.g. \code{c(-10, 10)}. Text columns are positioned relative to
#'   these limits — setting them too narrow will cause label overlap.
#' @param favors_positive Logical. Controls the direction labels printed above
#'   the x-axis. If \code{TRUE}, positive values favour early surgery and
#'   negative values favour later surgery. If \code{FALSE}, the labels are
#'   reversed.
#' @param subgroup_defs A list of named lists defining the row structure of
#'   the plot. Each entry must have three string fields: \code{var} (matches
#'   \code{var} column in \code{plot_data}), \code{label} (display label for
#'   the row), and \code{group} (section header; rows sharing a \code{group}
#'   value are rendered together under a bold header). Subgroups absent from
#'   \code{plot_data} are included as empty rows so layout remains consistent
#'   across outcomes.
#'
#' @return A \code{ggplot} object. Save with \code{ggplot2::ggsave()} or pass
#'   to \code{export_figure_docx()}.
#'
#' @details
#' Text columns are positioned in data units as fixed offsets from
#' \code{x_limits}: the subgroup label column at \code{-35\%}, the N column
#' at \code{-12\%}, and the CI string column at \code{+12\%}. If the plot is
#' saved at a non-standard width these offsets may need adjustment.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate if_else bind_rows
#' @importFrom purrr map_chr
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_vline geom_segment geom_point
#'   geom_text annotate scale_x_continuous coord_cartesian labs
#'   theme_classic theme element_blank element_line element_text margin
#' @export
make_forest_plot <- function(plot_data,
                             outcome_label,
                             horizon,
                             x_limits,
                             favors_positive,
                             subgroup_defs,
                             as_percentage = FALSE) {

  # Resolve direction labels before ggplot call
  label_left  <- if (favors_positive) "Favors later surgery" else "Favors early surgery"
  label_right <- if (favors_positive) "Favors early surgery" else "Favors later surgery"

  # Build ordered row structure with group headers inserted
  # Each group gets a header row (no data) + subgroup rows
  groups <- unique(purrr::map_chr(subgroup_defs, "group"))
  ordered_rows  <- list()
  row_counter   <- 0

  for (grp in groups) {
    # Header row
    row_counter <- row_counter + 1
    ordered_rows[[length(ordered_rows) + 1]] <- tibble::tibble(
      y_pos       = row_counter,
      var         = NA_character_,
      label       = grp,
      n           = NA_integer_,
      est         = NA_real_,
      ci_lo       = NA_real_,
      ci_hi       = NA_real_,
      is_header   = TRUE
    )

    # Subgroup rows for this group
    subs <- subgroup_defs[purrr::map_chr(subgroup_defs, "group") == grp]
    for (s in subs) {
      row_counter <- row_counter + 1
      d <- plot_data %>% filter(var == s$var)
      ordered_rows[[length(ordered_rows) + 1]] <- tibble::tibble(
        y_pos     = row_counter,
        var       = s$var,
        label     = s$label,
        n         = if (nrow(d) > 0) d$n       else NA_integer_,
        est       = if (nrow(d) > 0) d$est     else NA_real_,
        ci_lo     = if (nrow(d) > 0) d$ci_lo   else NA_real_,
        ci_hi     = if (nrow(d) > 0) d$ci_hi   else NA_real_,
        is_header = FALSE
      )
    }
  }

  df <- dplyr::bind_rows(ordered_rows) %>%
    # Reverse so first row plots at top
    mutate(y_pos = max(y_pos) - y_pos + 1)
  
  if (as_percentage) {
  df <- df %>%
    mutate(
      est   = est   * 100,
      ci_lo = ci_lo * 100,
      ci_hi = ci_hi * 100
    )
}

  # CI label for right-hand column
  df <- df %>%
    mutate(
      ci_label = dplyr::if_else(
        !is_header & !is.na(est),
        sprintf("%.3f (%.3f, %.3f)", est, ci_lo, ci_hi),
        NA_character_
      )
    )

  # x positions for text columns (in data units)
  x_left_label  <- x_limits[1] - 0.35 * diff(x_limits)
  x_left_n      <- x_limits[1] - 0.12 * diff(x_limits)
  x_right_ci    <- x_limits[2] + 0.12 * diff(x_limits)

  ggplot(df, aes(y = y_pos)) +

    # Zero line
    geom_vline(xintercept = 0, linetype = "solid", colour = "grey50", linewidth = 0.4) +

    # CI lines (non-header rows only)
    geom_segment(
      data = df %>% filter(!is_header & !is.na(est)),
      aes(x = ci_lo, xend = ci_hi, y = y_pos, yend = y_pos),
      linewidth = 0.5
    ) +

    # Whisker end caps
    geom_segment(
      data = df %>% filter(!is_header & !is.na(est)),
      aes(x = ci_lo, xend = ci_lo,
          y = y_pos - 0.2, yend = y_pos + 0.2),
      linewidth = 0.5
    ) +
    geom_segment(
      data = df %>% filter(!is_header & !is.na(est)),
      aes(x = ci_hi, xend = ci_hi,
          y = y_pos - 0.2, yend = y_pos + 0.2),
      linewidth = 0.5
    ) +

    # Point estimate diamond
    geom_point(
      data = df %>% filter(!is_header & !is.na(est)),
      aes(x = est),
      shape = 18, size = 3
    ) +

    # Subgroup labels (left, indented for non-headers)
    geom_text(
      data = df %>% filter(!is_header),
      aes(x = x_left_label, label = label),
      hjust = 0, size = 3, nudge_x = diff(x_limits) * 0.02
    ) +

    # Group header labels (bold)
    geom_text(
      data = df %>% filter(is_header),
      aes(x = x_left_label, label = label),
      hjust = 0, size = 3, fontface = "bold"
    ) +

    # N column
    geom_text(
      data = df %>% filter(!is_header & !is.na(n)),
      aes(x = x_left_n, label = formatC(n, format = "d", big.mark = ",")),
      hjust = 1, size = 3
    ) +

    # N column header
    annotate(
      "text",
      x     = x_left_n,
      y     = max(df$y_pos) + 1,
      label = "Sample size",
      hjust = 1, size = 3
    ) +

    # CI text column (right)
    geom_text(
      data = df %>% filter(!is.na(ci_label)),
      aes(x = x_right_ci, label = ci_label),
      hjust = 0, size = 3
    ) +

    # CI column header
    annotate(
      "text",
      x     = x_right_ci,
      y     = max(df$y_pos) + 1,
      label = "Mean Difference (95% CI)",
      hjust = 0, size = 3
    ) +
    
    # Direction labels — outcome-specific
    annotate(
      "text",
      x = (x_limits[1] / 2) - diff(x_limits) * 0.05, y = max(df$y_pos) + 2,
      label = label_left,
      hjust = 0, size = 3, fontface = "italic"
    ) +
    annotate(
      "text",
      x = (x_limits[2] / 2) + diff(x_limits) * 0.05, y = max(df$y_pos) + 2,
      label = label_right,
      hjust = 1, size = 3, fontface = "italic"
    ) +

    # Scales and theme
    scale_x_continuous(
      limits = c(x_left_label, x_right_ci + diff(x_limits) * 0.25),
      breaks = pretty(x_limits, n = 5),
      labels = function(x) ifelse(x >= x_limits[1] & x <= x_limits[2],
                                  as.character(x), "")
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = sprintf("%s (%d days)", outcome_label, horizon),
      x = NULL, y = NULL
    ) +
    theme_classic() +
    theme(
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      axis.line.x       = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.x      = element_line(colour = "black"),
      axis.text.x       = element_text(size = 9),
      plot.title        = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.margin       = margin(t = 20, r = 120, b = 10, l = 10)
    )
}
