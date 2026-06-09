# =============================================================================
# Two-panel forest plot:
# Non-elective population vs elective population
#
# Plot shows mean difference and 95% CI for:
#   1. Main analysis
#   2. Alternate analysis 1: DAOH measure scores deaths as zero
#   3. Alternate analysis 2: assumes no unmeasured confounding
#   4. Alternate analysis 3: shorter early-surgery window, elective only
# =============================================================================

library(tidyverse)
library(patchwork)
library(grid)

# -----------------------------------------------------------------------------
# 1. ENTER YOUR RESULTS HERE
# -----------------------------------------------------------------------------
# Fill in estimate, lower, and upper.
# For analyses that do not apply, leave estimate/lower/upper as NA and use
# custom_text = "Not applicable".

forest_df <- tribble(
  ~cohort,                  ~analysis,                                                              ~estimate, ~lower, ~upper, ~custom_text,
  "Non-elective population", "Main analysis",                                                         0.075,        -7.847,     7.997,     NA,
  "Non-elective population", "Alternate analysis 1\nDeaths scored as zero",                           -0.639,        -8.449,     7.171,     NA,
  "Non-elective population", "Alternate analysis 2\nAssumes no unmeasured confounding",               3.079,        1.898,     4.260,     NA,
  "Non-elective population", "Alternate analysis 3\nShorter early-surgery window (15–28 days)",      NA,        NA,     NA,     "Not applicable",

  "Elective population",     "Main analysis",                                                         -4.415,       -11.017,     2.188,     NA,
  "Elective population",     "Alternate analysis 1\nDeaths scored as zero",                           -3.060,        -10.728,     4.608,     NA,
  "Elective population",     "Alternate analysis 2\nAssumes no unmeasured confounding",               0.251,        -1.068,     1.569,     NA,
  "Elective population",     "Alternate analysis 3\nShorter early-surgery window (15–28 days)",      -4.342,        -11.112,     2.427,     NA
)

# -----------------------------------------------------------------------------
# 2. SETTINGS
# -----------------------------------------------------------------------------

plot_title <- "DAOH at 90 days"

output_path <- "Z:/PHP/HSR/ESORT-V/ESORT-V/forest_plot_non_elective_vs_elective.png"

analysis_levels <- c(
  "Main analysis",
  "Alternate analysis 1\nDeaths scored as zero",
  "Alternate analysis 2\nAssumes no unmeasured confounding",
  "Alternate analysis 3\nShorter early-surgery window (15–28 days)"
)

cohort_levels <- c("Non-elective population", "Elective population")

digits_est <- 3

# -----------------------------------------------------------------------------
# 3. PREPARE DATA
# -----------------------------------------------------------------------------

fmt_num <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

forest_df <- forest_df %>%
  mutate(
    cohort   = factor(cohort, levels = cohort_levels),
    analysis = factor(analysis, levels = analysis_levels)
  ) %>%
  arrange(cohort, analysis) %>%
  group_by(cohort) %>%
  mutate(
    y = rev(seq_along(analysis_levels))
  ) %>%
  ungroup() %>%
  mutate(
    md_ci_text = case_when(
      !is.na(custom_text) ~ custom_text,
      !is.na(estimate) & !is.na(lower) & !is.na(upper) ~
        paste0(
          fmt_num(estimate, digits_est), " (",
          fmt_num(lower, digits_est), ", ",
          fmt_num(upper, digits_est), ")"
        ),
      TRUE ~ ""
    )
  )

# -----------------------------------------------------------------------------
# 4. X-AXIS RANGE AND TEXT POSITIONS
# -----------------------------------------------------------------------------

max_abs_ci <- max(abs(c(forest_df$lower, forest_df$upper)), na.rm = TRUE)

if (!is.finite(max_abs_ci)) max_abs_ci <- 10

x_plot_lim <- max_abs_ci * 1.25
text_x     <- x_plot_lim * 1.55
x_rightlim <- x_plot_lim * 2.25

x_breaks <- pretty(c(-x_plot_lim, x_plot_lim), n = 5)

# -----------------------------------------------------------------------------
# 5. HELPER FUNCTION FOR ONE FOREST PLOT
# -----------------------------------------------------------------------------

make_single_forest_plot <- function(plot_df, cohort_title, show_y_axis = TRUE, plot_margin = margin(20, 40, 20, 20)) {
  
  header_df <- tibble(
    x     = text_x,
    y     = 4.7,
    label = "Mean difference (95% CI)"
  )
  
  # Move labels close to the zero line.
  # If still too far apart, reduce 0.12 to 0.08.
  ann_df <- tibble(
    x     = c(-x_plot_lim * 0.12, x_plot_lim * 0.12),
    y     = c(4.9, 4.9),
    label = c("Favours later surgery", "Favours early surgery"),
    hjust = c(1, 0)
  )
  
  y_scale <- if (show_y_axis) {
    scale_y_continuous(
      breaks = c(4, 3, 2, 1),
      labels = analysis_levels,
      expand = expansion(mult = c(0.06, 0.18))
    )
  } else {
    scale_y_continuous(
      breaks = c(4, 3, 2, 1),
      labels = rep("", 4),
      expand = expansion(mult = c(0.06, 0.18))
    )
  }
  
  ggplot(plot_df, aes(y = y)) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linewidth = 0.7
    ) +
    geom_segment(
      data = plot_df %>%
        filter(!is.na(estimate), !is.na(lower), !is.na(upper)),
      aes(x = lower, xend = upper, y = y, yend = y),
      linewidth = 0.8
    ) +
    geom_segment(
      data = plot_df %>%
        filter(!is.na(estimate), !is.na(lower), !is.na(upper)),
      aes(x = lower, xend = lower, y = y - 0.10, yend = y + 0.10),
      linewidth = 0.8
    ) +
    geom_segment(
      data = plot_df %>%
        filter(!is.na(estimate), !is.na(lower), !is.na(upper)),
      aes(x = upper, xend = upper, y = y - 0.10, yend = y + 0.10),
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_df %>% filter(!is.na(estimate)),
      aes(x = estimate),
      shape = 18,
      size = 3
    ) +
    geom_text(
      aes(x = text_x, label = md_ci_text),
      hjust = 0,
      size = 3.5
    ) +
    geom_text(
      data = header_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.7
    ) +
    geom_text(
      data = ann_df,
      aes(x = x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      fontface = "italic",
      size = 3.5
    ) +
    y_scale +
    scale_x_continuous(
      breaks = x_breaks,
      limits = c(-x_plot_lim, x_rightlim)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = cohort_title,
      x = "Days",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.y  = element_text(hjust = 0),
      axis.ticks.y = element_blank(),
      plot.margin  = plot_margin
    )
}

# -----------------------------------------------------------------------------
# 6. SPLIT DATA
# -----------------------------------------------------------------------------

forest_non_elective <- forest_df %>%
  filter(cohort == "Non-elective population")

forest_elective <- forest_df %>%
  filter(cohort == "Elective population")

# -----------------------------------------------------------------------------
# 7. CREATE PLOTS
# -----------------------------------------------------------------------------

plot_non_elective <- make_single_forest_plot(
  plot_df      = forest_non_elective,
  cohort_title = "Non-elective population",
  show_y_axis  = TRUE,
  plot_margin  = margin(20, 40, 20, 20)
)

plot_elective <- make_single_forest_plot(
  plot_df      = forest_elective,
  cohort_title = "Elective population",
  show_y_axis  = FALSE,
  plot_margin  = margin(20, 80, 20, 20) 
)

# -----------------------------------------------------------------------------
# 8. COMBINE PLOTS
# -----------------------------------------------------------------------------
# Increase the middle value if you want a larger gap between panels.

combined_forest_plot <-
  plot_non_elective + plot_spacer() + plot_elective +
  plot_layout(widths = c(1, 0.05, 1)) +
  plot_annotation(title = plot_title) &
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )

print(combined_forest_plot)

# -----------------------------------------------------------------------------
# 9. SAVE
# -----------------------------------------------------------------------------

ggsave(
  filename = output_path,
  plot     = combined_forest_plot,
  width    = 16,
  height   = 7,
  dpi      = 300
)

