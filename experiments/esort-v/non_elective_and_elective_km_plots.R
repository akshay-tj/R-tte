# =============================================================================
# KM curves — Non-elective vs elective
# Mortality (overall survival) and amputation-free survival (AFS)
#
# Assumes the following objects are already in the environment:
#
# From non_elective_create_analysis_df.R:
#   - non_elective_cohort   : tibble with columns STUDY_ID,
#                             NvrEpisode:ProcedureStartDate, early_surgery
#   - non_elective_outcomes : tibble with columns study_id, afs_days, afs_event
#
# From elective_create_analysis_df.R:
#   - elective_cohort       : tibble with columns STUDY_ID,
#                             NvrEpisode:ProcedureStartDate, early_surgery
#   - elective_outcomes     : tibble with columns study_id, afs_days, afs_event
#
# Shared:
#   - mortality_clean       : tibble with columns study_id, death_date
#
# Output saved to disk:
#   - combined_km_plot_path
# =============================================================================

# load packages ----------------------------------------------------------------
library(survival)
library(ggsurvfit)
library(dplyr)
library(lubridate)
library(patchwork)

# parameters -------------------------------------------------------------------
end_of_study <- as.Date("2024-03-31")

combined_km_plot_path <-
  "Z:/PHP/HSR/ESORT-V/ESORT-V/km_non_elective_vs_elective_combined.png"

# common plot settings ---------------------------------------------------------
x_axis_km <- scale_x_continuous(
  limits = c(0, 365),
  breaks = c(0, 90, 180, 365),
  expand = expansion(mult = c(0, 0.02))
)

y_axis_km <- scale_y_continuous(
  limits = c(0.50, 1.00),
  breaks = c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00),
  expand = expansion(mult = c(0, 0.02)),
  labels = scales::percent_format(accuracy = 1)
)

group_colours <- list(
  scale_color_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "Later surgery", "1" = "Early surgery")
  ),
  scale_fill_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "Later surgery", "1" = "Early surgery")
  )
)

km_theme <- theme_classic() +
  theme(
    legend.position = "top",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(size = 12, face = "bold", hjust = 0)
  )

# helper functions -------------------------------------------------------------
make_mortality_df <- function(cohort_df) {
  cohort_df %>%
    select(STUDY_ID, `NvrEpisode:ProcedureStartDate`, early_surgery) %>%
    mutate(study_id = as.character(STUDY_ID)) %>%
    left_join(mortality_clean, by = "study_id") %>%
    rename(follow_up_start_date = `NvrEpisode:ProcedureStartDate`) %>%
    mutate(
      follow_up_start_date = as.Date(follow_up_start_date),
      death_date           = as.Date(death_date),
      exit_date            = if_else(!is.na(death_date), death_date, end_of_study),
      time_to_event        = as.numeric(exit_date - follow_up_start_date),
      died                 = as.integer(!is.na(death_date))
    ) %>%
    filter(time_to_event >= 0)
}

make_afs_df <- function(cohort_df, outcomes_df) {
  cohort_df %>%
    select(STUDY_ID, early_surgery) %>%
    mutate(study_id = as.character(STUDY_ID)) %>%
    left_join(
      outcomes_df %>% select(study_id, afs_days, afs_event),
      by = "study_id"
    ) %>%
    filter(afs_days >= 0)
}

make_km_plot <- function(data, time_var, event_var, title_text,
                         y_label, x_label = "Days from surgery") {

  surv_formula <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ early_surgery")
  )

  survfit2(surv_formula, data = data) %>%
    ggsurvfit() +
    add_confidence_interval() +
    x_axis_km +
    y_axis_km +
    group_colours +
    labs(
      title = title_text,
      x     = x_label,
      y     = y_label,
      color = NULL,
      fill  = NULL
    ) +
    km_theme
}

# prepare data -----------------------------------------------------------------
df_mortality_non_elective <- make_mortality_df(non_elective_cohort)
df_afs_non_elective       <- make_afs_df(non_elective_cohort, non_elective_outcomes)

df_mortality_elective <- make_mortality_df(elective_cohort)
df_afs_elective       <- make_afs_df(elective_cohort, elective_outcomes)

# create plots -----------------------------------------------------------------
km_mortality_non_elective <- make_km_plot(
  data       = df_mortality_non_elective,
  time_var   = "time_to_event",
  event_var  = "died",
  title_text = "A. Overall survival: Non-elective population",
  y_label    = "Survival probability",
  x_label    = NULL
) +
  theme(axis.title.x = element_blank())

km_mortality_elective <- make_km_plot(
  data       = df_mortality_elective,
  time_var   = "time_to_event",
  event_var  = "died",
  title_text = "C. Overall survival: Elective population",
  y_label    = "Survival probability",
  x_label    = NULL
) +
  theme(axis.title.x = element_blank())

km_afs_non_elective <- make_km_plot(
  data       = df_afs_non_elective,
  time_var   = "afs_days",
  event_var  = "afs_event",
  title_text = "B. Amputation-free survival: Non-elective population",
  y_label    = "Amputation-free survival",
  x_label    = "Days from surgery"
)

km_afs_elective <- make_km_plot(
  data       = df_afs_elective,
  time_var   = "afs_days",
  event_var  = "afs_event",
  title_text = "D. Amputation-free survival: Elective population",
  y_label    = "Amputation-free survival",
  x_label    = "Days from surgery"
)

# combine plots ----------------------------------------------------------------
combined_km_plot <-
  (km_mortality_non_elective | km_mortality_elective) /
  (km_afs_non_elective | km_afs_elective) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top"
  )

# save output ------------------------------------------------------------------
ggsave(
  filename = combined_km_plot_path,
  plot     = combined_km_plot,
  width    = 14,
  height   = 10,
  dpi      = 300
)