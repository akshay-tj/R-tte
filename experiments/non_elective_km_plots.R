# =============================================================================
# KM curves — Mortality and Amputation-free survival (AFS)
#
# Assumes the following objects are already in the environment,
# produced by non_elective_create_analysis_df.R:
#   - non_elective_cohort   : cohort tibble with columns STUDY_ID,
#                             NvrEpisode.ProcedureStartDate, early_surgery
#   - mortality_clean       : tibble with columns study_id, death_date
#   - non_elective_outcomes : tibble with columns study_id, afs_days,
#                             afs_event (joined from compute_limb_outcomes())
#
# Outputs (saved to disk):
#   - mortality_km_plot_path
#   - afs_km_plot_path
# =============================================================================

# load packages
library(survival)
library(ggsurvfit)
library(dplyr)
library(lubridate)

# parameters
end_of_study           <- as.Date("2024-03-31")
mortality_km_plot_path <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/mortality_km_plot.png"
afs_km_plot_path       <- "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/afs_km_plot.png"

# =============================================================================
# KM curve — Mortality
# =============================================================================
df_mortality <- non_elective_cohort %>%
  select(STUDY_ID, `NvrEpisode:ProcedureStartDate`, early_surgery) %>%
  mutate(study_id = as.character(STUDY_ID)) %>%
  left_join(mortality_clean, by = "study_id") %>%
  rename(follow_up_start_date = `NvrEpisode:ProcedureStartDate`) %>%
  mutate(
    follow_up_start_date = as.Date(follow_up_start_date),
    exit_date            = if_else(!is.na(death_date), death_date, end_of_study),
    time_to_event        = as.numeric(exit_date - follow_up_start_date),
    died                 = as.integer(!is.na(death_date))
  ) %>%
  filter(time_to_event >= 0)

km_mortality <- survfit2(Surv(time_to_event, died) ~ early_surgery,
                         data = df_mortality) %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  scale_x_continuous(
    limits = c(0, 365),
    breaks = c(0, 90, 180, 365)
  ) +
  scale_y_continuous(
    limits = c(0.65, 1.0),
    breaks = c(0.65, 0.7, 0.8, 0.9, 1.0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "Late surgery", "1" = "Early surgery")
  ) +
  scale_fill_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "Late surgery", "1" = "Early surgery")
  ) +
  add_pvalue(location = "annotation") +
  labs(
    x     = "Days from surgery",
    y     = "Survival probability",
    color = NULL,
    fill  = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12)
  )

ggsave(mortality_km_plot_path, km_mortality, width = 8, height = 6, dpi = 300)

# =============================================================================
# KM curve — Amputation-free survival (AFS)
# =============================================================================
df_afs <- non_elective_cohort %>%
  select(STUDY_ID, early_surgery) %>%
  mutate(study_id = as.character(STUDY_ID)) %>%
  left_join(
    non_elective_outcomes %>% select(study_id, afs_days, afs_event),
    by = "study_id"
  ) %>%
  filter(afs_days >= 0)

km_afs <- survfit2(Surv(afs_days, afs_event) ~ early_surgery,
                   data = df_afs) %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  scale_x_continuous(
    limits = c(0, 365),
    breaks = c(0, 90, 180, 365)
  ) +
  scale_y_continuous(
    limits = c(0.65, 1.0),
    breaks = c(0.65, 0.7, 0.8, 0.9, 1.0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "Late surgery", "1" = "Early surgery")
  ) +
  scale_fill_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "Late surgery", "1" = "Early surgery")
  ) +
  add_pvalue(location = "annotation") +
  labs(
    x     = "Days from surgery",
    y     = "Amputation-free survival",
    color = NULL,
    fill  = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12)
  )

ggsave(afs_km_plot_path, km_afs, width = 8, height = 6, dpi = 300)