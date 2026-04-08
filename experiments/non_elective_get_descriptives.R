# =============================================================================
# SECTION X: Descriptive outcomes tables (Table 2)
# =============================================================================

cont_90 <- c("daoh_bypass_surg_90d", "total_los_no_90d",
             "bypass_surg_proc_los_no", "post_bypass_surg_los_no_90d",
             "bypass_surg_los_no")
cat_90  <- c("readmit_post_bypass_surg_90d", "died_post_bypass_surg_90d")
outcomes_90_days <- table2_outcomes(non_elective_outcomes,
                                    cont_vars = cont_90, cat_vars = cat_90,
                                    horizon_label = "90 days")

cont_180 <- c("daoh_bypass_surg_180d", "total_los_no_180d",
              "bypass_surg_proc_los_no", "post_bypass_surg_los_no_180d",
              "bypass_surg_los_no")
cat_180  <- c("readmit_post_bypass_surg_180d", "died_post_bypass_surg_180d")
outcomes_180_days <- table2_outcomes(non_elective_outcomes,
                                     cont_vars = cont_180, cat_vars = cat_180,
                                     horizon_label = "180 days")

cont_365 <- c("daoh_bypass_surg_365d", "total_los_no_365d",
              "bypass_surg_proc_los_no", "post_bypass_surg_los_no_365d",
              "bypass_surg_los_no")
cat_365  <- c("readmit_post_bypass_surg_365d", "died_post_bypass_surg_365d")
outcomes_365_days <- table2_outcomes(non_elective_outcomes,
                                     cont_vars = cont_365, cat_vars = cat_365,
                                     horizon_label = "365 days")
