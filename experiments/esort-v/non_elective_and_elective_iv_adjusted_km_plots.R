# =============================================================================
# IV-adjusted predicted AFS curves — Non-elective vs elective
#
# Requires model bundles created after running Steps 1–9 for each cohort:
#   - iv_adjusted_afs_non_elective_model_bundle.rds
#   - iv_adjusted_afs_elective_model_bundle.rds
#
# Output:
#   - IV-adjusted predicted AFS plot, non-elective and elective side-by-side
# =============================================================================


# load packages ----------------------------------------------------------------

library(survival)
library(tidyverse)
library(patchwork)
library(scales)


# paths ------------------------------------------------------------------------

non_elective_bundle_path <-
  "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/iv_adjusted_afs_non_elective_model_bundle.rds"

elective_bundle_path <-
  "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/iv_adjusted_afs_elective_model_bundle.rds"

iv_adjusted_afs_plot_path <-
  "Z:/PHP/HSR/ESORT-V/ESORT-V/iv_adjusted_afs_non_elective_vs_elective.png"


# common plot settings ---------------------------------------------------------

plot_times <- seq(0, 365, by = 5)

x_axis_iv_afs <- scale_x_continuous(
  limits = c(0, 365),
  breaks = c(0, 90, 180, 365),
  expand = expansion(mult = c(0, 0.02))
)

y_axis_iv_afs <- scale_y_continuous(
  limits = c(0.50, 1.00),
  breaks = c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00),
  expand = expansion(mult = c(0, 0.02)),
  labels = percent_format(accuracy = 1)
)

iv_afs_colours <- list(
  scale_color_manual(
    values = c(
      "Later surgery" = "steelblue",
      "Early surgery" = "firebrick"
    )
  ),
  scale_fill_manual(
    values = c(
      "Later surgery" = "steelblue",
      "Early surgery" = "firebrick"
    )
  )
)

iv_afs_theme <- theme_classic() +
  theme(
    legend.position = "top",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(size = 12, face = "bold", hjust = 0)
  )


# helper functions -------------------------------------------------------------

make_preddata <- function(data, d_value, lasso_obj) {
  
  pred <- data
  
  # Set treatment strategy counterfactually
  pred$early_surgery <- d_value
  
  # Identify true D*X terms, excluding internal residual placeholders
  dx_terms <- grep("_d$", lasso_obj$final_outcome_terms, value = TRUE)
  dx_terms <- setdiff(dx_terms, c(".ctrl_resid_x_d"))
  dx_terms <- intersect(dx_terms, names(pred))
  
  # Recalculate D*X interactions
  for (term in dx_terms) {
    
    base_var <- sub("_d$", "", term)
    
    if (!base_var %in% names(pred)) {
      stop(sprintf(
        "Cannot create '%s': base variable '%s' not found in prediction data.",
        term, base_var
      ))
    }
    
    pred[[term]] <- pred$early_surgery * pred[[base_var]]
  }
  
  # Recalculate residual-by-treatment interaction
  if (!"generalised_residual" %in% names(pred)) {
    stop("generalised_residual not found in prediction data.")
  }
  
  pred$generalised_residual_x_d <-
    pred$generalised_residual * pred$early_surgery
  
  pred
}


extract_hr <- function(model, term = "early_surgery") {
  
  c_idx <- which(names(coef(model)) == term)
  
  if (length(c_idx) != 1) {
    stop(sprintf("Term '%s' not found uniquely in model coefficients.", term))
  }
  
  beta <- coef(model)[c_idx]
  se   <- sqrt(vcov(model)[c_idx, c_idx])
  
  tibble(
    term = term,
    beta = beta,
    se   = se,
    HR   = exp(beta),
    LCI  = exp(beta - 1.96 * se),
    UCI  = exp(beta + 1.96 * se),
    p    = 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  )
}


average_survfit_at_times <- function(sf, arm_label, times = plot_times) {
  
  sf_sum <- summary(sf, times = times, extend = TRUE)
  
  n_times <- length(times)
  n_vals  <- length(sf_sum$surv)
  
  if (n_vals %% n_times != 0) {
    stop("Unexpected survfit summary structure: cannot reshape survival output.")
  }
  
  n_curves <- n_vals / n_times
  
  surv_mat <- matrix(
    sf_sum$surv,
    nrow = n_times,
    ncol = n_curves,
    byrow = FALSE
  )
  
  lower_mat <- matrix(
    sf_sum$lower,
    nrow = n_times,
    ncol = n_curves,
    byrow = FALSE
  )
  
  upper_mat <- matrix(
    sf_sum$upper,
    nrow = n_times,
    ncol = n_curves,
    byrow = FALSE
  )
  
  tibble(
    time  = times,
    surv  = rowMeans(surv_mat,  na.rm = TRUE),
    lower = rowMeans(lower_mat, na.rm = TRUE),
    upper = rowMeans(upper_mat, na.rm = TRUE),
    arm   = factor(
      arm_label,
      levels = c("Later surgery", "Early surgery")
    )
  )
}


make_iv_adjusted_afs_curve <- function(bundle) {
  
  afs_df          <- bundle$afs_df
  lasso_afs       <- bundle$lasso_afs
  cox_formula_str <- bundle$cox_formula_str
  cox_afs         <- bundle$cox_afs
  
  # Prediction-only model without frailty.
  # survfit.coxph() cannot use newdata with frailty terms.
  cox_afs_pred <- coxph(
    as.formula(cox_formula_str),
    data  = afs_df,
    x     = TRUE,
    model = TRUE
  )
  
  # Print HR comparison
  hr_compare <- bind_rows(
    extract_hr(cox_afs, "early_surgery") %>%
      mutate(model = "Cox 2SRI with patient frailty"),
    extract_hr(cox_afs_pred, "early_surgery") %>%
      mutate(model = "Cox 2SRI without frailty; prediction model")
  ) %>%
    mutate(
      cohort = bundle$cohort_label,
      HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCI, UCI),
      p_fmt = sprintf("%.4f", p)
    ) %>%
    select(cohort, model, HR_CI, p_fmt)
  
  print(hr_compare, n = Inf)
  
  # Counterfactual prediction datasets
  pred_later <- make_preddata(
    data      = afs_df,
    d_value   = 0,
    lasso_obj = lasso_afs
  )
  
  pred_early <- make_preddata(
    data      = afs_df,
    d_value   = 1,
    lasso_obj = lasso_afs
  )
  
  # Predict survival under each strategy
  sf_later <- survfit(
    cox_afs_pred,
    newdata = pred_later,
    se.fit  = TRUE
  )
  
  sf_early <- survfit(
    cox_afs_pred,
    newdata = pred_early,
    se.fit  = TRUE
  )
  
  # Average predicted survival across patients
  curve_later <- average_survfit_at_times(
    sf        = sf_later,
    arm_label = "Later surgery",
    times     = plot_times
  )
  
  curve_early <- average_survfit_at_times(
    sf        = sf_early,
    arm_label = "Early surgery",
    times     = plot_times
  )
  
  bind_rows(curve_later, curve_early) %>%
    mutate(
      lower  = pmax(lower, 0),
      upper  = pmin(upper, 1),
      cohort = bundle$cohort_label
    )
}


make_iv_afs_plot <- function(curve_df, title_text, show_y_title = TRUE) {
  
  y_label <- if (show_y_title) {
    "Predicted amputation-free survival"
  } else {
    NULL
  }
  
  ggplot(
    curve_df,
    aes(x = time, y = surv, colour = arm, fill = arm)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.18,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    x_axis_iv_afs +
    y_axis_iv_afs +
    iv_afs_colours +
    labs(
      title  = title_text,
      x      = "Days since surgery",
      y      = y_label,
      colour = NULL,
      fill   = NULL
    ) +
    iv_afs_theme
}


# read model bundles -----------------------------------------------------------

non_elective_bundle <- readRDS(non_elective_bundle_path)
elective_bundle     <- readRDS(elective_bundle_path)


# create curve data ------------------------------------------------------------

iv_afs_non_elective <- make_iv_adjusted_afs_curve(non_elective_bundle)
iv_afs_elective     <- make_iv_adjusted_afs_curve(elective_bundle)


# create plots -----------------------------------------------------------------

plot_iv_afs_non_elective <- make_iv_afs_plot(
  curve_df     = iv_afs_non_elective,
  title_text   = "A. Non-elective population",
  show_y_title = TRUE
)

plot_iv_afs_elective <- make_iv_afs_plot(
  curve_df     = iv_afs_elective,
  title_text   = "B. Elective population",
  show_y_title = FALSE
)


# combine plots ----------------------------------------------------------------

combined_iv_afs_plot <-
  plot_iv_afs_non_elective | plot_iv_afs_elective

combined_iv_afs_plot <-
  combined_iv_afs_plot +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top"
  )


# save output ------------------------------------------------------------------

ggsave(
  filename = iv_adjusted_afs_plot_path,
  plot     = combined_iv_afs_plot,
  width    = 14,
  height   = 5.5,
  dpi      = 300
)

print(combined_iv_afs_plot)