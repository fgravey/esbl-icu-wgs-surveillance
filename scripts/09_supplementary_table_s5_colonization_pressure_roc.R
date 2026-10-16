# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 09_supplementary_table_s5_colonization_pressure_roc.R
# Purpose: Dynamic computation of ward colonization pressure (CP), stratified
#          ROC curves with bootstrap-validated AUC, optimal Youden inflection
#          thresholds, and step-function odds ratios with exact Fisher
#          fallback (Supplementary Table S5; panels 3A/3B of Figure 3 are
#          also exported from this stratified analysis).
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Figure_3AB_ROC_Colonization_Pressure.pdf (.png)
#          output/tables/Supplementary_Table_S5_Logistic_Models.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_fig <- file.path(base_dir, "output", "figures")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. DATA PREPARATION: CLINICAL BURDEN AND TRANSMISSION STATUS
# ------------------------------------------------------------------------------
df_roc_prep <- df %>%
  mutate(
    # Native weekly clinical burden (continuous score)
    pressure_num  = as.numeric(nb_bmr_semaine),
    # Confirmed secondary genomic cross-transmission (<= 30 days, intra-ward)
    event_tx      = ifelse(is_cross_transmission_30d %in% c(TRUE, "TRUE", 1, "1"), 1, 0),
    ward_clean    = case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med", cores,  ignore.case = TRUE) ~ "Medical ICU",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pressure_num), !is.na(event_tx), !is.na(ward_clean))

cat(sprintf("--> Analyzable episodes for ROC analysis: %d\n", nrow(df_roc_prep)))

# ------------------------------------------------------------------------------
# 2. STRATIFIED ROC ANALYSIS WITH SEPARATION-PROOF ODDS RATIOS
# ------------------------------------------------------------------------------
compute_roc_by_ward <- function(ward_name, data) {
  sub_data <- data %>% filter(ward_clean == ward_name)

  # Direction: controls < cases
  roc_obj <- roc(sub_data$event_tx, sub_data$pressure_num, quiet = TRUE, direction = "<")

  # 2,000 stratified bootstrap replicates for the AUC 95% CI
  set.seed(42)
  ci_auc <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000, progress = "none")

  # Unconstrained optimal Youden index (J = Sensitivity + Specificity - 1)
  coords_opt <- coords(roc_obj, x = "best", method = "youden",
                       ret = c("threshold", "specificity", "sensitivity", "youden"))

  optimal_cutoff <- coords_opt$threshold[1]

  # Dichotomise exposure at the optimal cutoff
  sub_data <- sub_data %>% mutate(above_cutoff = ifelse(pressure_num >= optimal_cutoff, 1, 0))

  # 2x2 contingency table (Above cutoff x Transmission event)
  tab_2x2 <- table(Exposure = factor(sub_data$above_cutoff, levels = c(0, 1)),
                   Event    = factor(sub_data$event_tx,     levels = c(0, 1)))

  has_zero_cell <- any(tab_2x2 == 0)

  if (has_zero_cell) {
    # Exact median-unbiased Fisher estimation (avoids infinite OR)
    fisher_res <- fisher.test(tab_2x2)
    or_val     <- as.numeric(fisher_res$estimate)
    ci_or      <- as.numeric(fisher_res$conf.int)
    p_val      <- fisher_res$p.value
    or_method  <- "Fisher exact (Median-Unbiased)"
  } else {
    mod_glm <- glm(event_tx ~ above_cutoff, data = sub_data, family = binomial)
    or_val  <- exp(coef(mod_glm))[2]
    ci_or   <- exp(confint.default(mod_glm))[2, ]
    p_val   <- summary(mod_glm)$coefficients[2, 4]
    or_method <- "Logistic regression (Wald)"
  }

  list(
    roc_obj   = roc_obj,
    ci_auc    = ci_auc,
    threshold = optimal_cutoff,
    sens      = coords_opt$sensitivity[1],
    spec      = coords_opt$specificity[1],
    youden_j  = coords_opt$youden[1],
    or        = or_val,
    ci_or     = ci_or,
    p_val     = p_val,
    method    = or_method,
    sub_data  = sub_data,
    summary   = tibble(
      Ward             = ward_name,
      AUC              = round(as.numeric(auc(roc_obj)), 3),
      AUC_95_CI        = sprintf("[%.3f\u2013%.3f]", ci_auc[1], ci_auc[3]),
      Critical_Cutoff  = optimal_cutoff,
      Sensitivity_Pct  = round(coords_opt$sensitivity[1] * 100, 1),
      Specificity_Pct  = round(coords_opt$specificity[1] * 100, 1),
      Youden_Index     = round(coords_opt$youden[1], 3),
      Odds_Ratio       = round(or_val, 2),
      OR_95_CI         = sprintf("[%.2f\u2013%.2f]", ci_or[1], ci_or[2]),
      P_Value          = ifelse(p_val < 0.001, "< 0.001", sprintf("%.4f", p_val)),
      Estimation_Model = or_method
    )
  )
}

res_sicu <- compute_roc_by_ward("Surgical ICU", df_roc_prep)
res_micu <- compute_roc_by_ward("Medical ICU",  df_roc_prep)

roc_summary_table <- bind_rows(res_sicu$summary, res_micu$summary)

cat("\n=======================================================\n")
cat("DYNAMIC ROC & YOUDEN RESULTS SUMMARY:\n")
print(as.data.frame(roc_summary_table))
cat("=======================================================\n\n")

# ------------------------------------------------------------------------------
# 3. DYNAMIC PLOTTING: FIGURE 3A (ROC) AND FIGURE 3B (ODDS RATIOS)
# ------------------------------------------------------------------------------
roc_df_sicu <- tibble(FPR = 1 - res_sicu$roc_obj$specificities,
                      TPR = res_sicu$roc_obj$sensitivities, Unit = "Surgical ICU")
roc_df_micu <- tibble(FPR = 1 - res_micu$roc_obj$specificities,
                      TPR = res_micu$roc_obj$sensitivities, Unit = "Medical ICU")
roc_curves  <- bind_rows(roc_df_sicu, roc_df_micu)

youden_pts <- tibble(
  FPR  = c(1 - res_sicu$spec, 1 - res_micu$spec),
  TPR  = c(res_sicu$sens, res_micu$sens),
  Unit = c("Surgical ICU", "Medical ICU")
)

p3A <- ggplot(roc_curves, aes(x = FPR, y = TPR, color = Unit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.2) +
  geom_point(data = youden_pts, aes(x = FPR, y = TPR, color = Unit), size = 3.5, shape = 19) +
  annotate(
    "label", x = 0.97, y = 0.05, hjust = 1, vjust = 0, size = 2.8, fill = alpha("white", 0.95),
    label = sprintf(
      "Optimal Inflection Thresholds:\n\u2022 SICU: >= %.1f cases/wk (AUC = %.3f, Sens: %.1f%%, Spec: %.1f%%)\n\u2022 MICU: >= %.1f cases/wk (AUC = %.3f, Sens: %.1f%%, Spec: %.1f%%)",
      res_sicu$threshold, res_sicu$summary$AUC, res_sicu$summary$Sensitivity_Pct, res_sicu$summary$Specificity_Pct,
      res_micu$threshold, res_micu$summary$AUC, res_micu$summary$Sensitivity_Pct, res_micu$summary$Specificity_Pct
    )
  ) +
  scale_color_manual(values = c("Surgical ICU" = "#e41a1c", "Medical ICU" = "#377eb8")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title    = "A  ROC Curves: Colonization Pressure Thresholds",
    subtitle = "Youden index inflection points for ward cross-transmission risk",
    x        = "1 - Specificity (False Positive Rate)",
    y        = "Sensitivity (True Positive Rate)"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# PANEL 3B: STEP-FUNCTION ODDS RATIOS (FINITE AXIS AND Inf HANDLING)
# ------------------------------------------------------------------------------
df_or_plot <- tibble(
  Param = factor(c(sprintf("Surgical ICU\n(>= %.1f cases/wk)", res_sicu$threshold),
                   sprintf("Medical ICU\n(>= %.1f cases/wk)", res_micu$threshold)),
                 levels = c(sprintf("Surgical ICU\n(>= %.1f cases/wk)", res_sicu$threshold),
                            sprintf("Medical ICU\n(>= %.1f cases/wk)", res_micu$threshold))),
  OR    = c(res_sicu$or, res_micu$or),
  Lower = c(res_sicu$ci_or[1], res_micu$ci_or[1]),
  Upper = c(res_sicu$ci_or[2], res_micu$ci_or[2]),
  Color = c("#e41a1c", "#377eb8")
) %>%
  mutate(
    # Infinite-bound detection
    Is_Inf_Upper = !is.finite(Upper),
    # Label text (displays 'Inf' when the upper bound is infinite)
    Label_Text   = if_else(
      Is_Inf_Upper,
      sprintf("OR = %.2f\n[%.2f\u2013Inf]", OR, Lower),
      sprintf("OR = %.2f\n[%.2f\u2013%.2f]", OR, Lower, Upper)
    ),
    # Upper bound for the error bar (capped at 30 if infinite)
    Upper_Plot   = if_else(Is_Inf_Upper, 30, Upper)
  )

max_val_finite <- max(df_or_plot$Upper_Plot, na.rm = TRUE)
max_x_axis     <- max(35, ceiling(max_val_finite * 1.25))

p3B <- ggplot(df_or_plot, aes(x = OR, y = Param, color = Color)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper_Plot), height = 0.25, linewidth = 0.9) +
  geom_point(size = 3.8, shape = 18) +
  geom_text(aes(label = Label_Text), vjust = -0.65, size = 3.0, fontface = "bold") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, max_x_axis), breaks = seq(0, max_x_axis, by = 10)) +
  labs(
    title    = "B  Step-Function Increase in Transmission Odds",
    subtitle = "Odds ratios of active transmission upon exceeding critical threshold",
    x        = "Adjusted Odds Ratio (95% CI)",
    y        = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(plot.margin = margin(6, 30, 6, 6))

fig_3AB <- p3A | p3B

# ------------------------------------------------------------------------------
# 4. EXPORT TO STANDARDISED PUBLICATION PATHS
# ------------------------------------------------------------------------------
ggsave(file.path(out_dir_fig, "Figure_3AB_ROC_Colonization_Pressure.pdf"), fig_3AB, width = 12, height = 5.2, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_3AB_ROC_Colonization_Pressure.png"), fig_3AB, width = 12, height = 5.2, dpi = 300)

write_csv(roc_summary_table, file.path(out_dir_tab, "Supplementary_Table_S5_Logistic_Models.csv"))

cat(sprintf("\n--> Figure 3AB exported to: %s\n", out_dir_fig))
cat(sprintf("--> Supplementary Table S5 exported to: %s\n", out_dir_tab))
