# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 08_figure3_tipping_points_cox.R
# Purpose: Dynamic generation of Figure 3 (three panels):
#          - 3A: Empirical ROC curves and Youden inflection thresholds of
#                ward colonization pressure
#          - 3B: Step-function odds ratios upon exceeding pressure thresholds
#          - 3C: Multivariable Cox proportional hazards model of
#                cross-transmission velocity (forest plot)
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Figure_3_Tipping_Points_Cox.pdf (.png, .tiff)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(survival)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_fig <- file.path(base_dir, "output", "figures")
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# Cleaning and standardisation
df_prep <- df %>%
  mutate(
    ward_clean = case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med", cores,  ignore.case = TRUE) ~ "Medical ICU",
      TRUE ~ NA_character_
    ),
    is_trans  = ifelse(vraie_transmission_croisee %in% c(TRUE, "TRUE", 1, "1"), 1, 0),
    cp_num    = as.numeric(as.character(nb_bmr_semaine)),
    date_obj  = as.Date(date_sampling),
    is_covid  = ifelse(date_obj >= as.Date("2020-03-01") & date_obj <= as.Date("2021-06-30"), 1, 0),
    Taxon_factor = factor(case_when(
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE) ~ "Klebsiella pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)    ~ "Enterobacter cloacae complex",
      TRUE                                                 ~ "Escherichia coli"
    ), levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  ) %>%
  filter(!is.na(ward_clean) & !is.na(cp_num))

# ------------------------------------------------------------------------------
# 1. PANEL 3A: ROC ANALYSIS AND YOUDEN THRESHOLDS
# ------------------------------------------------------------------------------
calc_roc_data <- function(sub_df, unit_name) {
  roc_obj <- roc(sub_df$is_trans, sub_df$cp_num, quiet = TRUE, direction = "<")
  auc_val <- round(auc(roc_obj), 3)
  coords_best <- coords(roc_obj, x = "best", method = "youden",
                        ret = c("threshold", "specificity", "sensitivity"))

  roc_curve_df <- tibble(
    FPR  = 1 - roc_obj$specificities,
    TPR  = roc_obj$sensitivities,
    Unit = sprintf("%s (AUC = %.3f)", unit_name, auc_val)
  ) %>% arrange(FPR, TPR)

  youden_pt <- tibble(
    FPR       = 1 - coords_best$specificity,
    TPR       = coords_best$sensitivity,
    Unit      = sprintf("%s (AUC = %.3f)", unit_name, auc_val),
    Threshold = coords_best$threshold,
    Sens      = round(coords_best$sensitivity * 100, 1),
    Spec      = round(coords_best$specificity * 100, 1)
  )
  list(curve = roc_curve_df, point = youden_pt, auc = auc_val, threshold = coords_best$threshold)
}

roc_sicu <- calc_roc_data(filter(df_prep, ward_clean == "Surgical ICU"), "Surgical ICU")
roc_micu <- calc_roc_data(filter(df_prep, ward_clean == "Medical ICU"),  "Medical ICU")

df_roc_plot <- bind_rows(roc_sicu$curve, roc_micu$curve)
df_pts_plot <- bind_rows(roc_sicu$point, roc_micu$point)

unit_colors_roc <- setNames(c("#e41a1c", "#377eb8"), c(roc_sicu$point$Unit, roc_micu$point$Unit))

sicu_label_text <- sprintf("\u2022 SICU: >= %.1f cases/wk (Sens: %.1f%%, Spec: %.1f%%)",
                           roc_sicu$point$Threshold, roc_sicu$point$Sens, roc_sicu$point$Spec)
micu_label_text <- sprintf("\u2022 MICU: >= %.1f cases/wk (Sens: %.1f%%, Spec: %.1f%%)",
                           roc_micu$point$Threshold, roc_micu$point$Sens, roc_micu$point$Spec)

p3A <- ggplot(df_roc_plot, aes(x = FPR, y = TPR, color = Unit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.2) +
  geom_point(data = df_pts_plot, aes(x = FPR, y = TPR, color = Unit), size = 3.5, shape = 19) +
  annotate(
    "label", x = 0.97, y = 0.05,
    label = paste("Optimal Inflection Thresholds:\n", sicu_label_text, "\n", micu_label_text),
    hjust = 1, vjust = 0, size = 2.8, fill = alpha("white", 0.92), color = "black",
    fontface = "plain", label.size = 0.3
  ) +
  scale_color_manual(values = unit_colors_roc) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(
    title    = "A  ROC Curves: Colonization Pressure Thresholds",
    subtitle = "Youden index inflection points for ward cross-transmission risk",
    x        = "1 - Specificity (False Positive Rate)",
    y        = "Sensitivity (True Positive Rate)",
    color    = "Unit:"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = rel(1.02)),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35)
  )

# ------------------------------------------------------------------------------
# 2. PANEL 3B: STEP-FUNCTION ODDS RATIOS (LOGISTIC REGRESSION)
# ------------------------------------------------------------------------------
calc_or <- function(sub_df, thresh, label) {
  dat  <- sub_df %>% mutate(over = ifelse(cp_num >= thresh, 1, 0))
  mod  <- glm(is_trans ~ over, data = dat, family = binomial)
  est  <- exp(coef(mod)["over"])
  ci   <- exp(confint.default(mod)["over", ])
  pval <- summary(mod)$coefficients["over", "Pr(>|z|)"]
  p_txt <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.4f", pval)
  tibble(Param = label, OR = est, Lower = ci[1], Upper = ci[2], Pval = p_txt)
}

df_threshold_or <- bind_rows(
  calc_or(filter(df_prep, ward_clean == "Surgical ICU"), roc_sicu$point$Threshold,
          sprintf("Surgical ICU\n(>= %.1f cases/wk)", roc_sicu$point$Threshold)),
  calc_or(filter(df_prep, ward_clean == "Medical ICU"), roc_micu$point$Threshold,
          sprintf("Medical ICU\n(>= %.1f cases/wk)", roc_micu$point$Threshold))
) %>%
  mutate(Param = factor(Param, levels = Param))

p3B <- ggplot(df_threshold_or, aes(x = OR, y = Param, color = Param)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.25, linewidth = 0.9) +
  geom_point(size = 3.8, shape = 18) +
  geom_text(
    aes(label = sprintf("OR = %.2f\n[%.2f\u2013%.2f], %s", OR, Lower, Upper, Pval)),
    vjust = -0.65, size = 3.0, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  scale_x_continuous(limits = c(0, max(df_threshold_or$Upper) * 1.15)) +
  labs(
    title    = "B  Step-Function Increase in Transmission Odds",
    subtitle = "Odds ratios of active transmission upon exceeding critical threshold",
    x        = "Adjusted Odds Ratio (95% CI)",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = rel(1.02)),
    axis.text.y        = element_text(face = "bold", size = rel(0.9)),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35)
  )

# ------------------------------------------------------------------------------
# 3. PANEL 3C: MULTIVARIABLE COX MODEL (TRANSMISSION VELOCITY)
# ------------------------------------------------------------------------------
df_trans_net <- df_prep %>%
  filter(is_trans == 1) %>%
  mutate(
    delay_days = pmax(1, as.numeric(as.character(communoso))),
    delay_days = ifelse(delay_days > 30, 30, delay_days),
    event      = 1,
    ward_factor = factor(ward_clean, levels = c("Medical ICU", "Surgical ICU"))
  )

cox_fit <- coxph(Surv(delay_days, event) ~ is_covid + ward_factor + Taxon_factor, data = df_trans_net)
s_cox   <- summary(cox_fit)

df_cox <- tibble(
  Variable = factor(c(
    "COVID-19 Surge Periods",
    "Surgical ICU vs Medical ICU",
    "Enterobacter cloacae complex",
    "Klebsiella pneumoniae",
    "Escherichia coli (Reference)"
  ), levels = rev(c(
    "COVID-19 Surge Periods",
    "Surgical ICU vs Medical ICU",
    "Enterobacter cloacae complex",
    "Klebsiella pneumoniae",
    "Escherichia coli (Reference)"
  ))),
  aHR   = c(exp(coef(cox_fit)["is_covid"]),
            exp(coef(cox_fit)["ward_factorSurgical ICU"]),
            exp(coef(cox_fit)["Taxon_factorEnterobacter cloacae complex"]),
            exp(coef(cox_fit)["Taxon_factorKlebsiella pneumoniae"]),
            1.00),
  Lower = c(s_cox$conf.int["is_covid", "lower .95"],
            s_cox$conf.int["ward_factorSurgical ICU", "lower .95"],
            s_cox$conf.int["Taxon_factorEnterobacter cloacae complex", "lower .95"],
            s_cox$conf.int["Taxon_factorKlebsiella pneumoniae", "lower .95"],
            NA),
  Upper = c(s_cox$conf.int["is_covid", "upper .95"],
            s_cox$conf.int["ward_factorSurgical ICU", "upper .95"],
            s_cox$conf.int["Taxon_factorEnterobacter cloacae complex", "upper .95"],
            s_cox$conf.int["Taxon_factorKlebsiella pneumoniae", "upper .95"],
            NA),
  IsRef = c(FALSE, FALSE, FALSE, FALSE, TRUE)
) %>%
  mutate(
    Label = ifelse(IsRef, "Reference (1.00)",
                   sprintf("aHR = %.2f [%.2f\u2013%.2f]", aHR, Lower, Upper))
  )

cox_labels_expr <- c(
  "COVID-19 Surge Periods"       = "COVID-19\nSurge Periods",
  "Surgical ICU vs Medical ICU"  = "Surgical ICU\nvs Medical ICU",
  "Enterobacter cloacae complex" = expression(italic("E. cloacae")~"complex"),
  "Klebsiella pneumoniae"        = expression(italic("K. pneumoniae")),
  "Escherichia coli (Reference)" = expression(italic("E. coli")~"(Reference)")
)

p3C <- ggplot(df_cox, aes(x = aHR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(data = filter(df_cox, !IsRef), aes(xmin = Lower, xmax = Upper),
                 height = 0.22, linewidth = 0.85, color = "#084594") +
  geom_point(data = filter(df_cox, !IsRef), size = 3.6, shape = 15, color = "#084594") +
  geom_point(data = filter(df_cox, IsRef), size = 3.0, shape = 18, color = "grey40") +
  geom_text(aes(label = Label), vjust = -0.75, size = 3.0, fontface = "bold", color = "black") +
  scale_y_discrete(labels = cox_labels_expr) +
  scale_x_continuous(limits = c(0.4, max(df_cox$Upper, na.rm = TRUE) * 1.25)) +
  labs(
    title    = "C  Multivariable Cox Proportional Hazards Model",
    subtitle = "Independent predictors of cross-transmission velocity",
    x        = "Adjusted Hazard Ratio (95% CI)",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = rel(1.02)),
    axis.text.y        = element_text(face = "bold", size = rel(0.85), lineheight = 0.9),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35)
  )

# ------------------------------------------------------------------------------
# 4. ASSEMBLY AND EXPORT
# ------------------------------------------------------------------------------
fig3_complete <- ((p3A | p3B) / p3C + plot_layout(heights = c(1, 1.1))) &
  theme(plot.margin = margin(t = 10, r = 15, b = 10, l = 20, unit = "pt"))

ggsave(file.path(out_dir_fig, "Figure_3_Tipping_Points_Cox.pdf"), fig3_complete, width = 12.5, height = 9.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_3_Tipping_Points_Cox.png"), fig3_complete, width = 12.5, height = 9.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_3_Tipping_Points_Cox.tiff"), fig3_complete, width = 12.5, height = 9.5, dpi = 300, compression = "lzw")

cat("--> Figure 3 regenerated successfully in:", out_dir_fig, "\n")
