# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 10_cox_transmission_velocity_pairs.R
# Purpose: Pairwise transmission network reconstruction, multivariable Cox
#          proportional hazards modeling of acute cross-transmission velocity
#          (<= 30 days), and publication-grade Forest Plot (Figure 3C).
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Figure_3C_Cox_Forest.pdf (.png)
#          output/tables/Table_Cox_Proportional_Hazards.csv
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
  library(scales)
  library(readr)
})

base_dir    <- "."
data_dir    <- file.path(base_dir, "data")
out_dir_fig <- file.path(base_dir, "output", "figures")
out_dir_tab <- file.path(base_dir, "output", "tables")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
df <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. RECONSTRUCTION OF DIRECTIONAL INTRA-WARD TRANSMISSION PAIRS
# ------------------------------------------------------------------------------
# Index episode per patient within each genomic cluster
df_clones <- df %>%
  mutate(
    date_sampling = as.Date(date_sampling),
    ward_clean    = case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med", cores, ignore.case = TRUE)  ~ "Medical ICU",
      TRUE ~ NA_character_
    ),
    is_valid_cl = is_valid_cluster %in% c(TRUE, "TRUE", 1, "1")
  ) %>%
  filter(is_valid_cl == TRUE, !is.na(ward_clean), !is.na(date_sampling)) %>%
  group_by(Numero.patient, Taxon_group, snp_cluster_clean) %>%
  filter(date_sampling == min(date_sampling)) %>%
  slice(1) %>%
  ungroup()

# Directional self-join (predecessor -> secondary case)
df_historique <- df_clones %>%
  select(
    Patient_anc       = Numero.patient,
    Taxon_group,
    snp_cluster_clean,
    Date_anc          = date_sampling,
    ward_anc          = ward_clean
  )

# Table of acute transmission pairs (<= 30 days, same ward)
df_surv <- df_clones %>%
  inner_join(df_historique, by = c("Taxon_group", "snp_cluster_clean"), relationship = "many-to-many") %>%
  filter(
    Date_anc < date_sampling,
    Numero.patient != Patient_anc,
    ward_clean == ward_anc
  ) %>%
  mutate(
    delai_jours = as.numeric(difftime(date_sampling, Date_anc, units = "days")),
    Evenement   = 1 # Each pair documents a real transmission link
  ) %>%
  filter(delai_jours <= 30)

cat(sprintf("--> Acute transmission pairs analyzed (<= 30 days): %d\n", nrow(df_surv)))
print(table(df_surv$ward_clean, df_surv$Taxon_group))

# ------------------------------------------------------------------------------
# 2. COX PROPORTIONAL HAZARDS MODEL (E. COLI AS REFERENCE)
# ------------------------------------------------------------------------------
df_surv <- df_surv %>%
  mutate(
    # Ward reference: Surgical ICU (measures the ratio in Medical ICU)
    Ward_factor    = factor(ward_clean, levels = c("Surgical ICU", "Medical ICU")),
    # Species reference: Escherichia coli (community comparator)
    Species_factor = factor(Taxon_group, levels = c("Escherichia coli",
                                                    "Klebsiella pneumoniae",
                                                    "Enterobacter cloacae complex"))
  )

modele_cox <- coxph(Surv(delai_jours, Evenement) ~ Ward_factor + Species_factor, data = df_surv)
summary_cox <- summary(modele_cox)
print(summary_cox)

# Proportional hazards assumption validation
ph_test <- cox.zph(modele_cox)
print(ph_test)

# ------------------------------------------------------------------------------
# 3. DYNAMIC EXTRACTION OF COVARIATE METRICS FOR THE FOREST PLOT
# ------------------------------------------------------------------------------
coef_tab <- summary_cox$coefficients
conf_tab <- summary_cox$conf.int

get_cox_metric <- function(var_pattern, metric_col, target_matrix) {
  row_match <- grep(var_pattern, rownames(target_matrix))
  if (length(row_match) == 1) {
    return(as.numeric(target_matrix[row_match, metric_col]))
  }
  return(NA_real_)
}

df_forest <- tibble(
  Covariate = c(
    "Medical ICU vs Surgical ICU",
    "Klebsiella pneumoniae",
    "Enterobacter cloacae complex",
    "Escherichia coli (Reference)"
  ),
  aHR   = c(
    get_cox_metric("Medical", "exp(coef)", conf_tab),
    get_cox_metric("pneumoniae", "exp(coef)", conf_tab),
    get_cox_metric("cloacae", "exp(coef)", conf_tab),
    1.00
  ),
  Lower = c(
    get_cox_metric("Medical", "lower .95", conf_tab),
    get_cox_metric("pneumoniae", "lower .95", conf_tab),
    get_cox_metric("cloacae", "lower .95", conf_tab),
    NA_real_
  ),
  Upper = c(
    get_cox_metric("Medical", "upper .95", conf_tab),
    get_cox_metric("pneumoniae", "upper .95", conf_tab),
    get_cox_metric("cloacae", "upper .95", conf_tab),
    NA_real_
  ),
  P_Val = c(
    get_cox_metric("Medical", "Pr(>|z|)", coef_tab),
    get_cox_metric("pneumoniae", "Pr(>|z|)", coef_tab),
    get_cox_metric("cloacae", "Pr(>|z|)", coef_tab),
    NA_real_
  ),
  IsRef = c(FALSE, FALSE, FALSE, TRUE)
) %>%
  mutate(
    Covariate = factor(Covariate, levels = rev(c(
      "Medical ICU vs Surgical ICU",
      "Klebsiella pneumoniae",
      "Enterobacter cloacae complex",
      "Escherichia coli (Reference)"
    ))),
    Label_Text = ifelse(
      IsRef,
      "Reference (1.00)",
      sprintf("aHR = %.2f [%.2f-%.2f], %s",
              aHR, Lower, Upper,
              ifelse(P_Val < 0.001, "p < 0.001", sprintf("p = %.4f", P_Val)))
    )
  )

# Dynamic X-axis limit without truncation
max_upper_cox <- max(df_forest$Upper, na.rm = TRUE)
lim_x_cox     <- max(3.0, ceiling(max_upper_cox * 1.25))

# ------------------------------------------------------------------------------
# 4. FOREST PLOT RENDERING (PANEL 3C)
# ------------------------------------------------------------------------------
cox_labels_expr <- c(
  "Medical ICU vs Surgical ICU"  = "Medical ICU\nvs Surgical ICU",
  "Klebsiella pneumoniae"        = expression(italic("K. pneumoniae")),
  "Enterobacter cloacae complex" = expression(italic("E. cloacae")~"complex"),
  "Escherichia coli (Reference)" = expression(italic("E. coli")~"(Reference)")
)

p3C <- ggplot(df_forest, aes(x = aHR, y = Covariate)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    data = filter(df_forest, !IsRef),
    aes(xmin = Lower, xmax = Upper),
    height = 0.22, linewidth = 0.85, color = "#084594"
  ) +
  geom_point(data = filter(df_forest, !IsRef), size = 3.6, shape = 15, color = "#084594") +
  geom_point(data = filter(df_forest, IsRef), size = 3.0, shape = 18, color = "grey40") +
  geom_text(aes(label = Label_Text), vjust = -0.75, size = 3.1, fontface = "bold", color = "black") +
  scale_y_discrete(labels = cox_labels_expr) +
  scale_x_continuous(limits = c(0.7, lim_x_cox), breaks = pretty_breaks(n = 5)) +
  labs(
    title    = "C  Multivariable Cox Proportional Hazards Model",
    subtitle = "Independent predictors of acute cross-transmission velocity (<= 30-day intra-ward window)",
    x        = "Adjusted Hazard Ratio (95% CI)",
    y        = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y        = element_text(face = "bold", size = rel(0.92), lineheight = 0.9),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(6, 25, 6, 6)
  )

# Final exports
ggsave(file.path(out_dir_fig, "Figure_3C_Cox_Forest.pdf"), p3C, width = 10.0, height = 4.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_3C_Cox_Forest.png"), p3C, width = 10.0, height = 4.5, dpi = 300)

write_csv(df_forest %>% select(-IsRef), file.path(out_dir_tab, "Table_Cox_Proportional_Hazards.csv"))

cat(sprintf("\n--> Figure 3C exported to: %s\n", out_dir_fig))
cat(sprintf("--> Cox model table exported to: %s\n", out_dir_tab))
