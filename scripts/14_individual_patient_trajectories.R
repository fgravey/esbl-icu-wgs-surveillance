# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 14_individual_patient_trajectories.R
# Purpose: Comprehensive patient-level clinical and epidemiological evaluation:
#          1. Length of stay / acquisition delay distributions & Wilcoxon tests
#          2. Non-hierarchical anatomical infection sites breakdown
#          3. Carriage-to-infection conversion rates by species
#          4. Multivariable Cox models of incident acquisition (> 48h)
#          5. Multivariable logistic regression of secondary clinical infection
#          6. Native ggplot2 Kaplan-Meier curve with risk table (no survminer
#             dependency)
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Figure_Individual_Acquisition_Survival.pdf (.png)
#          output/tables/Table_Delais_Acquisition_Especes.csv
#          output/tables/Table_Sites_Infectieux_Details.csv
#          output/tables/Table_Taux_Conversion_Especes.csv
#          output/tables/Table_Individual_Cox_Acquisition.csv
#          output/tables/Table_Individual_Infection_Risk.csv
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(survival)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
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
# 1. COMPREHENSIVE CLINICAL PREPARATION OF THE INDIVIDUAL COHORT (N = 689)
# ------------------------------------------------------------------------------
# NOTE: several matching patterns below target French free-text labels of the
# source database (e.g. "Masculin"/"Feminin" for sex, French infection-site
# keywords). They are kept for reproducibility and must not be translated.
df_indiv <- df %>%
  mutate(
    # Hospital length-of-stay delays
    communoso_num = as.numeric(as.character(communoso)),
    los_days      = pmax(0, ifelse(is.na(communoso_num), 0, communoso_num)),

    # Classical binary 48h classification
    statut_48h    = ifelse(los_days <= 2, "Imported (<= 48h)", "Acquired (> 48h)"),
    is_acquired   = ifelse(los_days > 2, 1, 0),

    # Clinical infection status (0 = carriage only, >= 1 = infection)
    inf_count     = as.numeric(as.character(Infection.BLSE)),
    has_infection = ifelse(!is.na(inf_count) & inf_count > 0, 1, 0),

    # Robust sex coding (handles French labels "Masculin" / "Feminin")
    sexe_clean = case_when(
      grepl("masc|homme|^m$", Sexe, ignore.case = TRUE) ~ "M",
      grepl("fem|^f$",       Sexe, ignore.case = TRUE) ~ "F",
      TRUE                                             ~ NA_character_
    ),
    sexe_clean = factor(sexe_clean, levels = c("M", "F")),

    # Species and ward factors
    Species_factor = factor(case_when(
      Taxon_group == "Klebsiella pneumoniae"        ~ "Klebsiella pneumoniae",
      Taxon_group == "Enterobacter cloacae complex" ~ "Enterobacter cloacae complex",
      TRUE                                          ~ "Escherichia coli"
    ), levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex")),

    Ward_factor = factor(case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med", cores, ignore.case = TRUE)  ~ "Medical ICU",
      TRUE ~ NA_character_
    ), levels = c("Surgical ICU", "Medical ICU")),

    # Discretized global colonization pressure (threshold 4.5)
    cp_global_num = as.numeric(as.character(nb_bmr_semaine)),
    cp_high       = factor(ifelse(cp_global_num >= 4.5, "High (>= 4.5)", "Standard (< 4.5)"),
                           levels = c("Standard (< 4.5)", "High (>= 4.5)")),

    # Continuous species-specific colonization pressure
    cp_species_num = case_when(
      Taxon_group == "Escherichia coli"             ~ as.numeric(as.character(nb_ecoli_service_semaine)),
      Taxon_group == "Klebsiella pneumoniae"        ~ as.numeric(as.character(nb_kpn_service_semaine)),
      Taxon_group == "Enterobacter cloacae complex" ~ as.numeric(as.character(nb_ecc_service_semaine)),
      TRUE ~ 0
    ),

    # Clinical and baseline covariates
    age_num    = as.numeric(as.character(age)),
    effort_dep = ifelse(is.na(as.numeric(as.character(nb_depistage_patient))), 1,
                        as.numeric(as.character(nb_depistage_patient))),

    # Tripartite resistance block (blaCTX-M-15 + blaOXA-1 + blaTEM)
    amr_beta_str = paste(ifelse(is.na(AMR_BETA_LACTAM), "", AMR_BETA_LACTAM),
                         ifelse(is.na(Beta.lactam), "", Beta.lactam)),
    tripartite   = ifelse(str_detect(amr_beta_str, "CTX-M-15") &
                            str_detect(amr_beta_str, "OXA-1") &
                            str_detect(amr_beta_str, "TEM"), 1, 0)
  ) %>%
  filter(!is.na(Ward_factor), !is.na(sexe_clean)) %>%
  droplevels()

cat(sprintf("--> Complete individual cohort available for analysis: %d patients\n", nrow(df_indiv)))

# ------------------------------------------------------------------------------
# 2. HOSPITAL DELAYS AND NON-PARAMETRIC TESTS (FIG 1B LEGACY)
# ------------------------------------------------------------------------------
cat("\n=== 1. DISTRIBUTION OF ACQUISITION DELAYS BY SPECIES (DAYS) ===\n")
synthese_delais <- df_indiv %>%
  group_by(Species_factor) %>%
  summarise(
    N_Total       = n(),
    Mediane_LOS   = median(los_days),
    Q1_LOS        = quantile(los_days, 0.25),
    Q3_LOS        = quantile(los_days, 0.75),
    N_Importe     = sum(los_days <= 2),
    Pct_Importe   = round(100 * sum(los_days <= 2) / n(), 1),
    N_Acquis      = sum(los_days > 2),
    Pct_Acquis    = round(100 * sum(los_days > 2) / n(), 1),
    .groups       = "drop"
  )
print(as.data.frame(synthese_delais))
write_csv(synthese_delais, file.path(out_dir_tab, "Table_Delais_Acquisition_Especes.csv"))

kw_test <- kruskal.test(los_days ~ Species_factor, data = df_indiv)
cat(sprintf("Kruskal-Wallis p-value: %.4e\n", kw_test$p.value))

# ------------------------------------------------------------------------------
# 3. TEXT-MINING OF ANATOMICAL INFECTION SITES (FIG 4C / TABLE S6 LEGACY)
# ------------------------------------------------------------------------------
cat("\n=== 2. DETAILED BREAKDOWN OF ANATOMICAL INFECTION SITES ===\n")
df_sites_infectieux <- df_indiv %>%
  filter(has_infection == 1) %>%
  mutate(nature_brute = tolower(trimws(as.character(infection.nature)))) %>%
  filter(!nature_brute %in% c("-", "", "na", "nr", "unknown")) %>%
  separate_rows(nature_brute, sep = "[,;/]") %>%
  mutate(nature_brute = trimws(nature_brute)) %>%
  filter(nature_brute != "") %>%
  mutate(
    Site_Standard = case_when(
      str_detect(nature_brute, "bacterie|hemoc|sepsis")             ~ "Bacteremia / Sepsis",
      str_detect(nature_brute, "poumon|pneumo|asp|lba|trache|pavm") ~ "Pneumonia (VAP)",
      str_detect(nature_brute, "urin")                              ~ "Urinary Tract Infection",
      str_detect(nature_brute, "abdom|adbom|periton|bile")          ~ "Intra-Abdominal Infection",
      str_detect(nature_brute, "catheter|kt")                       ~ "Catheter-Related Infection",
      str_detect(nature_brute, "plaie|pus|paroi|os|artic|peau")     ~ "Skin, Soft Tissue and Bone",
      TRUE                                                          ~ "Other site"
    )
  )

synthese_sites <- df_sites_infectieux %>%
  group_by(Ward_factor, Species_factor, Site_Standard) %>%
  summarise(Nombre_Foyers = n(), .groups = "drop") %>%
  arrange(Ward_factor, Species_factor, desc(Nombre_Foyers))

print(head(as.data.frame(synthese_sites), 15))
write_csv(synthese_sites, file.path(out_dir_tab, "Table_Sites_Infectieux_Details.csv"))

# ------------------------------------------------------------------------------
# 4. CARRIAGE -> INFECTION SECONDARY CONVERSION RATES BY SPECIES
# ------------------------------------------------------------------------------
cat("\n=== 3. SECONDARY CONVERSION RATES BY SPECIES ===\n")
taux_conversion <- df_indiv %>%
  group_by(Species_factor) %>%
  summarise(
    Total_Porteurs   = n(),
    Infections_Obs   = sum(has_infection),
    Taux_Conversion  = round(100 * sum(has_infection) / n(), 2),
    .groups          = "drop"
  )
print(as.data.frame(taux_conversion))
write_csv(taux_conversion, file.path(out_dir_tab, "Table_Taux_Conversion_Especes.csv"))

# ------------------------------------------------------------------------------
# 5. MODEL 1: RISK OF INCIDENT NOSOCOMIAL ACQUISITION (> 48H)
#    (Strict exclusion of cases prevalent at admission)
# ------------------------------------------------------------------------------
df_incid <- df_indiv %>%
  filter(los_days > 2) %>%
  mutate(
    time_incident = los_days - 2,
    event_acq     = 1
  ) %>%
  droplevels()

cat(sprintf("\n--> Incident cohort for Cox regression: %d patients\n", nrow(df_incid)))

# Model 1A: discretized global pressure (>= 4.5 cases/week)
cox_acq_global <- coxph(
  Surv(time_incident, event_acq) ~ Species_factor + Ward_factor + age_num + sexe_clean + effort_dep + cp_high,
  data = df_incid
)
cat("\n--- MODEL 1A: COX ON ACQUISITION (GLOBAL PRESSURE) ---\n")
print(summary(cox_acq_global))

# Model 1B: species-specific pressure
cox_acq_species <- coxph(
  Surv(time_incident, event_acq) ~ Species_factor + Ward_factor + age_num + sexe_clean + effort_dep + cp_species_num,
  data = df_incid
)
cat("\n--- MODEL 1B: COX ON ACQUISITION (SPECIES-SPECIFIC PRESSURE) ---\n")
print(summary(cox_acq_species))

# Summary table of the Cox models
extract_cox_table <- function(mod, label) {
  s <- summary(mod)
  tibble(
    Model     = label,
    Covariate = rownames(s$coefficients),
    aHR       = round(s$conf.int[, "exp(coef)"], 2),
    Lower_95  = round(s$conf.int[, "lower .95"], 2),
    Upper_95  = round(s$conf.int[, "upper .95"], 2),
    P_Value   = ifelse(s$coefficients[, "Pr(>|z|)"] < 0.001, "< 0.001",
                       sprintf("%.4f", s$coefficients[, "Pr(>|z|)"]))
  )
}

table_acq_cox <- bind_rows(
  extract_cox_table(cox_acq_global, "Model 1A (Global pressure >= 4.5)"),
  extract_cox_table(cox_acq_species, "Model 1B (Species-specific pressure)")
)
write_csv(table_acq_cox, file.path(out_dir_tab, "Table_Individual_Cox_Acquisition.csv"))

# ------------------------------------------------------------------------------
# 6. MODEL 2: FACTORS ASSOCIATED WITH INFECTIOUS COMPLICATION IN CARRIERS
# ------------------------------------------------------------------------------
cat("\n=== 4. MULTIVARIABLE LOGISTIC MODEL: RISK OF CLINICAL INFECTION ===\n")
mod_infection_risk <- glm(
  has_infection ~ Species_factor + Ward_factor + age_num + sexe_clean + tripartite + cp_high,
  data = df_indiv,
  family = binomial
)
summary_inf <- summary(mod_infection_risk)
print(summary_inf)

ci_inf <- exp(confint.default(mod_infection_risk))
table_inf_risk <- tibble(
  Covariate  = names(coef(mod_infection_risk)),
  Odds_Ratio = round(exp(coef(mod_infection_risk)), 2),
  Lower_95   = round(ci_inf[, 1], 2),
  Upper_95   = round(ci_inf[, 2], 2),
  P_Value    = ifelse(summary_inf$coefficients[, 4] < 0.001, "< 0.001",
                      sprintf("%.4f", summary_inf$coefficients[, 4]))
) %>%
  filter(Covariate != "(Intercept)")

write_csv(table_inf_risk, file.path(out_dir_tab, "Table_Individual_Infection_Risk.csv"))

# ------------------------------------------------------------------------------
# 7. KAPLAN-MEIER SURVIVAL CURVES AND NATIVE RISK TABLE (WITHOUT SURVMINER)
# ------------------------------------------------------------------------------
# Non-parametric survival model
fit_km <- survfit(Surv(los_days, is_acquired) ~ Species_factor, data = df_indiv)

# Extract coordinates for ggplot
summary_km <- summary(fit_km)
df_km_plot <- tibble(
  time   = summary_km$time,
  surv   = summary_km$surv,
  lower  = summary_km$lower,
  upper  = summary_km$upper,
  strata = summary_km$strata
) %>%
  mutate(
    Species = gsub("Species_factor=", "", strata),
    Species = factor(Species, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  )

# Initial point at t = 0 for each species
df_km_init <- tibble(
  time   = rep(0, 3),
  surv   = rep(1, 3),
  lower  = rep(1, 3),
  upper  = rep(1, 3),
  strata = paste0("Species_factor=", c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex")),
  Species = factor(c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"),
                   levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
)

df_km_all <- bind_rows(df_km_init, df_km_plot) %>%
  arrange(Species, time)

# Global log-rank test
lr_test <- survdiff(Surv(los_days, is_acquired) ~ Species_factor, data = df_indiv)
lr_pval <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
lr_label <- ifelse(lr_pval < 0.001, "Log-rank p < 0.001", sprintf("Log-rank p = %.4f", lr_pval))

palette_species <- c(
  "Escherichia coli"             = "#1f78b4",
  "Klebsiella pneumoniae"        = "#e31a1c",
  "Enterobacter cloacae complex" = "#ff7f00"
)

# Main survival plot
p_km_main <- ggplot(df_km_all, aes(x = time, y = surv, color = Species)) +
  geom_step(linewidth = 1.1) +
  annotate(
    "label", x = 50, y = 0.85, label = lr_label,
    size = 3.5, fontface = "bold", fill = alpha("white", 0.9), color = "black"
  ) +
  scale_color_manual(
    values = palette_species,
    labels = c(expression(italic("E. coli")),
               expression(italic("K. pneumoniae")),
               expression(italic("E. cloacae")~"complex")),
    name = "Species:"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.02)) +
  scale_x_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  labs(
    title    = "Individual Kinetics of Hospital-Acquired ESBL-E Detection",
    subtitle = sprintf("Kaplan-Meier estimates of remaining uncolonized post-admission (N = %d episodes)", nrow(df_indiv)),
    x        = NULL,
    y        = "Uncolonized Probability"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = rel(1.05)),
    axis.text       = element_text(color = "black")
  )

# Risk table at regular intervals
time_ticks <- seq(0, 60, by = 15)

summary_risk <- summary(fit_km, times = time_ticks, extend = TRUE)
df_risk_table <- tibble(
  time    = summary_risk$time,
  n_risk  = summary_risk$n.risk,
  strata  = summary_risk$strata
) %>%
  mutate(
    Species = gsub("Species_factor=", "", strata),
    Species = factor(Species, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  )

p_km_risk <- ggplot(df_risk_table, aes(x = time, y = Species, label = n_risk)) +
  geom_text(size = 3.2, fontface = "bold", color = "black") +
  scale_y_discrete(
    labels = c(
      "Escherichia coli"             = expression(italic("E. coli")),
      "Klebsiella pneumoniae"        = expression(italic("K. pneumoniae")),
      "Enterobacter cloacae complex" = expression(italic("E. cloacae")~"c.")
    ),
    limits = rev(c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  ) +
  scale_x_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) +
  labs(
    x = "Hospital length of stay before detection (Days)",
    y = "Number at risk"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(face = "bold", size = rel(0.9)),
    axis.text.y  = element_text(face = "bold", size = rel(0.85), color = "black"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# Composition via patchwork
km_composed <- p_km_main / p_km_risk + plot_layout(heights = c(3.2, 1))

ggsave(file.path(out_dir_fig, "Figure_Individual_Acquisition_Survival.pdf"), km_composed, width = 9.0, height = 7.0, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_Individual_Acquisition_Survival.png"), km_composed, width = 9.0, height = 7.0, dpi = 300)

cat(sprintf("\n--> All individual-level analyses executed and exported to:\n%s\n", out_dir_fig))
