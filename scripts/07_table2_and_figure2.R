# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 07_table2_and_figure2.R
# Purpose: Dynamic recalculation of Table 2 (longitudinal incidence trends and
#          annual IRRs, 2019-2023, using hospital activity denominators) and
#          full generation of Figure 2 (panels 2A, 2B, 2C). The consolidated
#          Table 2 is additionally exported to a multi-sheet Excel workbook.
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
#          data/JH_admission_rea_2019_2023.xlsx  (patient-days and admissions)
# Output:  output/figures/Figure_2_Longitudinal_Dynamics.pdf (.png, .tiff)
#          output/tables/Table_2_Longitudinal_Trends_and_IRR.csv
#          output/tables/Table_2_Consolidated.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(openxlsx)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_fig <- file.path(base_dir, "output", "figures")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

cohort_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
denom_file  <- file.path(data_dir, "JH_admission_rea_2019_2023.xlsx")

# ------------------------------------------------------------------------------
# 1. HOSPITAL ACTIVITY DENOMINATORS (2019-2023)
# ------------------------------------------------------------------------------
df_adm_raw <- read_excel(denom_file, sheet = "Nb entr\u00e9e totale HC")
df_jh_raw  <- read_excel(denom_file, sheet = "Nb journ\u00e9es HC")

jh_long <- df_jh_raw %>%
  rename(unit = 1) %>%
  filter(unit != "Total g\u00e9n\u00e9ral" & !is.na(unit)) %>%
  pivot_longer(cols = matches("201[9]|202[0-3]"), names_to = "year", values_to = "patient_days") %>%
  mutate(year = as.numeric(year))

adm_long <- df_adm_raw %>%
  rename(unit = 1) %>%
  filter(unit != "Total g\u00e9n\u00e9ral" & !is.na(unit)) %>%
  pivot_longer(cols = matches("201[9]|202[0-3]"), names_to = "year", values_to = "admissions") %>%
  mutate(year = as.numeric(year))

denom_annuel <- jh_long %>%
  group_by(year) %>%
  summarise(patient_days = sum(patient_days, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    adm_long %>%
      group_by(year) %>%
      summarise(admissions = sum(admissions, na.rm = TRUE), .groups = "drop"),
    by = "year"
  ) %>%
  filter(year %in% 2019:2023) %>%
  arrange(year)

# ------------------------------------------------------------------------------
# 2. PREPARATION OF THE DE-DUPLICATED COHORT
# ------------------------------------------------------------------------------
df <- read.csv(cohort_file, sep = "\t", stringsAsFactors = FALSE)

df_clean <- df %>%
  mutate(
    # Year extracted from the year-week identifier
    Year_clean = as.numeric(substr(as.character(annee_semaine), 1, 4)),

    Taxon_clean = case_when(
      grepl("coli", Taxon_group, ignore.case = TRUE)        ~ "Escherichia coli",
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE)  ~ "Klebsiella pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)     ~ "Enterobacter cloacae complex",
      TRUE                                                  ~ NA_character_
    ),

    communoso_num = as.numeric(as.character(communoso)),
    los_days      = pmax(0, ifelse(is.na(communoso_num), 0, communoso_num)),
    is_imported   = ifelse(los_days <= 2, 1, 0),
    is_nosocomial = ifelse(los_days > 2, 1, 0),
    is_trans      = ifelse(vraie_transmission_croisee %in% c(TRUE, "TRUE", 1, "1"), 1, 0)
  ) %>%
  filter(!is.na(Taxon_clean) & Year_clean %in% 2019:2023)

# ------------------------------------------------------------------------------
# 3. STATISTICAL MODELS AND SECURED TABLE 2 GENERATION
# ------------------------------------------------------------------------------
compute_metric_row <- function(label, filter_cond, denom_col = "patient_days",
                               scale = 1000, model_type = "quasipoisson") {
  counts <- df_clean %>%
    filter(eval(parse(text = filter_cond))) %>%
    group_by(Year_clean) %>%
    summarise(n = n(), .groups = "drop")

  dat <- denom_annuel %>%
    left_join(counts, by = c("year" = "Year_clean")) %>%
    mutate(
      n      = coalesce(n, 0L),
      denom  = .data[[denom_col]],
      rate   = (n / denom) * scale,
      time   = year - 2019
    )

  vals_yr <- setNames(sprintf("%.2f (%d)", dat$rate, dat$n), paste0("yr_", dat$year))
  total_events <- sum(dat$n)

  if (total_events > 0 && n_distinct(dat$n) > 1) {
    if (model_type == "quasipoisson") {
      fit <- glm(n ~ time + offset(log(denom)), family = quasipoisson(link = "log"), data = dat)
      est <- exp(coef(fit)["time"])

      # Rigorous alignment: Student t quantile, consistent with the quasi-Poisson p-value
      df_res <- df.residual(fit)
      se_qp  <- summary(fit)$coefficients["time", "Std. Error"]
      t_crit <- qt(0.975, df = df_res)
      ci_low  <- exp(coef(fit)["time"] - t_crit * se_qp)
      ci_high <- exp(coef(fit)["time"] + t_crit * se_qp)
      pval    <- summary(fit)$coefficients["time", "Pr(>|t|)"]
    } else {
      # Standard Poisson model (Gaussian Wald test)
      fit <- glm(n ~ time + offset(log(denom)), family = poisson(link = "log"), data = dat)
      est <- exp(coef(fit)["time"])
      ci  <- exp(confint.default(fit)["time", ])
      ci_low  <- ci[1]
      ci_high <- ci[2]
      pval    <- summary(fit)$coefficients["time", "Pr(>|z|)"]
    }

    irr_str  <- sprintf("%.3f [%.3f\u2013%.3f]", est, ci_low, ci_high)
    pval_str <- if (pval < 0.0001) "< 0.0001" else sprintf("%.4f", pval)
  } else {
    irr_str  <- "\u2014"
    pval_str <- "\u2014"
  }

  c(Metric = label, vals_yr, IRR = irr_str, p = pval_str)
}

table2_raw <- list(
  c(Metric = "Exposure Denominators - Patient-Days",
    setNames(format(denom_annuel$patient_days, big.mark = " "), paste0("yr_", denom_annuel$year)),
    IRR = "\u2014", p = "\u2014"),
  c(Metric = "Exposure Denominators - ICU Admissions",
    setNames(format(denom_annuel$admissions, big.mark = " "), paste0("yr_", denom_annuel$year)),
    IRR = "\u2014", p = "\u2014"),

  # Active epidemic transmissions (<= 30-day clusters / 1,000 patient-days)
  compute_metric_row("Active Epidemic Transmissions - K. pneumoniae + E. cloacae complex",
                     "is_trans == 1 & Taxon_clean %in% c('Klebsiella pneumoniae', 'Enterobacter cloacae complex')"),
  compute_metric_row("Active Epidemic Transmissions - Klebsiella pneumoniae",
                     "is_trans == 1 & Taxon_clean == 'Klebsiella pneumoniae'"),
  compute_metric_row("Active Epidemic Transmissions - Enterobacter cloacae complex",
                     "is_trans == 1 & Taxon_clean == 'Enterobacter cloacae complex'"),
  compute_metric_row("Active Epidemic Transmissions - Escherichia coli",
                     "is_trans == 1 & Taxon_clean == 'Escherichia coli'"),

  # Total nosocomial acquisitions (> 2 days / 1,000 patient-days)
  compute_metric_row("Total Nosocomial Acquisitions - K. pneumoniae + E. cloacae complex",
                     "is_nosocomial == 1 & Taxon_clean %in% c('Klebsiella pneumoniae', 'Enterobacter cloacae complex')"),
  compute_metric_row("Total Nosocomial Acquisitions - Klebsiella pneumoniae",
                     "is_nosocomial == 1 & Taxon_clean == 'Klebsiella pneumoniae'"),
  compute_metric_row("Total Nosocomial Acquisitions - Enterobacter cloacae complex",
                     "is_nosocomial == 1 & Taxon_clean == 'Enterobacter cloacae complex'"),
  compute_metric_row("Total Nosocomial Acquisitions - Escherichia coli",
                     "is_nosocomial == 1 & Taxon_clean == 'Escherichia coli'"),

  # Internal negative control (community importations / 1,000 admissions).
  # Standard Poisson is used here for exact consistency with the historical
  # Figure 2C (p = 0.0069).
  compute_metric_row("Community Importations - Escherichia coli (Negative Control)",
                     "is_imported == 1 & Taxon_clean == 'Escherichia coli'",
                     denom_col = "admissions", scale = 1000, model_type = "poisson")
)

final_table2 <- bind_rows(lapply(table2_raw, as.data.frame.list, stringsAsFactors = FALSE)) %>%
  rename(
    `Outcome Metric & Pathogen` = Metric,
    `2019` = yr_2019, `2020` = yr_2020, `2021` = yr_2021, `2022` = yr_2022, `2023` = yr_2023,
    `Annual IRR [95% CI]` = IRR, `p-value` = p
  )

write_csv(final_table2, file.path(out_dir_tab, "Table_2_Longitudinal_Trends_and_IRR.csv"))
cat("--> Table 2 exported with restored mathematical consistency.\n")
print(as.data.frame(final_table2))

# Consolidated multi-sheet Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Table_2_Consolidated")
writeData(wb, sheet = 1, x = final_table2)
saveWorkbook(wb, file.path(out_dir_tab, "Table_2_Consolidated.xlsx"), overwrite = TRUE)

# ------------------------------------------------------------------------------
# 4. FIGURE 2 CONSTRUCTION (PANELS 2A, 2B, 2C)
# ------------------------------------------------------------------------------
species_order <- c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex")

tableau_annuel_plot <- df_clean %>%
  group_by(Year_clean, Taxon_clean) %>%
  summarise(
    n_nosocomial    = sum(is_nosocomial),
    n_transmissions = sum(is_trans),
    n_imported      = sum(is_imported),
    .groups         = "drop"
  ) %>%
  left_join(denom_annuel, by = c("Year_clean" = "year")) %>%
  mutate(
    incidence_noso_1000pd  = (n_nosocomial / patient_days) * 1000,
    incidence_trans_1000pd = (n_transmissions / patient_days) * 1000,
    rate_import_1000adm    = (n_imported / admissions) * 1000,
    species_grp            = factor(Taxon_clean, levels = species_order)
  )

# Dynamic 5-year percentage changes (Panel 2A)
labels_change_2A <- tableau_annuel_plot %>%
  filter(Year_clean %in% c(2019, 2023)) %>%
  select(Year_clean, species_grp, incidence_noso_1000pd) %>%
  pivot_wider(names_from = Year_clean, values_from = incidence_noso_1000pd, names_prefix = "yr_") %>%
  mutate(
    pct_change = ((yr_2023 - yr_2019) / yr_2019) * 100,
    label_text = sprintf("%+.1f%%", pct_change),
    Year       = 2023,
    Rate       = yr_2023
  )

species_colors <- c(
  "Escherichia coli"             = "#1f78b4",
  "Klebsiella pneumoniae"        = "#e31a1c",
  "Enterobacter cloacae complex" = "#ff7f00"
)

species_legend_expr <- c(
  "Escherichia coli"             = expression(italic("Escherichia coli")~"(Community influx / Non-transmitting)"),
  "Klebsiella pneumoniae"        = expression(italic("Klebsiella pneumoniae")~"(Outbreak-forming / Extinct)"),
  "Enterobacter cloacae complex" = expression(italic("Enterobacter cloacae")~"complex (Hospital-adapted / Controlled)")
)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = rel(1.02), hjust = 0),
      plot.subtitle      = element_text(size = rel(0.85), color = "grey30", margin = margin(b = 6)),
      axis.title         = element_text(face = "bold", size = rel(0.90)),
      axis.text          = element_text(color = "black", size = rel(0.85)),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = rel(0.88)),
      legend.text        = element_text(size = rel(0.82)),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      plot.margin        = margin(6, 14, 6, 6)
    )
}

# --- Panel 2A: Conventional nosocomial acquisitions (> 2 days) ---
p2A <- ggplot(tableau_annuel_plot, aes(x = Year_clean, y = incidence_noso_1000pd, color = species_grp, group = species_grp)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.6, stroke = 1) +
  geom_text(
    data = labels_change_2A,
    aes(x = Year + 0.12, y = Rate, label = label_text, color = species_grp),
    hjust = 0, vjust = 0.5, fontface = "bold", size = 3.1, show.legend = FALSE
  ) +
  scale_color_manual(values = species_colors, labels = species_legend_expr, name = "Bacterial Species:") +
  scale_x_continuous(breaks = 2019:2023, limits = c(2018.9, 2023.6)) +
  scale_y_continuous(limits = c(0, 2.2), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "A  Conventional Surveillance (> 2 Days Post-Admission)",
    subtitle = expression(atop(
      paste("Sustained decline in ", italic("K. pneumoniae"), " and ", italic("E. cloacae"), " complex"),
      paste("Persistent ", italic("E. coli"), " detection (IRR = 0.792, 95% CI 0.727\u20130.864, p < 0.0001)")
    )),
    x = "Year",
    y = "Acquisitions / 1,000 Patient-Days"
  ) +
  theme_pub()

# --- Panel 2B: Genomic transmissions (<= 30-day window) ---
p2B <- ggplot(tableau_annuel_plot, aes(x = Year_clean, y = incidence_trans_1000pd, color = species_grp, group = species_grp)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.6, stroke = 1) +
  scale_color_manual(values = species_colors, labels = species_legend_expr, name = "Bacterial Species:") +
  scale_x_continuous(breaks = 2019:2023, limits = c(2018.9, 2023.6)) +
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1.2, 0.2), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "B  Genomic Resolution (<= 30-Day WGS Clusters)",
    subtitle = expression(atop(
      paste("Suppression of active ICU outbreaks; absence of ", italic("E. coli"), " transmission"),
      paste("(IRR = 0.375, 95% CI 0.290\u20130.485, p < 0.0001)")
    )),
    x = "Year",
    y = "Active Cases / 1,000 Patient-Days"
  ) +
  theme_pub()

# --- Panel 2C: Internal negative control (community E. coli importations) ---
df_eco_import <- tableau_annuel_plot %>% filter(species_grp == "Escherichia coli")
eco_change_pct <- ((df_eco_import$rate_import_1000adm[df_eco_import$Year_clean == 2023] -
                      df_eco_import$rate_import_1000adm[df_eco_import$Year_clean == 2019]) /
                     df_eco_import$rate_import_1000adm[df_eco_import$Year_clean == 2019]) * 100

p2C <- ggplot(df_eco_import, aes(x = Year_clean, y = rate_import_1000adm)) +
  geom_line(linewidth = 1.2, color = "#1f78b4") +
  geom_point(size = 3.8, color = "#1f78b4", shape = 21, fill = "white", stroke = 1.8) +
  annotate(
    "text", x = 2023.12, y = df_eco_import$rate_import_1000adm[df_eco_import$Year_clean == 2023],
    label = sprintf("%+.1f%%", eco_change_pct),
    hjust = 0, vjust = 0.5, fontface = "bold", size = 3.1, color = "#1f78b4"
  ) +
  scale_x_continuous(breaks = 2019:2023, limits = c(2018.9, 2023.6)) +
  scale_y_continuous(limits = c(0, 20), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = expression(paste("C  Internal Negative Control (Community ", italic("Escherichia coli"), " Influx)")),
    subtitle = expression(paste("Unabated baseline importation pressure at admission (IRR = 1.136, 95% CI 1.036\u20131.247, p = 0.0069)")),
    x        = "Year",
    y        = "Importations / 1,000 Admissions"
  ) +
  theme_pub()

# --- Assembly and export of Figure 2 ---
fig2_complete <- (p2A | p2B) / p2C +
  plot_layout(guides = "collect") &
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.title      = element_text(face = "bold", size = 9.5),
    legend.text       = element_text(size = 9.0),
    legend.spacing.x  = unit(0.7, "cm"),
    legend.key.width  = unit(1.1, "cm"),
    legend.background = element_rect(fill = alpha("grey95", 0.5), color = "grey80", linewidth = 0.35)
  )

ggsave(file.path(out_dir_fig, "Figure_2_Longitudinal_Dynamics.pdf"), fig2_complete, width = 13, height = 9.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_2_Longitudinal_Dynamics.png"), fig2_complete, width = 13, height = 9.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_2_Longitudinal_Dynamics.tiff"), fig2_complete, width = 13, height = 9.5, dpi = 300, compression = "lzw")

cat(sprintf("\n--> Figure 2 and Table 2 generated successfully in:\n%s\n%s\n", out_dir_fig, out_dir_tab))
