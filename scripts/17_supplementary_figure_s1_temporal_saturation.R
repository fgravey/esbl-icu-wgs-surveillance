# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 17_supplementary_figure_s1_temporal_saturation.R
# Purpose: Dynamic non-linear temporal saturation modelling (NLS) and
#          multi-window sensitivity analysis justifying the <= 30-day
#          transmission window:
#          - Panel S1A: Empirical saturation curve of genomic linkage
#            (logistic NLS model across inter-isolate delays)
#          - Panel S1B: Sensitivity of transmission assignments across
#            14-, 21-, 30- and 60-day rolling windows
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Supplementary_Figure_S1_Temporal_Saturation.pdf (.png)
#          output/tables/Supplementary_Table_S1_Sensitivity.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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

# ------------------------------------------------------------------------------
# 1. DATA LOADING AND DYNAMIC DENOMINATORS
# ------------------------------------------------------------------------------
df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# Dynamic species denominators extracted directly from the cleaned cohort
species_denominators <- df %>% count(Taxon_group, name = "Denom_Total")

cat("--> Exact species denominators in cleaned cohort:\n")
print(species_denominators)

# ------------------------------------------------------------------------------
# 2. DYNAMIC NEAREST TEMPORAL NEIGHBOUR COMPUTATION
# ------------------------------------------------------------------------------
# Focus on valid genomic clusters for nosocomial opportunistic clades
# (K. pneumoniae and E. cloacae complex)
df_clustered <- df %>%
  filter(
    is_valid_cluster == TRUE,
    Taxon_group %in% c("Klebsiella pneumoniae", "Enterobacter cloacae complex"),
    !is.na(date_sampling)
  )

# Minimal pairwise delay to a co-hospitalised intra-cluster predecessor
df_nn <- df_clustered %>%
  group_by(ward_clean, snp_cluster_clean, Taxon_group) %>%
  group_modify(function(dat, key) {
    if (nrow(dat) < 2) {
      dat$min_delay <- NA_real_
      return(dat)
    }

    deltas <- sapply(seq_len(nrow(dat)), function(i) {
      cur_pat  <- dat$Correspondant[i]
      cur_date <- as.Date(dat$date_sampling[i])

      other_dates <- as.Date(dat$date_sampling[dat$Correspondant != cur_pat & !is.na(dat$date_sampling)])
      if (length(other_dates) == 0) return(NA_real_)

      min(abs(as.numeric(difftime(cur_date, other_dates, units = "days"))))
    })

    dat$min_delay <- deltas
    dat
  }) %>%
  ungroup() %>%
  filter(!is.na(min_delay))

cat(sprintf("\nDocumented nosocomial isolates with valid linkage delay: %d\n", nrow(df_nn)))

# ------------------------------------------------------------------------------
# 3. NON-LINEAR LEAST SQUARES (NLS) KINETIC MODELLING (0-60 DAYS)
# ------------------------------------------------------------------------------
n_acute_60d <- sum(df_nn$min_delay <= 60)

# Daily empirical cumulative series
df_sat_daily <- tibble(t = 1:60) %>%
  rowwise() %>%
  mutate(
    n_cum   = sum(df_nn$min_delay <= t),
    pct_obs = (n_cum / n_acute_60d) * 100
  ) %>%
  ungroup()

# Dynamic logistic model: pct_obs = L / (1 + exp(-k * (t - x0)))
fit_nls <- nls(
  pct_obs ~ L / (1 + exp(-k * (t - x0))),
  data = df_sat_daily,
  start = list(L = max(df_sat_daily$pct_obs), k = 0.15, x0 = 9)
)

params_fit <- coef(fit_nls)
L_est  <- as.numeric(params_fit["L"])
k_est  <- as.numeric(params_fit["k"])
x0_est <- as.numeric(params_fit["x0"])

val_30d          <- L_est / (1 + exp(-k_est * (30 - x0_est)))
sat_capacity_30d <- (val_30d / L_est) * 100

cat("\n=======================================================\n")
cat("DYNAMICALLY ESTIMATED NLS PARAMETERS (ACUTE CAPACITY):\n")
cat(sprintf("\u2022 Acute asymptotic plateau (L)      : %.2f %%\n", L_est))
cat(sprintf("\u2022 Propagation rate constant (k)     : %.4f\n", k_est))
cat(sprintf("\u2022 Half-saturation inflection (x0)   : %.2f days\n", x0_est))
cat(sprintf("\u2022 Absolute modelled level at 30 days: %.2f %%\n", val_30d))
cat(sprintf("\u2022 Acute capacity captured at 30 days: %.2f %%\n", sat_capacity_30d))
cat("=======================================================\n\n")

# ------------------------------------------------------------------------------
# 4. PANEL S1A: DYNAMIC EMPIRICAL SATURATION CURVE
# ------------------------------------------------------------------------------
df_curve <- tibble(
  Interval_Days = seq(0, 100, by = 0.5),
  Pct_Modeled   = L_est / (1 + exp(-k_est * (Interval_Days - x0_est)))
)

p_s1A <- ggplot(df_curve, aes(x = Interval_Days, y = Pct_Modeled)) +
  annotate("rect", xmin = 30, xmax = 100, ymin = 0, ymax = 100, fill = "grey92", alpha = 0.6) +
  geom_point(data = df_sat_daily, aes(x = t, y = pct_obs), color = "grey40", size = 1.2, alpha = 0.7) +
  geom_line(color = "#00468B", linewidth = 1.2) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "#ED0000", linewidth = 0.8) +
  geom_point(aes(x = 30, y = val_30d), color = "#ED0000", size = 3.5) +
  annotate(
    "label", x = 33, y = 35, hjust = 0, size = 3.0, fill = alpha("white", 0.95),
    label = sprintf(
      "Inflection Point:\n\u2022 Temporal cutoff: <= 30 days\n\u2022 Acute capacity captured: %.1f%%\n\u2022 Predicted cumulative linkage at day 30: %.1f%%\n\u2022 Censored outliers (>30 days)",
      sat_capacity_30d, val_30d
    )
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 102), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "A  Empirical Saturation of Genomic Linkage",
    subtitle = "Non-linear saturation kinetics (NLS) across inter-isolate delays",
    x        = "Inter-Isolate Isolation Interval (Days)",
    y        = "Cumulative Linkage (%)"
  ) +
  theme_classic(base_size = 10)

# ------------------------------------------------------------------------------
# 5. PANEL S1B: MULTI-WINDOW SENSITIVITY ANALYSIS (DYNAMIC RECALCULATION)
# ------------------------------------------------------------------------------
windows <- c(14, 21, 30, 60)

calc_window_cases <- function(w) {
  df %>%
    filter(is_valid_cluster == TRUE) %>%
    group_by(ward_clean, snp_cluster_clean, Taxon_group) %>%
    group_modify(function(dat, key) {
      if (nrow(dat) < 2) {
        dat$is_sec <- FALSE
        return(dat)
      }
      dat$is_sec <- sapply(seq_len(nrow(dat)), function(i) {
        p_i <- dat$Correspondant[i]
        d_i <- as.Date(dat$date_sampling[i])
        others <- as.Date(dat$date_sampling[dat$Correspondant != p_i & !is.na(dat$date_sampling)])
        if (length(others) == 0) return(FALSE)
        min(abs(as.numeric(difftime(d_i, others, units = "days")))) <= w
      })
      dat
    }) %>%
    ungroup() %>%
    group_by(Taxon_group) %>%
    summarise(Secondary_Cases = sum(is_sec), .groups = "drop") %>%
    mutate(Window = w)
}

sens_table <- bind_rows(lapply(windows, calc_window_cases)) %>%
  left_join(species_denominators, by = "Taxon_group") %>%
  mutate(
    Pct = 100 * Secondary_Cases / Denom_Total,
    Window_Label = factor(paste0(Window, " days", ifelse(Window == 30, " (Selected)", "")),
                          levels = c("14 days", "21 days", "30 days (Selected)", "60 days")),
    Taxon_group = factor(Taxon_group, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  )

cat("--> Sensitivity analysis across rolling windows:\n")
print(sens_table %>% select(Window_Label, Taxon_group, Secondary_Cases, Denom_Total, Pct))

species_colors <- c(
  "Escherichia coli"             = "#1f78b4",
  "Klebsiella pneumoniae"        = "#e31a1c",
  "Enterobacter cloacae complex" = "#ff7f00"
)

p_s1B <- ggplot(sens_table, aes(x = Window_Label, y = Pct, fill = Taxon_group)) +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(n=%d)", Pct, Secondary_Cases)),
    position = position_dodge(0.75),
    vjust = -0.2, size = 2.6, fontface = "bold"
  ) +
  scale_fill_manual(
    values = species_colors,
    labels = c(expression(italic("E. coli")),
               expression(italic("K. pneumoniae")),
               expression(italic("E. cloacae")~"complex")),
    name = "Species:"
  ) +
  scale_y_continuous(limits = c(0, max(sens_table$Pct) * 1.25), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "B  Sensitivity Analysis Across Rolling Temporal Windows",
    subtitle = "Stability of transmission assignments across 14-, 21-, 30-, and 60-day cutoffs",
    x        = "Rolling Time Window within Same Ward",
    y        = "Confirmed Secondary Acquisitions (%)"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.justification = "right")

# ------------------------------------------------------------------------------
# 6. COMPOSITION AND EXPORT
# ------------------------------------------------------------------------------
fig_s1 <- p_s1A / p_s1B + plot_layout(heights = c(1, 1))

ggsave(file.path(out_dir_fig, "Supplementary_Figure_S1_Temporal_Saturation.pdf"), fig_s1, width = 8.5, height = 9.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Supplementary_Figure_S1_Temporal_Saturation.png"), fig_s1, width = 8.5, height = 9.5, dpi = 300)
write.csv(sens_table, file.path(out_dir_tab, "Supplementary_Table_S1_Sensitivity.csv"), row.names = FALSE)

cat(sprintf("\n--> Supplementary Figure S1 exported to: %s\n", out_dir_fig))
cat(sprintf("--> Supplementary Table S1 exported to: %s\n", out_dir_tab))
