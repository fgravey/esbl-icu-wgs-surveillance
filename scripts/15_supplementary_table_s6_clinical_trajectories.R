# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 15_supplementary_table_s6_clinical_trajectories.R
# Purpose: Dynamic analysis of individual clinical trajectories and the
#          "Surgical Paradox" (isolated infection with negative screening),
#          with export of the trajectory distribution underlying
#          Supplementary Table S6 and Figure 4.
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/tables/Supplementary_Table_S6_Clinical_Trajectories.csv
#          output/figures/Figure_4_Clinical_Trajectories_Paradox.pdf (.png)
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
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
# 1. DYNAMIC CLINICAL HARMONIZATION OF TRAJECTORIES (N = 689)
# ------------------------------------------------------------------------------
# NOTE: the chronological screening status uses French labels in the source
# database ("toujours clinique" = always clinical, "apres"/"apr\u00e8s" = after,
# "avant" = before, "meme jour"/"m\u00eame jour" = same day). These literal
# strings are matched below and must not be translated.
df_prep <- df %>%
  mutate(
    ward_clean = case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med", cores, ignore.case = TRUE)  ~ "Medical ICU",
      TRUE ~ NA_character_
    ),
    inf_count  = as.numeric(as.character(Infection.BLSE)),
    has_inf    = ifelse(!is.na(inf_count) & inf_count > 0, 1, 0),
    chrono_raw = tolower(trimws(as.character(depistage.avant.apres))),

    # Dynamic classification of the 5 clinical trajectories
    trajectory = case_when(
      chrono_raw == "toujours clinique" ~ "4. Isolated Infection ('Surgical Paradox')",
      chrono_raw %in% c("apres", "apr\u00e8s") ~ "2. True Conversion (Infection after screen)",
      chrono_raw == "avant" ~ "3. Acute Inoculation (Infection before screen)",
      chrono_raw %in% c("meme jour", "m\u00eame jour") ~ "5. Same-Day Isolation",
      has_inf == 0 ~ "1. Silent Colonization",
      TRUE ~ "1. Silent Colonization"
    ),
    trajectory = factor(trajectory, levels = rev(c(
      "1. Silent Colonization",
      "2. True Conversion (Infection after screen)",
      "3. Acute Inoculation (Infection before screen)",
      "4. Isolated Infection ('Surgical Paradox')",
      "5. Same-Day Isolation"
    )))
  ) %>%
  filter(!is.na(ward_clean))

cat(sprintf("--> Cohort available for trajectory analysis: %d episodes\n", nrow(df_prep)))

# ------------------------------------------------------------------------------
# 2. PANEL 4A: DYNAMIC DISTRIBUTION OF CLINICAL TRAJECTORIES
# ------------------------------------------------------------------------------
# Strict per-ward denominators
denom_ward <- df_prep %>%
  count(ward_clean, name = "N_ward_total")

traj_summary <- df_prep %>%
  count(ward_clean, trajectory, name = "Count") %>%
  complete(ward_clean, trajectory, fill = list(Count = 0)) %>%
  left_join(denom_ward, by = "ward_clean") %>%
  mutate(
    Percentage = (Count / N_ward_total) * 100,
    Unit_Label = sprintf("%s (n=%d)", ward_clean, N_ward_total)
  )

p4A <- ggplot(traj_summary, aes(x = Percentage, y = trajectory, fill = ward_clean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(
    aes(label = sprintf("%.1f%% (%d)", Percentage, Count)),
    position = position_dodge(width = 0.8),
    hjust = -0.15, size = 3.0, fontface = "plain"
  ) +
  scale_fill_manual(values = c("Surgical ICU" = "#e41a1c", "Medical ICU" = "#377eb8"), name = "Unit:") +
  scale_x_continuous(limits = c(0, 95), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "A  Chronological Carriage-to-Infection Trajectories",
    subtitle = sprintf("Dynamic distribution across %d ESBL-E episodes by intensive care unit", nrow(df_prep)),
    x        = "Proportion of Unit Episodes (%)",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = rel(1.05)),
    axis.text          = element_text(color = "black"),
    panel.grid.major.y = element_blank()
  )

# ------------------------------------------------------------------------------
# 3. PANEL 4B: STATISTICAL FOCUS ON THE SURGICAL PARADOX
# ------------------------------------------------------------------------------
df_paradox <- df_prep %>%
  filter(trajectory == "4. Isolated Infection ('Surgical Paradox')")

total_paradox <- nrow(df_paradox)

paradox_summary <- df_paradox %>%
  count(ward_clean, name = "Count") %>%
  mutate(
    Percentage = (Count / total_paradox) * 100,
    Unit_Label = sprintf("%s\n(%.1f%%, n=%d)", ward_clean, Percentage, Count)
  )

# Over-representation test (conformity chi-square on the 2x2 table)
tab_paradox_2x2 <- table(df_prep$ward_clean, df_prep$trajectory == "4. Isolated Infection ('Surgical Paradox')")
p_val_paradox   <- chisq.test(tab_paradox_2x2)$p.value
p_str_paradox   <- ifelse(p_val_paradox < 0.001, "p < 0.001", sprintf("p = %.4f", p_val_paradox))

max_y_b <- max(paradox_summary$Count, na.rm = TRUE)

p4B <- ggplot(paradox_summary, aes(x = ward_clean, y = Count, fill = ward_clean)) +
  geom_bar(stat = "identity", width = 0.55, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%d cases\n(%.1f%%)", Count, Percentage)),
    vjust = -0.3, size = 3.2, fontface = "bold"
  ) +
  scale_fill_manual(values = c("Surgical ICU" = "#e41a1c", "Medical ICU" = "#377eb8")) +
  scale_y_continuous(limits = c(0, max_y_b * 1.3), expand = expansion(mult = c(0, 0.05))) +
  annotate(
    "label", x = 1.5, y = max_y_b * 1.15,
    label = sprintf("Chi-square test:\n%s", p_str_paradox),
    size = 3.0, fontface = "bold", fill = alpha("white", 0.9)
  ) +
  labs(
    title    = sprintf("B  The 'Surgical Paradox' Distribution (n = %d)", total_paradox),
    subtitle = "Severe infections with negative rectal screening",
    x        = NULL,
    y        = "Number of Patients"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = rel(1.05)),
    axis.text  = element_text(color = "black")
  )

# ------------------------------------------------------------------------------
# 4. PANEL 4C: DOCUMENTED CLINICAL INFECTION SITES (TEXT MINING)
# ------------------------------------------------------------------------------
# NOTE: site classification keywords match French free-text entries of the
# source database (e.g. "poumon" = lung, "urin" = urinary).
df_sites_raw <- df_prep %>%
  filter(has_inf == 1) %>%
  mutate(nature_brute = tolower(trimws(as.character(infection.nature)))) %>%
  filter(!nature_brute %in% c("-", "", "na", "nr", "unknown")) %>%
  separate_rows(nature_brute, sep = "[,;/]") %>%
  mutate(nature_brute = trimws(nature_brute)) %>%
  filter(nature_brute != "") %>%
  mutate(
    Site = case_when(
      str_detect(nature_brute, "urin")                              ~ "Urinary Tract (UTI)",
      str_detect(nature_brute, "poumon|pneumo|asp|lba|trache|pavm") ~ "Ventilator Pneumonia (VAP)",
      str_detect(nature_brute, "bacterie|hemoc|sepsis")             ~ "Bloodstream / Sepsis (BSI)",
      str_detect(nature_brute, "catheter|kt")                       ~ "Catheter-Related (CRI)",
      str_detect(nature_brute, "plaie|pus|paroi|os|artic|peau")     ~ "Skin, Soft Tissue & Bone",
      str_detect(nature_brute, "abdom|adbom|periton|bile")          ~ "Intra-Abdominal (IAI)",
      TRUE                                                          ~ "Other Sites"
    )
  )

total_sites <- nrow(df_sites_raw)

sites_summary <- df_sites_raw %>%
  count(Site, name = "Count") %>%
  mutate(
    Percentage = (Count / total_sites) * 100,
    Site = factor(Site, levels = rev(c(
      "Urinary Tract (UTI)",
      "Ventilator Pneumonia (VAP)",
      "Bloodstream / Sepsis (BSI)",
      "Catheter-Related (CRI)",
      "Skin, Soft Tissue & Bone",
      "Intra-Abdominal (IAI)",
      "Other Sites"
    )))
  ) %>%
  filter(!is.na(Site))

max_x_c <- max(sites_summary$Percentage, na.rm = TRUE)

p4C <- ggplot(sites_summary, aes(x = Percentage, y = Site)) +
  geom_bar(stat = "identity", fill = "#4daf4a", width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%% (n=%d)", Percentage, Count)),
    hjust = -0.1, size = 3.0, fontface = "plain"
  ) +
  scale_x_continuous(limits = c(0, max_x_c * 1.3), labels = function(x) paste0(x, "%")) +
  labs(
    title    = sprintf("C  Clinical Infection Sites (n = %d)", total_sites),
    subtitle = "Anatomical breakdown of all confirmed infections",
    x        = "Proportion of Clinical Infections (%)",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = rel(1.05)),
    axis.text          = element_text(color = "black"),
    panel.grid.major.y = element_blank()
  )

# ------------------------------------------------------------------------------
# 5. ASSEMBLY OF FIGURE 4 AND EXPORTS
# ------------------------------------------------------------------------------
fig4_complete <- p4A / (p4B | p4C) +
  plot_layout(heights = c(1.15, 1))

ggsave(file.path(out_dir_fig, "Figure_4_Clinical_Trajectories_Paradox.pdf"),
       fig4_complete, width = 12.0, height = 9.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_4_Clinical_Trajectories_Paradox.png"),
       fig4_complete, width = 12.0, height = 9.5, dpi = 300)

# Tabular export of Table S6
table_s6 <- traj_summary %>%
  select(Ward = ward_clean, Trajectory = trajectory, Count, Percentage) %>%
  arrange(Ward, desc(Count))

write_csv(table_s6, file.path(out_dir_tab, "Supplementary_Table_S6_Clinical_Trajectories.csv"))

cat(sprintf("\n--> Figure 4 saved to: %s\n", out_dir_fig))
cat(sprintf("--> Table S6 saved to: %s\n", out_dir_tab))
