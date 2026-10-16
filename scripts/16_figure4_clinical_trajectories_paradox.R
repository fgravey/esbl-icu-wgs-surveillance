# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 16_figure4_clinical_trajectories_paradox.R
# Purpose: Dynamic generation of Figure 4 (The Lancet Microbe):
#          - 4A: Chronological Carriage-to-Infection Trajectories by Unit
#          - 4B: Genomic Dissection of the Surgical Paradox (SICU)
#          - 4C: Clinical Infection Anatomical Breakdown (Text-mining)
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Figure_4_Clinical_Trajectories_Paradox.pdf (.png, .tiff)
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
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
df <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. CLINICAL CLEANING
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

    trajectory = case_when(
      chrono_raw == "toujours clinique"           ~ "4. Isolated Infection ('Surgical Paradox')",
      chrono_raw %in% c("apres", "apr\u00e8s")         ~ "2. True Conversion (Infection after screen)",
      chrono_raw == "avant"                       ~ "3. Acute Inoculation (Infection before screen)",
      chrono_raw %in% c("meme jour", "m\u00eame jour") ~ "5. Same-Day Isolation",
      TRUE                                        ~ "1. Silent Colonization"
    ),
    trajectory = factor(trajectory, levels = rev(c(
      "1. Silent Colonization",
      "2. True Conversion (Infection after screen)",
      "3. Acute Inoculation (Infection before screen)",
      "4. Isolated Infection ('Surgical Paradox')",
      "5. Same-Day Isolation"
    ))),

    Species_clean = case_when(
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE) ~ "Klebsiella pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)    ~ "Enterobacter cloacae complex",
      TRUE                                                 ~ "Escherichia coli"
    ),
    ST_clean = str_trim(as.character(ST)),
    ST_clean = ifelse(!grepl("^ST", ST_clean) & !ST_clean %in% c("Unknown", "ND", "", "NA"), paste0("ST", ST_clean), ST_clean)
  ) %>%
  filter(!is.na(ward_clean))

# ------------------------------------------------------------------------------
# 1. PANEL 4A: CLINICAL TRAJECTORIES BY UNIT
# ------------------------------------------------------------------------------
denom_ward <- df_prep %>% count(ward_clean, name = "N_ward_total")

traj_summary <- df_prep %>%
  count(ward_clean, trajectory, name = "Count") %>%
  complete(ward_clean, trajectory, fill = list(Count = 0)) %>%
  left_join(denom_ward, by = "ward_clean") %>%
  mutate(
    Percentage = (Count / N_ward_total) * 100,
    Unit_Label = sprintf("%s (n = %d)", ward_clean, N_ward_total)
  )

p4A <- ggplot(traj_summary, aes(x = Percentage, y = trajectory, fill = Unit_Label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%% (%d)", Percentage, Count)),
    position = position_dodge(width = 0.75),
    hjust = -0.1, size = 2.9, fontface = "bold"
  ) +
  scale_fill_manual(values = setNames(c("#e41a1c", "#377eb8"), unique(traj_summary$Unit_Label))) +
  scale_x_continuous(limits = c(0, 95), breaks = seq(0, 80, 20), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "A  Chronological Carriage-to-Infection Trajectories",
    subtitle = sprintf("Distribution across %d documented ESBL-E episodes by intensive care unit", nrow(df_prep)),
    x        = "Proportion of Unit Episodes (%)",
    y        = NULL,
    fill     = "Unit:"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = rel(1.02)),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35)
  )

# ------------------------------------------------------------------------------
# 2. PANEL 4B: GENOMIC DISSECTION OF THE SURGICAL PARADOX (SICU)
# ------------------------------------------------------------------------------
df_sicu_paradox <- df_prep %>%
  filter(ward_clean == "Surgical ICU" & trajectory == "4. Isolated Infection ('Surgical Paradox')")

total_sicu_paradox <- nrow(df_sicu_paradox)

# Dominant clone per species
dominant_clones <- list(
  "Enterobacter cloacae complex" = "ST66",
  "Escherichia coli"             = "ST131",
  "Klebsiella pneumoniae"        = "ST405"
)

paradox_genomics <- df_sicu_paradox %>%
  group_by(Species_clean) %>%
  summarise(
    Total_Species = n(),
    Dom_ST        = dominant_clones[[unique(Species_clean)]],
    N_Dom         = sum(ST_clean == Dom_ST),
    N_Other       = Total_Species - N_Dom,
    Pct_Total     = (Total_Species / total_sicu_paradox) * 100,
    Pct_Dom       = (N_Dom / total_sicu_paradox) * 100,
    Pct_Other     = (N_Other / total_sicu_paradox) * 100,
    .groups       = "drop"
  )

df_4B_segments <- paradox_genomics %>%
  pivot_longer(cols = c(Pct_Dom, Pct_Other), names_to = "SegType", values_to = "Pct") %>%
  mutate(
    Species   = factor(Species_clean, levels = c("Klebsiella pneumoniae", "Escherichia coli", "Enterobacter cloacae complex")),
    FillGroup = case_when(
      Species_clean == "Enterobacter cloacae complex" & SegType == "Pct_Dom"   ~ "Ecc_dom",
      Species_clean == "Enterobacter cloacae complex" & SegType == "Pct_Other" ~ "Ecc_oth",
      Species_clean == "Escherichia coli" & SegType == "Pct_Dom"               ~ "Eco_dom",
      Species_clean == "Escherichia coli" & SegType == "Pct_Other"             ~ "Eco_oth",
      Species_clean == "Klebsiella pneumoniae" & SegType == "Pct_Dom"          ~ "Kpn_dom",
      TRUE                                                                    ~ "Kpn_oth"
    )
  )

df_4B_labels <- paradox_genomics %>%
  mutate(
    Species = factor(Species_clean, levels = c("Klebsiella pneumoniae", "Escherichia coli", "Enterobacter cloacae complex")),
    Text_Desc = sprintf("Dominant clone: %s (n = %d, %.1f%%)  |  other STs (n = %d)", Dom_ST, N_Dom, (N_Dom/total_sicu_paradox)*100, N_Other)
  )

species_labels_expr <- c(
  "Enterobacter cloacae complex" = expression(italic("E. cloacae")~"complex"),
  "Escherichia coli"             = expression(italic("E. coli")),
  "Klebsiella pneumoniae"        = expression(italic("K. pneumoniae"))
)

palette_4B <- c(
  "Ecc_dom" = "#ff7f00", "Ecc_oth" = "#fdbf6f",
  "Eco_dom" = "#1f78b4", "Eco_oth" = "#a6cee3",
  "Kpn_dom" = "#e31a1c", "Kpn_oth" = "#fb9a99"
)

p4B <- ggplot() +
  geom_col(data = df_4B_segments, aes(x = Pct, y = Species, fill = FillGroup), position = "stack", width = 0.50, show.legend = FALSE) +
  geom_text(
    data = df_4B_labels,
    aes(x = Pct_Total + 1.2, y = Species, label = sprintf("%.1f%% (n = %d/%d)", Pct_Total, Total_Species, total_sicu_paradox)),
    hjust = 0, vjust = 0.5, fontface = "bold", size = 3.1, color = "black"
  ) +
  geom_text(
    data = df_4B_labels,
    aes(x = 0.5, y = as.numeric(Species) + 0.35, label = Text_Desc),
    hjust = 0, fontface = "italic", size = 2.9, color = "grey20"
  ) +
  scale_fill_manual(values = palette_4B) +
  scale_y_discrete(labels = species_labels_expr, expand = expansion(mult = c(0.15, 0.25))) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 15), labels = function(x) paste0(x, "%")) +
  labs(
    title    = sprintf("B  Genomic Dissection of the Surgical Paradox (SICU, n = %d)", total_sicu_paradox),
    subtitle = "Dark shading: dominant epidemic clone (ST >= 4 cases); Light shading: other STs",
    x        = "Proportion of Isolated Surgical Infections (%)",
    y        = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title  = element_text(face = "bold", size = rel(1.02)),
    axis.text.y = element_text(face = "bold", size = rel(0.92))
  )

# ------------------------------------------------------------------------------
# 3. PANEL 4C: ANATOMICAL INFECTION SITES (TEXT MINING)
# ------------------------------------------------------------------------------
# NOTE: site classification keywords match French free-text entries of the
# source database (e.g. "poumon" = lung, "urin" = urinary).
df_sites <- df_prep %>%
  filter(has_inf == 1) %>%
  mutate(nature_raw = tolower(trimws(as.character(infection.nature)))) %>%
  filter(!nature_raw %in% c("-", "", "na", "nr", "unknown")) %>%
  separate_rows(nature_raw, sep = "[,;/]") %>%
  mutate(nature_raw = trimws(nature_raw)) %>%
  filter(nature_raw != "") %>%
  mutate(
    Site = case_when(
      str_detect(nature_raw, "urin")                              ~ "Urinary Tract (UTI)",
      str_detect(nature_raw, "poumon|pneumo|asp|lba|trache|pavm") ~ "Ventilator Pneumonia (VAP)",
      str_detect(nature_raw, "bacterie|hemoc|sepsis")             ~ "Bloodstream/Sepsis (BSI)",
      str_detect(nature_raw, "catheter|kt")                       ~ "Catheter-Related (CRI)",
      str_detect(nature_raw, "plaie|pus|paroi|os|artic|peau")     ~ "Skin, Soft Tissue & Bone",
      str_detect(nature_raw, "abdom|adbom|periton|bile")          ~ "Intra-Abdominal (IAI)",
      TRUE                                                        ~ "Other (CNS & Deep Pus)"
    )
  )

total_sites_count <- nrow(df_sites)

site_summary <- df_sites %>%
  count(Site, name = "N") %>%
  mutate(
    Pct = (N / total_sites_count) * 100,
    Site = factor(Site, levels = rev(c(
      "Urinary Tract (UTI)",
      "Ventilator Pneumonia (VAP)",
      "Bloodstream/Sepsis (BSI)",
      "Catheter-Related (CRI)",
      "Skin, Soft Tissue & Bone",
      "Intra-Abdominal (IAI)",
      "Other (CNS & Deep Pus)"
    )))
  ) %>%
  filter(!is.na(Site))

p4C <- ggplot(site_summary, aes(x = Pct, y = Site)) +
  geom_col(fill = "#41b6c4", width = 0.58) +
  geom_text(
    aes(label = sprintf("%.1f%% (n = %d)", Pct, N)),
    hjust = -0.1, size = 3.0, fontface = "bold", color = "black"
  ) +
  scale_x_continuous(limits = c(0, max(site_summary$Pct) * 1.3), labels = function(x) paste0(x, "%")) +
  labs(
    title    = sprintf("C  Clinical Infection Sites (n = %d)", total_sites_count),
    subtitle = "Anatomical breakdown of all confirmed clinical infections",
    x        = "Proportion of Clinical Infections (%)",
    y        = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title         = element_text(face = "bold", size = rel(1.02)),
    axis.text.y        = element_text(size = rel(0.85)),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35)
  )

# ------------------------------------------------------------------------------
# 4. ASSEMBLY AND EXPORT
# ------------------------------------------------------------------------------
fig4_complete <- p4A / (p4B | p4C) +
  plot_layout(heights = c(1, 1.15), widths = c(1, 1))

ggsave(file.path(out_dir_fig, "Figure_4_Clinical_Trajectories_Paradox.pdf"), fig4_complete, width = 13.5, height = 10, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_4_Clinical_Trajectories_Paradox.png"), fig4_complete, width = 13.5, height = 10, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_4_Clinical_Trajectories_Paradox.tiff"), fig4_complete, width = 13.5, height = 10, dpi = 300, compression = "lzw")

cat("--> Figure 4 successfully regenerated in:", out_dir_fig, "\n")
