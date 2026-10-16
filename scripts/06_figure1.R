# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 06_figure1.R
# Purpose: Dynamic creation of Figure 1 (three panels):
#          - Panel 1A: Acquisition origin (community-imported <= 2 days vs
#                      nosocomial > 2 days) by species
#          - Panel 1B: Tripartite resistome block (blaCTX-M-15 + blaOXA-1 +
#                      blaTEM-1) and phenotypic piperacillin-tazobactam
#                      resistance
#          - Panel 1C: Clonal population architecture and sequence types (ST)
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/figures/Figure_1_Species_Dichotomy_Resistome.pdf (.png, .tiff)
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

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_fig <- file.path(base_dir, "output", "figures")
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. CLEANING AND STANDARDISATION OF THE DE-DUPLICATED COHORT
# ------------------------------------------------------------------------------
df_clean <- df %>%
  mutate(
    # Species harmonisation
    Taxon_clean = case_when(
      grepl("coli", Taxon_group, ignore.case = TRUE)        ~ "Escherichia coli",
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE)  ~ "Klebsiella pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)     ~ "Enterobacter cloacae complex",
      TRUE                                                  ~ NA_character_
    ),
    # Days in hospital before first ESBL-E sampling
    communoso_num = as.numeric(as.character(communoso)),
    los_days      = pmax(0, ifelse(is.na(communoso_num), 0, communoso_num)),
    is_imported   = ifelse(los_days <= 2, 1, 0),
    is_nosocomial = ifelse(los_days > 2, 1, 0),

    # Phenotypic piperacillin-tazobactam resistance (R)
    tzp_r = toupper(trimws(as.character(piperacilline_tazobactam_Inter))) == "R",

    # AMRFinderPlus-based resistome detection (logic validated in script 05)
    amr_beta     = ifelse(is.na(AMR_BETA_LACTAM), "", as.character(AMR_BETA_LACTAM)),
    beta_nat     = ifelse(is.na(Beta.lactam), "", as.character(Beta.lactam)),
    has_ctx_m_15 = grepl("CTX-M-15", amr_beta, ignore.case = TRUE) | grepl("CTX-M-15", beta_nat, ignore.case = TRUE),
    has_oxa_1    = grepl("OXA-1($|[^0-9])", amr_beta, ignore.case = TRUE) | grepl("OXA-1($|[^0-9])", beta_nat, ignore.case = TRUE),
    has_tem_1    = grepl("TEM-1|TEM\\b", amr_beta, ignore.case = TRUE) | grepl("TEM", beta_nat, ignore.case = TRUE),
    tripartite   = has_ctx_m_15 & has_oxa_1 & has_tem_1,

    # Sequence type (ST) cleaning
    st_clean = str_trim(as.character(ST)),
    st_clean = ifelse(is.na(st_clean) | st_clean %in% c("", "NA", "-", "ND"), "Unknown", st_clean),
    st_clean = ifelse(!grepl("^ST", st_clean) & st_clean != "Unknown", paste0("ST", st_clean), st_clean)
  ) %>%
  filter(!is.na(Taxon_clean))

# Exact per-species denominators
species_counts <- df_clean %>% count(Taxon_clean, name = "N_Species")

cat("--> Strict per-species denominators:\n")
print(species_counts)

# Axis label expressions with actual isolate counts
n_eco <- species_counts$N_Species[species_counts$Taxon_clean == "Escherichia coli"]
n_kpn <- species_counts$N_Species[species_counts$Taxon_clean == "Klebsiella pneumoniae"]
n_ecc <- species_counts$N_Species[species_counts$Taxon_clean == "Enterobacter cloacae complex"]

species_labels_expr <- c(
  "Escherichia coli"             = parse(text = sprintf("atop(italic('E. coli'), '(n=%d)')", n_eco)),
  "Klebsiella pneumoniae"        = parse(text = sprintf("atop(italic('K. pneumoniae'), '(n=%d)')", n_kpn)),
  "Enterobacter cloacae complex" = parse(text = sprintf("atop(italic('E. cloacae')~'complex', '(n=%d)')", n_ecc))
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
      plot.margin        = margin(6, 12, 6, 6)
    )
}

# ------------------------------------------------------------------------------
# 1. PANEL 1A: ACQUISITION ORIGIN (COMMUNITY-IMPORTED VS NOSOCOMIAL)
# ------------------------------------------------------------------------------
df_1A <- df_clean %>%
  group_by(Taxon_clean) %>%
  summarise(
    N_Imported   = sum(is_imported),
    N_Nosocomial = sum(is_nosocomial),
    Total        = n(),
    .groups      = "drop"
  ) %>%
  pivot_longer(cols = c(N_Imported, N_Nosocomial), names_to = "Origin_Type", values_to = "Count") %>%
  mutate(
    Origin = factor(
      ifelse(Origin_Type == "N_Imported",
             "Community / Imported (<= 2 days)",
             "Nosocomial Acquisition (> 2 days)"),
      levels = c("Community / Imported (<= 2 days)", "Nosocomial Acquisition (> 2 days)")
    ),
    Percentage  = (Count / Total) * 100,
    Taxon_clean = factor(Taxon_clean, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  )

p1A <- ggplot(df_1A, aes(x = Taxon_clean, y = Percentage, fill = Origin)) +
  geom_bar(stat = "identity", position = "fill", width = 0.58) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(%d)", Percentage, Count)),
    position = position_fill(vjust = 0.5),
    size = 3.1, fontface = "bold", color = "white"
  ) +
  scale_x_discrete(labels = species_labels_expr) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      "Community / Imported (<= 2 days)"  = "#2ca02c",
      "Nosocomial Acquisition (> 2 days)" = "#756bb1"
    ),
    name = "Acquisition Origin:"
  ) +
  labs(
    title    = "A  Origin of Acquisition & Delays",
    subtitle = "\u2264 2 days (< 48 h) = Imported vs > 2 days (\u2265 48 h) = Nosocomial",
    x        = NULL,
    y        = "Proportion of Isolates (%)"
  ) +
  theme_pub() +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 8.2),
    legend.text       = element_text(size = 7.8),
    legend.box.margin = margin(t = 2, b = 2),
    plot.margin       = margin(t = 8, r = 18, b = 8, l = 8)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# ------------------------------------------------------------------------------
# 2. PANEL 1B: TRIPARTITE RESISTOME AND TZP RESISTANCE (R)
# ------------------------------------------------------------------------------
df_1B <- df_clean %>%
  group_by(Taxon_clean) %>%
  summarise(
    N_Tripartite = sum(tripartite, na.rm = TRUE),
    N_TZP_R      = sum(tzp_r, na.rm = TRUE),
    Total        = n(),
    .groups      = "drop"
  ) %>%
  pivot_longer(cols = c(N_Tripartite, N_TZP_R), names_to = "Feature_Type", values_to = "Count") %>%
  mutate(
    Feature = factor(
      ifelse(Feature_Type == "N_Tripartite",
             "Tripartite Block (CTX-M-15 + OXA-1 + TEM-1)",
             "Phenotypic TZP Resistance (R)"),
      levels = c("Tripartite Block (CTX-M-15 + OXA-1 + TEM-1)", "Phenotypic TZP Resistance (R)")
    ),
    Percentage  = (Count / Total) * 100,
    Taxon_clean = factor(Taxon_clean, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"))
  )

max_pct_1B <- max(df_1B$Percentage, na.rm = TRUE)
lim_y_1B   <- max(60, ceiling((max_pct_1B + 10) / 10) * 10)

resistome_labels <- c(
  expression("Tripartite Block ("*italic("bla")["CTX-M-15"]~"+"~italic("bla")["OXA-1"]~"+"~italic("bla")["TEM-1"]*")"),
  "Phenotypic TZP Resistance (R)"
)

p1B <- ggplot(df_1B, aes(x = Taxon_clean, y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.68) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(%d)", Percentage, Count)),
    position = position_dodge(width = 0.78),
    vjust = -0.2, size = 3.0, fontface = "bold"
  ) +
  scale_x_discrete(labels = species_labels_expr) +
  scale_y_continuous(limits = c(0, lim_y_1B), labels = function(x) paste0(x, "%")) +
  scale_fill_manual(
    values = c(
      "Tripartite Block (CTX-M-15 + OXA-1 + TEM-1)" = "#d95f02",
      "Phenotypic TZP Resistance (R)"               = "#386cb0"
    ),
    labels = resistome_labels,
    name   = "Resistance Determinant:"
  ) +
  labs(
    title    = "B  Tripartite Resistome & Piperacillin\u2013Tazobactam Resistance",
    subtitle = "Consolidated multiresistance modules enriched in nosocomial opportunistic clades",
    x        = NULL,
    y        = "Prevalence in Species (%)"
  ) +
  theme_pub() +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 8.2),
    legend.text       = element_text(size = 7.8),
    legend.box.margin = margin(t = 2, b = 2),
    plot.margin       = margin(t = 8, r = 8, b = 8, l = 18)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# ------------------------------------------------------------------------------
# 3. PANEL 1C: CLONAL ARCHITECTURE OF SEQUENCE TYPES (ST)
# ------------------------------------------------------------------------------
# Epidemiologically documented high-risk clones tracked in the cohort
target_sts <- list(
  "Escherichia coli"             = c("ST131", "ST69", "ST38", "ST10", "ST1193"),
  "Klebsiella pneumoniae"        = c("ST405", "ST307", "ST584"),
  "Enterobacter cloacae complex" = c("ST66", "ST78", "ST114", "ST145")
)

df_1C_list <- list()

for (sp in names(target_sts)) {
  sub_df <- df_clean %>% filter(Taxon_clean == sp)
  sp_total <- nrow(sub_df)
  sts_tracked <- target_sts[[sp]]

  # Per-ST counts
  st_counts <- sub_df %>%
    count(st_clean, name = "Count") %>%
    mutate(
      Is_Major  = st_clean %in% sts_tracked,
      ST_Group  = ifelse(Is_Major, st_clean, "Other STs")
    )

  # Aggregation
  st_agg <- st_counts %>%
    group_by(ST_Group, Is_Major) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    mutate(
      Percentage = (Count / sp_total) * 100,
      Species    = sp,
      Total_Sp   = sp_total
    )

  # Number of distinct STs aggregated into "Other STs"
  n_other_types <- sum(!st_counts$Is_Major)
  st_agg <- st_agg %>%
    mutate(
      ST_Label = ifelse(ST_Group == "Other STs",
                        sprintf("Other STs (n=%d)", n_other_types),
                        ST_Group)
    )

  df_1C_list[[sp]] <- st_agg
}

df_1C <- bind_rows(df_1C_list) %>%
  mutate(
    Species_Facet = factor(
      sprintf("%s (n=%d)",
              case_when(
                Species == "Escherichia coli"             ~ "E. coli",
                Species == "Klebsiella pneumoniae"        ~ "K. pneumoniae",
                Species == "Enterobacter cloacae complex" ~ "E. cloacae complex"
              ), Total_Sp),
      levels = c(sprintf("E. coli (n=%d)", n_eco),
                 sprintf("K. pneumoniae (n=%d)", n_kpn),
                 sprintf("E. cloacae complex (n=%d)", n_ecc))
    ),
    FillGroup = case_when(
      Species == "Escherichia coli" & Is_Major              ~ "E. coli Clones",
      Species == "Escherichia coli" & !Is_Major             ~ "Other E. coli",
      Species == "Klebsiella pneumoniae" & Is_Major         ~ "K. pneumoniae Clones",
      Species == "Klebsiella pneumoniae" & !Is_Major        ~ "Other K. pneumoniae",
      Species == "Enterobacter cloacae complex" & Is_Major  ~ "E. cloacae Clones",
      TRUE                                                  ~ "Other E. cloacae"
    )
  )

# Bar ordering: "Other STs" at the bottom of each facet
df_1C <- df_1C %>%
  group_by(Species) %>%
  arrange(Is_Major, desc(Percentage), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(ST_Ordered = factor(ST_Label, levels = unique(ST_Label)))

facet_labels <- as_labeller(c(
  setNames(sprintf("atop(italic('E. coli'), '(n=%d)')", n_eco), sprintf("E. coli (n=%d)", n_eco)),
  setNames(sprintf("atop(italic('K. pneumoniae'), '(n=%d)')", n_kpn), sprintf("K. pneumoniae (n=%d)", n_kpn)),
  setNames(sprintf("atop(italic('E. cloacae')~'complex', '(n=%d)')", n_ecc), sprintf("E. cloacae complex (n=%d)", n_ecc))
), default = label_parsed)

max_pct_1C <- max(df_1C$Percentage, na.rm = TRUE)
lim_x_1C   <- max(70, ceiling((max_pct_1C + 10) / 10) * 10)

p1C <- ggplot(df_1C, aes(x = Percentage, y = ST_Ordered, fill = FillGroup)) +
  facet_grid(Species_Facet ~ ., scales = "free_y", space = "free_y", labeller = facet_labels) +
  geom_bar(stat = "identity", width = 0.72) +
  geom_text(
    aes(label = sprintf("%.1f%% (n=%d)", Percentage, Count)),
    hjust = -0.08, size = 2.9, fontface = "plain"
  ) +
  scale_fill_manual(
    values = c(
      "E. coli Clones"        = "#1f78b4",
      "Other E. coli"         = "#a6cee3",
      "K. pneumoniae Clones"  = "#e31a1c",
      "Other K. pneumoniae"   = "#fb9a99",
      "E. cloacae Clones"     = "#ff7f00",
      "Other E. cloacae"      = "#fdbf6f"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, lim_x_1C), labels = function(x) paste0(x, "%")) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "C  Clonal Population Structure and Lineage Diversity",
    subtitle = expression(
      "Polyclonal " * italic("E. coli") * " vs oligoclonal " *
        italic("K. pneumoniae") * " and endemic " * italic("E. cloacae") * " complex ST66 monopoly"
    ),
    x = "Proportion of Species Isolates (%)",
    y = NULL
  ) +
  theme_pub() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y       = element_text(size = 9.5, face = "bold", lineheight = 1.1, margin = margin(l = 10, r = 10)),
    strip.background   = element_rect(fill = "grey92", color = "grey60", linewidth = 0.5),
    plot.margin        = margin(t = 10, r = 35, b = 10, l = 10)
  )

# ------------------------------------------------------------------------------
# 4. ASSEMBLY WITH PATCHWORK AND MULTI-FORMAT EXPORT
# ------------------------------------------------------------------------------
fig1_complete <- (p1A | p1B) / p1C + plot_layout(heights = c(1, 1.25))

fig1_complete

ggsave(file.path(out_dir_fig, "Figure_1_Species_Dichotomy_Resistome.pdf"),
       fig1_complete, width = 13.5, height = 11.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_1_Species_Dichotomy_Resistome.png"),
       fig1_complete, width = 13.5, height = 11.5, dpi = 300)
ggsave(file.path(out_dir_fig, "Figure_1_Species_Dichotomy_Resistome.tiff"),
       fig1_complete, width = 13.5, height = 11.5, dpi = 300, compression = "lzw")

cat(sprintf("\n--> Figure 1 regenerated dynamically in: %s\n", out_dir_fig))
