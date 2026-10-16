# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 13_figure_s3_plasmid_architecture.R
# Purpose: Dynamic analysis of PlasmidFinder replicon profiles and generation
#          of Supplementary Figure S3 with standardized publication colors,
#          plus the underlying replicon frequency table (Supplementary Table S3).
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
#          data/article-blse-plasmidfinder.csv
# Output:  output/figures/Supplementary_Figure_S3_Plasmid_Architecture.pdf (.png)
#          output/tables/Supplementary_Table_S3_Plasmids.csv
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(scales)
})

base_dir    <- "."
data_dir    <- file.path(base_dir, "data")
out_dir_fig <- file.path(base_dir, "output", "figures")
out_dir_tab <- file.path(base_dir, "output", "tables")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. LOADING OF THE FINAL CLEANED COHORT AND PLASMIDFINDER RESULTS
# ------------------------------------------------------------------------------
cohort_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
df_cohort <- read.csv(cohort_file, sep = "\t", stringsAsFactors = FALSE)

# PlasmidFinder export: the first 6 lines are metadata header lines and are skipped
pf_file <- file.path(data_dir, "article-blse-plasmidfinder.csv")
df_pf_raw <- read.csv(pf_file, sep = ";", skip = 6, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. CLEANING AND CLASSIFICATION OF PLASMID REPLICONS
# ------------------------------------------------------------------------------
# Extract base strain names for matching with the cohort
df_pf <- df_pf_raw %>%
  filter(!is.na(plasmids), plasmids != "") %>%
  mutate(
    # Clean isolate identifier (e.g. Ec20190504.scfd.fasta -> Ec20190504)
    isolate_id = str_remove(Souche, "\\..*$"),
    plasmid_raw = trimws(plasmids)
  )

# Harmonized classification of replicon families and specific types
df_pf_tagged <- df_pf %>%
  mutate(
    # Binary indicators of specific replicons
    is_IncFIB_K   = grepl("IncFIB\\(K\\)", plasmid_raw, ignore.case = TRUE),
    is_IncFII_K   = grepl("IncFII\\(K\\)", plasmid_raw, ignore.case = TRUE),
    is_IncFIB_AP  = grepl("IncFIB\\(AP001918\\)|IncFIB$", plasmid_raw, ignore.case = TRUE) & !is_IncFIB_K,
    is_IncFIA     = grepl("IncFIA", plasmid_raw, ignore.case = TRUE),
    is_IncFII_can = grepl("^IncFII$", plasmid_raw, ignore.case = TRUE) | grepl("IncFII\\(p", plasmid_raw, ignore.case = TRUE) & !is_IncFII_K,

    is_IncHI1     = grepl("IncHI1", plasmid_raw, ignore.case = TRUE),
    is_IncHI2     = grepl("IncHI2", plasmid_raw, ignore.case = TRUE),
    is_Col_pHAD28 = grepl("Col\\(pHAD28\\)", plasmid_raw, ignore.case = TRUE),
    is_Col_any    = grepl("^Col", plasmid_raw, ignore.case = TRUE),

    is_IncI1      = grepl("IncI1", plasmid_raw, ignore.case = TRUE),
    is_IncBOKZ    = grepl("IncB|IncO|IncK|IncZ", plasmid_raw, ignore.case = TRUE),
    is_IncX       = grepl("IncX", plasmid_raw, ignore.case = TRUE),
    is_IncR       = grepl("IncR", plasmid_raw, ignore.case = TRUE),

    # Large-family indicators
    is_Any_IncF   = grepl("^IncF", plasmid_raw, ignore.case = TRUE),
    is_Any_IncHI  = grepl("^IncHI", plasmid_raw, ignore.case = TRUE)
  )

# Aggregate per isolate to obtain one row per PlasmidFinder-typed isolate
df_pf_summary <- df_pf_tagged %>%
  group_by(isolate_id) %>%
  summarise(
    `Any IncF`              = any(is_Any_IncF),
    `IncFIB(K)`             = any(is_IncFIB_K),
    `IncFII(K)`             = any(is_IncFII_K),
    `IncFIB(AP001918)`      = any(is_IncFIB_AP),
    `IncFIA`                = any(is_IncFIA),
    `IncFII (canonical)`    = any(is_IncFII_can),
    `Any IncHI`             = any(is_Any_IncHI),
    `IncHI1`                = any(is_IncHI1),
    `IncHI2`                = any(is_IncHI2),
    `Any Col-type`          = any(is_Col_any),
    `Col(pHAD28)`           = any(is_Col_pHAD28),
    `IncI1-I(alpha)`        = any(is_IncI1),
    `IncB/O/K/Z`            = any(is_IncBOKZ),
    `IncX lineage`          = any(is_IncX),
    `IncR`                  = any(is_IncR),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 3. DYNAMIC MAPPING WITH THE FINAL CLEANED COHORT
# ------------------------------------------------------------------------------
df_cohort_clean <- df_cohort %>%
  mutate(
    isolate_match = case_when(
      "Nom.etude" %in% colnames(df_cohort) ~ str_remove(trimws(Nom.etude), "\\..*$"),
      TRUE ~ as.character(Numero.patient)
    ),
    Species = case_when(
      grepl("coli", Taxon_group, ignore.case = TRUE)        ~ "Escherichia coli",
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE)  ~ "Klebsiella pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)     ~ "Enterobacter cloacae complex",
      TRUE                                                  ~ NA_character_
    )
  ) %>%
  filter(!is.na(Species))

# Exact cohort denominators
denom_species <- df_cohort_clean %>%
  count(Species, name = "N_Total")

cat("--> Strict per-species denominators for PlasmidFinder:\n")
print(denom_species)

# Merge cohort + PlasmidFinder data
df_merged <- df_cohort_clean %>%
  left_join(df_pf_summary, by = c("isolate_match" = "isolate_id")) %>%
  mutate(across(where(is.logical), ~ coalesce(., FALSE)))

# ------------------------------------------------------------------------------
# 4. DYNAMIC COMPUTATION OF PER-REPLICON / PER-SPECIES FREQUENCIES
# ------------------------------------------------------------------------------
replicon_cols <- c(
  "Any IncF", "IncFIB(K)", "IncFII(K)", "IncFIB(AP001918)", "IncFIA", "IncFII (canonical)",
  "Any IncHI", "IncHI1", "IncHI2",
  "Any Col-type", "Col(pHAD28)",
  "IncI1-I(alpha)", "IncB/O/K/Z", "IncX lineage", "IncR"
)

results_list <- list()

for (rep_name in replicon_cols) {
  df_sub <- df_merged %>%
    group_by(Species) %>%
    summarise(
      Count = sum(.data[[rep_name]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(denom_species, by = "Species") %>%
    mutate(
      Replicon = rep_name,
      Percentage = (Count / N_Total) * 100
    )
  results_list[[rep_name]] <- df_sub
}

df_plot <- bind_rows(results_list) %>%
  mutate(
    # Strict replicon order (identical to the original figure)
    Replicon = factor(Replicon, levels = rev(replicon_cols)),
    Species  = factor(Species, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex")),
    Label_Text = sprintf("%.1f%%", Percentage)
  )

# Export of the underlying data table (Supplementary Table S3)
write_csv(df_plot %>% select(Species, Replicon, Count, Denominator = N_Total, Percentage),
          file.path(out_dir_tab, "Supplementary_Table_S3_Plasmids.csv"))

# ------------------------------------------------------------------------------
# 5. OFFICIAL THE LANCET MICROBE RENDERING (SUPPLEMENTARY FIGURE S3)
# ------------------------------------------------------------------------------
# Harmonized official palette
species_colors <- c(
  "Escherichia coli"             = "#1f78b4",
  "Klebsiella pneumoniae"        = "#e31a1c",
  "Enterobacter cloacae complex" = "#ff7f00"
)

p_s3 <- ggplot(df_plot, aes(x = Percentage, y = Replicon, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = Label_Text),
    position = position_dodge(width = 0.8),
    hjust = -0.15, size = 2.9, fontface = "bold", color = "black"
  ) +
  scale_fill_manual(values = species_colors, name = "Bacterial species:") +
  scale_x_continuous(
    limits = c(0, 106),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    title    = "Supplementary Figure S3. Comparative Plasmid Replicon Architecture Across Species",
    subtitle = sprintf("In silico PlasmidFinder profiling highlighting the IncF vs. IncHI / Col dichotomy (N = %d isolates)", nrow(df_cohort_clean)),
    x        = "Isolates harboring replicon (%)",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = rel(1.05), hjust = 0),
    plot.subtitle      = element_text(size = rel(0.88), color = "grey30", margin = margin(b = 10)),
    axis.text.y        = element_text(face = "bold", size = rel(0.9), color = "black"),
    axis.text.x        = element_text(size = rel(0.85), color = "black"),
    legend.position    = "top",
    legend.justification = "left",
    legend.title       = element_text(face = "bold", size = rel(0.9)),
    legend.text        = element_text(size = rel(0.85)),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35, linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(10, 25, 10, 10)
  )

# Vector and high-resolution exports
ggsave(file.path(out_dir_fig, "Supplementary_Figure_S3_Plasmid_Architecture.pdf"),
       p_s3, width = 12.0, height = 9.0, dpi = 300)
ggsave(file.path(out_dir_fig, "Supplementary_Figure_S3_Plasmid_Architecture.png"),
       p_s3, width = 12.0, height = 9.0, dpi = 300)

cat(sprintf("\n--> Figure S3 successfully generated in: %s\n", out_dir_fig))
cat(sprintf("--> Table S3 successfully generated in: %s\n", out_dir_tab))
