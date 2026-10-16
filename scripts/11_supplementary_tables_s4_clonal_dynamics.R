# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 11_supplementary_tables_s4_clonal_dynamics.R
# Purpose: Dynamic computation of Supplementary Tables S4a, S4b, S4c (and S4d):
#          - Table S4a: Global summary of longitudinal clonal profiles
#          - Table S4b: Persistent / endemic lineages (>= 4 years)
#          - Table S4c: Emerging lineages (appeared after 2019, present in 2023)
#          - Table S4d: Disappeared / extinct lineages (absent in 2022-2023)
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/tables/Supplementary_Table_S4a_Summary.csv
#          output/tables/Supplementary_Table_S4b_Persistent.csv
#          output/tables/Supplementary_Table_S4c_Emerging.csv
#          output/tables/Supplementary_Table_S4d_Disappeared.csv
#          output/tables/Supplementary_Table_S4_Clonal_Dynamics.xlsx
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(openxlsx)
})

base_dir    <- "."
data_dir    <- file.path(base_dir, "data")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
df <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. ISOLATE AND YEAR PREPARATION (2019-2023)
# ------------------------------------------------------------------------------
df_clean <- df %>%
  mutate(
    # Robust year extraction from annee_semaine or date_sampling
    annee = case_when(
      !is.na(annee_semaine) ~ as.numeric(substr(as.character(annee_semaine), 1, 4)),
      !is.na(date_sampling) ~ as.numeric(substr(as.character(date_sampling), 1, 4)),
      TRUE                  ~ NA_real_
    ),
    Taxon_group = case_when(
      grepl("coli", Taxon_group, ignore.case = TRUE)        ~ "E. coli",
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE)  ~ "K. pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)     ~ "E. cloacae complex",
      TRUE                                                  ~ NA_character_
    ),
    ST_clean = str_trim(as.character(ST)),
    ST_clean = ifelse(is.na(ST_clean) | ST_clean %in% c("", "NA", "-", "ND", "Unknown"), "Unknown", ST_clean),
    ST_clean = ifelse(!grepl("^ST", ST_clean) & ST_clean != "Unknown", paste0("ST", ST_clean), ST_clean)
  ) %>%
  filter(!is.na(Taxon_group) & ST_clean != "Unknown" & annee %in% 2019:2023)

cat(sprintf("--> MLST-typed isolates available for analysis (2019-2023): %d\n", nrow(df_clean)))

# ------------------------------------------------------------------------------
# 2. LONGITUDINAL CHARACTERIZATION AND CLONAL PROFILE MATRIX
# ------------------------------------------------------------------------------
st_dynamics <- df_clean %>%
  group_by(Taxon_group, ST = ST_clean) %>%
  summarise(
    total_strains            = n(),
    across_nb_years          = n_distinct(annee),
    First_year_of_appearance = min(annee, na.rm = TRUE),
    Last_year_of_appearance  = max(annee, na.rm = TRUE),
    years_of_presence        = paste(sort(unique(annee)), collapse = ", "),
    .groups                  = "drop"
  ) %>%
  mutate(
    temporal_profile = case_when(
      across_nb_years >= 4 ~ "Persistent_endemic (>= 4 years)",
      First_year_of_appearance > 2019 & Last_year_of_appearance == 2023 ~ "Emerging (appeared after 2019, present in 2023)",
      Last_year_of_appearance < 2022 ~ "Disappeared (absent in 2022-2023)",
      TRUE ~ "Sporadic"
    )
  )

# ------------------------------------------------------------------------------
# 3. GENERATION OF THE 4 TABLES (S4a, S4b, S4c, S4d)
# ------------------------------------------------------------------------------

# --- Table S4a: Global summary of temporal profiles ---
table_s4a <- st_dynamics %>%
  group_by(Taxon_group, temporal_profile) %>%
  summarise(
    nb_ST         = n(),
    total_strains = sum(total_strains),
    .groups       = "drop"
  ) %>%
  group_by(Taxon_group) %>%
  mutate(
    pct_strains = round((total_strains / sum(total_strains)) * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(Taxon_group, desc(total_strains))

# --- Table S4b: Persistent / endemic STs (>= 4 years) ---
table_s4b <- st_dynamics %>%
  filter(temporal_profile == "Persistent_endemic (>= 4 years)") %>%
  select(
    Taxon_group, ST, total_strains, across_nb_years,
    First_year_of_appearance, Last_year_of_appearance,
    years_of_presence, temporal_profile
  ) %>%
  arrange(Taxon_group, desc(total_strains))

# --- Table S4c: Emerging STs (appeared after 2019, present in 2023) ---
table_s4c <- st_dynamics %>%
  filter(temporal_profile == "Emerging (appeared after 2019, present in 2023)") %>%
  select(
    Taxon_group, ST, total_strains, across_nb_years,
    First_year_of_appearance, Last_year_of_appearance,
    years_of_presence, temporal_profile
  ) %>%
  arrange(Taxon_group, desc(total_strains), ST)

# --- Table S4d: Disappeared / extinct STs (absent in 2022-2023) ---
table_s4d <- st_dynamics %>%
  filter(temporal_profile == "Disappeared (absent in 2022-2023)") %>%
  select(
    Taxon_group, ST, total_strains, across_nb_years,
    First_year_of_appearance, Last_year_of_appearance,
    years_of_presence, temporal_profile
  ) %>%
  arrange(Taxon_group, desc(total_strains), ST)

# ------------------------------------------------------------------------------
# 4. CSV EXPORTS AND OFFICIAL MULTI-SHEET EXCEL WORKBOOK
# ------------------------------------------------------------------------------
# Individual CSV exports
write_csv(table_s4a, file.path(out_dir_tab, "Supplementary_Table_S4a_Summary.csv"))
write_csv(table_s4b, file.path(out_dir_tab, "Supplementary_Table_S4b_Persistent.csv"))
write_csv(table_s4c, file.path(out_dir_tab, "Supplementary_Table_S4c_Emerging.csv"))
write_csv(table_s4d, file.path(out_dir_tab, "Supplementary_Table_S4d_Disappeared.csv"))

# Formatted Excel workbook
wb <- createWorkbook()

header_style <- createStyle(
  fontSize = 11, fontColour = "#FFFFFF", fgFill = "#00468B",
  halign = "center", valign = "center", textDecoration = "bold"
)

add_formatted_sheet <- function(wb, sheet_name, df_data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df_data, headerStyle = header_style)
  setColWidths(wb, sheet = sheet_name, cols = 1:ncol(df_data), widths = "auto")
}

add_formatted_sheet(wb, "Table_S4a_Summary", table_s4a)
add_formatted_sheet(wb, "Table_S4b_Persistent", table_s4b)
add_formatted_sheet(wb, "Table_S4c_Emerging", table_s4c)
add_formatted_sheet(wb, "Table_S4d_Disappeared", table_s4d)

excel_output <- file.path(out_dir_tab, "Supplementary_Table_S4_Clonal_Dynamics.xlsx")
saveWorkbook(wb, excel_output, overwrite = TRUE)

cat(sprintf("\n--> Tables S4a, S4b, S4c and S4d successfully generated in:\n%s\n", out_dir_tab))
