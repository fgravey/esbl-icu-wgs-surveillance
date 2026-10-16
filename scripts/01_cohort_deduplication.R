# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 01_cohort_deduplication.R
# Purpose: Strict epidemiological de-duplication of the ICU cohort.
#          Only the chronologically earliest isolate per patient / stay /
#          species / sequence type (ST) combination is retained.
# Input :  data/BLSE_cohort_curated.tsv  (raw curated isolate-level database)
# Output:  data/BLSE_cohort_curated_clean.tsv
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")

cat("--> Step 1: Loading raw curated database...\n")
input_file <- file.path(data_dir, "BLSE_cohort_curated.tsv")
df_raw <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)

cat(sprintf("Initial isolate count: %d\n", nrow(df_raw)))
cat(sprintf("Initial unique patients: %d\n", n_distinct(df_raw$Numero.patient)))

# 1. Parse dates safely (tolerate heterogeneous date formats)
df_prep <- df_raw %>%
  mutate(
    date_sampling  = as.Date(parse_date_time(Date.de.prel., orders = c("dmy", "ymd", "d/m/Y", "Y-m-d"))),
    date_admission = as.Date(parse_date_time(Date.entree,  orders = c("dmy", "ymd", "d/m/Y", "Y-m-d")))
  ) %>%
  filter(!is.na(date_sampling))

# 2. Strict episode-level de-duplication:
#    Keep the earliest isolate per Patient / Stay / Species / ST combination
df_clean <- df_prep %>%
  arrange(Numero.patient, Numero.sejour, Taxon_group, ST, date_sampling) %>%
  group_by(Numero.patient, Numero.sejour, Taxon_group, ST) %>%
  slice(1) %>%
  ungroup()

cat(sprintf("De-duplicated isolate count: %d (removed %d redundant isolates)\n",
            nrow(df_clean), nrow(df_raw) - nrow(df_clean)))
cat(sprintf("Final unique patients: %d\n", n_distinct(df_clean$Numero.patient)))

# 3. Export the cleaned dataset
output_file <- file.path(data_dir, "BLSE_cohort_curated_clean.tsv")
write.table(df_clean, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Cleaned database exported successfully to: %s\n\n", output_file))
