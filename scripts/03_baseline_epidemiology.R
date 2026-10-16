# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 03_baseline_epidemiology.R
# Purpose: Full recomputation of the Results section "Overall Cohort Dynamics
#          and Baseline Epidemiology": cohort size, species distribution,
#          ward distribution, sample types, delays from admission to first
#          ESBL-E isolation, and multi-species co-colonisation.
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/tables/Results_Baseline_Species.csv
#          output/tables/Results_Baseline_Delays.csv
#          output/tables/Results_Baseline_Multispecies.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. IMPORT
# ------------------------------------------------------------------------------
df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. PREPARATION
# ------------------------------------------------------------------------------
df <- df %>%
  mutate(
    d_entree = as.Date(parse_date_time(Date.entree, orders = c("dmy", "ymd", "d/m/Y", "Y-m-d"))),
    d_prel   = as.Date(parse_date_time(Date.de.prel., orders = c("dmy", "ymd", "d/m/Y", "Y-m-d"))),

    los_before_isolation = pmax(0, as.numeric(difftime(d_prel, d_entree, units = "days"))),

    # Convention: isolates recovered <= 2 days after admission are considered
    # community-imported; isolates recovered > 2 days are hospital-acquired
    acquisition_status = ifelse(los_before_isolation <= 2,
                                 "Community-acquired", "Hospital-acquired")
  )

# ------------------------------------------------------------------------------
# 3. OVERALL COHORT
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("OVERALL COHORT\n")
cat("====================================================\n")

cat("Episodes:", nrow(df), "\n")
cat("Unique patients:", n_distinct(df$Numero.patient), "\n")

# ------------------------------------------------------------------------------
# 4. ISOLATE COUNTS BY SPECIES
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("MAIN SPECIES\n")
cat("====================================================\n")

species_table <- df %>%
  count(Taxon_group, name = "N") %>%
  mutate(Pct = round(100 * N / sum(N), 1))

print(species_table)

# ------------------------------------------------------------------------------
# 5. SURGICAL ICU / MEDICAL ICU DISTRIBUTION
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("DISTRIBUTION BY WARD\n")
cat("====================================================\n")

service_table <- df %>%
  count(cores, name = "N") %>%
  mutate(Pct = round(100 * N / sum(N), 1))

print(service_table)

# ------------------------------------------------------------------------------
# 6. SCREENING VS CLINICAL SAMPLES
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SAMPLE TYPE\n")
cat("====================================================\n")

df <- df %>%
  mutate(
    sample_type = case_when(
      grepl("SELLES", Nature, ignore.case = TRUE) ~ "Rectal screening",
      TRUE ~ "Clinical sample"
    )
  )

sample_table <- df %>%
  count(sample_type, name = "N") %>%
  mutate(Pct = round(100 * N / sum(N), 1))

print(sample_table)

# ------------------------------------------------------------------------------
# 7. DETAILS OF CLINICAL SAMPLES
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("CLINICAL SAMPLE DETAILS\n")
cat("====================================================\n")

clinical_sites <- df %>%
  filter(sample_type == "Clinical sample") %>%
  mutate(
    Site = case_when(
      grepl("URIN", Nature, ignore.case = TRUE) ~ "Urine",
      grepl("HEMOC", Nature, ignore.case = TRUE) ~ "Blood culture",
      grepl("BRONCH|LBA|TRACHE|PLEUR", Nature, ignore.case = TRUE) ~ "Respiratory",
      grepl("PUS|LIQ PERITO|PDP|BILE|ARTIC", Nature, ignore.case = TRUE) ~ "Wound / Deep collection",
      grepl("CATHETER", Nature, ignore.case = TRUE) ~ "Catheter tip",
      grepl("LCR", Nature, ignore.case = TRUE) ~ "CSF",
      TRUE ~ "Other"
    )
  ) %>%
  count(Site, sort = TRUE)

print(clinical_sites)

# ------------------------------------------------------------------------------
# 8. OVERALL DELAY BEFORE FIRST ESBL-E ISOLATION
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("OVERALL DELAY BEFORE ISOLATION\n")
cat("====================================================\n")

summary_los <- df %>%
  summarise(
    Median = median(los_before_isolation, na.rm = TRUE),
    Q1     = quantile(los_before_isolation, 0.25, na.rm = TRUE),
    Q3     = quantile(los_before_isolation, 0.75, na.rm = TRUE),
    Min    = min(los_before_isolation, na.rm = TRUE),
    Max    = max(los_before_isolation, na.rm = TRUE)
  )

print(summary_los)

# ------------------------------------------------------------------------------
# 9. DELAYS BY SPECIES
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("DELAYS BY SPECIES\n")
cat("====================================================\n")

species_delay <- df %>%
  group_by(Taxon_group) %>%
  summarise(
    N = n(),
    Median = median(los_before_isolation, na.rm = TRUE),
    Q1     = quantile(los_before_isolation, 0.25, na.rm = TRUE),
    Q3     = quantile(los_before_isolation, 0.75, na.rm = TRUE),
    Community_N   = sum(acquisition_status == "Community-acquired"),
    Community_Pct = round(100 * Community_N / N, 1),
    Hospital_N   = sum(acquisition_status == "Hospital-acquired"),
    Hospital_Pct = round(100 * Hospital_N / N, 1),
    .groups = "drop"
  )

print(species_delay)

# ------------------------------------------------------------------------------
# 10. MULTI-SPECIES CO-COLONISATION
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("MULTI-SPECIES CO-COLONISATION\n")
cat("====================================================\n")

patient_taxa <- df %>% distinct(Numero.patient, Taxon_group)

multi_species <- patient_taxa %>%
  count(Numero.patient) %>%
  filter(n > 1)

cat("Multi-species patients:", nrow(multi_species), "\n")

pairwise_species <- patient_taxa %>%
  filter(Numero.patient %in% multi_species$Numero.patient) %>%
  inner_join(patient_taxa, by = "Numero.patient", relationship = "many-to-many") %>%
  filter(Taxon_group.x < Taxon_group.y) %>%
  count(Taxon_group.x, Taxon_group.y, sort = TRUE)

print(pairwise_species)

# ------------------------------------------------------------------------------
# 11. EXPORT
# ------------------------------------------------------------------------------
write.csv(species_table, file.path(out_dir_tab, "Results_Baseline_Species.csv"), row.names = FALSE)
write.csv(species_delay, file.path(out_dir_tab, "Results_Baseline_Delays.csv"), row.names = FALSE)
write.csv(pairwise_species, file.path(out_dir_tab, "Results_Baseline_Multispecies.csv"), row.names = FALSE)

cat("\nBaseline epidemiology exports completed.\n")
