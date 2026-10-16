# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 02_genomic_transmission_linkage.R
# Purpose: Define intra-ward genomic cross-transmission events based on
#          core-genome SNP (cgSNP) clusters combined with a rolling <= 30-day
#          temporal window (nearest temporal neighbour from a different
#          patient within the same ward).
# Input :  data/BLSE_cohort_curated_clean.tsv  (output of script 01)
# Output:  data/BLSE_cohort_final_with_transmissions.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")

cat("--> Step 2: Defining secondary genomic cross-transmissions...\n")
df <- read.csv(file.path(data_dir, "BLSE_cohort_curated_clean.tsv"), sep = "\t", stringsAsFactors = FALSE)

# 1. Identify valid core-genome clusters
df_clustered <- df %>%
  mutate(
    snp_cluster_clean = str_trim(as.character(snp_clusters)),
    is_valid_cluster = case_when(
      is.na(snp_cluster_clean) ~ FALSE,
      snp_cluster_clean %in% c("", "non", "Non", "singleton", "Singleton",
                               "fastq perdu", "genome manquant", "NA", "0", "-") ~ FALSE,
      TRUE ~ TRUE
    ),
    ward_clean = case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med", cores,  ignore.case = TRUE) ~ "Medical ICU",
      TRUE ~ str_trim(as.character(cores))
    )
  )

# 2. Compute the minimal inter-isolate time delay within the same ward and
#    the same cgSNP cluster
df_linked <- df_clustered %>%
  group_by(ward_clean, snp_cluster_clean) %>%
  arrange(date_sampling) %>%
  mutate(
    ward_cluster_patient_count = n_distinct(Correspondant[is_valid_cluster]),
    # Secondary transmission flag:
    # requires a valid cluster, >= 2 distinct patients in the ward, and a
    # delay <= 30 days to another patient of the same cluster
    is_cross_transmission_30d = case_when(
      !is_valid_cluster ~ FALSE,
      ward_cluster_patient_count < 2 ~ FALSE,
      is.na(date_sampling) ~ FALSE,
      sapply(seq_along(date_sampling), function(i) {
        current_patient <- Correspondant[i]
        current_date    <- date_sampling[i]
        other_dates     <- date_sampling[Correspondant != current_patient & !is.na(date_sampling)]
        if (length(other_dates) == 0) return(FALSE)
        min_delta <- min(abs(as.numeric(difftime(current_date, other_dates, units = "days"))))
        return(min_delta <= 30)
      }) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

cat("\nSummary of confirmed genomic cross-transmissions (<= 30 days):\n")
print(table(df_linked$Taxon_group, df_linked$is_cross_transmission_30d))

output_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
write.table(df_linked, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nConsolidated transmission database exported to: %s\n\n", output_file))
