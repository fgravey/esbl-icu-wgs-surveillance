# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 04_population_structure.R
# Purpose: Full recomputation of the Results section "Clonal Architecture and
#          High-Risk Lineages": sequence type (ST) diversity per species,
#          ST131 / fimH subclade breakdown for E. coli, and dominant lineages
#          of K. pneumoniae and E. cloacae complex.
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/tables/Results_PopStructure_Global.csv
#          output/tables/Results_PopStructure_Ecoli.csv
#          output/tables/Results_PopStructure_Kpn.csv
#          output/tables/Results_PopStructure_Ecc.csv
#          output/tables/Results_ST131_fimH.csv
#          output/tables/Results_ST131_Clades.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# IMPORT
# ------------------------------------------------------------------------------
df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# COHORT PER SPECIES
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SEQUENCE TYPES PER SPECIES\n")
cat("====================================================\n")

st_summary <- df %>%
  group_by(Taxon_group) %>%
  summarise(
    Total_Isolates = n(),
    Distinct_STs   = n_distinct(ST),
    .groups        = "drop"
  )

print(st_summary)

# ------------------------------------------------------------------------------
# E. COLI
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("E. COLI\n")
cat("====================================================\n")

ecoli <- df %>% filter(Taxon_group == "Escherichia coli")

ecoli_st <- ecoli %>%
  count(ST, sort = TRUE) %>%
  mutate(Pct = round(100 * n / sum(n), 1))

print(head(ecoli_st, 20))

# ------------------------------------------------------------------------------
# ST131
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("ST131\n")
cat("====================================================\n")

st131 <- ecoli %>% filter(grepl("131", ST))

cat("ST131 count:", nrow(st131), "\n")
cat("ST131 percentage:", round(100 * nrow(st131) / nrow(ecoli), 1), "\n")

# Automatic detection of fimH / clade columns
fimh_cols <- grep("fimh|fimH|clade", colnames(st131), value = TRUE, ignore.case = TRUE)
print(fimh_cols)

if (length(fimh_cols) > 0) {
  for (col in fimh_cols) {
    cat("\n")
    cat("Variable:", col, "\n")
    print(table(st131[[col]], useNA = "ifany"))
  }
}

# ------------------------------------------------------------------------------
# ST131: FimH SUBTYPES AND HISTORICAL CLADES
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("ST131 FimH SUBTYPES\n")
cat("====================================================\n")

st131 <- ecoli %>% filter(grepl("131", ST))

cat("Total ST131:", nrow(st131), "\n")

fimh_table <- st131 %>%
  count(FimH, sort = TRUE) %>%
  mutate(Pct = round(100 * n / sum(n), 1))

print(fimh_table)

# Historical clades (fimH-based typing)
st131_clades <- st131 %>%
  mutate(
    ST131_Clade = case_when(
      grepl("30", FimH) ~ "Clade C / H30",
      grepl("41", FimH) ~ "Clade A / H41",
      grepl("22", FimH) ~ "Clade B / H22",
      TRUE ~ "Other fimH variants"
    )
  ) %>%
  count(ST131_Clade, sort = TRUE) %>%
  mutate(Pct = round(100 * n / sum(n), 1))

print(st131_clades)

write.csv(fimh_table, file.path(out_dir_tab, "Results_ST131_fimH.csv"), row.names = FALSE)
write.csv(st131_clades, file.path(out_dir_tab, "Results_ST131_Clades.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# KLEBSIELLA PNEUMONIAE
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("KLEBSIELLA PNEUMONIAE\n")
cat("====================================================\n")

kpn <- df %>% filter(Taxon_group == "Klebsiella pneumoniae")

kpn_st <- kpn %>%
  count(ST, sort = TRUE) %>%
  mutate(Pct = round(100 * n / sum(n), 1))

print(head(kpn_st, 20))

# Top 3 lineages
top3_kpn <- sum(head(kpn_st$n, 3))

cat("\nTop 3 Kpn lineages:", top3_kpn, "/", nrow(kpn), "=",
    round(100 * top3_kpn / nrow(kpn), 1), "%\n")

# ------------------------------------------------------------------------------
# ENTEROBACTER CLOACAE COMPLEX
# ------------------------------------------------------------------------------
cat("\n====================================================\n")
cat("ENTEROBACTER CLOACAE COMPLEX\n")
cat("====================================================\n")

ecc <- df %>% filter(Taxon_group == "Enterobacter cloacae complex")

ecc_st <- ecc %>%
  count(ST, sort = TRUE) %>%
  mutate(Pct = round(100 * n / sum(n), 1))

print(head(ecc_st, 20))

# ------------------------------------------------------------------------------
# ST66 (endemic E. cloacae complex lineage)
# ------------------------------------------------------------------------------
st66 <- ecc %>% filter(grepl("66", ST))

cat("\nST66:", nrow(st66), "/", nrow(ecc), "=",
    round(100 * nrow(st66) / nrow(ecc), 1), "%\n")

# ------------------------------------------------------------------------------
# EXPORT
# ------------------------------------------------------------------------------
write.csv(ecoli_st, file.path(out_dir_tab, "Results_PopStructure_Ecoli.csv"), row.names = FALSE)
write.csv(kpn_st,    file.path(out_dir_tab, "Results_PopStructure_Kpn.csv"),    row.names = FALSE)
write.csv(ecc_st,    file.path(out_dir_tab, "Results_PopStructure_Ecc.csv"),    row.names = FALSE)
write.csv(st_summary, file.path(out_dir_tab, "Results_PopStructure_Global.csv"), row.names = FALSE)

cat("\nPopulation structure exports completed.\n")
