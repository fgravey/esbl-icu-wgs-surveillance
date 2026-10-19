# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 05_resistome_table2.R
# Purpose: Dynamic computation of EUCAST phenotypes and AMRFinderPlus resistome
#          markers, generation of Table 2 (resistome & susceptibility summary)
#          and Supplementary Table S5 (full S/SFP/R susceptibility profiles
#          across 14 antimicrobials).
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/tables/Table_2_Resistome_and_Susceptibility.csv
#          output/tables/Supplementary_Table_S5_Full_Phenotypes.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. DYNAMIC PHENOTYPIC AND GENOTYPIC VARIABLE EXTRACTION
# ------------------------------------------------------------------------------
taxons <- df$Taxon_group
species_levels <- c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex")

# Helper: phenotypic resistance category (R) for a given antibiotic column
is_r <- function(col_name) {
  toupper(trimws(as.character(df[[col_name]]))) == "R"
}

# Resistome determinants (AMRFinderPlus & Beta.lactam fields)
amr_beta <- ifelse(is.na(df$AMR_BETA_LACTAM), "", as.character(df$AMR_BETA_LACTAM))
amr_quin <- ifelse(is.na(df$AMR_QUINOLONE), "", as.character(df$AMR_QUINOLONE))
amr_amg  <- ifelse(is.na(df$AMR_AMINOGLYCOSIDE), "", as.character(df$AMR_AMINOGLYCOSIDE))
beta_nat <- ifelse(is.na(df$Beta.lactam), "", as.character(df$Beta.lactam))

has_ctx_m_15 <- grepl("CTX-M-15", amr_beta, ignore.case = TRUE) | grepl("CTX-M-15", beta_nat, ignore.case = TRUE)
has_ctx_m_27 <- grepl("CTX-M-27", amr_beta, ignore.case = TRUE) | grepl("CTX-M-27", beta_nat, ignore.case = TRUE)
has_ctx_m_14 <- grepl("CTX-M-14", amr_beta, ignore.case = TRUE) | grepl("CTX-M-14", beta_nat, ignore.case = TRUE)
has_shv_12   <- grepl("SHV-12", amr_beta, ignore.case = TRUE)   | grepl("SHV-12", beta_nat, ignore.case = TRUE)
has_oxa_1     <- grepl("OXA-1($|[^0-9])", amr_beta, ignore.case = TRUE) | grepl("OXA-1($|[^0-9])", beta_nat, ignore.case = TRUE)
has_tem_1     <- grepl("TEM-1|TEM\\b", amr_beta, ignore.case = TRUE) | grepl("TEM", beta_nat, ignore.case = TRUE)
has_tripartite <- has_ctx_m_15 & has_oxa_1 & has_tem_1

# Carbapenemases
has_kpc   <- grepl("KPC", amr_beta, ignore.case = TRUE)
has_oxa48 <- grepl("OXA-48|OXA-181|OXA-232", amr_beta, ignore.case = TRUE)
has_ndm   <- grepl("NDM", amr_beta, ignore.case = TRUE)
has_cpe   <- has_kpc | has_oxa48 | has_ndm

# Quinolones
has_gyra <- grepl("gyrA", amr_quin, ignore.case = TRUE)
has_parc <- grepl("parC", amr_quin, ignore.case = TRUE)
has_qnr  <- grepl("qnr", amr_quin, ignore.case = TRUE)

# Aminoglycosides
has_aac6 <- grepl("aac\\(6'\\)|aac6", amr_amg, ignore.case = TRUE)
has_aac3 <- grepl("aac\\(3\\)|aac3", amr_amg, ignore.case = TRUE)
has_ant2 <- grepl("ant\\(2''\\)|ant2", amr_amg, ignore.case = TRUE)
has_arm  <- grepl("armA|rmt[A-H]", amr_amg, ignore.case = TRUE)

# ------------------------------------------------------------------------------
# 2. DYNAMIC TABLE 1 GENERATION
# ------------------------------------------------------------------------------
format_table2_row <- function(cond_vec) {
  eco_n <- sum(cond_vec[taxons == "Escherichia coli"], na.rm = TRUE)
  kpn_n <- sum(cond_vec[taxons == "Klebsiella pneumoniae"], na.rm = TRUE)
  ecc_n <- sum(cond_vec[taxons == "Enterobacter cloacae complex"], na.rm = TRUE)
  tot_n <- sum(cond_vec, na.rm = TRUE)

  d_eco <- sum(taxons == "Escherichia coli")
  d_kpn <- sum(taxons == "Klebsiella pneumoniae")
  d_ecc <- sum(taxons == "Enterobacter cloacae complex")
  d_tot <- nrow(df)

  # Chi-square across species, with Fisher fallback for sparse cells
  p_str <- "\u2014"
  if (tot_n > 0) {
    m <- matrix(c(eco_n, d_eco - eco_n, kpn_n, d_kpn - kpn_n, ecc_n, d_ecc - ecc_n),
                nrow = 3, byrow = TRUE)
    p <- tryCatch(suppressWarnings(chisq.test(m)$p.value), error = function(e) NA)
    if (is.na(p) || is.nan(p)) {
      p <- tryCatch(fisher.test(m, simulate.p.value = TRUE, B = 2000)$p.value, error = function(e) NA)
    }
    if (!is.na(p) && !is.nan(p)) {
      p_str <- if (p < 0.0001) "< 0.0001" else sprintf("%.4f", p)
    }
  }

  c(
    E_coli        = sprintf("%d (%.1f%%)", eco_n, 100 * eco_n / d_eco),
    K_pneumoniae  = sprintf("%d (%.1f%%)", kpn_n, 100 * kpn_n / d_kpn),
    E_cloacae     = sprintf("%d (%.1f%%)", ecc_n, 100 * ecc_n / d_ecc),
    Overall       = sprintf("%d (%.1f%%)", tot_n, 100 * tot_n / d_tot),
    p_value       = p_str
  )
}

rows_table2 <- list(
  "Number of isolates"                          = c(E_coli = as.character(sum(taxons == "Escherichia coli")),
                                                    K_pneumoniae = as.character(sum(taxons == "Klebsiella pneumoniae")),
                                                    E_cloacae = as.character(sum(taxons == "Enterobacter cloacae complex")),
                                                    Overall = as.character(nrow(df)), p_value = "\u2014"),
  "Cefepime (R)"                                = format_table2_row(is_r("cefepime_Inter")),
  "Piperacillin-tazobactam (R)"                 = format_table2_row(is_r("piperacilline_tazobactam_Inter")),
  "Ciprofloxacin (R)"                           = format_table2_row(is_r("ciprofloxacine_Inter")),
  "Gentamicin (R)"                              = format_table2_row(is_r("gentamicine_Inter")),
  "Amikacin (R)"                                = format_table2_row(is_r("amikacine_Inter")),
  "Imipenem (R)"                                = format_table2_row(is_r("imipeneme_Inter")),
  "blaCTX-M-15"                                 = format_table2_row(has_ctx_m_15),
  "blaCTX-M-27"                                 = format_table2_row(has_ctx_m_27),
  "blaCTX-M-14"                                 = format_table2_row(has_ctx_m_14),
  "blaSHV-12"                                   = format_table2_row(has_shv_12),
  "blaOXA-1"                                    = format_table2_row(has_oxa_1),
  "blaTEM-1 / blaTEM-1B"                        = format_table2_row(has_tem_1),
  "Tripartite block (CTX-M-15 + OXA-1 + TEM-1)" = format_table2_row(has_tripartite),
  "Carbapenemase-producing (CPE)"               = format_table2_row(has_cpe),
  "GyrA QRDR mutation"                          = format_table2_row(has_gyra),
  "ParC QRDR mutation"                          = format_table2_row(has_parc),
  "PMQR (qnr determinants)"                     = format_table2_row(has_qnr),
  "aac(6')-Ib family"                           = format_table2_row(has_aac6),
  "aac(3)-II family"                            = format_table2_row(has_aac3),
  "ant(2'')-Ia family"                          = format_table2_row(has_ant2),
  "16S rRNA methyltransferase (armA/rmt)"       = format_table2_row(has_arm)
)

table2_df <- bind_rows(lapply(names(rows_table2), function(k) {
  data.frame(Characteristic = k, t(rows_table2[[k]]), stringsAsFactors = FALSE)
}))

write_csv(table2_df, file.path(out_dir_tab, "Table_2_Resistome_and_Susceptibility.csv"))
cat("--> Table 2 generated and saved successfully.\n")

# ------------------------------------------------------------------------------
# 3. DYNAMIC SUPPLEMENTARY TABLE S5 (14 ANTIMICROBIALS, S/SFP/R CATEGORIES)
# ------------------------------------------------------------------------------
abx_cols <- c(
  "amoxicilline_ac_clavulanique_Inter" = "Amoxicillin-clavulanic acid",
  "piperacilline_tazobactam_Inter"     = "Piperacillin-tazobactam",
  "temocilline_Inter"                   = "Temocillin",
  "cefoxitine_Inter"                    = "Cefoxitin",
  "ceftriaxone_Inter"                   = "Ceftriaxone",
  "ceftazidime_Inter"                   = "Ceftazidime",
  "cefepime_Inter"                      = "Cefepime",
  "imipeneme_Inter"                     = "Imipenem",
  "meropeneme_Inter"                    = "Meropenem",
  "gentamicine_Inter"                   = "Gentamicin",
  "amikacine_Inter"                     = "Amikacin",
  "ciprofloxacine_Inter"                = "Ciprofloxacin",
  "fosfomycine_Inter"                   = "Fosfomycin",
  "sulfamides_trimetoprime_Inter"       = "Trimethoprim-sulfamethoxazole"
)

s5_list <- list()
for (col_abx in names(abx_cols)) {
  label_abx <- abx_cols[col_abx]
  raw_vals  <- toupper(trimws(as.character(df[[col_abx]])))

  for (cat_sir in c("S", "SFP", "R")) {
    match_v <- (raw_vals == cat_sir)
    eco_n   <- sum(match_v[taxons == "Escherichia coli"], na.rm = TRUE)
    kpn_n   <- sum(match_v[taxons == "Klebsiella pneumoniae"], na.rm = TRUE)
    ecc_n   <- sum(match_v[taxons == "Enterobacter cloacae complex"], na.rm = TRUE)
    tot_n   <- sum(match_v, na.rm = TRUE)

    d_eco <- sum(taxons == "Escherichia coli")
    d_kpn <- sum(taxons == "Klebsiella pneumoniae")
    d_ecc <- sum(taxons == "Enterobacter cloacae complex")
    d_tot <- nrow(df)

    s5_list[[paste(label_abx, cat_sir)]] <- data.frame(
      Antibiotic        = label_abx,
      Clinical_Category = cat_sir,
      E_coli       = sprintf("%d (%.1f%%)", eco_n, 100 * eco_n / d_eco),
      K_pneumoniae = sprintf("%d (%.1f%%)", kpn_n, 100 * kpn_n / d_kpn),
      E_cloacae    = sprintf("%d (%.1f%%)", ecc_n, 100 * ecc_n / d_ecc),
      Overall      = sprintf("%d (%.1f%%)", tot_n, 100 * tot_n / d_tot),
      stringsAsFactors = FALSE
    )
  }
}

table_s5_df <- bind_rows(s5_list)
write_csv(table_s5_df, file.path(out_dir_tab, "Supplementary_Table_S5_Full_Phenotypes.csv"))
cat("--> Supplementary Table S5 generated and saved successfully.\n")
