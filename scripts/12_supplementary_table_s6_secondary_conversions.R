# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 12_supplementary_table_s6_secondary_conversions.R
# Purpose: Dynamic calculation of secondary conversion rates (true secondary
#          infection after a documented positive screening) and their
#          anatomical breakdown (Supplementary Table S6, conversion panel).
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
# Output:  output/tables/Supplementary_Table_S6_Secondary_Conversions.csv
# Author:  Francois Gravey and co-authors (see README.md)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

base_dir    <- "."
data_dir    <- file.path(base_dir, "data")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv")
df <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. CLINICAL PREPARATION OF THE CLEANED COHORT (N = 689)
# ------------------------------------------------------------------------------
# NOTE: the source database uses French labels for chronological screening
# status ("avant" = before, "apres"/"apr\u00e8s" = after, "meme jour"/"m\u00eame
# jour" = same day). These literal strings are matched below for consistency
# with the original data and must not be translated.
df_prep <- df %>%
  mutate(
    Species = case_when(
      grepl("coli", Taxon_group, ignore.case = TRUE)        ~ "Escherichia coli",
      grepl("pneumoniae", Taxon_group, ignore.case = TRUE)  ~ "Klebsiella pneumoniae",
      grepl("cloacae", Taxon_group, ignore.case = TRUE)     ~ "Enterobacter cloacae complex",
      TRUE                                                  ~ NA_character_
    ),
    inf_count     = as.numeric(as.character(Infection.BLSE)),
    has_infection = ifelse(!is.na(inf_count) & inf_count > 0, 1, 0),
    chrono_clean  = tolower(trimws(as.character(depistage.avant.apres))),

    # Documented rectal carriage (excludes the pure "always clinical" cases
    # of the Surgical Paradox)
    has_rectal_colonization = ifelse(
      chrono_clean %in% c("avant", "apres", "apr\u00e8s", "meme jour", "m\u00eame jour") | has_infection == 0,
      1, 0
    ),

    # True secondary conversion: infection occurring AFTER a positive screening
    is_true_conversion = ifelse(chrono_clean %in% c("apres", "apr\u00e8s") & has_infection == 1, 1, 0),

    nature_clean = tolower(trimws(as.character(infection.nature))),
    nature_clean = ifelse(nature_clean %in% c("-", "", "na", "nr", "unknown"), NA_character_, nature_clean)
  ) %>%
  filter(!is.na(Species))

# ------------------------------------------------------------------------------
# 2. ENUMERATION OF CONVERSIONS AND ANATOMICAL SITES
# ------------------------------------------------------------------------------
# A. Digestive carriage denominators
denom_colonized <- df_prep %>%
  filter(has_rectal_colonization == 1) %>%
  count(Species, name = "N_Colonized")

total_colonized <- sum(denom_colonized$N_Colonized)

cat(sprintf("--> Patients with documented rectal carriage: %d\n", total_colonized))
print(denom_colonized)

# B. Numerator: patients with strict secondary conversion (chrono == "apres")
conv_patients <- df_prep %>%
  filter(has_rectal_colonization == 1) %>%
  group_by(Species) %>%
  summarise(
    N_Conv_Patients = sum(is_true_conversion),
    .groups = "drop"
  ) %>%
  left_join(denom_colonized, by = "Species") %>%
  mutate(
    Pct_Conv = round((N_Conv_Patients / N_Colonized) * 100, 2)
  )

cat("\n=== TRUE SECONDARY CONVERSIONS (CHRONO = AFTER) ===\n")
print(as.data.frame(conv_patients))

# C. Anatomical breakdown of secondary conversions (events)
df_conv_sites <- df_prep %>%
  filter(is_true_conversion == 1 & !is.na(nature_clean)) %>%
  separate_rows(nature_clean, sep = "[,;/]") %>%
  mutate(nature_clean = trimws(nature_clean)) %>%
  filter(nature_clean != "") %>%
  mutate(
    Site = case_when(
      str_detect(nature_clean, "poumon|pneumo|asp|lba|trache|pavm") ~ "Ventilator-associated pneumonia (VAP)",
      str_detect(nature_clean, "bacterie|hemoc|sepsis")             ~ "Secondary bloodstream infection (BSI)",
      str_detect(nature_clean, "catheter|kt")                       ~ "Catheter-related infection",
      str_detect(nature_clean, "urin")                              ~ "Urinary tract infection (UTI)",
      str_detect(nature_clean, "abdom|adbom|periton|bile")          ~ "Intra-abdominal infection",
      str_detect(nature_clean, "plaie|pus|paroi|os|artic|peau")     ~ "Surgical site / other deep-seated",
      TRUE                                                          ~ "Other"
    )
  )

# NOTE: anatomical site keywords are matched against French free-text
# entries of the source database (e.g. "poumon" = lung, "urin" = urinary).
sites_order <- c(
  "Ventilator-associated pneumonia (VAP)",
  "Secondary bloodstream infection (BSI)",
  "Catheter-related infection",
  "Urinary tract infection (UTI)",
  "Intra-abdominal infection",
  "Surgical site / other deep-seated"
)

tab_sites_conv <- df_conv_sites %>%
  count(Species, Site, name = "N_Events") %>%
  complete(Species = c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex"),
           Site = sites_order, fill = list(N_Events = 0)) %>%
  left_join(denom_colonized, by = "Species") %>%
  mutate(
    Pct = round((N_Events / N_Colonized) * 100, 2),
    Cell_Label = sprintf("%d (%.2f%%)", N_Events, Pct)
  )

# ------------------------------------------------------------------------------
# 3. COMPILATION OF THE REVISED TABLE S6
# ------------------------------------------------------------------------------
res_table <- tab_sites_conv %>%
  select(Species, Site, Cell_Label) %>%
  pivot_wider(names_from = Species, values_from = Cell_Label) %>%
  arrange(match(Site, sites_order))

write_csv(res_table, file.path(out_dir_tab, "Supplementary_Table_S6_Secondary_Conversions.csv"))
cat(sprintf("\n--> Table S6 exported to: %s\n", out_dir_tab))
