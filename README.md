# ESBL-Producing Enterobacterales ICU Genomic Surveillance — Analysis Pipeline

This repository contains the complete, reproducible R analysis pipeline for:

> **Prospective Weekly Genomic Surveillance alerted and controlled Nosocomial Transmission Chains of ESBL-Producing Enterobacterales in Intensive Care Units: A 5-Year Longitudinal Study**

François Gravey, Camille Jeanne-Leroyer, Coraline Bernachot, Killian Le Neindre, Christophe Isnard, Clément Gakuba, Damien Du Cheyron, Pascal Thibon, Olivier Join-Lambert, Simon Le Hello.

The pipeline covers a 5-year (2019–2023) prospective genomic surveillance program of ESBL-producing *Enterobacterales* (ESBL-E) in three intensive care units (ICUs) of a French university hospital, combining weekly systematic rectal screening, whole-genome sequencing (WGS), core-genome SNP typing, and longitudinal epidemiological modelling.

## Repository structure

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── scripts/          # 17 numbered R scripts (run in order)
│   └── *.R
├── data/             # Input files (see data/README.md)
└── output/           # Generated figures and tables (created automatically)
    ├── figures/
    └── tables/
```

## Requirements

- R >= 4.3
- R packages: `dplyr`, `tidyr`, `stringr`, `readr`, `readxl`, `openxlsx`, `lubridate`, `ggplot2`, `patchwork`, `scales`, `survival`, `pROC`

Install them all at once:

```r
install.packages(c("dplyr", "tidyr", "stringr", "readr", "readxl", "openxlsx",
                   "lubridate", "ggplot2", "patchwork", "scales", "survival", "pROC"))
```

## Usage

1. Place the three input files in the `data/` folder (see `data/README.md` for the exact expected file names and formats).
2. Run the scripts in numbered order from the repository root:

```r
setwd("<path/to/repository>")
source("scripts/01_cohort_deduplication.R")
source("scripts/02_genomic_transmission_linkage.R")
# ... and so on up to script 17
```

All figures (PDF/PNG/TIFF) and tables (CSV/XLSX) are written to `output/figures/` and `output/tables/`. Intermediate cohort files produced by scripts 01–02 are written back into `data/` and consumed by all downstream scripts.

## Script description

| # | Script | Purpose | Main outputs |
|---|--------|---------|--------------|
| 01 | `01_cohort_deduplication.R` | Strict epidemiological de-duplication of the ICU cohort: only the earliest isolate per patient/stay/species/ST combination is retained. | `data/BLSE_cohort_curated_clean.tsv` |
| 02 | `02_genomic_transmission_linkage.R` | Defines intra-ward genomic cross-transmission events based on cgSNP clusters and a rolling ≤ 30-day window (nearest temporal neighbour). | `data/BLSE_cohort_final_with_transmissions.tsv` |
| 03 | `03_baseline_epidemiology.R` | Full recomputation of the baseline epidemiology results (cohort dynamics, species distribution, imported vs acquired status). | `output/tables/Results_Baseline_*.csv` |
| 04 | `04_population_structure.R` | Clonal architecture and high-risk lineages (MLST distribution, temporal dynamics of sequence types). | `output/tables/Results_PopStructure_*.csv` |
| 05 | `05_resistome_table1.R` | EUCAST phenotypes and AMRFinderPlus resistome markers; builds Table 1 (resistome & susceptibility) and Supplementary Table S2. | `output/tables/Table_1_*.csv` |
| 06 | `06_figure1.R` | Figure 1: acquisition origin (imported ≤ 48 h vs nosocomial), tripartite resistome (blaCTX-M-15 + blaOXA-1 + blaTEM-1), piperacillin-tazobactam resistance. | `output/figures/Figure_1_*.pdf` |
| 07 | `07_table2_and_figure2.R` | Table 2 (longitudinal incidence trends, rates, IRR — requires the admission numerators file) and full Figure 2. | `output/tables/Table_2_*.csv`, `output/figures/Figure_2_*.pdf` |
| 08 | `08_figure3_tipping_points_cox.R` | Figure 3: empirical ROC curves & Youden thresholds, step-function odds ratios above colonization-pressure thresholds, and the multivariable Cox forest plot. | `output/figures/Figure_3_*.pdf` |
| 09 | `09_supplementary_table_s5_colonization_pressure_roc.R` | Ward colonization pressure computation, stratified ROC curves with bootstrap-validated AUC, Youden inflection thresholds, step-function ORs with Fisher fallback (Table S5; panels 3A/3B). | `output/tables/Supplementary_Table_S5_*.csv` |
| 10 | `10_cox_transmission_velocity_pairs.R` | Pairwise transmission network reconstruction and multivariable Cox model of acute cross-transmission velocity (≤ 30 days); Figure 3C forest plot. | `output/figures/Figure_3C_*.pdf`, `output/tables/Table_Cox_*.csv` |
| 11 | `11_supplementary_tables_s4_clonal_dynamics.R` | Supplementary Tables S4a–S4d: global clonal profile summary, persistent/endemic lineages (≥ 4 years), emerging lineages, disappeared lineages. | `output/tables/Supplementary_Table_S4*.csv/.xlsx` |
| 12 | `12_supplementary_table_s6_secondary_conversions.R` | Secondary conversion rates (infection after documented positive screening) and anatomical breakdown of conversion events (Table S6, conversion panel). | `output/tables/Supplementary_Table_S6_Secondary_Conversions.csv` |
| 13 | `13_figure_s3_plasmid_architecture.R` | PlasmidFinder replicon profiles and comparative plasmid architecture across species (Figure S3, Table S3). | `output/figures/Supplementary_Figure_S3_*.pdf`, `output/tables/Supplementary_Table_S3_*.csv` |
| 14 | `14_individual_patient_trajectories.R` | Patient-level evaluation: acquisition delays, infection-site text mining, conversion rates, incident-acquisition Cox models, logistic model of clinical infection, and a native Kaplan-Meier curve with risk table. | `output/figures/Figure_Individual_Acquisition_Survival.*`, several `output/tables/*.csv` |
| 15 | `15_supplementary_table_s6_clinical_trajectories.R` | Individual clinical trajectories (5 categories) and the "Surgical Paradox" analysis (Table S6; legacy Figure 4 rendering). | `output/tables/Supplementary_Table_S6_Clinical_Trajectories.csv`, `output/figures/Figure_4_*.pdf` |
| 16 | `16_figure4_clinical_trajectories_paradox.R` | Publication version of Figure 4: carriage-to-infection trajectories by unit, genomic dissection of the Surgical Paradox (SICU), anatomical infection sites. | `output/figures/Figure_4_Clinical_Trajectories_Paradox.pdf/.png/.tiff` |
| 17 | `17_supplementary_figure_s1_temporal_saturation.R` | Non-linear temporal saturation modelling (NLS logistic fit) of genomic linkage delays and multi-window (14/21/30/60-day) sensitivity analysis (Figure S1, Table S1). | `output/figures/Supplementary_Figure_S1_*.pdf`, `output/tables/Supplementary_Table_S1_Sensitivity.csv` |

## Data availability

The patient-level cohort data supporting the findings cannot be shared publicly because they contain protected health information (French GDPR / RGPD). They are available from the corresponding author upon reasonable academic request. The `data/README.md` file documents the exact expected structure and file names so that authorized researchers can run the full pipeline.

## Notes on language and data labels

All code, comments, and console output are in English. The source database, however, contains a small number of French free-text labels (e.g. chronological screening status "avant"/"apres", anatomical infection-site keywords). Where these strings are used for matching, they are kept verbatim — with Unicode escapes (`\u00e8`, …) — and documented in the code, so that the pipeline remains fully reproducible against the original data.

## Citation

If you use this code, please cite the article and this repository (see `CITATION.cff`):

```bibtex
@article{gravey2026esbl,
  title   = {Prospective Weekly Genomic Surveillance alerted and controlled Nosocomial Transmission Chains of ESBL-Producing Enterobacterales in Intensive Care Units: A 5-Year Longitudinal Study},
  author  = {Gravey, Fran{\c{c}}ois and Jeanne-Leroyer, Camille and Bernachot, Coraline and Le Neindre, Killian and Isnard, Christophe and Gakuba, Cl{\'e}ment and Du Cheyron, Damien and Thibon, Pascal and Join-Lambert, Olivier and Le Hello, Simon},
  journal = {The Lancet Microbe},
  year    = {2026},
  note    = {Under review}
}
```

## License

This code is released under the MIT License (see `LICENSE`). If you use or adapt this pipeline in a scientific publication, please cite the article above.
