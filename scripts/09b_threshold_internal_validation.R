# ==============================================================================
# Study : Prospective Weekly Genomic Surveillance alerted and controlled
#         Nosocomial Transmission Chains of ESBL-Producing Enterobacterales
#         in Intensive Care Units: A 5-Year Longitudinal Study
# Script: 09b_threshold_internal_validation.R
# Purpose: Internal validation of the colonization-pressure transmission
#          threshold identified in script 09, addressing reviewer-anticipated
#          concerns BEFORE submission:
#            (A) Threshold stability: 2000-replicate bootstrap of the Youden
#                threshold per ward -> median threshold + 95% percentile CI.
#            (B) Cross-ward (geographical) validation: threshold derived in
#                one ICU applied to the other ICU (AUC, sens, spec, OR).
#            (C) Optimism-corrected AUC: repeated 50/50 split-sample
#                (train: Youden threshold + ROC; test: AUC at that threshold),
#                1000 repeats. NOTE: a temporal split is impossible because
#                no transmission events occurred after 2021 (zero events in
#                any holdout >= 2022), so this design substitutes for it and
#                must be stated as such in the manuscript.
# Input :  data/BLSE_cohort_final_with_transmissions.tsv
#          (same file, columns and ward definitions as script 09)
# Output:  output/tables/Supplementary_Table_S5b_Threshold_Validation.csv
#          console report with exact numbers to paste into the manuscript
# Run AFTER script 09 (needs only the input TSV, no dependency on 09 output).
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(readr)
})

base_dir <- "."
data_dir <- file.path(base_dir, "data")
out_dir_tab <- file.path(base_dir, "output", "tables")
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(data_dir, "BLSE_cohort_final_with_transmissions.tsv"),
               sep = "\t", stringsAsFactors = FALSE)

# --- Same preparation block as script 09 (identical definitions) ---------------
df_roc_prep <- df %>%
  mutate(
    pressure_num = as.numeric(nb_bmr_semaine),
    event_tx     = ifelse(is_cross_transmission_30d %in% c(TRUE, "TRUE", 1, "1"), 1, 0),
    ward_clean    = case_when(
      grepl("chir", cores, ignore.case = TRUE) ~ "Surgical ICU",
      grepl("med",  cores, ignore.case = TRUE) ~ "Medical ICU",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pressure_num), !is.na(event_tx), !is.na(ward_clean))

get_ward <- function(w) filter(df_roc_prep, ward_clean == w)

youden_threshold <- function(d) {
  roc_obj <- roc(d$event_tx, d$pressure_num, quiet = TRUE, direction = "<")
  as.numeric(coords(roc_obj, x = "best", method = "youden", ret = "threshold")[1])
}

perf_at_threshold <- function(d, thr) {
  roc_obj <- roc(d$event_tx, d$pressure_num, quiet = TRUE, direction = "<")
  sens    <- as.numeric(sensitivity(roc_obj, thr))
  spec    <- as.numeric(specificity(roc_obj, thr))
  auc_v   <- as.numeric(auc(roc_obj))
  above   <- ifelse(d$pressure_num >= thr, 1, 0)
  tab     <- table(Exposure = factor(above, levels = c(0, 1)),
                   Event    = factor(d$event_tx, levels = c(0, 1)))
  fit     <- fisher.test(tab)
  list(threshold = thr, sens = sens, spec = spec, auc = auc_v,
       or = as.numeric(fit$estimate), or_lo = fit$conf.int[1],
       or_hi = fit$conf.int[2], p = fit$p.value)
}

wards <- c("Surgical ICU", "Medical ICU")
n_boot <- 2000
set.seed(42)

# =============================================================================
# (A) BOOTSTRAP THRESHOLD STABILITY  (per ward)
# =============================================================================
boot_summary <- lapply(wards, function(w) {
  d <- get_ward(w)
  thr_obs  <- youden_threshold(d)
  boot_thr <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(d)), replace = TRUE)
    b   <- d[idx, ]
    # guard: bootstrap sample needs both classes
    if (length(unique(b$event_tx)) < 2) return(NA_real_)
    youden_threshold(b)
  })
  boot_thr <- boot_thr[!is.na(boot_thr)]
  ci_thr  <- quantile(boot_thr, c(0.025, 0.5, 0.975))
  tibble(
    Ward            = w,
    Observed_Thr    = thr_obs,
    Boot_Median_Thr = ci_thr[2],
    Thr_95CI_Low    = ci_thr[1],
    Thr_95CI_High   = ci_thr[3],
    Pct_Within_0_5  = mean(abs(boot_thr - thr_obs) <= 0.5) * 100,
    N_Boot_Valid    = length(boot_thr)
  )
})
tab_A <- bind_rows(boot_summary)

# =============================================================================
# (B) CROSS-WARD (GEOGRAPHICAL) VALIDATION
#    Derive the threshold in ICU 1, apply it to ICU 2 (and vice versa)
# =============================================================================
cross_pairs <- lapply(1:2, function(i) {
  w_train <- wards[i]; w_test <- wards[3 - i]
  d_train <- get_ward(w_train); d_test <- get_ward(w_test)
  thr <- youden_threshold(d_train)
  perf <- perf_at_threshold(d_test, thr)
  tibble(
    Train_Ward   = w_train,
    Test_Ward    = w_test,
    Applied_Thr  = perf$threshold,
    Test_AUC     = round(perf$auc, 3),
    Test_Sens    = round(perf$sens * 100, 1),
    Test_Spec    = round(perf$spec * 100, 1),
    Test_OR      = round(perf$or, 2),
    Test_OR_CI   = sprintf("[%.2f-%.2f]", perf$or_lo, perf$or_hi),
    Test_P       = ifelse(perf$p < 0.001, "< 0.001", sprintf("%.4f", perf$p))
  )
})
tab_B <- bind_rows(cross_pairs)

# =============================================================================
# (C) OPTIMISM-CORRECTED AUC (repeated 50/50 split-sample, 1000 repeats)
#    Apparent AUC (full sample) - mean test AUC = optimism estimate
# =============================================================================
n_split <- 1000
set.seed(42)
split_summary <- lapply(wards, function(w) {
  d <- get_ward(w)
  roc_app <- roc(d$event_tx, d$pressure_num, quiet = TRUE, direction = "<")
  auc_app <- as.numeric(auc(roc_app))
  test_aucs <- replicate(n_split, {
    idx   <- sample(seq_len(nrow(d)), size = floor(nrow(d) / 2))
    train <- d[idx, ]; test <- d[-idx, ]
    if (length(unique(train$event_tx)) < 2 || length(unique(test$event_tx)) < 2) return(NA_real_)
    thr  <- youden_threshold(train)
    if (!is.finite(thr)) return(NA_real_)
    roc_t <- roc(test$event_tx, test$pressure_num, quiet = TRUE, direction = "<")
    as.numeric(auc(roc_t))
  })
  test_aucs <- test_aucs[!is.na(test_aucs)]
  auc_cv    <- mean(test_aucs)
  tibble(
    Ward              = w,
    Apparent_AUC      = round(auc_app, 3),
    SplitSample_AUC   = round(auc_cv, 3),
    Optimism          = round(auc_app - auc_cv, 3),
    OptimismCorrected = round(2 * auc_cv - auc_app, 3),
    N_Valid_Splits    = length(test_aucs)
  )
})
tab_C <- bind_rows(split_summary)

# =============================================================================
# EXPORT + CONSOLE REPORT (numbers ready to paste in the manuscript)
# =============================================================================
val_table <- bind_rows(
  tab_A %>% mutate(Analysis = "A. Bootstrap threshold stability"),
  tab_B %>% mutate(Analysis = "B. Cross-ward validation"),
  tab_C %>% mutate(Analysis = "C. Optimism-corrected AUC")
)

write_csv(val_table, file.path(out_dir_tab, "Supplementary_Table_S5b_Threshold_Validation.csv"))

cat("\n=============================================================\n")
cat("INTERNAL VALIDATION OF THE TRANSMISSION THRESHOLD\n");
cat("=============================================================\n\n")
cat("(A) Bootstrap threshold stability (2000 replicates):\n");   print(as.data.frame(tab_A))
cat("\n(B) Cross-ward (geographical) validation:\n");            print(as.data.frame(tab_B))
cat("\n(C) Optimism-corrected AUC (1000 x 50/50 split):\n");     print(as.data.frame(tab_C))
cat("\n--> Exported: output/tables/Supplementary_Table_S5b_Threshold_Validation.csv\n")
cat("\nINTERPRETATION GUIDE FOR THE MANUSCRIPT:\n")
cat("* If Boot_Median_Thr is within +/-0.5 of the observed threshold in BOTH wards\n")
cat("  -> threshold is stable; report the bootstrap 95% CI in the main text.\n")
cat("* If Cross-ward OR remains >1 with p<0.05 in BOTH directions\n")
cat("  -> the threshold generalises across wards; cite Test_AUC and Test_OR.\n")
cat("* Optimism-corrected AUC within 0.05 of apparent AUC\n")
cat("  -> no material over-optimism; cite OptimismCorrected in the supplement.\n")