#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("composition")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Seurat)
  library(ggprism)
})

if (!requireNamespace("speckle", quietly = TRUE)) {
  stop("Package 'speckle' is required for propeller statistics. Install/load speckle and rerun this add-on.")
}

# =========================================================
# CONFIG
# =========================================================
base_scRNA_dir <- study_path("base_scRNA_dir")
out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

figshare_dir <- file.path(out_dir, "figshare_datasets", "scRNA_Fig2", "propeller_Fig2E")
dir.create(figshare_dir, showWarnings = FALSE, recursive = TRUE)
figshare_file <- function(...) file.path(figshare_dir, ...)

plot_dir <- file.path(out_dir, "Fig2E_propeller_outputs")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
plot_file <- function(...) file.path(plot_dir, ...)

# Original manuscript group labels used by propeller output.
propeller_group_levels <- c("HCT PBL", "CAR Marrow", "CAR PBMC", "Control Marrow", "Control PBMC")

group_order_timepoint2 <- c(
  "Untreated Control PBL",
  "Untreated Control Marrow",
  "HCT PBL",
  "CAR PBL",
  "CAR Marrow"
)

celltype_levels_propeller <- c(
  "Mature B-cell", "Erythroid Cell", "Naive T-cell", "Plasma Cell",
  "Immature B-cell", "Stromal Cell", "Other T-cell",
  "Megakaryocyte and Platelet", "Myeloid Cell", "Progenitor Cell",
  "Effector T-cell", "NK Cell", "Memory T-cell", "Other"
)

cols_propeller_group <- c(
  "HCT PBL"        = "#8DA0CB",
  "CAR Marrow"     = "#A6D854",
  "CAR PBMC"       = "#E78AC3",
  "Control Marrow" = "#FC8D62",
  "Control PBMC"   = "#66C2A5"
)

manifest <- tibble::tibble(
  file = character(), rows = integer(), columns = integer(), description = character()
)

# =========================================================
# HELPERS
# =========================================================
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

first_existing <- function(candidates, nm) {
  hit <- intersect(candidates, nm)
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

remove_internal_id_columns <- function(df) {
  bad_name <- stringr::regex("(^|_)(patient|patient_id|sample|sample_id|orig.ident|barcode|raw|mrn|cct|sid|trial_id)(_|$)", ignore_case = TRUE)
  drop_cols <- names(df)[stringr::str_detect(names(df), bad_name)]
  drop_cols <- setdiff(drop_cols, names(df)[stringr::str_detect(names(df), "^figshare_")])
  df %>% dplyr::select(-dplyr::any_of(drop_cols))
}

audit_no_disallowed_ids <- function(df, label) {
  if (nrow(df) == 0 || ncol(df) == 0) return(invisible(TRUE))
  chr_df <- df %>% dplyr::select(where(~is.character(.x) || is.factor(.x)))
  if (ncol(chr_df) == 0) return(invisible(TRUE))
  vals <- unlist(lapply(chr_df, as.character), use.names = FALSE)
  vals <- vals[!is.na(vals)]
  # Only block true raw patient IDs like [study-specific]/[study-specific], not CCT-family gene symbols.
  bad_pattern <- "(\\bCCT[0-9]{3,}(?:[-_A-Za-z0-9]*)?\\b|\\bS[0-9]{4}(?:o[0-9]+|[-_][A-Za-z0-9]+)?\\b|\\b[0-9]{7,}\\b)"
  bad_vals <- unique(vals[stringr::str_detect(vals, bad_pattern)])
  if (length(bad_vals) > 0) {
    stop(
      "Disallowed raw identifier detected in ", label, ": ",
      paste(utils::head(bad_vals, 20), collapse = ", "),
      ifelse(length(bad_vals) > 20, " ...", "")
    )
  }
  invisible(TRUE)
}

write_figshare_csv <- function(df, filename, description) {
  df_out <- df %>% remove_internal_id_columns()
  audit_no_disallowed_ids(df_out, filename)
  readr::write_csv(df_out, figshare_file(filename), na = "")
  manifest <<- dplyr::bind_rows(
    manifest,
    tibble::tibble(file = filename, rows = nrow(df_out), columns = ncol(df_out), description = description)
  )
  invisible(df_out)
}

# =========================================================
# LOAD + LIGHT ANNOTATION
# =========================================================
setwd(base_scRNA_dir)
Marrow_merge_idents <- readRDS(study_path("input_final_scrnaseq"))

# Cohort harmonization
Marrow_merge_idents$timepoint_2 <- dplyr::case_when(
  Marrow_merge_idents$timepoint == "Control PBMC"   ~ "Untreated Control PBL",
  Marrow_merge_idents$timepoint == "Control Marrow" ~ "Untreated Control Marrow",
  Marrow_merge_idents$timepoint == "Auto PBMC"      ~ "HCT PBL",
  Marrow_merge_idents$timepoint == "CAR PBMC"       ~ "CAR PBL",
  Marrow_merge_idents$timepoint == "CAR Marrow"     ~ "CAR Marrow",
  TRUE                                              ~ as.character(Marrow_merge_idents$timepoint)
)
Marrow_merge_idents$timepoint_2 <- factor(Marrow_merge_idents$timepoint_2, levels = group_order_timepoint2)

# Level 2 and level 3 cell-type annotations from the Fig 2 script.
T_cell_l2 <- c("CD8 Effector_2", "CD8 Memory", "CD14 Mono", "Macrophage", "NK", "CD8 Effector_1", "MAIT", "CD4 Memory", "CD8 Naive", "ILC", "T Proliferating", "CD8 Effector_3", "CD4 Naive", "CD4 Effector")
Myeloid_l2 <- c("CD14 Mono", "pre-mDC", "Macrophage", "CD16 Mono", "cDC1", "cDC2", "pDC", "pre-pDC", "ASDC")
Progenitor_l2 <- c("HSC", "BaEoMa", "GMP", "EMP", "LMPP", "CLP")
NK_l2 <- c("NK", "ILC", "NK CD56+", "NK Proliferating")
Stromal_l2 <- c("Stromal")
Plasma_l2 <- c("Plasma")
MegaPlt_l2 <- c("Prog Mk", "Platelet")
B_l2 <- c("Memory B", "Naive B", "transitional B", "pre B", "pro B")
Eryth_l2 <- c("Late Eryth", "Early Eryth")

Marrow_merge_idents$predicted.celltype_2 <- "Other"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% B_l2] <- "B-cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% T_cell_l2] <- "T-cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% NK_l2] <- "NK Cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% Myeloid_l2] <- "Myeloid Cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% Progenitor_l2] <- "Progenitor Cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% Stromal_l2] <- "Stromal Cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% Plasma_l2] <- "Plasma Cell"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% MegaPlt_l2] <- "Megakaryocyte and Platelet"
Marrow_merge_idents$predicted.celltype_2[Marrow_merge_idents$predicted.celltype.l2 %in% Eryth_l2] <- "Erythroid Cell"

T_effector_l2 <- c("CD8 Effector_1", "CD8 Effector_2", "CD8 Effector_3", "CD4 Effector")
T_naive_l2 <- c("CD8 Naive", "CD4 Naive")
T_memory_l2 <- c("CD8 Memory", "CD4 Memory")
T_other_l2 <- c("MAIT", "T Proliferating")
B_mature_l2 <- c("Memory B", "Naive B")
B_immature_l2 <- c("transitional B", "pre B", "pro B")

Marrow_merge_idents$predicted.celltype_3 <- as.character(Marrow_merge_idents$predicted.celltype_2)
Marrow_merge_idents$predicted.celltype_3[Marrow_merge_idents$predicted.celltype.l2 %in% B_mature_l2] <- "Mature B-cell"
Marrow_merge_idents$predicted.celltype_3[Marrow_merge_idents$predicted.celltype.l2 %in% B_immature_l2] <- "Immature B-cell"
Marrow_merge_idents$predicted.celltype_3[Marrow_merge_idents$predicted.celltype.l2 %in% T_effector_l2] <- "Effector T-cell"
Marrow_merge_idents$predicted.celltype_3[Marrow_merge_idents$predicted.celltype.l2 %in% T_naive_l2] <- "Naive T-cell"
Marrow_merge_idents$predicted.celltype_3[Marrow_merge_idents$predicted.celltype.l2 %in% T_memory_l2] <- "Memory T-cell"
Marrow_merge_idents$predicted.celltype_3[Marrow_merge_idents$predicted.celltype.l2 %in% T_other_l2] <- "Other T-cell"
Marrow_merge_idents$predicted.celltype_3 <- factor(Marrow_merge_idents$predicted.celltype_3, levels = celltype_levels_propeller)

# Deidentified maps
meta_raw <- Marrow_merge_idents@meta.data %>%
  tibble::rownames_to_column("raw_cell_id") %>%
  dplyr::mutate(
    patient_id = as.character(patient_id),
    sample_id  = as.character(sample_id)
  )

patient_map <- meta_raw %>%
  dplyr::filter(!is.na(patient_id)) %>%
  dplyr::distinct(raw_patient_id = patient_id) %>%
  dplyr::arrange(raw_patient_id) %>%
  dplyr::mutate(figshare_patient_id = paste0("scRNA_Patient_", sprintf("%03d", dplyr::row_number())))

sample_map <- meta_raw %>%
  dplyr::filter(!is.na(sample_id)) %>%
  dplyr::distinct(raw_sample_id = sample_id) %>%
  dplyr::arrange(raw_sample_id) %>%
  dplyr::mutate(figshare_sample_id = paste0("scRNA_Sample_", sprintf("%03d", dplyr::row_number())))

cell_propeller_public <- meta_raw %>%
  dplyr::left_join(patient_map, by = c("patient_id" = "raw_patient_id")) %>%
  dplyr::left_join(sample_map, by = c("sample_id" = "raw_sample_id")) %>%
  dplyr::transmute(
    figshare_patient_id,
    figshare_sample_id,
    timepoint_2 = as.character(timepoint_2),
    propeller_group = dplyr::case_when(
      timepoint_2 == "Untreated Control PBL"    ~ "Control PBMC",
      timepoint_2 == "Untreated Control Marrow" ~ "Control Marrow",
      timepoint_2 == "HCT PBL"                  ~ "HCT PBL",
      timepoint_2 == "CAR PBL"                  ~ "CAR PBMC",
      timepoint_2 == "CAR Marrow"               ~ "CAR Marrow",
      TRUE                                       ~ NA_character_
    ),
    source_group = dplyr::case_when(
      timepoint_2 %in% c("Untreated Control Marrow", "CAR Marrow") ~ "Marrow",
      timepoint_2 %in% c("Untreated Control PBL", "CAR PBL", "HCT PBL") ~ "PBL",
      TRUE ~ NA_character_
    ),
    car_control_group = dplyr::case_when(
      timepoint_2 %in% c("CAR Marrow", "CAR PBL") ~ "CAR",
      timepoint_2 %in% c("Untreated Control Marrow", "Untreated Control PBL") ~ "Control",
      TRUE ~ NA_character_
    ),
    predicted.celltype_2 = as.character(predicted.celltype_2),
    predicted.celltype_3 = as.character(predicted.celltype_3)
  ) %>%
  dplyr::filter(!is.na(figshare_sample_id), !is.na(propeller_group), !is.na(predicted.celltype_3)) %>%
  dplyr::mutate(
    propeller_group = factor(propeller_group, levels = propeller_group_levels),
    predicted.celltype_3 = factor(predicted.celltype_3, levels = celltype_levels_propeller)
  ) %>%
  droplevels()

# =========================================================
# PROPELLER INPUTS, STATS, AND PLOTS
# =========================================================
propeller_sample_fractions <- cell_propeller_public %>%
  dplyr::count(figshare_patient_id, figshare_sample_id, timepoint_2, propeller_group, source_group, car_control_group, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(figshare_patient_id, figshare_sample_id, timepoint_2, propeller_group, source_group, car_control_group) %>%
  tidyr::complete(
    predicted.celltype_3 = celltype_levels_propeller,
    fill = list(n_cells = 0L)
  ) %>%
  dplyr::mutate(
    total_cells = sum(n_cells, na.rm = TRUE),
    fraction = ifelse(total_cells > 0, n_cells / total_cells, NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(propeller_group), !is.na(predicted.celltype_3))

propeller_sample_summary <- propeller_sample_fractions %>%
  dplyr::group_by(propeller_group, timepoint_2, predicted.celltype_3) %>%
  dplyr::summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    se_fraction = se(fraction),
    n_samples = dplyr::n_distinct(figshare_sample_id),
    .groups = "drop"
  )

# Level 2 fractions are used for the manuscript statement that CAR+ patients globally lacked total B cells.
celltype_levels_lvl2 <- c(
  "T-cell", "Myeloid Cell", "NK Cell", "Progenitor Cell", "B-cell",
  "Erythroid Cell", "Megakaryocyte and Platelet", "Plasma Cell",
  "Stromal Cell", "Other"
)

propeller_sample_fractions_lvl2 <- cell_propeller_public %>%
  dplyr::count(figshare_patient_id, figshare_sample_id, timepoint_2, propeller_group, source_group, car_control_group, predicted.celltype_2, name = "n_cells") %>%
  dplyr::group_by(figshare_patient_id, figshare_sample_id, timepoint_2, propeller_group, source_group, car_control_group) %>%
  tidyr::complete(
    predicted.celltype_2 = celltype_levels_lvl2,
    fill = list(n_cells = 0L)
  ) %>%
  dplyr::mutate(
    total_cells = sum(n_cells, na.rm = TRUE),
    fraction = ifelse(total_cells > 0, n_cells / total_cells, NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(propeller_group), !is.na(predicted.celltype_2))

propeller_sample_summary_lvl2 <- propeller_sample_fractions_lvl2 %>%
  dplyr::group_by(propeller_group, timepoint_2, predicted.celltype_2) %>%
  dplyr::summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    se_fraction = se(fraction),
    n_samples = dplyr::n_distinct(figshare_sample_id),
    .groups = "drop"
  )

# =========================================================
# MANUSCRIPT CLAIM-SPECIFIC BETA-BINOMIAL TESTS
# =========================================================
# These are the direct tests for the Fig 2E red/purple asterisks described in the text:
#   Red asterisk: marrow capture/source effect (Marrow vs PBL), especially erythroid,
#                 immature B-cell, and progenitor cells.
#   Purple asterisk: CAR-patient global composition effect (CAR vs untreated control,
#                    excluding HCT PBL), especially B-cell loss, reduced naive T cells,
#                    and expanded effector T cells.
# The global speckle::propeller ANOVA table is still exported above/below; these models
# provide targeted two-group beta-binomial tests with BH/FDR correction.

safe_coef <- function(fit, prefix = "test_group") {
  co <- tryCatch(summary(fit)$coefficients$cond, error = function(e) NULL)
  if (is.null(co)) return(tibble::tibble(log_odds_ratio = NA_real_, se = NA_real_, z = NA_real_, p_wald = NA_real_))
  rn <- rownames(co)
  ix <- grep(paste0("^", prefix), rn)
  if (length(ix) == 0) return(tibble::tibble(log_odds_ratio = NA_real_, se = NA_real_, z = NA_real_, p_wald = NA_real_))
  tibble::tibble(
    log_odds_ratio = as.numeric(co[ix[1], "Estimate"]),
    se = as.numeric(co[ix[1], "Std. Error"]),
    z = as.numeric(co[ix[1], "z value"]),
    p_wald = as.numeric(co[ix[1], "Pr(>|z|)"])
  )
}

run_beta_binomial_comparison <- function(fraction_df,
                                         celltype_col,
                                         grouping_col,
                                         ref_level,
                                         alt_level,
                                         analysis_name,
                                         analysis_claim,
                                         target_celltypes = NULL) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required for the manuscript claim beta-binomial tests. Install glmmTMB and rerun this add-on.")
  }

  dat0 <- fraction_df %>%
    dplyr::filter(.data[[grouping_col]] %in% c(ref_level, alt_level)) %>%
    dplyr::mutate(
      test_group = factor(as.character(.data[[grouping_col]]), levels = c(ref_level, alt_level)),
      celltype = as.character(.data[[celltype_col]]),
      failures = total_cells - n_cells
    ) %>%
    dplyr::filter(!is.na(test_group), !is.na(celltype), !is.na(n_cells), !is.na(total_cells), total_cells > 0, failures >= 0)

  if (!is.null(target_celltypes)) dat0 <- dat0 %>% dplyr::filter(celltype %in% target_celltypes)

  out <- purrr::map_dfr(sort(unique(dat0$celltype)), function(ct) {
    dd <- dat0 %>% dplyr::filter(celltype == ct)

    mean_tbl <- dd %>%
      dplyr::group_by(test_group) %>%
      dplyr::summarise(
        mean_fraction = mean(fraction, na.rm = TRUE),
        median_fraction = median(fraction, na.rm = TRUE),
        n_samples = dplyr::n_distinct(figshare_sample_id),
        total_cells = sum(total_cells, na.rm = TRUE),
        total_celltype_cells = sum(n_cells, na.rm = TRUE),
        .groups = "drop"
      )

    get_mean <- function(level, col) {
      val <- mean_tbl %>% dplyr::filter(test_group == level) %>% dplyr::pull({{ col }})
      if (length(val) == 0) NA_real_ else as.numeric(val[[1]])
    }

    base <- tibble::tibble(
      analysis_name = analysis_name,
      analysis_claim = analysis_claim,
      celltype_level = celltype_col,
      celltype = ct,
      ref_group = ref_level,
      alt_group = alt_level,
      n_samples_ref = get_mean(ref_level, n_samples),
      n_samples_alt = get_mean(alt_level, n_samples),
      mean_fraction_ref = get_mean(ref_level, mean_fraction),
      mean_fraction_alt = get_mean(alt_level, mean_fraction),
      median_fraction_ref = get_mean(ref_level, median_fraction),
      median_fraction_alt = get_mean(alt_level, median_fraction),
      delta_mean_alt_minus_ref = mean_fraction_alt - mean_fraction_ref
    )

    if (nrow(dd) < 4 || dplyr::n_distinct(dd$test_group) < 2) {
      return(base %>% dplyr::mutate(log_odds_ratio = NA_real_, odds_ratio = NA_real_, se = NA_real_, z = NA_real_, p_wald = NA_real_, p_value = NA_real_, model_status = "insufficient_samples"))
    }

    fit_res <- tryCatch({
      fit0 <- glmmTMB::glmmTMB(cbind(n_cells, failures) ~ 1, data = dd, family = glmmTMB::betabinomial(link = "logit"))
      fit1 <- glmmTMB::glmmTMB(cbind(n_cells, failures) ~ test_group, data = dd, family = glmmTMB::betabinomial(link = "logit"))
      lrt <- stats::anova(fit0, fit1)
      p_lrt <- suppressWarnings(as.numeric(lrt$`Pr(>Chisq)`[2]))
      safe_coef(fit1) %>%
        dplyr::mutate(p_value = p_lrt, model_status = "ok")
    }, error = function(e) {
      tibble::tibble(log_odds_ratio = NA_real_, se = NA_real_, z = NA_real_, p_wald = NA_real_, p_value = NA_real_, model_status = paste0("failed: ", conditionMessage(e)))
    })

    base %>%
      dplyr::bind_cols(fit_res) %>%
      dplyr::mutate(odds_ratio = exp(log_odds_ratio), .after = log_odds_ratio)
  }) %>%
    dplyr::group_by(analysis_name) %>%
    dplyr::mutate(
      FDR = stats::p.adjust(p_value, method = "BH"),
      p_label = fmt_p(p_value),
      FDR_label = fmt_p(FDR),
      direction = dplyr::case_when(
        is.na(delta_mean_alt_minus_ref) ~ NA_character_,
        delta_mean_alt_minus_ref > 0 ~ paste0("Higher in ", alt_level),
        delta_mean_alt_minus_ref < 0 ~ paste0("Lower in ", alt_level),
        TRUE ~ "No mean difference"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(analysis_name, FDR, p_value)

  out
}

# Red asterisk: tissue/source effect. Include all samples, plus control-only and CAR-only
# sensitivity analyses so the red marrow-capture claim is fully traceable.
red_source_all_lvl3 <- run_beta_binomial_comparison(
  propeller_sample_fractions,
  celltype_col = "predicted.celltype_3",
  grouping_col = "source_group",
  ref_level = "PBL",
  alt_level = "Marrow",
  analysis_name = "source_effect_all_samples_marrow_vs_PBL",
  analysis_claim = "Fig 2E red asterisk: marrow specimens enriched for marrow-resident populations",
  target_celltypes = c("Erythroid Cell", "Immature B-cell", "Progenitor Cell")
)

red_source_control_lvl3 <- run_beta_binomial_comparison(
  propeller_sample_fractions %>% dplyr::filter(car_control_group == "Control"),
  celltype_col = "predicted.celltype_3",
  grouping_col = "source_group",
  ref_level = "PBL",
  alt_level = "Marrow",
  analysis_name = "source_effect_control_samples_marrow_vs_PBL",
  analysis_claim = "Fig 2E red asterisk sensitivity: untreated control marrow versus untreated control PBL",
  target_celltypes = c("Erythroid Cell", "Immature B-cell", "Progenitor Cell")
)

red_source_car_lvl3 <- run_beta_binomial_comparison(
  propeller_sample_fractions %>% dplyr::filter(car_control_group == "CAR"),
  celltype_col = "predicted.celltype_3",
  grouping_col = "source_group",
  ref_level = "PBL",
  alt_level = "Marrow",
  analysis_name = "source_effect_CAR_samples_marrow_vs_PBL",
  analysis_claim = "Fig 2E red asterisk sensitivity: CAR marrow versus CAR PBL",
  target_celltypes = c("Erythroid Cell", "Immature B-cell", "Progenitor Cell")
)

# Purple asterisk: global CAR-patient effect, excluding HCT PBL. Level 2 supports total B-cell loss;
# Level 3 supports mature/immature B-cell, naive T-cell, and effector T-cell statements.
purple_car_control_lvl2 <- run_beta_binomial_comparison(
  propeller_sample_fractions_lvl2 %>% dplyr::filter(car_control_group %in% c("Control", "CAR")),
  celltype_col = "predicted.celltype_2",
  grouping_col = "car_control_group",
  ref_level = "Control",
  alt_level = "CAR",
  analysis_name = "CAR_vs_control_global_level2",
  analysis_claim = "Fig 2E purple asterisk: CAR patients globally lack total B cells",
  target_celltypes = c("B-cell")
)

purple_car_control_lvl3 <- run_beta_binomial_comparison(
  propeller_sample_fractions %>% dplyr::filter(car_control_group %in% c("Control", "CAR")),
  celltype_col = "predicted.celltype_3",
  grouping_col = "car_control_group",
  ref_level = "Control",
  alt_level = "CAR",
  analysis_name = "CAR_vs_control_global_level3",
  analysis_claim = "Fig 2E purple asterisk: CAR patients have reduced B/naive T cells and expanded effector T cells",
  target_celltypes = c("Mature B-cell", "Immature B-cell", "Naive T-cell", "Effector T-cell")
)

paper_claim_beta_binomial_stats <- dplyr::bind_rows(
  red_source_all_lvl3,
  red_source_control_lvl3,
  red_source_car_lvl3,
  purple_car_control_lvl2,
  purple_car_control_lvl3
) %>%
  dplyr::mutate(
    asterisk_color = dplyr::case_when(
      grepl("red asterisk", analysis_claim, ignore.case = TRUE) ~ "red",
      grepl("purple asterisk", analysis_claim, ignore.case = TRUE) ~ "purple",
      TRUE ~ NA_character_
    ),
    significant_FDR_0.05 = !is.na(FDR) & FDR < 0.05,
    significant_FDR_0.10 = !is.na(FDR) & FDR < 0.10
  ) %>%
  dplyr::select(
    asterisk_color, analysis_name, analysis_claim, celltype_level, celltype,
    ref_group, alt_group, n_samples_ref, n_samples_alt,
    mean_fraction_ref, mean_fraction_alt, delta_mean_alt_minus_ref,
    odds_ratio, log_odds_ratio, se, z, p_wald, p_value, FDR,
    p_label, FDR_label, direction, significant_FDR_0.05, significant_FDR_0.10,
    model_status
  )

# Compact manuscript table with only the primary tests named in the paper text.
paper_text_primary_stats <- paper_claim_beta_binomial_stats %>%
  dplyr::filter(
    analysis_name %in% c("source_effect_all_samples_marrow_vs_PBL", "CAR_vs_control_global_level2", "CAR_vs_control_global_level3")
  ) %>%
  dplyr::arrange(factor(asterisk_color, levels = c("red", "purple")), celltype)


propeller_statistics <- tryCatch({
  speckle_out <- speckle::propeller(
    clusters = cell_propeller_public$predicted.celltype_3,
    sample   = cell_propeller_public$figshare_sample_id,
    group    = cell_propeller_public$propeller_group
  )
  as.data.frame(speckle_out) %>%
    tibble::rownames_to_column(var = "predicted.celltype_3") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      P.Value = if ("P.Value" %in% names(.)) P.Value else NA_real_,
      FDR = if ("FDR" %in% names(.)) FDR else NA_real_,
      p_label = fmt_p(P.Value),
      FDR_label = fmt_p(FDR)
    )
}, error = function(e) {
  stop("speckle::propeller failed: ", conditionMessage(e))
})

propeller_input_marrow <- cell_propeller_public %>%
  dplyr::filter(propeller_group %in% c("Control Marrow", "CAR Marrow")) %>%
  dplyr::mutate(propeller_group = factor(propeller_group, levels = c("Control Marrow", "CAR Marrow"))) %>%
  droplevels()

propeller_statistics_marrow <- tryCatch({
  speckle_out <- speckle::propeller(
    clusters = propeller_input_marrow$predicted.celltype_3,
    sample   = propeller_input_marrow$figshare_sample_id,
    group    = propeller_input_marrow$propeller_group
  )
  as.data.frame(speckle_out) %>%
    tibble::rownames_to_column(var = "predicted.celltype_3") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      P.Value = if ("P.Value" %in% names(.)) P.Value else NA_real_,
      FDR = if ("FDR" %in% names(.)) FDR else NA_real_,
      p_label = fmt_p(P.Value),
      FDR_label = fmt_p(FDR)
    )
}, error = function(e) {
  warning("Marrow-only speckle::propeller failed: ", conditionMessage(e))
  tibble::tibble()
})

write_figshare_csv(
  propeller_sample_fractions,
  "figshare_scRNA_Fig2_propeller_per_sample_celltype_fractions.csv",
  "Per-sample cell-type fractions used as input for Fig 2E speckle propeller composition analysis."
)
write_figshare_csv(
  propeller_sample_summary,
  "figshare_scRNA_Fig2_propeller_celltype_fraction_summary.csv",
  "Mean, median, and SE cell-type fractions by cohort for Fig 2E composition plots."
)
write_figshare_csv(
  propeller_sample_fractions_lvl2,
  "figshare_scRNA_Fig2_propeller_per_sample_level2_celltype_fractions.csv",
  "Per-sample Level 2 cell-type fractions used for total B-cell Fig 2E CAR-versus-control beta-binomial tests."
)
write_figshare_csv(
  propeller_sample_summary_lvl2,
  "figshare_scRNA_Fig2_propeller_level2_celltype_fraction_summary.csv",
  "Mean, median, and SE Level 2 cell-type fractions by cohort for Fig 2E total B-cell composition analysis."
)
write_figshare_csv(
  paper_claim_beta_binomial_stats,
  "figshare_scRNA_Fig2_paper_text_red_purple_asterisk_beta_binomial_statistics.csv",
  "Targeted beta-binomial tests for the Figure 2E manuscript claims: red asterisk marrow-vs-PBL source effect and purple asterisk CAR-vs-control global composition effect, with BH/FDR correction."
)
write_figshare_csv(
  paper_text_primary_stats,
  "figshare_scRNA_Fig2_paper_text_primary_beta_binomial_statistics.csv",
  "Compact primary beta-binomial statistics table for the specific Figure 2E claims stated in the manuscript text."
)
write_figshare_csv(
  propeller_statistics,
  "figshare_scRNA_Fig2_propeller_statistics_all_cohorts.csv",
  "speckle::propeller statistics for Fig 2E cell-type composition across HCT PBL, CAR marrow/PBMC, and control marrow/PBMC cohorts."
)
write_figshare_csv(
  propeller_statistics_marrow,
  "figshare_scRNA_Fig2_propeller_statistics_marrow_only.csv",
  "speckle::propeller statistics for marrow-only CAR marrow versus untreated control marrow composition."
)


# =========================================================
# STACKED BAR PLOT DATA + GLOBAL MARROW-vs-PBL COMPARISONS
# =========================================================
# These outputs restore the stacked Fig 2E composition panels and add a global
# Marrow-vs-PBL comparison across all cell types, not only the red-asterisk targets.

cols_lvl3_stacked <- c(
  "Myeloid Cell"               = "#3288BD",
  "NK Cell"                    = "#D8A499",
  "Effector T-cell"            = "#9E0142",
  "Memory T-cell"              = "#FB9A99",
  "Naive T-cell"               = "#F98400",
  "Progenitor Cell"            = "#FFFFBF",
  "Immature B-cell"            = "#5E4FA2",
  "Mature B-cell"              = "#CAB2D6",
  "Other T-cell"               = "#D53E4F",
  "Megakaryocyte and Platelet" = "#ABDDA4",
  "Erythroid Cell"             = "#66C2A5",
  "Stromal Cell"               = "#FDAE61",
  "Plasma Cell"                = "#E6F598",
  "Other"                      = "grey80"
)

# Pooled all-cell stacked composition by cohort. This matches the classic
# stacked-bar visual where every cell contributes once to its cohort total.
stacked_pooled_cohort <- cell_propeller_public %>%
  dplyr::count(propeller_group, timepoint_2, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(propeller_group, timepoint_2) %>%
  tidyr::complete(predicted.celltype_3 = celltype_levels_propeller, fill = list(n_cells = 0L)) %>%
  dplyr::mutate(
    total_cells = sum(n_cells, na.rm = TRUE),
    fraction = ifelse(total_cells > 0, n_cells / total_cells, NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    propeller_group = factor(propeller_group, levels = propeller_group_levels),
    predicted.celltype_3 = factor(predicted.celltype_3, levels = celltype_levels_propeller)
  )

# Mean-of-samples stacked composition by cohort. This is the patient/sample-level
# complement to the pooled stacked bar and matches the denominator used by boxplots.
stacked_mean_sample_cohort <- propeller_sample_fractions %>%
  dplyr::group_by(propeller_group, timepoint_2, predicted.celltype_3) %>%
  dplyr::summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    se_fraction = se(fraction),
    n_samples = dplyr::n_distinct(figshare_sample_id),
    .groups = "drop"
  ) %>%
  dplyr::group_by(propeller_group, timepoint_2) %>%
  dplyr::mutate(mean_fraction_rescaled_to_sum_1 = mean_fraction / sum(mean_fraction, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    propeller_group = factor(propeller_group, levels = propeller_group_levels),
    predicted.celltype_3 = factor(predicted.celltype_3, levels = celltype_levels_propeller)
  )

# Pooled and mean-of-samples stacked composition by global source group.
stacked_pooled_source <- cell_propeller_public %>%
  dplyr::filter(source_group %in% c("PBL", "Marrow")) %>%
  dplyr::count(source_group, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(source_group) %>%
  tidyr::complete(predicted.celltype_3 = celltype_levels_propeller, fill = list(n_cells = 0L)) %>%
  dplyr::mutate(
    total_cells = sum(n_cells, na.rm = TRUE),
    fraction = ifelse(total_cells > 0, n_cells / total_cells, NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    source_group = factor(source_group, levels = c("PBL", "Marrow")),
    predicted.celltype_3 = factor(predicted.celltype_3, levels = celltype_levels_propeller)
  )

stacked_mean_sample_source <- propeller_sample_fractions %>%
  dplyr::filter(source_group %in% c("PBL", "Marrow")) %>%
  dplyr::group_by(source_group, predicted.celltype_3) %>%
  dplyr::summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    se_fraction = se(fraction),
    n_samples = dplyr::n_distinct(figshare_sample_id),
    .groups = "drop"
  ) %>%
  dplyr::group_by(source_group) %>%
  dplyr::mutate(mean_fraction_rescaled_to_sum_1 = mean_fraction / sum(mean_fraction, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    source_group = factor(source_group, levels = c("PBL", "Marrow")),
    predicted.celltype_3 = factor(predicted.celltype_3, levels = celltype_levels_propeller)
  )

# Global marrow-vs-PBL beta-binomial tests across every Level 3 and Level 2 cell type.
global_source_all_lvl3 <- run_beta_binomial_comparison(
  propeller_sample_fractions,
  celltype_col = "predicted.celltype_3",
  grouping_col = "source_group",
  ref_level = "PBL",
  alt_level = "Marrow",
  analysis_name = "global_source_effect_all_celltypes_level3_marrow_vs_PBL",
  analysis_claim = "Global Fig 2E source comparison: marrow versus PBL across all Level 3 cell types",
  target_celltypes = NULL
) %>%
  dplyr::mutate(significant_FDR_0.05 = !is.na(FDR) & FDR < 0.05)

global_source_all_lvl2 <- run_beta_binomial_comparison(
  propeller_sample_fractions_lvl2,
  celltype_col = "predicted.celltype_2",
  grouping_col = "source_group",
  ref_level = "PBL",
  alt_level = "Marrow",
  analysis_name = "global_source_effect_all_celltypes_level2_marrow_vs_PBL",
  analysis_claim = "Global Fig 2E source comparison: marrow versus PBL across all Level 2 cell types",
  target_celltypes = NULL
) %>%
  dplyr::mutate(significant_FDR_0.05 = !is.na(FDR) & FDR < 0.05)

propeller_input_source <- cell_propeller_public %>%
  dplyr::filter(source_group %in% c("PBL", "Marrow")) %>%
  dplyr::mutate(source_group = factor(source_group, levels = c("PBL", "Marrow"))) %>%
  droplevels()

propeller_statistics_global_source <- tryCatch({
  speckle_out <- speckle::propeller(
    clusters = propeller_input_source$predicted.celltype_3,
    sample   = propeller_input_source$figshare_sample_id,
    group    = propeller_input_source$source_group
  )
  as.data.frame(speckle_out) %>%
    tibble::rownames_to_column(var = "predicted.celltype_3") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      P.Value = if ("P.Value" %in% names(.)) P.Value else NA_real_,
      FDR = if ("FDR" %in% names(.)) FDR else NA_real_,
      p_label = fmt_p(P.Value),
      FDR_label = fmt_p(FDR)
    )
}, error = function(e) {
  warning("Global source speckle::propeller failed: ", conditionMessage(e))
  tibble::tibble()
})

write_figshare_csv(
  stacked_pooled_cohort,
  "figshare_scRNA_Fig2_stacked_bar_pooled_cohort_composition.csv",
  "Pooled all-cell Level 3 cell-type fractions by cohort for the Fig 2E stacked bar plot."
)
write_figshare_csv(
  stacked_mean_sample_cohort,
  "figshare_scRNA_Fig2_stacked_bar_mean_sample_cohort_composition.csv",
  "Mean-of-samples Level 3 cell-type fractions by cohort for the Fig 2E stacked bar plot."
)
write_figshare_csv(
  stacked_pooled_source,
  "figshare_scRNA_Fig2_stacked_bar_pooled_global_marrow_vs_PBL_composition.csv",
  "Pooled all-cell Level 3 cell-type fractions for the global marrow-versus-PBL stacked bar plot."
)
write_figshare_csv(
  stacked_mean_sample_source,
  "figshare_scRNA_Fig2_stacked_bar_mean_sample_global_marrow_vs_PBL_composition.csv",
  "Mean-of-samples Level 3 cell-type fractions for the global marrow-versus-PBL stacked bar plot."
)
write_figshare_csv(
  global_source_all_lvl3,
  "figshare_scRNA_Fig2_global_marrow_vs_PBL_beta_binomial_statistics_level3.csv",
  "Global beta-binomial marrow-versus-PBL comparison across all Level 3 cell types with BH/FDR correction."
)
write_figshare_csv(
  global_source_all_lvl2,
  "figshare_scRNA_Fig2_global_marrow_vs_PBL_beta_binomial_statistics_level2.csv",
  "Global beta-binomial marrow-versus-PBL comparison across all Level 2 cell types with BH/FDR correction."
)
write_figshare_csv(
  propeller_statistics_global_source,
  "figshare_scRNA_Fig2_propeller_statistics_global_marrow_vs_PBL.csv",
  "speckle::propeller statistics for the global marrow-versus-PBL source comparison across Level 3 cell types."
)

# Stacked bar plots. Both pooled and mean-of-samples versions are saved because the pooled
# plot is the classic visual, while mean-of-samples matches the sample-level boxplot denominator.
p_stacked_pooled_cohort <- ggplot(
  stacked_pooled_cohort %>% dplyr::filter(predicted.celltype_3 != "Other"),
  aes(x = propeller_group, y = fraction, fill = predicted.celltype_3)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols_lvl3_stacked, name = "Cell Type") +
  labs(x = NULL, y = "Fraction of cells", title = "Fig 2E stacked composition by cohort") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_stacked_mean_cohort <- ggplot(
  stacked_mean_sample_cohort %>% dplyr::filter(predicted.celltype_3 != "Other"),
  aes(x = propeller_group, y = mean_fraction_rescaled_to_sum_1, fill = predicted.celltype_3)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols_lvl3_stacked, name = "Cell Type") +
  labs(x = NULL, y = "Mean fraction of cells", title = "Fig 2E mean sample-level stacked composition by cohort") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_stacked_pooled_source <- ggplot(
  stacked_pooled_source %>% dplyr::filter(predicted.celltype_3 != "Other"),
  aes(x = source_group, y = fraction, fill = predicted.celltype_3)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols_lvl3_stacked, name = "Cell Type") +
  labs(x = NULL, y = "Fraction of cells", title = "Global source composition: marrow versus PBL") +
  theme_classic(base_size = 18)

p_stacked_mean_source <- ggplot(
  stacked_mean_sample_source %>% dplyr::filter(predicted.celltype_3 != "Other"),
  aes(x = source_group, y = mean_fraction_rescaled_to_sum_1, fill = predicted.celltype_3)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols_lvl3_stacked, name = "Cell Type") +
  labs(x = NULL, y = "Mean fraction of cells", title = "Mean sample-level source composition: marrow versus PBL") +
  theme_classic(base_size = 18)

for (fmt in c("pdf", "eps", "png")) {
  dev <- switch(fmt, pdf = cairo_pdf, eps = cairo_ps, png = NULL)
  dpi_arg <- if (fmt == "png") list(dpi = 400) else list(device = dev)
  do.call(ggsave, c(list(filename = plot_file(paste0("Fig2E_stacked_bar_pooled_cohort_composition.", fmt)), plot = p_stacked_pooled_cohort, width = 11, height = 7, units = "in"), dpi_arg))
  do.call(ggsave, c(list(filename = plot_file(paste0("Fig2E_stacked_bar_mean_sample_cohort_composition.", fmt)), plot = p_stacked_mean_cohort, width = 11, height = 7, units = "in"), dpi_arg))
  do.call(ggsave, c(list(filename = plot_file(paste0("Fig2E_stacked_bar_pooled_global_marrow_vs_PBL.", fmt)), plot = p_stacked_pooled_source, width = 6.5, height = 6, units = "in"), dpi_arg))
  do.call(ggsave, c(list(filename = plot_file(paste0("Fig2E_stacked_bar_mean_sample_global_marrow_vs_PBL.", fmt)), plot = p_stacked_mean_source, width = 6.5, height = 6, units = "in"), dpi_arg))
}

p_propeller_box <- ggplot(
  propeller_sample_fractions %>% dplyr::filter(predicted.celltype_3 != "Other"),
  aes(x = propeller_group, y = fraction, fill = propeller_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.2, alpha = 0.7) +
  facet_wrap(~ predicted.celltype_3, ncol = 4, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = cols_propeller_group, guide = "none") +
  labs(x = NULL, y = "Fraction of cells", title = "Fig 2E propeller composition input") +
  theme_prism(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)
  )

ggsave(plot_file("Fig2E_propeller_celltype_fraction_boxplots.pdf"), p_propeller_box, width = 14, height = 10, units = "in", device = cairo_pdf)
ggsave(plot_file("Fig2E_propeller_celltype_fraction_boxplots.eps"), p_propeller_box, width = 14, height = 10, units = "in", device = cairo_ps)
ggsave(plot_file("Fig2E_propeller_celltype_fraction_boxplots.png"), p_propeller_box, width = 14, height = 10, units = "in", dpi = 400)

# Targeted boxplots for the red/purple asterisk tests in the manuscript text.
p_red_source <- ggplot(
  propeller_sample_fractions %>%
    dplyr::filter(predicted.celltype_3 %in% c("Erythroid Cell", "Immature B-cell", "Progenitor Cell"), !is.na(source_group)),
  aes(x = source_group, y = fraction, fill = source_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.5, alpha = 0.75) +
  facet_wrap(~ predicted.celltype_3, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("PBL" = "grey70", "Marrow" = "#D55E00"), guide = "none") +
  labs(x = NULL, y = "Fraction of cells", title = "Fig 2E red asterisk: marrow capture/source effect") +
  theme_prism(base_size = 14)

ggsave(plot_file("Fig2E_red_asterisk_source_effect_beta_binomial_boxplots.pdf"), p_red_source, width = 11, height = 4, units = "in", device = cairo_pdf)
ggsave(plot_file("Fig2E_red_asterisk_source_effect_beta_binomial_boxplots.eps"), p_red_source, width = 11, height = 4, units = "in", device = cairo_ps)
ggsave(plot_file("Fig2E_red_asterisk_source_effect_beta_binomial_boxplots.png"), p_red_source, width = 11, height = 4, units = "in", dpi = 400)

p_purple_car <- ggplot(
  propeller_sample_fractions %>%
    dplyr::filter(predicted.celltype_3 %in% c("Mature B-cell", "Immature B-cell", "Naive T-cell", "Effector T-cell"), car_control_group %in% c("Control", "CAR")),
  aes(x = car_control_group, y = fraction, fill = car_control_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.5, alpha = 0.75) +
  facet_wrap(~ predicted.celltype_3, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Control" = "grey70", "CAR" = "#CC79A7"), guide = "none") +
  labs(x = NULL, y = "Fraction of cells", title = "Fig 2E purple asterisk: CAR versus untreated control effect") +
  theme_prism(base_size = 14)

ggsave(plot_file("Fig2E_purple_asterisk_CAR_vs_control_beta_binomial_boxplots.pdf"), p_purple_car, width = 13, height = 4, units = "in", device = cairo_pdf)
ggsave(plot_file("Fig2E_purple_asterisk_CAR_vs_control_beta_binomial_boxplots.eps"), p_purple_car, width = 13, height = 4, units = "in", device = cairo_ps)
ggsave(plot_file("Fig2E_purple_asterisk_CAR_vs_control_beta_binomial_boxplots.png"), p_purple_car, width = 13, height = 4, units = "in", dpi = 400)

p_purple_total_b <- ggplot(
  propeller_sample_fractions_lvl2 %>%
    dplyr::filter(predicted.celltype_2 == "B-cell", car_control_group %in% c("Control", "CAR")),
  aes(x = car_control_group, y = fraction, fill = car_control_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.8, alpha = 0.75) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Control" = "grey70", "CAR" = "#CC79A7"), guide = "none") +
  labs(x = NULL, y = "B-cell fraction", title = "Fig 2E purple asterisk: total B-cell loss in CAR patients") +
  theme_prism(base_size = 14)

ggsave(plot_file("Fig2E_purple_asterisk_total_Bcell_CAR_vs_control_beta_binomial_boxplot.pdf"), p_purple_total_b, width = 4.5, height = 4, units = "in", device = cairo_pdf)
ggsave(plot_file("Fig2E_purple_asterisk_total_Bcell_CAR_vs_control_beta_binomial_boxplot.eps"), p_purple_total_b, width = 4.5, height = 4, units = "in", device = cairo_ps)
ggsave(plot_file("Fig2E_purple_asterisk_total_Bcell_CAR_vs_control_beta_binomial_boxplot.png"), p_purple_total_b, width = 4.5, height = 4, units = "in", dpi = 400)

readr::write_csv(manifest, figshare_file("figshare_scRNA_Fig2_propeller_dataset_manifest.csv"), na = "")
writeLines(capture.output(sessionInfo()), con = out_file("scRNA_Fig2_propeller_composition_sessionInfo.txt"))

cat("\nSaved Fig 2E propeller Figshare datasets to:\n", figshare_dir, "\n", sep = "")
cat("Saved Fig 2E propeller plots to:\n", plot_dir, "\n", sep = "")
