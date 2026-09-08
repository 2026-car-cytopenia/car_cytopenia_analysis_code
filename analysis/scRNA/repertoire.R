#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("repertoire")

#!/usr/bin/env Rscript

# =========================================================
# Figshare-ready datasets and figure-only code: TCR / BCR
# Project: CAR cytopenia paper
# Purpose:
#   1) Recreate only the TCR diversity and BCR/Ig-count figure panels.
#   2) Export cleaned, figure-level datasets for Figshare upload.
#   3) Apply the same exclusions/collapsing rules used in the plotted data.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(scRepertoire)
  library(scales)
  library(purrr)
  library(tibble)
})

# ----------------------------
# Paths
# ----------------------------
setwd(study_path("study_path_001"))

out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

figshare_dir <- out_file("figshare_datasets")
dir.create(figshare_dir, showWarnings = FALSE, recursive = TRUE)
figshare_file <- function(...) file.path(figshare_dir, ...)

# ----------------------------
# Input helpers
# ----------------------------
read_table_auto <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE, guess_max = 100000, na = c("", "NA", "NaN"))
  } else {
    readr::read_tsv(path, show_col_types = FALSE, guess_max = 100000, na = c("", "NA", "NaN"))
  }
}

first_existing <- function(candidates, pool) {
  hit <- intersect(candidates, pool)
  if (length(hit) > 0) return(hit[1])

  hit_made <- intersect(unique(make.names(candidates)), pool)
  if (length(hit_made) > 0) return(hit_made[1])

  normalize_name <- function(x) {
    stringr::str_to_lower(stringr::str_replace_all(as.character(x), "[^A-Za-z0-9]+", ""))
  }
  cand_norm <- normalize_name(candidates)
  pool_norm <- normalize_name(pool)
  idx <- match(cand_norm, pool_norm, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) > 0) return(pool[idx[1]])

  NA_character_
}

`%notin%` <- Negate(`%in%`)

normalize_token <- function(x) {
  x <- as.character(x)

  tok <- dplyr::case_when(
    stringr::str_detect(x, "^Auto\\d+") ~ stringr::str_extract(x, "^Auto\\d+"),
    stringr::str_detect(x, "^Ctrl\\d+") ~ stringr::str_extract(x, "^Ctrl\\d+"),
    stringr::str_detect(x, "^S\\d{4}")  ~ stringr::str_extract(x, "^S\\d{4}(?:[-o]\\d+|\\d)?"),
    TRUE ~ x
  )

  ifelse(
    stringr::str_detect(tok, "^S\\d{4}"),
    {
      base <- substr(tok, 1, 5)
      suf  <- substr(tok, 6, nchar(tok))
      suf  <- stringr::str_replace(suf, "^o", "-")
      suf  <- ifelse(stringr::str_detect(suf, "^\\d+$"), paste0("-", suf), suf)
      ifelse(nchar(suf) == 0, base, paste0(base, suf))
    },
    tok
  )
}

bucket_tissue <- function(tissue) {
  t2 <- tolower(as.character(tissue))
  dplyr::case_when(
    is.na(t2) ~ NA_character_,
    stringr::str_detect(t2, "pbmc|pbl|blood") ~ "PBMC",
    stringr::str_detect(t2, "marrow|bone")    ~ "Marrow",
    TRUE                                      ~ "Other"
  )
}

canonical_bcr_type <- function(type) {
  type <- as.character(type)
  dplyr::case_when(
    type %in% c("Auto PBMC", "HCT PBMC", "HCT PBL") ~ "HCT PBMC",
    type %in% c("Control PBMC", "Control PBL")       ~ "Control PBMC",
    type %in% c("Control Marrow")                     ~ "Control Marrow",
    type %in% c("CAR PBMC", "CAR PBL")               ~ "CAR PBMC",
    type %in% c("CAR Marrow")                         ~ "CAR Marrow",
    TRUE                                               ~ type
  )
}

# TCR marrow script assumptions:
# - Auto* are PBMC/HCT samples.
# - All other TCR samples are treated as marrow for this analysis.
infer_tcr_tissue <- function(norm_sample) {
  dplyr::case_when(
    stringr::str_detect(norm_sample, "^Auto\\d+") ~ "PBMC",
    TRUE                                          ~ "Marrow"
  )
}

run_pairwise_wilcox <- function(dat, outcome, group_var, comparisons) {
  purrr::map_dfr(comparisons, function(comp) {
    d <- dat %>%
      dplyr::filter(.data[[group_var]] %in% comp) %>%
      dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[group_var]]))

    g1 <- comp[1]
    g2 <- comp[2]
    x1 <- d %>% dplyr::filter(.data[[group_var]] == g1) %>% dplyr::pull(.data[[outcome]])
    x2 <- d %>% dplyr::filter(.data[[group_var]] == g2) %>% dplyr::pull(.data[[outcome]])

    wt <- if (length(x1) > 0 && length(x2) > 0) {
      tryCatch(wilcox.test(x1, x2, exact = FALSE), error = function(e) NULL)
    } else {
      NULL
    }

    tibble::tibble(
      outcome       = outcome,
      group_var     = group_var,
      group1        = g1,
      group2        = g2,
      n_group1      = length(x1),
      n_group2      = length(x2),
      median_group1 = ifelse(length(x1) > 0, median(x1, na.rm = TRUE), NA_real_),
      median_group2 = ifelse(length(x2) > 0, median(x2, na.rm = TRUE), NA_real_),
      statistic     = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
      p_value       = if (!is.null(wt)) wt$p.value else NA_real_
    )
  })
}

# ----------------------------
# Inputs
# ----------------------------
combined <- readRDS(study_path("input_all_tcr_combined_contigs"))

# Supply the BCR/Ig count table in private settings.
vdj_b_input_candidates <- study_path("input_bcr_counts")
vdj_b_input_file <- vdj_b_input_candidates[file.exists(vdj_b_input_candidates)][1]
if (is.na(vdj_b_input_file)) {
  stop("Could not find VDJ-B input file. Expected one of: ", paste(vdj_b_input_candidates, collapse = ", "))
}

message("[INFO] Reading VDJ-B/BCR input file: ", vdj_b_input_file)
vdj_b_raw <- read_table_auto(vdj_b_input_file)

# ----------------------------
# Internal TCR-to-paper timepoint map
# ----------------------------
# Used only to map scRepertoire sample names onto paper IDs/timepoints.
# These internal tokens are NOT written to Figshare datasets.
tcr_bcr_internal_timepoint_map <- study_setting("tcr_bcr_internal_timepoint_map")

# ----------------------------
# Robust column detection for either:
#   A) deidentified source: figshare_id, paper_id, day_post_cct, cohort_type, tissue, VDJ_B_predicted_cells
#   B) raw mapping source: Samples, Sample_ID, SID, MRN, Paper_ID, Day_Pot_CCT, Type, Tissue, VDJ_B_Predicted_Cells
# ----------------------------
vdj_sample_col <- first_existing(
  c("Sample_ID", "Sample.ID", "SampleID", "sample_id", "sample.id", "sample", "Sample"),
  names(vdj_b_raw)
)
vdj_paper_col <- first_existing(
  c("paper_id", "Paper_ID", "Paper.ID", "patient_id", "Patient_ID", "Patient.ID"),
  names(vdj_b_raw)
)
vdj_day_col <- first_existing(
  c("day_post_cct", "Day_Pot_CCT", "Day.Post.CCT", "Day_Post_CCT", "Day.Pot.CCT", "Day", "sample_day", "sample.day"),
  names(vdj_b_raw)
)
vdj_tissue_col <- first_existing(
  c("tissue", "Tissue", "Sample_Tissue", "Sample.Tissue", "sample_tissue", "sample.tissue"),
  names(vdj_b_raw)
)
vdj_type_col <- first_existing(
  c("cohort_type", "Type", "type", "Cohort_Type", "cohort.type"),
  names(vdj_b_raw)
)
vdj_count_col <- first_existing(
  c("VDJ_B_predicted_cells", "VDJ_B_Predicted_Cells", "VDJ_B_predicted", "VDJ.B.predicted", "VDJ_B_predicted_cells"),
  names(vdj_b_raw)
)

if (is.na(vdj_paper_col))  stop("VDJ-B/BCR file is missing paper_id/Paper_ID.")
if (is.na(vdj_tissue_col)) stop("VDJ-B/BCR file is missing tissue/Tissue.")
if (is.na(vdj_type_col))   stop("VDJ-B/BCR file is missing cohort_type/Type.")
if (is.na(vdj_count_col))  stop("VDJ-B/BCR file is missing VDJ_B_predicted_cells/VDJ_B_Predicted_Cells.")
if (is.na(vdj_day_col)) {
  message("[INFO] VDJ-B/BCR file has no day column; day_post_cct will be NA for BCR output.")
}
if (is.na(vdj_sample_col)) {
  message("[INFO] VDJ-B/BCR file has no Sample_ID column. This is OK for the deidentified Figshare source. TCR mapping will use the internal map embedded in this script.")
}

# =========================================================
# PART 0. Deidentified copy of the supplied VDJ-B/BCR file
# =========================================================

vdj_b_deidentified <- vdj_b_raw %>%
  dplyr::transmute(
    figshare_id = if ("figshare_id" %in% names(.)) as.character(.data[["figshare_id"]]) else sprintf("BCR_SOURCE_%03d", dplyr::row_number()),
    paper_id    = as.character(.data[[vdj_paper_col]]),
    day_post_cct = if (!is.na(vdj_day_col)) suppressWarnings(as.numeric(.data[[vdj_day_col]])) else NA_real_,
    cohort_type = canonical_bcr_type(.data[[vdj_type_col]]),
    tissue      = bucket_tissue(.data[[vdj_tissue_col]]),
    VDJ_B_predicted_cells = suppressWarnings(as.numeric(.data[[vdj_count_col]]))
  )

# Add the deidentified source file to the output directory and Figshare dataset directory.
readr::write_csv(vdj_b_deidentified, out_file("BCR_productive_Ig_counts_deidentified_for_figshare.csv"))
readr::write_csv(vdj_b_deidentified, figshare_file("figshare_BCR_productive_Ig_counts_deidentified_source.csv"))

# =========================================================
# PART 1. TCR diversity figure dataset
# =========================================================

# Cohort definitions used internally for mapping. These IDs are not exported.
remove  <- study_setting("remove")
auto    <- c("Auto1", "Auto2", "Auto3", "Auto4")
control <- c("Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4", "Ctrl5")

fill_cols_tcr <- c(
  "HCT"     = "#66C2A5",
  "Control" = "#FC8D62",
  "CAR"     = "#E78AC3"
)

clonaldiversity <- scRepertoire::clonalDiversity(
  combined,
  cloneCall   = "gene",
  group       = "sample",
  exportTable = TRUE
)

if (!"sample" %in% names(clonaldiversity)) {
  stop("Expected a sample column in clonalDiversity export.")
}

vdj_tp_map <- tcr_bcr_internal_timepoint_map %>%
  dplyr::mutate(
    tp_id = paste(
      paper_id,
      ifelse(is.na(day_post_cct), "day_not_available", as.character(day_post_cct)),
      tissue_bucket,
      sep = "|"
    )
  ) %>%
  dplyr::distinct(norm_sample, tissue_bucket, tp_id, paper_id, day_post_cct)

ambig <- vdj_tp_map %>%
  dplyr::count(norm_sample, tissue_bucket, name = "n_tp") %>%
  dplyr::filter(n_tp > 1)
if (nrow(ambig) > 0) {
  cat("\n[ERROR] Ambiguous TCR timepoint mapping for internal tokens. These are not exported.\n")
  print(ambig)
  stop("Fix internal TCR/BCR timepoint map so each internal token+tissue maps to one timepoint.")
}

tcr_tagged <- clonaldiversity %>%
  dplyr::filter(sample %notin% remove) %>%
  dplyr::mutate(
    norm_sample = normalize_token(sample),
    dataset = dplyr::case_when(
      norm_sample %in% auto    ~ "HCT",
      norm_sample %in% control ~ "Control",
      TRUE                     ~ "CAR"
    ),
    tissue_bucket = infer_tcr_tissue(norm_sample)
  ) %>%
  dplyr::left_join(vdj_tp_map, by = c("norm_sample", "tissue_bucket")) %>%
  dplyr::mutate(mapped_to_vdj = !is.na(tp_id))

# Count-only audit file; no internal sample IDs.
tcr_exclusion_summary <- tcr_tagged %>%
  dplyr::filter(dataset == "CAR", !mapped_to_vdj) %>%
  dplyr::summarise(
    excluded_unmapped_CAR_TCR_rows = dplyr::n(),
    .groups = "drop"
  )
readr::write_csv(tcr_exclusion_summary, out_file("audit_TCR_exclusion_summary_no_internal_ids.csv"))

tcr_tagged <- tcr_tagged %>%
  dplyr::filter(!(dataset == "CAR" & !mapped_to_vdj))

metrics_all <- intersect(c("Shannon", "Inv.Simpson", "Chao", "ACE", "Inv.Pielou"), names(tcr_tagged))
metrics_plot <- intersect(c("Inv.Simpson", "Chao", "Shannon"), metrics_all)
if (!all(c("Inv.Simpson", "Chao", "Shannon") %in% metrics_plot)) {
  stop("One or more plotted TCR diversity metrics are missing: Inv.Simpson, Chao, Shannon.")
}

tcr_diversity_figshare_wide <- tcr_tagged %>%
  dplyr::group_by(dataset, tp_id, paper_id, day_post_cct, tissue_bucket) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(metrics_all), ~mean(.x, na.rm = TRUE)),
    source_sample_count = dplyr::n_distinct(norm_sample),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dataset = factor(dataset, levels = c("HCT", "Control", "CAR"))) %>%
  dplyr::arrange(dataset, paper_id, day_post_cct, tissue_bucket) %>%
  dplyr::mutate(figshare_id = sprintf("TCR_%03d", dplyr::row_number())) %>%
  dplyr::select(
    figshare_id,
    cohort = dataset,
    paper_id,
    day_post_cct,
    tissue = tissue_bucket,
    source_sample_count,
    dplyr::all_of(metrics_all)
  )

tcr_diversity_figshare_long <- tcr_diversity_figshare_wide %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(metrics_plot),
    names_to  = "metric",
    values_to = "estimate"
  ) %>%
  dplyr::mutate(metric = factor(metric, levels = c("Inv.Simpson", "Chao", "Shannon"))) %>%
  dplyr::arrange(metric, cohort, paper_id, day_post_cct)

readr::write_csv(tcr_diversity_figshare_wide, figshare_file("figshare_TCR_diversity_plot_data_wide.csv"))
readr::write_csv(tcr_diversity_figshare_long, figshare_file("figshare_TCR_diversity_plot_data_long.csv"))

tcr_comparisons <- list(c("Control", "CAR"), c("HCT", "Control"), c("HCT", "CAR"))
tcr_stats <- purrr::map_dfr(metrics_plot, function(metric_i) {
  run_pairwise_wilcox(
    dat         = tcr_diversity_figshare_wide,
    outcome     = metric_i,
    group_var   = "cohort",
    comparisons = tcr_comparisons
  )
})
readr::write_csv(tcr_stats, figshare_file("figshare_TCR_diversity_plot_statistics.csv"))

tcr_n_labs <- tcr_diversity_figshare_wide %>%
  dplyr::count(cohort) %>%
  dplyr::mutate(lab = paste0(as.character(cohort), "\n", "n = ", n))
tcr_x_labs <- stats::setNames(tcr_n_labs$lab, as.character(tcr_n_labs$cohort))

make_tcr_div_plot <- function(yvar, title_txt) {
  ggplot(tcr_diversity_figshare_wide, aes(x = cohort, y = .data[[yvar]], fill = cohort)) +
    geom_boxplot(outlier.shape = NA, width = 0.65) +
    geom_point(color = "black", position = position_jitter(width = 0.12, height = 0), alpha = 0.75, size = 2) +
    stat_compare_means(comparisons = tcr_comparisons, method = "wilcox.test", label = "p.format", hide.ns = FALSE) +
    scale_fill_manual(values = fill_cols_tcr, drop = FALSE) +
    theme_classic(base_size = 24, base_family = "Arial") +
    labs(x = "", y = "Estimate", title = title_txt) +
    theme(legend.position = "none") +
    scale_x_discrete(name = "Cohort", labels = tcr_x_labs, drop = FALSE)
}

p_tcr_inv  <- make_tcr_div_plot("Inv.Simpson", "Inverted Simpson")
p_tcr_chao <- make_tcr_div_plot("Chao",        "Chao")
p_tcr_shan <- make_tcr_div_plot("Shannon",     "Shannon")

fig_tcr <- ggpubr::ggarrange(p_tcr_inv, p_tcr_chao, p_tcr_shan, ncol = 3, nrow = 1, align = "hv")

ggsave(out_file("marrow_TCR_diversity_plot.pdf"), plot = fig_tcr, width = 14, height = 6, device = cairo_pdf)
ggsave(out_file("marrow_TCR_diversity_plot.eps"), plot = fig_tcr, width = 14, height = 6, device = cairo_ps)
ggsave(out_file("marrow_TCR_diversity_plot.png"), plot = fig_tcr, width = 14, height = 6, dpi = 300)

# =========================================================
# PART 2. BCR / productive Ig-chain count figure dataset
# =========================================================

bcr_type_levels <- c("Control PBMC", "Control Marrow", "HCT PBMC", "CAR PBMC", "CAR Marrow")
bcr_comparisons <- list(
  c("Control PBMC", "Control Marrow"),
  c("CAR PBMC",     "CAR Marrow"),
  c("CAR Marrow",   "Control Marrow"),
  c("HCT PBMC",     "CAR PBMC")       # AutoPBL vs CAR PBL comparison
)
fill_cols_bcr <- c(
  "Control PBMC"   = "#8DA0CB",
  "Control Marrow" = "#FC8D62",
  "HCT PBMC"       = "#66C2A5",
  "CAR PBMC"       = "#A6D854",
  "CAR Marrow"     = "#E78AC3"
)

bcr_figshare <- vdj_b_deidentified %>%
  dplyr::mutate(
    cohort_type = factor(cohort_type, levels = bcr_type_levels),
    log10_VDJ_B_predicted_plus_0p01 = log10(VDJ_B_predicted_cells + 0.01)
  ) %>%
  dplyr::filter(!is.na(cohort_type), !is.na(VDJ_B_predicted_cells)) %>%
  dplyr::arrange(cohort_type, paper_id, day_post_cct, tissue) %>%
  dplyr::mutate(figshare_id = sprintf("BCR_%03d", dplyr::row_number())) %>%
  dplyr::select(
    figshare_id,
    cohort_type,
    paper_id,
    day_post_cct,
    tissue,
    VDJ_B_predicted_cells,
    log10_VDJ_B_predicted_plus_0p01
  )

readr::write_csv(bcr_figshare, figshare_file("figshare_BCR_productive_Ig_count_plot_data.csv"))

bcr_stats <- run_pairwise_wilcox(
  dat         = bcr_figshare,
  outcome     = "log10_VDJ_B_predicted_plus_0p01",
  group_var   = "cohort_type",
  comparisons = bcr_comparisons
)
readr::write_csv(bcr_stats, figshare_file("figshare_BCR_productive_Ig_count_plot_statistics.csv"))

bcr_n_labs <- bcr_figshare %>%
  dplyr::count(cohort_type) %>%
  dplyr::mutate(lab = paste0(as.character(cohort_type), "\n", "n = ", n))
bcr_x_labs <- stats::setNames(bcr_n_labs$lab, as.character(bcr_n_labs$cohort_type))

fig_bcr <- ggplot(bcr_figshare, aes(y = log10_VDJ_B_predicted_plus_0p01, x = cohort_type, fill = cohort_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_point(color = "black", position = position_jitter(width = 0.12, height = 0), alpha = 0.75, size = 2) +
  geom_hline(yintercept = -2) +
  stat_compare_means(comparisons = bcr_comparisons, method = "wilcox.test", label = "p.format", hide.ns = FALSE) +
  scale_fill_manual(values = fill_cols_bcr, drop = FALSE) +
  theme_classic(base_size = 24, base_family = "Arial") +
  labs(x = "", y = "Cell Count (log10)", title = "Cells with productive Ig chain") +
  theme(legend.position = "none") +
  scale_x_discrete(name = "Cohort", labels = bcr_x_labs, drop = FALSE)

ggsave(out_file("marrow_Ig_counts_plot.pdf"), plot = fig_bcr, width = 10, height = 6, device = cairo_pdf)
ggsave(out_file("marrow_Ig_counts_plot.eps"), plot = fig_bcr, width = 10, height = 6, device = cairo_ps)
ggsave(out_file("marrow_Ig_counts_plot.png"), plot = fig_bcr, width = 10, height = 6, dpi = 300)

# =========================================================
# PART 3. Figshare README / manifest
# =========================================================

manifest <- tibble::tribble(
  ~file, ~description,
  "figshare_TCR_diversity_plot_data_wide.csv", "Final wide dataset underlying the TCR diversity plot after exclusions, VDJ-B timepoint matching, and duplicate timepoint collapsing. Contains Paper_ID only; raw CCT/S/sample identifiers are removed.",
  "figshare_TCR_diversity_plot_data_long.csv", "Long-format version of the TCR diversity plot dataset with one row per plotted metric per sample/timepoint. Contains Paper_ID only; raw CCT/S/sample identifiers are removed.",
  "figshare_TCR_diversity_plot_statistics.csv", "Wilcoxon rank-sum comparisons used for the TCR diversity plot.",
  "figshare_BCR_productive_Ig_counts_deidentified_source.csv", "Deidentified BCR/VDJ-B source file used for Figshare. Raw Samples, Sample_ID, SID, MRN, CCT IDs, and S IDs are removed.",
  "figshare_BCR_productive_Ig_count_plot_data.csv", "Final dataset underlying the productive Ig-chain BCR/VDJ-B count plot after non-plotted Type levels are removed. Contains Paper_ID only; raw CCT/S/sample identifiers are removed.",
  "figshare_BCR_productive_Ig_count_plot_statistics.csv", "Wilcoxon rank-sum comparisons used for the BCR/VDJ-B productive Ig-chain plot, including HCT PBMC versus CAR PBMC."
)
readr::write_csv(manifest, figshare_file("figshare_TCR_BCR_dataset_manifest.csv"))

readr::write_lines(
  c(
    "Figshare datasets for TCR/BCR figure panels",
    "============================================================",
    "",
    "Generated by: TCR_BCR_figshare_datasets_figure_only_v6_autoPBL_vs_CAR_PBL.R",
    "",
    "TCR diversity dataset:",
    "- Source object: 20240929_All_TCR_combined_contigs.rds.",
    "- Diversity metrics computed with scRepertoire::clonalDiversity(cloneCall = 'gene', group = 'sample').",
    study_setting("study_value_001"),
    "- CAR samples without a matching final-paper VDJ-B timepoint map are excluded.",
    "- Duplicate TCR rows mapping to the same paper timepoint are collapsed by averaging diversity metrics.",
    "- Final cohorts: HCT, Control, CAR.",
    "",
    "BCR/Ig-count dataset:",
    paste0("- Source file used by script: ", vdj_b_input_file, "."),
    "- The preferred source is BCR_productive_Ig_counts_deidentified_for_figshare.csv.",
    "- Rows are retained only if Type/cohort_type maps to one of the plotted figure levels:",
    paste0("  ", paste(bcr_type_levels, collapse = ", ")),
    "- Plotted value is log10(VDJ_B_predicted_cells + 0.01).
- Pairwise comparisons include HCT PBMC versus CAR PBMC (AutoPBL vs CAR PBL).",
    "",
    "Deidentification / minimization:",
    "- Final Figshare datasets retain Paper_ID, cohort/type, tissue, day_post_cct, and plotted numeric values.",
    "- Raw CCT IDs, S IDs, Sample_ID, SID, MRN, and local sample labels are not exported in Figshare upload datasets.",
    "- A deidentified copy of the supplied BCR/VDJ-B file is written to both script_output and figshare_datasets.",
    "",
    "Output directory:",
    figshare_dir
  ),
  figshare_file("figshare_TCR_BCR_README.txt")
)

writeLines(capture.output(sessionInfo()), con = out_file("TCR_BCR_figshare_sessionInfo.txt"))

cat("\n[DONE] TCR/BCR figure plots and Figshare datasets written to:\n")
cat("  Figures:  ", out_dir, "\n", sep = "")
cat("  Datasets: ", figshare_dir, "\n", sep = "")

print(fig_tcr)
print(fig_bcr)
