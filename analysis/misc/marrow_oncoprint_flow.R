#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("marrow_oncoprint_flow")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(reshape2)
  library(survival)
  library(survminer)
  library(tableone)
  library(ggprism)
  library(broom)
  library(grid)
  library(Cairo)
})

`%notin%` <- Negate(`%in%`)

# =========================================================
# CONFIG
# =========================================================
base_oncoprint_dir <- study_path("base_oncoprint_dir")
base_flow_dir      <- study_path("base_flow_dir")

out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

figshare_dir <- file.path(out_dir, "figshare_datasets")
dir.create(figshare_dir, showWarnings = FALSE, recursive = TRUE)
figshare_file <- function(...) file.path(figshare_dir, ...)

MIN_SAMPLES_PER_GENE <- 2
font_family <- "Arial"
theme_set(theme_classic(base_family = font_family, base_size = 18))

pal_variant <- c(
  "DNMT3A" = "#0072B2",
  "TET2"   = "#009E73",
  "TP53"   = "#C43A31",
  "PPM1D"  = "#CC79A7",
  "ASXL1"  = "#E69F00",
  "ATM"    = "#56B4E9",
  "EZH2"   = "#000000",
  "IDH1"   = "#F0E442"
)

manifest <- tibble::tibble(
  file = character(),
  rows = integer(),
  columns = integer(),
  description = character()
)

# =========================================================
# HELPERS
# =========================================================
safe_read_delim <- function(path, delim = NULL) {
  if (!file.exists(path)) stop("Missing required input file: ", path)
  if (is.null(delim)) {
    readr::read_delim(path, show_col_types = FALSE, guess_max = 100000, progress = FALSE)
  } else {
    readr::read_delim(path, delim = delim, show_col_types = FALSE, guess_max = 100000, progress = FALSE)
  }
}

safe_read_tsv <- function(path) safe_read_delim(path, delim = "\t")

first_existing <- function(candidates, nm) {
  hit <- intersect(candidates, nm)
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

extract_lymid <- function(x) {
  stringr::str_extract(as.character(x), "LYM[0-9]+")
}

make_public_subject_id <- function(raw_id, prefix = "Subject") {
  raw_chr <- as.character(raw_id)
  lym <- extract_lymid(raw_chr)
  fallback <- paste0(prefix, "_", sprintf("%03d", as.integer(factor(raw_chr))))
  ifelse(!is.na(lym) & lym != "", lym, fallback)
}

remove_internal_id_columns <- function(df) {
  bad_exact <- c(
    "MRN", "mrn", "SID", "sid",
    "Sample_ID", "sample_id", "SAMPLE_ID", "SampleID", "sampleID",
    "CCT_ID", "CCTN", "Trial_ID", "trial_id", "UNIQUE.ID",
    "Specimen_ID", "Specimen", "raw_sample", "raw_subject_id", "raw_pair_id"
  )
  df %>% dplyr::select(-dplyr::any_of(bad_exact))
}

audit_no_disallowed_ids <- function(df, label) {
  if (nrow(df) == 0 || ncol(df) == 0) return(invisible(TRUE))

  chr_df <- df %>% dplyr::select(where(~is.character(.x) || is.factor(.x)))
  if (ncol(chr_df) == 0) return(invisible(TRUE))

  vals <- unlist(lapply(chr_df, as.character), use.names = FALSE)
  vals <- vals[!is.na(vals)]

  # LYM/AAID identifiers are intentionally allowed.
  bad_pattern <- study_setting("study_value_001")
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
    tibble::tibble(
      file = filename,
      rows = nrow(df_out),
      columns = ncol(df_out),
      description = description
    )
  )

  invisible(df_out)
}

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

norm_col <- function(df, choices, to, required = TRUE) {
  hit <- intersect(choices, names(df))
  if (length(hit) == 0) {
    if (required) stop("Missing required column among: ", paste(choices, collapse = ", "))
    df[[to]] <- NA
    return(df)
  }
  if (hit[[1]] != to) df <- dplyr::rename(df, !!to := !!rlang::sym(hit[[1]]))
  df
}

find_sample_col <- function(df) {
  candidates <- c("Sample", "sample", "SAMPLE", "Sample_ID", "SampleID", "CCT_ID", "CCTN", "Patient", "Patient_ID", "ID", "Specimen", "Specimen_ID")
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0) return(hit[[1]])
  pat_hits <- grep("(sample|cct|patient|specimen)", names(df), ignore.case = TRUE, value = TRUE)
  if (length(pat_hits) > 0) return(pat_hits[[1]])
  if (ncol(df) >= 18) return(names(df)[18])
  NA_character_
}

recode_variant <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("stop[-_ ]gained|nonsense|stopgain", "stopgain", x)
  x <- gsub("splice[-_ ](donor|acceptor)|splice[_ ]site|splice", "splice", x)
  x <- gsub("^missense.*", "missense", x)
  x <- gsub("frame[-_ ]shift.*", "frameshift", x)
  x <- gsub("^del(?!etion).*", "deletion", x, perl = TRUE)
  x <- gsub("^dup(?!lication).*", "duplication", x, perl = TRUE)
  x
}

to01 <- function(x) {
  if (is.numeric(x)) {
    ux <- unique(stats::na.omit(x))
    if (length(ux) > 0 && all(ux %in% c(0, 1))) return(as.integer(x))
    return(as.integer(ifelse(is.na(x), NA, x)))
  }
  x <- toupper(stringr::str_trim(as.character(x)))
  dplyr::case_when(
    x %in% c("1", "Y", "YES", "TRUE", "T", "POS", "POSITIVE") ~ 1L,
    x %in% c("0", "N", "NO", "FALSE", "F", "NEG", "NEGATIVE") ~ 0L,
    TRUE ~ NA_integer_
  )
}

# =========================================================
# PART 1. INPUTS: BM ONCOPRINT / tMN CLINICAL
# =========================================================
setwd(base_oncoprint_dir)

SNV <- safe_read_tsv(study_path("input_bm_oncoprint_snvs_4_update_20241011_with_mm"))
MDS_AML <- safe_read_tsv(study_path("input_mds_aml_clinical"))
MDS_AML_auto_v_CAR <- safe_read_tsv(study_path("input_mds_aml_auto_vs_car"))
clin <- safe_read_delim(study_path("input_marrow_clinical_5"))

# Normalize SNV columns
sample_col <- find_sample_col(SNV)
if (is.na(sample_col)) stop("Could not find a Sample-like column in SNV.")
if (sample_col != "Sample") SNV <- SNV %>% dplyr::rename(Sample = !!rlang::sym(sample_col))
SNV <- norm_col(SNV, c("Gene", "Hugo_Symbol", "Symbol", "Gene_Symbol"), "Gene")
SNV <- norm_col(SNV, c("VariantType", "Variant_Type", "Consequence", "Variant_Classification", "VC", "var_type"), "VariantType")
SNV <- SNV %>%
  dplyr::mutate(
    raw_sample = as.character(Sample),
    LYMID = make_public_subject_id(raw_sample, prefix = "OncoprintSubject"),
    Gene = as.character(Gene),
    VariantType = recode_variant(VariantType)
  )

clin <- norm_col(clin, c("CCT_ID", "CCTN", "Sample", "Patient", "UNIQUE.ID"), "CCT_ID") %>%
  dplyr::mutate(
    raw_sample = as.character(CCT_ID),
    LYMID = make_public_subject_id(raw_sample, prefix = "OncoprintSubject")
  )

# =========================================================
# PART 2. AUTO vs CAR tMN OUTCOMES: SURVIVAL / INCIDENCE DATASETS
# =========================================================
MDS_AML_auto_v_CAR_f <- MDS_AML_auto_v_CAR %>% dplyr::filter(Double_Keep < 3)
MDS_AML_auto_v_CAR_CI <- MDS_AML_auto_v_CAR %>% dplyr::filter(Double_Keep < 3, `Prior CCT` < 1)
MDS_AML_auto_v_CAR_NoHCT <- MDS_AML_auto_v_CAR_f %>% dplyr::filter(`Prior CCT` == 0)

subject_col_auto <- first_existing(c("LYMID", "AAID", "AA_Lab", "CCTN", "CCT_ID", "Patient", "Patient_ID"), names(MDS_AML_auto_v_CAR_f))
if (is.na(subject_col_auto)) {
  MDS_AML_auto_v_CAR_f$raw_subject_id <- seq_len(nrow(MDS_AML_auto_v_CAR_f))
  MDS_AML_auto_v_CAR_CI$raw_subject_id <- seq_len(nrow(MDS_AML_auto_v_CAR_CI))
  MDS_AML_auto_v_CAR_NoHCT$raw_subject_id <- seq_len(nrow(MDS_AML_auto_v_CAR_NoHCT))
  subject_col_auto <- "raw_subject_id"
}

if (!"CT_to_tMN" %in% names(MDS_AML_auto_v_CAR_f)) MDS_AML_auto_v_CAR_f$CT_to_tMN <- NA_real_
if (!"Cencode" %in% names(MDS_AML_auto_v_CAR_f)) MDS_AML_auto_v_CAR_f$Cencode <- NA_integer_

survival_plot_data <- MDS_AML_auto_v_CAR_f %>%
  dplyr::mutate(
    figure_subject_id = make_public_subject_id(.data[[subject_col_auto]], prefix = "AutoCARSubject")
  ) %>%
  dplyr::transmute(
    figure_subject_id,
    Cohort,
    prior_cellular_therapy = `Prior CCT`,
    double_keep = Double_Keep,
    PFS_Days,
    PFS_Event,
    OS_Days,
    OS_Event,
    CT_to_tMN = as.numeric(CT_to_tMN),
    Cencode = as.integer(Cencode)
  )

write_figshare_csv(
  survival_plot_data,
  "figshare_Auto_CAR_tMN_survival_plot_data.csv",
  "Deidentified data used for Auto vs CAR tMN PFS/OS and time-to-tMN analyses."
)

p_pfs <- ggsurvplot(
  survfit(Surv(PFS_Days, PFS_Event) ~ Cohort, data = MDS_AML_auto_v_CAR_NoHCT),
  xlab = "Months Post-Infusion", ylab = "Progression-Free Survival",
  conf.int = FALSE, pval = TRUE, xscale = "d_m",
  surv.median.line = "hv", palette = c("#0072B2", "#C43A31"),
  risk.table = TRUE, break.time.by = 365.25/4,
  risk.table.height = 0.32, xlim = c(0, 1000), risk.table.y.text = TRUE
)
ggsave(out_file("Auto_CAR_tMN_PFS.pdf"), plot = p_pfs$plot, width = 6, height = 5, device = cairo_pdf)
ggsave(out_file("Auto_CAR_tMN_PFS.eps"), plot = p_pfs$plot, width = 6, height = 5, device = cairo_ps)

p_os <- ggsurvplot(
  survfit(Surv(OS_Days, OS_Event) ~ Cohort, data = MDS_AML_auto_v_CAR_f),
  xlab = "Months Post-Infusion", ylab = "Overall Survival",
  conf.int = FALSE, pval = TRUE, xscale = "d_m",
  surv.median.line = "hv", palette = c("#0072B2", "#C43A31"),
  risk.table = TRUE, break.time.by = 365.25/4,
  risk.table.height = 0.32, xlim = c(0, 1000), risk.table.y.text = TRUE
)
ggsave(out_file("Auto_CAR_tMN_OS.pdf"), plot = p_os$plot, width = 6, height = 5, device = cairo_pdf)
ggsave(out_file("Auto_CAR_tMN_OS.eps"), plot = p_os$plot, width = 6, height = 5, device = cairo_ps)

# =========================================================
# PART 3. GENE ODDS RATIO FOR tMN / MYELOID MALIGNANCY
# =========================================================
needed_genes <- c("TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1")
for (g in needed_genes) {
  if (!(g %in% names(clin))) clin[[g]] <- 0L
  clin[[g]] <- as.integer(clin[[g]] %in% c(1, "1", "Yes", "YES", "TRUE", "True", TRUE))
}
if (!"MM" %in% names(clin)) stop("Column 'MM' not found in clin.")
clin$MM <- as.integer(clin$MM %in% c(1, "1", "Yes", "YES", "TRUE", "True", TRUE))

clin_or_input <- clin %>%
  dplyr::transmute(
    LYMID,
    myeloid_malignancy = MM,
    TP53, DNMT3A, PPM1D, TET2, ASXL1
  )

write_figshare_csv(
  clin_or_input,
  "figshare_gene_odds_ratio_tMN_model_input.csv",
  "Deidentified gene mutation indicators and tMN/myeloid malignancy outcome used for gene odds-ratio analyses."
)

fit_gene_or <- NULL
or_tbl <- tibble::tibble(
  term = needed_genes,
  OR = NA_real_,
  LCL = NA_real_,
  UCL = NA_real_,
  p = NA_real_,
  p_label = NA_character_,
  method = NA_character_
)

if (requireNamespace("logistf", quietly = TRUE)) {
  fit_gene_or <- logistf::logistf(MM ~ TP53 + DNMT3A + PPM1D + TET2 + ASXL1, data = clin)
  ci_f <- confint(fit_gene_or)
  co_f <- coef(fit_gene_or)
  summ <- summary(fit_gene_or)
  p_vec <- NULL
  if (!is.null(summ$coefficients) && any(colnames(summ$coefficients) %in% c("p", "P", "Prob", "prob"))) {
    p_col <- intersect(c("p", "P", "Prob", "prob"), colnames(summ$coefficients))[[1]]
    p_vec <- summ$coefficients[, p_col, drop = TRUE]
  } else if (!is.null(summ$prob)) {
    p_vec <- summ$prob
  }
  if (is.null(p_vec) && !is.null(summ$coefficients) && "Chisq" %in% colnames(summ$coefficients)) {
    p_vec <- pchisq(summ$coefficients[, "Chisq", drop = TRUE], df = 1, lower.tail = FALSE)
  }
  if (is.null(p_vec)) p_vec <- rep(NA_real_, length(co_f))
  names(p_vec) <- names(co_f)

  or_tbl <- tibble::tibble(
    term = names(co_f),
    OR = exp(co_f),
    LCL = exp(ci_f[, 1]),
    UCL = exp(ci_f[, 2]),
    p = as.numeric(p_vec[names(co_f)]),
    method = "Firth logistic regression"
  ) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      term = factor(term, levels = needed_genes),
      p_label = fmt_p(p)
    )
} else {
  warning("Package 'logistf' is not installed; using standard glm for gene OR export.")
  fit_gene_or <- glm(MM ~ TP53 + DNMT3A + PPM1D + TET2 + ASXL1, data = clin, family = binomial)
  or_tbl <- broom::tidy(fit_gene_or, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::transmute(term, OR = estimate, LCL = conf.low, UCL = conf.high, p = p.value, p_label = fmt_p(p), method = "glm") %>%
    dplyr::mutate(term = factor(term, levels = needed_genes))
}

write_figshare_csv(
  or_tbl %>% dplyr::mutate(term = as.character(term)),
  "figshare_gene_odds_ratio_tMN_statistics.csv",
  "Gene-specific odds ratios for tMN/myeloid malignancy."
)

x_min <- max(1e-3, min(or_tbl$LCL, na.rm = TRUE) / 5)
x_max <- min(1e3, max(or_tbl$UCL, na.rm = TRUE) * 5)
if (!is.finite(x_min)) x_min <- 0.01
if (!is.finite(x_max)) x_max <- 100
plot_or <- or_tbl %>%
  dplyr::mutate(
    term = factor(as.character(term), levels = needed_genes),
    LCL_clip = pmax(LCL, x_min),
    UCL_clip = pmin(UCL, x_max)
  )

gene_cols <- c(
  TP53 = pal_variant[["TP53"]],
  DNMT3A = pal_variant[["DNMT3A"]],
  PPM1D = pal_variant[["PPM1D"]],
  TET2 = pal_variant[["TET2"]],
  ASXL1 = pal_variant[["ASXL1"]]
)

p_or <- ggplot(plot_or, aes(x = OR, y = term, color = term)) +
  geom_errorbarh(aes(xmin = LCL_clip, xmax = UCL_clip), height = 0.15, linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, color = "grey35") +
  scale_x_log10(limits = c(x_min, x_max)) +
  scale_color_manual(values = gene_cols, guide = "none") +
  labs(title = "Odds Ratio for Myeloid Malignancy by Gene Mutation", x = "Odds Ratio for tMN", y = NULL) +
  theme_classic(base_family = font_family, base_size = 18)

ggsave(out_file("Gene_OR_tMN_firth_colored.pdf"), p_or, width = 7, height = 5, device = cairo_pdf)
ggsave(out_file("Gene_OR_tMN_firth_colored.png"), p_or, width = 7, height = 5, dpi = 400)

# =========================================================
# PART 4. ONCOPRINT MATRIX AND ANNOTATIONS
# =========================================================
op_data_full <- reshape2::acast(
  SNV,
  Gene ~ raw_sample,
  value.var = "VariantType",
  fun.aggregate = function(v) paste(sort(unique(v[!is.na(v) & v != ""])), collapse = ",")
)

op_data_recurrent <- op_data_full[
  rowSums(!is.na(op_data_full) & op_data_full != "") >= MIN_SAMPLES_PER_GENE,
  , drop = FALSE
]
op_data_recurrent <- op_data_recurrent[
  !(rownames(op_data_recurrent) %in% c("none", "NA", "Unknown", "")),
  , drop = FALSE
]

common_samples <- intersect(colnames(op_data_recurrent), clin$raw_sample)
if (length(common_samples) == 0) stop("No overlap between oncoprint columns and clinical metadata sample IDs.")
op_data_recurrent <- op_data_recurrent[, common_samples, drop = FALSE]
clin_genotyping <- clin %>% dplyr::filter(raw_sample %in% common_samples)
clin_genotyping <- clin_genotyping[match(common_samples, clin_genotyping$raw_sample), , drop = FALSE]

patient_order_tx <- study_setting("patient_order_tx")
raw_to_public <- tibble::tibble(raw_sample = common_samples, LYMID = make_public_subject_id(common_samples, prefix = "OncoprintSubject"))
ordered_raw <- raw_to_public %>%
  dplyr::mutate(order_rank = match(LYMID, patient_order_tx)) %>%
  dplyr::arrange(is.na(order_rank), order_rank, LYMID) %>%
  dplyr::pull(raw_sample)

op_data_recurrent <- op_data_recurrent[, ordered_raw, drop = FALSE]
clin_genotyping <- clin_genotyping[match(ordered_raw, clin_genotyping$raw_sample), , drop = FALSE]

oncoprint_long <- as.data.frame(as.table(op_data_recurrent), stringsAsFactors = FALSE) %>%
  dplyr::rename(Gene = Var1, raw_sample = Var2, variant_type = Freq) %>%
  dplyr::filter(!is.na(variant_type), variant_type != "") %>%
  tidyr::separate_rows(variant_type, sep = ",") %>%
  dplyr::mutate(
    variant_type = stringr::str_trim(variant_type),
    LYMID = make_public_subject_id(raw_sample, prefix = "OncoprintSubject"),
    oncoprint_column_order = match(raw_sample, ordered_raw),
    oncoprint_gene_order = match(Gene, rownames(op_data_recurrent))
  ) %>%
  dplyr::select(LYMID, oncoprint_column_order, Gene, oncoprint_gene_order, variant_type)

write_figshare_csv(
  oncoprint_long,
  "figshare_BM_oncoprint_mutation_long.csv",
  "Long-format mutation data used to generate the recurrent-gene bone marrow oncoprint."
)

annotation_candidates <- intersect(c("Diagnosis", "Histology", "MM", "TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1"), names(clin_genotyping))
oncoprint_annotations <- clin_genotyping %>%
  dplyr::mutate(
    LYMID = make_public_subject_id(raw_sample, prefix = "OncoprintSubject"),
    oncoprint_column_order = match(raw_sample, ordered_raw)
  ) %>%
  dplyr::select(LYMID, oncoprint_column_order, dplyr::all_of(annotation_candidates))

write_figshare_csv(
  oncoprint_annotations,
  "figshare_BM_oncoprint_annotations.csv",
  "Clinical annotations aligned to the recurrent-gene bone marrow oncoprint."
)

# Draw oncoprint using raw matrix columns internally; display LYMID/public labels.
anno_df <- data.frame(
  Diagnosis = if ("Diagnosis" %in% names(clin_genotyping)) clin_genotyping$Diagnosis else factor(rep(NA, nrow(clin_genotyping))),
  Histology = if ("Histology" %in% names(clin_genotyping)) clin_genotyping$Histology else factor(rep(NA, nrow(clin_genotyping))),
  row.names = clin_genotyping$raw_sample,
  stringsAsFactors = TRUE
)
anno_cols <- list(
  Diagnosis = c("Cytopenia" = "#0072B2", "MDS/AML" = "#C43A31", "AML" = "#C43A31", "MDS" = "#E69F00"),
  Histology = c("LBCL" = "steelblue3", "FL" = "orchid4", "MCL" = "palegreen4")
)
top_anno <- ComplexHeatmap::HeatmapAnnotation(
  df = anno_df[colnames(op_data_recurrent), , drop = FALSE],
  col = anno_cols,
  annotation_name_gp = grid::gpar(fontsize = 12, fontfamily = font_family, col = "black"),
  simple_anno_size = grid::unit(4, "mm")
)

col_map <- c(
  "splice" = "dodgerblue",
  "stopgain" = "cyan3",
  "missense" = "deepskyblue4",
  "deletion" = "#2c7fb8",
  "duplication" = "#253494",
  "frameshift" = "#003366"
)

alter_fun <- function(x, y, w, h, v) {
  n <- sum(v)
  grid::grid.rect(x, y, w, h * 0.9, gp = grid::gpar(fill = "#DDDDDD", col = NA))
  if (n == 0) return(invisible(NULL))
  classes <- names(which(v))
  seg_h <- (h * 0.9) / n
  y_top <- y + (h * 0.45)
  for (i in seq_len(n)) {
    grid::grid.rect(x, y_top - (i - 1) * seg_h, w * 0.9, seg_h,
                    just = "top", gp = grid::gpar(fill = col_map[classes[i]], col = NA))
  }
}

column_labels <- make_public_subject_id(colnames(op_data_recurrent), prefix = "OncoprintSubject")
op <- ComplexHeatmap::oncoPrint(
  op_data_recurrent,
  alter_fun = alter_fun,
  col = col_map,
  show_column_names = TRUE,
  column_order = colnames(op_data_recurrent),
  column_labels = column_labels,
  top_annotation = top_anno,
  row_names_gp = grid::gpar(fontsize = 10, fontfamily = font_family, col = "black"),
  column_names_gp = grid::gpar(fontsize = 9, fontfamily = font_family, col = "black"),
  pct_gp = grid::gpar(fontsize = 9, fontfamily = font_family, col = "black")
)
Cairo::CairoPDF(out_file("BM_oncoprint_recurrent.pdf"), width = 10, height = 8, family = font_family)
ComplexHeatmap::draw(op, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

# =========================================================
# PART 5. MARROW FLOW: CD19 / CAR19 PBMC vs BM
# =========================================================
setwd(base_flow_dir)
flow <- read.delim(study_path("input_flow"), stringsAsFactors = FALSE)
car_cd19 <- read.delim(study_path("input_car_cd19"), stringsAsFactors = FALSE)

flow_pair_col <- first_existing(c("LYMID", "AAID", "AA_Lab", "Sample_ID", "CCT_ID", "CCTN", "Patient_ID"), names(flow))
if (is.na(flow_pair_col)) {
  flow$raw_pair_id <- seq_len(nrow(flow))
  flow_pair_col <- "raw_pair_id"
}

car_cd19_pair_col <- first_existing(c("LYMID", "AAID", "AA_Lab", "Sample_ID", "CCT_ID", "CCTN", "Patient_ID"), names(car_cd19))
if (is.na(car_cd19_pair_col)) {
  car_cd19$raw_pair_id <- seq_len(nrow(car_cd19))
  car_cd19_pair_col <- "raw_pair_id"
}

car_cd19_plot_data <- car_cd19 %>%
  dplyr::mutate(
    figure_pair_id = make_public_subject_id(.data[[car_cd19_pair_col]], prefix = "FlowPair"),
    log10_CD19 = log10(CD19_Flow_BM + 0.01),
    log10_CAR19 = log10(Marrow_Sum + 0.01)
  ) %>%
  dplyr::transmute(
    figure_pair_id,
    CD19_Flow_BM,
    Marrow_Sum,
    log10_CD19,
    log10_CAR19
  )
write_figshare_csv(
  car_cd19_plot_data,
  "figshare_CD19_vs_CAR19_correlation_plot_data.csv",
  "Data used for the marrow CD19 versus CAR19 correlation plot."
)

p_scatter <- ggpubr::ggscatter(
  car_cd19 %>% dplyr::mutate(log10_CD19 = log10(CD19_Flow_BM + 0.01), log10_CAR19 = log10(Marrow_Sum + 0.01)),
  x = "log10_CAR19", y = "log10_CD19",
  color = "black", shape = 21, size = 3,
  add = "reg.line", add.params = list(color = "grey30", fill = "grey85"),
  conf.int = TRUE, cor.coef = TRUE,
  cor.coeff.args = list(method = "pearson", label.sep = "\n", size = 6)
) +
  xlab("Log10 CAR19+ Bone Marrow Cells") +
  ylab("Log10 CD19+ Bone Marrow Cells")
ggsave(out_file("CD19_vs_CAR19.pdf"), p_scatter, width = 9, height = 7, device = cairo_pdf)
ggsave(out_file("CD19_vs_CAR19.eps"), p_scatter, width = 9, height = 7, device = cairo_ps)

flow_public <- flow %>%
  dplyr::mutate(
    figure_pair_id = make_public_subject_id(.data[[flow_pair_col]], prefix = "FlowPair")
  )

flow_cyto_wide <- flow_public %>%
  dplyr::filter(Category %in% c("Long_Term_Cytopenia", "AML_Post_Treatment")) %>%
  dplyr::select(figure_pair_id, Category, CD19_Flow_PBMC, CD19_Flow_BM) %>%
  tidyr::drop_na(CD19_Flow_PBMC, CD19_Flow_BM)

flow_cyto2_wide <- flow_public %>%
  dplyr::filter(Category %in% c("Long_Term_Cytopenia", "AML_Post_Treatment")) %>%
  dplyr::select(figure_pair_id, Category, PB_Sum, Marrow_Sum) %>%
  tidyr::drop_na(PB_Sum, Marrow_Sum)

flow_ctrl_wide <- flow_public %>%
  dplyr::filter(Category == "Healthy Control") %>%
  dplyr::select(figure_pair_id, Category, CD19_Flow_PBMC, CD19_Flow_BM) %>%
  tidyr::drop_na(CD19_Flow_PBMC, CD19_Flow_BM)

flow_ctrl_car_wide <- flow_public %>%
  dplyr::filter(Category == "Healthy Control") %>%
  dplyr::select(figure_pair_id, Category, PB_Sum, Marrow_Sum) %>%
  tidyr::drop_na(PB_Sum, Marrow_Sum)

flow_cd19_cases_long <- flow_cyto_wide %>%
  tidyr::pivot_longer(cols = c("CD19_Flow_BM", "CD19_Flow_PBMC"), names_to = "source", values_to = "percent_CD19_positive") %>%
  dplyr::mutate(source = factor(source, levels = c("CD19_Flow_PBMC", "CD19_Flow_BM"), labels = c("PBMC", "BM")))

flow_car_cases_long <- flow_cyto2_wide %>%
  tidyr::pivot_longer(cols = c("PB_Sum", "Marrow_Sum"), names_to = "source", values_to = "percent_CAR19_positive") %>%
  dplyr::mutate(source = factor(source, levels = c("PB_Sum", "Marrow_Sum"), labels = c("PBMC", "BM")))

flow_cd19_controls_long <- flow_ctrl_wide %>%
  tidyr::pivot_longer(cols = c("CD19_Flow_PBMC", "CD19_Flow_BM"), names_to = "source", values_to = "percent_CD19_positive") %>%
  dplyr::mutate(source = factor(source, levels = c("CD19_Flow_PBMC", "CD19_Flow_BM"), labels = c("PBMC", "BM")))

flow_car_controls_long <- flow_ctrl_car_wide %>%
  tidyr::pivot_longer(cols = c("PB_Sum", "Marrow_Sum"), names_to = "source", values_to = "percent_CAR19_positive") %>%
  dplyr::mutate(source = factor(source, levels = c("PB_Sum", "Marrow_Sum"), labels = c("PBMC", "BM")))

write_figshare_csv(flow_cd19_cases_long, "figshare_CD19_PBMC_vs_BM_cases_plot_data.csv", "Paired CD19+ cell percentage data for cytopenia/tMN case PBMC versus bone marrow plot.")
write_figshare_csv(flow_car_cases_long, "figshare_CAR19_PBMC_vs_BM_cases_plot_data.csv", "Paired CAR19+ T-cell percentage data for cytopenia/tMN case PBMC versus bone marrow plot.")
write_figshare_csv(flow_cd19_controls_long, "figshare_CD19_PBMC_vs_BM_controls_plot_data.csv", "Paired CD19+ cell percentage data for healthy control PBMC versus bone marrow plot.")
write_figshare_csv(flow_car_controls_long, "figshare_CAR19_PBMC_vs_BM_controls_plot_data.csv", "Paired CAR19+ T-cell percentage data for healthy control PBMC versus bone marrow plot.")

paired_wilcox <- function(df, x, y, label) {
  wt <- wilcox.test(df[[x]], df[[y]], paired = TRUE, conf.int = TRUE, exact = FALSE)
  tibble::tibble(
    comparison = label,
    n_pairs = nrow(df),
    median_first = median(df[[x]], na.rm = TRUE),
    median_second = median(df[[y]], na.rm = TRUE),
    p_value = wt$p.value,
    p_label = fmt_p(wt$p.value)
  )
}

flow_stats <- dplyr::bind_rows(
  paired_wilcox(flow_cyto_wide, "CD19_Flow_PBMC", "CD19_Flow_BM", "Cases: CD19 PBMC vs BM"),
  paired_wilcox(flow_cyto2_wide, "PB_Sum", "Marrow_Sum", "Cases: CAR19 PBMC vs BM"),
  paired_wilcox(flow_ctrl_wide, "CD19_Flow_PBMC", "CD19_Flow_BM", "Healthy controls: CD19 PBMC vs BM"),
  paired_wilcox(flow_ctrl_car_wide, "PB_Sum", "Marrow_Sum", "Healthy controls: CAR19 PBMC vs BM")
)
write_figshare_csv(flow_stats, "figshare_BM_vs_PBMC_flow_paired_wilcoxon_statistics.csv", "Paired Wilcoxon statistics for CD19/CAR19 PBMC versus bone marrow plots.")

# Plot linear versions used for editing/figures
stat_case_cd19 <- tibble::tibble(group1 = "PBMC", group2 = "BM", p.adj = flow_stats$p_value[flow_stats$comparison == "Cases: CD19 PBMC vs BM"])
stat_case_car  <- tibble::tibble(group1 = "PBMC", group2 = "BM", p.adj = flow_stats$p_value[flow_stats$comparison == "Cases: CAR19 PBMC vs BM"])
stat_ctrl_cd19 <- tibble::tibble(group1 = "PBMC", group2 = "BM", p.adj = flow_stats$p_value[flow_stats$comparison == "Healthy controls: CD19 PBMC vs BM"])
stat_ctrl_car  <- tibble::tibble(group1 = "PBMC", group2 = "BM", p.adj = flow_stats$p_value[flow_stats$comparison == "Healthy controls: CAR19 PBMC vs BM"])

lab_pct <- function(x) paste0(x, "%")
y_breaks_car19 <- c(0, 5, 10, 15, 25)
y_lim_car19 <- c(0, 25)
y_breaks_cd19 <- c(0, 10, 20, 30, 40, 50)
y_lim_cd19 <- c(0, 50)

plot_pair_lines <- function(dat, y_col, stat_dat, y_position, y_breaks, y_lim, y_lab, out_base) {
  p <- ggplot(dat, aes(x = source, y = .data[[y_col]], group = figure_pair_id)) +
    geom_line(aes(color = figure_pair_id), linewidth = 1) +
    geom_point(aes(color = figure_pair_id), size = 2) +
    ggpubr::stat_pvalue_manual(stat_dat, label = "p.adj", y.position = y_position, size = 6) +
    scale_color_prism("colors") +
    scale_y_continuous(breaks = y_breaks, labels = lab_pct) +
    coord_cartesian(ylim = y_lim) +
    theme_classic(base_size = 20, base_family = font_family) +
    theme(legend.position = "none") +
    ylab(y_lab) + xlab("")
  ggsave(out_file(paste0(out_base, ".pdf")), p, width = 9, height = 7, device = cairo_pdf)
  ggsave(out_file(paste0(out_base, ".png")), p, width = 9, height = 7, dpi = 400)
  p
}

p_slope_cd19_linear <- plot_pair_lines(flow_cd19_cases_long, "percent_CD19_positive", stat_case_cd19, 50, y_breaks_cd19, y_lim_cd19, "Percent CD19+ Cells", "CD19_PBMC_vs_marrow_LINEAR")
p_slope_car_linear  <- plot_pair_lines(flow_car_cases_long, "percent_CAR19_positive", stat_case_car, 25, y_breaks_car19, y_lim_car19, "Percent CAR19+ T cells", "CAR19_PB_vs_BM_LINEAR")
p_ctrl_cd19_linear  <- plot_pair_lines(flow_cd19_controls_long, "percent_CD19_positive", stat_ctrl_cd19, 50, y_breaks_cd19, y_lim_cd19, "Percent CD19+ Cells", "CD19_controls_LINEAR")
p_ctrl_car_linear   <- plot_pair_lines(flow_car_controls_long, "percent_CAR19_positive", stat_ctrl_car, 25, y_breaks_car19, y_lim_car19, "Percent CAR19+ T cells", "CAR19_controls_LINEAR")

# BM enrichment summaries
summ_enrich <- function(df, metric, pb_col, bm_col, cohort) {
  dat <- df %>%
    dplyr::mutate(
      ratio_BM_over_PB = (.data[[bm_col]] + 0.01) / (.data[[pb_col]] + 0.01),
      log10_BM = log10(.data[[bm_col]] + 0.01),
      log10_PB = log10(.data[[pb_col]] + 0.01)
    )
  wt <- wilcox.test(dat$log10_BM, dat$log10_PB, paired = TRUE, conf.int = TRUE, conf.level = 0.95, exact = FALSE)
  HL_logFC <- unname(wt$estimate)
  tibble::tibble(
    cohort = cohort,
    metric = metric,
    n_pairs = nrow(dat),
    median_BM = median(dat[[bm_col]], na.rm = TRUE),
    median_PBMC = median(dat[[pb_col]], na.rm = TRUE),
    median_ratio_BM_over_PB = median(dat$ratio_BM_over_PB, na.rm = TRUE),
    geom_mean_ratio_BM_over_PB = 10^mean(dat$log10_BM - dat$log10_PB, na.rm = TRUE),
    prop_enriched_BM = mean(dat$ratio_BM_over_PB > 1, na.rm = TRUE),
    HL_median_fold_change = 10^HL_logFC,
    HL_FC_CI_low = 10^wt$conf.int[[1]],
    HL_FC_CI_high = 10^wt$conf.int[[2]],
    wilcox_p = wt$p.value,
    p_label = fmt_p(wt$p.value)
  )
}

bm_enrichment_summary <- dplyr::bind_rows(
  summ_enrich(flow_cyto_wide, "CD19", "CD19_Flow_PBMC", "CD19_Flow_BM", "Cytopenia/tMN cases"),
  summ_enrich(flow_cyto2_wide, "CAR19", "PB_Sum", "Marrow_Sum", "Cytopenia/tMN cases"),
  summ_enrich(flow_ctrl_wide, "CD19", "CD19_Flow_PBMC", "CD19_Flow_BM", "Healthy controls"),
  summ_enrich(flow_ctrl_car_wide, "CAR19", "PB_Sum", "Marrow_Sum", "Healthy controls")
) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4)))

write_figshare_csv(bm_enrichment_summary, "figshare_BM_vs_PBMC_enrichment_summary.csv", "Summary statistics for BM enrichment of CD19+ cells and CAR19+ T cells.")
readr::write_tsv(bm_enrichment_summary, out_file("BM_vs_PBMC_enrichment_summary.tsv"))

# =========================================================
# PART 6. tMN CUMULATIVE INCIDENCE + RISK FACTORS
# =========================================================
setwd(base_oncoprint_dir)
tmn_risk <- safe_read_delim(study_path("input_tmn_risk_variables"))

tmn_subject_col <- first_existing(c("LYMID", "AAID", "AA_Lab", "CCTN", "CCT_ID", "Patient_ID", "Patient"), names(tmn_risk))
if (is.na(tmn_subject_col)) {
  tmn_risk$raw_subject_id <- seq_len(nrow(tmn_risk))
  tmn_subject_col <- "raw_subject_id"
}

tmn_risk_f <- tmn_risk %>%
  dplyr::distinct(.data[[tmn_subject_col]], .keep_all = TRUE) %>%
  dplyr::mutate(
    figure_subject_id = make_public_subject_id(.data[[tmn_subject_col]], prefix = "tMNRiskSubject"),
    Time_to_Event_or_follow_up = as.numeric(Time_to_Event_or_follow_up),
    cen_code_death_2 = ifelse(Time_to_Event_or_follow_up >= 1500, 0, cen_code_death),
    Event_time_limited = pmin(Time_to_Event_or_follow_up, 1500)
  )

if (!"tMN" %in% names(tmn_risk_f)) tmn_risk_f$tMN <- NA_integer_

ci_dataset <- tmn_risk_f %>%
  dplyr::transmute(
    figure_subject_id,
    time_to_event_or_follow_up_days = Time_to_Event_or_follow_up,
    event_time_limited_days = Event_time_limited,
    competing_event_code = cen_code_death,
    competing_event_code_limited = cen_code_death_2,
    tMN = to01(tMN)
  )
write_figshare_csv(ci_dataset, "figshare_tMN_cumulative_incidence_plot_data.csv", "Data used for cumulative incidence of tMN or death plot.")

# cmprsk2 does not export cuminc(); the cumulative-incidence estimator is from cmprsk.
# Use cmprsk::cuminc when available and otherwise skip only the plot, while still exporting
# figshare_tMN_cumulative_incidence_plot_data.csv above.
if (requireNamespace("cmprsk", quietly = TRUE)) {
  ci_tmn_2 <- cmprsk::cuminc(
    ftime   = tmn_risk_f$Event_time_limited,
    fstatus = tmn_risk_f$cen_code_death_2,
    cencode = 0
  )

  p_ci_all <- survminer::ggcompetingrisks(
    ci_tmn_2,
    multiple_panels = FALSE,
    legend = "right",
    risk.table = TRUE,
    conf.int = TRUE
  ) +
    scale_color_brewer(palette = "Paired") +
    ylim(0, 1) +
    scale_x_continuous(
      breaks = c(0, 365, 730, 1095, 1460),
      labels = c("0", "12m", "24m", "36m", "48m")
    ) +
    labs(y = "Cumulative Incidence", x = "Time from Infusion (Days)") +
    ggtitle("Cumulative Incidence of tMN or Death")

  ggsave(out_file("tMN_competing_risks_all.pdf"), p_ci_all, width = 8, height = 6, device = cairo_pdf)
  ggsave(out_file("tMN_competing_risks_all.png"), p_ci_all, width = 8, height = 6, dpi = 400)
} else {
  message("Package 'cmprsk' is not installed; skipping tMN_competing_risks_all plot. The Figshare cumulative-incidence dataset was still exported.")
}

has_brglm2 <- requireNamespace("brglm2", quietly = TRUE)
fit_logit_safe <- function(formula, data) {
  m1 <- suppressWarnings(glm(formula, data = data, family = binomial(link = "logit")))
  co <- coef(m1)
  v <- suppressWarnings(try(vcov(m1), silent = TRUE))
  bad <- any(!is.finite(co)) || inherits(v, "try-error") || any(!is.finite(diag(v))) || any(abs(co[is.finite(co)]) > 12)
  if (!bad) return(list(model = m1, method = "glm"))
  if (has_brglm2) {
    m2 <- brglm2::brglm(formula, data = data, family = binomial("logit"), type = "AS_mean")
    return(list(model = m2, method = "brglm2"))
  }
  list(model = m1, method = "glm_separated")
}

tidy_or_ci <- function(fit_obj) {
  m <- fit_obj$model
  td <- broom::tidy(m) %>% dplyr::filter(term != "(Intercept)")
  V <- vcov(m)
  se <- sqrt(diag(V))
  td %>%
    dplyr::mutate(
      se = se[match(term, names(se))],
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * se),
      CI_high = exp(estimate + 1.96 * se),
      method = fit_obj$method
    ) %>%
    dplyr::select(term, estimate, se, OR, CI_low, CI_high, p.value, method)
}

pretty_term <- function(var, term) {
  if (term == var) return(var)
  if (startsWith(term, var)) {
    lvl <- sub(paste0("^", var), "", term)
    return(paste0(var, ": ", lvl))
  }
  paste0(var, ": ", term)
}

tmn_risk_model <- tmn_risk_f %>%
  dplyr::mutate(
    tMN = to01(tMN),
    AGE = as.numeric(AGE),
    SEX = factor(SEX),
    multipe_infusion_all = to01(multipe_infusion_all),
    Prior_Auto_fast = to01(Prior_Auto_fast),
    car_after_car = to01(car_after_car),
    Tumor = factor(Tumor)
  ) %>%
  dplyr::filter(!is.na(Time_to_Event_or_follow_up), Time_to_Event_or_follow_up > 90, !is.na(tMN)) %>%
  dplyr::mutate(
    SEX = forcats::fct_relevel(SEX, "F"),
    Tumor = forcats::fct_relevel(Tumor, "LBCL")
  )

model_input <- tmn_risk_model %>%
  dplyr::transmute(
    figure_subject_id,
    follow_up_days = Time_to_Event_or_follow_up,
    tMN,
    AGE,
    SEX = as.character(SEX),
    Tumor = as.character(Tumor),
    Prior_Auto_fast,
    multipe_infusion_all,
    car_after_car
  )
write_figshare_csv(model_input, "figshare_tMN_logistic_risk_model_input_D90.csv", "Deidentified follow-up >90 day model input for tMN logistic risk-factor analysis.")

candidate_vars <- c("AGE", "SEX", "Prior_Auto_fast", "multipe_infusion_all", "car_after_car")
candidate_vars <- candidate_vars[candidate_vars %in% names(tmn_risk_model)]

uni_tbl <- dplyr::bind_rows(lapply(candidate_vars, function(var) {
  df <- tmn_risk_model %>%
    dplyr::select(tMN, dplyr::all_of(var)) %>%
    dplyr::filter(!is.na(.data[[var]]))
  if (dplyr::n_distinct(df[[var]]) < 2) return(NULL)
  fit <- fit_logit_safe(as.formula(paste0("tMN ~ ", var)), df)
  tidy_or_ci(fit) %>%
    dplyr::mutate(
      predictor = var,
      label = vapply(term, function(tt) pretty_term(var, tt), character(1)),
      N = nrow(df),
      events = sum(df$tMN == 1)
    ) %>%
    dplyr::select(predictor, label, term, OR, CI_low, CI_high, p.value, method, N, events)
})) %>%
  dplyr::mutate(p_adj = p.adjust(p.value, method = "BH"), p_label = fmt_p(p.value)) %>%
  dplyr::arrange(p.value)

write_figshare_csv(uni_tbl, "figshare_tMN_logistic_univariable_D90_statistics.csv", "Univariable logistic odds ratios for tMN risk factors among subjects with >90 days follow-up.")
write.csv(uni_tbl, out_file("tMN_logistic_univariable_D90.csv"), row.names = FALSE)

mv_vars <- c("AGE", "SEX", "Prior_Auto_fast")
mv_vars <- mv_vars[mv_vars %in% names(tmn_risk_model)]
if (length(mv_vars) > 0) {
  mv_fit <- fit_logit_safe(as.formula(paste0("tMN ~ ", paste(mv_vars, collapse = " + "))), tmn_risk_model)
  mv_tbl <- tidy_or_ci(mv_fit) %>%
    dplyr::mutate(predictor = "Multivariable", label = term, p_label = fmt_p(p.value)) %>%
    dplyr::select(predictor, label, term, OR, CI_low, CI_high, p.value, p_label, method)
  write_figshare_csv(mv_tbl, "figshare_tMN_logistic_multivariable_D90_statistics.csv", "Conservative multivariable logistic odds ratios for tMN risk factors among subjects with >90 days follow-up.")
  write.csv(mv_tbl, out_file("tMN_logistic_multivariable_D90.csv"), row.names = FALSE)
}

plot_df <- uni_tbl %>%
  dplyr::mutate(
    label_clean = dplyr::case_when(
      predictor == "AGE" ~ "Age",
      predictor == "SEX" & stringr::str_detect(term, "M") ~ "Male Sex",
      predictor == "car_after_car" ~ "CAR after CAR",
      predictor == "Prior_Auto_fast" ~ "CAR after Auto",
      predictor == "multipe_infusion_all" ~ "Multiple Cell Therapies",
      TRUE ~ label
    ),
    OR = pmax(OR, .Machine$double.xmin),
    CI_low = pmax(CI_low, .Machine$double.xmin),
    CI_high = pmax(CI_high, .Machine$double.xmin),
    ord = dplyr::case_when(
      label_clean == "Age" ~ 1,
      label_clean == "Male Sex" ~ 2,
      label_clean == "CAR after CAR" ~ 3,
      label_clean == "CAR after Auto" ~ 4,
      label_clean == "Multiple Cell Therapies" ~ 5,
      TRUE ~ 99
    )
  ) %>%
  dplyr::filter(is.finite(OR), is.finite(CI_low), is.finite(CI_high)) %>%
  dplyr::arrange(ord, label_clean) %>%
  dplyr::mutate(label_clean = forcats::fct_rev(factor(label_clean)))

if (nrow(plot_df) > 0) {
  p_forest <- ggplot(plot_df, aes(y = label_clean, x = OR)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
    geom_point(size = 2) +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32), minor_breaks = NULL) +
    coord_cartesian(xlim = c(0.1, 50)) +
    labs(x = "Odds Ratio (log scale)", y = NULL, title = "tMN risk factors (univariable logistic; follow-up >90 days)") +
    theme_classic(base_size = 18, base_family = font_family)
  ggsave(out_file("tMN_logistic_univariable_forest_D90.pdf"), p_forest, width = 10, height = 7, device = cairo_pdf)
  ggsave(out_file("tMN_logistic_univariable_forest_D90.png"), p_forest, width = 10, height = 7, dpi = 300)
}

# =========================================================
# MANIFEST / README / SESSION INFO
# =========================================================
readr::write_csv(manifest, figshare_file("figshare_BM_oncoprint_flow_tMN_dataset_manifest.csv"), na = "")

readr::write_lines(
  c(
    "BM oncoprint, marrow flow, and tMN Figshare datasets",
    "",
    "These files were generated from the cleaned figure-only export script for the CAR cytopenia paper.",
    "Final CSVs are restricted to data used for figure generation and associated statistics.",
    "LYMID/AAID-style analysis identifiers are retained when available.",
    "",
    "Dataset manifest: figshare_BM_oncoprint_flow_tMN_dataset_manifest.csv"
  ),
  figshare_file("figshare_BM_oncoprint_flow_tMN_README.txt")
)

writeLines(capture.output(sessionInfo()), out_file("sessionInfo_BM_oncoprint_flow_tMN_figshare.txt"))

cat("\n[DONE] Figshare datasets written to:\n", figshare_dir, "\n", sep = "")
cat("[DONE] Figure outputs/session info written to:\n", out_dir, "\n", sep = "")
