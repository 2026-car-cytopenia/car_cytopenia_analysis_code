#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("clinical_ch_figures")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(survival)
  library(survminer)
  library(broom)
})

# ---------------------- Directories ----------------------
analysis_dir <- study_path("analysis_dir")
setwd(analysis_dir)

out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

figshare_dir <- file.path(out_dir, "figshare_datasets")
dir.create(figshare_dir, showWarnings = FALSE, recursive = TRUE)
figshare_file <- function(...) file.path(figshare_dir, ...)

# ---------------------- Helpers --------------------------
`%notin%` <- Negate(`%in%`)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

first_existing <- function(candidates, x) {
  hit <- intersect(candidates, x)
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_read_delim <- function(path, delim = NULL) {
  if (!file.exists(path)) stop("Missing required input file: ", path)
  if (is.null(delim)) {
    readr::read_delim(path, show_col_types = FALSE, guess_max = 100000)
  } else {
    readr::read_delim(path, delim = delim, show_col_types = FALSE, guess_max = 100000)
  }
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input file: ", path)
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

# tMN_CAR.txt is a one-column AAID file. read_delim() may fail because there is
# effectively no delimiter to infer, so read it line-by-line and standardize to AAID.
read_tmn_car_aai_file <- function(path) {
  if (!file.exists(path)) stop("Missing required input file: ", path)

  x <- readr::read_lines(path, skip_empty_rows = TRUE)
  x <- stringr::str_trim(x)
  x <- x[!is.na(x) & x != ""]

  # Drop a header line if present.
  if (length(x) > 0 && stringr::str_to_upper(x[[1]]) %in% c("AAID", "AA_ID", "LYMID", "PATIENT_ID")) {
    x <- x[-1]
  }

  tibble::tibble(AAID = x) %>%
    dplyr::filter(!is.na(AAID), AAID != "")
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ NA_character_,
    p < 0.001   ~ "p < 0.001",
    TRUE        ~ paste0("p = ", signif(p, 3))
  )
}

raw_id_value_regex <- "\\b(CCT[0-9A-Za-z_-]*|S[0-9]{4}(?:[-o]?[0-9]+)?)\\b"
raw_id_name_regex  <- "(?i)(^|_)(mrn|sid|sample[._]?id|sample_id|trial_id|cct|cct_id)(_|$)"

strip_raw_id_columns <- function(dat) {
  raw_cols <- names(dat)[stringr::str_detect(names(dat), raw_id_name_regex)]
  raw_cols <- setdiff(raw_cols, names(dat)[stringr::str_detect(names(dat), "^figshare_")])
  if (length(raw_cols) == 0) return(dat)
  dat %>% dplyr::select(-dplyr::all_of(raw_cols))
}

assert_no_raw_ids <- function(dat, filename) {
  dat_chr <- dat %>% mutate(across(everything(), as.character))
  bad_cols <- names(dat_chr)[vapply(dat_chr, function(x) {
    any(str_detect(x[!is.na(x)], raw_id_value_regex))
  }, logical(1))]
  if (length(bad_cols) > 0) {
    stop(
      "Disallowed internal identifiers detected in Figshare export ", filename,
      " in column(s): ", paste(bad_cols, collapse = ", "),
      ". Remove or deidentify CCT/S/MRN/sample-level identifiers before export."
    )
  }
  invisible(TRUE)
}

write_figshare_csv <- function(dat, filename) {
  dat2 <- dat %>% strip_raw_id_columns()
  assert_no_raw_ids(dat2, filename)
  readr::write_csv(dat2, figshare_file(filename), na = "")
  invisible(dat2)
}

save_plot_all <- function(plot, stem, width, height, dpi = 400) {
  ggsave(out_file(paste0(stem, ".pdf")), plot = plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(out_file(paste0(stem, ".eps")), plot = plot, width = width, height = height, units = "in", device = cairo_ps)
  ggsave(out_file(paste0(stem, ".png")), plot = plot, width = width, height = height, units = "in", dpi = dpi)
}

# Create stable deidentified patient IDs within this script.
make_patient_map <- function(...) {
  ids <- unlist(list(...), use.names = FALSE)
  ids <- sort(unique(as.character(ids[!is.na(ids) & ids != ""])))
  tibble(
    raw_patient_id = ids,
    figshare_patient_id = sprintf("Patient_%03d", seq_along(ids))
  )
}

add_deid_patient <- function(dat, raw_col, patient_map) {
  if (!raw_col %in% names(dat)) return(dat)
  dat %>%
    left_join(patient_map, by = setNames("raw_patient_id", raw_col))
}

add_deid_variant <- function(dat, raw_col = "VARIANT_ID", prefix = "Variant") {
  if (!raw_col %in% names(dat)) return(dat)
  vmap <- dat %>%
    distinct(.raw_variant = .data[[raw_col]]) %>%
    filter(!is.na(.raw_variant), .raw_variant != "") %>%
    arrange(.raw_variant) %>%
    mutate(figshare_variant_id = sprintf(paste0(prefix, "_%03d"), row_number()))
  dat %>% left_join(vmap, by = setNames(".raw_variant", raw_col))
}

# ---------------------- Domain definitions ----------------
functional <- c(
  "missense", "stopgain", "startloss", "splicing", "splice donor",
  "stoploss", "frameshift", "duplication", "deletion", "splice acceptor"
)

CH <- c(
  "ASXL1", "CBL", "CHEK2", "DNMT3A", "GNAS", "GNB1", "IDH1", "IDH2",
  "JAK2", "MYD88", "PPM1D", "SRSF2", "STAT3", "TET2", "U2AF1", "TP53", "SF3B1"
)

favorite_CH_names <- c("DNMT3A", "TP53", "PPM1D", "TET2", "ASXL1", "IDH2", "ATM", "U2AF1", "EZH2", "SF3B1")

nc_types <- c(
  "ncRNA_exonic", "ncRNA_intronic", "intergenic", "intronic", "synonymous",
  "upstream", "upstream;downstream", "UTR5", "UTR3", "downstream"
)

# Nick/pre-infusion CH list used in the paper cytopenia figure.
nick_list <- study_setting("nick_list")

all_pre_ch <- study_setting("all_pre_ch")

# ---------------------- Inputs ---------------------------
# Note: tMN_CAR.txt is read with a custom one-column AAID reader because
# readr::read_delim() cannot always infer a delimiter for a single-column file.
all_variants_raw <- safe_read_delim(study_path("input_all_snv_indel_monitoring"))
pre_inf_variants_raw <- safe_read_delim(study_path("input_v16_raw_calls"))
just_CH_mutations_raw <- safe_read_delim(study_path("input_v16_justch_mutations_snvs_indels"))
bm_variants_raw <- safe_read_delim(study_path("input_bm_oncoprint_snvs_4_update_20241011_with_mm_format_for_merge"))
BM_vs_raw <- safe_read_delim(study_path("input_bm_vs_cfdna_vs_pbl"))
tMN_molecular_raw <- safe_read_delim(study_path("input_tmn_molecular_trace"))
preCH_raw <- safe_read_delim(study_path("input_v15_precar_ch_list"))
cytopenia_raw <- safe_read_delim(study_path("input_cytopenia_list"))
pre_CH_data_raw <- safe_read_delim(study_path("input_pre_car_ch_workbook"))
tMN_car_raw <- read_tmn_car_aai_file(study_path("input_tmn_car"))
expansion_data_raw <- safe_read_delim(study_path("input_ch_expand_eval"))
tmn_CAR_raw <- safe_read_delim(study_path("input_tmn_axicel_depth_workbook"))

# Optional; used only for matched pre/post CH summary if available.
match_var_raw <- if (file.exists(study_path("input_pretx_precar_match", optional = TRUE))) safe_read_delim(study_path("input_pretx_precar_match", optional = TRUE)) else tibble()

# Build deidentification map from every raw ID source used in joins/filtering.
tMN_car_id_col <- first_existing(c("AAID", "AA_ID", "AA_Lab", "patient_id", "LYMID"), names(tMN_car_raw))

patient_map <- make_patient_map(
  all_variants_raw$patient_id,
  bm_variants_raw$patient_id,
  pre_inf_variants_raw$patient_id,
  just_CH_mutations_raw$patient_id,
  tMN_molecular_raw$patient_id,
  preCH_raw$patient_id,
  cytopenia_raw$LYMID,
  pre_CH_data_raw$AA_Lab,
  if (!is.na(tMN_car_id_col)) tMN_car_raw[[tMN_car_id_col]],
  expansion_data_raw$Lab_ID,
  tmn_CAR_raw$Lab_ID,
  tmn_CAR_raw$patient_id,
  nick_list,
  all_pre_ch
)

# ---------------------- Core variant data ----------------
# Variant IDs are created for internal matching only and are not exported.
pre_inf_variants <- pre_inf_variants_raw %>%
  filter(sample_type != "cfDNA") %>%
  mutate(
    POS = as.character(POS),
    VARIANT_ID = paste(patient_id, GENE, CHR, POS, paste0(REF, ">", TUMOR), sep = ":")
  ) %>%
  filter(TYPE %in% functional, GENE %in% CH)
pre_vars <- pre_inf_variants$VARIANT_ID

just_CH_mutations <- just_CH_mutations_raw %>%
  mutate(
    POS = as.character(POS),
    VARIANT_ID = paste(patient_id, GENE, CHR, POS, paste0(REF, ">", TUMOR), sep = ":")
  )

all_variants <- all_variants_raw %>%
  mutate(POS = as.character(POS)) %>%
  bind_rows(bm_variants_raw %>% mutate(POS = as.character(POS))) %>%
  mutate(
    VARIANT = paste(GENE, CHR, POS, paste0(REF, ">", TUMOR), sep = ":"),
    VARIANT_ID = paste(patient_id, GENE, CHR, POS, paste0(REF, ">", TUMOR), sep = ":"),
    VARIANT_day = paste(patient_id, GENE, CHR, POS, paste0(REF, ">", TUMOR), sample_day, sep = ":"),
    Type = case_when(
      sample_type %in% c("PBMC", "Normal") ~ "PBL",
      TRUE ~ as.character(sample_type)
    )
  ) %>%
  filter(TYPE %in% functional, GENE %in% CH)

# =========================================================
# FIGURE DATASET 1: BM versus cfDNA/PBL VAF correlation
# =========================================================
BM_vs_export <- BM_vs_raw %>%
  mutate(data_point_id = sprintf("BMcorr_%03d", row_number())) %>%
  transmute(
    data_point_id,
    gene = if ("GENE" %in% names(.)) as.character(GENE) else NA_character_,
    bone_marrow_vaf_percent = as_num(AF),
    cfdna_vaf_percent = as_num(cfDNA),
    pbl_vaf_percent = as_num(PBL)
  )
write_figshare_csv(BM_vs_export, "figshare_BM_cfDNA_PBL_correlation_plot_data.csv")

p_BM_cfDNA <- ggpubr::ggscatter(
  BM_vs_export %>% drop_na(cfdna_vaf_percent),
  x = "bone_marrow_vaf_percent", y = "cfdna_vaf_percent",
  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
  xlab = "VAF in Bone Marrow", ylab = "VAF in cfDNA"
) + ggtitle("Bone marrow versus cfDNA") + theme_classic(base_size = 20, base_family = "Arial")

p_BM_PBL <- ggpubr::ggscatter(
  BM_vs_export %>% drop_na(pbl_vaf_percent),
  x = "bone_marrow_vaf_percent", y = "pbl_vaf_percent",
  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
  xlab = "VAF in Bone Marrow", ylab = "VAF in PBL"
) + ggtitle("Bone marrow versus peripheral blood leukocytes") + theme_classic(base_size = 20, base_family = "Arial")

save_plot_all(p_BM_cfDNA + p_BM_PBL, "BM_cfDNA_PBL_correlation", width = 12, height = 6)

# =========================================================
# FIGURE DATASET 2: tMN mutation change over time
# =========================================================
tMN_molecular <- tMN_molecular_raw %>%
  filter(patient_id != study_setting("study_value_001")) %>%
  add_deid_patient("patient_id", patient_map) %>%
  add_deid_variant("VARIANT_ID", prefix = "tMN_variant") %>%
  transmute(
    LYMID = as.character(patient_id),
    figshare_patient_id,
    figshare_variant_id,
    gene = as.character(GENE),
    timepoint = as.character(Timpoint),
    vaf_percent = as_num(firstAF),
    sample_type = as.character(Type)
  )
write_figshare_csv(tMN_molecular, "figshare_tMN_mutation_change_over_time_plot_data.csv")

gene_cols <- c(
  "DNMT3A" = "#0072B2", "TET2" = "#009E73", "TP53" = "#C43A31",
  "PPM1D" = "#CC79A7", "ASXL1" = "#999999", "ATM" = "#56B4E9",
  "EZH2" = "#000000", "IDH1" = "#F0E442", "JAK2" = "#A6761D"
)

p_tMN_trace <- ggplot(tMN_molecular, aes(x = timepoint, y = vaf_percent, group = figshare_variant_id, color = gene)) +
  geom_line(linewidth = 1.2, alpha = 0.65) +
  geom_point(size = 2.5, aes(shape = sample_type)) +
  scale_color_manual(values = gene_cols, na.value = "grey50") +
  theme_classic(base_family = "Arial", base_size = 20) +
  xlab("Timepoint") + ylab("VAF") +
  ggtitle("tMN mutation change over time")
save_plot_all(p_tMN_trace, "tMN_mutation_change_over_time", width = 7, height = 5)

# =========================================================
# FIGURE DATASET 3: pre-infusion cfDNA versus PBL CH correlation
# =========================================================
pre_CH_cfdna_pbl <- all_variants %>%
  filter(VARIANT_ID %in% pre_vars, GENE %in% CH, TYPE %notin% nc_types, sample_day < 1) %>%
  filter(Type %in% c("PBL", "cfDNA")) %>%
  add_deid_patient("patient_id", patient_map) %>%
  add_deid_variant("VARIANT_ID", prefix = "preCH_variant") %>%
  select(LYMID = patient_id, figshare_patient_id, figshare_variant_id, gene = GENE, VARIANT_day, Type, AF) %>%
  distinct()

cfDNA_pre <- pre_CH_cfdna_pbl %>%
  filter(Type == "cfDNA") %>%
  select(VARIANT_day, LYMID, figshare_patient_id, figshare_variant_id, gene, cfdna_vaf_percent = AF)
PBL_pre <- pre_CH_cfdna_pbl %>%
  filter(Type == "PBL") %>%
  select(VARIANT_day, pbl_vaf_percent = AF)

cfDNA_PBL_pre_merge <- cfDNA_pre %>%
  inner_join(PBL_pre, by = "VARIANT_day") %>%
  mutate(
    log10_cfdna_vaf_plus_0.01 = log10(as_num(cfdna_vaf_percent) + 0.01),
    log10_pbl_vaf_plus_0.01   = log10(as_num(pbl_vaf_percent) + 0.01)
  ) %>%
  select(-VARIANT_day)
write_figshare_csv(cfDNA_PBL_pre_merge, "figshare_preinfusion_cfDNA_PBL_CH_correlation_plot_data.csv")

p_pre_corr <- ggpubr::ggscatter(
  cfDNA_PBL_pre_merge,
  x = "log10_cfdna_vaf_plus_0.01", y = "log10_pbl_vaf_plus_0.01",
  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
  add.params = list(color = "black", fill = "lightgray"),
  xlab = "Log10 VAF in cfDNA", ylab = "Log10 VAF in PBL"
) +
  theme_classic(base_size = 24, base_family = "Arial") +
  geom_hline(yintercept = -2, linetype = "dashed") +
  geom_vline(xintercept = -2, linetype = "dashed") +
  ggtitle("PBL vs cfDNA pre-infusion\nCH mutation correlation")
save_plot_all(p_pre_corr, "preinfusion_cfDNA_PBL_CH_correlation", width = 7, height = 6)

# =========================================================
# FIGURE DATASET 4: cfDNA CH first-to-last and delta VAF
# =========================================================
full_dataset <- all_variants %>%
  filter(sample_type == "cfDNA", TYPE %notin% nc_types) %>%
  arrange(VARIANT_ID, sample_day) %>%
  group_by(patient_id, VARIANT, VARIANT_ID, GENE) %>%
  mutate(first_day = first(sample_day), last_day = last(sample_day)) %>%
  filter(sample_day == first_day | sample_day == last_day) %>%
  filter(first_day < 1, last_day > 25) %>%
  mutate(last_VAF = last(AF), first_VAF = first(AF)) %>%
  distinct(VARIANT_ID, sample_day, .keep_all = TRUE) %>%
  select(patient_id, VARIANT_ID, GENE, first_day, last_day, first_VAF, last_VAF) %>%
  ungroup() %>%
  distinct() %>%
  filter(VARIANT_ID %in% just_CH_mutations$VARIANT_ID) %>%
  mutate(
    delta_vaf_percent = as_num(last_VAF) - as_num(first_VAF),
    log10_fold_change_plus_0.01 = log10((as_num(last_VAF) + 0.01) / (as_num(first_VAF) + 0.01))
  ) %>%
  add_deid_patient("patient_id", patient_map) %>%
  add_deid_variant("VARIANT_ID", prefix = "cfDNA_CH_variant")

first_last_wide <- full_dataset %>%
  transmute(
    LYMID = patient_id,
    figshare_patient_id,
    figshare_variant_id,
    gene = GENE,
    first_day_post_infusion = as_num(first_day),
    last_day_post_infusion  = as_num(last_day),
    first_vaf_percent = as_num(first_VAF),
    last_vaf_percent  = as_num(last_VAF),
    delta_vaf_percent,
    log10_fold_change_plus_0.01
  )
write_figshare_csv(first_last_wide, "figshare_cfDNA_first_last_CH_variant_plot_data_wide.csv")

first_last_long <- first_last_wide %>%
  select(LYMID, figshare_patient_id, figshare_variant_id, gene, first_vaf_percent, last_vaf_percent) %>%
  pivot_longer(c(first_vaf_percent, last_vaf_percent), names_to = "timepoint", values_to = "vaf_percent") %>%
  mutate(timepoint = recode(timepoint, first_vaf_percent = "First", last_vaf_percent = "Last"))
write_figshare_csv(first_last_long, "figshare_cfDNA_first_last_CH_variant_plot_data_long.csv")

zero_delta_long <- first_last_wide %>%
  select(LYMID, figshare_patient_id, figshare_variant_id, gene, delta_vaf_percent) %>%
  mutate(zero = 0) %>%
  pivot_longer(c(zero, delta_vaf_percent), names_to = "timepoint", values_to = "vaf_percent") %>%
  mutate(timepoint = recode(timepoint, zero = "Pre-infusion baseline", delta_vaf_percent = "Post-infusion delta"))
write_figshare_csv(zero_delta_long, "figshare_cfDNA_delta_CH_variant_plot_data_long.csv")

make_delta_plot <- function(gene_name, y_lim = NULL) {
  p <- ggplot(zero_delta_long %>% filter(gene == gene_name), aes(x = timepoint, y = vaf_percent, group = figshare_variant_id)) +
    geom_line(linewidth = 1.1, alpha = 0.55, color = gene_cols[[gene_name]] %||% "grey40") +
    geom_point(size = 2.5, color = gene_cols[[gene_name]] %||% "grey40") +
    theme_classic(base_family = "Arial", base_size = 20) +
    xlab("") + ylab("cfDNA delta VAF") + ggtitle(gene_name) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
  if (!is.null(y_lim)) p <- p + coord_cartesian(ylim = y_lim)
  p
}

p_delta <- make_delta_plot("TP53", c(-50, 100)) + make_delta_plot("PPM1D", c(-50, 100)) + make_delta_plot("DNMT3A", c(-50, 100))
save_plot_all(p_delta, "cfDNA_delta_CH_variant_panel", width = 15, height = 5)

# =========================================================
# FIGURE DATASET 5: TP53 transformation example traces
# =========================================================
# These are the TP53 example patients used in the legacy figure code.  Raw IDs are
# retained only for filtering and are exported as deidentified case labels.
tp53_example_ids <- study_setting("tp53_example_ids")
tp53_case_map <- tibble(
  patient_id = tp53_example_ids,
  case_id = sprintf("TP53_case_%02d", seq_along(tp53_example_ids))
)

tp53_examples <- all_variants %>%
  filter(patient_id %in% tp53_example_ids, GENE == "TP53") %>%
  left_join(tp53_case_map, by = "patient_id") %>%
  add_deid_variant("VARIANT_ID", prefix = "TP53_variant") %>%
  transmute(
    LYMID = patient_id,
    case_id,
    figshare_variant_id,
    gene = GENE,
    day_post_infusion = as_num(sample_day),
    vaf_percent = as_num(AF),
    sample_type = Type
  ) %>%
  distinct()
write_figshare_csv(tp53_examples, "figshare_TP53_transformation_example_trace_plot_data.csv")

p_tp53_examples <- ggplot(tp53_examples, aes(x = day_post_infusion, y = vaf_percent, group = figshare_variant_id, color = figshare_variant_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = sample_type), size = 1.8) +
  facet_wrap(~ case_id, scales = "free", ncol = 2) +
  theme_classic(base_family = "Arial", base_size = 16) +
  xlab("Day Post Infusion") + ylab("CH VAF (percent)") +
  guides(color = "none")
save_plot_all(p_tp53_examples, "TP53_transformation_example_traces", width = 10, height = 8)

# =========================================================
# FIGURE DATASET 6: PFS by pre-infusion CH
# =========================================================
pre_CH_data <- pre_CH_data_raw %>%
  filter(Protocol %notin% c("5908", "5910", "5022", "NA"), AA_Lab != "NA", AA_Lab %in% all_pre_ch) %>%
  mutate(
    preinfusion_CH = if_else(AA_Lab %in% nick_list, "Yes", "No"),
    preinfusion_CH = factor(preinfusion_CH, levels = c("No", "Yes"))
  ) %>%
  add_deid_patient("AA_Lab", patient_map)

survival_export <- pre_CH_data %>%
  transmute(
    LYMID = as.character(AA_Lab),
    figshare_patient_id,
    preinfusion_CH = as.character(preinfusion_CH),
    pfs_day = as_num(PFS_day),
    pfs_event = as.integer(PFS),
    age = as_num(Age),
    sex = as.character(Sex),
    prior_lines = as_num(Prior_Lines)
  )
write_figshare_csv(survival_export, "figshare_preinfusion_CH_PFS_plot_data.csv")

p_surv <- survminer::ggsurvplot(
  fit = survival::survfit(Surv(pfs_day, pfs_event) ~ preinfusion_CH, data = survival_export),
  data = survival_export,
  xlab = "Months Post-Infusion",
  ylab = "Progression Free Survival",
  conf.int = FALSE,
  pval = TRUE,
  xscale = "d_m",
  surv.median.line = "hv",
  palette = c("#7FC97F", "#BEAED4"),
  risk.table = TRUE,
  risk.table.col = "strata",
  break.time.by = 365.25 / 4,
  risk.table.height = 0.35,
  xlim = c(0, 890),
  risk.table.y.text = FALSE,
  legend = "top",
  legend.labs = c("No", "Yes"),
  ggtheme = theme_classic(base_size = 20, base_family = "Arial")
)
# ggsurvplot object needs arranged output for reliable saving.
ggsurv <- survminer::arrange_ggsurvplots(list(p_surv), print = FALSE, ncol = 1, nrow = 1)
ggsave(out_file("preinfusion_CH_PFS.pdf"), plot = ggsurv, width = 8, height = 10, units = "in", device = cairo_pdf)
ggsave(out_file("preinfusion_CH_PFS.eps"), plot = ggsurv, width = 8, height = 10, units = "in", device = cairo_ps)
ggsave(out_file("preinfusion_CH_PFS.png"), plot = ggsurv, width = 8, height = 10, units = "in", dpi = 400)

# =========================================================
# FIGURE DATASET 7: cytopenia recovery by pre-infusion CH
# =========================================================
evaluable_list <- unique(preCH_raw$patient_id)

preCH_CH <- preCH_raw %>%
  filter(GENE %in% CH, TYPE %notin% nc_types, sample_day < 1, AF > 0.5, sample_type %in% c("Normal", "PBMC"))
preCH_CH_TP53 <- preCH_raw %>%
  filter(GENE == "TP53", TYPE %notin% c("intronic", "synonymous", "UTR3", "UTR5"), sample_day < 1, sample_type == "Normal")
preCH_fav <- preCH_raw %>%
  filter(GENE %in% favorite_CH_names, TYPE %notin% nc_types, sample_day < 1, sample_type %in% c("Normal", "PBMC"))
preCH_CH_high <- preCH_raw %>%
  filter(GENE %in% CH, TYPE %notin% nc_types, sample_day < 1, AF >= 2, sample_type %in% c("Normal", "PBMC"))

any_pre_CH <- unique(preCH_CH$patient_id)
any_pre_TP53 <- unique(preCH_CH_TP53$patient_id)
any_pre_favs <- unique(preCH_fav$patient_id)
any_high_CH <- unique(preCH_CH_high$patient_id)

cytopenia_eval <- cytopenia_raw %>%
  filter(LYMID %in% evaluable_list) %>%
  mutate(
    CH = if_else(LYMID %in% any_pre_CH, 1L, 0L),
    preTP53 = if_else(LYMID %in% any_pre_TP53, 1L, 0L),
    pre_favs = if_else(LYMID %in% any_pre_favs, 1L, 0L),
    high_CH = if_else(LYMID %in% any_high_CH, 1L, 0L),
    nick_list = if_else(LYMID %in% nick_list, 1L, 0L),
    Day = as.character(Day)
  ) %>%
  filter(Day %notin% c("365", "500", "800"), Protocol == "5906") %>%
  add_deid_patient("LYMID", patient_map)

cytopenia_individual_export <- cytopenia_eval %>%
  transmute(
    LYMID = as.character(LYMID),
    figshare_patient_id,
    day_post_infusion = as_num(Day),
    preinfusion_CH = if_else(nick_list == 1L, "Yes", "No"),
    ANC = as_num(ANC),
    ALC = as_num(ALC),
    hemoglobin = as_num(HG),
    platelets = as_num(PLT)
  )
write_figshare_csv(cytopenia_individual_export, "figshare_cytopenia_recovery_individual_plot_data.csv")

cytopenia_long <- cytopenia_individual_export %>%
  pivot_longer(c(ANC, ALC, hemoglobin, platelets), names_to = "lineage", values_to = "value") %>%
  drop_na(value)

cytopenia_summary <- cytopenia_long %>%
  group_by(day_post_infusion, preinfusion_CH, lineage) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )
write_figshare_csv(cytopenia_summary, "figshare_cytopenia_recovery_summary_plot_data.csv")

cytopenia_stats <- cytopenia_individual_export %>%
  filter(day_post_infusion %in% c(60, 90, 180)) %>%
  pivot_longer(c(ANC, ALC, hemoglobin, platelets), names_to = "lineage", values_to = "value") %>%
  drop_na(value) %>%
  group_by(day_post_infusion, lineage) %>%
  summarise(
    p_value = tryCatch(wilcox.test(value ~ preinfusion_CH)$p.value, error = function(e) NA_real_),
    p_label = fmt_p(p_value),
    .groups = "drop"
  )
write_figshare_csv(cytopenia_stats, "figshare_cytopenia_recovery_wilcoxon_statistics.csv")

lineage_titles <- c(
  ANC = "Neutrophils (10^9/L)", ALC = "Lymphocytes (10^9/L)",
  hemoglobin = "Hemoglobin (g/dL)", platelets = "Platelets (10^9/L)"
)
ref_lines <- c(ANC = 1.5, ALC = 1, hemoglobin = 12.6, platelets = 150)

make_cyto_plot <- function(lineage_name) {
  d <- cytopenia_summary %>% filter(lineage == lineage_name)
  ggplot(d, aes(x = factor(day_post_infusion), y = median, group = preinfusion_CH, color = preinfusion_CH)) +
    geom_vline(xintercept = which(sort(unique(d$day_post_infusion)) == 28), linetype = "dotted") +
    geom_hline(yintercept = ref_lines[[lineage_name]], linetype = "dotted") +
    geom_point(size = 2.2) +
    geom_line(linewidth = 1.1) +
    geom_errorbar(aes(ymin = median - se, ymax = median + se), width = 0.1, linewidth = 0.8) +
    scale_color_manual(values = c("No" = "#7FC97F", "Yes" = "#BEAED4"), name = "Pre-Infusion CH") +
    ylab("Median") + xlab("Day") + ggtitle(lineage_titles[[lineage_name]]) +
    theme_classic(base_size = 20, base_family = "Arial")
}

p_cyto <- make_cyto_plot("ANC") + make_cyto_plot("hemoglobin") + make_cyto_plot("platelets") + make_cyto_plot("ALC") + plot_layout(ncol = 4, guides = "collect")
save_plot_all(p_cyto, "cytopenia_recovery_preinfusion_CH", width = 28, height = 6)

# =========================================================
# FIGURE DATASET 8: CH expansion and D7 CAR expansion
# =========================================================
cfDNA_first_last_for_expansion <- all_variants %>%
  filter(sample_type == "cfDNA", TYPE %notin% nc_types) %>%
  arrange(VARIANT_ID, sample_day) %>%
  group_by(patient_id, VARIANT, VARIANT_ID, GENE) %>%
  mutate(first_day = first(sample_day), last_day = last(sample_day)) %>%
  filter(sample_day == first_day | sample_day == last_day) %>%
  filter(first_day < 1, last_day > 25) %>%
  mutate(last_VAF = last(AF), first_VAF = first(AF)) %>%
  distinct(VARIANT_ID, sample_day, .keep_all = TRUE) %>%
  ungroup() %>%
  filter(VARIANT_ID %in% just_CH_mutations$VARIANT_ID) %>%
  mutate(change = log10((as_num(last_VAF) + 0.01) / (as_num(first_VAF) + 0.01)))

df_max <- cfDNA_first_last_for_expansion %>%
  group_by(patient_id) %>%
  summarise(Max_expansion = max(change, na.rm = TRUE), .groups = "drop") %>%
  rename(Lab_ID = patient_id)

expansion_data_eval <- expansion_data_raw %>%
  left_join(df_max, by = "Lab_ID") %>%
  filter(Expand_eval == "Yes", Progression == "No", Sample_Type == "cfDNA", Day == 7) %>%
  mutate(axicel_log10_hge_per_ml = log10(as_num(Axicel) + 1)) %>%
  add_deid_patient("Lab_ID", patient_map)

expansion_export <- expansion_data_eval %>%
  transmute(
    LYMID = as.character(Lab_ID),
    figshare_patient_id,
    day_post_infusion = as_num(Day),
    sample_type = as.character(Sample_Type),
    axicel_hge_per_ml = as_num(Axicel),
    axicel_log10_hge_per_ml,
    max_CH_log10_fold_change_plus_0.01 = as_num(Max_expansion),
    CH_expansion_group = as.character(CH_expansion),
    progression = as.character(Progression)
  )
write_figshare_csv(expansion_export, "figshare_CH_expansion_vs_D7_CAR_expansion_plot_data.csv")

p_expansion_scatter <- ggpubr::ggscatter(
  expansion_export,
  x = "axicel_log10_hge_per_ml", y = "max_CH_log10_fold_change_plus_0.01",
  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
  xlab = "Axi-cel Expansion on D7 (hGE/mL)",
  ylab = "Maximum CH Fold Change",
  title = ""
) + geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_classic(base_size = 20, base_family = "Arial")
save_plot_all(p_expansion_scatter, "CH_expansion_vs_D7_CAR_expansion", width = 7, height = 6)

if ("CH_expansion" %in% names(expansion_data_eval)) {
  p_expansion_box <- ggplot(expansion_export, aes(y = axicel_hge_per_ml, x = CH_expansion_group, fill = CH_expansion_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.5, dotsize = 0.5) +
    ylab("CAR Expansion Day 7 (hGE/mL)") + xlab("") +
    theme_classic(base_size = 20, base_family = "Arial") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    theme(legend.position = "none")
  save_plot_all(p_expansion_box, "D7_CAR_expansion_by_CH_expansion_group", width = 6, height = 6)
}

# =========================================================
# FIGURE DATASET 9: Axi-cel detection around tMN
# =========================================================
tmn_CAR <- tmn_CAR_raw %>%
  mutate(
    tMN_timepoint = factor(tMN_timepoint, levels = c("Pre", "D7_14", "Other", "D85plus", "Y")),
    Axicel_hge = as_num(Axicel_hge),
    axicel_depth = as_num(axicel_depth)
  )

# Add LYMID/AAID-style ID and deidentified patient ID if an eligible raw ID column is present.
tmn_id_col <- NA_character_
for (id_col in c("Lab_ID", "patient_id", "AA_Lab", "AAID", "LYMID")) {
  if (id_col %in% names(tmn_CAR)) {
    tmn_id_col <- id_col
    tmn_CAR <- tmn_CAR %>% mutate(LYMID = as.character(.data[[id_col]]))
    tmn_CAR <- add_deid_patient(tmn_CAR, id_col, patient_map)
    break
  }
}
if (!"LYMID" %in% names(tmn_CAR)) {
  tmn_CAR <- tmn_CAR %>% mutate(LYMID = NA_character_)
}
if (!"figshare_patient_id" %in% names(tmn_CAR)) {
  tmn_CAR <- tmn_CAR %>% mutate(figshare_patient_id = NA_character_)
}

tmn_car_export <- tmn_CAR %>%
  filter(tMN_timepoint != "Other", Sample_Type %in% c("cfDNA", "PBL")) %>%
  transmute(
    LYMID,
    figshare_patient_id,
    sample_type = as.character(Sample_Type),
    tMN_timepoint = as.character(tMN_timepoint),
    axicel_hge_per_ml = Axicel_hge,
    axicel_depth = axicel_depth,
    axicel_hge_per_ml_plot_value = Axicel_hge + 0.001,
    axicel_depth_plot_value = axicel_depth + 0.1
  )
write_figshare_csv(tmn_car_export, "figshare_tMN_axicel_detection_plot_data.csv")

p_tmn_cfdna <- ggplot(tmn_car_export %>% filter(sample_type == "cfDNA"), aes(fill = tMN_timepoint, x = tMN_timepoint, y = axicel_hge_per_ml_plot_value)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Axi-cel Detection (hGE/mL)") + xlab("") +
  theme_classic(base_size = 20, base_family = "Arial") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "gray") +
  theme(legend.position = "none") +
  ggtitle("Axi-cel Detection in cfDNA")

p_tmn_pbl <- ggplot(tmn_car_export %>% filter(sample_type == "PBL"), aes(fill = tMN_timepoint, x = tMN_timepoint, y = axicel_depth_plot_value)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Depth of Axi-cel Coverage") + xlab("") +
  theme_classic(base_size = 20, base_family = "Arial") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray") +
  theme(legend.position = "none") +
  ggtitle("Axi-cel Detection in PBL")

save_plot_all(p_tmn_cfdna + p_tmn_pbl, "tMN_axicel_detection", width = 18, height = 8)

# =========================================================
# FIGURE DATASET 10: pre-chemo versus post-chemo CH frequencies
# =========================================================
# These are the paper counts hard-coded in the legacy script for the unmatched
# pre-chemo and post-chemo comparison shown in Figure 3f/g.
pre_post_counts <- tribble(
  ~analysis_set, ~gene,    ~cohort,       ~CH_no, ~CH_yes,
  "unmatched",  "Any",    "Pre-Chemo",      78,      46,
  "unmatched",  "Any",    "Post-Chemo",     39,      62,
  "unmatched",  "ASXL1",  "Pre-Chemo",     124,       0,
  "unmatched",  "ASXL1",  "Post-Chemo",     98,       3,
  "unmatched",  "DNMT3A", "Pre-Chemo",      97,      27,
  "unmatched",  "DNMT3A", "Post-Chemo",     71,      30,
  "unmatched",  "PPM1D",  "Pre-Chemo",     119,       5,
  "unmatched",  "PPM1D",  "Post-Chemo",     64,      37,
  "unmatched",  "TP53",   "Pre-Chemo",     115,       9,
  "unmatched",  "TP53",   "Post-Chemo",     82,      19,
  "unmatched",  "TET2",   "Pre-Chemo",     109,      15,
  "unmatched",  "TET2",   "Post-Chemo",     88,      13,
  "matched",    "Any",    "Pre-Chemo",      58,      43,
  "matched",    "Any",    "Post-Chemo",     39,      62,
  "matched",    "PPM1D",  "Pre-Chemo",      96,       5,
  "matched",    "PPM1D",  "Post-Chemo",     67,      34,
  "matched",    "TP53",   "Pre-Chemo",      92,       9,
  "matched",    "TP53",   "Post-Chemo",     84,      17
) %>%
  mutate(
    n_total = CH_no + CH_yes,
    frequency_percent = 100 * CH_yes / n_total
  )

pre_post_stats <- pre_post_counts %>%
  group_by(analysis_set, gene) %>%
  group_modify(function(.x, .y) {
    tab <- as.matrix(.x[, c("CH_no", "CH_yes")])
    rownames(tab) <- .x$cohort
    tibble(p_value = tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_))
  }) %>%
  ungroup() %>%
  mutate(p_label = fmt_p(p_value))
write_figshare_csv(pre_post_counts, "figshare_pre_post_chemo_CH_frequency_plot_data.csv")
write_figshare_csv(pre_post_stats, "figshare_pre_post_chemo_CH_frequency_fisher_statistics.csv")

p_freq <- pre_post_counts %>%
  filter(analysis_set == "unmatched", gene != "Any") %>%
  mutate(gene = factor(gene, levels = c("DNMT3A", "PPM1D", "TP53", "TET2", "ASXL1"))) %>%
  ggplot(aes(x = frequency_percent, y = gene, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  theme_classic(base_size = 16, base_family = "Arial") +
  xlab("Frequency (% of patients with mutation)") + ylab("CH Gene Mutations in LBCL") +
  scale_fill_manual(values = c("Pre-Chemo" = "grey75", "Post-Chemo" = "#CDB4DB")) +
  ggtitle("More CH mutations after chemo")
save_plot_all(p_freq, "pre_post_chemo_CH_frequency", width = 7, height = 5)

# =========================================================
# FIGURE DATASET 11: tMN odds ratios by gene mutation
# =========================================================
if (!is.na(tMN_car_id_col)) {
  tMN_ids <- unique(as.character(tMN_car_raw[[tMN_car_id_col]]))
} else {
  tMN_ids <- character(0)
}

evaluable <- tibble(raw_patient_id = unique(preCH_raw$patient_id)) %>%
  mutate(tMN = if_else(raw_patient_id %in% tMN_ids, 1L, 0L))

preCH_gene_status <- preCH_raw %>%
  filter(GENE %in% c("TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1"), TYPE %notin% nc_types, sample_day < 1, sample_type %in% c("Normal", "PBMC")) %>%
  distinct(raw_patient_id = patient_id, GENE) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0L)

for (g in c("TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1")) {
  if (!g %in% names(preCH_gene_status)) preCH_gene_status[[g]] <- 0L
}

gene_or_input <- evaluable %>%
  left_join(preCH_gene_status, by = "raw_patient_id") %>%
  mutate(across(c("TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1"), ~replace_na(.x, 0L)))

gene_or_stats <- purrr::map_dfr(c("TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1"), function(g) {
  f <- as.formula(paste("tMN ~", g))

  fit <- tryCatch(
    glm(f, data = gene_or_input, family = binomial),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      gene = g, odds_ratio = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p_value = NA_real_
    ))
  }

  tid <- tryCatch(
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL
  )

  if (is.null(tid) || !"p.value" %in% names(tid)) {
    return(tibble(
      gene = g, odds_ratio = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p_value = NA_real_
    ))
  }

  tid_g <- tid %>% dplyr::filter(term == g)

  if (nrow(tid_g) == 0) {
    return(tibble(
      gene = g, odds_ratio = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p_value = NA_real_
    ))
  }

  tid_g %>%
    dplyr::transmute(
      gene = g,
      odds_ratio = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p_value = p.value
    )
}) %>%
  dplyr::mutate(p_label = fmt_p(p_value))
write_figshare_csv(gene_or_stats, "figshare_tMN_gene_odds_ratio_statistics.csv")

p_gene_or <- ggplot(gene_or_stats, aes(x = odds_ratio, y = factor(gene, levels = rev(c("TP53", "DNMT3A", "PPM1D", "TET2", "ASXL1"))), color = gene)) +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.4) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.06) +
  scale_x_log10() +
  scale_color_manual(values = gene_cols, na.value = "grey50") +
  theme_classic(base_size = 16, base_family = "Arial") +
  xlab("Odds Ratio for tMN") + ylab("") +
  ggtitle("Odds ratio for myeloid malignancy by gene mutation") +
  theme(legend.position = "none")
save_plot_all(p_gene_or, "tMN_gene_odds_ratio", width = 6.5, height = 4.5)

# =========================================================
# Manifest + README
# =========================================================
manifest <- tibble(
  file = list.files(figshare_dir, pattern = "^figshare_.*\\.csv$", full.names = FALSE),
  description = case_when(
    str_detect(file, "BM_cfDNA_PBL") ~ "Bone marrow versus cfDNA/PBL VAF correlation plot data.",
    str_detect(file, "tMN_mutation_change") ~ "tMN mutation change over time plot data.",
    str_detect(file, "preinfusion_cfDNA_PBL") ~ "Pre-infusion cfDNA versus PBL CH VAF correlation plot data.",
    str_detect(file, "cfDNA_first_last") ~ "First-to-last cfDNA CH variant VAF plot data.",
    str_detect(file, "cfDNA_delta") ~ "Zeroed delta cfDNA CH variant VAF plot data.",
    str_detect(file, "TP53_transformation") ~ "TP53 malignant transformation example trace plot data.",
    str_detect(file, "preinfusion_CH_PFS") ~ "PFS by pre-infusion CH status plot data.",
    str_detect(file, "cytopenia_recovery_individual") ~ "Individual cytopenia recovery values used to generate longitudinal summaries.",
    str_detect(file, "cytopenia_recovery_summary") ~ "Median/SE cytopenia recovery summary plotted by CH group.",
    str_detect(file, "cytopenia_recovery_wilcoxon") ~ "Wilcoxon statistics for cytopenia recovery comparisons.",
    str_detect(file, "CH_expansion_vs_D7") ~ "D7 CAR expansion versus maximum CH expansion plot data.",
    str_detect(file, "tMN_axicel") ~ "Axi-cel detection around tMN plot data.",
    str_detect(file, "pre_post_chemo_CH_frequency_plot") ~ "Pre-chemo versus post-chemo CH frequency plot data.",
    str_detect(file, "pre_post_chemo_CH_frequency_fisher") ~ "Fisher exact statistics for pre/post chemo CH frequency comparisons.",
    str_detect(file, "tMN_gene_odds_ratio") ~ "Univariable tMN odds ratios by CH gene.",
    TRUE ~ "Figshare-ready deidentified figure dataset."
  )
)
readr::write_csv(manifest, figshare_file("figshare_PBMC_cfDNA_CH_tMN_dataset_manifest.csv"), na = "")

readme <- c(
  "PBMC/cfDNA CH, tMN, cytopenia, and CAR-expansion Figshare datasets",
  "===================================================================",
  "",
  "This directory contains figure-specific, deidentified datasets generated from the PBMC/cfDNA CH analysis script.",
  "LYMID/AAID-style analysis identifiers are retained when useful for linking figure datasets. CCT IDs, S IDs, MRNs, SID, Sample_ID, and Trial_ID are not exported.",
  "Longitudinal rows use LYMID when available, plus figshare_patient_id and figshare_variant_id labels generated within this script.",
  "",
  "Datasets include only rows retained by the analysis filters used for figure generation. Values excluded in code are not retained in the final upload files.",
  "",
)
writeLines(readme, con = figshare_file("figshare_PBMC_cfDNA_CH_tMN_README.txt"))

writeLines(capture.output(sessionInfo()), con = out_file("PBMC_cfDNA_CH_tMN_figshare_sessionInfo.txt"))

cat("\nSaved PBMC/cfDNA CH/tMN figure datasets to:\n", figshare_dir, "\n", sep = "")
cat("Saved figure outputs to:\n", out_dir, "\n", sep = "")
