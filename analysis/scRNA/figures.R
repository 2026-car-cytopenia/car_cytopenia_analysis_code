#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("figures")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(Seurat)
  library(patchwork)
  library(ggh4x)
  library(ggrepel)
  library(ggpubr)
  library(ggprism)
  library(fgsea)
  library(msigdbr)
})

# ggrastr is optional but recommended for hybrid vector/raster UMAPs.
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
if (!has_ggrastr) {
  message("Package 'ggrastr' is not installed. UMAPs will still save, but points will not be rasterized. Install with install.packages('ggrastr').")
}

# =========================================================
# CONFIG
# =========================================================
base_scRNA_dir <- study_path("base_scRNA_dir")
out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

figshare_dir <- file.path(out_dir, "figshare_datasets", "scRNA_Fig2")
dir.create(figshare_dir, showWarnings = FALSE, recursive = TRUE)
figshare_file <- function(...) file.path(figshare_dir, ...)

hybrid_umap_dir <- file.path(out_dir, "umap_vector_raster_points", "scRNA_Fig2")
dir.create(hybrid_umap_dir, showWarnings = FALSE, recursive = TRUE)
hybrid_file <- function(...) file.path(hybrid_umap_dir, ...)

font_family <- "Arial"

# Timepoint / cohort order used everywhere
group_order_timepoint2 <- c(
  "Untreated Control PBL",
  "Untreated Control Marrow",
  "HCT PBL",
  "CAR PBL",
  "CAR Marrow"
)

group_order_marrow <- c(
  "Untreated Control Marrow",
  "CAR Marrow"
)

axis_trunc <- ggh4x::guide_axis_truncated(
  trunc_lower = grid::unit(0, "npc"),
  trunc_upper = grid::unit(3, "cm")
)

cols_lvl2 <- c(
  "T-cell"                     = "#B40F20",
  "Myeloid Cell"               = "#3288BD",
  "NK Cell"                    = "#D8A499",
  "Progenitor Cell"            = "#FFFFBF",
  "B-cell"                     = "#5E4FA2",
  "Megakaryocyte and Platelet" = "#ABDDA4",
  "Erythroid Cell"             = "#66C2A5",
  "Stromal Cell"               = "#FDAE61",
  "Plasma Cell"                = "#E6F598",
  "Other"                      = "grey80"
)

cols_lvl3 <- c(
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

timepoint_cols <- c(
  "Untreated Control PBL"    = "#66C2A5",
  "Untreated Control Marrow" = "#FC8D62",
  "HCT PBL"                  = "#8DA0CB",
  "CAR PBL"                  = "#E78AC3",
  "CAR Marrow"               = "#A6D854"
)

colorblind_vector <- grDevices::colorRampPalette(rev(c(
  "#0D0887FF", "#47039FFF", "#7301A8FF", "#9C179EFF",
  "#BD3786FF", "#D8576BFF", "#ED7953FF", "#FA9E3BFF",
  "#FDC926FF", "#F0F921FF"
)))

manifest <- tibble::tibble(
  file = character(),
  rows = integer(),
  columns = integer(),
  description = character()
)

# =========================================================
# HELPERS
# =========================================================
first_existing <- function(candidates, nm) {
  hit <- intersect(candidates, nm)
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

safe_wilcox_2group <- function(dat, value_col, group_col, g1, g2, paired = FALSE, id_col = NULL) {
  out_empty <- tibble::tibble(
    group1 = g1,
    group2 = g2,
    n_group1 = NA_integer_,
    n_group2 = NA_integer_,
    median_group1 = NA_real_,
    median_group2 = NA_real_,
    p_value = NA_real_,
    p_label = NA_character_,
    test = ifelse(paired, "paired Wilcoxon", "Wilcoxon rank-sum")
  )

  dat <- dat %>% dplyr::filter(.data[[group_col]] %in% c(g1, g2), !is.na(.data[[value_col]]))
  if (paired && !is.null(id_col)) {
    wide <- dat %>%
      dplyr::select(dplyr::all_of(c(id_col, group_col, value_col))) %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(group_col), values_from = dplyr::all_of(value_col)) %>%
      tidyr::drop_na(dplyr::all_of(c(g1, g2)))
    if (nrow(wide) < 2) return(out_empty)
    wt <- tryCatch(stats::wilcox.test(wide[[g1]], wide[[g2]], paired = TRUE, exact = FALSE), error = function(e) NULL)
    if (is.null(wt)) return(out_empty)
    return(tibble::tibble(
      group1 = g1,
      group2 = g2,
      n_group1 = nrow(wide),
      n_group2 = nrow(wide),
      median_group1 = stats::median(wide[[g1]], na.rm = TRUE),
      median_group2 = stats::median(wide[[g2]], na.rm = TRUE),
      p_value = wt$p.value,
      p_label = fmt_p(wt$p.value),
      test = "paired Wilcoxon"
    ))
  }

  x <- dat %>% dplyr::filter(.data[[group_col]] == g1) %>% dplyr::pull(dplyr::all_of(value_col))
  y <- dat %>% dplyr::filter(.data[[group_col]] == g2) %>% dplyr::pull(dplyr::all_of(value_col))
  if (length(x) < 1 || length(y) < 1) return(out_empty)
  wt <- tryCatch(stats::wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
  if (is.null(wt)) return(out_empty)
  tibble::tibble(
    group1 = g1,
    group2 = g2,
    n_group1 = length(x),
    n_group2 = length(y),
    median_group1 = stats::median(x, na.rm = TRUE),
    median_group2 = stats::median(y, na.rm = TRUE),
    p_value = wt$p.value,
    p_label = fmt_p(wt$p.value),
    test = "Wilcoxon rank-sum"
  )
}

remove_internal_id_columns <- function(df) {
  # Keep figshare_* identifiers, remove raw patient/sample/barcode/clinical IDs.
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
  # Blocks raw CCT IDs, S#### IDs, and MRN-like long numeric strings.
  bad_pattern <- "(CCT[0-9A-Za-z_-]+|\\bS[0-9]{4}(?:o[0-9]+|[-_][A-Za-z0-9]+)?\\b|\\b[0-9]{7,}\\b)"
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

write_figshare_csv_gz <- function(df, filename, description) {
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

get_umap_df <- function(obj, reduction = "umap") {
  emb <- tryCatch(Seurat::Embeddings(obj, reduction = reduction), error = function(e) NULL)
  if (is.null(emb)) {
    redn <- names(obj@reductions)
    hit <- redn[stringr::str_detect(redn, stringr::regex("umap", ignore_case = TRUE))]
    if (length(hit) == 0) stop("Could not find a UMAP reduction in object.")
    emb <- Seurat::Embeddings(obj, reduction = hit[[1]])
  }
  tibble::tibble(
    raw_cell_id = rownames(emb),
    UMAP1 = as.numeric(emb[, 1]),
    UMAP2 = as.numeric(emb[, 2])
  )
}

get_feature_vec <- function(obj, feature) {
  md <- obj@meta.data
  if (feature %in% names(md)) return(as.numeric(md[[feature]]))
  val <- tryCatch(
    Seurat::FetchData(obj, vars = feature)[[feature]],
    error = function(e) rep(NA_real_, ncol(obj))
  )
  as.numeric(val)
}

geom_point_raster_if_available <- function(mapping = NULL, data = NULL, ..., raster.dpi = 450, inherit.aes = TRUE) {
  if (has_ggrastr) {
    ggrastr::geom_point_rast(mapping = mapping, data = data, ..., raster.dpi = raster.dpi, inherit.aes = inherit.aes)
  } else {
    ggplot2::geom_point(mapping = mapping, data = data, ..., inherit.aes = inherit.aes)
  }
}

theme_umap <- function() {
  ggplot2::theme_classic(base_size = 20, base_family = font_family) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(0.3, "cm"), ends = "last")
      ),
      axis.line.x = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(0.3, "cm"), ends = "last")
      ),
      axis.title = ggplot2::element_text(hjust = 0)
    )
}

save_hybrid_plot <- function(plot, basename, width, height) {
  pdf_path <- hybrid_file(paste0(basename, "_hybrid_raster_points.pdf"))
  eps_path <- hybrid_file(paste0(basename, "_hybrid_raster_points.eps"))
  png_path <- hybrid_file(paste0(basename, "_hybrid_raster_points.png"))
  tryCatch(ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf), error = function(e) message("PDF save failed for ", basename, ": ", e$message))
  tryCatch(ggplot2::ggsave(eps_path, plot, width = width, height = height, units = "in", device = cairo_ps),  error = function(e) message("EPS save failed for ", basename, ": ", e$message))
  tryCatch(ggplot2::ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 400),       error = function(e) message("PNG save failed for ", basename, ": ", e$message))
  invisible(c(pdf_path, eps_path, png_path))
}


read_rds_optional <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    message("Optional RDS not found for ", label, ": ", path)
    return(NULL)
  }
  readRDS(path)
}

# =========================================================
# LOAD + ANNOTATE OBJECTS
# =========================================================
setwd(base_scRNA_dir)

combined <- readRDS(study_path("input_all_tcr_combined_contigs"))
Marrow_merge_idents <- readRDS(study_path("input_final_scrnaseq"))

# These marrow-only RDS files are optional. If one is absent, it is rebuilt
# below from the fully annotated Marrow_merge_idents object.
CAR_marrow_only  <- read_rds_optional(study_path("input_car_marrow_only", optional = TRUE),  "CAR_marrow_only")
CTRL_marrow_only <- read_rds_optional(study_path("input_ctrl_marrow_only", optional = TRUE), "CTRL_marrow_only")

mds_pt <- study_setting("mds_pt")
Marrow_merge_idents$MDS <- ifelse(Marrow_merge_idents$patient_id %in% mds_pt, "Yes", "No")

Marrow_merge_idents$timepoint_2 <- dplyr::case_when(
  Marrow_merge_idents$timepoint == "Control PBMC"   ~ "Untreated Control PBL",
  Marrow_merge_idents$timepoint == "Control Marrow" ~ "Untreated Control Marrow",
  Marrow_merge_idents$timepoint == "Auto PBMC"      ~ "HCT PBL",
  Marrow_merge_idents$timepoint == "CAR PBMC"       ~ "CAR PBL",
  Marrow_merge_idents$timepoint == "CAR Marrow"     ~ "CAR Marrow",
  TRUE                                              ~ as.character(Marrow_merge_idents$timepoint)
)
Marrow_merge_idents$timepoint_2 <- factor(Marrow_merge_idents$timepoint_2, levels = group_order_timepoint2)
Marrow_merge_idents$source <- ifelse(Marrow_merge_idents$timepoint_2 %in% c("Untreated Control Marrow", "CAR Marrow"), "Marrow", "PBL")

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
Marrow_merge_idents$predicted.celltype_2 <- factor(
  Marrow_merge_idents$predicted.celltype_2,
  levels = c("T-cell", "Myeloid Cell", "NK Cell", "Progenitor Cell", "B-cell", "Erythroid Cell", "Megakaryocyte and Platelet", "Plasma Cell", "Stromal Cell", "Other")
)

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
Marrow_merge_idents$predicted.celltype_3 <- factor(
  Marrow_merge_idents$predicted.celltype_3,
  levels = c("Naive T-cell", "Memory T-cell", "Effector T-cell", "Other T-cell", "Mature B-cell", "Immature B-cell", "Myeloid Cell", "NK Cell", "Progenitor Cell", "Erythroid Cell", "Megakaryocyte and Platelet", "Plasma Cell", "Stromal Cell", "Other")
)

DefaultAssay(Marrow_merge_idents) <- "SCT"
Marrow_merge_idents$axicel_expr <- get_feature_vec(Marrow_merge_idents, "axicel")

Marrow_merge_car_t <- subset(
  Marrow_merge_idents,
  subset = timepoint_2 %in% c("CAR PBL", "CAR Marrow") & predicted.celltype_2 == "T-cell"
)
Marrow_merge_car_t$axicel_expr <- get_feature_vec(Marrow_merge_car_t, "axicel")
Marrow_merge_car_t$CAR_group <- ifelse(Marrow_merge_car_t$axicel_expr > 0, "CAR+", "CAR-")
Marrow_merge_car_t$CAR_group <- factor(Marrow_merge_car_t$CAR_group, levels = c("CAR-", "CAR+"))

CAR_barcodes <- rownames(Marrow_merge_car_t)[Marrow_merge_car_t$CAR_group == "CAR+"]
Marrow_merge_idents$predicted.celltype_4 <- as.character(Marrow_merge_idents$predicted.celltype_2)
Marrow_merge_idents$predicted.celltype_4[rownames(Marrow_merge_idents) %in% CAR_barcodes] <- "CAR T-cell"
Marrow_merge_idents$predicted.celltype_4 <- factor(
  Marrow_merge_idents$predicted.celltype_4,
  levels = c("CAR T-cell", "T-cell", "Myeloid Cell", "NK Cell", "Progenitor Cell", "B-cell", "Erythroid Cell", "Megakaryocyte and Platelet", "Plasma Cell", "Stromal Cell", "Other")
)

# Rebuild missing marrow-only objects from the main annotated object.
# This avoids requiring 20241216_ctrl_marrow_only.rds to exist on disk.
if (is.null(CAR_marrow_only)) {
  message("Rebuilding CAR_marrow_only from Marrow_merge_idents")
  CAR_marrow_only <- subset(Marrow_merge_idents, subset = timepoint_2 == "CAR Marrow")
}
if (is.null(CTRL_marrow_only)) {
  message("Rebuilding CTRL_marrow_only from Marrow_merge_idents")
  CTRL_marrow_only <- subset(Marrow_merge_idents, subset = timepoint_2 == "Untreated Control Marrow")
}

# Carry the harmonized annotations into optional pre-saved subset objects when possible.
sync_subset_metadata <- function(sub_obj, main_obj) {
  common_cells <- intersect(rownames(sub_obj@meta.data), rownames(main_obj@meta.data))
  cols_to_copy <- intersect(
    c("timepoint_2", "predicted.celltype_2", "predicted.celltype_3", "predicted.celltype_4", "MDS", "source", "axicel_expr", "cloneType"),
    colnames(main_obj@meta.data)
  )
  if (length(common_cells) > 0 && length(cols_to_copy) > 0) {
    sub_obj@meta.data[common_cells, cols_to_copy] <- main_obj@meta.data[common_cells, cols_to_copy, drop = FALSE]
  }
  sub_obj
}
CAR_marrow_only <- sync_subset_metadata(CAR_marrow_only, Marrow_merge_idents)
CTRL_marrow_only <- sync_subset_metadata(CTRL_marrow_only, Marrow_merge_idents)

# =========================================================
# DEIDENTIFIED METADATA TABLES
# =========================================================
meta_raw <- Marrow_merge_idents@meta.data %>%
  tibble::rownames_to_column("raw_cell_id") %>%
  dplyr::mutate(
    patient_id = as.character(patient_id),
    sample_id = as.character(sample_id)
  ) %>%
  dplyr::left_join(get_umap_df(Marrow_merge_idents), by = "raw_cell_id")

patient_map <- meta_raw %>%
  dplyr::filter(!is.na(patient_id)) %>%
  dplyr::distinct(raw_patient_id = as.character(patient_id)) %>%
  dplyr::arrange(raw_patient_id) %>%
  dplyr::mutate(figshare_patient_id = paste0("scRNA_Patient_", sprintf("%03d", dplyr::row_number())))

sample_map <- meta_raw %>%
  dplyr::filter(!is.na(sample_id)) %>%
  dplyr::distinct(raw_sample_id = as.character(sample_id)) %>%
  dplyr::arrange(raw_sample_id) %>%
  dplyr::mutate(figshare_sample_id = paste0("scRNA_Sample_", sprintf("%03d", dplyr::row_number())))

cell_meta_public <- meta_raw %>%
  dplyr::left_join(patient_map, by = c("patient_id" = "raw_patient_id")) %>%
  dplyr::left_join(sample_map, by = c("sample_id" = "raw_sample_id")) %>%
  dplyr::mutate(figshare_cell_id = paste0("scRNA_Cell_", sprintf("%07d", dplyr::row_number()))) %>%
  dplyr::transmute(
    figshare_cell_id,
    figshare_patient_id,
    figshare_sample_id,
    UMAP1, UMAP2,
    timepoint_2 = as.character(timepoint_2),
    source,
    MDS,
    predicted.celltype.l2 = as.character(predicted.celltype.l2),
    predicted.celltype_2 = as.character(predicted.celltype_2),
    predicted.celltype_3 = as.character(predicted.celltype_3),
    predicted.celltype_4 = as.character(predicted.celltype_4),
    cloneType = if ("cloneType" %in% names(meta_raw)) as.character(cloneType) else NA_character_,
    CAR19_RNA_expression = axicel_expr,
    CAR_T_cell_status = ifelse(predicted.celltype_4 == "CAR T-cell", "CAR T-cell", "Non-CAR/other")
  )

write_figshare_csv_gz(
  cell_meta_public,
  "figshare_scRNA_Fig2_full_single_cell_metadata.csv.gz",
  "Deidentified cell-level metadata for the Fig 2 scRNA object, including UMAP coordinates and annotated cell types."
)

write_figshare_csv_gz(
  cell_meta_public %>% dplyr::select(figshare_cell_id, UMAP1, UMAP2, predicted.celltype_2),
  "figshare_scRNA_Fig2_celltype_umap_plot_data.csv.gz",
  "UMAP plot data for Fig 2 cell-type panel."
)

write_figshare_csv_gz(
  cell_meta_public %>% dplyr::select(figshare_cell_id, UMAP1, UMAP2, timepoint_2),
  "figshare_scRNA_Fig2_source_umap_plot_data.csv.gz",
  "UMAP plot data for Fig 2 sample-source/cohort panel."
)

write_figshare_csv_gz(
  cell_meta_public %>% dplyr::select(figshare_cell_id, UMAP1, UMAP2, CAR19_RNA_expression),
  "figshare_scRNA_Fig2_CAR19_expression_umap_plot_data.csv.gz",
  "UMAP plot data for CAR19 RNA expression feature panel."
)

# =========================================================
# HYBRID VECTOR/RASTER UMAPS
# =========================================================
p_celltype_hybrid <- ggplot(cell_meta_public, aes(UMAP1, UMAP2, color = predicted.celltype_2)) +
  geom_point_raster_if_available(size = 0.12, alpha = 0.85) +
  scale_color_manual(values = cols_lvl2, name = "Cell Type", na.value = "grey80") +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("Hematopoietic Cell Types") +
  theme_umap() +
  guides(x = axis_trunc, y = axis_trunc, color = guide_legend(override.aes = list(size = 4, alpha = 1)))
save_hybrid_plot(p_celltype_hybrid, "Fig2B_celltype_umap", width = 9, height = 7.06)

p_source_hybrid <- ggplot(cell_meta_public, aes(UMAP1, UMAP2, color = timepoint_2)) +
  geom_point_raster_if_available(size = 0.12, alpha = 0.85) +
  scale_color_manual(values = timepoint_cols, name = "Cohort", na.value = "grey80") +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("Sample Source") +
  theme_umap() +
  guides(x = axis_trunc, y = axis_trunc, color = guide_legend(override.aes = list(size = 4, alpha = 1)))

p_car19_hybrid <- ggplot(cell_meta_public, aes(UMAP1, UMAP2, color = CAR19_RNA_expression)) +
  geom_point_raster_if_available(size = 0.15, alpha = 0.9) +
  scale_color_gradient(low = "lightgrey", high = "darkgreen", name = "axicel") +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("CAR19 (axi-cel/brexu-cel) RNA Expression") +
  theme_umap() + guides(x = axis_trunc, y = axis_trunc)

save_hybrid_plot(p_source_hybrid | p_car19_hybrid, "Fig2C_source_and_CAR19", width = 12, height = 5.06)

# TCR clonality UMAPs from pre-saved marrow-only subsets
make_subset_umap_public <- function(obj, cohort_label) {
  md <- obj@meta.data %>% tibble::rownames_to_column("raw_cell_id") %>% dplyr::left_join(get_umap_df(obj), by = "raw_cell_id")
  pt_col <- first_existing(c("patient_id", "Patient_ID", "orig.ident"), names(md))
  sm_col <- first_existing(c("sample_id", "Sample_ID", "orig.ident"), names(md))
  md %>%
    dplyr::mutate(
      raw_patient_tmp = if (!is.na(pt_col)) as.character(.data[[pt_col]]) else NA_character_,
      raw_sample_tmp = if (!is.na(sm_col)) as.character(.data[[sm_col]]) else NA_character_
    ) %>%
    dplyr::left_join(patient_map, by = c("raw_patient_tmp" = "raw_patient_id")) %>%
    dplyr::left_join(sample_map, by = c("raw_sample_tmp" = "raw_sample_id")) %>%
    dplyr::mutate(
      figshare_cell_id = paste0("scRNA_", stringr::str_replace_all(cohort_label, "[^A-Za-z0-9]", "_"), "_Cell_", sprintf("%07d", dplyr::row_number())),
      cohort = cohort_label
    ) %>%
    dplyr::transmute(
      figshare_cell_id,
      figshare_patient_id,
      figshare_sample_id,
      UMAP1, UMAP2,
      cohort,
      cloneType = if ("cloneType" %in% names(md)) as.character(cloneType) else NA_character_,
      predicted.celltype_2 = if ("predicted.celltype_2" %in% names(md)) as.character(predicted.celltype_2) else NA_character_,
      predicted.celltype_3 = if ("predicted.celltype_3" %in% names(md)) as.character(predicted.celltype_3) else NA_character_
    )
}

clon_umap_public <- dplyr::bind_rows(
  make_subset_umap_public(CTRL_marrow_only, "Untreated Control Marrow"),
  make_subset_umap_public(CAR_marrow_only,  "CAR Marrow")
)

write_figshare_csv_gz(
  clon_umap_public,
  "figshare_scRNA_Fig2_TCR_clonality_marrow_umap_plot_data.csv.gz",
  "UMAP plot data for untreated control marrow versus post-CAR marrow clonotype category panels."
)

clone_cols <- stats::setNames(colorblind_vector(5), c(
  "Hyperexpanded (100 < X <= 500)",
  "Large (20 < X <= 100)",
  "Medium (5 < X <= 20)",
  "Small (1 < X <= 5)",
  "Single (0 < X <= 1)"
))

p_clon_hybrid <- ggplot(clon_umap_public, aes(UMAP1, UMAP2, color = cloneType)) +
  geom_point_raster_if_available(size = 0.12, alpha = 0.85) +
  facet_wrap(~cohort, nrow = 1) +
  scale_color_manual(values = clone_cols, na.value = "grey", name = "Clone type") +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("TCR clonality in marrow") +
  theme_umap() +
  guides(x = axis_trunc, y = axis_trunc, color = guide_legend(override.aes = list(size = 4, alpha = 1)))
save_hybrid_plot(p_clon_hybrid, "Fig2D_TCR_clonality_marrow", width = 14, height = 6.06)

# =========================================================
# COMPOSITION DATASETS + STATS
# =========================================================
cell_type_props_global <- cell_meta_public %>%
  dplyr::count(timepoint_2, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(timepoint_2) %>%
  dplyr::mutate(total_cells = sum(n_cells), fraction = n_cells / total_cells) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(timepoint_2 = factor(timepoint_2, levels = group_order_timepoint2))

write_figshare_csv(
  cell_type_props_global,
  "figshare_scRNA_Fig2_global_celltype_composition_all_cells.csv",
  "Global cell-type composition by cohort/timepoint used for stacked bar plot."
)

cell_type_props_sample <- cell_meta_public %>%
  dplyr::count(figshare_patient_id, figshare_sample_id, timepoint_2, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(figshare_sample_id) %>%
  dplyr::mutate(total_cells = sum(n_cells), fraction = n_cells / total_cells) %>%
  dplyr::ungroup() %>%
  dplyr::filter(timepoint_2 %in% group_order_marrow)

cell_type_props_summary <- cell_type_props_sample %>%
  dplyr::group_by(timepoint_2, predicted.celltype_3) %>%
  dplyr::summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    se_fraction = se(fraction),
    n_samples = dplyr::n_distinct(figshare_sample_id),
    .groups = "drop"
  )

cell_type_props_stats <- cell_type_props_sample %>%
  dplyr::group_by(predicted.celltype_3) %>%
  dplyr::group_modify(~safe_wilcox_2group(.x, "fraction", "timepoint_2", "Untreated Control Marrow", "CAR Marrow")) %>%
  dplyr::ungroup()

write_figshare_csv(
  cell_type_props_sample,
  "figshare_scRNA_Fig2_marrow_celltype_composition_per_sample.csv",
  "Per-sample marrow cell-type fractions used for mean±SE marrow composition plot."
)
write_figshare_csv(
  cell_type_props_summary,
  "figshare_scRNA_Fig2_marrow_celltype_composition_summary_mean_se.csv",
  "Mean and standard error of per-sample marrow cell-type fractions."
)
write_figshare_csv(
  cell_type_props_stats,
  "figshare_scRNA_Fig2_marrow_celltype_composition_wilcoxon_statistics.csv",
  "Wilcoxon statistics for per-sample marrow cell-type fractions, CAR marrow versus untreated control marrow."
)

# B-cell total fraction and immature fraction
b_total_fraction <- cell_meta_public %>%
  dplyr::filter(timepoint_2 %in% group_order_marrow) %>%
  dplyr::group_by(figshare_patient_id, figshare_sample_id, timepoint_2) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_B_cell = sum(predicted.celltype_3 %in% c("Mature B-cell", "Immature B-cell"), na.rm = TRUE),
    B_cell_fraction = n_B_cell / n_total,
    .groups = "drop"
  )

b_total_stats <- safe_wilcox_2group(b_total_fraction, "B_cell_fraction", "timepoint_2", "Untreated Control Marrow", "CAR Marrow")

b_comp_fraction <- cell_meta_public %>%
  dplyr::filter(timepoint_2 %in% group_order_marrow, predicted.celltype_3 %in% c("Mature B-cell", "Immature B-cell")) %>%
  dplyr::count(figshare_patient_id, figshare_sample_id, timepoint_2, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(figshare_sample_id) %>%
  dplyr::mutate(n_B_cells = sum(n_cells), B_subset_fraction = n_cells / n_B_cells) %>%
  dplyr::ungroup()

b_immature_stats <- b_comp_fraction %>%
  dplyr::filter(predicted.celltype_3 == "Immature B-cell") %>%
  safe_wilcox_2group("B_subset_fraction", "timepoint_2", "Untreated Control Marrow", "CAR Marrow")

b_pooled_counts <- cell_meta_public %>%
  dplyr::filter(timepoint_2 %in% group_order_marrow, predicted.celltype_3 %in% c("Mature B-cell", "Immature B-cell")) %>%
  dplyr::count(timepoint_2, predicted.celltype_3, name = "n_cells") %>%
  tidyr::complete(timepoint_2 = group_order_marrow, predicted.celltype_3 = c("Mature B-cell", "Immature B-cell"), fill = list(n_cells = 0L))

b_fisher_p <- tryCatch({
  mat <- b_pooled_counts %>% tidyr::pivot_wider(names_from = timepoint_2, values_from = n_cells) %>% tibble::column_to_rownames("predicted.celltype_3") %>% as.matrix()
  stats::fisher.test(mat)$p.value
}, error = function(e) NA_real_)
b_fisher_stats <- tibble::tibble(test = "Fisher exact test", p_value = b_fisher_p, p_label = fmt_p(b_fisher_p))

write_figshare_csv(b_total_fraction, "figshare_scRNA_Fig2_Bcell_total_fraction_per_sample.csv", "Per-sample total B-cell fraction among marrow cells.")
write_figshare_csv(b_total_stats, "figshare_scRNA_Fig2_Bcell_total_fraction_wilcoxon_statistics.csv", "Wilcoxon statistics for total B-cell fraction in marrow.")
write_figshare_csv(b_comp_fraction, "figshare_scRNA_Fig2_Bcell_subset_fraction_per_sample.csv", "Per-sample mature/immature B-cell fractions among marrow B cells.")
write_figshare_csv(b_immature_stats, "figshare_scRNA_Fig2_immature_Bcell_fraction_wilcoxon_statistics.csv", "Wilcoxon statistics for immature B-cell fraction among marrow B cells.")
write_figshare_csv(b_pooled_counts, "figshare_scRNA_Fig2_Bcell_subset_pooled_counts.csv", "Pooled mature/immature B-cell counts by marrow cohort.")
write_figshare_csv(b_fisher_stats, "figshare_scRNA_Fig2_Bcell_subset_fisher_statistics.csv", "Fisher exact statistics for pooled mature/immature B-cell counts by marrow cohort.")

# =========================================================
# PSEUDOBULK DE + GSEA
# =========================================================
DefaultAssay(Marrow_merge_idents) <- "RNA"
Marrow_marrow_only <- subset(Marrow_merge_idents, subset = timepoint_2 %in% group_order_marrow)

bulk_myeloid_de <- tibble::tibble()
bulk_prog_de <- tibble::tibble()
fgsea_myeloid <- tibble::tibble()
fgsea_prog <- tibble::tibble()
fgsea_sel <- tibble::tibble()

tryCatch({
  pseudo_marrow <- AggregateExpression(
    Marrow_marrow_only,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = c("timepoint_2", "patient_id", "predicted.celltype_4")
  )
  pseudo_marrow$celltype_timepoint <- paste(pseudo_marrow$predicted.celltype_4, pseudo_marrow$timepoint_2, sep = "_")
  Idents(pseudo_marrow) <- "celltype_timepoint"

  bulk_myeloid_de <- FindMarkers(
    object = pseudo_marrow,
    ident.1 = "Myeloid Cell_CAR Marrow",
    ident.2 = "Myeloid Cell_Untreated Control Marrow",
    test.use = "DESeq2"
  ) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      comparison = "Myeloid Cell: CAR Marrow vs Untreated Control Marrow",
      stat = sign(avg_log2FC) * (-log10(p_val + 1e-300)),
      updown = dplyr::case_when(
        p_val_adj < 0.05 & avg_log2FC > 0 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < 0 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )

  bulk_prog_de <- FindMarkers(
    object = pseudo_marrow,
    ident.1 = "Progenitor Cell_CAR Marrow",
    ident.2 = "Progenitor Cell_Untreated Control Marrow",
    test.use = "DESeq2"
  ) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      comparison = "Progenitor Cell: CAR Marrow vs Untreated Control Marrow",
      stat = sign(avg_log2FC) * (-log10(p_val + 1e-300)),
      updown = dplyr::case_when(
        p_val_adj < 0.05 & avg_log2FC > 0 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < 0 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )

  msig_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  hallmark_sets <- msig_h %>% split(x = .$gene_symbol, f = .$gs_name)

  ranks_myeloid <- bulk_myeloid_de$stat
  names(ranks_myeloid) <- bulk_myeloid_de$gene
  ranks_myeloid <- sort(ranks_myeloid[is.finite(ranks_myeloid)], decreasing = TRUE)

  ranks_prog <- bulk_prog_de$stat
  names(ranks_prog) <- bulk_prog_de$gene
  ranks_prog <- sort(ranks_prog[is.finite(ranks_prog)], decreasing = TRUE)

  fgsea_myeloid <- fgsea::fgsea(pathways = hallmark_sets, stats = ranks_myeloid, minSize = 15, maxSize = 500) %>%
    as_tibble() %>%
    dplyr::mutate(celltype = "Myeloid Cell", leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) %>%
    dplyr::arrange(padj)

  fgsea_prog <- fgsea::fgsea(pathways = hallmark_sets, stats = ranks_prog, minSize = 15, maxSize = 500) %>%
    as_tibble() %>%
    dplyr::mutate(celltype = "Progenitor Cell", leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) %>%
    dplyr::arrange(padj)

  pathways_of_interest <- c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING"
  )

  fgsea_sel <- bind_rows(fgsea_myeloid, fgsea_prog) %>%
    dplyr::filter(pathway %in% pathways_of_interest)
}, error = function(e) {
  message("Pseudobulk DE/GSEA block skipped or failed: ", e$message)
})

write_figshare_csv(bulk_myeloid_de, "figshare_scRNA_Fig2_myeloid_pseudobulk_DE.csv", "Pseudobulk differential expression results for myeloid cells, CAR marrow versus untreated control marrow.")
write_figshare_csv(bulk_prog_de, "figshare_scRNA_Fig2_progenitor_pseudobulk_DE.csv", "Pseudobulk differential expression results for progenitor cells, CAR marrow versus untreated control marrow.")
write_figshare_csv(fgsea_myeloid, "figshare_scRNA_Fig2_myeloid_hallmark_fgsea.csv", "Hallmark GSEA results for myeloid pseudobulk DE ranking.")
write_figshare_csv(fgsea_prog, "figshare_scRNA_Fig2_progenitor_hallmark_fgsea.csv", "Hallmark GSEA results for progenitor-cell pseudobulk DE ranking.")
write_figshare_csv(fgsea_sel, "figshare_scRNA_Fig2_selected_inflammatory_IFN_GSEA_plot_data.csv", "Selected inflammatory, IFN, NF-kB, and IL6/JAK/STAT Hallmark GSEA plot data.")

# =========================================================
# MODULE SCORES
# =========================================================
DefaultAssay(Marrow_merge_idents) <- "SCT"

AP1_genes <- c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "BATF", "BATF3")
NFkB_genes <- c("NFKB1", "NFKB2", "RELA", "RELB", "REL", "NFKBIA", "TNFAIP3", "TNFAIP2")
IFNa_response <- list(c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "FAM46A", "GBP2", "GBP4", "GMPR", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LY6E", "MX1", "NMI", "OAS1", "OASL", "PARP12", "PARP14", "PLSCR1", "PSMB8", "PSMB9", "PSME1", "PSME2", "RSAD2", "SAMD9L", "SELL", "STAT2", "TAP1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS"))
IFNy_response <- list(c("ADAR", "APOL6", "BATF2", "B2M", "C1S", "CASP1", "CCL2", "CCL5", "CD274", "CD40", "CD69", "CD74", "CD86", "CIITA", "CMPK2", "CXCL9", "CXCL10", "CXCL11", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EPSTI1", "FAS", "GBP4", "GBP6", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20", "MX1", "MX2", "NMI", "NLRC5", "OAS2", "OAS3", "OASL", "PARP12", "PARP14", "PLSCR1", "PSMB8", "PSMB9", "PSME1", "PSME2", "RSAD2", "RTP4", "SAMHD1", "SOCS1", "STAT1", "STAT2", "STAT4", "TAP1", "TAPBP", "TNFSF10", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TXNIP", "UBE2L6", "USP18", "VCAM1", "WARS", "XAF1", "ZBP1"))

if (!"AP1_score1" %in% colnames(Marrow_merge_idents@meta.data)) Marrow_merge_idents <- AddModuleScore(Marrow_merge_idents, features = list(AP1_genes), name = "AP1_score", random.seed = 1)
if (!"NFkB_score1" %in% colnames(Marrow_merge_idents@meta.data)) Marrow_merge_idents <- AddModuleScore(Marrow_merge_idents, features = list(NFkB_genes), name = "NFkB_score", random.seed = 1)
if (!"IFNa_response1" %in% colnames(Marrow_merge_idents@meta.data)) Marrow_merge_idents <- AddModuleScore(Marrow_merge_idents, features = IFNa_response, name = "IFNa_response", random.seed = 1)
if (!"IFNy_response1" %in% colnames(Marrow_merge_idents@meta.data)) Marrow_merge_idents <- AddModuleScore(Marrow_merge_idents, features = IFNy_response, name = "IFNy_response", random.seed = 1)

module_scores_raw <- Marrow_merge_idents@meta.data %>%
  tibble::rownames_to_column("raw_cell_id") %>%
  dplyr::mutate(
    patient_id = as.character(patient_id),
    sample_id = as.character(sample_id)
  ) %>%
  dplyr::left_join(patient_map, by = c("patient_id" = "raw_patient_id")) %>%
  dplyr::left_join(sample_map, by = c("sample_id" = "raw_sample_id")) %>%
  dplyr::mutate(figshare_cell_id = paste0("scRNA_Cell_", sprintf("%07d", match(raw_cell_id, meta_raw$raw_cell_id)))) %>%
  dplyr::select(figshare_cell_id, figshare_patient_id, figshare_sample_id, timepoint_2, predicted.celltype_2, AP1_score1, NFkB_score1, IFNa_response1, IFNy_response1) %>%
  dplyr::filter(timepoint_2 %in% group_order_marrow, predicted.celltype_2 %in% c("Myeloid Cell", "Progenitor Cell")) %>%
  tidyr::pivot_longer(cols = c(AP1_score1, NFkB_score1, IFNa_response1, IFNy_response1), names_to = "module", values_to = "score")

module_scores_sample <- module_scores_raw %>%
  dplyr::group_by(figshare_patient_id, figshare_sample_id, timepoint_2, predicted.celltype_2, module) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE), median_score = median(score, na.rm = TRUE), n_cells = dplyr::n(), .groups = "drop")

module_scores_stats <- module_scores_sample %>%
  dplyr::group_by(predicted.celltype_2, module) %>%
  dplyr::group_modify(~safe_wilcox_2group(.x, "mean_score", "timepoint_2", "Untreated Control Marrow", "CAR Marrow")) %>%
  dplyr::ungroup()

write_figshare_csv_gz(module_scores_raw, "figshare_scRNA_Fig2_module_scores_cell_level_long.csv.gz", "Cell-level AP1, NFkB, IFN-alpha, and IFN-gamma module scores in marrow myeloid/progenitor cells.")
write_figshare_csv(module_scores_sample, "figshare_scRNA_Fig2_module_scores_sample_summary.csv", "Per-sample module score summaries for marrow myeloid/progenitor cells.")
write_figshare_csv(module_scores_stats, "figshare_scRNA_Fig2_module_scores_wilcoxon_statistics.csv", "Wilcoxon statistics for module scores, CAR marrow versus untreated control marrow.")

# =========================================================
# T-CELL SUBSET ANALYSES: NON-CAR MARROW T CELLS
# =========================================================
t_subset_counts <- cell_meta_public %>%
  dplyr::filter(predicted.celltype_2 == "T-cell", predicted.celltype_4 != "CAR T-cell", timepoint_2 %in% group_order_marrow) %>%
  dplyr::filter(predicted.celltype_3 %in% c("Naive T-cell", "Memory T-cell", "Effector T-cell", "Other T-cell")) %>%
  dplyr::count(figshare_patient_id, figshare_sample_id, timepoint_2, predicted.celltype_3, name = "n_cells") %>%
  dplyr::group_by(figshare_sample_id) %>%
  dplyr::mutate(total_nonCAR_T_cells = sum(n_cells), T_subset_fraction = n_cells / total_nonCAR_T_cells) %>%
  dplyr::ungroup()

t_subset_stats <- t_subset_counts %>%
  dplyr::group_by(predicted.celltype_3) %>%
  dplyr::group_modify(~safe_wilcox_2group(.x, "T_subset_fraction", "timepoint_2", "Untreated Control Marrow", "CAR Marrow")) %>%
  dplyr::ungroup()

write_figshare_csv(t_subset_counts, "figshare_scRNA_Fig2_nonCAR_Tcell_subset_fraction_per_sample.csv", "Per-sample non-CAR marrow T-cell subset fractions.")
write_figshare_csv(t_subset_stats, "figshare_scRNA_Fig2_nonCAR_Tcell_subset_wilcoxon_statistics.csv", "Wilcoxon statistics for non-CAR marrow T-cell subset fractions.")

# =========================================================
# EXTERNAL TSVs: TCR DIVERSITY + FLOW IMMATURE B-CELL FRACTIONS
# =========================================================
if (file.exists(study_path("input_tcr_diversity", optional = TRUE))) {
  TCR_div <- read.delim(study_path("input_tcr_diversity", optional = TRUE), stringsAsFactors = FALSE)
  TCR_div_shannon <- TCR_div %>% dplyr::filter(Measure == "Shannon")
  tcr_comparisons <- list(c("Auto PBMC", "CAR PBMC"), c("CAR Marrow", "Control Marrow"))
  tcr_stats <- purrr::map_dfr(tcr_comparisons, function(comp) {
    dat <- TCR_div_shannon %>% dplyr::filter(Timepoint %in% comp)
    x <- dat %>% dplyr::filter(Timepoint == comp[[1]]) %>% dplyr::pull(Value)
    y <- dat %>% dplyr::filter(Timepoint == comp[[2]]) %>% dplyr::pull(Value)
    p <- tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    tibble::tibble(group1 = comp[[1]], group2 = comp[[2]], n_group1 = length(x), n_group2 = length(y), median_group1 = median(x, na.rm = TRUE), median_group2 = median(y, na.rm = TRUE), p_value = p, p_label = fmt_p(p), test = "Wilcoxon rank-sum")
  })
  write_figshare_csv(TCR_div_shannon, "figshare_scRNA_Fig2_TCR_diversity_shannon_plot_data.csv", "TCR Shannon diversity values used for supplementary TCR diversity plot.")
  write_figshare_csv(tcr_stats, "figshare_scRNA_Fig2_TCR_diversity_shannon_wilcoxon_statistics.csv", "Wilcoxon statistics for TCR Shannon diversity comparisons.")
}

if (file.exists(study_path("input_b_cell_proportions", optional = TRUE))) {
  b_cell_prop <- read.delim(study_path("input_b_cell_proportions", optional = TRUE), stringsAsFactors = FALSE)
  b_cell_prop <- b_cell_prop %>%
    dplyr::mutate(Type_label = dplyr::recode(Type, "CAR Marrow" = "CAR Marrow", "Control Marrow" = "Untreated Control Marrow"))
  p <- tryCatch({
    x <- b_cell_prop %>% dplyr::filter(Type_label == "CAR Marrow") %>% dplyr::pull(Immature)
    y <- b_cell_prop %>% dplyr::filter(Type_label == "Untreated Control Marrow") %>% dplyr::pull(Immature)
    stats::wilcox.test(x, y, exact = FALSE)$p.value
  }, error = function(e) NA_real_)
  b_cell_prop_stats <- tibble::tibble(group1 = "CAR Marrow", group2 = "Untreated Control Marrow", p_value = p, p_label = fmt_p(p), test = "Wilcoxon rank-sum")
  write_figshare_csv(b_cell_prop, "figshare_scRNA_Fig2_flow_immature_Bcell_fraction_plot_data.csv", "External flow-derived immature B-cell fraction plot data.")
  write_figshare_csv(b_cell_prop_stats, "figshare_scRNA_Fig2_flow_immature_Bcell_fraction_wilcoxon_statistics.csv", "Wilcoxon statistics for external flow immature B-cell fraction plot.")
}

# =========================================================
# MANIFEST + README
# =========================================================
readr::write_csv(manifest, figshare_file("figshare_scRNA_Fig2_dataset_manifest.csv"), na = "")

readme <- study_setting("readme")
writeLines(readme, con = figshare_file("figshare_scRNA_Fig2_README.txt"))
writeLines(capture.output(sessionInfo()), con = out_file("scRNA_Fig2_figshare_hybrid_sessionInfo.txt"))

cat("\nSaved scRNA Fig 2 Figshare datasets to:\n", figshare_dir, "\n", sep = "")
cat("Saved hybrid vector/rastered-point UMAP outputs to:\n", hybrid_umap_dir, "\n", sep = "")
