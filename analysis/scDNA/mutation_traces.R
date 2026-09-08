#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("mutation_traces")

# Hybrid vector/raster t-SNE export version generated for Figshare/AI editing
#!/usr/bin/env Rscript

setwd(study_path("study_path_001"))

## Final output directory for paper revision figures/tables/logs
out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

############################################################
## HYBRID VECTOR t-SNE/UMAP EXPORTS WITH RASTERED POINTS
## Creates Illustrator-friendly PDF/EPS files: axes/text/legends remain vector,
## while dense point layers are rasterized to keep file size manageable.
############################################################
tsne_vector_raster_dir <- out_file("tsne_vector_raster_points")
dir.create(tsne_vector_raster_dir, showWarnings = FALSE, recursive = TRUE)
tsne_vector_raster_file <- function(...) file.path(tsne_vector_raster_dir, ...)

## Use ggrastr only for the duplicate hybrid exports. Original outputs are unchanged.
if (!requireNamespace("ggrastr", quietly = TRUE)) {
  stop("Package 'ggrastr' is required for hybrid vector/raster t-SNE exports. Install with install.packages('ggrastr') and rerun.")
}

## Filenames matching these patterns will be duplicated into tsne_vector_raster_points/.
## This catches the saved t-SNE/UMAP embedding panels and marker-overlay embedding panels.
tsne_vector_raster_patterns <- study_setting("tsne_vector_raster_patterns")

should_save_tsne_hybrid <- function(filename) {
  fn <- basename(as.character(filename))
  ext <- tolower(tools::file_ext(fn))
  if (!ext %in% c("pdf", "eps", "svg")) return(FALSE)
  any(stringr::str_detect(fn, stringr::regex(paste(tsne_vector_raster_patterns, collapse = "|"), ignore_case = TRUE)))
}

geom_point_ai_rast <- function(..., raster.dpi = 600) {
  ggrastr::geom_point_rast(..., raster.dpi = raster.dpi)
}

rasterize_point_layers_for_ai <- function(plot, dpi = 600) {
  tryCatch(
    ggrastr::rasterise(plot, layers = "Point", dpi = dpi),
    error = function(e) {
      message("[WARN] Could not rasterize point layers with ggrastr; saving original plot for hybrid export: ", conditionMessage(e))
      plot
    }
  )
}

## Wrapper around ggplot2::ggsave. It still writes the original file, then writes a
## duplicate hybrid PDF/EPS/SVG in tsne_vector_raster_points/ when the filename is a t-SNE/UMAP panel.
ggsave <- function(filename, plot = ggplot2::last_plot(), ..., device = NULL) {
  out <- ggplot2::ggsave(filename = filename, plot = plot, ..., device = device)

  if (should_save_tsne_hybrid(filename)) {
    fn <- basename(as.character(filename))
    ext <- tolower(tools::file_ext(fn))
    stem <- tools::file_path_sans_ext(fn)
    hybrid_filename <- tsne_vector_raster_file(paste0(stem, "_hybrid_raster_points.", ext))
    hybrid_plot <- rasterize_point_layers_for_ai(plot, dpi = 600)
    ggplot2::ggsave(filename = hybrid_filename, plot = hybrid_plot, ..., device = device)
    message("[HYBRID t-SNE EXPORT] ", hybrid_filename)
  }

  invisible(out)
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(broom)
  library(ggh4x)
  library(grid)
  library(patchwork)
  library(stringr)
})

############################################################
## Shared helpers, palettes, themes
############################################################

## CAR colors – match original script
CAR_GREEN <- "darkgreen"
pal_car   <- c("Absent" = "grey80", "Present" = CAR_GREEN)

## collapse cell types for scDNA
collapse_celltype_scDNA <- function(x) {
  case_when(
    x %in% c("CD4_Tcell", "CD4_T-cell", "CD4 T-cell") ~ "CD4 T-cell",
    x %in% c("CD8_Tcell", "CD8_T-cell", "CD8 T-cell") ~ "CD8 T-cell",
    x %in% c("NK_cell", "NK_Cell", "NK cell")         ~ "NK cell",

    ## Treat all of these as myeloid for clonal tracing
    x %in% c(
      "Myeloid",
      "Myeloid_NOS",
      "CD34_progenitor",   # now myeloid
      "CD71_erythroid"     # optional but usually myeloid-lineage
    ) ~ "Myeloid",

    ## If you still want “Stem Cell” separated:
    x %in% c("Stem Cell") ~ "Myeloid Progenitor",

    TRUE ~ "Other/Unknown"
  )
}

## publication theme (no explicit font to avoid PDF font errors)
theme_pub <- theme_bw(base_size = 16) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(hjust = 0)
  )

pal_variant <- c(
  "DNMT3A" = "#0072B2",  # Okabe-Ito blue: good / less toxic
  "TET2"   = "#009E73",  # bluish green
  "TP53"   = "#C43A31",  # red-vermilion/crimson: bad / more toxic
  "PPM1D"  = "#CC79A7",  # reddish purple
  "BCOR"   = "#E69F00",  # orange
  "ATM"    = "#56B4E9",  # sky blue
  "EZH2"   = "#000000",  # black
  "IDH1"   = "#F0E442"   # yellow
)

## Robustly recover the CH gene from full Variant_Call/Base_Variant strings.
## This prevents PPM1D or other mutations from falling to ggplot's NA/grey
## when the plotted label is a full variant string rather than exactly "PPM1D".
gene_from_variant <- function(x) {
  x <- as.character(x)
  gene_regex <- paste(names(pal_variant), collapse = "|")
  stringr::str_extract(stringr::str_to_upper(x), gene_regex)
}

## CH palette for t-SNE (adds Other + WT grey)
pal_CH_category <- c(
  pal_variant[c("DNMT3A","TET2","PPM1D","TP53","BCOR","ATM","EZH2")],
  "Other" = "grey50",
  "WT"    = "grey85"
)

## For clonal tracing plots where the color aesthetic is Base_Variant
## rather than Gene, force each variant label to inherit its gene color.
variant_gene_palette <- function(labels, genes = NULL, default = "grey50") {
  labels <- as.character(labels)
  genes <- if (is.null(genes)) gene_from_variant(labels) else gene_from_variant(genes)

  keep <- !is.na(labels) & labels != ""
  labels <- labels[keep]
  genes   <- genes[keep]

  dedup <- !duplicated(labels)
  labels <- labels[dedup]
  genes   <- genes[dedup]

  cols  <- unname(pal_variant[genes])
  cols[is.na(cols)] <- default

  stats::setNames(cols, labels)
}

## Lineage palette refinements requested for color-blind contrast.
## NK is a saturated purple to separate it from gray Other/Unknown cells.
## Myeloid progenitor-like cells are dark teal to separate from orange myeloid cells.
NK_CELL_COLOR <- "#AA4499"
MYELOID_PROGENITOR_LIKE_COLOR <- "#00A6D6"

pal_celltype <- c(
  "CD4 T-cell"    = "#C6CDF7",
  "CD8 T-cell"    = "#7294D4",
  "NK cell"       = NK_CELL_COLOR,
  "Myeloid"       = "#F98400",
  "Myeloid Progenitor"      = MYELOID_PROGENITOR_LIKE_COLOR,
  "Myeloid Progenitor-Like" = MYELOID_PROGENITOR_LIKE_COLOR,
  "Unknown"       = "grey70",
  "Other/Unknown" = "grey70"
)

## Backward-compatible alias used by later lineage-composition plots
pal_lineage5 <- pal_celltype

## CCT ID -> Paper ID key
id_key <- study_setting("id_key")

## CAR palette for all new plots
pal_CAR <- c("Present" = CAR_GREEN, "Absent" = "lightgray")

axis_trunc <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

theme_tsne <- function() {
  theme_classic(base_size = 20) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.line.y  = element_line(
        arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")
      ),
      axis.line.x  = element_line(
        arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")
      ),
      axis.title      = element_text(hjust = 0),
      legend.position = "right"
    )
}

############################################################
## PART A – All cells (df_all): t-SNEs, CAR enrichment, CD4/CD8 stats
############################################################
input_file_all <- study_path("input_file_all")

df_all <- readr::read_csv(input_file_all, guess_max = 100000, show_col_types = FALSE)

required_cols_all <- c("Patient_ID", "OmiqFilter", "CAR",
                       "Variant_Call", "optsne_1", "optsne_2")
missing_cols_all <- setdiff(required_cols_all, names(df_all))
if (length(missing_cols_all) > 0) {
  stop(
    "Missing expected columns in df_all input: ",
    paste(missing_cols_all, collapse = ", ")
  )
}

df_all <- df_all %>%
  mutate(
    ## cell type
    Cell_Type5 = case_when(
      OmiqFilter %in% c("CD4_Tcell", "CD4_T") ~ "CD4 T-cell",
      OmiqFilter %in% c("CD8_Tcell", "CD8_T") ~ "CD8 T-cell",
      OmiqFilter %in% c("NK Cell","NK cell","NK_Cell","NK_cell","NKcell","NK")
        ~ "NK cell",
      OmiqFilter %in% c(
        "Myeloid","Myeloid_NOS","Myeloid-Monocyte-Neutrophil-Macrophage",
        "CD71+Erythroid","CD71_erythroid"
      ) ~ "Myeloid",
      OmiqFilter %in% c("Myeloid Progenitor","CD34+Progenitor","CD34_progenitor")
        ~ "Myeloid Progenitor",
      TRUE ~ "Unknown"
    ),
    Cell_Type5 = factor(
      Cell_Type5,
      levels = c("CD4 T-cell","CD8 T-cell","NK cell",
                 "Myeloid","Myeloid Progenitor","Unknown")
    ),

    ## CAR detection from df_all$CAR
    CAR_detected = if_else(CAR == "CAR", "Present", "Absent"),
    CAR_detected = factor(CAR_detected, levels = c("Present","Absent")),
    is_CAR       = CAR_detected == "Present",

    ## patient + timepoint (CTN ID here)
    patient_core     = sub("_[^_]+$", "", Patient_ID),
    sample_timepoint = sub("^[^_]+_", "", Patient_ID),

    ## T-cell flags
    is_Tcell = Cell_Type5 %in% c("CD4 T-cell","CD8 T-cell"),
    is_CD4   = Cell_Type5 == "CD4 T-cell",
    is_CD8   = Cell_Type5 == "CD8 T-cell",

    ## CH annotations
    is_CH   = !is.na(Variant_Call) & Variant_Call != "WT",
    CH_gene = if_else(is_CH, sub(":.*$", "", Variant_Call), NA_character_),

    ## collapse TP53 flavours
    CH_gene_simple = case_when(
      grepl("^TP53", CH_gene) ~ "TP53",
      TRUE                    ~ CH_gene
    ),

    CH_category = case_when(
      !is_CH ~ "WT",
      CH_gene_simple %in% c("DNMT3A","TET2","PPM1D",
                            "TP53","BCOR","ATM","EZH2") ~ CH_gene_simple,
      TRUE ~ "Other"
    ),

    is_CAR_CH        = is_CAR & is_CH,
    CH_gene_for_plot = if_else(
      is_CAR_CH & !is.na(CH_gene_simple),
      CH_gene_simple,
      NA_character_
    )
  ) %>%
  ## bring in Paper_ID and create patient_paper used for all plots/statistics
  left_join(id_key, by = c("patient_core" = "CTN_ID")) %>%
  mutate(
    patient_paper = dplyr::coalesce(Paper_ID, patient_core)
  )

cat("\n[PART A] Columns after annotation (df_all):\n")
print(colnames(df_all))
## 1) t-SNE: lineage
p_tsne_celltype <- ggplot(
  df_all %>% arrange(Cell_Type5),
  aes(x = optsne_1, y = optsne_2, color = Cell_Type5)
) +
  geom_point_ai_rast(size = 0.3, alpha = 0.8) +
  scale_colour_manual(
    values = pal_celltype,
    na.value = "grey80",
    name = "Cell Lineage"
  ) +
  xlab("t-SNE1") + ylab("t-SNE2") +
  ggtitle("Lineage Distribution of PBL Cells") +
  theme_tsne() +
  guides(
    x = axis_trunc,
    y = axis_trunc,
    color = guide_legend(override.aes = list(size = 4, alpha = 1))
  )

## 2) t-SNE: CH mutations (no CH in grey)
p_tsne_CH <- ggplot() +
  geom_point_ai_rast(
    data = df_all %>% filter(CH_category == "WT"),
    aes(x = optsne_1, y = optsne_2),
    color = "grey80",
    size = 0.3, alpha = 0.4
  ) +
  geom_point_ai_rast(
    data = df_all %>% filter(CH_category != "WT"),
    aes(x = optsne_1, y = optsne_2, color = CH_category),
    size = 0.3, alpha = 0.9
  ) +
  scale_colour_manual(
    values = pal_CH_category,
    na.value = "grey80",
    name = "CH Mutation"
  ) +
  xlab("t-SNE1") + ylab("t-SNE2") +
  ggtitle("All CH Mutant Cells") +
  theme_tsne() +
  guides(
    x = axis_trunc,
    y = axis_trunc,
    color = guide_legend(override.aes = list(size = 4, alpha = 1))
  )

## 3) t-SNE: CAR detection (CAR– grey background, CAR+ darkgreen on top)
p_tsne_CAR <- ggplot() +
  geom_point_ai_rast(
    data = df_all %>% filter(!is_CAR),
    aes(x = optsne_1, y = optsne_2, color = "Absent"),
    size = 0.3, alpha = 0.4
  ) +
  geom_point_ai_rast(
    data = df_all %>% filter(is_CAR),
    aes(x = optsne_1, y = optsne_2, color = "Present"),
    size = 0.3, alpha = 0.9
  ) +
  scale_colour_manual(
    values = pal_CAR,
    breaks = c("Present","Absent"),
    name   = "CAR Detected"
  ) +
  xlab("t-SNE1") + ylab("t-SNE2") +
  ggtitle("All CAR+ Cells") +
  theme_tsne() +
  guides(
    x = axis_trunc,
    y = axis_trunc,
    color = guide_legend(override.aes = list(size = 4, alpha = 1))
  )

## 4) t-SNE: CAR+ CH+ cells highlighted
## - make CAR+CH+ points larger + on top
## - italicize gene names (kept) BUT remove legend entirely (redundant)

gene_labeller_expr <- function(x) {
  lapply(x, function(g) parse(text = paste0("italic(", g, ")")))
}

car_ch_genes <- intersect(
  unique(na.omit(df_all$CH_gene_for_plot)),
  names(pal_CH_category)[names(pal_CH_category) != "WT"]
)

p_tsne_CAR_CH <- ggplot() +
  geom_point_ai_rast(
    data = df_all,
    aes(x = optsne_1, y = optsne_2),
    color = "grey80",
    size = 0.25, alpha = 0.35
  ) +
  geom_point_ai_rast(
    data = df_all %>% filter(is_CAR_CH, !is.na(CH_gene_for_plot)),
    aes(x = optsne_1, y = optsne_2, color = CH_gene_for_plot),
    size = 1.0, alpha = 0.98
  ) +
  scale_colour_manual(
    values = pal_CH_category[names(pal_CH_category) != "WT"],
    breaks = car_ch_genes,
    labels = gene_labeller_expr(car_ch_genes),
    name   = "CH Mutation"   # label here doesn't matter since legend is removed
  ) +
  xlab("t-SNE1") + ylab("t-SNE2") +
  ggtitle("All CAR+ CH Mutant Cells") +
  theme_tsne() +
  guides(
    x = axis_trunc,
    y = axis_trunc,
    color = "none"  # REMOVE legend for panel 4
  )

## Row of 4 t-SNEs
p_tsne_row <- p_tsne_celltype + p_tsne_CH + p_tsne_CAR + p_tsne_CAR_CH +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "right")

print(p_tsne_row)

## Save: PDF + EPS + PNG
out_base <- out_file("row_tsne_pbl_expansion")

ggsave(
  paste0(out_base, ".pdf"),
  plot  = p_tsne_row,
  width = 26, height = 8, units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(out_base, ".eps"),
  plot  = p_tsne_row,
  width = 26, height = 8, units = "in",
  device = cairo_ps
)

# PNG (raster) - use DPI; do NOT use cairo_ps here
ggsave(
  paste0(out_base, ".png"),
  plot  = p_tsne_row,
  width = 26, height = 8, units = "in",
  dpi   = 400
)


############################################################
## COUNTS to annotate the t-SNE panels (df_all)
## - total cells
## - CH mutant cells
## - CAR+ cells
## - CAR+ CH mutant cells
############################################################

stopifnot(all(c("is_CAR","is_CH","is_CAR_CH","CH_category") %in% colnames(df_all)))

# Core counts
n_total         <- nrow(df_all)
n_CH_mut        <- sum(df_all$is_CH, na.rm = TRUE)                          # Variant_Call != WT (per your is_CH)
n_CAR_pos       <- sum(df_all$is_CAR, na.rm = TRUE)                         # CAR_detected == Present
n_CAR_CH_mut    <- sum(df_all$is_CAR_CH & df_all$is_CH, na.rm = TRUE)       # CAR+ AND CH mutant

cat("\n=== [t-SNE COUNTS: df_all] ===\n")
cat("Total cells:                 ", n_total, "\n", sep = "")
cat("CH mutant cells:             ", n_CH_mut, "\n", sep = "")
cat("CAR+ cells:                  ", n_CAR_pos, "\n", sep = "")
cat("CAR+ CH mutant cells:        ", n_CAR_CH_mut, "\n", sep = "")

# Optional sanity checks (helpful if anything seems off)
cat("\n[CHECK] is_CAR_CH vs (is_CAR & is_CH):\n")
print(table(df_all$is_CAR_CH, (df_all$is_CAR & df_all$is_CH), useNA = "ifany"))

cat("\n[CHECK] CH_category breakdown:\n")
print(sort(table(df_all$CH_category, useNA = "ifany"), decreasing = TRUE))

# Ready-to-use panel subtitles (paste into labs(subtitle=...) or ggtitle/subtitle)
sub_lineage <- paste0("n=", n_total,
                      " | CH mutant=", n_CH_mut,
                      " | CAR+=", n_CAR_pos,
                      " | CAR+CH mutant=", n_CAR_CH_mut)

sub_CH      <- paste0("CH mutant n=", n_CH_mut, " (of ", n_total, ")")
sub_CAR     <- paste0("CAR+ n=", n_CAR_pos, " (of ", n_total, ")")
sub_CAR_CH  <- paste0("CAR+CH mutant n=", n_CAR_CH_mut, " (of ", n_total, ")")

cat("\nSuggested subtitles:\n")
cat("Lineage panel: ", sub_lineage, "\n", sep = "")
cat("CH panel:      ", sub_CH, "\n", sep = "")
cat("CAR panel:     ", sub_CAR, "\n", sep = "")
cat("CAR+CH panel:  ", sub_CAR_CH, "\n", sep = "")

# If you want, store these in a named list for reuse
tsne_counts <- list(
  n_total      = n_total,
  n_CH_mut     = n_CH_mut,
  n_CAR_pos    = n_CAR_pos,
  n_CAR_CH_mut = n_CAR_CH_mut,
  subtitle = list(
    lineage = sub_lineage,
    CH      = sub_CH,
    CAR     = sub_CAR,
    CAR_CH  = sub_CAR_CH
  )
)


## Peak timepoint lookup based on CAR expansion (by Paper_ID)
peak_meta <- df_all %>%
  group_by(patient_paper, sample_timepoint) %>%
  summarise(
    n_cells      = n(),
    frac_CAR_all = mean(is_CAR),
    n_T          = sum(is_Tcell),
    frac_CAR_T   = if_else(
      n_T > 0,
      sum(is_CAR & is_Tcell) / n_T,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    peak_metric = if_else(
      !is.na(frac_CAR_T) & frac_CAR_T > 0,
      frac_CAR_T,
      frac_CAR_all
    )
  )

peak_lookup <- peak_meta %>%
  group_by(patient_paper) %>%
  slice_max(order_by = peak_metric, n = 1, with_ties = FALSE) %>%
  ungroup()

df_peak <- df_all %>%
  inner_join(
    peak_lookup %>% select(patient_paper, sample_timepoint),
    by = c("patient_paper","sample_timepoint")
  )

cat("\n=== [PART A] Peak timepoint per patient (by Paper_ID) ===\n")
print(peak_lookup, n = Inf)

## CD4 vs CD8 CAR+ T cells at peak (per-patient, by Paper_ID)
peak_car_cd4_cd8 <- df_peak %>%
  filter(is_CAR, is_Tcell, Cell_Type5 %in% c("CD4 T-cell","CD8 T-cell")) %>%
  count(patient_paper, Cell_Type5, name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from  = Cell_Type5,
    values_from = n_cells,
    values_fill = 0
  )

cat("\n=== [PART A] PEAK: per-patient CAR+ CD4 vs CD8 T-cell counts ===\n")
print(peak_car_cd4_cd8, n = Inf)

if (nrow(peak_car_cd4_cd8) > 0) {
  w_cd4_cd8_peak <- with(
    peak_car_cd4_cd8,
    wilcox.test(`CD4 T-cell`, `CD8 T-cell`, paired = TRUE)
  )

  cat("\nWilcoxon (paired) CD4 vs CD8 CAR+ T cells at PEAK (per-patient counts):\n")
  print(w_cd4_cd8_peak)

  total_cd4 <- sum(peak_car_cd4_cd8$`CD4 T-cell`)
  total_cd8 <- sum(peak_car_cd4_cd8$`CD8 T-cell`)
  b_cd4_cd8_peak <- binom.test(total_cd4, total_cd4 + total_cd8, p = 0.5)

  cat("\nBinomial test (pooled) CD4 vs CD8 CAR+ T cells at PEAK:\n")
  print(b_cd4_cd8_peak)
}

## CD4 vs CD8 CAR+ CH+ T cells at peak (Paper_ID)
peak_car_ch_cd4_cd8 <- df_peak %>%
  filter(is_CAR, is_CH, is_Tcell,
         Cell_Type5 %in% c("CD4 T-cell","CD8 T-cell")) %>%
  count(patient_paper, Cell_Type5, name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from  = Cell_Type5,
    values_from = n_cells,
    values_fill = 0
  )

cat("\n=== [PART A] PEAK: per-patient CAR+ CH+ CD4 vs CD8 T-cell counts ===\n")
print(peak_car_ch_cd4_cd8, n = Inf)

if (nrow(peak_car_ch_cd4_cd8) > 0) {
  w_cd4_cd8_peak_ch <- with(
    peak_car_ch_cd4_cd8,
    wilcox.test(`CD4 T-cell`, `CD8 T-cell`, paired = TRUE)
  )

  cat("\nWilcoxon (paired) CD4 vs CD8 CAR+ CH+ T cells at PEAK (per-patient counts):\n")
  print(w_cd4_cd8_peak_ch)

  total_cd4_ch <- sum(peak_car_ch_cd4_cd8$`CD4 T-cell`)
  total_cd8_ch <- sum(peak_car_ch_cd4_cd8$`CD8 T-cell`)
  b_cd4_cd8_peak_ch <- binom.test(total_cd4_ch, total_cd4_ch + total_cd8_ch, p = 0.5)

  cat("\nBinomial test (pooled) CD4 vs CD8 CAR+ CH+ T cells at PEAK:\n")
  print(b_cd4_cd8_peak_ch)
}

## Enrichment / depletion of CH genes in CAR+ cells (all CH+ cells)
df_ch_enrich <- df_all %>%
  filter(is_CH, !is.na(CH_gene_simple))

gene_counts <- df_ch_enrich %>%
  count(CH_gene_simple, sort = TRUE)

cat("\n=== [PART A] CH+ cells per gene (collapsed TP53) ===\n")
print(gene_counts, n = Inf)

genes_to_test <- gene_counts %>%
  filter(n >= 20) %>%
  pull(CH_gene_simple)

enrich_list <- lapply(genes_to_test, function(g) {
  a <- sum(df_ch_enrich$CH_gene_simple == g & df_ch_enrich$is_CAR)      # gene g, CAR+
  b <- sum(df_ch_enrich$CH_gene_simple == g & !df_ch_enrich$is_CAR)     # gene g, CAR-
  c <- sum(df_ch_enrich$CH_gene_simple != g & df_ch_enrich$is_CAR)      # others, CAR+
  d <- sum(df_ch_enrich$CH_gene_simple != g & !df_ch_enrich$is_CAR)     # others, CAR-

  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      gene = c(g, "Other"),
      CAR  = c("CAR+","CAR-")
    )
  )

  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(tibble(
      CH_gene      = g,
      log2_OR      = NA_real_,
      log2_OR_low  = NA_real_,
      log2_OR_high = NA_real_,
      p_value      = NA_real_
    ))
  }

  ft <- fisher.test(mat)

  tibble(
    CH_gene      = g,
    log2_OR      = log2(unname(ft$estimate)),
    log2_OR_low  = log2(ft$conf.int[1]),
    log2_OR_high = log2(ft$conf.int[2]),
    p_value      = ft$p.value
  )
})

enrich_df <- bind_rows(enrich_list) %>%
  arrange(desc(log2_OR))

cat("\n=== [PART A] Enrichment of CH genes in CAR+ vs CAR- (log2 OR) ===\n")
print(enrich_df, n = Inf)

enrich_df_plot <- enrich_df %>% filter(!is.na(log2_OR))

p_enrich <- ggplot(
  enrich_df_plot,
  aes(
    x = log2_OR,
    y = forcats::fct_reorder(CH_gene, log2_OR, .na_rm = TRUE),
    color = CH_gene
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = log2_OR_low, xmax = log2_OR_high),
    height = 0.2
  ) +
  geom_point(size = 3) +
  scale_color_manual(values = pal_CH_category, guide = "none") +
  xlab("log2 odds ratio for CAR+ vs CAR- (CH+ cells)") +
  ylab("CH gene") +
  ggtitle("Enrichment / depletion of CH in CAR+ cells") +
  theme_classic(base_size = 18)

print(p_enrich)

ggsave(
  out_file("enrichmentplot.eps"),
  plot  = p_enrich,
  width = 8, height = 5, units = "in",
  device = cairo_ps
)

############################################################
## PART B – Clonal tracing of myeloid CH (df3)
############################################################

## 1) Load & annotate ----------------------------------------

df3 <- readr::read_delim(
  study_path("input_scdna_pbmc_bm_merged"),
  show_col_types = FALSE
)

df3 <- df3 %>%
  mutate(
    Cell_Type5 = collapse_celltype_scDNA(OmiqFilter),
    Cell_Type5 = factor(
      Cell_Type5,
      levels = c("CD4 T-cell","CD8 T-cell","NK cell",
                 "Myeloid","Myeloid Progenitor","Other/Unknown")
    ),
    is_myeloid = Cell_Type5 == "Myeloid",
    is_CH      = Variant_Call != "WT",

    ## identify point mutations
    is_point_mut = str_detect(Variant_Call, "^[^:]+:p\\."),

    ## founder (Base_Variant): first "GENE:p.xxx" if present
    Base_Variant = if_else(
      is_point_mut,
      str_extract(Variant_Call, "^[^:]+:p\\.[^:_]+"),
      Variant_Call
    ),

    ## gene from founder
    Gene = gene_from_variant(Base_Variant),

    ## clone label for subclone plot:
    ##  - founder-only cells: just Base_Variant
    ##  - non-founder point-mutation cells: full Variant_Call
    ##  - non point-mutations: full Variant_Call
    Clone_Label = case_when(
      is_point_mut & Variant_Call == Base_Variant ~ Base_Variant,
      is_point_mut                                ~ Variant_Call,
      TRUE                                        ~ Variant_Call
    )
  ) %>%
  ## add Paper_ID and plotting ID
  left_join(id_key, by = c("Patient" = "CTN_ID")) %>%
  mutate(
    Patient_Paper = dplyr::coalesce(Paper_ID, Patient)
  )

cat("\n=== Collapsed cell-type counts (df3) ===\n")
df3 %>% count(Cell_Type5) %>% print(n = Inf)

## 2) Restrict to patients with >1 timepoint -----------------

pt_time_summary <- df3 %>%
  distinct(Patient_Paper, Day) %>%
  count(Patient_Paper, name = "n_timepoints")

multi_pt <- pt_time_summary %>% filter(n_timepoints > 1)

cat("\n=== Patients (Paper_ID) with >1 timepoint in df3 ===\n")
multi_pt %>% arrange(Patient_Paper) %>% print(n = Inf)

df3_multi <- df3 %>% semi_join(multi_pt, by = "Patient_Paper")

all_days   <- sort(unique(df3_multi$Day))
day_labels <- paste0("D", all_days)

## 3) Myeloid totals & CH+ fraction over time -----------------

myeloid_totals <- df3_multi %>%
  filter(is_myeloid) %>%
  group_by(Patient_Paper, Day) %>%
  summarise(
    n_myeloid       = n(),
    n_myeloid_CH    = sum(is_CH),
    frac_CH_myeloid = n_myeloid_CH / pmax(n_myeloid, 1L),
    .groups         = "drop"
  )

cat("\n=== Myeloid totals & CH+ fractions (multi-timepoint) ===\n")
myeloid_totals %>%
  arrange(Patient_Paper, Day) %>%
  print(n = Inf)

## Plot: CH+ fraction in myeloid over time (per patient; Paper_ID facets)
p_myeloid_CH_time <- ggplot(
  myeloid_totals,
  aes(x = Day, y = frac_CH_myeloid, group = Patient_Paper)
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = all_days, labels = day_labels) +
  labs(
    title = "Fraction of myeloid cells that are CH+ over time",
    x = "Timepoint (Day)",
    y = "Fraction of myeloid cells that are CH+"
  ) +
  theme_pub

## 4) Founder-level myeloid clone trajectories ---------------

## observed founder clone counts in myeloid
myeloid_clone_obs <- df3_multi %>%
  filter(is_myeloid, is_CH) %>%
  group_by(Patient_Paper, Day, Base_Variant, Gene) %>%
  summarise(
    n_cells_clone = n(),
    .groups       = "drop"
  )

## clone list per patient
clone_meta <- myeloid_clone_obs %>%
  distinct(Patient_Paper, Base_Variant, Gene)

## full grid of (Patient_Paper, Day, Base_Variant, Gene) with myeloid totals
myeloid_clone_time <- myeloid_totals %>%
  select(Patient_Paper, Day, n_myeloid) %>%
  left_join(
    clone_meta,
    by = "Patient_Paper",
    relationship = "many-to-many"
  ) %>%
  left_join(
    myeloid_clone_obs,
    by = c("Patient_Paper","Day","Base_Variant","Gene")
  ) %>%
  mutate(
    n_cells_clone   = coalesce(n_cells_clone, 0L),
    frac_of_myeloid = n_cells_clone / pmax(n_myeloid, 1L)
  )

cat("\n=== Founder-level myeloid clones over time (first 50 rows) ===\n")
myeloid_clone_time %>%
  arrange(Patient_Paper, Base_Variant, Day) %>%
  print(n = 50)

## keep founders that ever reach >= 1% of myeloid
clone_keep <- myeloid_clone_time %>%
  group_by(Patient_Paper, Base_Variant, Gene) %>%
  summarise(
    max_frac = max(frac_of_myeloid, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  filter(max_frac >= 0.01)

cat("\n=== Founder clones (Paper_ID) reaching ≥1% of myeloid ===\n")
clone_keep %>%
  arrange(Patient_Paper, desc(max_frac)) %>%
  print(n = Inf)

myeloid_clone_time_filt <- myeloid_clone_time %>%
  inner_join(
    clone_keep %>% select(Patient_Paper, Base_Variant),
    by = c("Patient_Paper","Base_Variant")
  )

## drop [study-specific] ([study-specific]) from this figure
myeloid_clone_time_filt <- myeloid_clone_time_filt %>%
  filter(Patient_Paper != study_setting("study_value_001"))

cat("\n[COLOR CHECK] Genes in main myeloid clone plot not in pal_variant:\n")
print(setdiff(unique(na.omit(myeloid_clone_time_filt$Gene)), names(pal_variant)))

## Plot: founder clone trajectories, colored by CHIP gene palette
p_myeloid_clone_traj <- ggplot(
  myeloid_clone_time_filt,
  aes(
    x     = Day,
    y     = frac_of_myeloid,
    group = Base_Variant,
    color = Gene
  )
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(
    labels = function(x) paste0("D", x)
  ) +
  scale_color_manual(
    values = pal_variant,
    breaks = c("DNMT3A","TET2","TP53","PPM1D","BCOR","ATM","EZH2","IDH1"),
    name   = "Gene",
    na.value = "grey50"
  ) +
  labs(
    title = "CH gene mutations in myeloid cells over time measured by scDNA",
    x = "Timepoint (Day)",
    y = "Fraction of myeloid compartment"
  ) +
  theme_pub +
  theme(legend.position = "right")

ggsave(
  out_file("myeloid_cells_with_CH_together.eps"),
  plot  = p_myeloid_clone_traj,
  width = 10, height = 5, units = "in",
  device = cairo_ps
)

## 5) Subclone trajectories – point mutations only -----------

myeloid_subclone_obs <- df3_multi %>%
  filter(is_myeloid, is_CH, is_point_mut) %>%
  group_by(Patient_Paper, Day, Clone_Label, Base_Variant, Gene) %>%
  summarise(
    n_cells_clone = n(),
    .groups       = "drop"
  )

subclone_meta <- myeloid_subclone_obs %>%
  distinct(Patient_Paper, Clone_Label, Base_Variant, Gene)

subclone_time <- myeloid_totals %>%
  select(Patient_Paper, Day, n_myeloid) %>%
  left_join(
    subclone_meta,
    by = "Patient_Paper",
    relationship = "many-to-many"
  ) %>%
  left_join(
    myeloid_subclone_obs,
    by = c("Patient_Paper","Day","Clone_Label","Base_Variant","Gene")
  ) %>%
  mutate(
    n_cells_clone   = coalesce(n_cells_clone, 0L),
    frac_of_myeloid = n_cells_clone / pmax(n_myeloid, 1L),
    founder_status  = if_else(Clone_Label == Base_Variant,
                              "Founder", "Non-founder")
  )

## keep subclones that ever reach >=1% of myeloid
subclone_keep <- subclone_time %>%
  group_by(Patient_Paper, Clone_Label, Base_Variant, Gene) %>%
  summarise(
    max_frac = max(frac_of_myeloid, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  filter(max_frac >= 0.01)

cat("\n=== Subclones (Paper_ID) reaching ≥1% of myeloid (point mutations) ===\n")
subclone_keep %>%
  arrange(Patient_Paper, desc(max_frac)) %>%
  print(n = Inf)

subclone_time_filt <- subclone_time %>%
  inner_join(
    subclone_keep %>% select(Patient_Paper, Clone_Label),
    by = c("Patient_Paper","Clone_Label")
  )

p_myeloid_subclone_traj <- ggplot(
  subclone_time_filt,
  aes(
    x        = Day,
    y        = frac_of_myeloid,
    group    = Clone_Label,
    color    = Base_Variant,
    linetype = founder_status
  )
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.0) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = all_days, labels = day_labels) +
  scale_color_manual(
    values = variant_gene_palette(
      subclone_time_filt$Base_Variant,
      subclone_time_filt$Gene
    ),
    name = "Founder (Base_Variant)"
  ) +
  scale_linetype_manual(
    values = c("Founder" = "solid", "Non-founder" = "dashed"),
    name   = "Clone type"
  ) +
  labs(
    title = "Myeloid point-mutation founders vs subclones over time",
    x = "Timepoint (Day)",
    y = "Fraction of myeloid compartment"
  ) +
  theme_pub +
  theme(legend.position = "right")

## 6) CH+ myeloid fraction over time – single black line ------
## Patients with >1 timepoint (Paper_ID)
multi_pt_patients <- myeloid_totals %>%
  count(Patient_Paper) %>%
  filter(n > 1) %>%
  pull(Patient_Paper)

## Keep multi-timepoint patients, drop [study-specific], and build per-day labels
myeloid_totals_multi <- myeloid_totals %>%
  filter(Patient_Paper %in% multi_pt_patients,
         Patient_Paper != study_setting("study_value_001")) %>%
  mutate(
    Day_factor = factor(
      Day,
      levels = sort(unique(Day)),
      labels = paste0("D", sort(unique(Day)))
    )
  )

p_myeloid_state_lines <- ggplot(
  myeloid_totals_multi,
  aes(x = Day_factor, y = frac_CH_myeloid, group = 1)
) +
  geom_line(color = "black", size = 0.9) +
  geom_point(color = "black", size = 2.2) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Fraction of myeloid cells that are CH+ over time",
    x = "Timepoint (Day)",
    y = "Fraction of myeloid compartment that is CH+"
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_myeloid_state_lines

ggsave(
  out_file("myeloid_cells_with_CH.eps"),
  plot  = p_myeloid_state_lines,
  width = 8, height = 5, units = "in",
  device = cairo_ps
)

## 7) Per-patient founder clone plots (printed & saved individually) ---

patients_for_clones <- myeloid_clone_time_filt %>%
  group_by(Patient_Paper) %>%
  filter(Patient_Paper != study_setting("study_value_001")) %>%
  summarise(
    total_clone_cells = sum(n_cells_clone, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_clone_cells > 0) %>%
  pull(Patient_Paper) %>%
  sort()

cat("\n=== Making per-patient founder clone plots for (Paper_ID): ===\n")
print(patients_for_clones)

for (pt in patients_for_clones) {
  dat_pt <- myeloid_clone_time_filt %>%
    filter(Patient_Paper == pt)

  pt_days <- sort(unique(dat_pt$Day))

  p_pt <- ggplot(
    dat_pt,
    aes(
      x     = Day,
      y     = frac_of_myeloid,
      group = Base_Variant,
      color = Base_Variant
    )
  ) +
    geom_line(size = 0.9) +
    geom_point(size = 2.2) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(
      breaks = pt_days,
      labels = paste0("D", pt_days)
    ) +
    scale_color_manual(
      values = variant_gene_palette(dat_pt$Base_Variant, dat_pt$Gene),
      name   = "Founder mutation\n(Base_Variant)"
    ) +
    labs(
      title = paste0("Myeloid founder clones over time: ", pt),
      x     = "Timepoint (Day)",
      y     = "Fraction of myeloid compartment",
      color = "Founder mutation\n(Base_Variant)"
    ) +
    theme_pub +
    theme(
      legend.position = "right",
      plot.title      = element_text(hjust = 0)
    )

  print(p_pt)

  outfile_pdf <- out_file(paste0("Myeloid_founder_clones_", pt, ".pdf"))
  outfile_eps <- out_file(paste0("Myeloid_founder_clones_", pt, ".eps"))

  ggsave(
    filename     = outfile_pdf,
    plot         = p_pt,
    width        = 7,
    height       = 5,
    dpi          = 300,
    useDingbats  = FALSE
  )

  ggsave(
    filename = outfile_eps,
    plot     = p_pt,
    width    = 7,
    height   = 5,
    dpi      = 300,
    device   = cairo_ps
  )
}


############################################################
## PART B.2 – [study-specific]: DNMT3A+ T cells (CAR+ vs CAR−) over time
############################################################

CAR_GREEN <- "darkgreen"
patient_of_interest <- study_setting("study_value_002")

paper_of_interest <- id_key %>%
  filter(CTN_ID == patient_of_interest) %>%
  pull(Paper_ID) %>%
  { if (length(.) == 0 || is.na(.)) patient_of_interest else . }

## Days present for this patient (e.g. -5, 7, 1043)
days_0114 <- df3 %>%
  dplyr::filter(Patient == patient_of_interest) %>%
  dplyr::pull(Day) %>%
  unique() %>%
  sort()

## T cells only for this patient, + CH + CAR annotations
df_0114_T <- df3 %>%
  dplyr::filter(Patient == patient_of_interest) %>%
  dplyr::mutate(
    ## treat these as T-cell compartment
    is_Tcell_0114 = OmiqFilter %in% c("CD4_Tcell", "CD8_Tcell", "Lymphoid")
  ) %>%
  dplyr::filter(is_Tcell_0114) %>%
  dplyr::mutate(
    ## CH gene & DNMT3A flag from Variant_Call
    CH_gene = dplyr::if_else(
      Variant_Call == "WT",
      NA_character_,
      sub(":.*$", "", Variant_Call)
    ),
    is_DNMT3A = CH_gene == "DNMT3A",

    ## CAR compartment from CAR column
    CAR_status = dplyr::if_else(CAR == "CAR", "CAR+", "CAR-")
  )

cat(study_setting("study_value_003"))
df_0114_T %>%
  dplyr::count(Day, CAR_status) %>%
  print(n = Inf)

## Summarise: fraction of T cells that are DNMT3A+ in each CAR compartment
summary_dnmt3a <- df_0114_T %>%
  dplyr::group_by(Day, CAR_status) %>%
  dplyr::summarise(
    n_Tcells      = dplyr::n(),
    n_dnmt3a_T    = sum(is_DNMT3A, na.rm = TRUE),
    frac_dnmt3a_T = dplyr::if_else(
      n_Tcells > 0,
      n_dnmt3a_T / n_Tcells,
      0
    ),
    .groups = "drop"
  )

## Fill in missing Day × CAR_status combos with 0
grid <- tidyr::expand_grid(
  Day        = days_0114,
  CAR_status = c("CAR+", "CAR-")
)

summ_dnmt3a_grid <- grid %>%
  dplyr::left_join(summary_dnmt3a, by = c("Day", "CAR_status")) %>%
  dplyr::mutate(
    n_Tcells      = tidyr::replace_na(n_Tcells, 0L),
    n_dnmt3a_T    = tidyr::replace_na(n_dnmt3a_T, 0L),
    frac_dnmt3a_T = tidyr::replace_na(frac_dnmt3a_T, 0),
    time_label    = factor(
      Day,
      levels = days_0114,
      labels = paste0("D", days_0114)
    ),
    group_label   = dplyr::if_else(
      CAR_status == "CAR+",
      "DNMT3A+ CAR+",
      "DNMT3A+ CAR-"
    )
  )

cat(study_setting("study_value_004"))
summ_dnmt3a_grid %>%
  dplyr::arrange(Day, CAR_status) %>%
  dplyr::mutate(percent_dnmt3a_T = 100 * frac_dnmt3a_T) %>%
  dplyr::select(time_label, CAR_status, n_Tcells, n_dnmt3a_T, percent_dnmt3a_T) %>%
  print(n = Inf)

## Plot: DNMT3A+ CAR+ vs CAR- over time (title by Paper_ID)
p_dnmt3_0114 <- ggplot(
  summ_dnmt3a_grid,
  aes(
    x     = time_label,
    y     = frac_dnmt3a_T,
    group = group_label,
    color = group_label
  )
) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    values = c(
      "DNMT3A+ CAR+" = unname(pal_variant["DNMT3A"]),
      "DNMT3A+ CAR-" = "grey40"
    ),
    breaks = c("DNMT3A+ CAR+", "DNMT3A+ CAR-"),
    name   = NULL
  ) +
  labs(
    title = paste0(paper_of_interest, ": DNMT3A mutant T cells over time"),
    x     = "Timepoint",
    y     = "Percent of T cells that are DNMT3A+"
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

print(p_dnmt3_0114)

ggsave(
  out_file(study_setting("study_value_005")),
  plot   = p_dnmt3_0114,
  width  = 8,
  height = 5,
  units  = "in"
)

##### test if there is CD4 predominance in CH+ cells at expansion #####
############################################################
## CH+CAR+ enrichment in CD4 vs CD8 (at PEAK, by Paper_ID)
############################################################

## Restrict to PEAK, CAR+ T cells (CD4/CD8 only)
df_peak_car_T <- df_peak %>%
  filter(
    is_Tcell,
    is_CAR,
    Cell_Type5 %in% c("CD4 T-cell", "CD8 T-cell")
  )

cat("\n=== [PART A] PEAK: CAR+ T cells (CD4/CD8 only, by Paper_ID) ===\n")
df_peak_car_T %>%
  count(patient_paper, Cell_Type5, is_CH) %>%
  print(n = Inf)

## 1) Global enrichment test: CH+ vs CH- in CD4 vs CD8
counts_CH_lineage <- df_peak_car_T %>%
  mutate(
    lineage = if_else(Cell_Type5 == "CD4 T-cell", "CD4", "CD8"),
    CH_flag = if_else(is_CH, "CH+", "CH-")
  ) %>%
  count(lineage, CH_flag) %>%
  tidyr::pivot_wider(
    names_from  = CH_flag,
    values_from = n,
    values_fill = 0
  )

cat("\n=== [PART A] PEAK: pooled counts of CAR+ T cells by lineage × CH status ===\n")
print(counts_CH_lineage, n = Inf)

if (all(c("CD4","CD8") %in% counts_CH_lineage$lineage)) {
  get_val <- function(lin, ch) {
    counts_CH_lineage %>%
      filter(lineage == lin) %>%
      pull(!!sym(ch)) %>%
      { if (length(.) == 0) 0L else . }
  }

  mat_CH_enrich <- matrix(
    c(
      get_val("CD4","CH+"),
      get_val("CD4","CH-"),
      get_val("CD8","CH+"),
      get_val("CD8","CH-")
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Lineage = c("CD4","CD8"),
      CH      = c("CH+","CH-")
    )
  )

  cat("\nContingency table: Lineage (rows) × CH status (cols) within CAR+ T cells at PEAK:\n")
  print(mat_CH_enrich)

  ft_CH_enrich <- fisher.test(mat_CH_enrich)

  cat("\nFisher test: enrichment of CH+ in CD4 vs CD8 (within CAR+ T cells at PEAK)\n")
  print(ft_CH_enrich)
  cat(
    sprintf(
      "log2(OR) = %.3f\n\n",
      log2(unname(ft_CH_enrich$estimate))
    )
  )
} else {
  cat("\n[WARN] Missing either CD4 or CD8 in PEAK CAR+ T cells – skipping global enrichment test.\n")
}

## 2) Per-patient enrichment: fraction of CH+ among CAR+,
##    CD4 vs CD8 (Paper_ID-based)
per_pt_lineage_CHfrac <- df_peak_car_T %>%
  mutate(
    lineage = if_else(Cell_Type5 == "CD4 T-cell", "CD4", "CD8")
  ) %>%
  group_by(patient_paper, lineage) %>%
  summarise(
    n_car_T  = n(),
    n_car_CH = sum(is_CH),
    frac_CH  = if_else(n_car_T > 0, n_car_CH / n_car_T, NA_real_),
    .groups  = "drop"
  )

cat("\n=== [PART A] PEAK: per-patient fraction CH+ among CAR+ T cells (CD4 vs CD8) ===\n")
print(per_pt_lineage_CHfrac, n = Inf)

paired_CHfrac <- per_pt_lineage_CHfrac %>%
  select(patient_paper, lineage, frac_CH) %>%
  tidyr::pivot_wider(
    names_from  = lineage,
    values_from = frac_CH
  ) %>%
  filter(!is.na(CD4), !is.na(CD8))

cat("\nPatients (Paper_ID) with both CD4 and CD8 CAR+ T cells at PEAK (for paired test):\n")
print(paired_CHfrac, n = Inf)

if (nrow(paired_CHfrac) > 0) {
  w_CHfrac_cd4_vs_cd8 <- wilcox.test(
    paired_CHfrac$CD4,
    paired_CHfrac$CD8,
    paired      = TRUE,
    alternative = "greater"  # test CD4 > CD8
  )

  cat("\nWilcoxon (paired) test: fraction CH+ among CAR+ CD4 vs CD8 T cells at PEAK\n")
  print(w_CHfrac_cd4_vs_cd8)
} else {
  cat("\n[WARN] No patients with both CD4 and CD8 CAR+ T cells at PEAK – skipping paired test.\n")
}


############### Part C – CAR summaries & PEAK T-cell CH analyses ##################

## 1) CAR overall summaries ------------------------------------

## (a) Per-sample summary
sample_summary <- df_all %>%
  group_by(Patient_ID, patient_paper, sample_timepoint) %>%
  summarise(
    n_cells  = n(),
    n_CAR    = sum(CAR_detected == "Present"),
    frac_CAR = n_CAR / pmax(n_cells, 1L),
    .groups  = "drop"
  )

cat("\n=== Per-sample CAR summary (all samples) ===\n")
sample_summary %>%
  arrange(patient_paper, sample_timepoint) %>%
  tibble::as_tibble() %>%
  print(n = Inf)

## (b) Timepoint-level summary
timepoint_car_summary <- sample_summary %>%
  group_by(sample_timepoint) %>%
  summarise(
    n_samples          = n(),
    n_samples_with_CAR = sum(n_CAR > 0),
    total_cells        = sum(n_cells),
    total_CAR          = sum(n_CAR),
    frac_cells_CAR     = total_CAR / pmax(total_cells, 1L),
    .groups            = "drop"
  ) %>%
  arrange(sample_timepoint)

cat("\n=== CAR presence by sample timepoint ===\n")
timepoint_car_summary %>%
  tibble::as_tibble() %>%
  print(n = Inf)

## (c) Pre vs non-pre (based on "pre" in sample_timepoint)
df_all <- df_all %>%
  mutate(pre_flag = str_detect(tolower(sample_timepoint), "pre"))

pre_vs_post_car <- df_all %>%
  group_by(pre_flag) %>%
  summarise(
    n_cells  = n(),
    n_CAR    = sum(CAR_detected == "Present"),
    frac_CAR = n_CAR / pmax(n_cells, 1L),
    .groups  = "drop"
  )

cat("\n=== CAR in pre vs non-pre samples ===\n")
pre_vs_post_car %>%
  tibble::as_tibble() %>%
  print(n = Inf)


## 2) Per-sample lineage × CAR composition ---------------------

sample_lineage_car <- df_all %>%
  group_by(Patient_ID, patient_paper, sample_timepoint,
           Cell_Type5, CAR_detected) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Patient_ID, patient_paper, sample_timepoint, Cell_Type5) %>%
  mutate(
    n_total_lineage     = sum(n),
    frac_within_lineage = n / pmax(n_total_lineage, 1L)
  ) %>%
  ungroup()

cat("\n=== Per-sample lineage x CAR counts (first 50 rows) ===\n")
sample_lineage_car %>%
  arrange(patient_paper, sample_timepoint, Cell_Type5, CAR_detected) %>%
  tibble::as_tibble() %>%
  print(n = 50)

p_sample_comp <- ggplot(
  sample_lineage_car %>% filter(Cell_Type5 != "Unknown"),
  aes(
    x     = sample_timepoint,
    y     = n,
    fill  = Cell_Type5,
    alpha = CAR_detected
  )
) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = pal_lineage5[c("CD4 T-cell","CD8 T-cell","NK cell","Myeloid","Myeloid Progenitor")],
    name   = "Lineage"
  ) +
  scale_alpha_manual(
    values = c("Present" = 0.95, "Absent" = 0.4),
    name   = "CAR status"
  ) +
  facet_wrap(~ patient_paper, scales = "free_y") +
  labs(
    title = "Per-sample lineage composition and CAR status",
    x = "Sample timepoint",
    y = "Cell count"
  ) +
  theme_pub +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom",
    legend.direction = "horizontal"
  )


## 3) CAR fraction by timepoint & lineage ----------------------

car_by_time_lineage <- df_all %>%
  group_by(sample_timepoint, Cell_Type5) %>%
  summarise(
    n_cells  = n(),
    n_CAR    = sum(CAR_detected == "Present"),
    frac_CAR = n_CAR / pmax(n_cells, 1L),
    .groups  = "drop"
  )

cat("\n=== CAR fraction by timepoint and lineage ===\n")
car_by_time_lineage %>%
  arrange(sample_timepoint, Cell_Type5) %>%
  tibble::as_tibble() %>%
  print(n = Inf)

p_car_by_time_lineage <- ggplot(
  car_by_time_lineage %>% filter(Cell_Type5 != "Unknown"),
  aes(
    x     = sample_timepoint,
    y     = frac_CAR,
    group = Cell_Type5,
    color = Cell_Type5
  )
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = pal_lineage5[c("CD4 T-cell","CD8 T-cell","NK cell","Myeloid","Myeloid Progenitor")],
    name   = "Lineage"
  ) +
  labs(
    title = "CAR+ fraction by lineage over time",
    x = "Sample timepoint",
    y = "Fraction of cells that are CAR+"
  ) +
  theme_pub +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


## 4) PEAK-only CAR summaries (T-cell & compartments) ----------

## Use df_peak defined earlier from peak_lookup (Paper_ID-based)
## df_peak has only the peak timepoint per patient_paper

## (a) Fraction of T cells that are CAR+ vs non-CAR (per patient)
tcell_fraction_peak <- df_peak %>%
  filter(is_Tcell) %>%
  group_by(patient_paper) %>%
  summarise(
    n_T        = n(),
    n_T_CAR    = sum(CAR_detected == "Present"),
    frac_T_CAR = n_T_CAR / pmax(n_T, 1L),
    .groups    = "drop"
  ) %>%
  mutate(frac_T_nonCAR = 1 - frac_T_CAR) %>%
  pivot_longer(
    cols      = c(frac_T_CAR, frac_T_nonCAR),
    names_to  = "CAR_status",
    values_to = "fraction"
  ) %>%
  mutate(
    CAR_status = recode(
      CAR_status,
      frac_T_CAR    = "CAR+ T cells",
      frac_T_nonCAR = "Non-CAR T cells"
    ),
    CAR_status = factor(CAR_status,
                        levels = c("Non-CAR T cells","CAR+ T cells"))
  )

p_tcell_CAR_box <- ggplot(
  tcell_fraction_peak,
  aes(x = CAR_status, y = fraction)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    fill  = "grey90",
    color = "black"
  ) +
  geom_jitter(
    width = 0.1, height = 0,
    size  = 1.6, alpha = 0.7
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "PEAK: CAR+ vs non-CAR T cells per patient",
    x = NULL,
    y = "Fraction of total T cells"
  ) +
  theme_pub

## (b) Fraction of CAR+ cells per compartment (PEAK)
compartment_CAR_frac_peak <- df_peak %>%
  filter(Cell_Type5 != "Unknown") %>%
  group_by(patient_paper, Cell_Type5) %>%
  summarise(
    n_comp     = n(),
    n_comp_CAR = sum(CAR_detected == "Present"),
    frac_CAR   = n_comp_CAR / pmax(n_comp, 1L),
    .groups    = "drop"
  )

p_compartment_CAR_box <- ggplot(
  compartment_CAR_frac_peak,
  aes(x = Cell_Type5, y = frac_CAR, fill = Cell_Type5)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.9
  ) +
  geom_jitter(
    width = 0.15, height = 0,
    size  = 1.6, alpha = 0.7
  ) +
  scale_fill_manual(
    values = pal_lineage5[c("CD4 T-cell","CD8 T-cell","NK cell","Myeloid","Myeloid Progenitor")],
    guide  = "none"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "PEAK: fraction of cells that are CAR+ in each compartment",
    x = "Compartment",
    y = "Fraction CAR+"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  out_file("CAR_composition_Boxplot.eps"),
  plot  = p_compartment_CAR_box,
  width = 8, height = 5, units = "in",
  device = cairo_ps
)


## 5) PEAK T-cell CH gene enrichment (CAR+ vs CAR−) ------------

df_peak_T_CH <- df_peak %>%
  filter(is_Tcell, is_CH, !is.na(CH_gene_simple))

gene_counts_peak_T <- df_peak_T_CH %>%
  count(CH_gene_simple, sort = TRUE)

cat("\n=== [PART C] PEAK T cells: CH+ counts per gene ===\n")
print(gene_counts_peak_T, n = Inf)

genes_to_test_peak_T <- gene_counts_peak_T %>%
  filter(n >= 20) %>%
  pull(CH_gene_simple)

enrich_peak_T_list <- lapply(genes_to_test_peak_T, function(g) {
  a <- sum(df_peak_T_CH$CH_gene_simple == g & df_peak_T_CH$is_CAR)      # gene g, CAR+
  b <- sum(df_peak_T_CH$CH_gene_simple == g & !df_peak_T_CH$is_CAR)     # gene g, CAR-
  c <- sum(df_peak_T_CH$CH_gene_simple != g & df_peak_T_CH$is_CAR)      # other genes, CAR+
  d <- sum(df_peak_T_CH$CH_gene_simple != g & !df_peak_T_CH$is_CAR)     # other genes, CAR-

  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      gene = c(g, "Other"),
      CAR  = c("CAR+","CAR-")
    )
  )

  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(tibble(
      CH_gene      = g,
      log2_OR      = NA_real_,
      log2_OR_low  = NA_real_,
      log2_OR_high = NA_real_,
      p_value      = NA_real_
    ))
  }

  ft <- fisher.test(mat)

  tibble(
    CH_gene      = g,
    log2_OR      = log2(unname(ft$estimate)),
    log2_OR_low  = log2(ft$conf.int[1]),
    log2_OR_high = log2(ft$conf.int[2]),
    p_value      = ft$p.value
  )
})

enrich_peak_T_df <- bind_rows(enrich_peak_T_list) %>%
  arrange(desc(log2_OR))

cat("\n=== [PART C] PEAK T cells: enrichment of CH genes in CAR+ vs CAR- (log2 OR) ===\n")
print(enrich_peak_T_df, n = Inf)

## Optional forest plot for PEAK T-cell enrichment
enrich_peak_T_df_plot <- enrich_peak_T_df %>%
  filter(!is.na(log2_OR))

p_enrich_peak_T <- ggplot(
  enrich_peak_T_df_plot,
  aes(
    x = log2_OR,
    y = forcats::fct_reorder(CH_gene, log2_OR, .na_rm = TRUE),
    color = CH_gene
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = log2_OR_low, xmax = log2_OR_high), height = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(values = pal_CH_category, guide = "none") +
  xlab("log2 odds ratio for CAR+ vs CAR- (PEAK CH+ T cells)") +
  ylab("CH gene") +
  ggtitle("PEAK T cells: enrichment / depletion of CH in CAR+ cells") +
  theme_classic(base_size = 18)


## 6) PEAK T-cell CAR enrichment within each CH gene ----------

car_within_CH_peak_list <- lapply(genes_to_test_peak_T, function(g) {
  dat_g <- df_peak_T_CH %>% filter(CH_gene_simple == g)
  n_CAR    <- sum(dat_g$is_CAR)
  n_nonCAR <- sum(!dat_g$is_CAR)
  n_tot    <- n_CAR + n_nonCAR

  if (n_tot == 0) {
    return(tibble(
      CH_gene        = g,
      n_CH_T_cells   = 0L,
      n_CAR_CH       = 0L,
      frac_CAR_CH    = NA_real_,
      p_value_binom  = NA_real_
    ))
  }

  bt <- binom.test(n_CAR, n_tot, p = 0.5)  # null: 50% CAR+

  tibble(
    CH_gene        = g,
    n_CH_T_cells   = n_tot,
    n_CAR_CH       = n_CAR,
    frac_CAR_CH    = n_CAR / n_tot,
    p_value_binom  = bt$p.value
  )
})

car_within_CH_peak_df <- bind_rows(car_within_CH_peak_list) %>%
  arrange(desc(frac_CAR_CH))

cat("\n=== [PART C] PEAK T cells: CAR+ fraction within CH+ cells for each gene ===\n")
print(car_within_CH_peak_df, n = Inf)


## 7) PEAK per-patient % of CAR+ T cells that are CH+ ----------

car_CH_T_per_patient <- df_peak %>%
  filter(is_Tcell, is_CAR) %>%
  group_by(patient_paper) %>%
  summarise(
    n_CAR_T      = n(),
    n_CAR_T_CH   = sum(is_CH),
    frac_CAR_T_CH = n_CAR_T_CH / pmax(n_CAR_T, 1L),
    .groups      = "drop"
  )

cat("\n=== [PART C] PEAK: per-patient fraction of CAR+ T cells that are CH+ ===\n")
print(car_CH_T_per_patient, n = Inf)

## 8) PEAK per-patient enrichment of CH+ in CAR+ vs CAR- T cells ----

## All PEAK T cells (CH+ and CH-)
df_peak_T_all <- df_peak %>%
  filter(is_Tcell) %>%
  mutate(
    CAR_flag = if_else(CAR_detected == "Present", "CAR+", "CAR-"),
    CH_flag  = if_else(is_CH, "CH+", "CH-")
  )

## (a) Per-patient CH+ fraction among CAR+ vs CAR- T cells
per_pt_CH_frac_CAR <- df_peak_T_all %>%
  count(patient_paper, CAR_flag, CH_flag) %>%
  tidyr::complete(
    patient_paper,
    CAR_flag = c("CAR+","CAR-"),
    CH_flag  = c("CH+","CH-"),
    fill = list(n = 0L)
  ) %>%
  group_by(patient_paper) %>%
  summarise(
    n_CAR_T        = sum(n[CAR_flag == "CAR+"]),
    n_nonCAR_T     = sum(n[CAR_flag == "CAR-"]),
    n_CAR_CH       = n[CAR_flag == "CAR+" & CH_flag == "CH+"],
    n_nonCAR_CH    = n[CAR_flag == "CAR-" & CH_flag == "CH+"],
    frac_CH_in_CAR    = if_else(n_CAR_T > 0,    n_CAR_CH    / n_CAR_T,    NA_real_),
    frac_CH_in_nonCAR = if_else(n_nonCAR_T > 0, n_nonCAR_CH / n_nonCAR_T, NA_real_),
    .groups = "drop"
  )

cat("\n=== [PART C] PEAK: per-patient CH+ fraction among CAR+ vs CAR- T cells ===\n")
print(per_pt_CH_frac_CAR, n = Inf)

## (b) Per-patient Fisher test: is CH+ enriched in CAR+ vs CAR-?
per_pt_counts <- df_peak_T_all %>%
  count(patient_paper, CAR_flag, CH_flag) %>%
  tidyr::complete(
    patient_paper,
    CAR_flag = c("CAR+","CAR-"),
    CH_flag  = c("CH+","CH-"),
    fill = list(n = 0L)
  )

per_patient_list <- per_pt_counts %>%
  group_split(patient_paper)

per_patient_enrichment <- lapply(per_patient_list, function(dat) {
  pt <- unique(dat$patient_paper)

  tab <- dat %>%
    select(CAR_flag, CH_flag, n) %>%
    tidyr::pivot_wider(
      names_from  = CH_flag,
      values_from = n
    ) %>%
    tibble::column_to_rownames("CAR_flag") %>%
    as.matrix()

  ## enforce consistent ordering
  tab <- tab[c("CAR+","CAR-"), c("CH+","CH-"), drop = FALSE]

  if (any(is.na(tab)) ||
      any(rowSums(tab) == 0) ||
      any(colSums(tab) == 0)) {
    tibble(
      patient_paper = pt,
      n_total_T     = sum(tab, na.rm = TRUE),
      log2_OR       = NA_real_,
      log2_OR_low   = NA_real_,
      log2_OR_high  = NA_real_,
      p_value       = NA_real_
    )
  } else {
    ft <- fisher.test(tab)
    tibble(
      patient_paper = pt,
      n_total_T     = sum(tab),
      log2_OR       = log2(unname(ft$estimate)),
      log2_OR_low   = log2(ft$conf.int[1]),
      log2_OR_high  = log2(ft$conf.int[2]),
      p_value       = ft$p.value
    )
  }
})

per_patient_enrichment <- bind_rows(per_patient_enrichment) %>%
  arrange(desc(log2_OR))

cat("\n=== [PART C] PEAK: per-patient enrichment of CH+ in CAR+ vs CAR- T cells ===\n")
print(per_patient_enrichment, n = Inf)

## (c) Across-patient test: are log2 ORs > 0 on average?
valid_enrich <- per_patient_enrichment %>% filter(!is.na(log2_OR))

if (nrow(valid_enrich) > 0) {
  w_per_patient <- wilcox.test(
    valid_enrich$log2_OR,
    mu          = 0,
    alternative = "greater"  # test enrichment of CH+ in CAR+ vs CAR-
  )

  cat("\nWilcoxon one-sample test (log2 OR > 0):\n")
  print(w_per_patient)
} else {
  cat("\n[PART C] No patients with sufficient counts for per-patient enrichment test.\n")
}

## === 9) Plot per-patient enrichment of CH+ in CAR+ vs CAR- T cells ===

per_patient_enrichment_plot_df <- per_patient_enrichment %>%
  # drop patients where we couldn't estimate OR
  filter(!is.na(log2_OR), is.finite(log2_OR)) %>%
  mutate(
    direction = case_when(
      log2_OR > 0 & p_value < 0.05 ~ "Enriched in CAR+",
      log2_OR < 0 & p_value < 0.05 ~ "Depleted in CAR+",
      TRUE                         ~ "NS"
    ),
    direction = factor(
      direction,
      levels = c("Enriched in CAR+","Depleted in CAR+","NS")
    )
  )

cat("\n=== [PART C] Plotting per-patient enrichment (log2 OR) ===\n")
print(per_patient_enrichment_plot_df, n = Inf)

p_per_patient_enrichment <- ggplot(
  per_patient_enrichment_plot_df,
  aes(
    x     = log2_OR,
    y     = forcats::fct_reorder(patient_paper, log2_OR),
    xmin  = log2_OR_low,
    xmax  = log2_OR_high,
    color = direction
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(height = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Enriched in CAR+" = "#1b9e77",
      "Depleted in CAR+" = "#d95f02",
      "NS"               = "grey50"
    ),
    name = NULL
  ) +
  xlab("log2 odds ratio (CH+ in CAR+ vs CAR- T cells)") +
  ylab("Patient") +
  ggtitle("PEAK T cells: per-patient enrichment of CH+ in CAR+ vs CAR-") +
  theme_pub +
  theme(legend.position = "bottom")

print(p_per_patient_enrichment)

ggsave(
  out_file("per_patient_CH_enrichment.eps"),
  plot  = p_per_patient_enrichment,
  width = 8, height = 5, units = "in",
  device = cairo_ps
)


############################################################
## PART D – Surface marker overlays on t-SNE (old-style overlays)
## Replicates your previous marker highlighting panels using df_all.
## Saves into the same folder as your other new plots.
############################################################

# --- output directory (match your other new plots) ---
OUT_DIR_MARKERS <- out_dir
dir.create(OUT_DIR_MARKERS, recursive = TRUE, showWarnings = FALSE)

# --- match old palette + limits ---
marker_cols_pal    <- c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07")
marker_limits_hard <- c(-5, 5)

# --- helper: normalize strings for tolerant matching (handles weird symbols like ε) ---
.norm_str <- function(x) {
  x %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "")
}

# --- helper: find a marker column from a list of candidate names ---
find_marker_col <- function(df, candidates) {
  cn <- colnames(df)

  # 1) exact match
  hit <- candidates[candidates %in% cn]
  if (length(hit) > 0) return(hit[[1]])

  # 2) case-insensitive exact match
  cn_low <- tolower(cn)
  cand_low <- tolower(candidates)
  hit2 <- cn[match(cand_low, cn_low)]
  hit2 <- hit2[!is.na(hit2)]
  if (length(hit2) > 0) return(hit2[[1]])

  # 3) normalized match (strip punctuation/special chars)
  cn_norm <- .norm_str(cn)
  cand_norm <- .norm_str(candidates)
  m <- match(cand_norm, cn_norm)
  hit3 <- cn[m[!is.na(m)]]
  if (length(hit3) > 0) return(hit3[[1]])

  NA_character_
}

# --- helper: placeholder plot when marker is missing ---
missing_marker_plot <- function(title_txt) {
  ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label = paste0(title_txt, "\nMISSING"),
      size = 6, fontface = "bold"
    ) +
    theme_void() +
    ggtitle("")
}

# --- helper: single marker t-SNE plot (old-style gradient + your axis arrows theme) ---
plot_marker_tsne <- function(df, marker_col, title_txt = NULL,
                             limits = marker_limits_hard) {
  if (is.null(title_txt)) title_txt <- marker_col

  # ensure numeric
  vals <- suppressWarnings(as.numeric(df[[marker_col]]))
  dfp  <- df %>% mutate(.marker_val = vals)

  ggplot(dfp, aes(x = optsne_1, y = optsne_2, color = .marker_val)) +
    geom_point_ai_rast(size = 0.2, alpha = 0.9) +
    xlab("t-SNE1") + ylab("t-SNE2") +
    scale_color_gradientn(
      colors = marker_cols_pal,
      limits = limits,
      oob    = scales::squish,
      name   = title_txt
    ) +
    theme_tsne() +
    guides(
      x     = axis_trunc,
      y     = axis_trunc,
      color = guide_colorbar(barheight = unit(3.0, "cm"))
    )
}

# --- helper: build a row/panel from a named marker spec ---
make_marker_row <- function(df, marker_spec_named, ncol = 4) {
  plots <- lapply(names(marker_spec_named), function(nm) {
    col_found <- find_marker_col(df, marker_spec_named[[nm]])
    if (is.na(col_found)) {
      message(
        "[Marker plot] Missing marker column for: ", nm,
        " (candidates: ", paste(marker_spec_named[[nm]], collapse = ", "), ")"
      )
      missing_marker_plot(nm)
    } else {
      message("[Marker plot] Using column '", col_found, "' for marker: ", nm)
      plot_marker_tsne(df, col_found, title_txt = nm)
    }
  })
  wrap_plots(plots, ncol = ncol)
}

# ------------------------------------------------------------
# Marker sets (rational, and mirrors your previous intent)
# ------------------------------------------------------------

# Progenitor / stem-ish
markers_progenitor <- list(
  "CD34"  = c("CD34"),
  "CD117" = c("CD117", "KIT"),
  "CD38"  = c("CD38"),
  "CD90"  = c("CD90", "THY1")
)

# Myeloid
markers_myeloid <- list(
  "CD16"     = c("CD16", "FCGR3A"),
  "CD62L"    = c("CD62L", "SELL"),
  "CD14"     = c("CD14"),
  "FcεRIα"   = c("FcεRIα", "FcERIa", "FceRIa", "FCER1A", "FceR1a", "Fcer1a",
                "FcER1A", "FCeRIa", "FceRIalpha", "FcεRIa")
)

# B cell
markers_bcell <- list(
  "CD10" = c("CD10", "MME"),
  "CD19" = c("CD19")
)

# T cell (CD45RA/RO will take first available of those)
markers_tcell <- list(
  "CD3"       = c("CD3", "CD3E", "CD3e", "CD3_Total"),
  "CD4"       = c("CD4"),
  "CD8"       = c("CD8", "CD8A", "CD8a"),
  "CD45RA/RO" = c("CD45RA", "CD45RO", "PTPRC", "CD45")
)

# NK cell
markers_nk <- list(
  "CD56"       = c("CD56", "NCAM1", "NCAM-1"),
  "CD7"       = c("CD7")
)


# ------------------------------------------------------------
# Build marker panels
# ------------------------------------------------------------
p_markers_progenitor <- make_marker_row(df_all, markers_progenitor, ncol = 4)
p_markers_myeloid    <- make_marker_row(df_all, markers_myeloid,    ncol = 4)
p_markers_tcell      <- make_marker_row(df_all, markers_tcell,      ncol = 4)
p_markers_nk         <- make_marker_row(df_all, markers_nk,         ncol = 2)
p_markers_bcell      <- make_marker_row(df_all, markers_bcell,      ncol = 2)

# ------------------------------------------------------------
# Save panels (PNG + EPS), same style/dimensions as your old code
# ------------------------------------------------------------

# Progenitor
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "myeloid_progenitor_markers.png"),
  plot     = p_markers_progenitor,
  width    = 21, height = 4.06, units = "in", dpi = 300
)
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "myeloid_progenitor_markers.eps"),
  plot     = p_markers_progenitor,
  width    = 21, height = 4.06, units = "in",
  device   = cairo_ps
)

# Myeloid
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "myeloid_markers.png"),
  plot     = p_markers_myeloid,
  width    = 21, height = 4.06, units = "in", dpi = 300
)
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "myeloid_markers.eps"),
  plot     = p_markers_myeloid,
  width    = 21, height = 4.06, units = "in",
  device   = cairo_ps
)

# T cell
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "t_cell_markers.png"),
  plot     = p_markers_tcell,
  width    = 21, height = 4.06, units = "in", dpi = 300
)
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "t_cell_markers.eps"),
  plot     = p_markers_tcell,
  width    = 21, height = 4.06, units = "in",
  device   = cairo_ps
)

# NK cell
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "nk_cell_markers.png"),
  plot     = p_markers_nk,
  width    = 10.5, height = 4.06, units = "in", dpi = 300
)
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "nk_cell_markers.eps"),
  plot     = p_markers_nk,
  width    = 21, height = 4.06, units = "in",
  device   = cairo_ps
)

# B cell
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "b_cell_markers.png"),
  plot     = p_markers_bcell,
  width    = 10.5, height = 4.06, units = "in", dpi = 300
)
ggsave(
  filename = file.path(OUT_DIR_MARKERS, "b_cell_markers.eps"),
  plot     = p_markers_bcell,
  width    = 10.5, height = 4.06, units = "in",
  device   = cairo_ps
)

# Optional: print (useful if running interactively)
print(p_markers_progenitor)
print(p_markers_myeloid)
print(p_markers_tcell)
print(p_markers_nk)
print(p_markers_bcell)


############################################################
## PART E – Variant printout by patient, lineage, clone (adjacent columns)
############################################################

# Where to save (same folder as other outputs)
OUT_DIR_VARIANTS <- out_dir
dir.create(OUT_DIR_VARIANTS, recursive = TRUE, showWarnings = FALSE)

# Optional: set to >1 if you want to drop ultra-rare clones
MIN_CLONE_CELLS <- 1

# ---- Parse variants into "clone" + "member variants" (like your df3 logic) ----
variants_by_patient_lineage <- df_peak %>%
  dplyr::filter(!is.na(Variant_Call), Variant_Call != "WT") %>%
  dplyr::mutate(
    # point mutations typically look like GENE:p.X123Y (keep founder compact)
    is_point_mut = stringr::str_detect(Variant_Call, "^[^:]+:p\\."),

    # founder label: first "GENE:p.xxx" token if present; otherwise keep whole call
    Base_Variant = dplyr::if_else(
      is_point_mut,
      stringr::str_extract(Variant_Call, "^[^:]+:p\\.[^:_]+"),
      Variant_Call
    ),

    # gene for ordering/annotation
    Gene = gene_from_variant(Base_Variant),

    # clone member label: founder-only stays compact; subclones keep full Variant_Call
    Clone_Variant = dplyr::case_when(
      is_point_mut & Variant_Call == Base_Variant ~ Base_Variant,
      is_point_mut                                ~ Variant_Call,
      TRUE                                        ~ Variant_Call
    )
  ) %>%
  # collapse at patient + lineage + founder clone
  dplyr::group_by(patient_paper, Cell_Type5, Base_Variant, Gene) %>%
  dplyr::summarise(
    n_cells_clone = dplyr::n(),
    variants      = list(sort(unique(Clone_Variant))),
    .groups       = "drop"
  ) %>%
  dplyr::filter(n_cells_clone >= MIN_CLONE_CELLS) %>%
  dplyr::mutate(
    # order variants with founder first, then the rest (so they land adjacent as Variant_1, Variant_2...)
    variants_ordered = purrr::map2(Base_Variant, variants, function(base, v) {
      v <- unique(v)
      c(base, setdiff(v, base))
    })
  ) %>%
  dplyr::select(patient_paper, Cell_Type5, Gene, Base_Variant, n_cells_clone, variants_ordered) %>%
  tidyr::unnest_longer(variants_ordered, values_to = "Variant") %>%
  dplyr::group_by(patient_paper, Cell_Type5, Gene, Base_Variant, n_cells_clone) %>%
  dplyr::mutate(var_idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    names_from  = var_idx,
    values_from = Variant,
    names_prefix = "Variant_"
  ) %>%
  # order by lineage (Cell_Type5 factor order), then patient, then gene/founder
  dplyr::arrange(Cell_Type5, patient_paper, Gene, Base_Variant, dplyr::desc(n_cells_clone))

# ---- Print to console (full) ----
cat("\n=== Variants by patient × lineage × clone (adjacent columns) ===\n")
print(variants_by_patient_lineage, n = Inf, width = Inf)

# ---- Save as TSV ----
variant_outfile <- file.path(OUT_DIR_VARIANTS, "variants_by_patient_lineage_clone.tsv")
readr::write_tsv(variants_by_patient_lineage, variant_outfile)

cat("\n[WROTE] ", variant_outfile, "\n", sep = "")

# ---- Optional: a compact per-patient summary (one line per patient) ----
variants_by_patient_compact <- variants_by_patient_lineage %>%
  dplyr::mutate(
    # collapse the Variant_* columns into a single string per clone row
    variants_concat = apply(
      dplyr::select(., dplyr::starts_with("Variant_")),
      1,
      function(x) paste(na.omit(x), collapse = "; ")
    ),
    clone_label = paste0(Cell_Type5, " | ", variants_concat)
  ) %>%
  dplyr::group_by(patient_paper) %>%
  dplyr::summarise(
    clones = paste(clone_label, collapse = "  ||  "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(patient_paper)

cat("\n=== Compact per-patient variant summary (ordered by lineage within each patient) ===\n")
print(variants_by_patient_compact, n = Inf, width = Inf)

compact_outfile <- file.path(OUT_DIR_VARIANTS, "variants_by_patient_compact.txt")
readr::write_lines(
  paste0(variants_by_patient_compact$patient_paper, "\t", variants_by_patient_compact$clones),
  compact_outfile
)
cat("\n[WROTE] ", compact_outfile, "\n", sep = "")


############################################################
## PART F – Variant lists (per patient), per-patient t-SNE, and variant×cell-type heatmap
############################################################

OUT_DIR_SUMMARY <- out_dir
dir.create(OUT_DIR_SUMMARY, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Gene priority for ordering variants ("type" order)
# Edit this list to your preferred presentation order.
# ----------------------------
GENE_PRIORITY <- c(
  "DNMT3A","TET2","PPM1D","TP53","BCOR","ATM","EZH2","ASXL1","JAK2",
  "SF3B1","SRSF2","U2AF1","IDH1","IDH2","RUNX1","CBL","KRAS","NRAS",
  "BRCA2","BRCA1","CHEK2","ATR","PALB2","RAD51","FANCA","NBN","MRE11","RAD50",
  "MSH2","MSH6","MLH1","PMS2","POLE","POLD1"
)

gene_rank <- function(g) {
  r <- match(g, GENE_PRIORITY)
  ifelse(is.na(r), length(GENE_PRIORITY) + 999, r)
}


# ----------------------------
# 1) Build a simple variant table (patient × variant)
# ----------------------------

# Parse variants into a compact "Base_Variant" (e.g. GENE:p.X123Y) when possible
df_var <- df_all %>%
  dplyr::filter(!is.na(Variant_Call), Variant_Call != "WT") %>%
  dplyr::mutate(
    is_point_mut = stringr::str_detect(Variant_Call, "^[^:]+:p\\."),

    # Founder/compact label when point mutation is present:
    Base_Variant = dplyr::if_else(
      is_point_mut,
      stringr::str_extract(Variant_Call, "^[^:]+:p\\.[^:_]+"),
      Variant_Call
    ),

    Gene = gene_from_variant(Base_Variant),

    # Optional: collapse TP53 flavors to TP53 at the gene level only
    Gene_simple = dplyr::case_when(
      stringr::str_detect(Gene, "^TP53") ~ "TP53",
      TRUE ~ Gene
    )
  )

# One row per patient × variant (deduplicated), ordered by gene type
variant_patient_tbl <- df_var %>%
  dplyr::group_by(patient_paper, Gene_simple, Base_Variant) %>%
  dplyr::summarise(
    n_cells_variant = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(gene_rank = gene_rank(Gene_simple)) %>%
  dplyr::arrange(patient_paper, gene_rank, Gene_simple, Base_Variant)

cat("\n=== Variants observed per patient (deduplicated; ordered by gene type) ===\n")
print(variant_patient_tbl, n = Inf, width = Inf)

readr::write_tsv(
  variant_patient_tbl %>% dplyr::select(-gene_rank),
  file.path(OUT_DIR_SUMMARY, "variants_by_patient.tsv")
)

# Compact “Patient -> variant list”, but *within each patient* it respects gene priority
patient_variant_list <- variant_patient_tbl %>%
  dplyr::group_by(patient_paper) %>%
  dplyr::arrange(gene_rank, Gene_simple, Base_Variant, .by_group = TRUE) %>%
  dplyr::summarise(
    n_variants = dplyr::n(),
    variants   = paste0(Base_Variant, collapse = "; "),
    .groups    = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_variants), patient_paper)

cat("\n=== Patient-level variant list (within-patient ordered by gene type) ===\n")
print(patient_variant_list, n = Inf, width = Inf)

readr::write_tsv(
  patient_variant_list,
  file.path(OUT_DIR_SUMMARY, "patient_variant_list.tsv")
)


# ----------------------------
# 2) DDR vs age-related / CH co-occurrence summary (for your narrative)
# ----------------------------

# You can edit these lists to match your exact “DDR” and “age-related/CH” definitions
AGE_RELATED_GENES <- c(
  "DNMT3A","TET2","ASXL1","PPM1D","TP53","BCOR","JAK2","SF3B1","SRSF2","U2AF1",
  "IDH1","IDH2","CBL","RUNX1","EZH2","ATM","KRAS","NRAS","GNAS","SETD2","ETV6"
)

DDR_GENES <- c(
  "TP53","ATM","ATR","CHEK2","BRCA1","BRCA2","PALB2","RAD51","RAD51C","RAD51D",
  "FANCA","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI","FANCJ","FANCL","FANCM",
  "BARD1","BRIP1","MRE11","NBN","RAD50","MSH2","MSH6","MLH1","PMS2","POLE","POLD1"
)

pt_gene_flags <- df_var %>%
  dplyr::distinct(patient_paper, Gene_simple) %>%
  dplyr::group_by(patient_paper) %>%
  dplyr::summarise(
    has_age_related = any(Gene_simple %in% AGE_RELATED_GENES),
    has_DDR         = any(Gene_simple %in% DDR_GENES),
    n_genes_total   = dplyr::n(),
    genes           = paste(sort(unique(Gene_simple)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(patient_paper)

cat("\n=== Per-patient gene flags (DDR vs age-related/CH) ===\n")
print(pt_gene_flags, n = Inf, width = Inf)

cat("\n=== Cross-tab: has_DDR vs has_age_related ===\n")
print(with(pt_gene_flags, table(has_DDR, has_age_related)))

readr::write_tsv(
  pt_gene_flags,
  file.path(OUT_DIR_SUMMARY, "patient_DDR_age_related_flags.tsv")
)

# ----------------------------
# 3) t-SNE highlighted per patient (paged grid)
# ----------------------------

# Make per-patient panels: each panel shows only that patient’s cells (colored by Cell_Type5)
# This avoids the unreadable “all patients in one legend” problem.
patients <- sort(unique(df_all$patient_paper))
PATIENTS_PER_PAGE <- 12  # 3x4 grid is usually readable
N_COL <- 4

# Ensure consistent axis limits across panels
x_rng <- range(df_all$optsne_1, na.rm = TRUE)
y_rng <- range(df_all$optsne_2, na.rm = TRUE)

make_patient_panel <- function(pt) {
  dat <- df_all %>% dplyr::filter(patient_paper == pt)

  ggplot(dat, aes(x = optsne_1, y = optsne_2, color = Cell_Type5)) +
    geom_point_ai_rast(size = 0.25, alpha = 0.85) +
    scale_colour_manual(values = pal_celltype, na.value = "grey80") +
    coord_cartesian(xlim = x_rng, ylim = y_rng) +
    ggtitle(pt) +
    theme_tsne() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

if (length(patients) > 0) {
  pages <- split(patients, ceiling(seq_along(patients) / PATIENTS_PER_PAGE))

  for (i in seq_along(pages)) {
    pts <- pages[[i]]
    plots <- lapply(pts, make_patient_panel)

    p_page <- patchwork::wrap_plots(plots, ncol = N_COL)

    out_png <- file.path(OUT_DIR_SUMMARY, sprintf("tsne_by_patient_page%02d.png", i))
    out_eps <- file.path(OUT_DIR_SUMMARY, sprintf("tsne_by_patient_page%02d.eps", i))

    ggsave(out_png, plot = p_page, width = 16, height = 12, units = "in", dpi = 300)
    ggsave(out_eps, plot = p_page, width = 16, height = 12, units = "in", device = cairo_ps)

    cat("[WROTE] ", out_png, "\n", sep = "")
  }
}

# ----------------------------
# 4) Heatmap: cell-type specificity per variant
#    (fraction of each variant’s cells that fall in each Cell_Type5)
# ----------------------------

# Choose whether to plot ALL variants or top N for readability
TOP_VARIANTS_TO_PLOT <- 80  # set to Inf if you truly want everything in one plot

variant_celltype <- df_var %>%
  dplyr::mutate(
    Variant = Base_Variant,
    Cell_Type5 = as.character(Cell_Type5)
  ) %>%
  dplyr::count(Variant, Gene_simple, Cell_Type5, name = "n_cells") %>%
  dplyr::group_by(Variant) %>%
  dplyr::mutate(
    n_total = sum(n_cells),
    frac_in_celltype = n_cells / pmax(n_total, 1L)
  ) %>%
  dplyr::ungroup()

# Order variants primarily by gene type (GENE_PRIORITY), then by abundance
variant_order_tbl <- variant_celltype %>%
  dplyr::distinct(Variant, Gene_simple, n_total) %>%
  dplyr::mutate(
    gene_rank = gene_rank(Gene_simple)
  ) %>%
  dplyr::arrange(gene_rank, dplyr::desc(n_total), Gene_simple, Variant)

variants_to_plot <- variant_order_tbl$Variant
if (is.finite(TOP_VARIANTS_TO_PLOT) && length(variants_to_plot) > TOP_VARIANTS_TO_PLOT) {
  variants_to_plot <- variants_to_plot[seq_len(TOP_VARIANTS_TO_PLOT)]
}

# keep Cell_Type5 ordering stable
ct_levels <- if (is.factor(df_all$Cell_Type5)) levels(df_all$Cell_Type5) else
  c("CD4 T-cell","CD8 T-cell","NK cell","Myeloid","Myeloid Progenitor","Unknown")

variant_celltype_plot_df <- variant_celltype %>%
  dplyr::filter(Variant %in% variants_to_plot) %>%
  dplyr::mutate(
    Variant    = factor(Variant, levels = variants_to_plot),
    Cell_Type5 = factor(Cell_Type5, levels = ct_levels)
  )

variants_to_plot <- variant_order_tbl$Variant
if (is.finite(TOP_VARIANTS_TO_PLOT) && length(variants_to_plot) > TOP_VARIANTS_TO_PLOT) {
  variants_to_plot <- variants_to_plot[seq_len(TOP_VARIANTS_TO_PLOT)]
}

variant_celltype_plot_df <- variant_celltype %>%
  dplyr::filter(Variant %in% variants_to_plot) %>%
  dplyr::mutate(
    Variant = factor(Variant, levels = variants_to_plot),
    # keep your Cell_Type5 ordering if already a factor in df_all; otherwise use a sensible order
    Cell_Type5 = factor(
      Cell_Type5,
      levels = levels(df_all$Cell_Type5) %||% c("CD4 T-cell","CD8 T-cell","NK cell","Myeloid","Myeloid Progenitor","Unknown")
    )
  )

p_variant_heatmap <- ggplot(
  variant_celltype_plot_df,
  aes(x = Cell_Type5, y = Variant, fill = frac_in_celltype)
) +
  geom_tile() +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    labels = scales::percent_format(accuracy = 1),
    name = "Fraction\nof variant"
  ) +
  labs(
    title = "Cell-type specificity per variant",
    x = "Cell type",
    y = "Variant"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

print(p_variant_heatmap)

ggsave(
  filename = file.path(OUT_DIR_SUMMARY, "variant_celltype_heatmap.png"),
  plot     = p_variant_heatmap,
  width    = 10,
  height   = max(6, 0.12 * length(variants_to_plot) + 3),
  units    = "in",
  dpi      = 300
)
ggsave(
  filename = file.path(OUT_DIR_SUMMARY, "variant_celltype_heatmap.eps"),
  plot     = p_variant_heatmap,
  width    = 10,
  height   = max(6, 0.12 * length(variants_to_plot) + 3),
  units    = "in",
  device   = cairo_ps
)

# Save the underlying heatmap table too (for supplements)
readr::write_tsv(
  variant_celltype %>% dplyr::arrange(dplyr::desc(n_total), Gene_simple, Variant, Cell_Type5),
  file.path(OUT_DIR_SUMMARY, "variant_celltype_fractions.tsv")
)

cat("\n[WROTE] variants_by_patient.tsv\n")
cat("[WROTE] patient_variant_list.tsv\n")
cat("[WROTE] patient_DDR_age_related_flags.tsv\n")
cat("[WROTE] variant_celltype_fractions.tsv\n")
cat("[WROTE] variant_celltype_heatmap.(png/eps)\n\n")


############################################################
## PART G – Single t-SNE colored by patient (all cells)
############################################################

OUT_DIR_SUMMARY <- out_dir
dir.create(OUT_DIR_SUMMARY, recursive = TRUE, showWarnings = FALSE)

# Build a stable palette for many patients
patients <- sort(unique(df_all$patient_paper))
n_pat <- length(patients)

# Use a qualitative palette that scales (falls back gracefully if many patients)
patient_cols <- scales::hue_pal()(n_pat)
names(patient_cols) <- patients

p_tsne_by_patient_all <- ggplot(
  df_all %>% dplyr::arrange(patient_paper),
  aes(x = optsne_1, y = optsne_2, color = patient_paper)
) +
  geom_point_ai_rast(size = 0.25, alpha = 0.85) +
  scale_colour_manual(values = patient_cols, name = "Patient") +
  xlab("t-SNE1") + ylab("t-SNE2") +
  ggtitle("t-SNE colored by patient") +
  theme_tsne() +
  guides(
    x = axis_trunc,
    y = axis_trunc,
    color = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.35, "cm"),
    legend.text = element_text(size = 9)
  )

print(p_tsne_by_patient_all)

ggsave(
  filename = file.path(OUT_DIR_SUMMARY, "tsne_colored_by_patient.png"),
  plot     = p_tsne_by_patient_all,
  width    = 10, height = 8, units = "in", dpi = 300
)
ggsave(
  filename = file.path(OUT_DIR_SUMMARY, "tsne_colored_by_patient.eps"),
  plot     = p_tsne_by_patient_all,
  width    = 10, height = 8, units = "in",
  device   = cairo_ps
)

############################################################
## PART Z – DIAGNOSTICS + FIXED PLOTS (variant heatmap + pies)
## Changes:
##  - Row labels on LEFT
##  - n_cells annotation on RIGHT only (drop n_patients)
##  - Black box around each gene block (row_split slice)
##  - Extra exports to trace TP53/PPM1D sources
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## ---------------------------
## Output directory
## ---------------------------
OUT_DIR_SUMMARY <- out_dir
dir.create(OUT_DIR_SUMMARY, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## Gene priority (ordering)
## ---------------------------
GENE_PRIORITY <- c(
  "DNMT3A","TET2","PPM1D","TP53","BCOR","ATM","EZH2","ASXL1",
  "JAK2","SF3B1","SRSF2","U2AF1","IDH1","IDH2","RUNX1"
)
gene_rank <- function(g) {
  r <- match(g, GENE_PRIORITY)
  ifelse(is.na(r), length(GENE_PRIORITY) + 999, r)
}

## ---------------------------
## Cell type order + colors
## ---------------------------
ct_levels <- c("CD4 T-cell","CD8 T-cell","NK cell","Myeloid","Myeloid Progenitor","Unknown")

pal_celltype <- c(
  "CD4 T-cell" = "#C6CDF7",
  "CD8 T-cell" = "#7294D4",
  "NK cell"    = NK_CELL_COLOR,
  "Myeloid"    = "#F98400",
  "Myeloid Progenitor"   = MYELOID_PROGENITOR_LIKE_COLOR,
  "Unknown"    = "grey80"
)

## ---------------------------
## Pretty gene palette
## ---------------------------
pal_gene <- c(
  pal_variant[c("DNMT3A","TET2","PPM1D","TP53","BCOR","ATM","EZH2","IDH1")],
  "ASXL1"  = "#6F7C91",
  "Other"  = "grey80"
)

## ---------------------------
## Robust weighting knobs
## ---------------------------
MIN_MUT_CELLS_PER_PATIENT_VARIANT  <- 20
MIN_MUT_CELLS_PER_PATIENT_CELLTYPE <- 25
MIN_MUT_CELLS_PER_PATIENT_GENE     <- 25
WEIGHT_CAP_Q <- 0.90
TOP_GENES_FOR_PIES <- 8

cat("\n============================================================\n")
cat("PART Z: DIAGNOSTICS + FIXED PLOTS\n")
cat("============================================================\n")

############################################################
## 0) Ensure df_var exists (mutated cells only)
############################################################
if (!exists("df_var")) {
  cat("\n[INFO] df_var not found; constructing df_var from df_all (mutated cells only)\n")
  df_var <- df_all %>%
    dplyr::filter(!is.na(Variant_Call), Variant_Call != "WT") %>%
    dplyr::mutate(
      Cell_Type5 = factor(as.character(Cell_Type5), levels = ct_levels),
      is_point_mut = stringr::str_detect(Variant_Call, "^[^:]+:p\\."),
      Base_Variant = dplyr::if_else(
        is_point_mut,
        stringr::str_extract(Variant_Call, "^[^:]+:p\\.[^:_]+"),
        Variant_Call
      ),
      Gene_simple = gene_from_variant(Base_Variant)
    ) %>%
    dplyr::select(patient_paper, Cell_Type5, Base_Variant, Gene_simple) %>%
    dplyr::filter(!is.na(patient_paper), !is.na(Cell_Type5), !is.na(Base_Variant), !is.na(Gene_simple))
} else {
  df_var <- df_var %>%
    mutate(
      Cell_Type5 = factor(as.character(Cell_Type5), levels = ct_levels),
      Gene_simple = as.character(Gene_simple),
      Base_Variant = as.character(Base_Variant)
    )
}

cat("\n[DIAG] df_var summary:\n")
cat("  Mutated cells (rows): ", nrow(df_var), "\n", sep = "")
cat("  Patients:            ", dplyr::n_distinct(df_var$patient_paper), "\n", sep = "")
cat("  Variants:            ", dplyr::n_distinct(df_var$Base_Variant), "\n", sep = "")
cat("  Genes:               ", dplyr::n_distinct(df_var$Gene_simple), "\n", sep = "")

## pooled counts by gene × cell type (NOT weighted; just raw composition)
diag_gene_ct_pooled <- df_var %>%
  count(Gene_simple, Cell_Type5, name = "n_cells") %>%
  group_by(Gene_simple) %>%
  mutate(frac_within_gene = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(gene_rank(Gene_simple), desc(n_cells))

cat("\n[DIAG] Pooled mutated-cell counts by gene × cell type:\n")
print(diag_gene_ct_pooled, n = Inf)

readr::write_tsv(
  diag_gene_ct_pooled,
  file.path(OUT_DIR_SUMMARY, "DIAG_pooled_gene_by_celltype.tsv")
)

############################################################
## 1) Diagnostics: where TP53/PPM1D are coming from
############################################################
diag_key_genes <- c("TP53","PPM1D","DNMT3A","TET2")

diag_key_by_patient <- df_var %>%
  filter(Gene_simple %in% diag_key_genes) %>%
  count(Gene_simple, patient_paper, Cell_Type5, Base_Variant, name = "n_cells") %>%
  arrange(gene_rank(Gene_simple), patient_paper, desc(n_cells))

cat("\n[DIAG] Key genes: per patient × cell type × variant (first 200 rows):\n")
print(diag_key_by_patient, n = 200)

readr::write_tsv(
  diag_key_by_patient,
  file.path(OUT_DIR_SUMMARY, "DIAG_key_genes_patient_celltype_variant.tsv")
)

diag_tp53_ct <- df_var %>%
  filter(Gene_simple == "TP53") %>%
  count(patient_paper, Cell_Type5, name = "n_cells") %>%
  arrange(desc(n_cells))

cat("\n[DIAG] TP53: mutated-cell counts per patient × cell type:\n")
print(diag_tp53_ct, n = Inf)

readr::write_tsv(
  diag_tp53_ct,
  file.path(OUT_DIR_SUMMARY, "DIAG_TP53_patient_celltype_counts.tsv")
)

diag_tp53_variants <- df_var %>%
  filter(Gene_simple == "TP53") %>%
  count(Base_Variant, Cell_Type5, name = "n_cells") %>%
  arrange(desc(n_cells))

cat("\n[DIAG] TP53: mutated-cell counts per TP53 variant × cell type:\n")
print(diag_tp53_variants, n = Inf)

readr::write_tsv(
  diag_tp53_variants,
  file.path(OUT_DIR_SUMMARY, "DIAG_TP53_variant_celltype_counts.tsv")
)

############################################################
## 2) FIXED VARIANT × CELLTYPE heatmap (patient-weighted; zero-completed)
############################################################

## 2a) Stable variant → gene mapping
variant_gene_map <- df_var %>%
  distinct(Base_Variant, Gene_simple) %>%
  group_by(Base_Variant) %>%
  slice_head(n = 1) %>%
  ungroup()

## 2b) Patient × variant × celltype counts
variant_ct_counts <- df_var %>%
  count(patient_paper, Base_Variant, Cell_Type5, name = "n_cells") %>%
  left_join(variant_gene_map, by = "Base_Variant") %>%
  mutate(Cell_Type5 = factor(as.character(Cell_Type5), levels = ct_levels))

## 2c) COMPLETE missing celltypes as zeros per patient×variant (CRITICAL FIX)
variant_ct_complete <- variant_ct_counts %>%
  group_by(patient_paper, Base_Variant, Gene_simple) %>%
  tidyr::complete(
    Cell_Type5 = factor(ct_levels, levels = ct_levels),
    fill = list(n_cells = 0L)
  ) %>%
  ungroup()

## export full patient×variant×celltype table (small: ~patients×variants×celltypes)
readr::write_tsv(
  variant_ct_complete,
  file.path(OUT_DIR_SUMMARY, "DIAG_patient_variant_celltype_counts_ZERO_COMPLETED.tsv")
)

## 2d) Per patient×variant totals + filters
variant_pt_totals <- variant_ct_complete %>%
  group_by(patient_paper, Base_Variant, Gene_simple) %>%
  summarise(n_variant_cells = sum(n_cells), .groups = "drop")

cat("\n[DIAG] Patient×variant totals (before filtering):\n")
print(summary(variant_pt_totals$n_variant_cells))

readr::write_tsv(
  variant_pt_totals,
  file.path(OUT_DIR_SUMMARY, "DIAG_patient_variant_totals.tsv")
)

variant_ct_complete2 <- variant_ct_complete %>%
  left_join(variant_pt_totals, by = c("patient_paper","Base_Variant","Gene_simple")) %>%
  filter(n_variant_cells >= MIN_MUT_CELLS_PER_PATIENT_VARIANT) %>%
  mutate(frac_pt = n_cells / pmax(n_variant_cells, 1L))

cat("\n[DIAG] Patient×variant totals (after filtering):\n")
print(summary((variant_ct_complete2 %>% distinct(patient_paper, Base_Variant, n_variant_cells))$n_variant_cells))

## Weight cap
w_cap_variant <- stats::quantile(
  (variant_ct_complete2 %>% distinct(patient_paper, Base_Variant, n_variant_cells))$n_variant_cells,
  probs = WEIGHT_CAP_Q, na.rm = TRUE
)

variant_ct_complete2 <- variant_ct_complete2 %>%
  mutate(w = pmin(n_variant_cells, w_cap_variant))

## 2e) Patient-weighted mean fractions per variant×celltype + counts used
variant_celltype_avg <- variant_ct_complete2 %>%
  group_by(Base_Variant, Gene_simple, Cell_Type5) %>%
  summarise(
    frac = sum(w * frac_pt, na.rm = TRUE) / pmax(sum(w, na.rm = TRUE), 1e-9),
    n_cells_used = sum(n_cells),
    n_patients_used = n_distinct(patient_paper),
    .groups = "drop"
  )

## 2f) Sanity: fractions should sum to ~1 per variant
variant_sumcheck <- variant_celltype_avg %>%
  group_by(Base_Variant, Gene_simple) %>%
  summarise(
    sum_frac = sum(frac, na.rm = TRUE),
    n_cells_used = sum(n_cells_used),
    n_patients_used = max(n_patients_used),
    .groups = "drop"
  ) %>%
  mutate(delta = abs(sum_frac - 1)) %>%
  arrange(desc(delta))

cat("\n[DIAG] Sum-to-1 check (largest deviations):\n")
print(head(variant_sumcheck, 25))

readr::write_tsv(
  variant_sumcheck,
  file.path(OUT_DIR_SUMMARY, "DIAG_variant_sum_to_1_check.tsv")
)

############################################################
## 2g–2j) REBUILD heatmap matrix + splits (fix orphan gene gaps)
############################################################

## 2g) Build ordered variant list with gene splits
variant_meta <- variant_celltype_avg %>%
  distinct(Base_Variant, Gene_simple) %>%
  mutate(gene_rank = gene_rank(Gene_simple)) %>%
  arrange(gene_rank, Gene_simple, Base_Variant)

variant_levels <- variant_meta$Base_Variant
gene_levels_present <- unique(variant_meta$Gene_simple[order(gene_rank(variant_meta$Gene_simple))])

## Named map: Base_Variant -> Gene_simple (THIS is the key to perfect alignment)
gene_by_variant <- variant_meta %>%
  distinct(Base_Variant, Gene_simple) %>%
  tibble::deframe()

## 2h) Zero-complete (already done) but ensure factor levels are correct
variant_celltype_wide <- variant_celltype_avg %>%
  mutate(
    Base_Variant = factor(Base_Variant, levels = variant_levels),
    Cell_Type5   = factor(as.character(Cell_Type5), levels = ct_levels)
  ) %>%
  select(-Gene_simple) %>%
  tidyr::complete(
    Base_Variant, Cell_Type5,
    fill = list(frac = 0, n_patients_used = 0, n_cells_used = 0)
  ) %>%
  mutate(Base_Variant = as.character(Base_Variant)) %>%
  left_join(variant_meta %>% select(Base_Variant, Gene_simple), by = "Base_Variant") %>%
  select(Base_Variant, Gene_simple, Cell_Type5, frac, n_cells_used)

## 2h2) REBUILD heatmap_mat from scratch (avoid stale heatmap_mat)
heatmap_mat <- variant_celltype_wide %>%
  mutate(
    Base_Variant = factor(Base_Variant, levels = variant_levels),
    Cell_Type5   = factor(as.character(Cell_Type5), levels = ct_levels)
  ) %>%
  arrange(Base_Variant, Cell_Type5) %>%
  select(Base_Variant, Cell_Type5, frac) %>%
  tidyr::pivot_wider(names_from = Cell_Type5, values_from = frac, values_fill = 0) %>%
  arrange(Base_Variant) %>%
  tibble::column_to_rownames("Base_Variant") %>%
  as.matrix()

heatmap_mat <- heatmap_mat[variant_levels, ct_levels, drop = FALSE]
storage.mode(heatmap_mat) <- "numeric"

## 2h3) row_split MUST be computed from rownames(heatmap_mat) (prevents orphan slices)
row_split <- factor(
  gene_by_variant[rownames(heatmap_mat)],
  levels = gene_levels_present
)

## --- Diagnostics to prove splits are correct ---
cat("\n[DIAG] Heatmap row_split counts by gene:\n")
print(table(row_split, useNA = "ifany"))

diag_orphan <- tibble(
  Base_Variant = rownames(heatmap_mat),
  Gene_split   = as.character(row_split)
) %>%
  mutate(Gene_from_map = gene_by_variant[Base_Variant]) %>%
  filter(is.na(Gene_split) | Gene_split != Gene_from_map)

if (nrow(diag_orphan) > 0) {
  cat("\n[WARN] Orphan / mismatch rows detected (these would create weird gaps):\n")
  print(diag_orphan, n = Inf)
  readr::write_tsv(diag_orphan, file.path(OUT_DIR_SUMMARY, "DIAG_orphan_rows_heatmap.tsv"))
}

## 2i) Color function (same as before)
col_fun <- circlize::colorRamp2(
  c(0, 0.05, 0.15, 0.50, 1.00),
  c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#b2182b")
)

## 2j) Row annotation: n_cells only (RIGHT)
row_anno_df <- variant_sumcheck %>%
  right_join(variant_meta, by = c("Base_Variant","Gene_simple")) %>%
  mutate(Base_Variant = factor(Base_Variant, levels = variant_levels)) %>%
  arrange(Base_Variant)

row_anno <- rowAnnotation(
  `n_cells` = anno_barplot(row_anno_df$n_cells_used, border = FALSE, gp = gpar(fill = "grey35")),
  width = unit(1.7, "cm"),
  annotation_name_gp = gpar(fontsize = 10)
)

## Heatmap: labels LEFT, gene splits, slice gaps
ht <- Heatmap(
  heatmap_mat,
  name = "Fraction",
  col = col_fun,
  row_split = row_split,
  row_gap = unit(2.5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 45,
  rect_gp = gpar(col = "grey92", lwd = 0.5),
  column_title = "Cell type specificity of variants (patient-weighted; zero-completed)",
  heatmap_legend_param = list(
    at = c(0, 0.5, 1),
    labels = c("0%", "50%", "100%")
  )
)

## Save heatmap with black box around each gene block
save_heatmap_with_gene_boxes <- function(filename_pdf, filename_png,
                                         ht_obj, anno_obj, n_slices,
                                         width_pdf = 8.5, height_pdf = 7.5,
                                         png_w = 2400, png_h = 1800, png_res = 300) {

  pdf(filename_pdf, width = width_pdf, height = height_pdf)
  draw(ht_obj + anno_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
  for (i in seq_len(n_slices)) {
    decorate_heatmap_body("Fraction", slice = i, {
      grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1.6))
    })
  }
  dev.off()

  png(filename_png, width = png_w, height = png_h, res = png_res)
  draw(ht_obj + anno_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
  for (i in seq_len(n_slices)) {
    decorate_heatmap_body("Fraction", slice = i, {
      grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1.6))
    })
  }
  dev.off()
}

save_heatmap_with_gene_boxes(
  file.path(OUT_DIR_SUMMARY, "variant_celltype_heatmap_FIXED_BOXED.pdf"),
  file.path(OUT_DIR_SUMMARY, "variant_celltype_heatmap_FIXED_BOXED.png"),
  ht, row_anno,
  n_slices = length(levels(row_split))
)

## Export heatmap values (no n_patients)
heatmap_export <- variant_celltype_wide %>%
  select(Base_Variant, Gene_simple, Cell_Type5, frac, n_cells_used) %>%
  arrange(gene_rank(Gene_simple), Gene_simple, Base_Variant, Cell_Type5)

readr::write_tsv(
  heatmap_export,
  file.path(OUT_DIR_SUMMARY, "EXPORT_variant_celltype_heatmap_values.tsv")
)

############################################################
## 3) Pie set A – Within each cell type: which genes?
############################################################

gene_abundance <- df_var %>%
  count(Gene_simple, name = "n_cells") %>%
  arrange(gene_rank(Gene_simple), desc(n_cells))

top_genes <- gene_abundance$Gene_simple
if (length(top_genes) > TOP_GENES_FOR_PIES) top_genes <- top_genes[seq_len(TOP_GENES_FOR_PIES)]

df_var_pies <- df_var %>%
  mutate(
    Gene_pie = if_else(Gene_simple %in% top_genes, Gene_simple, "Other"),
    Cell_Type5 = factor(as.character(Cell_Type5), levels = ct_levels)
  )

for (g in unique(df_var_pies$Gene_pie)) if (!g %in% names(pal_gene)) pal_gene[g] <- "grey80"

celltype_gene_counts <- df_var_pies %>%
  count(patient_paper, Cell_Type5, Gene_pie, name = "n_cells") %>%
  group_by(patient_paper, Cell_Type5) %>%
  tidyr::complete(
    Gene_pie = unique(c(top_genes, "Other")),
    fill = list(n_cells = 0L)
  ) %>%
  ungroup()

celltype_gene_totals <- celltype_gene_counts %>%
  group_by(patient_paper, Cell_Type5) %>%
  summarise(n_total = sum(n_cells), .groups = "drop")

celltype_gene_pt <- celltype_gene_counts %>%
  left_join(celltype_gene_totals, by = c("patient_paper","Cell_Type5")) %>%
  filter(n_total >= MIN_MUT_CELLS_PER_PATIENT_CELLTYPE) %>%
  mutate(frac_pt = n_cells / pmax(n_total, 1L))

w_cap_ct <- stats::quantile(
  (celltype_gene_pt %>% distinct(patient_paper, Cell_Type5, n_total))$n_total,
  probs = WEIGHT_CAP_Q, na.rm = TRUE
)

celltype_gene_avg <- celltype_gene_pt %>%
  mutate(w = pmin(n_total, w_cap_ct)) %>%
  group_by(Cell_Type5, Gene_pie) %>%
  summarise(
    frac = sum(w * frac_pt, na.rm = TRUE) / pmax(sum(w, na.rm = TRUE), 1e-9),
    n_cells_used = sum(n_cells),
    n_patients_used = n_distinct(patient_paper),
    .groups = "drop"
  ) %>%
  group_by(Cell_Type5) %>%
  mutate(frac = frac / pmax(sum(frac), 1e-9)) %>%
  ungroup()

cat("\n[DIAG] Pie A (cell type -> genes):\n")
print(celltype_gene_avg %>% arrange(Cell_Type5, desc(frac)), n = Inf)

readr::write_tsv(
  celltype_gene_avg,
  file.path(OUT_DIR_SUMMARY, "EXPORT_pieA_celltype_to_genes_values.tsv")
)

p_pies_by_celltype <- ggplot(celltype_gene_avg, aes(x = "", y = frac, fill = Gene_pie)) +
  geom_col(width = 1, color = "white", linewidth = 0.35) +
  coord_polar(theta = "y") +
  facet_wrap(~ Cell_Type5, nrow = 1) +
  scale_fill_manual(values = pal_gene, name = "Gene") +
  theme_void(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Gene mutation composition within each cell type (patient-weighted)")

ggsave(
  file.path(OUT_DIR_SUMMARY, "pies_gene_by_celltype_FIXED.pdf"),
  p_pies_by_celltype,
  width = 14, height = 3.2, units = "in", dpi = 300
)

############################################################
## 4) Pie set B – Within each gene: which cell types?
############################################################

gene_celltype_counts <- df_var_pies %>%
  filter(Gene_pie != "Other") %>%
  count(patient_paper, Gene_pie, Cell_Type5, name = "n_cells") %>%
  group_by(patient_paper, Gene_pie) %>%
  tidyr::complete(
    Cell_Type5 = factor(ct_levels, levels = ct_levels),
    fill = list(n_cells = 0L)
  ) %>%
  ungroup()

gene_celltype_totals <- gene_celltype_counts %>%
  group_by(patient_paper, Gene_pie) %>%
  summarise(n_total = sum(n_cells), .groups = "drop")

gene_celltype_pt <- gene_celltype_counts %>%
  left_join(gene_celltype_totals, by = c("patient_paper","Gene_pie")) %>%
  filter(n_total >= MIN_MUT_CELLS_PER_PATIENT_GENE) %>%
  mutate(frac_pt = n_cells / pmax(n_total, 1L))

w_cap_g <- stats::quantile(
  (gene_celltype_pt %>% distinct(patient_paper, Gene_pie, n_total))$n_total,
  probs = WEIGHT_CAP_Q, na.rm = TRUE
)

gene_celltype_avg <- gene_celltype_pt %>%
  mutate(w = pmin(n_total, w_cap_g)) %>%
  group_by(Gene_pie, Cell_Type5) %>%
  summarise(
    frac = sum(w * frac_pt, na.rm = TRUE) / pmax(sum(w, na.rm = TRUE), 1e-9),
    n_cells_used = sum(n_cells),
    n_patients_used = n_distinct(patient_paper),
    .groups = "drop"
  ) %>%
  group_by(Gene_pie) %>%
  mutate(frac = frac / pmax(sum(frac), 1e-9)) %>%
  ungroup()

## enforce gene order
gene_levels_present2 <- unique(as.character(gene_celltype_avg$Gene_pie))
gene_levels_present2 <- gene_levels_present2[order(gene_rank(gene_levels_present2))]
gene_celltype_avg <- gene_celltype_avg %>%
  mutate(Gene_pie = factor(Gene_pie, levels = gene_levels_present2))

cat("\n[DIAG] Pie B (gene -> cell types):\n")
print(gene_celltype_avg %>% arrange(Gene_pie, desc(frac)), n = Inf)

readr::write_tsv(
  gene_celltype_avg,
  file.path(OUT_DIR_SUMMARY, "EXPORT_pieB_gene_to_celltypes_values.tsv")
)

## Extra: per-patient TP53 composition (this is the “where is TP53 coming from” smoking gun)
tp53_patient_breakdown <- gene_celltype_pt %>%
  filter(Gene_pie == "TP53") %>%
  select(patient_paper, Gene_pie, Cell_Type5, n_cells, n_total, frac_pt) %>%
  arrange(desc(n_total), desc(n_cells))

readr::write_tsv(
  tp53_patient_breakdown,
  file.path(OUT_DIR_SUMMARY, "DIAG_TP53_per_patient_celltype_frac_pt.tsv")
)

p_pies_by_gene <- ggplot(gene_celltype_avg, aes(x = "", y = frac, fill = Cell_Type5)) +
  geom_col(width = 1, color = "white", linewidth = 0.35) +
  coord_polar(theta = "y") +
  facet_wrap(~ Gene_pie, nrow = 1) +
  scale_fill_manual(values = pal_celltype, name = "Cell type") +
  theme_void(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Cell-type composition within each gene mutation (patient-weighted)")

ggsave(
  file.path(OUT_DIR_SUMMARY, "pies_celltype_by_gene_FIXED.pdf"),
  p_pies_by_gene,
  width = 14, height = 3.2, units = "in", dpi = 300
)

cat("\n[DONE] Wrote outputs to:\n")
cat("  ", OUT_DIR_SUMMARY, "\n", sep = "")
cat("  - variant_celltype_heatmap_FIXED_BOXED.(pdf/png)\n")
cat("  - pies_gene_by_celltype_FIXED.png\n")
cat("  - pies_celltype_by_gene_FIXED.png\n")
cat("  - DIAG_*.tsv and EXPORT_*.tsv for tracing sources\n")


############ lymphoid CH clonal tracing #########
############################################################
## PART B2 – Clonal tracing of LYMPHOID CH (df3_multi)
## Goal: identical to myeloid clonal tracing plot, but
##       denominator = LYMPHOID compartment (PBMC CD4/CD8/NK
##       collapsed to "Lymphoid" to match BM "Lymphoid").
############################################################

## --- 0) Define a 2-level lineage that harmonizes PBMC vs BM ---
## PBMC: CD4_Tcell / CD8_Tcell / NK_cell -> Lymphoid
## BM:   Lymphoid already -> Lymphoid
## Myeloid: Myeloid + CD34_progenitor -> Myeloid
df3_multi <- df3_multi %>%
  mutate(
    lineage_group2 = case_when(
      OmiqFilter %in% c("CD4_Tcell","CD4_T-cell","CD4_T","CD4 T-cell",
                        "CD8_Tcell","CD8_T-cell","CD8_T","CD8 T-cell",
                        "NK_cell","NK_Cell","NK cell","NK Cell","NKcell","NK",
                        "Lymphoid") ~ "Lymphoid",
      OmiqFilter %in% c("Myeloid","Myeloid_NOS","CD34_progenitor","CD34+Progenitor",
                        "CD71_erythroid","CD71+Erythroid") ~ "Myeloid",
      TRUE ~ NA_character_
    ),
    lineage_group2 = factor(lineage_group2, levels = c("Lymphoid","Myeloid")),

    is_lymphoid = lineage_group2 == "Lymphoid",
    is_myeloid2 = lineage_group2 == "Myeloid",

    ## robust CH flag
    is_CH = !is.na(Variant_Call) & Variant_Call != "WT"
  )

cat("\n=== [LYMPHOID] OmiqFilter -> lineage_group2 mapping ===\n")
df3_multi %>% count(OmiqFilter, lineage_group2) %>% arrange(OmiqFilter, desc(n)) %>% print(n = Inf)

## --- 1) Lymphoid totals & CH+ fraction over time (per patient) ---
lymphoid_totals <- df3_multi %>%
  filter(is_lymphoid) %>%
  group_by(Patient_Paper, Day) %>%
  summarise(
    n_lymphoid       = n(),
    n_lymphoid_CH    = sum(is_CH),
    frac_CH_lymphoid = n_lymphoid_CH / pmax(n_lymphoid, 1L),
    .groups          = "drop"
  )

## Ensure: if lymphoid exists at a day but CH=0 -> fraction 0 (already)
## If lymphoid absent at a day -> it won't appear (which is correct; NA if desired)

cat("\n=== [LYMPHOID] totals & CH+ fractions (multi-timepoint) ===\n")
lymphoid_totals %>% arrange(Patient_Paper, Day) %>% print(n = Inf)

## Plot: CH+ fraction in lymphoid over time (same style as myeloid)
all_days_ly <- sort(unique(df3_multi$Day))
day_labels_ly <- paste0("D", all_days_ly)

p_lymphoid_CH_time <- ggplot(
  lymphoid_totals,
  aes(x = Day, y = frac_CH_lymphoid, group = Patient_Paper)
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = all_days_ly, labels = day_labels_ly) +
  labs(
    title = "Fraction of lymphoid cells that are CH+ over time (PBMC CD4/CD8/NK collapsed to Lymphoid)",
    x = "Timepoint (Day)",
    y = "Fraction of lymphoid cells that are CH+"
  ) +
  theme_pub

print(p_lymphoid_CH_time)

ggsave(
  out_file("lymphoid_cells_with_CH_fraction_over_time.eps"),
  plot  = p_lymphoid_CH_time,
  width = 10, height = 5, units = "in",
  device = cairo_ps
)

## --- 2) Founder-level LYMPHOID clone trajectories (identical pattern to myeloid) ---

## observed founder clone counts in lymphoid
lymphoid_clone_obs <- df3_multi %>%
  filter(is_lymphoid, is_CH) %>%
  group_by(Patient_Paper, Day, Base_Variant, Gene) %>%
  summarise(
    n_cells_clone = n(),
    .groups       = "drop"
  )

## clone list per patient
clone_meta_ly <- lymphoid_clone_obs %>%
  distinct(Patient_Paper, Base_Variant, Gene)

## full grid of (Patient_Paper, Day, Base_Variant, Gene) with lymphoid totals
lymphoid_clone_time <- lymphoid_totals %>%
  select(Patient_Paper, Day, n_lymphoid) %>%
  left_join(
    clone_meta_ly,
    by = "Patient_Paper",
    relationship = "many-to-many"
  ) %>%
  left_join(
    lymphoid_clone_obs,
    by = c("Patient_Paper","Day","Base_Variant","Gene")
  ) %>%
  mutate(
    n_cells_clone    = coalesce(n_cells_clone, 0L),
    frac_of_lymphoid = n_cells_clone / pmax(n_lymphoid, 1L)
  )

cat("\n=== [LYMPHOID] Founder-level clones over time (first 50 rows) ===\n")
lymphoid_clone_time %>%
  arrange(Patient_Paper, Base_Variant, Day) %>%
  print(n = 50)

## keep founders that ever reach >= 1% of lymphoid
clone_keep_ly <- lymphoid_clone_time %>%
  group_by(Patient_Paper, Base_Variant, Gene) %>%
  summarise(
    max_frac = max(frac_of_lymphoid, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  filter(max_frac >= 0.01)

cat("\n=== [LYMPHOID] Founder clones reaching ≥1% of lymphoid ===\n")
clone_keep_ly %>%
  arrange(Patient_Paper, desc(max_frac)) %>%
  print(n = Inf)

lymphoid_clone_time_filt <- lymphoid_clone_time %>%
  inner_join(
    clone_keep_ly %>% select(Patient_Paper, Base_Variant),
    by = c("Patient_Paper","Base_Variant")
  )

## optional: mirror your prior exclusion
lymphoid_clone_time_filt <- lymphoid_clone_time_filt %>%
  filter(Patient_Paper != study_setting("study_value_001"))

cat("\n[COLOR CHECK] Genes in lymphoid clone plot not in pal_variant:\n")
print(setdiff(unique(na.omit(lymphoid_clone_time_filt$Gene)), names(pal_variant)))

## Plot: founder clone trajectories in LYMPHOID, colored by CHIP gene palette
p_lymphoid_clone_traj <- ggplot(
  lymphoid_clone_time_filt,
  aes(
    x     = Day,
    y     = frac_of_lymphoid,
    group = Base_Variant,
    color = Gene
  )
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(labels = function(x) paste0("D", x)) +
  scale_color_manual(
    values = pal_variant,
    breaks = c("DNMT3A","TET2","TP53","PPM1D","BCOR","ATM","EZH2","IDH1"),
    name   = "Gene",
    na.value = "grey50"
  ) +
  labs(
    title = "CH gene mutations in lymphoid cells over time measured by scDNA",
    x = "Timepoint (Day)",
    y = "Fraction of lymphoid compartment"
  ) +
  theme_pub +
  theme(legend.position = "right")

print(p_lymphoid_clone_traj)

ggsave(
  out_file("lymphoid_cells_with_CH_together.eps"),
  plot  = p_lymphoid_clone_traj,
  width = 10, height = 5, units = "in",
  device = cairo_ps
)

## --- 3) Subclone trajectories – point mutations only (LYMPHOID) ---

lymphoid_subclone_obs <- df3_multi %>%
  filter(is_lymphoid, is_CH, is_point_mut) %>%
  group_by(Patient_Paper, Day, Clone_Label, Base_Variant, Gene) %>%
  summarise(
    n_cells_clone = n(),
    .groups       = "drop"
  )

subclone_meta_ly <- lymphoid_subclone_obs %>%
  distinct(Patient_Paper, Clone_Label, Base_Variant, Gene)

subclone_time_ly <- lymphoid_totals %>%
  select(Patient_Paper, Day, n_lymphoid) %>%
  left_join(
    subclone_meta_ly,
    by = "Patient_Paper",
    relationship = "many-to-many"
  ) %>%
  left_join(
    lymphoid_subclone_obs,
    by = c("Patient_Paper","Day","Clone_Label","Base_Variant","Gene")
  ) %>%
  mutate(
    n_cells_clone    = coalesce(n_cells_clone, 0L),
    frac_of_lymphoid = n_cells_clone / pmax(n_lymphoid, 1L),
    founder_status   = if_else(Clone_Label == Base_Variant,
                               "Founder", "Non-founder")
  )

## keep subclones that ever reach >=1% of lymphoid
subclone_keep_ly <- subclone_time_ly %>%
  group_by(Patient_Paper, Clone_Label, Base_Variant, Gene) %>%
  summarise(
    max_frac = max(frac_of_lymphoid, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  filter(max_frac >= 0.01)

cat("\n=== [LYMPHOID] Subclones reaching ≥1% of lymphoid (point mutations) ===\n")
subclone_keep_ly %>%
  arrange(Patient_Paper, desc(max_frac)) %>%
  print(n = Inf)

subclone_time_filt_ly <- subclone_time_ly %>%
  inner_join(
    subclone_keep_ly %>% select(Patient_Paper, Clone_Label),
    by = c("Patient_Paper","Clone_Label")
  ) %>%
  filter(Patient_Paper != study_setting("study_value_001"))

p_lymphoid_subclone_traj <- ggplot(
  subclone_time_filt_ly,
  aes(
    x        = Day,
    y        = frac_of_lymphoid,
    group    = Clone_Label,
    color    = Base_Variant,
    linetype = founder_status
  )
) +
  geom_line(size = 0.8) +
  geom_point(size = 2.0) +
  facet_wrap(~ Patient_Paper, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = all_days_ly, labels = day_labels_ly) +
  scale_color_manual(
    values = variant_gene_palette(
      subclone_time_filt_ly$Base_Variant,
      subclone_time_filt_ly$Gene
    ),
    name = "Founder (Base_Variant)"
  ) +
  scale_linetype_manual(
    values = c("Founder" = "solid", "Non-founder" = "dashed"),
    name   = "Clone type"
  ) +
  labs(
    title = "Lymphoid point-mutation founders vs subclones over time",
    x = "Timepoint (Day)",
    y = "Fraction of lymphoid compartment"
  ) +
  theme_pub +
  theme(legend.position = "right")

print(p_lymphoid_subclone_traj)


## Save session information for reproducibility
writeLines(
  capture.output(sessionInfo()),
  con = out_file("scDNA_mutation_traces_sessionInfo.txt")
)

cat("\nSaved final scDNA mutation-tracing outputs to:\n", out_dir, "\n", sep = "")


############################################################
## Briefing figure:
## Frequent detection of CAR+ CH mutant cells in expanding peripheral blood
## - Background: low-alpha cell-type context
## - Foreground: CAR+ CH mutant cells colored by CH mutation
## - Myeloid progenitors are collapsed into Myeloid
############################################################

## Optional: define output helper if this block is run standalone
if (!exists("out_dir")) {
  out_dir <- study_path("out_dir")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("out_file")) {
  out_file <- function(...) file.path(out_dir, ...)
}

## Optional: italic gene labels if not already defined
if (!exists("gene_labeller_expr")) {
  gene_labeller_expr <- function(x) {
    lapply(x, function(g) parse(text = paste0("italic('", g, "')")))
  }
}

## Muted cell-type context palette
## Intentionally subdued so CH mutation colors remain the visual focus
pal_cell_context <- c(
  "CD4 T-cell" = "#C5BAD8",  # soft lavender
  "CD8 T-cell" = "#D8B0B0",  # dusty rose
  "NK cell"    = "#AEBEAA",  # muted sage
  "Myeloid"    = "#CDBEAA"   # warm stone/tan
)

## Build background context data
## Collapse Myeloid Progenitor into Myeloid for this figure
df_context <- df_all %>%
  dplyr::mutate(
    Cell_Type_Context = dplyr::case_when(
      Cell_Type5 == "CD4 T-cell" ~ "CD4 T-cell",
      Cell_Type5 == "CD8 T-cell" ~ "CD8 T-cell",
      Cell_Type5 == "NK cell" ~ "NK cell",
      Cell_Type5 %in% c("Myeloid", "Myeloid Progenitor") ~ "Myeloid",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Cell_Type_Context)) %>%
  dplyr::mutate(
    Cell_Type_Context = factor(
      Cell_Type_Context,
      levels = c("CD4 T-cell", "CD8 T-cell", "NK cell", "Myeloid")
    )
  )

## Foreground CAR+ CH mutant cells
df_car_ch <- df_all %>%
  dplyr::filter(is_CAR_CH, !is.na(CH_gene_for_plot)) %>%
  dplyr::mutate(
    CH_gene_for_plot = as.character(CH_gene_for_plot)
  )

## Order mutation legend by mutations actually present in the plotted cells
car_ch_genes_briefing <- intersect(
  names(pal_CH_category)[names(pal_CH_category) != "WT"],
  unique(df_car_ch$CH_gene_for_plot)
)

## Main merged briefing plot
p_tsne_CAR_CH_cellcontext <- ggplot() +
  ## Very light all-cell embedding background
  geom_point_ai_rast(
    data = df_all,
    aes(x = optsne_1, y = optsne_2),
    color = "grey92",
    size = 0.20,
    alpha = 0.10
  ) +

  ## Low-alpha cell-type context layer
  geom_point_ai_rast(
    data = df_context,
    aes(x = optsne_1, y = optsne_2, color = Cell_Type_Context),
    size = 0.34,
    alpha = 0.18
  ) +

  ## Foreground CAR+ CH+ cells
  geom_point_ai_rast(
    data = df_car_ch,
    aes(x = optsne_1, y = optsne_2, fill = CH_gene_for_plot),
    shape = 21,
    color = "black",
    stroke = 0.20,
    size = 1.30,
    alpha = 0.98
  ) +

  scale_color_manual(
    values = pal_cell_context,
    name   = "Cell-type context",
    drop   = FALSE
  ) +

  scale_fill_manual(
    values   = pal_CH_category[names(pal_CH_category) != "WT"],
    breaks   = car_ch_genes_briefing,
    labels   = gene_labeller_expr(car_ch_genes_briefing),
    name     = "CAR+ CH mutation",
    na.value = "grey50"
  ) +

  xlab("t-SNE1") +
  ylab("t-SNE2") +
  ggtitle("Frequent Detection of CAR+ CH Mutant Cells in Expanding Peripheral Blood") +
  theme_tsne() +
  guides(
    x = axis_trunc,
    y = axis_trunc,
    color = guide_legend(
      override.aes = list(size = 4, alpha = 0.7)
    ),
    fill = guide_legend(
      override.aes = list(size = 4, alpha = 1, color = "black")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 15, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.key.size = grid::unit(0.55, "cm"),
    plot.title      = element_text(size = 18, face = "bold")
  )

print(p_tsne_CAR_CH_cellcontext)

## Save outputs
ggsave(
  out_file("CAR_CH_mutant_tsne_celltype_context_briefing.pdf"),
  plot = p_tsne_CAR_CH_cellcontext,
  width = 10.0, height = 7.4, units = "in",
  device = cairo_pdf
)

ggsave(
  out_file("CAR_CH_mutant_tsne_celltype_context_briefing.eps"),
  plot = p_tsne_CAR_CH_cellcontext,
  width = 10.0, height = 7.4, units = "in",
  device = cairo_ps
)

ggsave(
  out_file("CAR_CH_mutant_tsne_celltype_context_briefing.png"),
  plot = p_tsne_CAR_CH_cellcontext,
  width = 10.0, height = 7.4, units = "in",
  dpi = 400
)
