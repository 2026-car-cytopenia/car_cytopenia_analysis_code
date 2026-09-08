#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("after_omiq")

## This is a script file designed for analyzing scDNA seq data.
# %>%
library(ggplot2)  # graphics
library(survival)
library(patchwork)
library(ggprism)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(tableone)

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


## Optional execution log: console output is still printed and also written to this file.
log_con <- file(out_file("scDNA_after_omiq_session_output.txt"), open = "wt")
sink(log_con, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

################################ Define color palette

wes_palettes <- list(
  BottleRocket1 = c("#A42820", "#5F5647", "#9B110E", "#3F5151", "#4E2A1E", "#550307", "#0C1707"),
  BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"),
  Rushmore1 = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" ,"#F2300F"),
  Rushmore = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" ,"#F2300F"),
  Royal1 = c("#899DA4", "#C93312", "#FAEFD1", "#DC863B"),
  Royal2 = c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089"),
  Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"),
  Darjeeling1 = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"),
  Darjeeling2 = c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000"),
  Chevalier1 = c("#446455", "#FDD262", "#D3DDDC", "#C7B19C"),
  FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"),
  Moonrise1 = c("#F3DF6C", "#CEAB07", "#D5D5D3", "#24281A"),
  Moonrise2 = c("#798E87", "#C27D38", "#CCC591", "#29211F"),
  Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"),
  Cavalcanti1 = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15"),
  GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236"),
  GrandBudapest2 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4"),
  IsleofDogs1 = c("#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680"),
  IsleofDogs2 = c("#EAD3BF", "#AA9486", "#B6854D", "#39312F", "#1C1718"),
  FrenchDispatch = c("#90D4CC", "#BD3027", "#B0AFA2", "#7FC0C6", "#9D9C85")
)

set1 <- c("#C43A31", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#F0E442", "#000000", "#56B4E9")
dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

## Shared color-blind-friendly mutation palette for all paper figures.
## TP53 is intentionally a deeper crimson than standard Okabe-Ito vermillion.
pal_variant <- c(
  "DNMT3A" = "#0072B2",
  "TET2"   = "#009E73",
  "TP53"   = "#C43A31",
  "PPM1D"  = "#CC79A7",
  "BCOR"   = "#E69F00",
  "ATM"    = "#56B4E9",
  "EZH2"   = "#000000",
  "IDH1"   = "#F0E442",
  "ASXL1"  = "#999999",
  "MYD88"  = "#6A3D9A",
  "SF3B1"  = "#8C564B",
  "U2AF1"  = "#CC6677"
)
pal_variant_with_wt <- c(pal_variant, "WT" = "lightgray")

## Lineage palette refinements requested for color-blind contrast.
## NK is a saturated purple to separate it from gray Other/Unknown cells.
## Myeloid progenitor-like cells are dark teal to separate from orange myeloid cells.
NK_CELL_COLOR <- "#AA4499"
MYELOID_PROGENITOR_LIKE_COLOR <- "#00A6D6"

set1 <- unname(pal_variant[c("TP53", "DNMT3A", "PPM1D", "TET2", "BCOR", "ATM", "EZH2", "IDH1")])
#python

# define tiny axis things because Ash said so.  These will be imported into the umap graphs through theme to make silly tiny axis labels.
library(ggh4x)
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
    )

###
#unix
#cd [local path]
#awk 'FNR==1 && NR!=1{next;}{print}' *.csv > merged.tsne.csv
##

library(readxl)
setwd(study_path("study_path_001"))

pct_vs_stamp <- read_delim(study_path("input_percent_vs_stamp_scdna"))
merged_scDNA <- read_delim(study_path("input_merged_tsne"))


## ---- Robust Variant_Call / clone-label recovery ----
## The current merged.tsne.csv export does not contain a literal Variant_Call
## column.  In this file structure:
##   - Sample_ID contains the [study-specific]/[study-specific]/... sample code that the original
##     script used as Variant_Call for patient/sample remapping.
##   - Original_clone_ID contains the per-cell clone/mutation label that the
##     original script used as label_final.
## We recreate those expected column names before running the rest of the
## figure code, and then re-check them before the CreateTableOne sections.
first_existing_col <- function(candidates, nm) {
  hit <- intersect(candidates, nm)
  if (length(hit) > 0) return(hit[[1]])

  norm <- function(x) stringr::str_to_lower(stringr::str_replace_all(as.character(x), "[^A-Za-z0-9]+", ""))
  idx <- match(norm(candidates), norm(nm), nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) > 0) return(nm[[idx[[1]]]])

  NA_character_
}

ensure_scDNA_expected_columns <- function(dat) {
  ## 1) Clone/mutation label expected downstream as label_final.
  if (!"label_final" %in% names(dat)) {
    label_col <- first_existing_col(
      c("label_final", "Label_Final", "label.final", "Original_clone_ID", "Original Clone ID",
        "Original.clone.ID", "Clone_ID", "clone_id", "clone", "Original_Clone"),
      names(dat)
    )
    if (!is.na(label_col)) {
      dat$label_final <- as.character(dat[[label_col]])
      message("[INFO] label_final was not present; using column '", label_col, "' as label_final.")
    } else {
      dat$label_final <- "WT"
      warning("label_final was not present and no clone-label column was found; using 'WT'.")
    }
  }

  ## 2) Sample-code column expected downstream as Variant_Call.
  ## Prefer a true Variant_Call if present. Otherwise, use Sample_ID when it
  ## carries the [study-specific]/[study-specific]/[study-specific]/CAR1/CAR2 sample-code structure in this file.
  if (!"Variant_Call" %in% names(dat)) {
    vc_col <- first_existing_col(
      c("Variant_Call", "Variant Call", "Variant.Call", "variant_call", "variant.call",
        "VariantCall", "variantcall", "Variant_Code", "Variant Code", "Variant.Code"),
      names(dat)
    )

    if (!is.na(vc_col)) {
      dat$Variant_Call <- as.character(dat[[vc_col]])
      message("[INFO] Variant_Call was not present; using column '", vc_col, "' as Variant_Call.")
    } else if ("Sample_ID" %in% names(dat)) {
      sample_vals <- as.character(dat$Sample_ID)
      sample_code_fraction <- mean(stringr::str_detect(sample_vals, "^([A-Z][0-9][A-Z][0-9]|CAR[0-9]+|G[0-9]H[0-9])$"), na.rm = TRUE)
      if (is.nan(sample_code_fraction)) sample_code_fraction <- 0

      if (sample_code_fraction > 0.25) {
        dat$Variant_Call <- sample_vals
        message(study_setting("study_value_001"))
      } else {
        ## Last-resort fallback for summary tables only.
        dat$Variant_Call <- as.character(dat$label_final)
        message("[WARN] Variant_Call was not present and Sample_ID did not look like sample codes; using label_final as a fallback.")
      }
    } else {
      dat$Variant_Call <- as.character(dat$label_final)
      message("[WARN] Variant_Call was not present and Sample_ID was unavailable; using label_final as a fallback.")
    }
  }

  dat
}

## Backward-compatible wrapper used later in the script.
ensure_variant_call_for_scDNA <- ensure_scDNA_expected_columns

merged_scDNA <- ensure_scDNA_expected_columns(merged_scDNA)


# first show that the scDNA vaf and stamp vaf are the same


ggscatter(pct_vs_stamp %>% drop_na(STAMP), x = "scDNA", y = "STAMP",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          add.params = list(color = "black", fill = "lightgray"),
          xlab = "VAF in scDNA BM", ylab = "VAF in Clinical BM") +
          theme_classic(base_size = 28, base_family = "Arial") + geom_hline(yintercept = -2, linetype = "dashed") + geom_vline(xintercept = -2, linetype = "dashed") + ggtitle("Variant allele frequency in scDNA vs clinical")

# n = 25 shared variants

# Substantial relabeling necessary

# Relabel a single cluster that was not CAR
merged_scDNA$Sample_ID <- ifelse(merged_scDNA$Sample_ID == study_setting("study_value_002"), "Non-CAR", merged_scDNA$Sample_ID)
merged_scDNA$Sample_ID <- ifelse(merged_scDNA$Sample_ID == study_setting("study_value_003"), "CAR", merged_scDNA$Sample_ID)

# At patient IDs
merged_scDNA$patient_id <- study_setting("study_value_004")
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_005"), "Control", merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_006"), "Product_1", merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_007"), "Product_2", merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_008"), study_setting("study_value_009"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_010"), study_setting("study_value_011"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_012"), study_setting("study_value_013"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_014"), study_setting("study_value_015"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_016"), study_setting("study_value_017"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_018"), study_setting("study_value_019"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_020"), study_setting("study_value_021"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_022"), study_setting("study_value_004"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_023"), study_setting("study_value_024"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_025"), study_setting("study_value_026"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_027"), study_setting("study_value_004"), merged_scDNA$patient_id)
merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_028"), study_setting("study_value_029"), merged_scDNA$patient_id)

# Define cases and controls
merged_scDNA$Sample_Type <- "Case"
merged_scDNA$Sample_Type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_005"), "Healthy Control", merged_scDNA$Sample_Type)
merged_scDNA$Sample_Type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_006"), "Product", merged_scDNA$Sample_Type)
merged_scDNA$Sample_Type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_007"), "Product", merged_scDNA$Sample_Type)


# Define disease states
merged_scDNA$disease_type <- "CCUS"
merged_scDNA$disease_type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_005"), "Control", merged_scDNA$disease_type)
merged_scDNA$disease_type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_006"), "Product", merged_scDNA$disease_type)
merged_scDNA$disease_type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_007"), "Product", merged_scDNA$disease_type)
merged_scDNA$disease_type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_022"), "MDS", merged_scDNA$disease_type)
merged_scDNA$disease_type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_012"), "MDS", merged_scDNA$disease_type)
merged_scDNA$disease_type <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_027"), "Early_MDS", merged_scDNA$disease_type)

# Define mutation types.  Grouped mutations are always defined as the DDR variant.
merged_scDNA$Variant_Type <- "WT"
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "TP53_het", "TP53", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "TP53_EZH2_hom", "TP53", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "TP53_hom", "TP53", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_Het", "DNMT3A", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "PPM1D_Het", "PPM1D", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_2", "DNMT3A", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_1_ATM_het", "ATM", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "BCOR_Het", "BCOR", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "TP53_double", "TP53", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_2_TP53", "TP53", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "PPM1D:p.Q524*_Het", "PPM1D", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "PPM1D:p.L450*_Het", "PPM1D", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "PPM1D:p.S570*_Het", "PPM1D", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "EZH2_Het", "EZH2", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_R551Nfs", "PPM1D", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "TET2_Het", "TET2", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R635W_Het", "DNMT3A", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R882H_Het", "DNMT3A", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_1", "DNMT3A", merged_scDNA$Variant_Type)
merged_scDNA$Variant_Type <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_E540Dfs", "DNMT3A", merged_scDNA$Variant_Type)

# Define mutation types.  Grouped mutations are always defined as the DDR variant.  Differentiate mono/bialleic tp53
merged_scDNA$Variant_Type2 <- "WT"
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "TP53_het", "Monoallelic TP53", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "TP53_EZH2_hom", "Biallelic TP53", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "TP53_hom", "Biallelic TP53", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_Het", "DNMT3A", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "PPM1D_Het", "PPM1D", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_2", "DNMT3A", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_1_ATM_het", "ATM", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "BCOR_Het", "BCOR", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "TP53_double", "Biallelic TP53", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_2_TP53", "Monoallelic TP53", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "PPM1D:p.Q524*_Het", "PPM1D", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "PPM1D:p.L450*_Het", "PPM1D", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "PPM1D:p.S570*_Het", "PPM1D", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "EZH2_Het", "EZH2", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_R551Nfs", "PPM1D", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "TET2_Het", "TET2", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R635W_Het", "DNMT3A", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R882H_Het", "DNMT3A", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_1", "DNMT3A", merged_scDNA$Variant_Type2)
merged_scDNA$Variant_Type2 <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_E540Dfs", "DNMT3A", merged_scDNA$Variant_Type2)


# DDR or not
# Define mutation types.  Grouped mutations are always defined as the DDR variant.
merged_scDNA$Variant_Origin <- "WT"
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "TP53_het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "TP53_EZH2_hom", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "TP53_hom", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_Het", "Age Related", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "PPM1D_Het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_2", "Age Related", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_1_ATM_het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "BCOR_Het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "TP53_double", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_2_TP53", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "PPM1D:p.Q524*_Het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "PPM1D:p.L450*_Het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "PPM1D:p.S570*_Het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "EZH2_Het", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_R551Nfs", "DDR", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "TET2_Het", "Age Related", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R635W_Het", "Age Related", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R882H_Het", "Age Related", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_1", "Age Related", merged_scDNA$Variant_Origin)
merged_scDNA$Variant_Origin <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_E540Dfs", "DDR", merged_scDNA$Variant_Origin)


# Mark intermediary variant
merged_scDNA$Intermediate <- paste(merged_scDNA$Variant_Call, ":", merged_scDNA$label_final, sep = "")

# Define sepcific mutations/clones
merged_scDNA$Variant_Specific <- "WT"
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_030"), "TP53_G266V", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_031"), "TP53_F270S", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_032"), "TP53_F270S", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_033"), "DNMT3A_R729W", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_034"), "DNMT3A_splicing", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_035"), "DNMT3A_L547H", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_036"), "DNMT3A_M801V", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_037"), "DNMT3A_C497R", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_038"), "TET2_H1676fs", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_039"), "TET2_W1003*", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_040"), "DNMT3A_splicing", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_041"), "DNMT3A_T862I", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_041"), "DNMT3A_T862I", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "TP53_EZH2_hom", "TP53_F270S_TP53_del_EZH2_E249K", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "TP53_hom", "TP53_F270S_TP53_del", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "PPM1D_Het", "PPM1D_S412fs", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A_1_ATM_het", "DNMT3A_F752del_ATM_G2020V", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "BCOR_Het", "BCOR_V896fs", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "TP53_double", "TP53_Y205*_TP53_R282W", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A_2_TP53", "DNMT3A_T862I_TP53_splicing", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "PPM1D:p.Q524*_Het", "PPM1D_Q524*", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "PPM1D:p.L450*_Het", "PPM1D_L450*", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "PPM1D:p.S570*_Het", "PPM1D_S570*", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "EZH2_Het", "EZH2_Y741Lfs*22", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_R551Nfs", "DNMT3A_R882H_PPM1D_R551Nfs", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R635W_Het", "DNMT3A_R635W", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A:p.R882H_Het", "DNMT3A_R882H", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A_1", "DNMT3A_E733*", merged_scDNA$Variant_Specific)
merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$label_final == "DNMT3A_R882H_PPM1D_E540Dfs", "DNMT3A_R882H_PPM1D_E540Dfs", merged_scDNA$Variant_Specific)

# Broadly define cell types
merged_scDNA$Cell_Family <- "Myeloid"
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "myeloid NOS", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "pDC", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "neutrophil_cd10", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "erythroid progenitors", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "late CMP", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "monocytes", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "neutrophils", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "basophils", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "early CMP", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd1cDC", "Myeloid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "HSC", "Stem Cell", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "unclassifiable", "Unknown", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "NK" , "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25-" , "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd4-cd8- Tcells", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25+", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd8+cd25+", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd8+cd25-", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "early lymphocytes", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25+", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "mature Bcells", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "pre/pro Bcells", "Lymphoid", merged_scDNA$Cell_Family)
merged_scDNA$Cell_Family <- ifelse(merged_scDNA$OmiqFilter == "CD4+cd8+", "Lymphoid", merged_scDNA$Cell_Family)

# Specifically Designate cell types
merged_scDNA$Cell_Type <- "Myeloid"
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "myeloid NOS", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "pDC", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "neutrophil_cd10", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "erythroid progenitors", "Erythroid Progenitor", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "late CMP", "Myeloid Progenitor", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "monocytes", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "neutrophils", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "basophils", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "early CMP", "Myeloid Progenitor", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd1cDC", "Myeloid Effector", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "HSC", "Hematopoietic Stem Cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "unclassifiable", "Unknown", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "NK" , "NK Cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25-" , "CD4+ T-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd4-cd8- Tcells", "Double Negative T-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25+", "CD4+ T-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd8+cd25+", "CD8+ T-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd8+cd25-", "CD8+ T-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "early lymphocytes", "Lymphoid Progenitor", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25+", "CD4+ T-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "mature Bcells", "B-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "pre/pro Bcells", "B-cell", merged_scDNA$Cell_Type)
merged_scDNA$Cell_Type <- ifelse(merged_scDNA$OmiqFilter == "CD4+cd8+", "Double Positive T-cell", merged_scDNA$Cell_Type)

# Define T-cell populations

# Broadly define cell types
merged_scDNA$Tcells <- "Other"
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "myeloid NOS", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "pDC", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "neutrophil_cd10", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "erythroid progenitors", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "late CMP", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "monocytes", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "neutrophils", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "basophils", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "early CMP", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd1cDC", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "HSC", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "unclassifiable", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "NK" , "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25-" , "T-cell", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd4-cd8- Tcells", "T-cell", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25+", "T-cell", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd8+cd25+", "T-cell", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd8+cd25-", "T-cell", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "early lymphocytes", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "cd4+cd25+", "T-cell", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "mature Bcells", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "pre/pro Bcells", "Other", merged_scDNA$Tcells)
merged_scDNA$Tcells <- ifelse(merged_scDNA$OmiqFilter == "CD4+cd8+", "T-cell", merged_scDNA$Tcells)


# rename columns due to wonky naming
colnames(merged_scDNA)[5] <- "barcode_again"
colnames(merged_scDNA)[6] <- "CAR"
colnames(merged_scDNA)[7] <- "Sample_ID"
colnames(merged_scDNA)[8] <- "barcode_again_again"
colnames(merged_scDNA)[9] <- "Original_clone_ID"
colnames(merged_scDNA)[10] <- "Sample_ID_again"
colnames(merged_scDNA)[10] <- "barcode_again_again_again"

# make column of unique barcodes
merged_scDNA$uniq_barcode <- paste(merged_scDNA$Sample_ID, ":", merged_scDNA$barcode, sep = "")

# make a column of CAR positive mutant cells
# Make prettier version of car names
merged_scDNA$CAR_detected <- "Absent"
merged_scDNA$CAR_detected <- ifelse(merged_scDNA$CAR == "CAR", "Present", merged_scDNA$CAR_detected)

####
merged_scDNA %>% subset(CAR_detected != "Absent" & Variant_Specific != "WT") %>% subset(Cell_Family == "Lymphoid") -> car_pos_mutants
car_pos_mutants$uniq_barcode -> car_pos_mutants
merged_scDNA$CAR_mutant <- ifelse(merged_scDNA$uniq_barcode %in% car_pos_mutants, "Yes", "No")


# create a version that is without product
merged_scDNA_noprod <- subset(merged_scDNA, Sample_Type != "Product")

# Make graph of cell types
# Highlight cell type
cols <- c("Myeloid Effector" = "#F98400", "NK Cell" = NK_CELL_COLOR, "CD4+ T-cell" = "#C6CDF7", "CD8+ T-cell" = "#7294D4", "Double Negative T-cell" = "#046C9A", "Double Positive T-cell" = "#ABDDDE", "Myeloid Progenitor" = MYELOID_PROGENITOR_LIKE_COLOR, "B-cell" = "#E6A0C4", "Unknown" = "lightgray", "Erythroid Progenitor" = "#E1BD6D", "Lymphoid Progenitor" = "#35274A", "Hematopoietic Stem Cell" = "black")

q <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = Cell_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Cell Type") +
    xlab("t-SNE1") + ylab("t-SNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Cell Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

# Highlight mutations
cols <- pal_variant_with_wt

r <- ggplot(merged_scDNA_noprod %>% arrange(desc(Variant_Type)), aes(x = optsne_1, y = optsne_2, color = Variant_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Variant") +
    xlab("t-SNE1") + ylab("t-SNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Variant Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r + plot_layout(ncol = 2)
ggsave(out_file("Cell_type_variants.png"), width = 15, height = 4.06, units = "in")
ggsave(out_file("Cell_type_variants.pdf"), width = 15, height = 4.06, units = "in", device = cairo_pdf)
ggsave(out_file("Cell_type_variants.eps"), width = 15, height = 4.06, units = "in")

############################################################
## Diagnostics for Cell_type_variants plots (merged_scDNA_noprod)
## - total cells
## - CH mutant cells (non-WT Variant_Type)
############################################################

stopifnot(all(c("Variant_Type","optsne_1","optsne_2") %in% colnames(merged_scDNA_noprod)))

n_total_plot <- nrow(merged_scDNA_noprod)

# Define CH-mutant as anything not WT (and not NA)
n_CH_mut_plot <- merged_scDNA_noprod %>%
  dplyr::filter(!is.na(Variant_Type), Variant_Type != "WT") %>%
  nrow()

cat("\n=== [DIAGNOSTIC: merged_scDNA_noprod] ===\n")
cat("Total cells in plot dataset: ", n_total_plot, "\n", sep = "")
cat("CH mutant cells (Variant_Type != WT): ", n_CH_mut_plot, "\n", sep = "")

cat("\n[CHECK] Variant_Type breakdown:\n")
print(sort(table(merged_scDNA_noprod$Variant_Type, useNA = "ifany"), decreasing = TRUE))

# Optional: if you want to annotate the plot titles/subtitles
sub_diag <- paste0("n=", n_total_plot, " | CH mutant=", n_CH_mut_plot)


cols <- c("Absent" = "lightgray", "Present" = "darkgreen")
s <- ggplot(merged_scDNA_noprod %>% arrange(desc(CAR)), aes(x = optsne_1, y = optsne_2, color = CAR_detected)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "CAR DNA", breaks = c("Present", "Absent")) +
    xlab("t-SNE1") + ylab("t-SNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Axi-cel Detection") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))
ggsave(out_file("axice_positive.png"), width = 15, height = 4.06, units = "in")
ggsave(out_file("axice_positive.pdf"), width = 15, height = 4.06, units = "in", device = cairo_pdf)


# highlight axi-cel -> this is the group used in figure 6

merged_scDNA_case <- subset(merged_scDNA , Sample_Type == "Case")
merged_scDNA_control <- subset(merged_scDNA , Sample_Type == "Healthy Control")
merged_scDNA_product <- subset(merged_scDNA , Sample_Type == "Product")


cols <- c("Absent" = "lightgray", "Present" = "darkgreen")
s <- ggplot(merged_scDNA_case %>% arrange(desc(CAR)), aes(x = optsne_1, y = optsne_2, color = CAR_detected)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "CAR DNA") +
    xlab("t-SNE1") + ylab("t-SNE2") + xlim(-60,60) + ylim(-60,60) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Post-CAR Patients") +  theme(legend.position="none") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

cols <- c("Absent" = "lightgray", "Present" = "darkgreen")
t <- ggplot(merged_scDNA_control %>% arrange(desc(CAR)), aes(x = optsne_1, y = optsne_2, color = CAR_detected)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "CAR DNA") +
    xlab("t-SNE1") + ylab("t-SNE2") + xlim(-60,60) + ylim(-60,60) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Healthy Control") + theme(legend.position="none") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

cols <- c("Absent" = "lightgray", "Present" = "darkgreen")
u <- ggplot(merged_scDNA_product %>% arrange(desc(CAR)), aes(x = optsne_1, y = optsne_2, color = CAR_detected)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "CAR DNA") +
    xlab("t-SNE1") + ylab("t-SNE2") + xlim(-60,60) + ylim(-60,60) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Product") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


s+t+u + plot_layout(ncol = 3)

ggsave(out_file("axicel_detection.png"), width = 16, height = 4.06, units = "in")
ggsave(out_file("axicel_detection.pdf"), width = 16, height = 4.06, units = "in", device = cairo_pdf)


############################################################
## DIAGNOSTICS: n cells + n CAR+ for the final three graphs
############################################################

# Robust CAR+ flag (handles if CAR_detected is missing/dirty)
flag_car_present <- function(df) {
  if ("CAR_detected" %in% names(df)) {
    return(!is.na(df$CAR_detected) & df$CAR_detected == "Present")
  }
  if ("CAR" %in% names(df)) {
    # in your file CAR seems to be "CAR" vs something else
    return(!is.na(df$CAR) & df$CAR == "CAR")
  }
  stop("No CAR_detected or CAR column found in the provided data frame.")
}

count_cells_and_car <- function(df, label) {
  car_flag <- flag_car_present(df)
  tibble(
    Panel      = label,
    n_cells    = nrow(df),
    n_CAR_pos  = sum(car_flag, na.rm = TRUE),
    frac_CAR   = if_else(nrow(df) > 0, sum(car_flag, na.rm = TRUE) / nrow(df), NA_real_)
  )
}

diag_three <- bind_rows(
  count_cells_and_car(merged_scDNA_case,    "Post-CAR Patients"),
  count_cells_and_car(merged_scDNA_control, "Healthy Control"),
  count_cells_and_car(merged_scDNA_product, "Product")
)

cat("\n=== [DIAGNOSTIC] Final three t-SNE panels: cell counts + CAR+ counts ===\n")
print(diag_three)

cat("\n=== [DIAGNOSTIC] CAR labeling sanity checks ===\n")
if ("CAR_detected" %in% names(merged_scDNA_case)) {
  cat("\nPost-CAR Patients: table(CAR_detected)\n")
  print(table(merged_scDNA_case$CAR_detected, useNA = "ifany"))
}
if ("CAR_detected" %in% names(merged_scDNA_control)) {
  cat("\nHealthy Control: table(CAR_detected)\n")
  print(table(merged_scDNA_control$CAR_detected, useNA = "ifany"))
}
if ("CAR_detected" %in% names(merged_scDNA_product)) {
  cat("\nProduct: table(CAR_detected)\n")
  print(table(merged_scDNA_product$CAR_detected, useNA = "ifany"))
}

# Optional: also show totals for the non-product combined plot you saved earlier
if (exists("merged_scDNA_noprod")) {
  diag_noprod <- count_cells_and_car(merged_scDNA_noprod, "All (No Product)")
  cat("\n=== [DIAGNOSTIC] All (No Product) panel ===\n")
  print(diag_noprod)
}

# highlight axicel inside a mutation
cols <- c("No" = "lightgray", "Yes" = "Dodgerblue")
y <- ggplot(merged_scDNA_noprod  %>% arrange(CAR_mutant), aes(x = optsne_1, y = optsne_2, color = CAR_mutant)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "CAR with CH Mutation", breaks = c("Yes", "No")) +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Axi-cel within CH Mutant Cell") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


################# Highlight CAR CH mutations ##########
cols <- c("No" = "lightgray", "Yes" = "Dodgerblue")
y <- ggplot(merged_scDNA_noprod  %>% arrange(CAR_mutant), aes(x = optsne_1, y = optsne_2, color = CAR_mutant)) + geom_point_ai_rast(size = 0.5) +
    scale_colour_manual(values = cols, name = "CAR with CH Mutation", breaks = c("Yes", "No")) +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Axi-cel within CH Mutant Cell") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

cols <- pal_variant_with_wt


z <- ggplot(merged_scDNA_noprod  %>% subset(CAR_detected == "Present"), aes(x = optsne_1, y = optsne_2, color = Variant_Type)) + geom_point_ai_rast(size = 0.5) +
    scale_colour_manual(values = cols, name = "Variant") +
    xlab("TSNE1") + ylab("TSNE2") + xlim(-60,60) + ylim(-60,60) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("CAR+ Cells") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

y+z

ggsave(out_file("axicel_CH_scDNA.png"), width = 16, height = 5.06, units = "in")
ggsave(out_file("axicel_CH_scDNA.pdf"), width = 16, height = 5.06, units = "in", device = cairo_pdf)
############################################################
## DIAGNOSTIC: number of CAR+CH+ cells (CAR_mutant == "Yes")
############################################################

# ensure CAR_mutant exists
stopifnot("CAR_mutant" %in% colnames(merged_scDNA_noprod))

# count CAR+CH+ (as defined by your CAR_mutant logic)
n_car_ch_pos <- sum(merged_scDNA_noprod$CAR_mutant == "Yes", na.rm = TRUE)
n_total      <- nrow(merged_scDNA_noprod)

cat("\n=== [DIAGNOSTIC] CAR+CH+ cells (CAR_mutant == 'Yes') ===\n")
cat("Total cells (noprod): ", n_total, "\n", sep = "")
cat("CAR+CH+ cells:        ", n_car_ch_pos, "\n", sep = "")
cat("Fraction CAR+CH+:     ", round(n_car_ch_pos / max(n_total, 1), 4), "\n", sep = "")

cat("\nTable of CAR_mutant:\n")
print(table(merged_scDNA_noprod$CAR_mutant, useNA = "ifany"))

# Optional: among CAR+ cells only, how many are CAR+CH+?
if ("CAR_detected" %in% colnames(merged_scDNA_noprod)) {
  n_car_total <- sum(merged_scDNA_noprod$CAR_detected == "Present", na.rm = TRUE)
  n_car_ch_in_car <- sum(merged_scDNA_noprod$CAR_detected == "Present" &
                           merged_scDNA_noprod$CAR_mutant == "Yes", na.rm = TRUE)

  cat("\nAmong CAR+ cells:\n")
  cat("CAR+ total:          ", n_car_total, "\n", sep = "")
  cat("CAR+CH+ (subset):    ", n_car_ch_in_car, "\n", sep = "")
  cat("Fraction of CAR+ cells that are CH+: ",
      round(n_car_ch_in_car / max(n_car_total, 1), 4), "\n", sep = "")
}

# write CAR+ table for manipulation
CAR_pos <- as.data.frame(table(merged_scDNA_noprod$CAR, merged_scDNA_noprod$Cell_Type, merged_scDNA_noprod$Sample_ID))

# Isolate CAR+ samples
CAR_pos_pos <- subset(CAR_pos, Var1 == "CAR")

# Get fraction of CAR in each group
CAR_pos_pos %>%
dplyr::group_by(Var3) %>%
  dplyr::mutate(Freq / sum(Freq)) -> CAR_pos_fraction
  
# Sum all CAR+ cells
agg_df <- aggregate(CAR_pos_pos$Freq, by=list(CAR_pos_pos$Var2), FUN=sum)

cols <- c("Myeloid Effector" = "#F98400", "NK Cell" = NK_CELL_COLOR, "CD4+ T-cell" = "#C6CDF7", "CD8+ T-cell" = "#7294D4", "Double Negative T-cell" = "#046C9A", "Double Positive T-cell" = "#ABDDDE", "Myeloid Progenitor" = MYELOID_PROGENITOR_LIKE_COLOR, "B-cell" = "#E6A0C4", "Unknown" = "lightgray", "Erythroid Progenitor" = "#E1BD6D", "Lymphoid Progenitor" = "#35274A", "Hematopoietic Stem Cell" = "black")

p <- ggplot(CAR_pos_fraction, aes(y = Var2, x = `Freq/sum(Freq)`, fill = Var2)) +
    geom_boxplot() +
    scale_fill_manual(values = cols, name = "CAR+ Cell Fraction") +
    theme_classic(base_size=20, base_family = "Arial") +
    xlab("Fraction of CAR+ Cells in Group") + ylab("Cell Type") +
    theme(legend.position = "none")

library(dplyr)
library(tidyr)

# helper: make safe column keys from Cell_Type labels
key_from_label <- function(lbl) {
  x <- gsub("[^A-Za-z0-9]+", "_", lbl)
  gsub("_$", "", x)
}

# Build per-sample CAR+ totals + fractions, wide format
car_wide <- CAR_pos_pos %>%
  mutate(
    Sample_ID = Var3,
    Cell_Type = Var2,
    CAR_n     = Freq,
    Cell_Key  = key_from_label(Cell_Type)
  ) %>%
  group_by(Sample_ID) %>%
  mutate(
    CAR_total = sum(CAR_n),
    CAR_frac  = ifelse(CAR_total > 0, CAR_n / CAR_total, NA_real_)
  ) %>%
  ungroup() %>%
  select(Sample_ID, Cell_Key, CAR_n, CAR_frac, CAR_total) %>%
  pivot_wider(
    names_from   = Cell_Key,
    values_from  = c(CAR_n, CAR_frac),
    values_fill  = list(CAR_n = 0, CAR_frac = 0),
    names_sep    = "__"
  )

# Keys/columns for CD4 vs CD8
cd4_key <- key_from_label("CD4+ T-cell")
cd8_key <- key_from_label("CD8+ T-cell")

frac_cd4_col <- paste0("CAR_frac__", cd4_key)
frac_cd8_col <- paste0("CAR_frac__", cd8_key)

n_cd4_col <- paste0("CAR_n__", cd4_key)
n_cd8_col <- paste0("CAR_n__", cd8_key)

# Safe paired Wilcoxon runner (avoids crashing when all diffs are 0)
run_paired_wilcox <- function(df, x, y, label) {
  df2 <- df %>%
    filter(CAR_total > 0) %>%                      # only samples with any CAR+ cells
    mutate(x = .data[[x]], y = .data[[y]]) %>%
    filter(!is.na(x) & !is.na(y))

  n_pairs_all <- nrow(df2)
  n_pairs_informative <- sum((df2$x + df2$y) > 0)  # at least one of the two is >0
  df_inf <- df2 %>% filter((x + y) > 0)

  cat("\n============================\n")
  cat(label, "\n")
  cat("N pairs (CAR_total>0): ", n_pairs_all, "\n", sep = "")
  cat("N pairs (informative; x+y>0): ", n_pairs_informative, "\n", sep = "")

  if (nrow(df_inf) < 2) {
    cat("Not enough informative pairs to test.\n")
    return(invisible(NULL))
  }
  if (all(df_inf$x == df_inf$y)) {
    cat("All paired values are identical (all differences = 0); Wilcoxon not informative.\n")
    return(invisible(NULL))
  }

  wt <- wilcox.test(df_inf$x, df_inf$y, paired = TRUE, exact = FALSE)

  cat("Median(x): ", signif(median(df_inf$x), 4),
      " | Median(y): ", signif(median(df_inf$y), 4), "\n", sep = "")
  cat("Median(x - y): ", signif(median(df_inf$x - df_inf$y), 4), "\n", sep = "")
  cat("Wilcoxon V: ", unname(wt$statistic),
      " | p: ", format.pval(wt$p.value, digits = 3, eps = 1e-300), "\n", sep = "")

  invisible(list(df = df_inf, test = wt))
}

# --- PRIMARY: fractions (matches your plot concept) ---
res_frac <- run_paired_wilcox(
  car_wide,
  x = frac_cd4_col,
  y = frac_cd8_col,
  label = "[Paired Wilcoxon] CAR+ FRACTION: CD4 vs CD8"
)

# --- OPTIONAL: counts (raw CAR+ CD4 vs CAR+ CD8 counts per sample) ---
res_n <- run_paired_wilcox(
  car_wide,
  x = n_cd4_col,
  y = n_cd8_col,
  label = "[Paired Wilcoxon] CAR+ COUNTS: CD4 vs CD8"
)


CAR_pos_by_sample <- as.data.frame(table(merged_scDNA_noprod$CAR, merged_scDNA_noprod$Sample_ID))
CAR_pos_by_sample_full <- as.data.frame(table(merged_scDNA$CAR, merged_scDNA$Sample_ID))

# Isolate DDR mutations in lymphoid

merged_scDNA_noprod_lymphoid <- subset(merged_scDNA_noprod, Cell_Family == "Lymphoid")
merged_scDNA_noprod_lymphoid_DDR <- subset(merged_scDNA_noprod_lymphoid, Variant_Origin == "DDR")
merged_scDNA_noprod_lymphoid_Age <- subset(merged_scDNA_noprod_lymphoid, Variant_Origin == "Age Related")

# Table
DDR_lymphoid <- as.data.frame(table(merged_scDNA_noprod_lymphoid_DDR$Cell_Type))
DDR_lymphoid %>%
  dplyr::mutate(Freq / sum(Freq)) -> DDR_lymphoid

Age_lymphoid <- as.data.frame(table(merged_scDNA_noprod_lymphoid_Age$Cell_Type))
Age_lymphoid %>%
  dplyr::mutate(Freq / sum(Freq)) -> Age_lymphoid
  
  
All_lymphoid <- as.data.frame(table(merged_scDNA_noprod_lymphoid$Cell_Type, merged_scDNA_noprod_lymphoid$Variant_Origin))
All_lymphoid_mut <- subset(All_lymphoid, Var2 != "WT")

All_lymphoid %>%
  dplyr::mutate(Freq / sum(Freq)) -> All_lymphoid
  
  
# determine the enrichment of axi-cel in all CH variants

merged_scDNA_noprod_tcell <- subset(merged_scDNA_noprod, Tcells == "T-cell")
T_cell_variants <- as.data.frame(table(merged_scDNA_noprod_tcell$Intermediate, merged_scDNA_noprod_tcell$CAR))
T_cell_variants_2 <- as.data.frame(table(merged_scDNA_noprod_tcell$Intermediate))
T_cell_variants_3 <- as.data.frame(table(merged_scDNA_noprod_tcell$CAR))


# Get ATM mutation as anecdote

merged_scDNA_noprod_ATM <- subset(merged_scDNA_noprod, Variant_Type == "ATM")

# Load up pie plots
library(scales)

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

cols <- c("Myeloid Effector" = "#F98400", "NK Cell" = NK_CELL_COLOR, "CD4+ T-cell" = "#C6CDF7", "CD8+ T-cell" = "#7294D4", "Double Negative T-cell" = "#046C9A", "Double Positive T-cell" = "#ABDDDE", "Myeloid Progenitor" = MYELOID_PROGENITOR_LIKE_COLOR, "B-cell" = "#E6A0C4", "Unknown" = "lightgray", "Erythroid Progenitor" = "#E1BD6D", "Lymphoid Progenitor" = "#35274A", "Hematopoietic Stem Cell" = "black")


p <- ggplot(DDR_lymphoid, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = cols, name = "Lymphoid Cell Type") + blank_theme +
    theme(axis.text.x=element_blank()) + ggtitle("DDR Mutations in Lymphoid Cells")
    
q <- ggplot(Age_lymphoid, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = cols, name = "Lymphoid Cell Type") + blank_theme +
    theme(axis.text.x=element_blank()) + ggtitle("Age Related Mutations in Lymphoid Cells")

############# Stacked bar NK plots ########

ggplot(All_lymphoid_mut, aes(fill=Var1, y=Freq, x=Var2)) +
    geom_bar(position="fill", stat="identity", color = "black") +
    theme_classic(base_family = "Arial", base_size = 24) +
    scale_fill_manual(values = cols, name = "Cell Type") +
    ylab("Fraction of Cells") + xlab("Mutation Type") + ggtitle("Total cell fraction by mutation type")

dat <- data.frame(
  "NK_no" = c(1984, 880),
  "NK_yes" = c(540, 1261),
  row.names = c("immatureb", "matureb"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Age", "DDR")

dat

test <- fisher.test(dat)
test


# Myeloid vs lymphoid heatmap
keep <- c("Myeloid", "Lymphoid")
merged_scDNA_noprod_muts <- subset(merged_scDNA_noprod, Variant_Type != "WT")
merged_scDNA_noprod_muts_ML <- subset(merged_scDNA_noprod_muts, Cell_Family %in% keep)

type_vs_table <- as.data.frame(table(merged_scDNA_noprod_muts_ML$Intermediate, merged_scDNA_noprod_muts_ML$Cell_Family))

type_vs_table %>%
dplyr::group_by(Var1) %>%
  dplyr::mutate(Freq / sum(Freq)) -> type_vs_table

write_delim(type_vs_table, out_file("type_vs_table.txt"))


# Get CAR-CH patient data

merged_scDNA_selected_sample <- subset(merged_scDNA_noprod, Sample_ID == study_setting("study_value_028"))
merged_scDNA_selected_sample_CAR <- subset(merged_scDNA_noprod, CAR == "CAR")

# Highlight cell type
cols <- c("Myeloid Effector" = "#F98400", "NK Cell" = NK_CELL_COLOR, "CD4+ T-cell" = "#C6CDF7", "CD8+ T-cell" = "#7294D4", "Double Negative T-cell" = "#046C9A", "Double Positive T-cell" = "#ABDDDE", "Myeloid Progenitor" = MYELOID_PROGENITOR_LIKE_COLOR, "B-cell" = "#E6A0C4", "Unknown" = "lightgray", "Erythroid Progenitor" = "#E1BD6D", "Lymphoid Progenitor" = "#35274A", "Hematopoietic Stem Cell" = "black")

q <- ggplot(merged_scDNA_selected_sample, aes(x = optsne_1, y = optsne_2, color = Cell_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Cell Type") +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Cell Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

# Highlight mutations
cols <- pal_variant_with_wt

r <- ggplot(merged_scDNA_selected_sample %>% arrange(desc(Variant_Type)), aes(x = optsne_1, y = optsne_2, color = Variant_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Variant") +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Variant Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


# highlight axi-cel

cols <- c("Absent" = "lightgray", "Present" = "darkgreen")
s <- ggplot(merged_scDNA_selected_sample %>% arrange(desc(CAR)), aes(x = optsne_1, y = optsne_2, color = CAR_detected)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "CAR DNA") +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Axi-cel Detection") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r+s + plot_layout(ncol = 3)

ggsave(out_file(study_setting("study_value_042")), width = 16, height = 4.06, units = "in")
ggsave(out_file(study_setting("study_value_043")), width = 16, height = 4.06, units = "in", device = cairo_pdf)


# Highlight cell type
cols <- c("Myeloid Effector" = "#F98400", "NK Cell" = NK_CELL_COLOR, "CD4+ T-cell" = "#C6CDF7", "CD8+ T-cell" = "#7294D4", "Double Negative T-cell" = "#046C9A", "Double Positive T-cell" = "#ABDDDE", "Myeloid Progenitor" = MYELOID_PROGENITOR_LIKE_COLOR, "B-cell" = "#E6A0C4", "Unknown" = "lightgray", "Erythroid Progenitor" = "#E1BD6D", "Lymphoid Progenitor" = "#35274A", "Hematopoietic Stem Cell" = "black")

q <- ggplot(merged_scDNA_selected_sample_CAR, aes(x = optsne_1, y = optsne_2, color = Cell_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Cell Type") +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Cell Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

# Highlight mutations
cols <- pal_variant_with_wt

r <- ggplot(merged_scDNA_selected_sample_CAR %>% arrange(desc(Variant_Type)), aes(x = optsne_1, y = optsne_2, color = Variant_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Variant") +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Variant Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


# Highlight expansion of bialleic TP53


merged_scDNA_noprod_TP53 <- subset(merged_scDNA_noprod, Variant_Type == "TP53")

# Highlight cell type
cols <- c("Myeloid Effector" = "#F98400", "NK Cell" = NK_CELL_COLOR, "CD4+ T-cell" = "#C6CDF7", "CD8+ T-cell" = "#7294D4", "Double Negative T-cell" = "#046C9A", "Double Positive T-cell" = "#ABDDDE", "Myeloid Progenitor" = MYELOID_PROGENITOR_LIKE_COLOR, "B-cell" = "#E6A0C4", "Unknown" = "lightgray", "Erythroid Progenitor" = "#E1BD6D", "Lymphoid Progenitor" = "#35274A", "Hematopoietic Stem Cell" = "black")

q <- ggplot(merged_scDNA_noprod_TP53, aes(x = optsne_1, y = optsne_2, color = Cell_Type)) + geom_point_ai_rast(size = 0.2) +
    scale_colour_manual(values = cols, name = "Cell Type") +
    xlab("TSNE1") + ylab("TSNE2") + #xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("TP53+ Variant Localization") + facet_grid(~factor(Variant_Type2, levels = c("Monoallelic TP53", "Biallelic TP53"))) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


# Highlight specific proteins
# Highlight the proteins T cells
q <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD7)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD4)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

s <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD8)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

t <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD56)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r+s+t + plot_layout(ncol = 4)
ggsave(out_file("Tcell_markers.png"), width = 21, height = 4.06, units = "in")
ggsave(out_file("Tcell_markers.pdf"), width = 21, height = 4.06, units = "in", device = cairo_pdf)

# Highlight the proteins progenitor
q <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD34)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD117)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

s <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD38)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

t <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD90)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r+s+t + plot_layout(ncol = 4)
ggsave(out_file("myeloid_progenitor_markers.png"), width = 21, height = 4.06, units = "in")
ggsave(out_file("myeloid_progenitor_markers.pdf"), width = 21, height = 4.06, units = "in", device = cairo_pdf)


# Highlight the proteins myeloid
q <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD16)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD62L)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

s <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD14)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

t <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = FcεRIα)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r+s+t + plot_layout(ncol = 4)
ggsave(out_file("myeloid_markers.png"), width = 21, height = 4.06, units = "in")
ggsave(out_file("myeloid_markers.pdf"), width = 21, height = 4.06, units = "in", device = cairo_pdf)


# Highlight the proteins myeloid
q <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD10)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA_noprod, aes(x = optsne_1, y = optsne_2, color = CD19)) + geom_point_ai_rast(size = 0.2) +
    xlab("TSNE1") + ylab("TSNE2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


q+r + plot_layout(ncol = 2)
ggsave(out_file("b_cell_markers.png"), width = 10.5, height = 4.06, units = "in")
ggsave(out_file("b_cell_markers.pdf"), width = 10.5, height = 4.06, units = "in", device = cairo_pdf)


############################
############################
############################
############################
##### subset axicel positive cells
merged_scDNA_axicel <- subset(merged_scDNA, CAR == "CAR")

table(merged_scDNA_axicel$Variant_Call)
#DNMT3A_2       WT
#       9       35
        
# remake plots with axicel only

# Highlight the mutation call

cols <- pal_variant_with_wt

r <- ggplot(merged_scDNA_axicel, aes(x = UMAP1, y = UMAP2, color = Variant_Type)) + geom_point_ai_rast() +
    scale_colour_manual(values = cols, name = "Detected Mutation") +
    xlab("UMAP1") + ylab("UMAP2") + xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Variant Identification") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

 #   labels = c("Myeloid cell", "NK cell", "CD4 T cell", "CD8 T cell")
# Highlight the cell type
cols <- c("Myeloid cell" = "#E6A0C4", "NK cell" = NK_CELL_COLOR, "CD4 T cell" = "#C6CDF7", "CD8 T cell" = "#7294D4")

q <- ggplot(merged_scDNA_axicel, aes(x = UMAP1, y = UMAP2, color = Cell_Type)) + geom_point_ai_rast() +
    scale_colour_manual(values = cols, name = "Cell Type") +
    xlab("UMAP1") + ylab("UMAP2") + xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Cell Type") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

# Highlight axi-cel
merged_scDNA_axicel

colors <- c("Present" = "darkgreen", "Absent" = "lightgray")
s <- ggplot(merged_scDNA_axicel %>% arrange(CAR_detected), aes(x = UMAP1, y = UMAP2, color = Axicel)) + geom_point_ai_rast() +
    scale_color_manual(values = colors) +
    xlab("UMAP1") + ylab("UMAP2") + labs(color = "Axi-cel") + xlim(-8,8) + ylim(-8,8) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("Axi-cel Positive Cells") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))
    
    
######


# Highlight the proteins T cells
q <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD3)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD4)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

s <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD8)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r+s + plot_layout(ncol = 3)
ggsave(out_file(study_setting("study_value_044")), width = 15, height = 4.06, units = "in")
ggsave(out_file(study_setting("study_value_045")), width = 15, height = 4.06, units = "in", device = cairo_pdf)


### NK markers ###

r <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD45RA)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

s <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD56)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r+s + plot_layout(ncol = 2)
ggsave(out_file(study_setting("study_value_046")), width = 7.5, height = 4.06, units = "in")
ggsave(out_file(study_setting("study_value_047")), width = 7.5, height = 4.06, units = "in", device = cairo_pdf)


# Highlight the proteins

o <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD11b)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

p <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD14)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))


q <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD16)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 6), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD62L)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

s <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD33)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

t <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD34)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

o+p+q+r+s+t + plot_layout(ncol = 3)
ggsave(out_file(study_setting("study_value_048")), width = 15, height = 8.06, units = "in")
ggsave(out_file(study_setting("study_value_049")), width = 15, height = 8.06, units = "in", device = cairo_pdf)


#### B cell

q <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD10)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

r <- ggplot(merged_scDNA, aes(x = UMAP1, y = UMAP2, color = CD19)) + geom_point_ai_rast() +
    xlab("UMAP1") + ylab("UMAP2") +
      scale_color_gradientn(colors = c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07"), limits = c(-5, 5), oob = scales::squish) +
    theme_classic(base_size=20, base_family = "Arial") +
    ggtitle("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last")), axis.title = element_text(hjust = 0)) + guides(x = axis, y = axis, color = guide_legend(override.aes = list(size=7)))

q+r + plot_layout(ncol = 2)
ggsave(out_file(study_setting("study_value_050")), width = 7.5, height = 4.06, units = "in")
ggsave(out_file(study_setting("study_value_051")), width = 7.5, height = 4.06, units = "in", device = cairo_pdf)

# Summarize data without stratification
merged_scDNA <- ensure_variant_call_for_scDNA(merged_scDNA)
myVars1 <- intersect(c("CAR", "Cell_Type", "Variant_Call"), names(merged_scDNA))
All_CAR_Table1 <- CreateTableOne(data = merged_scDNA, vars = myVars1)
print(All_CAR_Table1, quote = TRUE, noSpaces = TRUE, nonnormal = TRUE)

# stratify by variant call; skip cleanly if the source file truly lacks a usable column
if ("Variant_Call" %in% names(merged_scDNA) && dplyr::n_distinct(merged_scDNA$Variant_Call, na.rm = TRUE) > 1) {
  myVars1 <- intersect(c("CAR", "Cell_Type", "Variant_Call"), names(merged_scDNA))
  All_CAR_Table1 <- CreateTableOne(data = merged_scDNA, vars = myVars1, strata = "Variant_Call")
  print(All_CAR_Table1, quote = TRUE, noSpaces = TRUE, nonnormal = TRUE)
} else {
  message("[INFO] Skipping CreateTableOne(..., strata = 'Variant_Call') because Variant_Call is unavailable or has only one level.")
}

# stratify by CAR+
myVars1 <- intersect(c("CAR", "Cell_Type", "Variant_Call"), names(merged_scDNA))
All_CAR_Table1 <- CreateTableOne(data = merged_scDNA, vars = myVars1, strata = "CAR")
print(All_CAR_Table1, quote = TRUE, noSpaces = TRUE, nonnormal = TRUE)


###### Determine amount of cells that are CAR+

# ---- packages ----
library(dplyr)
library(tidyr)

# --- set your object name here if different ---
df <- merged_scDNA_case

# ---- harmonize CAR status & define CD subset ----
dat <- df %>%
  mutate(
    CD = case_when(
      grepl("^CD4\\+", Cell_Type, ignore.case = TRUE) ~ "CD4",
      grepl("^CD8\\+", Cell_Type, ignore.case = TRUE) ~ "CD8",
      TRUE ~ NA_character_
    ),
    # Prefer explicit CAR_detected if present; fall back to CAR column
    CAR_status = case_when(
      CAR_detected %in% c("Present","Absent") ~ ifelse(CAR_detected == "Present", "CAR+", "CAR-"),
      CAR %in% c("CAR","Non-CAR") ~ ifelse(CAR == "CAR", "CAR+", "CAR-"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CD), !is.na(CAR_status))  # keep only CD4/CD8 with known CAR status

# ---- tables for CAR+ and CAR- separately ----
car_pos_tbl <- dat %>%
  filter(CAR_status == "CAR+") %>%
  count(CD, name = "n") %>%
  pivot_wider(names_from = CD, values_from = n, values_fill = 0)

car_neg_tbl <- dat %>%
  filter(CAR_status == "CAR-") %>%
  count(CD, name = "n") %>%
  pivot_wider(names_from = CD, values_from = n, values_fill = 0)

car_pos_tbl
car_neg_tbl

# ---- combined 2x2 contingency table & Fisher's exact test ----
contingency <- dat %>%
  count(CAR_status, CD) %>%
  pivot_wider(names_from = CD, values_from = n, values_fill = 0) %>%
  arrange(match(CAR_status, c("CAR+","CAR-"))) %>%  # order rows: CAR+, CAR-
  as.data.frame()

rownames(contingency) <- contingency$CAR_status
cont_2x2 <- as.matrix(contingency[, c("CD4","CD8")])

contingency
cont_2x2

fisher_res <- fisher.test(cont_2x2)

# ---- quick summary ----
cat("\nCounts (CAR+ vs CAR- by CD subset):\n")
print(cont_2x2)

cat("\nFisher's exact test:\n")
print(fisher_res)

# Helpful derived ratios:
ratio_car_pos <- cont_2x2["CAR+","CD4"] / max(1, cont_2x2["CAR+","CD8"])
ratio_car_neg <- cont_2x2["CAR-","CD4"] / max(1, cont_2x2["CAR-","CD8"])

cat(sprintf("\nCD4:CD8 ratio  (CAR+): %.3f\n", ratio_car_pos))
cat(sprintf("CD4:CD8 ratio  (CAR-): %.3f\n", ratio_car_neg))
if (!is.nan(fisher_res$estimate)) {
  cat(sprintf("Odds ratio (CD4 enrichment in CAR+ vs CAR-): %.3f\n", fisher_res$estimate))
}

####### now perform for CAR+CH+ vs CAR+CH- #########
# ---- classify CH status ----
dat_ch <- dat %>%
  filter(CAR_status == "CAR+") %>%   # restrict to CAR+ cells
  mutate(
    CH_status = ifelse(Variant_Type != "WT", "CH+", "CH-")
  ) %>%
  filter(CH_status %in% c("CH+","CH-"))

# ---- contingency table: CH_status x CD subset ----
contingency_ch <- dat_ch %>%
  count(CH_status, CD) %>%
  tidyr::pivot_wider(names_from = CD, values_from = n, values_fill = 0) %>%
  arrange(match(CH_status, c("CH+","CH-"))) %>%
  as.data.frame()

rownames(contingency_ch) <- contingency_ch$CH_status
cont_2x2_ch <- as.matrix(contingency_ch[, c("CD4","CD8")])

contingency_ch
cont_2x2_ch

# ---- fisher exact test ----
fisher_res_ch <- fisher.test(cont_2x2_ch)

# ---- ratios ----
ratio_ch_pos <- cont_2x2_ch["CH+","CD4"] / max(1, cont_2x2_ch["CH+","CD8"])
ratio_ch_neg <- cont_2x2_ch["CH-","CD4"] / max(1, cont_2x2_ch["CH-","CD8"])

cat("\nCounts (CAR+CH+ vs CAR+CH- by CD subset):\n")
print(cont_2x2_ch)

cat("\nFisher's exact test (CD4:CD8 ratio difference in CAR+ cells by CH status):\n")
print(fisher_res_ch)

cat(sprintf("\nCD4:CD8 ratio  (CAR+CH+): %.3f\n", ratio_ch_pos))
cat(sprintf("CD4:CD8 ratio  (CAR+CH-): %.3f\n", ratio_ch_neg))
if (!is.nan(fisher_res_ch$estimate)) {
  cat(sprintf("Odds ratio (CD4 enrichment in CAR+CH+ vs CAR+CH-): %.3f\n", fisher_res_ch$estimate))
}


###### compare CAR-CH+ vs CAR-CH-

library(dplyr)
library(tidyr)

# assume `dat` already created with CD (CD4/CD8) and CAR_status as before

# Restrict to CAR- cells only
dat_car_neg <- dat %>%
  filter(CAR_status == "CAR-") %>%
  mutate(
    CH_status = ifelse(Variant_Type != "WT", "CH+", "CH-")
  ) %>%
  filter(CH_status %in% c("CH+","CH-"))

# Build contingency table (CH x CD subset)
contingency_car_neg <- dat_car_neg %>%
  count(CH_status, CD) %>%
  pivot_wider(names_from = CD, values_from = n, values_fill = 0) %>%
  arrange(match(CH_status, c("CH+","CH-"))) %>%
  as.data.frame()

rownames(contingency_car_neg) <- contingency_car_neg$CH_status
cont_2x2_car_neg <- as.matrix(contingency_car_neg[, c("CD4","CD8")])

contingency_car_neg
cont_2x2_car_neg

# Fisher's exact test
fisher_res_car_neg <- fisher.test(cont_2x2_car_neg)

# Ratios
ratio_ch_pos <- cont_2x2_car_neg["CH+","CD4"] / max(1, cont_2x2_car_neg["CH+","CD8"])
ratio_ch_neg <- cont_2x2_car_neg["CH-","CD4"] / max(1, cont_2x2_car_neg["CH-","CD8"])

cat("\nCounts (CAR-CH+ vs CAR-CH- by CD subset):\n")
print(cont_2x2_car_neg)

cat("\nFisher's exact test (CD4:CD8 ratio difference in CAR- cells by CH status):\n")
print(fisher_res_car_neg)

cat(sprintf("\nCD4:CD8 ratio  (CAR-CH+): %.3f\n", ratio_ch_pos))
cat(sprintf("CD4:CD8 ratio  (CAR-CH-): %.3f\n", ratio_ch_neg))
if (!is.nan(fisher_res_car_neg$estimate)) {
  cat(sprintf("Odds ratio (CD4 enrichment in CAR-CH+ vs CAR-CH-): %.3f\n", fisher_res_car_neg$estimate))
}


###########

per_patient <- dat %>%
  filter(CAR_status == "CAR-") %>%
  mutate(CH_status = ifelse(Variant_Type != "WT", "CH+", "CH-")) %>%
  group_by(patient_id, CH_status, CD) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = CD, values_from = n, values_fill = 0) %>%
  mutate(ratio = CD4 / pmax(1, CD8))

# compare ratios by CH status
wilcox.test(ratio ~ CH_status, data = per_patient)

############

library(lme4)

# only CD4/CD8 cells
dat_car_neg <- dat %>%
  filter(CAR_status == "CAR-") %>%
  mutate(
    CH_status = ifelse(Variant_Type != "WT", "CH+", "CH-"),
    outcome_CD4 = ifelse(CD == "CD4", 1, 0)
  )

model <- glmer(outcome_CD4 ~ CH_status + (1|patient_id),
               data = dat_car_neg, family = binomial)

summary(model)
