# Copy outside the public repo and fill the sections for the scripts you need.
# Paths resolve relative to this file. NULL means unconfigured; optional paths may stay NULL.
# Use character(0) for an intentionally empty selection.
list(
  repertoire = list(
    # tcr_bcr_internal_timepoint_map <- study_setting("tcr_bcr_internal_timepoint_map")
    "tcr_bcr_internal_timepoint_map" = NULL,
    # remove  <- study_setting("remove")
    "remove" = NULL,
    # setwd(study_path("study_path_001"))
    "study_path_001" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # combined <- readRDS(study_path("input_all_tcr_combined_contigs"))
    "input_all_tcr_combined_contigs" = NULL,
    # study_setting("study_value_001"),
    "study_value_001" = NULL,
    # vdj_b_input_candidates <- study_path("input_bcr_counts")
    "input_bcr_counts" = NULL
  ),
  figures = list(
    # mds_pt <- study_setting("mds_pt")
    "mds_pt" = NULL,
    # readme <- study_setting("readme")
    "readme" = NULL,
    # base_scRNA_dir <- study_path("base_scRNA_dir")
    "base_scRNA_dir" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # combined <- readRDS(study_path("input_all_tcr_combined_contigs"))
    "input_all_tcr_combined_contigs" = NULL,
    # Marrow_merge_idents <- readRDS(study_path("input_final_scrnaseq"))
    "input_final_scrnaseq" = NULL,
    # CAR_marrow_only  <- read_rds_optional(study_path("input_car_marrow_only", optional = TRUE),  "CAR_marrow_only")
    "input_car_marrow_only" = NULL,
    # CTRL_marrow_only <- read_rds_optional(study_path("input_ctrl_marrow_only", optional = TRUE), "CTRL_marrow_only")
    "input_ctrl_marrow_only" = NULL,
    # if (file.exists(study_path("input_tcr_diversity", optional = TRUE))) {
    "input_tcr_diversity" = NULL,
    # if (file.exists(study_path("input_b_cell_proportions", optional = TRUE))) {
    "input_b_cell_proportions" = NULL
  ),
  composition = list(
    # base_scRNA_dir <- study_path("base_scRNA_dir")
    "base_scRNA_dir" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # Marrow_merge_idents <- readRDS(study_path("input_final_scrnaseq"))
    "input_final_scrnaseq" = NULL
  )
)
