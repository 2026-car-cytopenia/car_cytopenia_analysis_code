# Copy outside the public repo and fill the sections for the scripts you need.
# Paths resolve relative to this file. NULL means unconfigured; optional paths may stay NULL.
# Use character(0) for an intentionally empty selection.
list(
  after_omiq = list(
    # tsne_vector_raster_patterns <- study_setting("tsne_vector_raster_patterns")
    "tsne_vector_raster_patterns" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # setwd(study_path("study_path_001"))
    "study_path_001" = NULL,
    # pct_vs_stamp <- read_delim(study_path("input_percent_vs_stamp_scdna"))
    "input_percent_vs_stamp_scdna" = NULL,
    # merged_scDNA <- read_delim(study_path("input_merged_tsne"))
    "input_merged_tsne" = NULL,
    # message(study_setting("study_value_001"))
    "study_value_001" = NULL,
    # merged_scDNA$Sample_ID <- ifelse(merged_scDNA$Sample_ID == study_setting("study_value_002"), "Non-CAR", merged_scDNA$Sample_ID)
    "study_value_002" = NULL,
    # merged_scDNA$Sample_ID <- ifelse(merged_scDNA$Sample_ID == study_setting("study_value_003"), "CAR", merged_scDNA$Sample_ID)
    "study_value_003" = NULL,
    # merged_scDNA$patient_id <- study_setting("study_value_004")
    "study_value_004" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_005"), "Control", merged_scDNA$patient_id)
    "study_value_005" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_006"), "Product_1", merged_scDNA$patient_id)
    "study_value_006" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_007"), "Product_2", merged_scDNA$patient_id)
    "study_value_007" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_008"), study_setting("study_value_009"), merged_scDNA$patient_id)
    "study_value_008" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_008"), study_setting("study_value_009"), merged_scDNA$patient_id)
    "study_value_009" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_010"), study_setting("study_value_011"), merged_scDNA$patient_id)
    "study_value_010" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_010"), study_setting("study_value_011"), merged_scDNA$patient_id)
    "study_value_011" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_012"), study_setting("study_value_013"), merged_scDNA$patient_id)
    "study_value_012" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_012"), study_setting("study_value_013"), merged_scDNA$patient_id)
    "study_value_013" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_014"), study_setting("study_value_015"), merged_scDNA$patient_id)
    "study_value_014" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_014"), study_setting("study_value_015"), merged_scDNA$patient_id)
    "study_value_015" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_016"), study_setting("study_value_017"), merged_scDNA$patient_id)
    "study_value_016" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_016"), study_setting("study_value_017"), merged_scDNA$patient_id)
    "study_value_017" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_018"), study_setting("study_value_019"), merged_scDNA$patient_id)
    "study_value_018" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_018"), study_setting("study_value_019"), merged_scDNA$patient_id)
    "study_value_019" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_020"), study_setting("study_value_021"), merged_scDNA$patient_id)
    "study_value_020" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_020"), study_setting("study_value_021"), merged_scDNA$patient_id)
    "study_value_021" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_022"), study_setting("study_value_004"), merged_scDNA$patient_id)
    "study_value_022" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_023"), study_setting("study_value_024"), merged_scDNA$patient_id)
    "study_value_023" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_023"), study_setting("study_value_024"), merged_scDNA$patient_id)
    "study_value_024" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_025"), study_setting("study_value_026"), merged_scDNA$patient_id)
    "study_value_025" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_025"), study_setting("study_value_026"), merged_scDNA$patient_id)
    "study_value_026" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_027"), study_setting("study_value_004"), merged_scDNA$patient_id)
    "study_value_027" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_028"), study_setting("study_value_029"), merged_scDNA$patient_id)
    "study_value_028" = NULL,
    # merged_scDNA$patient_id <- ifelse(merged_scDNA$Variant_Call == study_setting("study_value_028"), study_setting("study_value_029"), merged_scDNA$patient_id)
    "study_value_029" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_030"), "TP53_G266V", merged_scDNA$Variant_Specific)
    "study_value_030" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_031"), "TP53_F270S", merged_scDNA$Variant_Specific)
    "study_value_031" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_032"), "TP53_F270S", merged_scDNA$Variant_Specific)
    "study_value_032" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_033"), "DNMT3A_R729W", merged_scDNA$Variant_Specific)
    "study_value_033" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_034"), "DNMT3A_splicing", merged_scDNA$Variant_Specific)
    "study_value_034" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_035"), "DNMT3A_L547H", merged_scDNA$Variant_Specific)
    "study_value_035" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_036"), "DNMT3A_M801V", merged_scDNA$Variant_Specific)
    "study_value_036" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_037"), "DNMT3A_C497R", merged_scDNA$Variant_Specific)
    "study_value_037" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_038"), "TET2_H1676fs", merged_scDNA$Variant_Specific)
    "study_value_038" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_039"), "TET2_W1003*", merged_scDNA$Variant_Specific)
    "study_value_039" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_040"), "DNMT3A_splicing", merged_scDNA$Variant_Specific)
    "study_value_040" = NULL,
    # merged_scDNA$Variant_Specific <- ifelse(merged_scDNA$Intermediate == study_setting("study_value_041"), "DNMT3A_T862I", merged_scDNA$Variant_Specific)
    "study_value_041" = NULL,
    # ggsave(out_file(study_setting("study_value_042")), width = 16, height = 4.06, units = "in")
    "study_value_042" = NULL,
    # ggsave(out_file(study_setting("study_value_043")), width = 16, height = 4.06, units = "in", device = cairo_pdf)
    "study_value_043" = NULL,
    # ggsave(out_file(study_setting("study_value_044")), width = 15, height = 4.06, units = "in")
    "study_value_044" = NULL,
    # ggsave(out_file(study_setting("study_value_045")), width = 15, height = 4.06, units = "in", device = cairo_pdf)
    "study_value_045" = NULL,
    # ggsave(out_file(study_setting("study_value_046")), width = 7.5, height = 4.06, units = "in")
    "study_value_046" = NULL,
    # ggsave(out_file(study_setting("study_value_047")), width = 7.5, height = 4.06, units = "in", device = cairo_pdf)
    "study_value_047" = NULL,
    # ggsave(out_file(study_setting("study_value_048")), width = 15, height = 8.06, units = "in")
    "study_value_048" = NULL,
    # ggsave(out_file(study_setting("study_value_049")), width = 15, height = 8.06, units = "in", device = cairo_pdf)
    "study_value_049" = NULL,
    # ggsave(out_file(study_setting("study_value_050")), width = 7.5, height = 4.06, units = "in")
    "study_value_050" = NULL,
    # ggsave(out_file(study_setting("study_value_051")), width = 7.5, height = 4.06, units = "in", device = cairo_pdf)
    "study_value_051" = NULL
  ),
  mutation_traces = list(
    # tsne_vector_raster_patterns <- study_setting("tsne_vector_raster_patterns")
    "tsne_vector_raster_patterns" = NULL,
    # id_key <- study_setting("id_key")
    "id_key" = NULL,
    # setwd(study_path("study_path_001"))
    "study_path_001" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # input_file_all <- study_path("input_file_all")
    "input_file_all" = NULL,
    # study_path("input_scdna_pbmc_bm_merged"),
    "input_scdna_pbmc_bm_merged" = NULL,
    # filter(Patient_Paper != study_setting("study_value_001"))
    "study_value_001" = NULL,
    # patient_of_interest <- study_setting("study_value_002")
    "study_value_002" = NULL,
    # cat(study_setting("study_value_003"))
    "study_value_003" = NULL,
    # cat(study_setting("study_value_004"))
    "study_value_004" = NULL,
    # out_file(study_setting("study_value_005")),
    "study_value_005" = NULL
  )
)
