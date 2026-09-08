# Copy outside the public repo and fill the sections for the scripts you need.
# Paths resolve relative to this file. NULL means unconfigured; optional paths may stay NULL.
# Use character(0) for an intentionally empty selection.
list(
  clinical_ch_figures = list(
    # nick_list <- study_setting("nick_list")
    "nick_list" = NULL,
    # all_pre_ch <- study_setting("all_pre_ch")
    "all_pre_ch" = NULL,
    # tp53_example_ids <- study_setting("tp53_example_ids")
    "tp53_example_ids" = NULL,
    # analysis_dir <- study_path("analysis_dir")
    "analysis_dir" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # all_variants_raw <- safe_read_delim(study_path("input_all_snv_indel_monitoring"))
    "input_all_snv_indel_monitoring" = NULL,
    # pre_inf_variants_raw <- safe_read_delim(study_path("input_v16_raw_calls"))
    "input_v16_raw_calls" = NULL,
    # just_CH_mutations_raw <- safe_read_delim(study_path("input_v16_justch_mutations_snvs_indels"))
    "input_v16_justch_mutations_snvs_indels" = NULL,
    # bm_variants_raw <- safe_read_delim(study_path("input_bm_oncoprint_snvs_4_update_20241011_with_mm_format_for_merge"))
    "input_bm_oncoprint_snvs_4_update_20241011_with_mm_format_for_merge" = NULL,
    # BM_vs_raw <- safe_read_delim(study_path("input_bm_vs_cfdna_vs_pbl"))
    "input_bm_vs_cfdna_vs_pbl" = NULL,
    # tMN_molecular_raw <- safe_read_delim(study_path("input_tmn_molecular_trace"))
    "input_tmn_molecular_trace" = NULL,
    # preCH_raw <- safe_read_delim(study_path("input_v15_precar_ch_list"))
    "input_v15_precar_ch_list" = NULL,
    # cytopenia_raw <- safe_read_delim(study_path("input_cytopenia_list"))
    "input_cytopenia_list" = NULL,
    # pre_CH_data_raw <- safe_read_delim(study_path("input_pre_car_ch_workbook"))
    "input_pre_car_ch_workbook" = NULL,
    # tMN_car_raw <- read_tmn_car_aai_file(study_path("input_tmn_car"))
    "input_tmn_car" = NULL,
    # expansion_data_raw <- safe_read_delim(study_path("input_ch_expand_eval"))
    "input_ch_expand_eval" = NULL,
    # tmn_CAR_raw <- safe_read_delim(study_path("input_tmn_axicel_depth_workbook"))
    "input_tmn_axicel_depth_workbook" = NULL,
    # match_var_raw <- if (file.exists(study_path("input_pretx_precar_match", optional = TRUE))) safe_read_delim(study_path("input_pretx_precar_match", optional = TRUE)) else tibble()
    "input_pretx_precar_match" = NULL,
    # filter(patient_id != study_setting("study_value_001")) %>%
    "study_value_001" = NULL
  )
)
