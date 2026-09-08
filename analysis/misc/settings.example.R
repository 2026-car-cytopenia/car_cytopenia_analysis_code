# Copy outside the public repo and fill the sections for the scripts you need.
# Paths resolve relative to this file. NULL means unconfigured; optional paths may stay NULL.
# Use character(0) for an intentionally empty selection.
list(
  marrow_oncoprint_flow = list(
    # patient_order_tx <- study_setting("patient_order_tx")
    "patient_order_tx" = NULL,
    # base_oncoprint_dir <- study_path("base_oncoprint_dir")
    "base_oncoprint_dir" = NULL,
    # base_flow_dir      <- study_path("base_flow_dir")
    "base_flow_dir" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # bad_pattern <- study_setting("study_value_001")
    "study_value_001" = NULL,
    # SNV <- safe_read_tsv(study_path("input_bm_oncoprint_snvs_4_update_20241011_with_mm"))
    "input_bm_oncoprint_snvs_4_update_20241011_with_mm" = NULL,
    # MDS_AML <- safe_read_tsv(study_path("input_mds_aml_clinical"))
    "input_mds_aml_clinical" = NULL,
    # MDS_AML_auto_v_CAR <- safe_read_tsv(study_path("input_mds_aml_auto_vs_car"))
    "input_mds_aml_auto_vs_car" = NULL,
    # clin <- safe_read_delim(study_path("input_marrow_clinical_5"))
    "input_marrow_clinical_5" = NULL,
    # tmn_risk <- safe_read_delim(study_path("input_tmn_risk_variables"))
    "input_tmn_risk_variables" = NULL,
    # flow <- read.delim(study_path("input_flow"), stringsAsFactors = FALSE)
    "input_flow" = NULL,
    # car_cd19 <- read.delim(study_path("input_car_cd19"), stringsAsFactors = FALSE)
    "input_car_cd19" = NULL
  )
)
