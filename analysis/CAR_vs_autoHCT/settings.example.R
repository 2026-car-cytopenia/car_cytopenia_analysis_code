# Copy outside the public repo and fill the sections for the scripts you need.
# Paths resolve relative to this file. NULL means unconfigured; optional paths may stay NULL.
# Use character(0) for an intentionally empty selection.
list(
  infection = list(
    # PATH_INF    <- study_path("PATH_INF")
    "PATH_INF" = NULL,
    # PATH_MATCH  <- study_path("PATH_MATCH")
    "PATH_MATCH" = NULL,
    # PATH_CYTO   <- study_path("PATH_CYTO")
    "PATH_CYTO" = NULL,
    # out_dir <- study_path("out_dir")
    "out_dir" = NULL,
    # bad_pattern <- study_setting("study_value_001")
    "study_value_001" = NULL
  )
)
