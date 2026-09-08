# Marrow, oncoprint and flow analyses

Input paths, cohort selections and sample crosswalks are supplied through `settings.example.R`.
Copy that file to a private location and fill the section for the chosen script.
All paths resolve against your settings file. Required `NULL` entries fail with
a named error; optional paths may remain `NULL`. Use strings for paths/selections,
character vectors for participant lists, and data frames for crosswalks. Comments
in the template show the consuming source expression.

### `marrow_oncoprint_flow.R`

Marrow mutation tables, clinical/outcome tables, flow measurements and CAR/CD19 comparisons.

```bash
Rscript analysis/misc/marrow_oncoprint_flow.R private/settings.R
```

Packages used (including conditional analyses): `Cairo`, `ComplexHeatmap`, `brglm2`, `broom`, `circlize`, `cmprsk`, `dplyr`, `forcats`, `ggplot2`, `ggprism`, `ggpubr`, `logistf`, `readr`, `reshape2`, `rlang`, `stringr`, `survival`, `survminer`, `tableone`, `tibble`, `tidyr`.

File settings: `input_bm_oncoprint_snvs_4_update_20241011_with_mm`, `input_mds_aml_clinical`, `input_mds_aml_auto_vs_car`, `input_marrow_clinical_5`, `input_tmn_risk_variables`, `input_flow`, `input_car_cd19`.

Install the packages required by the selected script in your R environment.
Some are Bioconductor packages; the source retains its original optional-package
checks and Arial/rasterization settings. Numerical outputs could not be tested
without these inputs; R syntax and preservation of the calculation expressions
were checked. These scripts write local figures/tables and do not upload data.
