# CH calls, correlations and figures

Run each script independently from the repository root. Python scripts use pandas
and NumPy. Inputs beyond the five preprocessing examples are supplied by the user.

## CH discovery

```bash
python3 analysis/ch_plots/prepare_ch.py \
  --monitoring outputs/cohort/processed/monitoring.csv \
  --snv-depth outputs/cohort/processed/snv_depth.csv \
  --indel-depth outputs/cohort/processed/indel_depth.csv \
  --calls private/de_novo_calls.csv --out outputs/cohort/ch
```

The calls file contains caller-positive `CHR,POS,TUMOR,REF,sample_id` entries;
monitoring measurements alone do not establish these calls. Optional `--cohort
pretx` selects untreated discovery. Optional `--members` accepts
`sample_id,patient_id`; `--bed` restricts baseline CH to shared assay territory.
`--position-base 0` reproduces the original notebook's BED query, or use 1 for
explicitly one-based coordinates. These extra inputs are not bundled.

Writes `raw_calls.csv`, `ch_monitoring.csv`, `ch_baseline_all_regions.csv` and
`ch_baseline.csv`. Uses the original functional/CH-gene rules and at least five
supporting reads. Run the prevalence script in `analysis/preTx_vs_preCAR/` next
when a complete cohort membership file is available.

## Longitudinal correlations

```bash
python3 analysis/ch_plots/correlations.py \
  --monitoring outputs/cohort/ch/ch_monitoring.csv --out outputs/cohort/correlations
```

Requires `detailed_sample_type` in the monitoring table, or an optional
`--metadata private/sample_metadata.csv` with `sample_id,detailed_sample_type`.
Labels distinguish cfDNA, PBMC and PDWB. The script retains the original seven-gene
subset, exact-day pairing (at least three shared days) and PDWB preference for PBL.
Writes four specimen-pair correlation CSVs. The input defines the variants studied;
use CH monitoring when reproducing the CH correlations.

## Summary plots

```bash
Rscript analysis/ch_plots/plot_summaries.R outputs/cohort/correlations outputs/cohort/plots
```

Reads correlation CSVs, `prevalence.csv` and/or `mutation_counts.csv` present in
the supplied directory; writes PNG/PDFs. To compare cohorts, concatenate their
prevalence/count tables while retaining the `source` labels. Requires `ggplot2`,
`readr`, `tidyr`, `dplyr`, `colorspace` and `scales`.

## Original R scripts

These scripts retain the source calculations and figure routines. Input paths,
cohort selections and sample crosswalks are supplied through `settings.example.R`.
Copy that file to a private location and fill the section for the chosen script.
All paths resolve against your settings file. Required `NULL` entries fail with
a named error; optional paths may remain `NULL`. Use strings for paths/selections,
character vectors for participant lists, and data frames for crosswalks. Comments
in the template show the consuming source expression. No study inputs are bundled.

### `clinical_ch_figures.R`

CH monitoring/call tables, marrow comparisons, molecular tracing, cytopenia and outcome tables/workbooks, and optional matching data.

```bash
Rscript analysis/ch_plots/clinical_ch_figures.R private/settings.R
```

Packages used (including conditional analyses): `broom`, `dplyr`, `ggplot2`, `ggpubr`, `patchwork`, `purrr`, `readr`, `scales`, `stringr`, `survival`, `survminer`, `tibble`, `tidyr`, `tidyverse`.

File settings: `input_all_snv_indel_monitoring`, `input_v16_raw_calls`, `input_v16_justch_mutations_snvs_indels`, `input_bm_oncoprint_snvs_4_update_20241011_with_mm_format_for_merge`, `input_bm_vs_cfdna_vs_pbl`, `input_tmn_molecular_trace`, `input_v15_precar_ch_list`, `input_cytopenia_list`, `input_pre_car_ch_workbook`, `input_tmn_car`, `input_ch_expand_eval`, `input_tmn_axicel_depth_workbook`, `input_pretx_precar_match`.

Install the packages required by the selected script in your R environment.
Some are Bioconductor packages; the source retains its original optional-package
checks and Arial/rasterization settings. Numerical outputs could not be tested
without these inputs; R syntax and preservation of the calculation expressions
were checked. These scripts write local figures/tables and do not upload data.
