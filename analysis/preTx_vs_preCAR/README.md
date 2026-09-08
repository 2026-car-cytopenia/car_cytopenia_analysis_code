# Baseline CH prevalence and sample depth

These are standalone extractions of the original prevalence and depth notebooks.
Run from the repository root with Python, pandas and NumPy.

```bash
python3 analysis/preTx_vs_preCAR/prevalence.py \
  --ch outputs/cohort/ch/ch_baseline.csv \
  --members private/cohort_members.csv --source 'Pre-CAR-T' --out outputs/cohort/prevalence
```

The CH input comes from `analysis/ch_plots/prepare_ch.py`. Membership requires
`sample_id,patient_id`, including eligible participants with no CH. These extra
inputs are not bundled. Use the selected/matched samples for your intended cohort.

Outputs: `prevalence.csv`, `patient_gene_depth.csv`, `mutation_counts.csv`.
Prevalence uses the original <2% and >=2% bins and earliest CH-positive baseline
PBL day. Per-gene bins are counted independently; `Any` gives high VAF precedence.
Denominators come from the membership file, including any additional controls.
For plots, pass the output directory to `analysis/ch_plots/plot_summaries.R`.

Optional coverage summary:

```bash
python3 analysis/preTx_vs_preCAR/sample_depths.py \
  --monitoring outputs/cohort/processed/monitoring.csv \
  --coverage private/coverage_percentiles.csv --members private/cohort_members.csv \
  --out outputs/cohort/sample_depths.csv
```

Coverage input requires `sample_id,type,percentile,depth`. The script selects
`type == "all"`, percentile 50, and baseline PBL samples. Variant-site depth cannot
substitute for coverage percentiles across assay positions. This input is not bundled.
