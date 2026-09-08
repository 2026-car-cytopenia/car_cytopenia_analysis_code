# CAR cytopenia analysis code

Analysis and plotting code for "CAR19 therapy drives expansion of clonal hematopoiesis and associated cytopenias" (Hamilton et al.).

```text
preprocessing/
  example_inputs/                   Five small, bundled example tables
  01_filter_variants.py             Filter monitoring and matching depth tables
  02_prepare_analysis_tables.py     Combine monitoring and prepare enrichment input
  common.py                        Shared table helpers
analysis/
  variant_enrichment/              Enrichment statistics and volcano plots
  ch_plots/                        CH calls, correlations and clinical CH figures
  preTx_vs_preCAR/                  Prevalence and baseline sample-depth summaries
  CAR_vs_autoHCT/                   Infection and clinical comparisons
  scRNA/                           Expression, composition and repertoire figures
  scDNA/                           Mutation tracing and embedding figures
  misc/                            Marrow, oncoprint and flow figures
  _shared/                         External-settings helper for original R scripts
run_preprocessing.sh               Runs the two preprocessing steps only
```

## Preprocess the example

Requires Python 3.9+ with pandas and NumPy:

```bash
python3 -m pip install -r requirements.txt
./run_preprocessing.sh --example
```

This reads the five files in `preprocessing/example_inputs/` and writes
`outputs/example/filtered/` and `outputs/example/processed/`.
The example includes synonymous annotations, repeated annotations, a gene outside
the supplied list, and a common indel so the filtering steps are visible.

For your own concatenated source files, supply five paths and an output directory:

```bash
./run_preprocessing.sh snv.csv indel.csv snv_depth.csv indel_depth.csv genes.tsv outputs/cohort
```

See [preprocessing/README.md](preprocessing/README.md) for columns and output files.

## Run an analysis separately

Each analysis directory has its own README and commands. For example:

```bash
Rscript -e 'install.packages(c("dplyr", "readr", "ggplot2", "ggrepel"))'
Rscript analysis/variant_enrichment/enrichment.R \
  outputs/example/processed/enrichment_input.csv outputs/example/enrichment
```

CH discovery requires caller-positive variants; prevalence additionally requires
cohort membership including participants without CH. Clinical and single-cell
analyses require their own study inputs.

The original clinical and single-cell R scripts accept a private settings file:

```bash
Rscript analysis/scRNA/figures.R private/settings.R
```

Copy the `settings.example.R` from that analysis directory and fill the selected
script's section. Paths, sample selections and crosswalks are supplied there.