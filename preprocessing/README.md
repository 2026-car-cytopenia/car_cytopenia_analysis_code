# Preprocessing

Start with concatenated monitoring and depth tables. Alignment, variant calling
and concatenation from individual caller files happen upstream of this workflow.
Run `../run_preprocessing.sh` from this directory, or `./run_preprocessing.sh` from
the repository root. `--example` uses the five bundled synthetic inputs.

## Inputs

| File | Required columns |
|---|---|
| `snv_monitoring.csv` | `CHR, POS, TUMOR, REF, GENE, TYPE, SAMPLE, AF, sample_id, patient_id, sample_day, sample_type, selector, progression` |
| `indel_monitoring.csv` | Same columns as SNV monitoring |
| `snv_depth.csv` | `CHR, POS, TUMOR, REF, sample_id, Depth` |
| `indel_depth.csv` | Same columns as SNV depth |
| `protein_coding_genes.tsv` | `symbol` |

Monitoring files contain sample–variant–annotation measurements. `TUMOR` means
alternate allele. `AF` is a percentage: `0.5` means 0.5%. `SAMPLE` is the original
measurement identifier; `sample_id` links sample metadata. Preserve their mapping
across the input files. `sample_day` is relative to treatment; `progression` is
0/1 (`Progression` is also accepted). Use the same coordinate/allele conventions
across all tables. SNV `Depth` is total depth; indel `Depth` is alternate-supporting
depth, following the supplied source workflow.

The example also includes `detailed_sample_type` for specimen comparisons. Every
record and sample label is invented. Gene symbols form an illustrative reference
subset; coordinates, gene assignments and measurements are fabricated, and the
five-gene list is not a complete protein-coding reference for real analyses.

## Ordered steps

1. **`01_filter_variants.py`** retains supplied protein-coding genes, removes
   ncRNA/intergenic annotations, and removes indel annotations with AF >0.5% in
   >80% of participants. It removes identical SNV rows and restricts depth tables
   to retained sample–allele keys, preserving the notebook's behavior.
2. **`02_prepare_analysis_tables.py`** concatenates filtered SNVs/indels and
   prepares the enrichment subset: excludes MRD selectors, takes the first
   semicolon-separated gene, removes synonymous annotations, then retains the
   first row per `CHR, POS, TUMOR, REF, SAMPLE`. The final enrichment subset also
   excludes genes containing `LINC`.

These operations retain the source's filtering order. The enrichment synonymous
filter does not remove intronic/UTR annotations; the separate CH classification
uses its original functional-consequence filter. The >0.5% endpoint threshold is
applied later by the enrichment R script.

## Outputs

| Directory | Files |
|---|---|
| `filtered/` | `snv_filtered.csv`, `indel_filtered.csv`, matching `snv_depth_filtered.csv` and `indel_depth_filtered.csv`, `excluded_common_indels.csv`, `filter_counts.json` |
| `processed/` | `monitoring.csv`, `snv_depth.csv`, `indel_depth.csv`, `enrichment_input.csv` |

`monitoring.csv` preserves the filtered annotations for downstream CH processing.
`enrichment_input.csv` applies the enrichment-specific annotation collapse. The
analysis scripts consume these tables independently; preprocessing does not run
plots or invent CH calls, cohort membership, clinical outcomes or capture regions.

Each numbered script has `--help` and can also be run directly.
