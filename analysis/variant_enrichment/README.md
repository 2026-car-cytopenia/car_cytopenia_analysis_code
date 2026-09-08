# Variant enrichment

Uses the corrected author R calculation: positive VAF at both endpoints, at least
one VAF strictly >0.5%, within-participant sample-SD z-scores, Wilcoxon tests and BH
correction. Compares earliest and latest responder cfDNA, with earliest day <1 and
latest day >=25. The Python preprocessing handles annotation filtering/collapse.

From the repository root:

```bash
Rscript analysis/variant_enrichment/enrichment.R \
  outputs/example/processed/enrichment_input.csv outputs/example/enrichment
```

Requires R with `dplyr`, `readr`, `ggplot2` and `ggrepel` (tested with R 4.4.1).
Optional third argument: minimum endpoint VAF in percent; default `0.5`.

Writes gene statistics, paired measurements and PNG/PDF volcano plots for combined
variants, SNVs and indels separately. Empty/undefined comparisons do not produce
a volcano plot. The bundled example demonstrates execution, not study results.
