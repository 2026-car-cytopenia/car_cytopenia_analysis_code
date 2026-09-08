# Single-cell DNA analyses

Input paths, cohort selections and sample crosswalks are supplied through `settings.example.R`.
Copy that file to a private location and fill the section for the chosen script.
All paths resolve against your settings file. Required `NULL` entries fail with
a named error; optional paths may remain `NULL`. Use strings for paths/selections,
character vectors for participant lists, and data frames for crosswalks. Comments
in the template show the consuming source expression.

### `after_omiq.R`

Merged OMIQ/t-SNE table and mutation-frequency comparison table, plus configured specimen selections.

```bash
Rscript analysis/scDNA/after_omiq.R private/settings.R
```

Packages used (including conditional analyses): `dplyr`, `ggh4x`, `ggplot2`, `ggprism`, `ggpubr`, `ggrastr`, `lme4`, `patchwork`, `readxl`, `scales`, `stringr`, `survival`, `tableone`, `tidyr`, `tidyverse`.

File settings: `input_percent_vs_stamp_scdna`, `input_merged_tsne`.

### `mutation_traces.R`

Merged scDNA and PBMC/marrow tables with molecular/clinical context and an explicit sample crosswalk.

```bash
Rscript analysis/scDNA/mutation_traces.R private/settings.R
```

Packages used (including conditional analyses): `ComplexHeatmap`, `broom`, `circlize`, `dplyr`, `forcats`, `ggh4x`, `ggplot2`, `ggrastr`, `patchwork`, `purrr`, `readr`, `scales`, `stringr`, `tibble`, `tidyr`, `tidyverse`.

File settings: `input_file_all`, `input_scdna_pbmc_bm_merged`.

Install the packages required by the selected script in your R environment.
Some are Bioconductor packages; the source retains its original optional-package
checks and Arial/rasterization settings. Numerical outputs could not be tested
without these inputs; R syntax and preservation of the calculation expressions
were checked. These scripts write local figures/tables and do not upload data.
