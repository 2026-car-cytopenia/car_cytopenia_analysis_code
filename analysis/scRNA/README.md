# Single-cell RNA and receptor analyses

Input paths, cohort selections and sample crosswalks are supplied through `settings.example.R`.
Copy that file to a private location and fill the section for the chosen script.
All paths resolve against your settings file. Required `NULL` entries fail with
a named error; optional paths may remain `NULL`. Use strings for paths/selections,
character vectors for participant lists, and data frames for crosswalks. Comments
in the template show the consuming source expression.

### `repertoire.R`

Combined receptor-contig RDS, BCR/Ig count table and sample/timepoint crosswalk.

```bash
Rscript analysis/scRNA/repertoire.R private/settings.R
```

Packages used (including conditional analyses): `dplyr`, `ggplot2`, `ggpubr`, `purrr`, `readr`, `scRepertoire`, `scales`, `stringr`, `tibble`, `tidyr`.

File settings: `input_all_tcr_combined_contigs`, `input_bcr_counts`.

### `figures.R`

Annotated Seurat RDS and combined receptor-contig RDS; optional marrow RDS objects, TCR diversity and B-cell proportion tables.

```bash
Rscript analysis/scRNA/figures.R private/settings.R
```

Packages used (including conditional analyses): `Seurat`, `dplyr`, `fgsea`, `ggh4x`, `ggplot2`, `ggprism`, `ggpubr`, `ggrastr`, `ggrepel`, `msigdbr`, `patchwork`, `purrr`, `readr`, `scales`, `stringr`, `tibble`, `tidyr`.

File settings: `input_all_tcr_combined_contigs`, `input_final_scrnaseq`, `input_car_marrow_only`, `input_ctrl_marrow_only`, `input_tcr_diversity`, `input_b_cell_proportions`.

### `composition.R`

Annotated Seurat RDS with the sample, cohort and cell-type metadata consumed by the source script.

```bash
Rscript analysis/scRNA/composition.R private/settings.R
```

Packages used (including conditional analyses): `Seurat`, `dplyr`, `ggplot2`, `ggprism`, `glmmTMB`, `purrr`, `readr`, `scales`, `speckle`, `stringr`, `tibble`, `tidyr`.

File settings: `input_final_scrnaseq`.

Install the packages required by the selected script in your R environment.
Some are Bioconductor packages; the source retains its original optional-package
checks and Arial/rasterization settings. Numerical outputs could not be tested
without these inputs; R syntax and preservation of the calculation expressions
were checked. These scripts write local figures/tables and do not upload data.
