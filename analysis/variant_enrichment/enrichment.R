#!/usr/bin/env Rscript
# Gene-level expansion using the original R Wilcoxon/scaling implementation.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

get_var_fold_change <- function(first_patient_df, last_patient_df, min_vaf = 0.5) {
  intersect <- suppressWarnings(inner_join(first_patient_df, last_patient_df,
    by = c("CHR", "POS", "TUMOR", "REF", "GENE"), suffix = c("_first", "_last")))
  intersect <- intersect %>% filter(AF_first > 0, AF_last > 0,
                                   AF_first > min_vaf | AF_last > min_vaf)
  intersect$logFC <- log(intersect$AF_last / intersect$AF_first)
  intersect$zscore <- as.numeric(scale(intersect$logFC))
  intersect %>% filter(is.finite(zscore))
}

compute_enrichment <- function(df, responders, min_vaf = 0.5) {
  pairs <- list()
  for (patient in responders) {
    d <- df %>% filter(patient_id == patient, sample_type == "cfDNA")
    if (!nrow(d) || anyNA(d$sample_day)) next
    first_day <- min(d$sample_day)
    last_day <- max(d$sample_day)
    if (first_day < 1 && last_day >= 25) {
      pairs[[length(pairs) + 1]] <- get_var_fold_change(
        filter(d, sample_day == first_day), filter(d, sample_day == last_day), min_vaf)
    }
  }
  pair_data <- bind_rows(pairs)
  if (!nrow(pair_data)) {
    return(list(pairs = pair_data, stats = tibble(gene = character(), meanzscore = numeric(),
      stat = numeric(), pval = numeric(), pval_corrected = numeric())))
  }
  result <- bind_rows(lapply(unique(pair_data$GENE), function(gene) {
    x <- pair_data$zscore[pair_data$GENE == gene]
    y <- pair_data$zscore[pair_data$GENE != gene]
    if (!length(y)) return(tibble(gene = gene, meanzscore = mean(x), stat = NA_real_, pval = NA_real_))
    test <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided"))
    tibble(gene = gene, meanzscore = mean(x), stat = unname(test$statistic), pval = test$p.value)
  }))
  result$pval_corrected <- p.adjust(result$pval, method = "fdr")
  list(pairs = pair_data, stats = result)
}

plot_volcano <- function(data, title) {
  cols <- c(DNMT3A = "#0072B2", TET2 = "#009E73", TP53 = "#C43A31", PPM1D = "#CC79A7",
            BCOR = "#E69F00", ATM = "#56B4E9", EZH2 = "#000000", IDH1 = "#F0E442")
  data <- data %>% mutate(color_label = ifelse(pval_corrected < 0.05 & gene %in% names(cols), gene, "Other"))
  ggplot(data, aes(meanzscore, -log10(pval_corrected))) +
    geom_point(aes(color = color_label)) + scale_color_manual(values = c(cols, Other = "black")) +
    theme_classic(base_family = "sans") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    labs(title = title, x = "Mean Z-Score", y = expression(-log[10](FDR~Corrected~P~Value))) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_text_repel(data = data %>% filter(is.finite(pval_corrected), pval_corrected > 0) %>%
                      slice_max(-log10(pval_corrected), n = 8),
                    aes(label = gene), box.padding = 0.6, point.padding = 0.5,
                    segment.color = "grey50", seed = 1) +
    scale_x_continuous(limits = c(-2.5, 2.5))
}

run_enrichment <- function(input, out, min_vaf = 0.5) {
  if (!is.finite(min_vaf) || min_vaf < 0 || min_vaf >= 100) stop("min_vaf must be a percentage in [0,100)")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  df <- read_csv(input, show_col_types = FALSE)
  required <- c("CHR", "POS", "TUMOR", "REF", "GENE", "SAMPLE", "AF", "patient_id",
                "sample_day", "sample_type", "progression", "variant_class")
  stopifnot(all(required %in% names(df)))
  if (anyDuplicated(df[c("CHR", "POS", "TUMOR", "REF", "SAMPLE")])) {
    stop("Input contains repeated sample–allele entries; run Python preprocessing first")
  }
  responders <- unique(df$patient_id[df$progression == 0 & !is.na(df$progression)])
  for (kind in c("all", "snv", "indel")) {
    unlink(file.path(out, c(paste0("volcano_", kind, c(".png", ".pdf")), paste0("variant_pairs_", kind, ".csv"))))
    d <- if (kind == "all") df else filter(df, variant_class == kind)
    result <- compute_enrichment(d, responders, min_vaf)
    write_csv(result$stats, file.path(out, paste0("enrichment_", kind, ".csv")))
    if (ncol(result$pairs)) write_csv(result$pairs, file.path(out, paste0("variant_pairs_", kind, ".csv")))
    if (nrow(result$stats)) {
      p <- plot_volcano(result$stats, if (kind == "all") "SNV and Indels" else toupper(kind))
      ggsave(file.path(out, paste0("volcano_", kind, ".png")), p, width = 5, height = 4, dpi = 300)
      ggsave(file.path(out, paste0("volcano_", kind, ".pdf")), p, width = 5, height = 4)
    }
    message(kind, ": ", nrow(result$pairs), " paired measurements; ", nrow(result$stats), " genes")
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript enrichment.R enrichment_input.csv output_dir [min_vaf_percent]")
  run_enrichment(args[1], args[2], if (length(args) >= 3) as.numeric(args[3]) else 0.5)
}
