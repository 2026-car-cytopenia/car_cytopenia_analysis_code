#!/usr/bin/env Rscript
# Adapted from the supplied prevalence, CH distribution and correlation figures.
suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(tidyr); library(dplyr)
  library(colorspace); library(scales)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript plot_summaries.R processed_dir output_dir")
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
cols <- c(DNMT3A = "#0072B2", TET2 = "#009E73", TP53 = "#C43A31", PPM1D = "#CC79A7",
          BCOR = "#E69F00", ATM = "#56B4E9", EZH2 = "#000000", IDH1 = "#F0E442")
save_plot <- function(p, name, width = 10, height = 6) {
  for (ext in c("png", "pdf")) {
    ggsave(file.path(args[2], paste0(name, ".", ext)), p, width = width, height = height, dpi = 300)
  }
}

file <- file.path(args[1], "prevalence.csv")
if (file.exists(file)) {
  d <- read_csv(file, show_col_types = FALSE) %>% filter(!gene %in% c("U2AF1", "IDH1", "GNAS"))
  if (nrow(d)) {
    sources <- unique(d$source)
    # Cohort order follows configuration, with the original labels ordered first.
    sources <- unique(c(intersect(c("Pre-Treatment", "Pre-CAR-T"), sources), sources))
    labels <- d %>% distinct(source, n_total) %>% mutate(label = paste0(source, " (n = ", n_total, ")"))
    facet_labels <- setNames(labels$label, labels$source)
    sort_source <- if ("Pre-CAR-T" %in% sources) "Pre-CAR-T" else tail(sources, 1)
    gene_order <- d %>% filter(source == sort_source) %>% arrange((low_af + high_af) / n_total) %>% pull(gene)
    d <- d %>% pivot_longer(c(low_af, high_af), names_to = "mutation", values_to = "n") %>%
      mutate(percent = 100 * n / n_total, source = factor(source, levels = sources),
             gene = factor(gene, levels = gene_order),
             mutation_type = ifelse(mutation == "low_af", "AF < 2%", "AF >= 2%"),
             legend_label = ifelse(gene %in% names(cols), paste(gene, mutation_type), mutation_type))
    palette <- bind_rows(
      data.frame(legend_label = paste(names(cols), "AF < 2%"), color = lighten(unname(cols), .3)),
      data.frame(legend_label = paste(names(cols), "AF >= 2%"), color = darken(unname(cols), .3)),
      data.frame(legend_label = c("AF < 2%", "AF >= 2%"), color = c("lightgray", "darkgray")))
    d$legend_label <- factor(d$legend_label, levels = palette$legend_label)
    p <- ggplot(d, aes(y = gene, x = percent, fill = legend_label)) + geom_col() +
      facet_grid(. ~ source, labeller = labeller(source = facet_labels)) +
      theme_classic(base_family = "sans", base_size = 20) +
      labs(y = "Gene", x = "Frequency (%)", fill = NULL) +
      theme(axis.text.y = element_text(size = 14, face = "italic"),
            legend.key.size = grid::unit(.5, "lines"), legend.text = element_text(size = 12),
            strip.background = element_blank()) +
      scale_fill_manual(values = setNames(palette$color, palette$legend_label))
    save_plot(p, "prevalence", 12, 8)
  }
}

file <- file.path(args[1], "mutation_counts.csv")
if (file.exists(file)) {
  d <- read_csv(file, show_col_types = FALSE)
  if (nrow(d)) {
    sources <- unique(d$source)
    d$source <- factor(d$source, levels = unique(c(intersect(c("Pre-Treatment", "Pre-CAR-T"), sources), sources)))
    breaks <- sort(unique(c(0, 1, 2, 3, 4, 5, 6, 8, 10, 15, max(d$unique_mutations))))
    fills <- setNames(rep(c("lightblue", "lightgrey"), length.out = length(sources)), sources)
    if ("Pre-Treatment" %in% sources) fills["Pre-Treatment"] <- "lightblue"
    if ("Pre-CAR-T" %in% sources) fills["Pre-CAR-T"] <- "lightgrey"
    p <- ggplot(d, aes(source, unique_mutations, fill = source)) +
      geom_hline(yintercept = breaks, linetype = "dashed", color = "grey") +
      geom_violin(alpha = .5) +
      geom_point(position = position_jitter(width = .3, height = 0, seed = 1), size = 1.5) +
      labs(x = "", y = "Unique CH Mutations per Participant", fill = NULL) +
      theme_classic(base_size = 20, base_family = "sans") +
      scale_y_continuous(breaks = breaks, limits = c(0, max(breaks)), transform = pseudo_log_trans(base = 10)) +
      scale_fill_manual(values = fills)
    save_plot(p, "mutation_counts")
  }
}

files <- list.files(args[1], pattern = "^(cfDNA_vs_|PBMC_vs_).*csv$", full.names = TRUE)
for (file in files) {
  d <- read_csv(file, show_col_types = FALSE)
  name <- sub("\\.csv$", "", basename(file))
  # The original figure displays these four genes; all seven stay in the CSVs.
  d <- d %>% filter(gene %in% c("TP53", "DNMT3A", "PPM1D", "TET2"), is.finite(corr))
  if (!nrow(d)) next
  counts <- d %>% count(gene) %>% mutate(label = paste0(gene, "\n(n=", n, ")"))
  d$gene <- factor(d$gene, levels = unique(d$gene))
  p <- ggplot(d, aes(gene, corr)) +
    geom_boxplot(aes(fill = gene), outlier.shape = NA) +
    geom_point(aes(color = gene), position = position_jitter(width = .2, height = 0, seed = 1), size = 2) +
    scale_fill_manual(values = setNames(lighten(unname(cols), .4), names(cols))) +
    scale_color_manual(values = setNames(darken(unname(cols), .2), names(cols))) +
    scale_x_discrete(labels = setNames(counts$label, counts$gene)) +
    theme_classic(base_family = "sans", base_size = 20) + theme(legend.position = "none") +
    labs(x = "Gene", y = "Correlation", title = gsub("_", " ", name))
  if (name %in% c("cfDNA_vs_PDWB", "cfDNA_vs_PBMC")) p <- p + coord_cartesian(ylim = c(-1, 1))
  save_plot(p, name)
}
