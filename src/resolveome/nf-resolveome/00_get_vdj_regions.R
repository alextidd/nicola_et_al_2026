#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
dir.create("out/vdj_coverage/regions/", showWarnings = FALSE, recursive = TRUE)

# load data
df <-
  readr::read_tsv("data/vdj_coverage/ig_tcr_genes_pseudogenes.tsv") %>%
  janitor::clean_names() %>%
  dplyr::select(chr = chromosome_scaffold_name, start = gene_start_bp,
                end = gene_end_bp, gene = gene_name, starts_with("gene_")) %>%
  dplyr::filter(chr %in% c("2", "7", "14", "22")) %>%
  dplyr::mutate(class = dplyr::case_when(grepl("^IG", gene) ~ "BCR",
                                         grepl("^TR", gene) ~ "TCR",
                                         TRUE ~ NA))

# write tsv of genes
df %>%
  readr::write_tsv("out/vdj_coverage/regions/ig_tcr_genes.tsv")

# write bed of genes
df %>%
  dplyr::select(chr, start, end, gene) %>%
  readr::write_tsv("out/vdj_coverage/regions/ig_tcr_genes.bed",
                   col_names = FALSE)

# summarise BCR/TCR regions per chr
df_regions <-
  df %>%
  dplyr::group_by(chr, class) %>%
  dplyr::summarise(start = min(start), end = max(end)) %>%
  dplyr::select(chr, start, end, class)

# write bed of regions
df_regions %>%
  dplyr::select(chr, start, end, class) %>%
  readr::write_tsv("out/vdj_coverage/regions/ig_tcr_regions.bed",
                   col_names = FALSE)
