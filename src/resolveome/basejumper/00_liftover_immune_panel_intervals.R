#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)

# dirs
dir.create("out/twist")

# load interval list
intervals <-
  readr::read_csv("data/twist/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19_gene_info.csv")

# get gene strandedness
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
strandedness <-
  getBM(
    attributes = c("hgnc_symbol", "strand"),
    filters = "hgnc_symbol",
    values = intervals$gene,
    mart = mart) %>%
  dplyr::transmute(gene = hgnc_symbol, strand) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(strand = ifelse(sum(strand) < 0, "-", "+"))

# convert to grch38, add strandedness
hg38_intervals <-
  intervals %>%
  {GRanges(
    seqnames = paste0("chr", .$chr),
    ranges = IRanges(start = .$bed_start, end = .$bed_end),
    gene = .$gene, target_type = .$target_type)} %>%
  liftOver(import.chain("../../reference/liftover/hg19ToHg38.over.chain")) %>%
  unlist() %>%
  tibble::as_tibble() %>%
  dplyr::select(seqnames, start, end, gene) %>%
  dplyr::left_join(strandedness, by = "gene") %>%
  # make missing genes + stranded
  dplyr::mutate(strand = ifelse(is.na(strand), "+", strand)) %>%
  dplyr::transmute(chr = seqnames, start, end, strand, gene)

# convert to picard-style .interval_list
# Picard-style interval files have a SAM-like header that includes a sequence
# dictionary. The intervals are given in the form
# <chr>:<start>-<stop> + <target_name>, with fields separated by tabs, and the
# coordinates are 1-based (first position in the genome is position 1, not
# position 0).
dict_header <-
  "/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa.dict" %>%
  readr::read_lines()

# save to picard interval_list format
interval_list <- "out/twist/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg38.interval_list"
readr::write_lines(dict_header, interval_list)
readr::write_tsv(hg38_intervals, interval_list,
                 append = TRUE, col_names = FALSE)

# save to bed format
hg38_intervals %>%
  dplyr::select(chr, start, end) %>%
  readr::write_tsv("out/twist/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg38.bed",
                   col_names = FALSE)
