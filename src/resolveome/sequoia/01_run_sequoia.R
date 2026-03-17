#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(dplyr)
library(tibble)

# dirs
bj_dir <- "out/resolveome/basejumper/bj-somatic-variantcalling/dna/PD63118/"
seq_dir <- "out/resolveome/sequoia/"
dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

# get ids in tables
ids <-
  file.path(bj_dir, "/PD63118_run/CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS/Mat_NV_null.tsv") %>%
  read.table() %>%
  colnames() %>%
  {tibble(id = ., cell_id = gsub("_dna.*", "", .))}

# load manual inspection
man_insp <-
  readr::read_tsv("data/resolveome/manual_inspection/PD63118.tsv") %>%
  full_join(ids) %>%
  filter(!is.na(id))

# get dropout / doublet cells
bad_cell_ids <-
  man_insp %>%
  filter((suspected_doublet == TRUE | chr_dropout == TRUE)) %>%
  pull(id)
good_cells <-
  man_insp %>%
  dplyr::filter(!(id %in% bad_cell_ids))

# get cnv cells
loh_cells <-
  good_cells %>%
  filter((loh_1p == 1 | loh_9p == 1)) %>%
  left_join(
    file.path(bj_dir, "samplesheet.csv") %>%
      readr::read_csv() %>%
      transmute(id = biosampleName, cell_id = gsub("_dna.*", "", biosampleName))
  ) %>%
  pull(id)

# load and lift over caveman snps
snps <-
  "out/resolveome/nf-resolveome/muts_and_snps/PD63118/caveman_snps.tsv" %>%
  readr::read_tsv() %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  alexr::lift_over(
    lift_over_chain = "data/reference/liftover/hg19ToHg38.over.chain") %>%
  dplyr::mutate(mut_id = paste(chr, pos, ref, alt, sep = "_"))

# load and lift over common snps
common_snps <-
  "data/reference/nanoseq/SNP_GRCh37.wgns.bed.gz" %>%
  gzfile() %>%
  readr::read_tsv(col_names = c("#CHROM", "START", "POS")) %>%
  dplyr::transmute(chr = paste0("chr", `#CHROM`), pos = POS) %>%
  alexr::lift_over(
    lift_over_chain = "data/reference/liftover/hg19ToHg38.over.chain") %>%
  dplyr::mutate(mut_id = paste(chr, pos, sep = "_"))

# mask sites
c("NV", "NR") %>%
  purrr::set_names() %>%
  purrr::walk(function(i) {

    # load matrix
    mat <-
      paste0(bj_dir, "/PD63118_run/CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS/Mat_", i,
             "_null.tsv") %>%
      data.table::fread() %>%
      as.matrix(rownames = 1)

    # remove caveman snps
    filt_mat <- mat[!(rownames(mat) %in% snps$mut_id), ]

    # remove common snps
    mut_ids <- sub("^((?:[^_]+_){1}[^_]+)_.*", "\\1", rownames(filt_mat))
    filt_mat <- filt_mat[!(mut_ids %in% common_snps$mut_id), ]

    # save
    filt_mat %>%
      write.table(file = paste0(seq_dir, "/Mat_", i, "_snpmasked.tsv"),
                  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  })

# run sequoia
system(paste0("Rscript bin/build_phylogeny.R",
              " --input_nv ", seq_dir, "/Mat_NV_snpmasked.tsv",
              " --input_nr ", seq_dir, "/Mat_NR_snpmasked.tsv",
              " --output_dir ", seq_dir,
              " --snv_rho 0.4",
              " --indel_rho 0.4",
              " --germline_cutoff -10",
              " --min_cov 10",
              " --max_cov 500",
              " --cnv_samples ", paste(loh_cells, collapse = ","),
              " --exclude_samples ", paste(bad_cell_ids, collapse = ","),
              " --mpboot_path /nfs/casm/team268im/at31/bin/"))