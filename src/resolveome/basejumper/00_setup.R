#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(purrr)
library(dplyr)

# dirs
data_dir <- file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis/data/")
fastq_dir <- file.path(data_dir, "/fastqs/reads")
out_dir <- file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis/out/resolveome/basejumper")
dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)

# read samplesheet
ss <- readr::read_csv(file.path(data_dir, "bams/samplesheet_local.csv"))

# write samplesheet for fastqs
ss %>%
  dplyr::transmute(sample_id = id, mapped = bam, index = paste0(bam, ".bai"),
                   file_type = "bam") %>%
  # check if bam exists
  dplyr::filter(file.exists(mapped)) %>%
  readr::write_csv(file.path(data_dir, "fastqs/samplesheet.csv"))

# create dnahyb/dna bj-somatic-variantcalling samplesheets
# columns: biosampleName,read1,read2,groups,isbulk,bam
ss_fastq <-
  ss %>%
  transmute(
    donor_id, biosampleName = id,
    read1 = file.path(fastq_dir, paste0(id, "_1.merged.fastq.gz")),
    read2 = file.path(fastq_dir, paste0(id, "_2.merged.fastq.gz")),
    groups = donor_id, isbulk = FALSE, bam = "", seq_type) %>%
  {split(., .$seq_type)} %>%
  map(function(df) {
    split(df %>% select(-donor_id, -seq_type), df$donor_id)
  })

# set up all runs
bj_runs <-
  list(list(pipeline = "bj-dna-qc", seq_type = "dna"),
       list(pipeline = "bj-expression", seq_type = "rna"),
       list(pipeline = "bj-somatic-variantcalling", seq_type = "dna"),
       list(pipeline = "bj-somatic-variantcalling", seq_type = "dnahyb"))

# save samplesheets
bj_runs %>%
  walk(function(run) {
    lis <- ss_fastq[[run$seq_type]]
    walk2(names(lis), lis, function(i, df) {
      out_dir_i <- file.path(out_dir, run$pipeline, run$seq_type, i)
      dir.create(out_dir_i, recursive = TRUE, showWarnings = FALSE)
      df %>%
        select(biosampleName, read1, read2) %>%
        readr::write_csv(file.path(out_dir_i, "samplesheet.csv"))
    })
  })
