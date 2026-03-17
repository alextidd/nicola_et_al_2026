#!/bin/bash

# dirs
wd=$(pwd)

# run
(
  cd out/resolveome/nf-resolveome/dnahyb/
  nextflow run nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --seq_type dnahyb \
    --bait_set_hyb $wd/data/twist/Sanger_Immune-v1_TE-91661256_hg19_reformatted_220.bed \
    --bait_set_vdj $wd/out/vdj_coverage/regions/ig_tcr_genes.bed \
    --fasta data/reference/1kgp/GRCh37/hs37d5.fa \
    --location local \
    --baf_chrs 1,4,9 \
    --out_dir ./ \
    -w $wd/work/nf-resolveome/dnahyb/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)
