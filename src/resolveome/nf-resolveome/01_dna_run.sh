#!/bin/bash

# dirs
wd=$(pwd)

# run
(
  cd out/resolveome/nf-resolveome/dna/
  nextflow run $NFS_TEAM/nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --seq_type dna \
    --bait_set_hyb $wd/out/twist/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19.bed \
    --bait_set_vdj $wd/out/vdj_coverage/regions/ig_tcr_genes.bed \
    --fasta data/reference/1kgp/GRCh37/hs37d5.fa \
    --location local \
    --baf_chrs 1,4,9 \
    --baf_genes TNFRSF14,CD274,TET2 \
    --out_dir ./ \
    -w $wd/work/nf-resolveome/dna/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)