#!/bin/bash

# dirs
wd=$(pwd)

# run
(
  cd out/resolveome/nf-resolveome/dnahyb/
  nextflow run $NFS_TEAM/nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --seq_type dnahyb \
    --bait_set_hyb $wd/data/twist/Sanger_Immune-v1_TE-91661256_hg19_reformatted_220.bed \
    --bait_set_vdj $wd/out/vdj_coverage/regions/ig_tcr_genes.bed \
    --fasta /lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh37d5/genome.fa \
    --location local \
    --baf_chrs 1,4,9 \
    --out_dir ./ \
    -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/nf-resolveome/dnahyb/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)
