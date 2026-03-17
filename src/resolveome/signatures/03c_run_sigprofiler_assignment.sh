#!/bin/bash

# modules
module load sigprofiler/extractor-1.1.24-GRCh38-GRCh37

# all signatures (cosmic + custom)
python3 bin/run_Sigprofiler_Assignment.py \
  --samples out/resolveome/signatures/matrices/trinuc_mut_mat_sigpro.txt \
  --output_dir out/resolveome/signatures/sigprofiler/assignment/all/ \
  --reference_genome GRCh38 \
  --signatures_file out/resolveome/signatures/ref_sigs/ref_sigs.tsv

# all signatures (cosmic only)
python3 bin/run_Sigprofiler_Assignment.py \
  --samples out/resolveome/signatures/matrices/trinuc_mut_mat_sigpro.txt \
  --output_dir out/resolveome/signatures/sigprofiler/assignment/cosmic/ \
  --reference_genome GRCh38 \
  --signatures_file ../../reference/cosmic/COSMIC_v3.4_SBS_GRCh38.txt

# consensus signatures (cosmic only)
python3 bin/run_Sigprofiler_Assignment.py \
  --samples out/resolveome/signatures/matrices/trinuc_mut_mat_sigpro.txt \
  --output_dir out/resolveome/signatures/sigprofiler/assignment/consensus_cosmic/ \
  --reference_genome GRCh38 \
  --include_sigs SBS1,SBS5,SBS9,SBS17a,SBS17b

# consensus signatures (cosmic + custom)
python3 bin/run_Sigprofiler_Assignment.py \
  --samples out/resolveome/signatures/matrices/trinuc_mut_mat_sigpro.txt \
  --output_dir out/resolveome/signatures/sigprofiler/assignment/consensus/ \
  --reference_genome GRCh38 \
  --signatures_file out/resolveome/signatures/ref_sigs/ref_sigs.tsv \
  --include_sigs SBS1,SBS5,SBS9,SBS17a,SBS17b,machado_2022_SBSblood,luquette_2022_PTA_artefact