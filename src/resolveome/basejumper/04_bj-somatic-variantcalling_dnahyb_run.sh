#!/bin/bash
# donor_id=PD63118
# donor_id=PD66718
# donor_id=PD63121
# donor_id=PD63126
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q week -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J resolveome_basejumper_04_bj-somatic-variantcalling_dnahyb_run -o log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_basejumper_04_bj-somatic-variantcalling_dnahyb_run.out -e log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_basejumper_04_bj-somatic-variantcalling_dnahyb_run.err "bash src/resolveome/basejumper/04_bj-somatic-variantcalling_dnahyb_run.sh ${donor_id}"

# parameters
donor_id=$1

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
run_dir=$lustre_dir/out/resolveome/basejumper/bj-somatic-variantcalling/dnahyb/$donor_id/

# modules
module load singularityce-4.1.0/python-3.11.6

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/basejumper/bj-somatic-variantcalling/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd $run_dir
  nextflow run $NFS_TEAM/nextflow/external/basejumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir $donor_id \
    --timestamp run \
    --variant_workflow_type somatic_heuristic_filter \
    --dnascope_model_selection bioskryb129 \
    --skip_variant_annotation false \
    --skip_sigprofile false \
    -c $wd/config/bj-somatic-variantcalling_dnahyb.config \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $lustre_dir/work/basejumper/bj-somatic-variantcalling/dnahyb/$donor_id/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume \
    -with-tower \
    -N at31@sanger.ac.uk
)