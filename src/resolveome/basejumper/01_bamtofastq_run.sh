
#!/bin/bash
# runsub src/resolveome/basejumper/01_bamtofastq_run.sh

# modules
module load singularity/3.11.4

# run
(
  cd $LUSTRE_125/projects/hashimoto_thyroiditis/data/fastqs/
  nextflow run nf-core/bamtofastq \
    -profile singularity,sanger \
    --input samplesheet.csv \
    --outdir . \
    -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/bamtofastq/ \
    -resume \
    -N at31@sanger.ac.uk
)