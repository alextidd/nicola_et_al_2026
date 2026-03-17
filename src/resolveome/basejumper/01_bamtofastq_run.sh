
#!/bin/bash
# runsub src/resolveome/basejumper/01_bamtofastq_run.sh

# modules
module load singularity/3.11.4

# dirs
wd=$(pwd)

# run
(
  cd data/fastqs/
  nextflow run nf-core/bamtofastq \
    -profile singularity,sanger \
    --input samplesheet.csv \
    --outdir . \
    -w $wd/work/bamtofastq/ \
    -resume \
    -N at31@sanger.ac.uk
)