
#!/bin/bash

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