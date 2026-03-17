
genome.fa

mpboot


# Polyclonal selection of immune checkpoint mutations in thyroid autoimmunity

> **Paper:** Polyclonal selection of immune checkpoint mutations in thyroid autoimmunity

> **Authors:** Nicola, Lawson, et al.  

> **Year:** 2026

This repository contains all analyses that were run on the PTA sequencing data
from donors with Hashimoto thyroiditis. This includes somatic variant calling 
with `BaseJumper`, somatic variant genotyping with `nf-resolveome`, phylogenetic 
analysis with `Sequoia`, and signature analysis with `HDP`, `sigfit` and 
`SigProfiler`. 

---

## Table of Contents

- [Requirements](#requirements)
- [Data](#data)
- [Repository structure](#repository-structure)
- [How to run](#how-to-run)
- [Citation](#citation)
- [License](#license)

---

## Requirements

### Languages

The scripts in this repository use R, Python and Nextflow.

Singularity v3.11.6 was used to run the BaseJumper Nextflow pipelines. The 
BaseJumper pipelines additionally require a Sentieon license in order to run 
DNAscope.

[MPBoot](https://github.com/diepthihoang/mpboot) must be downloaded in order to
run [Sequoia](https://github.com/TimCoorens/Sequoia) within the
`01_run_sequoia.R` script. Once the MPBoot binary has been downloaded, edit this
script to add the correct `--mpboot_path` parameter.

### Nextflow pipelines

- [bamtofastq](https://github.com/nf-core/bamtofastq)
- [nf-resolveome](https://github.com/alextidd/nf-resolveome)
- [bj-dna-qc](https://github.com/alextidd/bj-dna-qc)
- [bj-somatic-variantcalling](https://github.com/alextidd/bj-somatic-variantcalling)

### Python packages

- SigProfilerExtractor
- SigProfilerMatrixGenerator
- SigProfilerAssignment
- pandas
- argparse

### R packages

- magrittr
- tidyverse
- data.table
- ape
- patchwork
- RColorBrewer
- lsa
- slider
- ggh4x
- janitor
- knitr
- seqinr
- VGAM
- MASS
- devtools
- optparse
- hdp
- sigfit
- GenomicRanges
- rtracklayer
- biomaRt
- Rsamtools
- ggtree
- BiocManager
- treemut

Some scripts also use helper functions from the R package `alexr`. You can 
install this by running `devtools::install_github("alextidd/alexr")`.

---

## Data

### External data (must be downloaded)

The following datasets are not included in this repository and must be downloaded separately before running any analysis.

data/reference/liftover/hg19ToHg38.over.chain.gz
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

data/reference/nanoseq/SNP_GRCh37.wgns.bed.gz
https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL

data/reference/cosmic/COSMIC_v3.5_SBS_GRCh38.txt
https://cancer.sanger.ac.uk/signatures/downloads/

data/reference/COSMIC_v3.4_SBS_GRCh38.txt
https://cancer.sanger.ac.uk/signatures/downloads/

data/reference/COSMIC_v2_SBS_GRCh38.txt
https://cancer.sanger.ac.uk/signatures/downloads/

data/signatures/machado_2022/S8_finalsignaturetable.tsv
- upload to github

data/resolveome/shared_clades/shared_clades.tsv
- upload to github

data/signatures/petljak_2019/mmc1.tsv
- upload to github

data/signatures/lodato_2018/Lodato2018_SignatureData_Aging.csv
- upload to github

data/signatures/luquette_2022/snv.artifact.signature.v3.rda
- upload to github

data/resolveome/manual_inspection/20250902_pta_additional_annotation_H1.tsv
data/resolveome/manual_inspection/H1_PD63118_pta_additional_annotation.tsv
- upload to github
  

| Dataset | Source | Version / Build | Download link | Destination path |
|---------|--------|----------------|--------------|-----------------|
| Dataset A | Source name | vX.X | [Link](https://example.com) | `data/raw/dataset_a/` |
| Reference genome | Ensembl | GRCh38.XX | [Link](https://example.com) | `data/reference/` |
| Dataset B | Source name | — | [Link](https://example.com) | `data/raw/dataset_b/` |

After downloading, run the setup script to verify file integrity:

```bash
bash scripts/check_data.sh
```

### Processed / intermediate data

Processed data files that are too large to host on GitHub but are required to run downstream steps can be downloaded from [Zenodo (DOI: 10.5281/zenodo.XXXXXXX)](https://zenodo.org). Place them in `data/processed/` before proceeding.

---

## Repository structure

```
.
├── bin
│   ├── build_phylogeny.R
│   ├── run_Sigprofiler_Assignment.py
│   ├── run_Sigprofiler_Decompose.py
│   └── run_Sigprofiler_Extractor.py
├── config
│   ├── basejumper.config
│   ├── bj-somatic-variantcalling.config
│   └── bj-somatic-variantcalling_dnahyb.config
├── data
│   ├── nanoseq
│   │   ├── hashimoto_exome_targeted_combined_muts.tsv
│   │   └── metadata.yaml
│   └── resolveome
│       ├── manual_inspection
│       │   ├── metadata.yaml
│       │   └── PD63118.tsv
│       ├── twist
│       │   ├── metadata.yaml
│       │   ├── Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19_gene_info.csv
│       │   └── Sanger_Immune-v1_TE-91661256_hg19_reformatted_220.bed
│       └── vdj_coverage
│           ├── ig_tcr_genes_pseudogenes.tsv
│           └── metadata.yaml
├── README.md
└── src
    └── resolveome
        ├── basejumper
        │   ├── 00_liftover_immune_panel_intervals.R
        │   ├── 00_setup.R
        │   ├── 01_bamtofastq_run.sh
        │   ├── 02_bj-dna-qc_dna_run.sh
        │   ├── 03_bj-somatic-variantcalling_dna_run.sh
        │   └── 04_bj-somatic-variantcalling_dnahyb_run.sh
        ├── nf-resolveome
        │   ├── 00_get_vdj_regions.R
        │   ├── 00_setup.R
        │   ├── 01_dna_run.sh
        │   ├── 02_dnahyb_run.sh
        │   └── 03_phase_snps.Rmd
        ├── sequoia
        │   └── 01_run_sequoia.R
        └── signatures
            ├── 00_get_ref_signatures.R
            ├── 01_generate_matrices.R
            ├── 02a_run_hdp.R
            ├── 02b_run_sigfit.R
            ├── 03a_run_sigprofiler_extractor.sh
            ├── 03b_run_sigprofiler_decomposition.py
            └── 03c_run_sigprofiler_assignment.sh
```

---

## How to run

### 1. Clone the repository

```bash
git clone https://github.com/username/repo-name.git
cd repo-name
```

### 2. Install dependencies

```r
# In R
install.packages(c("tidyverse", "ggplot2", "patchwork"))
# Add any Bioconductor packages:
BiocManager::install(c("GenomicRanges", "deepSNV"))
```

### 3. Download data

Follow the instructions in the [Data](#data) section above, then verify:

```bash
bash scripts/check_data.sh
```

### 4. Run the analysis

Scripts are numbered and intended to be run in order:

```bash
Rscript scripts/01_preprocess.R
Rscript scripts/02_analysis.R
Rscript scripts/03_figures.R
```

Alternatively, run the full pipeline end-to-end:

```bash
bash run_all.sh
```

Output figures will be written to `results/figures/` and tables to `results/tables/`.

> **Note:** Approximate runtime and memory requirements on a standard workstation: ~X hours, ~X GB RAM.

---

## License

This project is licensed under the [MIT License](LICENSE).