
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

List software, languages, and key packages needed.

R, Python, Nextflow
mpboot

- **R** >= 4.x.x
  - `tidyverse`, `ggplot2`, `patchwork`, `<other packages>`
- **Python** >= 3.x (if applicable)
  - `pandas`, `numpy`, `<other packages>`
- **Other tools:** `samtools >= 1.x`, `bedtools >= 2.x`, etc.

---

## Data

### External data (must be downloaded)

hg19ToHg38.over.chain.gz
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

SNP_GRCh37.wgns.bed.gz
https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL

The following datasets are not included in this repository and must be downloaded separately before running any analysis.

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