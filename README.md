# Polyclonal selection of immune checkpoint mutations in thyroid autoimmunity

> **Paper:** Polyclonal selection of immune checkpoint mutations in thyroid autoimmunity

> **Authors:** Nicola, Lawson, et al.  

> **Year:** 2026

---

## Table of Contents

- [Paper analysis](#paper-analysis)
  - [Data](#data)
  - [Dependencies](#dependencies)
  - [How to run](#how-to-run)
- [PTA analysis](#pta-analysis)
  - [Data](#data-1)
  - [Dependencies](#dependencies-1)
  - [How to run](#how-to-run-1)
- [License](#license)

---

## Paper analysis

### Data

#### Internal data (already in repo)

Many of the inputs for `HashimotoAnalysis.Rmd` are included in this repository
in the zipped directory `data/rmd_input.zip`, which must be unzipped.

```bash
unzip data/rmd_input.zip -d data/
```

The directory should look like this when unzipped.

```
data/rmd_input/
├── 2025-07-18_TNFRSF14_Uniprot_Q92956_Domains.json
├── 2025-07-18_TNFRSF14_Uniprot_Q92956_PTM.json
├── 2025-07-22_CBL_Uniprot_P22681_Domains.json
├── 2025-07-22_CD274_Uniprot_Q9NZQ7_PTM.json
├── combined_B_mem_1_final.tsv
├── discarded_variants_high_vaf_drivers.tsv
├── hashimoto_emseq_cell_fractions_updated_atlas_epidish.tsv
├── hashimoto_exome_targeted_combined_muts.tsv
├── hashimoto_exome_targeted_combined_per_sample_cov.tsv
├── hashimoto_exome_targeted_combined_summary_stats.tsv
├── noshm_exome_5e-07_B_mem_v3.Rdat
├── RefCDS_GRCh37_v1_NSX_codon.Rdat
├── RefCDS_GRCh37_v1_NSX_noshm_codon.Rdat
├── RefCDS_GRCh37_v1.NSXupdate.Rdat
├── Sanger_Immune-v1_TE-91661256_hg19_gene_list.tsv
├── shm_exome_5e-07_B_mem_v3.Rdat
├── per_gene_per_site_cov/
├── antibody_synthesis/
├── interpro_domain/
├── pta/
├── spatial_mapping/
└── spatial_mapping_batch_2/

6 directories, 16 files
```

#### External data (must be downloaded)

The following dataset is not included in this repository and must be downloaded
before running the analysis.

- `data/reference/1kgp/GRCh37/hs37d5.fa` - https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/

The following datasets are not included in this repository and must be requested
and downloaded before running the analysis.

- `data/rmd_input/normal_lymphocytes/` - Data being released for manuscript in 
preparation. Please contact corresponding authors if access is required sooner.
- `data/rmd_input/other_normal_samples/` - Data being released for manuscript in 
preparation. Please contact corresponding authors if access is required sooner.
- `data/rmd_input/cosmic_data/` - Data must be downloaded from the COSMIC 
database at [www.cosmickb.org](https://www.cosmickb.org/). Only download data from samples that have had whole genome or exome sequencing. COSMIC data were re-annotated with 
`data/rmd_input/RefCDS_GRCh37_v1.NSXupdate.Rdat` using `dndscv` and split by
chromosome.

The `data/rmd_input/cosmic_data/` directory should contain the following files.

```bash
ls data/rmd_input/cosmic_data/
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr10.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr22.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr11.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr2.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr12.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr3.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr13.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr4.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr14.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr5.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr15.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr6.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr16.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr7.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr17.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr8.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr18.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chr9.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr19.tsv  COSMIC_v99_WholeGenomeExome_FullAnnotated_chrX.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr1.tsv   COSMIC_v99_WholeGenomeExome_FullAnnotated_chrY.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr20.tsv  COSMIC_v99_WholeGenomeExome_SamplePhenotypeInfo.tsv
COSMIC_v99_WholeGenomeExome_FullAnnotated_chr21.tsv

head -2 data/rmd_input/cosmic_data/COSMIC_v99_WholeGenomeExome_FullAnnotated_chr1.tsv
```
|mstr                  |sampleID    | chr|   pos|ref |mut |gene  | strand|ref_cod |mut_cod |ref3_cod |mut3_cod |aachange |ntchange |codonsub |impact     |pid             |
|:---------------------|:-----------|---:|-----:|:---|:---|:-----|------:|:-------|:-------|:--------|:--------|:--------|:--------|:--------|:----------|:---------------|
|COSS2120562:1:69224:C |COSS2120562 |   1| 69224|A   |C   |OR4F5 |      1|A       |C       |GAC      |GCC      |D45A     |A134C    |GAC>GCC  |Missense   |ENSP00000334393 |
|COSS2120551:1:69230:C |COSS2120551 |   1| 69230|A   |C   |OR4F5 |      1|A       |C       |CAC      |CCC      |H47P     |A140C    |CAC>CCC  |Missense   |ENSP00000334393 |
|COSS2120562:1:69236:C |COSS2120562 |   1| 69236|A   |C   |OR4F5 |      1|A       |C       |CAC      |CCC      |H49P     |A146C    |CAC>CCC  |Missense   |ENSP00000334393 |
|COSS2339604:1:69270:G |COSS2339604 |   1| 69270|A   |G   |OR4F5 |      1|A       |G       |CAC      |CGC      |S60S     |A180G    |TCA>TCG  |Synonymous |ENSP00000334393 |

---

### Dependencies

#### Languages

The scripts in this analysis use R (v4.5.0).

#### R packages

- `Rsamtools` (v2.26.0)
- `Biostrings` (v2.78.0)
- `XVector` (v0.50.0)
- `ggpubr` (v0.6.2)
- `ggh4x` (v0.3.1)
- `ggtree` (v4.0.4)
- `ape` (v5.8-1)
- `pander` (v0.6.6)
- `drc` (v3.0-1)
- `gtools` (v3.9.5)
- `stringi` (v1.8.7)
- `XML` (v3.99-0.20)
- `ggforce` (v0.5.0)
- `jsonlite` (v2.0.0)
- `MASS` (v7.3-65)
- `vcfR` (v1.15.0)
- `knitr` (v1.51)
- `latticeExtra` (v0.6-31)
- `lattice` (v0.22-6)
- `RColorBrewer` (v1.1-3)
- `viridis` (v0.6.5)
- `viridisLite` (v0.4.3)
- `patchwork` (v1.3.2)
- `scales` (v1.4.0)
- `dndscv` (v0.0.1.0)
- `lubridate` (v1.9.4)
- `forcats` (v1.0.1)
- `stringr` (v1.6.0)
- `dplyr` (v1.2.0)
- `purrr` (v1.2.1)
- `readr` (v2.1.6)
- `tidyr` (v1.3.2)
- `tibble` (v3.3.1)
- `ggplot2` (v4.0.2)
- `tidyverse` (v2.0.0)
- `GenomicRanges` (v1.62.1)
- `Seqinfo` (v1.0.0)
- `IRanges` (v2.44.0)
- `S4Vectors` (v0.48.0)
- `BiocGenerics` (v0.56.0)
- `generics` (v0.1.4)
- `readxl` (v1.4.5)

---

### How to run

Clone the repository.

# TODO: update to zenodo link

```bash
git clone https://github.com/alextidd/nicola_et_al_2026/
cd nicola_et_al_2026
```

Install all dependencies described in the [Dependencies](#dependencies) section
above. 

Download all data listed in the [Data](#data) section above.

#### Running the analysis

The `src/rmd/dNdS_shm_RefCDS_creation.R` script demonstrates how the RefCDS 
objects can be generated for running dNdSshm. These files are already available 
in `data/rmd_input/`.

All analyses are contained within the `src/rmd/HashimotoAnalysis.Rmd` 
Rmarkdown script. Render the report with the following command.

```r
rmarkdown::render('src/rmd/HashimotoAnalysis.Rmd')
```

This will save all outputs to `output/` and will regenerate the report at 
`src/rmd/HashimotoAnalysis.html`.

---

## PTA analysis

This section contains all analyses that were run on the PTA sequencing data
from donors with Hashimoto thyroiditis. This includes somatic variant calling 
with `BaseJumper`, somatic variant genotyping with `nf-resolveome`, phylogenetic 
analysis with `Sequoia`, and signature analysis with `HDP`, `sigfit` and 
`SigProfiler`. 

### Data

#### Internal data (already in repo)

The following datasets are already included in this repository.

```
data/
├── nanoseq
│   ├── hashimoto_exome_targeted_combined_muts.tsv
│   └── metadata.yaml
├── reference
│   └── gatk
│       └── GRCh38
│           └── genome.fa.dict
├── resolveome
│   ├── manual_inspection
│   │   ├── 20250902_pta_additional_annotation_H1.tsv
│   │   ├── H1_PD63118_pta_additional_annotation.tsv
│   │   ├── metadata.yaml
│   │   └── PD63118.tsv
│   └── shared_clades
│       ├── metadata.yaml
│       ├── shared_clades.tsv
│       └── shared_clades.txt
├── signatures
│   ├── lodato_2018
│   │   ├── Lodato2018_SignatureData_Aging.csv
│   │   └── metadata.yaml
│   ├── luquette_2022
│   │   ├── metadata.yaml
│   │   └── snv.artifact.signature.v3.rda
│   ├── machado_2022
│   │   ├── 41586_2022_5072_MOESM4_ESM.xlsx
│   │   ├── metadata.yaml
│   │   │   └── S8_finalsignaturetable.tsv
│   └── petljak_2019
│       ├── metadata.yaml
│       └── mmc1.tsv
├── twist
│   ├── metadata.yaml
│   ├── Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19_gene_info.csv
│   └── Sanger_Immune-v1_TE-91661256_hg19_reformatted_220.bed
└── vdj_coverage
    ├── ig_tcr_genes_pseudogenes.tsv
    └── metadata.yaml
```

#### External data (must be downloaded)

The following datasets are not included in this repository and must be 
downloaded separately before running any analysis. Please save them to the paths 
in the `file` column.

| file | link |
| --- | --- |
| data/reference/liftover/hg19ToHg38.over.chain.gz | http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz |
| data/reference/nanoseq/SNP_GRCh37.wgns.bed.gz | https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL |
| data/reference/cosmic/COSMIC_v3.5_SBS_GRCh38.txt | https://cancer.sanger.ac.uk/signatures/downloads/ |
| data/reference/gatk/GRCh38/Homo_sapiens_assembly38.fasta.gz | [ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/](ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/) |
| data/reference/1kgp/GRCh37/hs37d5.fa | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz |

##### BAMs and SNPs

The PTA BAMs and CaVEMan SNP calls must be downloaded from EGA. Please consult
the **Data availability** statement in the paper for information on how to 
download these.

Once the BAMs are downloaded, please update the `bam` column of the samplesheet 
at `data/resolveome/bams/samplesheet.csv` with full paths to their locations.

| id | donor_id | seq_type | bam |
| --- | --- | --- | --- |
| plate10_wellA10_dna_run50382 | PD63118 | dna | plate10_wellA10_dna_run50382.bam |
| plate10_wellA10_dnahyb_run50227 | PD63118 | dnahyb | plate10_wellA10_dnahyb_run50227.bam |
| plate10_wellA11_dna_run50382 | PD63118 | dna | plate10_wellA11_dna_run50382.bam |
| plate10_wellA11_dnahyb_run50227 | PD63118 | dnahyb | plate10_wellA11_dnahyb_run50227.bam |

Once the CaVEMan SNPs are downloaded, please save them to `data/caveman/`.

```
data/caveman/
├── PD66718b_lo0041.caveman_c.snps.vcf.gz
├── PD63118b_lo0044.caveman_c.snps.vcf.gz
├── PD63121d_lo0022.caveman_c.snps.vcf.gz
└── PD63126b_lo0010.caveman_c.snps.vcf.gz
```

---

### Dependencies

#### Languages

The scripts in this repository use R, Python and Nextflow.

Singularity v3.11.6 was used to run the BaseJumper Nextflow pipelines. The 
BaseJumper pipelines additionally require a Sentieon license in order to run 
DNAscope.

[MPBoot](https://github.com/diepthihoang/mpboot) must be downloaded in order to
run [Sequoia](https://github.com/TimCoorens/Sequoia) within the
`01_run_sequoia.R` script. Once the MPBoot binary has been downloaded, edit this
script to add the correct `--mpboot_path` parameter.

#### Python packages

- `SigProfilerExtractor`
- `SigProfilerMatrixGenerator`
- `SigProfilerAssignment`
- `pandas`
- `argparse`

#### R packages

- `magrittr` (v2.0.4)
- `tidyverse` (v2.0.0)
- `data.table` (v1.18.2.1)
- `ape` (v5.8)
- `patchwork` (v1.2.0)
- `RColorBrewer` (v1.1-3)
- `lsa` (v0.73.3)
- `slider` (v0.3.2)
- `ggh4x` (v0.2.8)
- `janitor` (v2.2.0)
- `knitr` (v1.51)
- `seqinr` (v4.2-36)
- `VGAM` (v1.1-12)
- `MASS` (v7.3-60.2)
- `devtools` (v2.4.6)
- `optparse` (v1.7.5)
- `hdp` (v0.1.5)
- `sigfit` (v2.2)
- `GenomicRanges` (v1.56.1)
- `rtracklayer` (v1.64.0)
- `biomaRt` (v2.60.1)
- `Rsamtools` (v2.20.0)
- `ggtree` (v3.12.0)
- `BiocManager` (v1.30.26)
- `treemut` (v1.1)

Some scripts also use helper functions from the R package
[alexr](https://github.com/alextidd/alexr). You can install this by running 
`devtools::install_github("alextidd/alexr@adea7780946218ac51211a61ab450a665c2f3cd1")`.

---

### How to run

Clone the repository.

```bash
git clone https://github.com/alextidd/nicola_et_al_2026/
cd nicola_et_al_2026
```

Install all dependencies described in the [Dependencies](#dependencies-1) section
above. 

Download all data listed in the [Data](#data-1) section above.

#### Nextflow pipelines

The analysis depends on the following Nextflow pipelines.

- [nf-core/bamtofastq](https://github.com/nf-core/bamtofastq) (commit hash: 8698321)
- [alextidd/nf-resolveome](https://github.com/alextidd/nf-resolveome) (commit hash: a785010)
- [alextidd/bj-dna-qc](https://github.com/alextidd/bj-dna-qc) (commit hash: 3149537)
- [alextidd/bj-somatic-variantcalling](https://github.com/alextidd/bj-somatic-variantcalling) (commit hash: ed6a84b)

Please download these into the `nextflow/` subdirectory.

```bash
$ tree -d nextflow/
nextflow/                   
├── bamtofastq
├── nf-resolveome
├── bj-dna-qc
└── bj-somatic-variantcalling
```

#### Running the analysis

Scripts are numbered and intended to be run in the following order.

```
src/resolveome/
│
│   # 1. DNA QC and somatic variant calling
├── basejumper                   
│   ├── 00_liftover_immune_panel_intervals.R
│   ├── 00_setup.R
│   ├── 01_bamtofastq_run.sh
│   ├── 02_bj-dna-qc_dna_run.sh
│   ├── 03_bj-somatic-variantcalling_dna_run.sh
│   └── 04_bj-somatic-variantcalling_dnahyb_run.sh
│
│   # 2. somatic variant and SNP genotyping
├── nf-resolveome                
│   ├── 00_get_vdj_regions.R
│   ├── 00_setup.R
│   ├── 01_dna_run.sh
│   ├── 02_dnahyb_run.sh
│   └── 03_phase_snps.Rmd
│
│   # 3. build phylogeny
├── sequoia                      
│   └── 01_run_sequoia.R
│
│   # 4. signature analysis
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

## License

This project is licensed under the [MIT License](LICENSE).