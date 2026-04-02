library(dndscv)
library(GenomicRanges)

source("/Users/al28/Projects/SortedBlood/refcds_masker_for_dndscv_clean_v3.R")

many_burdens = read.table("/Users/al28/Projects/Thyroid_final/input/combined_B_mem_1_final.tsv",sep="\t",header=T)

rows_to_remove = which(many_burdens$top_lci_burden > 5*10^-7)
exons_to_remove =  many_burdens$region_name[rows_to_remove]
shm_genes <- unique(gsub("_.*","",exons_to_remove))

# create GRanges and add 3 bp to cover splice sites:
exons_gr = GRanges(many_burdens[rows_to_remove,"chr"], IRanges(start=many_burdens[rows_to_remove,"start_pos"]-10,end=many_burdens[rows_to_remove,"end_pos"]+10))

subsetrefcds(regions_to_be_masked=exons_gr,
             genomefile="/Users/al28/Projects/Thyroid_final/input/hs37d5.fa",
             refcds="/Users/al28/Projects/Thyroid_final/input/RefCDS_GRCh37_v1.NSXupdate.Rdat",
             outfile="noshm_exome_5e-07_B_mem_v3.Rdat",
             invert=FALSE,
             subtractBed_path="/Users/al28/.homebrew/bin/subtractBed")

subsetrefcds(regions_to_be_masked=exons_gr,
             genomefile="/Users/al28/Projects/Thyroid_final/input/hs37d5.fa",
             refcds="/Users/al28/Projects/Thyroid_final/input/RefCDS_GRCh37_v1.NSXupdate.Rdat",
             outfile="shm_exome_5e-07_B_mem_v3.Rdat",
             invert=TRUE,
             subtractBed_path="/Users/al28/.homebrew/bin/subtractBed")
