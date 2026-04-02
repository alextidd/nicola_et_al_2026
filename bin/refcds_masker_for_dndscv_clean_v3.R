#' subsetrefcds
#'
#' Function to remove regions from a RefCDS object. These regions may be anything from a set of SNP sites, to a set of regions under a specific mutational process such as somatic hypermutation. 
#'
#' This function relies on bedtools - it needs to be installed and on a findable path
#'
#' @author fa8 (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#'
#' @param regions_to_be_masked GenomicRanges object listing the regions to be masked.
#' @param genomefile Path to the indexed reference genome file.
#' @param refcds Path to the RefCDS to be modified.
#' @param outfile Output file name .
#' @param invert, either TRUE or FALSE. If TRUE, for genes overlapping regions_to_be_masked the inverse will be masked. The rest of genes will be left unmodified.
#' @param subtractBed_path, path and command to bedtool's subtractBed (e.g. /Users/fa8/homebrew/bin/subtractBed)
#'
#' @export

subsetrefcds = function(regions_to_be_masked, genomefile, refcds, outfile, invert = FALSE, subtractBed_path="/Users/fa8/homebrew/bin/subtractBed") {

	library(GenomicRanges)
	library(Rsamtools    )
	library(dndscv       )

    message("[1/7] ",length(regions_to_be_masked)," regions to be masked, spanning ",
            sum(width(regions_to_be_masked))," bps")

	load(refcds) # this loads two objects RefCDS and gr_genes
	
	# Find which genes are affected (needed if invert is TRUE):
	ovs       = findOverlaps(regions_to_be_masked,gr_genes)
	ovs_df    = as.data.frame(ovs)
	df_genes  = as.data.frame(gr_genes)
	genes_with_overlaps = unique(df_genes[ovs_df$subjectHits,"names"])  # I will need those to prepare the inverse
	regions_to_be_masked.bck = regions_to_be_masked
	regions_to_be_masked = regions_to_be_masked[unique(ovs_df$queryHits)]
	message("[2/7] ",length(genes_with_overlaps), 
	        " genes overlapping the regions to be masked ")
	message("[3/7] ",length(regions_to_be_masked), 
	        " regions to be masked overlapping genes, spanning ",
            sum(width(regions_to_be_masked))," bps")
	
	bedTools.2in<-function(functionstring="bedIntersect",bed1,bed2,opt.string="")
	{
	  #create temp files
	  a.file=tempfile()
	  b.file=tempfile()
	  out   =tempfile()
	  options(scipen =99) # not to use scientific notation when writing out
	 
	  #write bed formatted dataframes to tempfile
	  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
	  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
	 
	  # create the command string and call the command using system()
	  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
	  cat(command,"\n")
	  try(system(command))
	 
	  res=read.table(out,header=F)
	  unlink(a.file);unlink(b.file);unlink(out)
	  options(scipen = 0) 
	  return(res)
	}

	if(invert == TRUE) {
		# Need to invert regions, but only for the genes with overlaps	
		# This is the only missing bit	
		new_gr_genes = as.data.frame(gr_genes)
		new_gr_genes = new_gr_genes[which(new_gr_genes$names %in% genes_with_overlaps),]
		# regions_to_filter = as.data.frame(setdiff(new_gr_genes,ranges_gr))
		# Bedtools:
			df_genes           = new_gr_genes
			df_genes           = df_genes[,c("seqnames","start","end","names")]
			colnames(df_genes) = c("chr","start","end","gene")
			df_genes$start     = df_genes$start-1 # to bed format
			
			df_regions_to_be_masked           = as.data.frame(regions_to_be_masked)
			df_regions_to_be_masked           = df_regions_to_be_masked[,c("seqnames","start","end")]
			colnames(df_regions_to_be_masked) = c("chr","start","end")
			df_regions_to_be_masked$start     = df_regions_to_be_masked$start-1 # to bed format
			
			resultado                         = bedTools.2in(functionstring=subtractBed_path,df_genes,df_regions_to_be_masked)
			colnames(resultado)               = c("chr","start","end","gene")
			resultado$start                   = resultado$start+1
			df_genes$start                    = df_genes$start+1
			df_regions_to_be_masked$start     = df_regions_to_be_masked$start+1
			
			regions_to_be_masked = GRanges(resultado$chr, IRanges(start=resultado$start,end=resultado$end))
			GenomicRanges::mcols(regions_to_be_masked)$names = resultado$gene


		## done
	}
	
	ranges_df = as.data.frame(regions_to_be_masked)
	chrs = unlist(sapply(1:nrow(ranges_df),function(x) rep(ranges_df[x,"seqnames"],ranges_df[x,"end"]-ranges_df[x,"start"]+1)))
	poss = unlist(sapply(1:nrow(ranges_df),function(x) ranges_df[x,"start"]:ranges_df[x,"end"]))
	to_filter = data.frame(chr=chrs,pos=poss)
	# For each mask site, get the base at that site:
	to_filter$ref = as.vector(scanFa(genomefile,GRanges(to_filter$chr, IRanges(start=to_filter$pos, end=to_filter$pos))))
	
	# Expand to all possible subs (and then remove cases with N in ref or with ref = mut)
	muts = rep(c("A","C","G","T"),nrow(to_filter))
	mut_matrix = data.frame(sampleID="none",chr=rep(to_filter$chr,each=4),pos=rep(to_filter$pos,each=4),ref=rep(to_filter$ref,each=4),mut=muts)
	mut_matrix = mut_matrix[which(mut_matrix$ref %in% c("A","C","G","T")),]
	mut_matrix = mut_matrix[which(mut_matrix$ref != mut_matrix$mut),]
	message("[4/7] Ready to mask ",nrow(mut_matrix), " subs for ",nrow(to_filter)," sites...")
	# Then run dndscv with that mutation matrix
	dndsout = dndscv(mut_matrix,refdb=refcds,cv=NULL,max_coding_muts_per_sample=Inf,max_muts_per_gene_per_sample=Inf,outmats=T,outp=1)
	
	# Modify the RefCDS and save it
	load(refcds) # this loads two objects RefCDS and gr_genes
	for (j in 1:length(RefCDS)) {
	  RefCDS[[j]]$L = RefCDS[[j]]$L - dndsout$N[,,j]
	}

	# Modify gr_genes:
	# Change this to avoid "reduce" in setdiff, doing one gene at a time
	# Now relying on bedtools
	
	gr_genes.bck       = gr_genes
	
	df_genes           = as.data.frame(gr_genes)
	df_genes           = df_genes[,c("seqnames","start","end","names")]
	colnames(df_genes) = c("chr","start","end","gene")
	df_genes$start     = df_genes$start-1 # to bed format
	
	df_regions_to_be_masked           = as.data.frame(regions_to_be_masked)
	df_regions_to_be_masked           = df_regions_to_be_masked[,c("seqnames","start","end")]
	colnames(df_regions_to_be_masked) = c("chr","start","end")
	df_regions_to_be_masked$start     = df_regions_to_be_masked$start-1 # to bed format
	
	resultado                         = bedTools.2in(functionstring=subtractBed_path,df_genes,df_regions_to_be_masked)
	colnames(resultado)               = c("chr","start","end","gene")
	resultado$start                   = resultado$start+1
	df_genes$start                    = df_genes$start+1
	df_regions_to_be_masked$start     = df_regions_to_be_masked$start+1
	
	gr_genes = GRanges(resultado$chr, IRanges(start=resultado$start,end=resultado$end))
	GenomicRanges::mcols(gr_genes)$names = resultado$gene
	
	if(invert == TRUE) {
		# Filter gr_genes
		new_gr_genes = as.data.frame(gr_genes)
		new_gr_genes = new_gr_genes[which(new_gr_genes$names %in% genes_with_overlaps),]
		gr_genes = GRanges(new_gr_genes$seqnames, IRanges(start=new_gr_genes$start,end=new_gr_genes$end))
		GenomicRanges::mcols(gr_genes)$names = new_gr_genes$names
		# Filter RefCDS:
		genes_in_refcds = sapply(1:length(RefCDS),function(x) RefCDS[[x]]$gene_name)
		keep_these = which(genes_in_refcds %in% genes_with_overlaps)
		newRefCDS = list()
		for(idx in keep_these) {
			newRefCDS[[length(newRefCDS)+1]] = RefCDS[[idx]]
		}
		RefCDS = newRefCDS
	}
	
	##### df_genes = as.data.frame(gr_genes)
	##### setdiffs = list()
	##### for(gene in genes_with_overlaps) {
	##### 	gr_gene = GRanges(df_genes[which(df_genes$names==gene),"seqnames"], 
	##### 	                  IRanges(start=df_genes[which(df_genes$names==gene),"start"], 
	##### 	                          end = df_genes[which(df_genes$names==gene),"end"]))
	##### 	setdiffs[[gene]] = as.data.frame(setdiff(gr_gene,regions_to_be_masked))
	##### 	setdiffs[[gene]][,"names"] = gene
	##### }
	##### untouched_genes = setdiff(unique(df_genes$names),genes_with_overlaps)
	##### for(gene in untouched_genes) {
	##### 	setdiffs[[gene]] = df_genes[which(df_genes[which(df_genes$names == gene),]),]
	##### }
	# new_genes_gr = setdiff(gr_genes,regions_to_be_masked)	

	# Add gene names back to new_genes_gr:
	# ovs = as.data.frame(findOverlaps(new_genes_gr,gr_genes))
	# new_genes_df = as.data.frame(new_genes_gr)
	# refcds_df    = as.data.frame(gr_genes)
	# new_genes_df[ovs$queryHits,"gene_name"] = refcds_df[ovs$subjectHits,"names"]
	# GenomicRanges::mcols(new_genes_gr)$names = new_genes_df$gene_name
	# gr_genes = new_genes_gr

	# Save new refcds objects:
	save(RefCDS,gr_genes,file=outfile)
	
}

