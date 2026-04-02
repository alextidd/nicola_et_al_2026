# Wrapper function to run the split-dnds model to handle somatic hypermutation in B cells
# WARNING: Still to implement:
#   1. Combining p-values using psubs when pglobal not available due to lack of indels.
#   2. Calculating a merged (average) dN/dS for split genes.
#   3. Possible changes to dndscv.R: kc="qsubs"

dnds_shm = function(mutations, refdb_noshm, refdb_shm, dc_noshm = NULL, gene_list = NULL, excl_shm = c("IGLL5","MYO1E","ZNF595"), excl_noshm = c("IGLL5"), sm = "192r_3w", kc = "qsubs", cv = "hg19", max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, numcode = 1, mingenecovs = 500, onesided = F) {
  
  mutations = unique(mutations[,1:5]) # Unique nutations per sampleID
  
  ## 1. Gene lists from RefCDS objects
  load(refdb_noshm)
  Lpergene = vector()
  for(i in c(1:length(RefCDS))) { 
    Lpergene[RefCDS[[i]]$gene_name] = sum(RefCDS[[i]]$L);
  }
  excl_noshm <- c(excl_noshm,names(Lpergene[which(Lpergene==0)]))
  gene_list_noshm = sapply(RefCDS, function(x) x$gene_name)
  gene_list_noshm = setdiff(gene_list_noshm, excl_noshm)
  
  load(refdb_shm)
  gene_list_shm = sapply(RefCDS, function(x) x$gene_name)
  gene_list_shm = setdiff(gene_list_shm, excl_shm)
  
  if (!is.null(gene_list)) {
    gene_list_noshm = intersect(gene_list_noshm, gene_list)
    gene_list_shm = intersect(gene_list_shm, gene_list)
    gene_list = unique(c(gene_list_noshm,gene_list_shm)) # Restricting gene_list to genes in the input RefCDSs and excluding excl_shm and excl_noshm
  }
  
  ## 2. Known cancer genes to exclude from the indel model
  if (kc[1] == "cgc81") {
    
    data(list=sprintf("cancergenes_%s",kc), package="dndscv") # dNdScv default
    
  } else if (kc[1] == "qsubs") { # Running dndscv to define kc based on substitution drivers
    
    dndsout = dndscv(mutations = mutations, gene_list = gene_list_noshm, max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, onesided = onesided, dc = dc_noshm, refdb = refdb_noshm, sm = sm, cv = cv, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, numcode = numcode, mingenecovs = mingenecovs)
    known_cancergenes = as.vector(dndsout$sel_cv$gene_name[which(dndsout$sel_cv$qallsubs_cv<0.05)])
    
  } else {
    known_cancergenes = kc
  }
  
  ## 3. Running a separate dNdS model for non-SHM regions and SHM regions
  dndsout = dndscv(mutations = mutations, gene_list = gene_list_noshm, max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, onesided = onesided, dc = dc_noshm, refdb = refdb_noshm, kc = known_cancergenes, sm = sm, cv = cv, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, numcode = numcode, mingenecovs = mingenecovs, outmats = T)
  dndsout_shm = dndscv(mutations = mutations, gene_list = gene_list_shm, max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, onesided = onesided, refdb = refdb_shm, kc = known_cancergenes, sm = sm, cv = cv, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, numcode = numcode, mingenecovs = mingenecovs)
  
  ## 4. Merging the outputs
  dndsout$sel_loc_shm = dndsout_shm$sel_loc
  dndsout$genemuts_shm = dndsout_shm$genemuts
  dndsout$sel_merged = merge(dndsout$sel_cv, dndsout$sel_loc_shm[,-(2:5)], by = "gene_name", all = T) # Merging the two outputs
  
  # Combining number of mutations
  if("n_ind" %in% colnames(dndsout_shm$sel_cv)){
    nmuts = as.matrix(dndsout_shm$sel_cv[,c("n_syn","n_mis","n_non","n_spl","n_ind")])
  }else{
    nmuts = as.matrix(dndsout_shm$sel_cv[,c("n_syn","n_mis","n_non","n_spl")])
  }
  
  rownames(nmuts) = dndsout_shm$sel_cv$gene_name
  natozero = function(x) { x[is.na(x)] = 0; return(x) } # Function to convert NAs to 0s from a vector
  dndsout$sel_merged$n_syn = natozero(dndsout$sel_merged$n_syn) + natozero(nmuts[,"n_syn"][dndsout$sel_merged$gene_name])
  dndsout$sel_merged$n_mis = natozero(dndsout$sel_merged$n_mis) + natozero(nmuts[,"n_mis"][dndsout$sel_merged$gene_name])
  dndsout$sel_merged$n_non = natozero(dndsout$sel_merged$n_non) + natozero(nmuts[,"n_non"][dndsout$sel_merged$gene_name])
  dndsout$sel_merged$n_spl = natozero(dndsout$sel_merged$n_spl) + natozero(nmuts[,"n_spl"][dndsout$sel_merged$gene_name])
  if("n_ind" %in% colnames(dndsout$sel_merged) & "n_ind" %in% colnames(nmuts)){
    dndsout$sel_merged$n_ind = natozero(dndsout$sel_merged$n_ind) + natozero(nmuts[,"n_ind"][dndsout$sel_merged$gene_name])
  }else if("n_ind" %in% colnames(nmuts)){
    dndsout$sel_merged$n_ind = natozero(nmuts[,"n_ind"][dndsout$sel_merged$gene_name])
  }
  
  # Combining p-values
  combpval = function(dndsout, p1, p2, pcomb) {
    dndsout$sel_merged[,pcomb] = 1 - pchisq(-2 * (log(dndsout$sel_merged[,p1]) + log(dndsout$sel_merged[,p2])), df = 4)
    nas = which(is.na(dndsout$sel_merged[,p1]))
    dndsout$sel_merged[nas,pcomb] = dndsout$sel_merged[nas,p2]
    nas = which(is.na(dndsout$sel_merged[,p2]))
    dndsout$sel_merged[nas,pcomb] = dndsout$sel_merged[nas,p1]
    qname = sprintf("q%s",substr(pcomb,2,nchar(pcomb)))
    dndsout$sel_merged[,qname] = p.adjust(dndsout$sel_merged[,pcomb], method = "BH")
    return(dndsout)
  }
  
  if (onesided == T) {
    if (!is.null(dndsout$sel_merged$pglobalpos_cv)) { # With indels
      dndsout = combpval(dndsout, p1="pglobalpos_cv", p2="ppos_loc", pcomb="pglobalpos_m")
      dndsout$sel_merged = dndsout$sel_merged[order(dndsout$sel_merged[,"pglobalpos_m"]), ]
    } else { # Without indels
      dndsout = combpval(dndsout, p1="ppos_cv", p2="ppos_loc", pcomb="psubpos_m")
      dndsout$sel_merged = dndsout$sel_merged[order(dndsout$sel_merged[,"psubpos_m"]), ]
    }
  } else {
    if (!is.null(dndsout$sel_merged$pglobal_cv)) { # With indels
      dndsout = combpval(dndsout, p1="pglobal_cv", p2="pall_loc", pcomb="pglobal_m")
      dndsout$sel_merged = dndsout$sel_merged[order(dndsout$sel_merged[,"pglobal_m"]), ]
    } else { # Without indels
      dndsout = combpval(dndsout, p1="pallsubs_cv", p2="pall_loc", pcomb="pallsubs_m")
      dndsout$sel_merged = dndsout$sel_merged[order(dndsout$sel_merged[,"pallsubs_m"]), ]
    }
  }
  
  return(dndsout)
}

# geneciloc
# Function to calculate confidence intervals for dN/dS values per gene under the dNdSloc model using a Poisson ratio.

geneciloc = function(genemuts, gene_list = NULL, level = 0.95) {
  
  # Ensuring valid level value
  if (level > 1) {
    warning("Confidence level must be lower than 1, using 0.95 as default")
    level = 0.95
  }
  
  # Subfunction to run the poisson.test function on each gene (r = expected neutral ratio given the substitution model and the gene sequence)
  poissonci = function(ns, s, r) {
    return(t(rbind((ns/s)/r, sapply(1:length(ns), function(j) poisson.test(x=c(ns[j],s[j]), conf.level = level)$conf.int / r[j]))))
  }
  
  # Calculating CI95% across all genes
  if (!is.null(gene_list)) {
    genemuts = genemuts[genemuts$gene_name %in% gene_list,]
  }
  ci95df = data.frame(gene_name = genemuts$gene_name, mis_mle=NA, tru_mle=NA, non_mle=NA, spl_mle = NA, mis_low = NA, tru_low=NA, non_low = NA, spl_low = NA, mis_high = NA, tru_high=NA, non_high = NA, spl_high = NA)
  
  # Missense
  ci95df[,c("mis_mle","mis_low","mis_high")] = poissonci(genemuts$n_mis, genemuts$n_syn, genemuts$exp_mis/genemuts$exp_syn)
  # Truncating (nonsense + essential splice sites)
  ci95df[,c("tru_mle","tru_low","tru_high")] = poissonci(genemuts$n_non+genemuts$n_spl, genemuts$n_syn, (genemuts$exp_non+genemuts$exp_spl)/genemuts$exp_syn)
  # Nonsense
  ci95df[,c("non_mle","non_low","non_high")] = poissonci(genemuts$n_non, genemuts$n_syn, genemuts$exp_non/genemuts$exp_syn)
  # Essential splice
  ci95df[,c("spl_mle","spl_low","spl_high")] = poissonci(genemuts$n_spl, genemuts$n_syn, genemuts$exp_spl/genemuts$exp_syn)
  
  return(ci95df)
}

# nb_profile_ci
# Function to calculate profile likelihood confidence intervals for indels

geneindelci <- function(geneindels, gene_list = NULL, theta, level = 0.95) {
  # Ensuring valid level value
  if (level > 1) {
    warning("Confidence level must be lower than 1, using 0.95 as default")
    level = 0.95
  }
  
  # Subfunction to calculate profile likelihood CI
  nb_profile_ci = function(y, theta, conf_level) {
    
    cutoff = dnbinom(y, size = theta, mu = y, log = TRUE) - qchisq(conf_level, 1) / 2
    obj = function(mu) dnbinom(y, size = theta, mu = mu, log = TRUE) - cutoff
    
    # Search for lower bound: Only for y>0 and ensuring that the lower limit is slightly above 0
    if (y == 0) {
      lower = 0
    } else {
      lower = uniroot(obj, interval = c(1e-10, y), extendInt = "no")$root
    }
    
    # Search for upper bound: Start from y and extend outward (bounded to avoid convergence problems for low theta)
    outofbounds = (dnbinom(y, size = theta, mu = y, log = TRUE) - dnbinom(y, size = theta, mu = y * 200 + 10, log = TRUE)) < (qchisq(conf_level, 1) / 2)
    if (!outofbounds) {
      upper = uniroot(obj, interval = c(y, y * 200 + 10), extendInt = "no")$root
    } else {
      upper = Inf
    }
    
    return(c(estimate = y, lower = lower, upper = upper))
  }
  
  # Restrict to gene list
  if (!is.null(gene_list)) {
    geneindels = geneindels[geneindels$gene_name %in% gene_list,]
  }
  ci95df = data.frame(gene_name = geneindels$gene_name, ind_mle=NA, ind_low = NA, ind_high=NA)
  
  for(x in 1:nrow(ci95df)){
    ci95df[x,c("ind_mle","ind_low","ind_high")] = nb_profile_ci(y = geneindels$n_indused[x], theta = theta, conf_level = level) / geneindels$exp_indcv[x]
  }
  
  return(ci95df)
}

# Likelihood-Ratio Test to compare dN/dS ratios per gene between two datasets.
# dnds1 and dnds2 are dndsout objects generated by running dndscv on a list of known cancer genes in two datasets.
#
# In a previous implementation, we used the assumptions of the uniform dNdS model to compare dN/dS ratios between two datasets. This assumed that
# there are no differences in the mutability or duplex coverage of a gene (relative to all other genes in the objects) between both datasets,
# but it corrects for differences in mutation spectra. In the implementation below, we instead use the uniform neutral expectation from the
# sel_cv model, which will include a correction for duplex coverage if this was used in the generation of the dnds1 and dnds2 input objects.

pairwise_dNdScv = function(dnds1, dnds2, genestotest) {
  pvec = rmisvec = rtruvec = setNames(rep(NA,length(genestotest)), genestotest) # Initialising vectors for p-values and for the ratios of wmis and wtru between dataset 1 and 2
  
  genemuts1 = dnds1$genemuts[,-1]; rownames(genemuts1) = dnds1$genemuts$gene_name
  rel_mr = genemuts1$exp_syn_cv/genemuts1$exp_syn
  genemuts1[,c(5:8)] = genemuts1[,c(5:8)] * rel_mr
  
  genemuts2 = dnds2$genemuts[,-1]; rownames(genemuts2) = dnds2$genemuts$gene_name
  rel_mr = genemuts2$exp_syn_cv/genemuts2$exp_syn
  genemuts2[,c(5:8)] = genemuts2[,c(5:8)] * rel_mr
  
  for (g in 1:length(genestotest)) {
    
    # LRT comparing the fit between: (H0) the same dN/dS ratio for a gene in both datasets, and (2) a different dN/dS ratio for each gene in both datasets.
    # H0: wmis1==wmis2 & wtru1==wtru2
    # H1: wmis1!=wmis2 & wtru1!=wtru2
    # This is simply done using obs1, exp1, obs2, exp2 (y1 and y2 vectors below), with exp1 and exp2 being the neutral expected value from the sel_cv models in the input objects.
    
    y1 = genemuts1[genestotest[g],]
    y2 = genemuts2[genestotest[g],]
    
    if (!is.na(y1$exp_syn)) {
      # MLE dN/dS ratios using the uniform model under H0 and H1
      wmis_mle0 = (y1$n_mis+y2$n_mis) / (y1$exp_mis+y2$exp_mis)
      wtru_mle0 = (y1$n_non+y2$n_non+y1$n_spl+y2$n_spl) / (y1$exp_non+y2$exp_non+y1$exp_spl+y2$exp_spl)
      wmis_mle1 = c(y1$n_mis, y2$n_mis) / c(y1$exp_mis, y2$exp_mis)
      wtru_mle1 = c(y1$n_non+y1$n_spl, y2$n_non+y2$n_spl) / c(y1$exp_non+y1$exp_spl, y2$exp_non+y2$exp_spl)
      
      # Observed and predicted counts under H0 and H1
      obs = c(y1$n_mis, y1$n_non+y1$n_spl, y2$n_mis, y2$n_non+y2$n_spl)
      exp0 = c(y1$exp_mis*wmis_mle0, (y1$exp_non+y1$exp_spl)*wtru_mle0, y2$exp_mis*wmis_mle0, (y2$exp_non+y2$exp_spl)*wtru_mle0)
      exp1 = c(y1$exp_mis*wmis_mle1[1], (y1$exp_non+y1$exp_spl)*wtru_mle1[1], y2$exp_mis*wmis_mle1[2], (y2$exp_non+y2$exp_spl)*wtru_mle1[2])
      
      ll0 = c(sum(dpois(x=obs[c(1,3)], lambda=exp0[c(1,3)], log=T)), sum(dpois(x=obs[c(2,4)], lambda=exp0[c(2,4)], log=T)))
      ll1 = c(sum(dpois(x=obs[c(1,3)], lambda=exp1[c(1,3)], log=T)), sum(dpois(x=obs[c(2,4)], lambda=exp1[c(2,4)], log=T)))
      
      # One-sided p-values and then Fisher's method to combine p_mis and p_trunc
      pvals = (1-pchisq(2*(ll1-ll0), df=1))
      if (wmis_mle1[1]<wmis_mle1[2]) { pvals[1] = 1 } else { pvals[1] = pvals[1]/2 }
      if (wtru_mle1[1]<wtru_mle1[2]) { pvals[2] = 1 } else { pvals[2] = pvals[2]/2 }
      
      # Saving the results
      pvec[g] = 1 - pchisq(-2 * sum(log(pvals)), df = 4) # Fisher combined p-value
      rmisvec[g] = wmis_mle1[1]/wmis_mle1[2]
      rtruvec[g] = wtru_mle1[1]/wtru_mle1[2]
    }
  }
  out = data.frame(gene_name=genestotest, rmis=rmisvec, rtru=rtruvec, pval=pvec)
  out$qval = p.adjust(out$pval, method="BH")
  out = out[order(out$pval),]
  return(out)
}

# Driver plot
# The driverdensity_sampleIDs is a vector with the IDs of the samples to include for the estimation of driver % (e.g. individuals of a certain age and/or meeting a minimum duplex depth)
# gene2dc allows to input a different coverage vector (such as gene2dc_blood for the blood dataset)
# logscale refers to the driver density plot

driverplot = function(genes2plot, muts, exome_dnds_sel_cv, exome_dnds_sel_cv_gene_ci = NULL, combined_targeted_only_dnds_sel_cv, combined_targeted_only_dnds_sel_cv_gene_ci = NULL, targeted_dnds_sel_cv, targeted_dnds_sel_cv_gene_ci = NULL, exome_dnds_loc, exome_dnds_loc_gene_ci = NULL, combined_targeted_only_dnds_loc, combined_targeted_only_dnds_loc_gene_ci = NULL, targeted_dnds_loc, targeted_dnds_loc_gene_ci = NULL, min_syn_loc = Inf, gene2dc_array, max_genes = Inf, sortbyfreq = T, onlysignifdnds = T, vafdnds = T, driverdensity_sampleIDs = NULL, full_sampleIDs = NULL, dnds_plot_cap = Inf, nmuts_cap = Inf, nonsyn_mcf_cap = Inf, max_duplex_cov, plotfilename = "Driver_plot.pdf", logscale_dnds = F, logscale_mcf = T, runman = T, plotpanels, plottitles = "", highlighted_genes = NULL) {
  
  if (runman) { dev.new() }
  
  plot_rows <- length(plotpanels)
  plot_cols <- length(genes2plot)
  
  layout(matrix(c(1:(plot_rows*plot_cols)), nrow = plot_rows, ncol = plot_cols, byrow = F), widths = (lengths(genes2plot) + 7.5) / sum(lengths(genes2plot)))
  par(mar = c(5,5.5,2,2))
  
  for(i in 1:length(genes2plot)){
    sample_muts <- muts[muts$sampleID %in% full_sampleIDs[[i]],]
    sample_muts_unique <- sample_muts[which(!(duplicated(paste(sample_muts$donor,sample_muts$chr,sample_muts$pos,sample_muts$ref,sample_muts$mut,sep = "_")))),]
    nmuts <- table(sample_muts_unique$gene,sample_muts_unique$impact)
    
    if (sortbyfreq) {
      genes2plot[[i]] = names(sort(rowSums(nmuts)[genes2plot[[i]]],decreasing=T)) # Sorting the genes in decreasing order of their number of mutations
    }
    if (max_genes<length(genes2plot[[i]])) {
      genes2plot[[i]] = head(genes2plot[[i]],max_genes) # We select the first N genes
    }
    
    nmuts = nmuts[genes2plot[[i]],,drop = F]
    
    # a. Mutations observed in combined exome and targeted data
    if("a" %in% plotpanels){
      subs_per_gene_init = t(nmuts)
      
      subs_per_gene = array(data = 0, dim = c(7,ncol(subs_per_gene_init)))
      rownames(subs_per_gene) <- c("Synonymous","Missense","Start_loss","Stop_loss","Nonsense","Essential_Splice","no-SNV")  
      colnames(subs_per_gene) <- colnames(subs_per_gene_init)
      
      subs_per_gene[rownames(subs_per_gene_init),] <- subs_per_gene_init
      
      rownames(subs_per_gene)[rownames(subs_per_gene) == "Essential_Splice"] <- "Splice"
      rownames(subs_per_gene)[rownames(subs_per_gene) == "no-SNV"] <- "Indels"
      
      colvec = c("Synonymous"="grey70","Missense"="cadetblue","Start_loss"="firebrick","Stop_loss"="red","Nonsense"="darkorchid4","Splice"="darkorchid2","Indels"="chocolate3")
      
      text_labels <- colSums(subs_per_gene)
      text_labels[which(as.numeric(text_labels) < nmuts_cap)] <- ""
      subs_per_gene_scaled <- subs_per_gene
      for(x in 1:ncol(subs_per_gene_scaled)){
        subs_per_gene_scaled[,x] <- subs_per_gene_scaled[,x] * min(sum(subs_per_gene_scaled[,x]),nmuts_cap) / sum(subs_per_gene_scaled[,x])
      }
      if(is.null(highlighted_genes[[i]])){
        pos = barplot(subs_per_gene_scaled, las=2, col=colvec, border=NA, ylim=c(0,if (nmuts_cap == Inf) max(colSums(subs_per_gene))*1.1 else min((1.1 * nmuts_cap),(max(colSums(subs_per_gene))*1.1))), ylab="Total mutations\n(combined unique per donor)", main = plottitles[i])
      }else{
        pos = barplot(subs_per_gene_scaled, las=2, col=colvec, border=NA, ylim=c(0,if (nmuts_cap == Inf) max(colSums(subs_per_gene))*1.1 else min((1.1 * nmuts_cap),(max(colSums(subs_per_gene))*1.1))), ylab="Total mutations\n(combined unique per donor)", main = plottitles[i], names.arg = rep("",times = length(genes2plot[[i]])))
        for(y in 1:length(highlighted_genes[[i]])){
          axis(1, las = 2, at = pos[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], labels = colnames(subs_per_gene_scaled)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
        }
        axis(1, las = 2, at = pos[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], labels = colnames(subs_per_gene_scaled)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
      }
      if(nmuts_cap != Inf){
        text(x = pos, y = nmuts_cap * 1.05, labels = text_labels)
      }
      if(i == 1){
        legend("topright",y=max(apply(subs_per_gene,2,sum))*1.09,legend=rownames(subs_per_gene),fill=colvec,border=NA,box.col=NA,bg="transparent")
      }
    }
    
    # b. Exome data dN/dS ratios from the dNdScv model
    exome_obsw_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),3)))
    colnames(exome_obsw_full) <- c("Missense","Nonsense+splice","Indels")
    rownames(exome_obsw_full) <- genes2plot[[i]]
    
    if(is.null(exome_dnds_sel_cv[[i]])){
      if("b" %in% plotpanels){
        plot.new()
      }
    }else{
      exome_obsw = as.matrix(exome_dnds_sel_cv[[i]][,c("wmis_cv","wnon_cv","wind_cv")])
      rownames(exome_obsw) = exome_dnds_sel_cv[[i]]$gene_name
      colnames(exome_obsw) = c("Missense","Nonsense+splice","Indels")
      exome_obsw_full = exome_obsw[genes2plot[[i]],,drop = F]
      if(!(is.null(exome_dnds_sel_cv_gene_ci[[i]]))){
        rownames(exome_dnds_sel_cv_gene_ci[[i]]) <- exome_dnds_sel_cv_gene_ci[[i]]$gene
        exome_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),6)))
        rownames(exome_ci) <- genes2plot[[i]]
        colnames(exome_ci) <- c("Missense_Low","Missense_High","Nonsense+splice_Low","Nonsense+splice_High","Indels_Low","Indels_High")
        exome_ci <- exome_dnds_sel_cv_gene_ci[[i]][rownames(exome_obsw_full),c("mis_low","mis_high","tru_low","tru_high","ind_low","ind_high")]
      }
      if (onlysignifdnds) {
        obsp = as.matrix(exome_dnds_sel_cv[[i]][,c("pmis_cv","ptrunc_cv","pindpos_cv")]); rownames(obsp) = exome_dnds_sel_cv[[i]]$gene_name; obsp = obsp[genes2plot[[i]],]; exome_obsw_full[obsp>0.05] = NA # Masking out P>0.05
        if(exists("exome_ci")){
          for(j in 1:ncol(obsp)){
            exome_ci[which(obsp[,j]>0.05),c(((2*j)-1):(2*j))] <- NA
          }
        }
      }
      if("b" %in% plotpanels){
        text_labels <- t(exome_obsw_full)
        text_labels <- round(text_labels)
        text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
        if(is.null(highlighted_genes[[i]])){
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(exome_obsw_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(exome_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Exome-wide\ndN/dS sel_cv ratios", log = "y")
          } else{
            pos = barplot(pmin(t(exome_obsw_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(exome_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Exome-wide\ndN/dS sel_cv ratios")
          }
        } else{
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(exome_obsw_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(exome_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Exome-wide\ndN/dS sel_cv ratios", names.arg = rep("",times = 3 * length(genes2plot[[i]])), log = "y")
          } else{
            pos = barplot(pmin(t(exome_obsw_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(exome_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Exome-wide\ndN/dS sel_cv ratios", names.arg = rep("",times = 3 * length(genes2plot[[i]])))
          }
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[2,which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], labels = rownames(exome_obsw_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[2,which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], labels = rownames(exome_obsw_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
        if(dnds_plot_cap != Inf){
          text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3","chocolate3"))
        }
        if(i == 1){
          legend("topright",y=max(apply(exome_obsw_full,2,sum))*1.09,legend=colnames(exome_obsw_full),fill=c("cadetblue","darkorchid3","chocolate3"),border=NA,box.col=NA,bg="transparent")
        }
        if(!(is.null(exome_dnds_sel_cv_gene_ci[[i]]))){
          if(logscale_dnds == T){
            segments(x0=pos, y0=pmax(t(exome_ci[,c("mis_low","tru_low","ind_low")]),0.1), y1=t(exome_ci[,c("mis_high","tru_high","ind_high")]))
          }else{
            segments(x0=pos, y0=t(exome_ci[,c("mis_low","tru_low","ind_low")]), y1=t(exome_ci[,c("mis_high","tru_high","ind_high")])) 
          }
        }
        abline(h=1, col="grey")
      }
    }
    
    # c. Exome data dN/dS ratios from the dNdSloc model
    exome_obsw_loc_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),2)))
    colnames(exome_obsw_loc_full) <- c("Missense","Nonsense+splice")
    rownames(exome_obsw_loc_full) <- genes2plot[[i]]
    
    if(is.null(exome_dnds_loc[[i]])){
      if("c" %in% plotpanels){
        plot.new()
      }
    }else{
      exome_obsw_loc = as.matrix(exome_dnds_loc[[i]][,c("wmis_loc","wnon_loc")])
      rownames(exome_obsw_loc) = exome_dnds_loc[[i]]$gene_name
      colnames(exome_obsw_loc) = c("Missense","Nonsense+splice")
      exome_obsw_loc_full = exome_obsw_loc[genes2plot[[i]],,drop = F]
      if(!(is.null(exome_dnds_loc_gene_ci[[i]]))){
        rownames(exome_dnds_loc_gene_ci[[i]]) <- exome_dnds_loc_gene_ci[[i]]$gene_name
        exome_loc_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),4)))
        rownames(exome_loc_ci) <- genes2plot[[i]]
        colnames(exome_loc_ci) <- c("Missense_Low","Missense_High","Nonsense+splice_Low","Nonsense+splice_High")
        exome_loc_ci <- exome_dnds_loc_gene_ci[[i]][rownames(exome_obsw_loc_full),c("mis_low","mis_high","tru_low","tru_high")]
      }
      if (onlysignifdnds) {
        obsp = as.matrix(exome_dnds_loc[[i]][,c("pmis_loc","ptrunc_loc")]); rownames(obsp) = exome_dnds_loc[[i]]$gene_name; obsp = obsp[genes2plot[[i]],]; exome_obsw_loc_full[obsp>0.05] = NA # Masking out P>0.05
        if(exists("exome_loc_ci")){
          for(j in 1:ncol(obsp)){
            exome_loc_ci[which(obsp[,j]>0.05),c(((2*j)-1):(2*j))] <- NA
          }
        }
      }
      
      syn_count = exome_dnds_loc[[i]][,c("n_syn")]; names(syn_count) = exome_dnds_loc[[i]]$gene_name
      syn_count_full <- rep(TRUE, times = c(length(genes2plot[[i]])))
      names(syn_count_full) <- genes2plot[[i]]
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% names(syn_count)){
          syn_count_full[j] <- syn_count[which(names(syn_count) == genes2plot[[i]][j])]
        }
      }
      exome_obsw_loc_full[which(syn_count_full< min_syn_loc),] = NA
      if(exists("exome_loc_ci")){
        exome_loc_ci[which(syn_count_full < min_syn_loc),] = NA
      }
      if("c" %in% plotpanels){
        if(length(which(!(is.na(exome_obsw_loc_full)))) > 0){
          text_labels <- t(exome_obsw_loc_full)
          text_labels <- round(text_labels)
          text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
          if(is.null(highlighted_genes[[i]])){
            if(logscale_dnds == T){
              pos = barplot(pmax(pmin(t(exome_obsw_loc_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(exome_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Exome-wide\ndN/dS loc (min ",min_syn_loc," syn)"), log = "y")
            }else{
              pos = barplot(pmin(t(exome_obsw_loc_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(exome_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Exome-wide\ndN/dS loc (min ",min_syn_loc," syn)"))
            }
          } else{
            if(logscale_dnds == T){
              pos = barplot(pmax(pmin(t(exome_obsw_loc_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(exome_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Exome-wide\ndN/dS loc (min ",min_syn_loc," syn)"), names.arg = rep("",times = 2 * length(genes2plot[[i]])), log = "y")
            } else{
              pos = barplot(pmin(t(exome_obsw_loc_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(exome_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Exome-wide\ndN/dS loc (min ",min_syn_loc," syn)"), names.arg = rep("",times = 2 * length(genes2plot[[i]])))
            }
            for(y in 1:length(highlighted_genes[[i]])){
              axis(1, las = 2, at = pos[1,which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))] + 0.5, labels = rownames(exome_obsw_loc_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
            }
            axis(1, las = 2, at = pos[1,which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))] + 0.5, labels = rownames(exome_obsw_loc_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
          }
          if(dnds_plot_cap != Inf){
            text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3"))
          }
          if(i == 1){
            legend("topright",y=max(apply(exome_obsw_loc_full,2,sum))*1.09,legend=colnames(exome_obsw_loc_full),fill=c("cadetblue","darkorchid3"),border=NA,box.col=NA,bg="transparent")
          }
          if(!(is.null(exome_dnds_loc_gene_ci[[i]]))){
            if(logscale_dnds == T){
              segments(x0=pos, y0=pmax(t(exome_loc_ci[,c("mis_low","tru_low")]),0.1), y1=t(exome_loc_ci[,c("mis_high","tru_high")]))
            }else{
              segments(x0=pos, y0=t(exome_loc_ci[,c("mis_low","tru_low")]), y1=t(exome_loc_ci[,c("mis_high","tru_high")]))
            }
          }
          abline(h=1, col="grey")
        }else{
          plot.new() 
        }
      }
    }
    
    # d. Targeted data dN/dS ratios from the dNdScv model
    targeted_obsw_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),3)))
    colnames(targeted_obsw_full) <- c("Missense","Nonsense+splice","Indels")
    rownames(targeted_obsw_full) <- genes2plot[[i]]
    
    if(is.null(targeted_dnds_sel_cv[[i]])){
      if("d" %in% plotpanels){
        plot.new()
      }
    }else{
      targeted_obsw = as.matrix(targeted_dnds_sel_cv[[i]][,c("wmis_cv","wnon_cv","wind_cv")])
      rownames(targeted_obsw) = targeted_dnds_sel_cv[[i]]$gene_name
      colnames(targeted_obsw) = c("Missense","Nonsense+splice","Indels")
      
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% rownames(targeted_obsw)){
          targeted_obsw_full[j,] <- targeted_obsw[which(rownames(targeted_obsw) == genes2plot[[i]][j]),,drop = F]
        }
      }
      if(!(is.null(targeted_dnds_sel_cv_gene_ci[[i]]))){
        rownames(targeted_dnds_sel_cv_gene_ci[[i]]) <- targeted_dnds_sel_cv_gene_ci[[i]]$gene
        targeted_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),6)))
        rownames(targeted_ci) <- genes2plot[[i]]
        colnames(targeted_ci) <- c("Missense_Low","Missense_High","Nonsense+splice_Low","Nonsense+splice_High","Indels_Low","Indels_High")
        targeted_ci <- targeted_dnds_sel_cv_gene_ci[[i]][rownames(targeted_obsw_full),c("mis_low","mis_high","tru_low","tru_high","ind_low","ind_high")]
      }
      if (onlysignifdnds) {
        obsp = as.matrix(targeted_dnds_sel_cv[[i]][,c("pmis_cv","ptrunc_cv","pindpos_cv")]); rownames(obsp) = targeted_dnds_sel_cv[[i]]$gene_name; 
        obsp_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),3)))
        colnames(obsp_full) <- colnames(obsp)
        rownames(obsp_full) <- genes2plot[[i]]
        
        for(j in 1:length(genes2plot[[i]])){
          if(genes2plot[[i]][j] %in% rownames(obsp)){
            obsp_full[j,] <- obsp[which(rownames(obsp) == genes2plot[[i]][j]),]
          }
        }
        targeted_obsw_full[obsp_full>0.05] = NA # Masking out P>0.05
        if(exists("targeted_ci")){
          for(j in 1:ncol(obsp)){
            targeted_ci[which(obsp_full[,j]>0.05),c(((2*j)-1):(2*j))] <- NA
          }
        }
      }
      if("d" %in% plotpanels){
        text_labels <- t(targeted_obsw_full)
        text_labels <- round(text_labels)
        text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
        if(is.null(highlighted_genes[[i]])){
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(targeted_obsw_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Targeted panel\ndN/dS sel_cv ratios", log = "y")
          }else{
            pos = barplot(pmin(t(targeted_obsw_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Targeted panel\ndN/dS sel_cv ratios")
          }
        }else{
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(targeted_obsw_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Targeted panel\ndN/dS sel_cv ratios", names.arg = rep("",times = 3 * length(genes2plot[[i]])), log = "y")
          }else{
            pos = barplot(pmin(t(targeted_obsw_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Targeted panel\ndN/dS sel_cv ratios", names.arg = rep("",times = 3 * length(genes2plot[[i]])))
          }
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[2,which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], labels = rownames(targeted_obsw_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[2,which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], labels = rownames(targeted_obsw_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
        if(dnds_plot_cap != Inf){
          text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3","chocolate3"))
        }
        if(i == 1){
          legend("topright",y=max(apply(targeted_obsw_full,2,sum))*1.09,legend=colnames(targeted_obsw_full),fill=c("cadetblue","darkorchid3","chocolate3"),border=NA,box.col=NA,bg="transparent")
        }
        if(!(is.null(targeted_dnds_sel_cv_gene_ci[[i]]))){
          if(logscale_dnds == T){
            segments(x0=pos, y0=pmax(t(targeted_ci[,c("mis_low","tru_low","ind_low")]),0.1), y1=t(targeted_ci[,c("mis_high","tru_high","ind_high")]))
          }else{
            segments(x0=pos, y0=t(targeted_ci[,c("mis_low","tru_low","ind_low")]), y1=t(targeted_ci[,c("mis_high","tru_high","ind_high")]))
          }
        }
        abline(h=1, col="grey")
      }
    }
    
    # e. Targeted data dN/dS ratios from dNdSloc (targeted genes only)
    targeted_obsw_loc_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),2)))
    colnames(targeted_obsw_loc_full) <- c("Missense","Nonsense+splice")
    rownames(targeted_obsw_loc_full) <- genes2plot[[i]]
    
    if(is.null(targeted_dnds_loc[[i]])){
      if("e" %in% plotpanels){
        plot.new()
      }
    }else{
      targeted_obsw_loc = as.matrix(targeted_dnds_loc[[i]][,c("wmis_loc","wnon_loc")])
      rownames(targeted_obsw_loc) = targeted_dnds_loc[[i]]$gene_name
      colnames(targeted_obsw_loc) = c("Missense","Nonsense+splice")
      
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% rownames(targeted_obsw_loc)){
          targeted_obsw_loc_full[j,] <- targeted_obsw_loc[which(rownames(targeted_obsw_loc) == genes2plot[[i]][j]),,drop = F]
        }
      }
      if(!(is.null(targeted_dnds_loc_gene_ci[[i]]))){
        rownames(targeted_dnds_loc_gene_ci[[i]]) <- targeted_dnds_loc_gene_ci[[i]]$gene_name
        targeted_loc_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),4)))
        rownames(targeted_loc_ci) <- genes2plot[[i]]
        colnames(targeted_loc_ci) <- c("Missense_Low","Missense_High","Nonsense+splice_Low","Nonsense+splice_High")
        targeted_loc_ci <- targeted_dnds_loc_gene_ci[[i]][rownames(targeted_obsw_loc_full),c("mis_low","mis_high","tru_low","tru_high")]
      }
      if (onlysignifdnds) {
        obsp = as.matrix(targeted_dnds_loc[[i]][,c("pmis_loc","ptrunc_loc")]); rownames(obsp) = targeted_dnds_loc[[i]]$gene_name; 
        obsp_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),2)))
        colnames(obsp_full) <- colnames(obsp)
        rownames(obsp_full) <- genes2plot[[i]]
        
        for(j in 1:length(genes2plot[[i]])){
          if(genes2plot[[i]][j] %in% rownames(obsp)){
            obsp_full[j,] <- obsp[which(rownames(obsp) == genes2plot[[i]][j]),,drop = F]
          }
        }
        targeted_obsw_loc_full[obsp_full>0.05] = NA # Masking out P>0.05
        if(exists("targeted_loc_ci")){
          for(j in 1:ncol(obsp)){
            targeted_loc_ci[which(obsp_full[,j]>0.05),c(((2*j)-1):(2*j))] <- NA
          }
        }
      }
      
      syn_count = targeted_dnds_loc[[i]][,c("n_syn")]; names(syn_count) = targeted_dnds_loc[[i]]$gene_name
      syn_count_full <- rep(TRUE, times = c(length(genes2plot[[i]])))
      names(syn_count_full) <- genes2plot[[i]]
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% names(syn_count)){
          syn_count_full[j] <- syn_count[which(names(syn_count) == genes2plot[[i]][j])]
        }
      }
      targeted_obsw_loc_full[which(syn_count_full< min_syn_loc),] = NA
      if(exists("targeted_loc_ci")){
        targeted_loc_ci[which(syn_count_full< min_syn_loc),] = NA
      }
      if("e" %in% plotpanels){
        text_labels <- t(targeted_obsw_loc_full)
        text_labels <- round(text_labels)
        text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
        if(is.null(highlighted_genes[[i]])){
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(targeted_obsw_loc_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Targeted panel\ndN/dS loc (min ",min_syn_loc," syn)"), log = "y")
          }else{
            pos = barplot(pmin(t(targeted_obsw_loc_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Targeted panel\ndN/dS loc (min ",min_syn_loc," syn)"))
          }
        }else{
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(targeted_obsw_loc_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Targeted panel\ndN/dS loc (min ",min_syn_loc," syn)"), names.arg = rep("",times = 2 * length(genes2plot[[i]])), log = "y")
          }else{
            pos = barplot(pmin(t(targeted_obsw_loc_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Targeted panel\ndN/dS loc (min ",min_syn_loc," syn)"), names.arg = rep("",times = 2 * length(genes2plot[[i]])))
          }
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[1,which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))] + 0.5, labels = rownames(targeted_obsw_loc_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[1,which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))] + 0.5, labels = rownames(targeted_obsw_loc_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
        if(dnds_plot_cap != Inf){
          text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3"))
        }
        if(i == 1){
          legend("topright",y=max(apply(targeted_obsw_loc_full,2,sum))*1.09,legend=colnames(targeted_obsw_loc_full),fill=c("cadetblue","darkorchid3"),border=NA,box.col=NA,bg="transparent") 
        }
        if(!(is.null(targeted_dnds_loc_gene_ci[[i]]))){
          if(logscale_dnds == T){
            segments(x0=pos, y0=t(targeted_loc_ci[,c("mis_low","tru_low")]), y1=t(targeted_loc_ci[,c("mis_high","tru_high")]))
          }else{
            segments(x0=pos, y0=pmax(t(targeted_loc_ci[,c("mis_low","tru_low")])), y1=t(targeted_loc_ci[,c("mis_high","tru_high")]))
          }
        }
        abline(h=1, col="grey")
      }
    }
    
    # f. Combined exome and targeted data dN/dS ratios from the dNdScv model (targeted genes only)
    combined_targeted_obsw_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),3)))
    colnames(combined_targeted_obsw_full) <- c("Missense","Nonsense+splice","Indels")
    rownames(combined_targeted_obsw_full) <- genes2plot[[i]]
    
    if(is.null(combined_targeted_only_dnds_sel_cv[[i]])){
      if("f" %in% plotpanels){
        plot.new()
      }
    }else{
      combined_targeted_obsw = as.matrix(combined_targeted_only_dnds_sel_cv[[i]][,c("wmis_cv","wnon_cv","wind_cv")])
      rownames(combined_targeted_obsw) = combined_targeted_only_dnds_sel_cv[[i]]$gene_name
      colnames(combined_targeted_obsw) = c("Missense","Nonsense+splice","Indels")
      
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% rownames(combined_targeted_obsw)){
          combined_targeted_obsw_full[j,] <- combined_targeted_obsw[which(rownames(combined_targeted_obsw) == genes2plot[[i]][j]),]
        }
      }
      if(!(is.null(combined_targeted_only_dnds_sel_cv_gene_ci[[i]]))){
        rownames(combined_targeted_only_dnds_sel_cv_gene_ci[[i]]) <- combined_targeted_only_dnds_sel_cv_gene_ci[[i]]$gene
        combined_targeted_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),6)))
        rownames(combined_targeted_ci) <- genes2plot[[i]]
        colnames(combined_targeted_ci) <- c("Missense_Low","Missense_High","Nonsense+splice_Low","Nonsense+splice_High","Indels_Low","Indels_High")
        combined_targeted_ci <- combined_targeted_only_dnds_sel_cv_gene_ci[[i]][rownames(combined_targeted_obsw_full),c("mis_low","mis_high","tru_low","tru_high","ind_low","ind_high")]
      }
      if (onlysignifdnds) {
        obsp = as.matrix(combined_targeted_only_dnds_sel_cv[[i]][,c("pmis_cv","ptrunc_cv","pindpos_cv")]); rownames(obsp) = combined_targeted_only_dnds_sel_cv[[i]]$gene_name; 
        obsp_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),3)))
        colnames(obsp_full) <- colnames(obsp)
        rownames(obsp_full) <- genes2plot[[i]]
        
        for(j in 1:length(genes2plot[[i]])){
          if(genes2plot[[i]][j] %in% rownames(obsp)){
            obsp_full[j,] <- obsp[which(rownames(obsp) == genes2plot[[i]][j]),]
          }
        }
        combined_targeted_obsw_full[obsp_full>0.05] = NA # Masking out P>0.05
        if(exists("combined_targeted_ci")){
          for(j in 1:ncol(obsp)){
            combined_targeted_ci[which(obsp_full[,j]>0.05),c(((2*j)-1):(2*j))] <- NA
          }
        }
      }
      if("f" %in% plotpanels){
        text_labels <- t(combined_targeted_obsw_full)
        text_labels <- round(text_labels)
        text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
        if(is.null(highlighted_genes[[i]])){
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(combined_targeted_obsw_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Combined (targeted panel)\ndN/dS sel_cv ratios", log = "y")           
          }else{
            pos = barplot(pmin(t(combined_targeted_obsw_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Combined (targeted panel)\ndN/dS sel_cv ratios")
          }
        }else{
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(combined_targeted_obsw_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Combined (targeted panel)\ndN/dS sel_cv ratios", names.arg = rep("",times = 3 * length(genes2plot[[i]])), log = "y")
          }else{
            pos = barplot(pmin(t(combined_targeted_obsw_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Combined (targeted panel)\ndN/dS sel_cv ratios", names.arg = rep("",times = 3 * length(genes2plot[[i]])))
          }
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[2,which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], labels = rownames(combined_targeted_obsw_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[2,which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], labels = rownames(combined_targeted_obsw_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
        if(dnds_plot_cap != Inf){
          text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3","chocolate3"))
        }
        if(i == 1){
          legend("topright",y=max(apply(combined_targeted_obsw_full,2,sum))*1.09,legend=colnames(combined_targeted_obsw_full),fill=c("cadetblue","darkorchid3","chocolate3"),border=NA,box.col=NA,bg="transparent")
        }
        if(!(is.null(combined_targeted_only_dnds_sel_cv_gene_ci[[i]]))){
          if(logscale_dnds == T){
            segments(x0=pos, y0=pmax(t(combined_targeted_ci[,c("mis_low","tru_low","ind_low")])), y1=t(combined_targeted_ci[,c("mis_high","tru_high","ind_high")]))
          }else{
            segments(x0=pos, y0=t(combined_targeted_ci[,c("mis_low","tru_low","ind_low")]), y1=t(combined_targeted_ci[,c("mis_high","tru_high","ind_high")]))
          }
        }
        abline(h=1, col="grey")
      }
    }
    
    # g. Combined exome and targeted data dN/dS ratios from dNdSloc (targeted genes only)
    combined_targeted_obsw_loc_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),2)))
    colnames(combined_targeted_obsw_loc_full) <- c("Missense","Nonsense+splice")
    rownames(combined_targeted_obsw_loc_full) <- genes2plot[[i]]
    
    if(is.null(combined_targeted_only_dnds_loc[[i]])){
      if("g" %in% plotpanels){
        plot.new()
      }
    }else{
      combined_targeted_obsw_loc = as.matrix(combined_targeted_only_dnds_loc[[i]][,c("wmis_loc","wnon_loc")])
      rownames(combined_targeted_obsw_loc) = combined_targeted_only_dnds_loc[[i]]$gene_name
      colnames(combined_targeted_obsw_loc) = c("Missense","Nonsense+splice")
      
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% rownames(combined_targeted_obsw_loc)){
          combined_targeted_obsw_loc_full[j,] <- combined_targeted_obsw_loc[which(rownames(combined_targeted_obsw_loc) == genes2plot[[i]][j]),]
        }
      }
      if(!(is.null(combined_targeted_only_dnds_loc_gene_ci[[i]]))){
        rownames(combined_targeted_only_dnds_loc_gene_ci[[i]]) <- combined_targeted_only_dnds_loc_gene_ci[[i]]$gene_name
        combined_targeted_loc_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),4)))
        rownames(combined_targeted_loc_ci) <- genes2plot[[i]]
        colnames(combined_targeted_loc_ci) <- c("Missense_Low","Missense_High","Nonsense+splice_Low","Nonsense+splice_High")
        combined_targeted_loc_ci <- combined_targeted_only_dnds_loc_gene_ci[[i]][rownames(combined_targeted_obsw_loc_full),c("mis_low","mis_high","tru_low","tru_high")]
      }
      if (onlysignifdnds) {
        obsp = as.matrix(combined_targeted_only_dnds_loc[[i]][,c("pmis_loc","ptrunc_loc")]); rownames(obsp) = combined_targeted_only_dnds_loc[[i]]$gene_name; 
        obsp_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),2)))
        colnames(obsp_full) <- colnames(obsp)
        rownames(obsp_full) <- genes2plot[[i]]
        
        for(j in 1:length(genes2plot[[i]])){
          if(genes2plot[[i]][j] %in% rownames(obsp)){
            obsp_full[j,] <- obsp[which(rownames(obsp) == genes2plot[[i]][j]),]
          }
        }
        combined_targeted_obsw_loc_full[obsp_full>0.05] = NA # Masking out P>0.05
        if(exists("combined_targeted_loc_ci")){
          for(j in 1:ncol(obsp)){
            combined_targeted_loc_ci[which(obsp_full[,j]>0.05),c(((2*j)-1):(2*j))] <- NA
          }
        }
      }
      
      syn_count = combined_targeted_only_dnds_loc[[i]][,c("n_syn")]; names(syn_count) = combined_targeted_only_dnds_loc[[i]]$gene_name
      syn_count_full <- rep(TRUE, times = c(length(genes2plot[[i]])))
      names(syn_count_full) <- genes2plot[[i]]
      for(j in 1:length(genes2plot[[i]])){
        if(genes2plot[[i]][j] %in% names(syn_count)){
          syn_count_full[j] <- syn_count[which(names(syn_count) == genes2plot[[i]][j])]
        }
      }
      combined_targeted_obsw_loc_full[which(syn_count_full< min_syn_loc),] = NA
      if(exists("combined_targeted_loc_ci")){
        combined_targeted_loc_ci[which(syn_count_full< min_syn_loc),] = NA
      }
      if("g" %in% plotpanels){
        text_labels <- t(combined_targeted_obsw_loc_full)
        text_labels <- round(text_labels)
        text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
        if(is.null(highlighted_genes[[i]])){
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(combined_targeted_obsw_loc_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Combined (targeted panel)\ndN/dS loc (min ",min_syn_loc," syn)"), log = "y")
          }else{
            pos = barplot(pmin(t(combined_targeted_obsw_loc_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Combined (targeted panel)\ndN/dS loc (min ",min_syn_loc," syn)"))
          }
        } else{
          if(logscale_dnds == T){
            pos = barplot(pmax(pmin(t(combined_targeted_obsw_loc_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Combined (targeted panel)\ndN/dS loc (min ",min_syn_loc," syn)"), names.arg = rep("",times = 2 * length(genes2plot[[i]])), log = "y")
          }else{
            pos = barplot(pmin(t(combined_targeted_obsw_loc_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(combined_targeted_obsw_loc_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab=paste0("Combined (targeted panel)\ndN/dS loc (min ",min_syn_loc," syn)"), names.arg = rep("",times = 2 * length(genes2plot[[i]])))
          }
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[1,which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))] + 0.5, labels = rownames(combined_targeted_obsw_loc_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[1,which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))] + 0.5, labels = rownames(combined_targeted_obsw_loc_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
        if(dnds_plot_cap != Inf){
          text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3"))
        }
        if(i == 1){
          legend("topright",y=max(apply(combined_targeted_obsw_loc_full,2,sum))*1.09,legend=colnames(combined_targeted_obsw_loc_full),fill=c("cadetblue","darkorchid3"),border=NA,box.col=NA,bg="transparent") 
        }
        if(!(is.null(combined_targeted_only_dnds_loc_gene_ci[[i]]))){
          if(logscale_dnds == T){
            segments(x0=pos, y0=pmax(t(combined_targeted_loc_ci[,c("mis_low","tru_low")]),0.1), y1=t(combined_targeted_loc_ci[,c("mis_high","tru_high")]))
          }else{
            segments(x0=pos, y0=t(combined_targeted_loc_ci[,c("mis_low","tru_low")]), y1=t(combined_targeted_loc_ci[,c("mis_high","tru_high")]))
          }
        }
        abline(h=1, col="grey")
      }
    }
    
    # h. Combined dN/dS ratios (sel_loc for subs for genes with more than threshold number of synonymous mutations & sel_cv for other genes)
    if("h" %in% plotpanels){
      merged_dnds_values_full <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),3)))
      colnames(merged_dnds_values_full) <- c("Missense","Nonsense+splice","Indels")
      rownames(merged_dnds_values_full) <- genes2plot[[i]]
      
      if(!(onlysignifdnds)){
        merge_order <- c("combined_targeted_obsw_loc_full","exome_obsw_loc_full","targeted_obsw_loc_full","combined_targeted_obsw_full","exome_obsw_full","targeted_obsw_full")
      }else if(!(all(is.na(combined_targeted_obsw_full)))){
        merge_order <- c("combined_targeted_obsw_loc_full","combined_targeted_obsw_full")
      }else if(!(all(is.na(exome_obsw_full)))){
        merge_order <- c("exome_obsw_loc_full","exome_obsw_full")
      }else if(!(all(is.na(targeted_obsw_full)))){
        merge_order <- c("targeted_obsw_loc_full","targeted_obsw_full")
      }
            
      for(x in 1:length(merge_order)){
        for(y in 1:nrow(merged_dnds_values_full)){
          for(z in 1:ncol(merged_dnds_values_full)){
            if(is.na(merged_dnds_values_full[y,z])){
              if(colnames(merged_dnds_values_full)[z] %in% colnames(eval(parse(text = merge_order[x])))){
                if(!(is.na(eval(parse(text = merge_order[x]))[y,z]))){
                  merged_dnds_values_full[y,z] <- eval(parse(text = merge_order[x]))[y,z]
                }
              }
            }
          }
        }
      }
      
      if(any(!(is.null(exome_dnds_sel_cv_gene_ci)),!(is.null(exome_dnds_loc_gene_ci)),!(is.null(targeted_dnds_sel_cv_gene_ci)),!(is.null(targeted_dnds_loc_gene_ci)),!(is.null(combined_targeted_only_dnds_sel_cv_gene_ci)),!(is.null(combined_targeted_only_dnds_loc_gene_ci)))){
        if(!(onlysignifdnds)){
          merge_ci_order <- c("combined_targeted_loc_ci","exome_loc_ci","targeted_loc_ci","combined_targeted_ci","exome_ci","targeted_ci")
        }else if(exists("combined_targeted_ci")){
          merge_ci_order <- c("combined_targeted_loc_ci","combined_targeted_ci")
        }else if(exists("exome_ci")){
          merge_ci_order <- c("exome_loc_ci","exome_ci")
        }else if(exists("targeted_ci")){
          merge_ci_order <- c("targeted_loc_ci","targeted_ci")
        }
        
        merged_ci <- as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]),6)))
        rownames(merged_ci) <- genes2plot[[i]]
        colnames(merged_ci) <- c("mis_low","mis_high","tru_low","tru_high","ind_low","ind_high")
        
        for(x in 1:length(merge_ci_order)){
          for(y in 1:nrow(merged_ci)){
            for(z in 1:ncol(merged_ci)){
              if(is.na(merged_ci[y,z])){
                if(exists(merge_ci_order[x])){
                  if(colnames(merged_ci)[z] %in% colnames(eval(parse(text = merge_ci_order[x])))){
                    if(!(is.na(eval(parse(text = merge_ci_order[x]))[y,z]))){
                      merged_ci[y,z] <- eval(parse(text = merge_ci_order[x]))[y,z]
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      text_labels <- t(merged_dnds_values_full)
      text_labels <- round(text_labels)
      text_labels[which(as.numeric(text_labels) < dnds_plot_cap)] <- ""
      if(is.null(highlighted_genes[[i]])){
        if(logscale_dnds == T){
          pos = barplot(pmax(pmin(t(merged_dnds_values_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(merged_dnds_values_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Merged\ndN/dS values", log = "y")
        } else{
          pos = barplot(pmin(t(merged_dnds_values_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(merged_dnds_values_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Merged\ndN/dS values")
        }
      } else{
        if(logscale_dnds == T){
          pos = barplot(pmax(pmin(t(merged_dnds_values_full), dnds_plot_cap),0.1), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0.1,if (dnds_plot_cap == Inf) max(merged_dnds_values_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Merged\ndN/dS values", names.arg = rep("",times = 3 * length(genes2plot[[i]])), log = "y")
        } else{
          pos = barplot(pmin(t(merged_dnds_values_full), dnds_plot_cap), beside=T, las=2, col=c("cadetblue","darkorchid3","chocolate3"), border=NA, ylim=c(0,if (dnds_plot_cap == Inf) max(merged_dnds_values_full,na.rm=T)+10 else 1.1 * dnds_plot_cap), ylab="Merged\ndN/dS values", names.arg = rep("",times = 3 * length(genes2plot[[i]])))
        }
        for(y in 1:length(highlighted_genes[[i]])){
          axis(1, las = 2, at = pos[(which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y])) * 3)] - 1, labels = rownames(merged_dnds_values_full)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
        }
        axis(1, las = 2, at = pos[(which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))* 3)] - 1, labels = rownames(merged_dnds_values_full)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
      }
      if(dnds_plot_cap != Inf){
        text(x = pos, y = dnds_plot_cap * 1.05, labels = text_labels, col=c("cadetblue","darkorchid3","chocolate3"))
      }
      if(i == 1){
        legend("topright",y=max(apply(merged_dnds_values_full,2,sum))*1.09,legend=colnames(merged_dnds_values_full),fill=c("cadetblue","darkorchid3","chocolate3"),border=NA,box.col=NA,bg="transparent")
      }
      if(exists("merged_ci")){
        if(logscale_dnds == T){
          segments(x0=pos, y0=pmax(t(merged_ci[,c("mis_low","tru_low","ind_low")]),0.1), y1=t(merged_ci[,c("mis_high","tru_high","ind_high")]))
        }else{
          segments(x0=pos, y0=t(merged_ci[,c("mis_low","tru_low","ind_low")]), y1=t(merged_ci[,c("mis_high","tru_high","ind_high")]))
        }
      }
      abline(h=1, col="grey")
    }
    
    # i. % of mutant cells with non-synonymous mutations
    maux = muts
    maux = maux[which(maux$sampleID %in% driverdensity_sampleIDs[[i]]), ,drop = F]
    cellfraction = as.data.frame(array(NA, dim = c(length(genes2plot[[i]]),8), dimnames = list(genes2plot[[i]], c("mis","non","spl","ind","mislow","nonlow","spllow","indlow"))))
    numsamples = length(unique(maux$sampleID))
    
    for (j in 1:length(genes2plot[[i]])) {
      aux = maux[which(maux$gene==genes2plot[[i]][j]),c("sampleID","impact","cellfraction","duplex_vaf"), drop = F]
      per_sample_gene_cov <- gene2dc_array[which(gene2dc_array$gene == genes2plot[[i]][j]),which(colnames(gene2dc_array) %in% driverdensity_sampleIDs[[i]]),drop = F]
      aux$sample_gene_cov <- NA
      aux$cov_weighting <- NA
      for(k in 1:length(per_sample_gene_cov)){
        aux$sample_gene_cov[which(aux$sampleID == names(per_sample_gene_cov)[k])] <- as.numeric(per_sample_gene_cov[k])
      }
      aux$cov_weighting <- aux$sample_gene_cov / sum(per_sample_gene_cov)
      
      aux$weighted_cellfraction <- aux$cellfraction * aux$cov_weighting
      aux$weighted_duplex_vaf <- aux$duplex_vaf * aux$cov_weighting
      
      cellfraction[j,c("mis","mislow")] = colSums(aux[aux$impact %in% c("Missense","Start_loss","Stop_loss"),c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
      cellfraction[j,c("non","nonlow")] = colSums(aux[aux$impact=="Nonsense",c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
      cellfraction[j,c("spl","spllow")] = colSums(aux[aux$impact=="Essential_Splice",c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
      cellfraction[j,c("ind","indlow")] = colSums(aux[aux$impact=="no-SNV",c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
    }
    
    nonsyndens = cellfraction * 100
    nonsyndens_bounds = cbind(up = rowSums(nonsyndens[,1:4]), low=rowSums(nonsyndens[,5:8]))
    aux = t(cbind(nonsyndens_bounds[,2],nonsyndens_bounds[,1]-nonsyndens_bounds[,2]))
    colnames(aux) <- rownames(nonsyndens_bounds)
    
    aux[which(aux == 0)] <- 0.001
    
    if("i" %in% plotpanels){
      if (logscale_mcf == T) {
        if(is.null(highlighted_genes[[i]])){
          pos = barplot(aux[2,], las=2, col=c("indianred3"), border=c("indianred3"), ylab="% cells non-synonymous mutation", ylim=c(0.01,if (nonsyn_mcf_cap == Inf) 2*max(aux,na.rm=T)+2 else 1.1 * nonsyn_mcf_cap), log="y", offset = aux[1,])
        } else{
          pos = barplot(aux[2,], las=2, col=c("indianred3"), border=c("indianred3"), ylab="% cells non-synonymous mutation", ylim=c(0.01,if (nonsyn_mcf_cap == Inf) 2*max(aux,na.rm=T)+2 else 1.1 * nonsyn_mcf_cap), log="y", offset = aux[1,],names.arg = rep("",times = length(genes2plot[[i]])))
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], labels = colnames(aux)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], labels = colnames(aux)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
      } else {
        if(is.null(highlighted_genes[[i]])){
          pos = barplot(aux[2,], las=2, col=c("indianred3"), border=c("indianred3"), ylab="% cells non-synonymous mutation", ylim=c(0,if (nonsyn_mcf_cap == Inf) 2*max(aux,na.rm=T)+2 else 1.1 * nonsyn_mcf_cap), offset = aux[1,])
        } else{
          pos = barplot(aux[2,], las=2, col=c("indianred3"), border=c("indianred3"), ylab="% cells non-synonymous mutation", ylim=c(0,if (nonsyn_mcf_cap == Inf) 2*max(aux,na.rm=T)+2 else 1.1 * nonsyn_mcf_cap), offset = aux[1,],names.arg = rep("",times = length(genes2plot[[i]])))
          for(y in 1:length(highlighted_genes[[i]])){
            axis(1, las = 2, at = pos[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], labels = colnames(aux)[which(genes2plot[[i]] %in% names(highlighted_genes[[i]][y]))], tick = F, col.axis = highlighted_genes[[i]][y])
          }
          axis(1, las = 2, at = pos[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], labels = colnames(aux)[which(!(genes2plot[[i]] %in% names(highlighted_genes[[i]])))], tick = F, col.axis = "black")
        }
      }
    }
    
    # j. additional information about each gene
    add_info <- as.data.frame(array(data = NA, dim = c(7,length(genes2plot[[i]]))))
    colnames(add_info) <- genes2plot[[i]]
    rownames(add_info) <- c("ExomeDuplexCov","ExomeWideSig","ExomeWideNom","TargetedPanel","TargetedDuplexCov","CombinedExomeTargetedSig","CombinedExomeTargetedNom")
    add_info["TargetedPanel",] <- genes2plot[[i]] %in% immune_genes$gene
    for(x in 1:ncol(add_info)){
      add_info["ExomeDuplexCov",x] <- sum(gene2dc_array[which(gene2dc_array$gene == colnames(add_info)[x]),which(colnames(gene2dc_array) %in% full_sampleIDs[[i]] & grepl("exome",colnames(gene2dc_array)))])
      add_info["TargetedDuplexCov",x] <- sum(gene2dc_array[which(gene2dc_array$gene == colnames(add_info)[x]),which(colnames(gene2dc_array) %in% full_sampleIDs[[i]] & grepl("targeted",colnames(gene2dc_array)))])
      if(!(is.null(exome_dnds_sel_cv[[i]]))){
        if(length(which(colnames(exome_dnds_sel_cv[[i]]) == "qglobalpos_m")) == 1){
          add_info["ExomeWideSig",x] <- exome_dnds_sel_cv[[i]]$qglobalpos_m[which(exome_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["ExomeWideNom",x] <- exome_dnds_sel_cv[[i]]$pglobalpos_m[which(exome_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }else if(length(which(colnames(exome_dnds_sel_cv[[i]]) == "qglobal_m")) == 1){
          add_info["ExomeWideSig",x] <- exome_dnds_sel_cv[[i]]$qglobal_m[which(exome_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["ExomeWideNom",x] <- exome_dnds_sel_cv[[i]]$pglobal_m[which(exome_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }else{
          add_info["ExomeWideSig",x] <- exome_dnds_sel_cv[[i]]$qglobalpos_cv[which(exome_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["ExomeWideNom",x] <- exome_dnds_sel_cv[[i]]$pglobalpos_cv[which(exome_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }
      }
      if(add_info["TargetedPanel",x] == 1 & !(is.null(targeted_dnds_sel_cv[[i]]))){
        if(length(which(colnames(targeted_dnds_sel_cv[[i]]) == "qglobalpos_m")) == 1){
          add_info["TargetedSig",x] <- targeted_dnds_sel_cv[[i]]$qglobalpos_m[which(targeted_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["TargetedNom",x] <- targeted_dnds_sel_cv[[i]]$pglobalpos_m[which(targeted_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }else if(length(which(colnames(targeted_dnds_sel_cv[[i]]) == "qglobal_m")) == 1){
          add_info["TargetedSig",x] <- targeted_dnds_sel_cv[[i]]$qglobal_m[which(targeted_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["TargetedNom",x] <- targeted_dnds_sel_cv[[i]]$pglobal_m[which(targeted_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }else{
          add_info["TargetedSig",x] <- targeted_dnds_sel_cv[[i]]$qglobalpos_cv[which(targeted_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["TargetedNom",x] <- targeted_dnds_sel_cv[[i]]$pglobalpos_cv[which(targeted_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }
      }
      if(add_info["TargetedPanel",x] == 1 & !(is.null(combined_targeted_only_dnds_sel_cv[[i]]))){
        if(length(which(colnames(combined_targeted_only_dnds_sel_cv[[i]]) == "qglobalpos_m")) == 1){
          add_info["CombinedExomeTargetedSig",x] <- combined_targeted_only_dnds_sel_cv[[i]]$qglobalpos_m[which(combined_targeted_only_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["CombinedExomeTargetedNom",x] <- combined_targeted_only_dnds_sel_cv[[i]]$pglobalpos_m[which(combined_targeted_only_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }else if(length(which(colnames(combined_targeted_only_dnds_sel_cv[[i]]) == "qglobal_m")) == 1){
          add_info["CombinedExomeTargetedSig",x] <- combined_targeted_only_dnds_sel_cv[[i]]$qglobal_m[which(combined_targeted_only_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["CombinedExomeTargetedNom",x] <- combined_targeted_only_dnds_sel_cv[[i]]$pglobal_m[which(combined_targeted_only_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }else{
          add_info["CombinedExomeTargetedSig",x] <- combined_targeted_only_dnds_sel_cv[[i]]$qglobalpos_cv[which(combined_targeted_only_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
          add_info["CombinedExomeTargetedNom",x] <- combined_targeted_only_dnds_sel_cv[[i]]$pglobalpos_cv[which(combined_targeted_only_dnds_sel_cv[[i]]$gene_name == colnames(add_info)[x])]
        }
      }
    }
    
    if("j" %in% plotpanels){
      pal_cov = mako(n = 100, direction = -1)
      pal_sig <- brewer.pal(n = 5, name = "YlOrRd")
      plot(x = c(0, length(genes2plot[[i]])), y = c(0, 5), type= "n", xlab = "", ylab = "", axes = F)
      
      for(x in 1:ncol(add_info)){
        if(add_info["ExomeDuplexCov",x] == 0){
          color = "grey"
        }else{
          color <- pal_cov[round((add_info["ExomeDuplexCov",x] / max_duplex_cov) * 100)]
        }
        rect(xleft = (x - 1), ybottom = 4.5, xright = x, ytop = 5, border = "white", col=color)
      }
      if(length(which(add_info["ExomeWideSig",] <= 0.0001) > 0)){
        rect(xleft = which(add_info["ExomeWideSig",] <= 0.0001) - 1, ybottom = 4, xright = which(add_info["ExomeWideSig",] <= 0.0001), ytop = 4.5, border = "white", col=pal_sig[5])
        text(x = rep(which(add_info["ExomeWideSig",] <= 0.0001), each = 4) - rep(rep(c(2/3,1/3), each = 2), times = length(which(add_info["ExomeWideSig",] < 0.0001))), y = rep(rep(c((4.5 - 1/6),(4.5 - 1/3)), times = 2), times = length(which(add_info["ExomeWideSig",] < 0.0001))), labels = "*")
      }
      if(length(which(add_info["ExomeWideSig",] <= 0.001 & add_info["ExomeWideSig",] > 0.0001)) > 0){
        rect(xleft = which(add_info["ExomeWideSig",] <= 0.001 & add_info["ExomeWideSig",] > 0.0001) - 1, ybottom = 4, xright = which(add_info["ExomeWideSig",] <= 0.001 & add_info["ExomeWideSig",] > 0.0001), ytop = 4.5, border = "white", col=pal_sig[4])
        text(x = rep(which(add_info["ExomeWideSig",] <= 0.001 & add_info["ExomeWideSig",] > 0.0001), each = 3) - rep(c(1/2,2/3,1/3), times = length(which(add_info["ExomeWideSig",] <= 0.001 & add_info["ExomeWideSig",] > 0.0001))), y = rep(c((4.5 - 1/6),(4.5 - 1/3),(4.5 - 1/3)), times = length(which(add_info["ExomeWideSig",] <= 0.001 & add_info["ExomeWideSig",] > 0.0001))), labels = "*")
      }
      if(length(which(add_info["ExomeWideSig",] <= 0.01 & add_info["ExomeWideSig",] > 0.001)) > 0){
        rect(xleft = which(add_info["ExomeWideSig",] <= 0.01 & add_info["ExomeWideSig",] > 0.001) - 1, ybottom = 4, xright = which(add_info["ExomeWideSig",] <= 0.01 & add_info["ExomeWideSig",] > 0.001), ytop = 4.5, border = "white", col=pal_sig[3])
        text(x = rep(which(add_info["ExomeWideSig",] <= 0.01 & add_info["ExomeWideSig",] > 0.001), each = 2) - rep(c(2/3,1/3), times = length(which(add_info["ExomeWideSig",] <= 0.01 & add_info["ExomeWideSig",] > 0.001))), y = rep(c((4.5 - 1/4),(4.5 - 1/4)), times = length(which(add_info["ExomeWideSig",] <= 0.01 & add_info["ExomeWideSig",] > 0.001))), labels = "*")
      }
      if(length(which(add_info["ExomeWideSig",] <= 0.1 & add_info["ExomeWideSig",] > 0.01)) > 0){
        rect(xleft = which(add_info["ExomeWideSig",] <= 0.1 & add_info["ExomeWideSig",] > 0.01) - 1, ybottom = 4, xright = which(add_info["ExomeWideSig",] <= 0.1 & add_info["ExomeWideSig",] > 0.01), ytop = 4.5, border = "white", col=pal_sig[2])
        text(x = which(add_info["ExomeWideSig",] <= 0.1 & add_info["ExomeWideSig",] > 0.01) - 0.5, y = rep((4.5 - 1/4), times = length(which(add_info["ExomeWideSig",] <= 0.1 & add_info["ExomeWideSig",] > 0.01))), labels = "*")
      }
      if(length(which(add_info["ExomeWideSig",] > 0.1 & add_info["ExomeWideNom",] <= 0.01)) > 0){
        rect(xleft = which(add_info["ExomeWideSig",] > 0.1 & add_info["ExomeWideNom",] <= 0.01) - 1, ybottom = 4, xright = which(add_info["ExomeWideSig",] > 0.1 & add_info["ExomeWideNom",] <= 0.01), ytop = 4.5, border = "white", col=pal_sig[1])
        text(x = which(add_info["ExomeWideSig",] > 0.1 & add_info["ExomeWideNom",] <= 0.01) - 0.5, y = rep((4.5 - 1/4), times = length(which(add_info["ExomeWideSig",] > 0.1 & add_info["ExomeWideNom",] <= 0.01))), labels = ".")
      }
      if(length(which(is.na(add_info["ExomeWideSig",])))){
        rect(xleft = which(is.na(add_info["ExomeWideSig",])) - 1, ybottom = 4, xright = which(is.na(add_info["ExomeWideSig",])), ytop = 4.5, border = "white", col="grey")
      }
      
      for(x in 1:ncol(add_info)){
        if(add_info["TargetedDuplexCov",x] == 0){
          color = "grey"
        }else{
          color <- pal_cov[round((add_info["TargetedDuplexCov",x] / max_duplex_cov) * 100)]
        }
        rect(xleft = (x - 1), ybottom = 3, xright = x, ytop = 3.5, border = "white", col=color)
      }
      
      if(length(which(add_info["TargetedSig",] <= 0.0001)) > 0){
        rect(xleft = which(add_info["TargetedSig",] <= 0.0001) - 1, ybottom = 2.5, xright = which(add_info["TargetedSig",] <= 0.0001), ytop = 3, border = "white", col=pal_sig[5])
        text(x = rep(which(add_info["TargetedSig",] <= 0.0001), each = 4) - rep(rep(c(2/3,1/3), each = 2), times = length(which(add_info["TargetedSig",] < 0.0001))), y = rep(rep(c((3 - 1/6),(3 - 1/3)), times = 2), times = length(which(add_info["TargetedSig",] < 0.0001))), labels = "*")
      }
      if(length(which(add_info["TargetedSig",] <= 0.001 & add_info["TargetedSig",] > 0.0001)) > 0){
        rect(xleft = which(add_info["TargetedSig",] <= 0.001 & add_info["TargetedSig",] > 0.0001) - 1, ybottom = 2.5, xright = which(add_info["TargetedSig",] <= 0.001 & add_info["TargetedSig",] > 0.0001), ytop = 3, border = "white", col=pal_sig[4])
        text(x = rep(which(add_info["TargetedSig",] <= 0.001 & add_info["TargetedSig",] > 0.0001), each = 3) - rep(c(1/2,2/3,1/3), times = length(which(add_info["TargetedSig",] <= 0.001 & add_info["TargetedSig",] > 0.0001))), y = rep(c((3 - 1/6),(3 - 1/3),(3 - 1/3)), times = length(which(add_info["TargetedSig",] <= 0.001 & add_info["TargetedSig",] > 0.0001))), labels = "*")
      }
      if(length(which(add_info["TargetedSig",] <= 0.01 & add_info["TargetedSig",] > 0.001)) > 0){
        rect(xleft = which(add_info["TargetedSig",] <= 0.01 & add_info["TargetedSig",] > 0.001) - 1, ybottom = 2.5, xright = which(add_info["TargetedSig",] <= 0.01 & add_info["TargetedSig",] > 0.001), ytop = 3, border = "white", col=pal_sig[3])
        text(x = rep(which(add_info["TargetedSig",] <= 0.01 & add_info["TargetedSig",] > 0.001), each = 2) - rep(c(2/3,1/3), times = length(which(add_info["TargetedSig",] <= 0.01 & add_info["TargetedSig",] > 0.001))), y = rep(c((3 - 1/4),(3 - 1/4)), times = length(which(add_info["TargetedSig",] <= 0.01 & add_info["TargetedSig",] > 0.001))), labels = "*")
      }
      if(length(which(add_info["TargetedSig",] <= 0.1 & add_info["TargetedSig",] > 0.01)) > 0){
        rect(xleft = which(add_info["TargetedSig",] <= 0.1 & add_info["TargetedSig",] > 0.01) - 1, ybottom = 2.5, xright = which(add_info["TargetedSig",] <= 0.1 & add_info["TargetedSig",] > 0.01), ytop = 3, border = "white", col=pal_sig[2])
        text(x = which(add_info["TargetedSig",] <= 0.1 & add_info["TargetedSig",] > 0.01) - 0.5, y = rep((3 - 1/4), times = length(which(add_info["TargetedSig",] <= 0.1 & add_info["TargetedSig",] > 0.01))), labels = "*")
      }
      if(length(which(add_info["TargetedSig",] > 0.1 & add_info["TargetedNom",] <= 0.01)) > 0){
        rect(xleft = which(add_info["TargetedSig",] > 0.1 & add_info["TargetedNom",] <= 0.01) - 1, ybottom = 2.5, xright = which(add_info["TargetedSig",] > 0.1 & add_info["TargetedNom",] <= 0.01), ytop = 3, border = "white", col=pal_sig[1])
        text(x = which(add_info["TargetedSig",] > 0.1 & add_info["TargetedNom",] <= 0.01) - 0.5, y = rep((3 - 1/4), times = length(which(add_info["TargetedSig",] > 0.1 & add_info["TargetedNom",] <= 0.01))), labels = ".")
      }
      if(length(which(is.na(add_info["TargetedSig",])))){
        rect(xleft = which(is.na(add_info["TargetedSig",])) - 1, ybottom = 2.5, xright = which(is.na(add_info["TargetedSig",])), ytop = 3, border = "white", col="grey")
      }
      
      if(length(which(add_info["CombinedExomeTargetedSig",] <= 0.0001) > 0)){
        rect(xleft = which(add_info["CombinedExomeTargetedSig",] <= 0.0001) - 1, ybottom = 1.5, xright = which(add_info["CombinedExomeTargetedSig",] <= 0.0001), ytop = 2, border = "white", col=pal_sig[5])
        text(x = rep(which(add_info["CombinedExomeTargetedSig",] <= 0.0001), each = 4) - rep(rep(c(2/3,1/3), each = 2), times = length(which(add_info["CombinedExomeTargetedSig",] < 0.0001))), y = rep(rep(c((2 - 1/6),(2 - 1/3)), times = 2), times = length(which(add_info["CombinedExomeTargetedSig",] < 0.0001))), labels = "*")
      }
      if(length(which(add_info["CombinedExomeTargetedSig",] <= 0.001 & add_info["CombinedExomeTargetedSig",] > 0.0001)) > 0){
        rect(xleft = which(add_info["CombinedExomeTargetedSig",] <= 0.001 & add_info["CombinedExomeTargetedSig",] > 0.0001) - 1, ybottom = 1.5, xright = which(add_info["CombinedExomeTargetedSig",] <= 0.001 & add_info["CombinedExomeTargetedSig",] > 0.0001), ytop = 2, border = "white", col=pal_sig[4])
        text(x = rep(which(add_info["CombinedExomeTargetedSig",] <= 0.001 & add_info["CombinedExomeTargetedSig",] > 0.0001), each = 3) - rep(c(1/2,2/3,1/3), times = length(which(add_info["CombinedExomeTargetedSig",] <= 0.001 & add_info["CombinedExomeTargetedSig",] > 0.0001))), y = rep(c((2 - 1/6),(2 - 1/3),(2 - 1/3)), times = length(which(add_info["CombinedExomeTargetedSig",] <= 0.001 & add_info["CombinedExomeTargetedSig",] > 0.0001))), labels = "*")
      }
      if(length(which(add_info["CombinedExomeTargetedSig",] <= 0.01 & add_info["CombinedExomeTargetedSig",] > 0.001)) > 0){
        rect(xleft = which(add_info["CombinedExomeTargetedSig",] <= 0.01 & add_info["CombinedExomeTargetedSig",] > 0.001) - 1, ybottom = 1.5, xright = which(add_info["CombinedExomeTargetedSig",] <= 0.01 & add_info["CombinedExomeTargetedSig",] > 0.001), ytop = 2, border = "white", col=pal_sig[3])
        text(x = rep(which(add_info["CombinedExomeTargetedSig",] <= 0.01 & add_info["CombinedExomeTargetedSig",] > 0.001), each = 2) - rep(c(2/3,1/3), times = length(which(add_info["CombinedExomeTargetedSig",] <= 0.01 & add_info["CombinedExomeTargetedSig",] > 0.001))), y = rep(c((2 - 1/4),(2 - 1/4)), times = length(which(add_info["CombinedExomeTargetedSig",] <= 0.01 & add_info["CombinedExomeTargetedSig",] > 0.001))), labels = "*")
      }
      if(length(which(add_info["CombinedExomeTargetedSig",] <= 0.1 & add_info["CombinedExomeTargetedSig",] > 0.01)) > 0){
        rect(xleft = which(add_info["CombinedExomeTargetedSig",] <= 0.1 & add_info["CombinedExomeTargetedSig",] > 0.01) - 1, ybottom = 1.5, xright = which(add_info["CombinedExomeTargetedSig",] <= 0.1 & add_info["CombinedExomeTargetedSig",] > 0.01), ytop = 2, border = "white", col=pal_sig[2])
        text(x = which(add_info["CombinedExomeTargetedSig",] <= 0.1 & add_info["CombinedExomeTargetedSig",] > 0.01) - 0.5, y = rep((2 - 1/4), times = length(which(add_info["CombinedExomeTargetedSig",] <= 0.1 & add_info["CombinedExomeTargetedSig",] > 0.01))), labels = "*")
      }
      if(length(which(add_info["CombinedExomeTargetedSig",] > 0.1 & add_info["CombinedExomeTargetedNom",] <= 0.01)) > 0){
        rect(xleft = which(add_info["CombinedExomeTargetedSig",] > 0.1 & add_info["CombinedExomeTargetedNom",] <= 0.01) - 1, ybottom = 1.5, xright = which(add_info["CombinedExomeTargetedSig",] > 0.1 & add_info["CombinedExomeTargetedNom",] <= 0.01), ytop = 2, border = "white", col=pal_sig[1])
        text(x = which(add_info["CombinedExomeTargetedSig",] > 0.1 & add_info["CombinedExomeTargetedNom",] <= 0.01) - 0.5, y = rep((2 - 1/4), times = length(which(add_info["CombinedExomeTargetedSig",] > 0.1 & add_info["CombinedExomeTargetedNom",] <= 0.01))), labels = ".")
      }
      if(length(which(is.na(add_info["CombinedExomeTargetedSig",])))){
        rect(xleft = which(is.na(add_info["CombinedExomeTargetedSig",])) - 1, ybottom = 1.5, xright = which(is.na(add_info["CombinedExomeTargetedSig",])), ytop = 2, border = "white", col="grey")
      }
      if(i == 1){
        mtext(text = "Exome dx", side = 2, at = 4.75, las = 2, cex = 0.5, adj = 1)
        mtext(text = "Exome dN/dS", side = 2, at = 4.25, las = 2, cex = 0.5, adj = 1)
        mtext(text = "Targeted dx", side = 2, at = 3.25, las = 2, cex = 0.5, adj = 1)
        mtext(text = "Targeted dN/dS", side = 2, at = 2.75, las = 2, cex = 0.5, adj = 1)
        mtext(text = "Combined dN/dS", side = 2, at = 1.75, las = 2, cex = 0.5, adj = 1)
        rasterImage(as.raster(t(matrix(pal_cov))), ybottom = 0.75, ytop = 1, xleft = 0, xright = ncol(add_info))
        lbsq <- seq(from = 0, to = ncol(add_info), length.out = 9)
        segments(x0 = lbsq, x1 = lbsq, y0 = rep(0.65, times = length(lbsq)), y1 = rep(0.75, times = length(lbsq)))
        text(x = lbsq, y = 0.5, labels = as.character(seq(from = 0, to = max_duplex_cov, length.out = length(lbsq))), cex = 0.75)
        text(x = lbsq[5], y = 0.1, labels = "Duplex coverage", cex = 1)
      }
    }
  }
  if (runman) { dev.copy(pdf, plotfilename, width=plotwidth, height=plotheight); dev.off() }
}

# Comparison non-synonymous mutant cell fraction plot
compplot = function(comp_list, genes2plot, driverdensity_sampleIDs, full_sampleIDs, muts, combined_targeted_only_dnds_sel_cv, targeted_dnds_sel_cv, gene2dc_array, max_duplex_cov, highlight_colour, pairwise_dnds = NULL, plotfilename = "Comparison_plot.pdf", plotwidth, plotheight, plottitles = plottitles, min_mcf_upper = 0, runman = T) {
  
  if (runman) { dev.new() }
  
  layout(matrix(c(1:(5*length(comp_list))), nrow = length(comp_list), ncol = 5, byrow = T), widths = c(5,0.75,1.6,0.75,5), heights = (lengths(genes2plot) + 5) / sum(lengths(genes2plot)))
  par(mar = c(1,0.4,6,0.4))
  
  for(i in 1:length(genes2plot)){
    maux = vector("list", 2)
    cellfraction = vector("list", 2)
    numsamples = vector("list", 2)
    nonsyndens = vector("list", 2)
    
    for(j in 1:2){
      maux[[j]] = muts[which(muts$sampleID %in% unlist(driverdensity_sampleIDs[[i]][j])),]
      cellfraction[[j]] = array(NA, dim = c(length(genes2plot[[i]]),8), dimnames = list(genes2plot[[i]], c("mis","non","spl","ind","mislow","nonlow","spllow","indlow")))
      numsamples[[j]] = length(unique(maux[[j]]$sampleID))
    }
    
    for (x in 1:length(genes2plot[[i]])) {
      for(j in 1:2){
        aux = maux[[j]][which(maux[[j]]$gene==genes2plot[[i]][x]),c("sampleID","impact","cellfraction","duplex_vaf")]
        per_sample_gene_cov <- gene2dc_array[which(gene2dc_array$gene == genes2plot[[i]][x]),which(colnames(gene2dc_array) %in% unlist(driverdensity_sampleIDs[[i]][j]))]
        if(nrow(aux) > 0){
          aux$sample_gene_cov <- NA
          aux$cov_weighting <- NA
          for(k in 1:length(per_sample_gene_cov)){
            aux$sample_gene_cov[which(aux$sampleID == names(per_sample_gene_cov)[k])] <- as.numeric(per_sample_gene_cov[k])
          }
          aux$cov_weighting <- aux$sample_gene_cov / sum(per_sample_gene_cov)
          
          aux$weighted_cellfraction <- aux$cellfraction * aux$cov_weighting
          aux$weighted_duplex_vaf <- aux$duplex_vaf * aux$cov_weighting
          
          cellfraction[[j]][x,c("mis","mislow")] = colSums(aux[aux$impact %in% c("Missense","Start_loss","Stop_loss"),c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
          cellfraction[[j]][x,c("non","nonlow")] = colSums(aux[aux$impact=="Nonsense",c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
          cellfraction[[j]][x,c("spl","spllow")] = colSums(aux[aux$impact=="Essential_Splice",c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
          cellfraction[[j]][x,c("ind","indlow")] = colSums(aux[aux$impact=="no-SNV",c("weighted_cellfraction","weighted_duplex_vaf")], na.rm = T)
        }else{
          cellfraction[[j]][x,c("mis","mislow")] = 0
          cellfraction[[j]][x,c("non","nonlow")] = 0
          cellfraction[[j]][x,c("spl","spllow")] = 0
          cellfraction[[j]][x,c("ind","indlow")] = 0
        }
      }
    }
    
    aux = as.data.frame(array(data = NA, dim = c(length(genes2plot[[i]]) * 2,12)))
    colnames(aux) <- c("gene","idx","cell_type","lower_bound","upper_bound","bound_diff","colour","combined_sig","combined_nom","combined_cov","pairwise_sig","pairwise_nom")
    aux$gene <- rep(genes2plot[[i]], times = 2)
    if(comp_list[[i]][[1]] != comp_list[[i]][[2]]){
      aux$cell_type <- rep(comp_list[[i]], each = length(genes2plot[[i]]))
    }else{
      aux$cell_type <- rep(paste(comp_list[[i]],c(1,2),sep = "_"), each = length(genes2plot[[i]]))
      comp_list[[i]][[1]] <- paste(comp_list[[i]][[1]],1, sep = "_")
      comp_list[[i]][[2]] <- paste(comp_list[[i]][[2]],2, sep = "_")
    }
    
    for(j in 1:2){
      nonsyndens[[j]] = cellfraction[[j]] * 100
      aux$upper_bound[(((j - 1) * length(genes2plot[[i]])) + 1):(j * length(genes2plot[[i]]))] <- rowSums(nonsyndens[[j]][,1:4])
      aux$lower_bound[(((j - 1) * length(genes2plot[[i]])) + 1):(j * length(genes2plot[[i]]))] <- rowSums(nonsyndens[[j]][,5:8])
    }
    
    aux$bound_diff <- aux$upper_bound - aux$lower_bound
    
    for(j in 1:length(genes2plot[[i]])){
      if(aux$upper_bound[j] > aux$upper_bound[j + length(genes2plot[[i]])]){
        aux$colour[j] <- "skyblue3"
        aux$colour[j +  length(genes2plot[[i]])] <- "lightskyblue1" 
      }else if(aux$upper_bound[j] < aux$upper_bound[j + length(genes2plot[[i]])]){
        aux$colour[j] <- "lightskyblue1"
        aux$colour[j +  length(genes2plot[[i]])] <- "skyblue3" 
      }
    }
    
    for(x in 1:length(genes2plot[[i]])){
      for(j in 1:2){
        aux$combined_cov[x + (j - 1) * (length(genes2plot[[i]]))] <- sum(gene2dc_array[which(gene2dc_array$gene == genes2plot[[i]][x]),which(colnames(gene2dc_array) %in% unlist(full_sampleIDs[[i]][[j]]))])
        if(!(is.null(combined_targeted_only_dnds_sel_cv[[i]][[j]]))){
          if(length(which(colnames(combined_targeted_only_dnds_sel_cv[[i]][[j]]) == "qglobalpos_m")) == 1){
            aux$combined_sig[x + (j - 1) * (length(genes2plot[[i]]))] <- combined_targeted_only_dnds_sel_cv[[i]][[j]]$qglobalpos_m[which(combined_targeted_only_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
            aux$combined_nom[x + (j - 1) * (length(genes2plot[[i]]))] <- combined_targeted_only_dnds_sel_cv[[i]][[j]]$pglobalpos_m[which(combined_targeted_only_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
          }else{
            aux$combined_sig[x + (j - 1) * (length(genes2plot[[i]]))] <- combined_targeted_only_dnds_sel_cv[[i]][[j]]$qglobalpos_cv[which(combined_targeted_only_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
            aux$combined_nom[x + (j - 1) * (length(genes2plot[[i]]))] <- combined_targeted_only_dnds_sel_cv[[i]][[j]]$pglobalpos_cv[which(combined_targeted_only_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
          }
        }else{
          if(length(which(colnames(targeted_dnds_sel_cv[[i]][[j]]) == "qglobalpos_m")) == 1){
            aux$combined_sig[x + (j - 1) * (length(genes2plot[[i]]))] <- targeted_dnds_sel_cv[[i]][[j]]$qglobalpos_m[which(targeted_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
            aux$combined_nom[x + (j - 1) * (length(genes2plot[[i]]))] <- targeted_dnds_sel_cv[[i]][[j]]$pglobalpos_m[which(targeted_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
          }else{
            aux$combined_sig[x + (j - 1) * (length(genes2plot[[i]]))] <- targeted_dnds_sel_cv[[i]][[j]]$qglobalpos_cv[which(targeted_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
            aux$combined_nom[x + (j - 1) * (length(genes2plot[[i]]))] <- targeted_dnds_sel_cv[[i]][[j]]$pglobalpos_cv[which(targeted_dnds_sel_cv[[i]][[j]]$gene_name == genes2plot[[i]][x])]
          }
        }
      }
    }
    
    if(!(is.null(pairwise_dnds))){
      for(x in 1:length(genes2plot[[i]])){
        for(j in 1:2){
          aux$pairwise_sig[x + (j - 1) * (length(genes2plot[[i]]))] <- pairwise_dnds[[i]][[j]]$qval[which(pairwise_dnds[[i]][[j]]$gene_name == genes2plot[[i]][x])]
          aux$pairwise_nom[x + (j - 1) * (length(genes2plot[[i]]))] <- pairwise_dnds[[i]][[j]]$pval[which(pairwise_dnds[[i]][[j]]$gene_name == genes2plot[[i]][x])]
        }
      }
    }
    
    aux$idx[order(aux$upper_bound[1:length(genes2plot[[i]])], decreasing = T)] <- 1:length(genes2plot[[i]])
    aux$idx[order(aux$upper_bound[(1:length(genes2plot[[i]]))], decreasing = T) + length(genes2plot[[i]])] <- 1:length(genes2plot[[i]])
    aux <- aux[order(aux$idx, decreasing = T),]
    
    plot_breaks <- pretty_breaks()(0:max(round(max(aux$upper_bound * 1.1)),min_mcf_upper))
    
    pos = barplot(height = -(aux[which(aux$cell_type == comp_list[[i]][1]),"bound_diff"]), offset = -(aux[which(aux$cell_type == comp_list[[i]][1]),"lower_bound"]), las = 2, horiz = T, col = aux[which(aux$cell_type == comp_list[[i]][1]),"colour"], border = aux[which(aux$cell_type == comp_list[[i]][1]),"colour"], ylab = "", xlim = c(-(max(round(max(aux$upper_bound * 1.1)),min_mcf_upper)),0.02), axes = F, axisnames = F)
    title(main = as.character(plottitles[[i]][1]), line = 2.5)
    axis(side = 3, at = -(rev(plot_breaks)), labels = rev(plot_breaks))
    abline(v = -(rev(setdiff(plot_breaks,0))), col = "gray90")
    abline(v = 0, col = "black")
    
    pal_cov = mako(n = 100, direction = -1)
    pal_sig <- brewer.pal(n = 5, name = "YlOrRd")
    
    pos_diff <- pos[2] - pos[1]
    
    xlim = max(round(max(aux$upper_bound * 1.1)),min_mcf_upper)
    
    if(!(is.null(pairwise_dnds))){
      if(length(which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])) > 0){
        rect(xleft = -(aux$upper_bound[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), xright = -(aux$upper_bound[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/90), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[5])
        text(y = rep(pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), each = 4) + rep(rep(c(2/3,1/3) * pos_diff, each = 2), times = length(which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1]))), x = rep(-(aux$upper_bound[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), each = 4) + rep(rep(c(2/3,1/3) * (xlim/18), times = 2), times = length(which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
      }
      if(length(which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])) > 0){
        rect(xleft = -(aux$upper_bound[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), xright = -(aux$upper_bound[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/90), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[4])
        text(y = rep(pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), each = 3) + rep(c(1/3,2/3,1/3) * pos_diff, times = length(which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1]))), x = rep(-(aux$upper_bound[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), each = 3) + rep(c(1/3,1/2,2/3) * (xlim/18), times = length(which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
      }
      if(length(which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])) > 0){
        rect(xleft = -(aux$upper_bound[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), xright = -(aux$upper_bound[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/90), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[3])
        text(y = rep(pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), each = 2) + rep(c(1/2,1/2) * pos_diff, times = length(which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1]))), x = rep(-(aux$upper_bound[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), each = 2) + rep(c(1/3,2/3) * (xlim/18), times = length(which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
      }
      if(length(which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])) > 0){
        rect(xleft = -(aux$upper_bound[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), xright = -(aux$upper_bound[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])]) - (xlim/90), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[2])
        text(y = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)) + rep(c(1/2) * pos_diff, times = length(which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1]))), x = -(aux$upper_bound[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15) + rep(1/2 * (xlim/18), times = length(which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
      }
      if(length(which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])) > 0){
        rect(xleft = -(aux$upper_bound[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15), xright = -(aux$upper_bound[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])]) - (xlim/90), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[1])
        text(y = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])] + 1] - ((pos_diff/2) - (pos_diff/20)) + rep(c(1/2) * pos_diff, times = length(which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1]))), x = -(aux$upper_bound[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1])]) - (xlim/15) + rep(1/2 * (xlim/18), times = length(which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][1]))), labels = ".")
      }
    }
    
    plot(x = c(0, 2), y = c(0, length(genes2plot[[i]])), type= "n", xlab = "", ylab = "", axes = F)
    
    for(x in 1:length(genes2plot[[i]])){
      color <- pal_cov[round((aux$combined_cov[which(aux$cell_type == comp_list[[i]][1] & aux$gene == genes2plot[[i]][x])] / max_duplex_cov) * 100)]
      rect(xleft = 1, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$gene == genes2plot[[i]][x] & aux$cell_type == comp_list[[i]][1])], xright = 2, ytop = length(genes2plot[[i]]) - aux$idx[which(aux$gene == genes2plot[[i]][x] & aux$cell_type == comp_list[[i]][1])] + 1, border = "white", col=color)
    }
    
    if(length(which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])) > 0){
      rect(xleft = 0, xright = 1, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1, border = "white", col=pal_sig[5])
      text(y = rep(length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1])], each = 4) + rep(rep(c(2/3,1/3), each = 2), times = length(which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1]))), x = rep(rep(c(1/3,2/3),times = 2), times = length(which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
    }
    
    if(length(which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])) > 0){
      rect(xleft = 0, xright = 1, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])] + 1, border = "white", col=pal_sig[4])
      text(y = rep(length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][1])], each = 3) + rep(c(1/3,2/3,1/3), times = length(which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][1]))), x = rep(c(1/3,1/2,2/3), times = length(which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
    }
    
    if(length(which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][1])) > 0){
      rect(xleft = 0, xright = 1, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][1])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][1])] + 1, border = "white", col=pal_sig[3])
      text(y = rep(length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][1])], each = 2) + rep(c(1/2,1/2), times = length(which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][1]))), x = rep(c(1/3,2/3), times = length(which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
    }
    
    if(length(which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][1])) > 0){
      rect(xleft = 0, xright = 1, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][1])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][1])] + 1, border = "white", col=pal_sig[2])
      text(y = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][1])] + rep(c(1/2), times = length(which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][1]))), x = rep(c(1/2), times = length(which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][1]))), labels = "*")
    }
    
    if(length(which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][1])) > 0){
      rect(xleft = 0, xright = 1, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][1])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][1])] + 1, border = "white", col=pal_sig[1])
      text(y = length(genes2plot[[i]]) - aux$idx[which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][1])] + rep(c(1/2), times = length(which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][1]))), x = rep(c(1/2), times = length(which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][1]))), labels = ".")
    }
    
    gene_names <- lapply(aux$gene[which(aux$cell_type == comp_list[[i]][1])], function(x) bquote(italic(.(x))))
    plot(x = c(-3, 3), y = c(0, length(genes2plot[[i]])), type= "n", xlab = "", ylab = "", axes = F)
    text(y = rev(c(aux$idx[which(aux$cell_type == comp_list[[i]][1])])) - 0.5, x = 0, labels = as.expression(gene_names), cex = 1.5, col = highlight_colour[[i]][aux$gene[which(aux$cell_type == comp_list[[i]][1])]])
    
    plot(x = c(0, 2), y = c(0, length(genes2plot[[i]])), type= "n", xlab = "", ylab = "", axes = F)
    
    for(x in 1:length(genes2plot[[i]])){
      color <- pal_cov[round((aux$combined_cov[which(aux$cell_type == comp_list[[i]][2] & aux$gene == genes2plot[[i]][x])] / max_duplex_cov) * 100)]
      rect(xleft = 0, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$gene == genes2plot[[i]][x] & aux$cell_type == comp_list[[i]][2])], xright = 1, ytop = length(genes2plot[[i]]) - aux$idx[which(aux$gene == genes2plot[[i]][x] & aux$cell_type == comp_list[[i]][2])] + 1, border = "white", col=color)
    }
    
    if(length(which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])) > 0){
      rect(xleft = 1, xright = 2, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1, border = "white", col=pal_sig[5])
      text(y = rep(length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])], each = 4) + rep(rep(c(2/3,1/3), each = 2), times = length(which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2]))), x = rep(rep(c(4/3,5/3),times = 2), times = length(which(aux$combined_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
    }
    
    if(length(which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])) > 0){
      rect(xleft = 1, xright = 2, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1, border = "white", col=pal_sig[4])
      text(y = rep(length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])], each = 3) + rep(c(1/3,2/3,1/3), times = length(which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][2]))), x = rep(c(4/3,3/2,5/3), times = length(which(aux$combined_sig <= 0.001 & aux$combined_sig > 0.0001 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
    }
    
    if(length(which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][2])) > 0){
      rect(xleft = 1, xright = 2, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][2])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][2])] + 1, border = "white", col=pal_sig[3])
      text(y = rep(length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][2])], each = 2) + rep(c(1/2,1/2), times = length(which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][2]))), x = rep(c(4/3,5/3), times = length(which(aux$combined_sig <= 0.01 & aux$combined_sig > 0.001 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
    }
    
    if(length(which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][2])) > 0){
      rect(xleft = 1, xright = 2, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][2])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][2])] + 1, border = "white", col=pal_sig[2])
      text(y = length(genes2plot[[i]]) - aux$idx[which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][2])] + rep(c(1/2), times = length(which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][2]))), x = rep(c(3/2), times = length(which(aux$combined_sig <= 0.1 & aux$combined_sig > 0.01 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
    }
    
    if(length(which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][2])) > 0){
      rect(xleft = 1, xright = 2, ybottom = length(genes2plot[[i]]) - aux$idx[which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][2])], ytop = length(genes2plot[[i]]) - aux$idx[which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][2])] + 1, border = "white", col=pal_sig[1])
      text(y = length(genes2plot[[i]]) - aux$idx[which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][2])] + rep(c(1/2), times = length(which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][2]))), x = rep(c(3/2), times = length(which(aux$combined_nom <= 0.01 & aux$combined_sig > 0.1 & aux$cell_type == comp_list[[i]][2]))), labels = ".")
    }
    
    barplot(height = aux[which(aux$cell_type == comp_list[[i]][2]),c("bound_diff")], offset = aux[which(aux$cell_type == comp_list[[i]][2]),"lower_bound"], las = 2, horiz = T, col=aux[which(aux$cell_type == comp_list[[i]][2]),"colour"], border = NA, ylab = "", xlim = c(-0.02,max(round(max(aux$upper_bound * 1.1)),min_mcf_upper)), axes = F, axisnames = F)
    title(as.character(plottitles[[i]][2]), line = 2.5)
    axis(side = 3, at = plot_breaks, labels = plot_breaks)
    abline(v = setdiff(plot_breaks,0), col = "gray90")
    abline(v = 0, col = "black")
    
    xlim = max(round(max(aux$upper_bound * 1.1)),min_mcf_upper)
    
    if(!(is.null(pairwise_dnds))){
      if(length(which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])) > 0){
        rect(xleft = aux$upper_bound[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])] + (xlim/90), xright = +(aux$upper_bound[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/15), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[5])
        text(y = rep(pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), each = 4) + rep(rep(c(2/3,1/3) * pos_diff, each = 2), times = length(which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2]))), x = rep(+(aux$upper_bound[which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), each = 4) + rep(rep(c(2/3,1/3) * (xlim/18), times = 2), times = length(which(aux$pairwise_sig <= 0.0001 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
      }
      if(length(which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])) > 0){
        rect(xleft = +(aux$upper_bound[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), xright = +(aux$upper_bound[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/15), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[4])
        text(y = rep(pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), each = 3) + rep(c(1/3,2/3,1/3) * pos_diff, times = length(which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2]))), x = rep(+(aux$upper_bound[which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), each = 3) + rep(c(1/3,1/2,2/3) * (xlim/18), times = length(which(aux$pairwise_sig <= 0.001 & aux$pairwise_sig > 0.0001 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
      }
      if(length(which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])) > 0){
        rect(xleft = +(aux$upper_bound[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), xright = +(aux$upper_bound[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/15), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[3])
        text(y = rep(pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), each = 2) + rep(c(1/2,1/2) * pos_diff, times = length(which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2]))), x = rep(+(aux$upper_bound[which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), each = 2) + rep(c(1/3,2/3) * (xlim/18), times = length(which(aux$pairwise_sig <= 0.01 & aux$pairwise_sig > 0.001 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
      }
      if(length(which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])) > 0){
        rect(xleft = +(aux$upper_bound[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), xright = +(aux$upper_bound[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])]) + (xlim/15), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[2])
        text(y = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)) + rep(c(1/2) * pos_diff, times = length(which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2]))), x = +(aux$upper_bound[which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90) + rep(1/2 * (xlim/18), times = length(which(aux$pairwise_sig <= 0.1 & aux$pairwise_sig > 0.01 & aux$cell_type == comp_list[[i]][2]))), labels = "*")
      }
      if(length(which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])) > 0){
        rect(xleft = +(aux$upper_bound[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90), xright = +(aux$upper_bound[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])]) + (xlim/15), ybottom = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)), ytop = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])] + 1] + ((pos_diff/2) - (pos_diff/20)), border = "white", col=pal_sig[1])
        text(y = pos[length(genes2plot[[i]]) - aux$idx[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])] + 1] - ((pos_diff/2) - (pos_diff/20)) + rep(c(1/2) * pos_diff, times = length(which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2]))), x = +(aux$upper_bound[which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2])]) + (xlim/90) + rep(1/2 * (xlim/18), times = length(which(aux$pairwise_nom <= 0.01 & aux$pairwise_sig > 0.1 & aux$cell_type == comp_list[[i]][2]))), labels = ".")
      }
    }
  }
  if (runman) { dev.copy(pdf, plotfilename, width=plotwidth, height=plotheight); dev.off() }
}

# Mutation distribution plot across the entire gene (i.e. including introns and UTRs)
full_gene_mutation_plot = function(gene, mutations, sample_list, title, gene_coverage_path, combined_per_sample_cov, immune_genes, plotfilename, runman = T) {
  muts <- mutations[which(mutations$sampleID %in% sample_list),]
  muts <- muts[which(!(duplicated(paste(muts$donor, muts$chr, muts$pos, muts$ref, muts$mut, sep = "_")))),]
  
  gene_coverage <- read.table(gene_coverage_path, header = T, stringsAsFactors = F)
  gene_coverage <- gene_coverage[,c(1:4,which(colnames(gene_coverage) %in% sample_list)),]
  
  data(list = "refcds_hg19", package = "dndscv")
  
  ref_cds_gene_list <- NULL
  
  for(i in 1:length(RefCDS)){
    ref_cds_gene_list <- c(ref_cds_gene_list,RefCDS[[i]]$gene_name)
  }
  
  gene_coverage$idx <- 1:nrow(gene_coverage)
  gene_coverage$total <- NA
  for(i in 1:nrow(gene_coverage)){
    gene_coverage$total[i] <- sum(gene_coverage[i,which(substr(colnames(gene_coverage),1,2) %in% c("PD","EL"))])
  }
  
  gene_coverage$coding_transitions <- NA
  for(i in 1:nrow(gene_coverage)){
    if(i == 1){
      gene_coverage$coding_transitions[i] <- 1
    }else if(gene_coverage$coding[i] == gene_coverage$coding[i - 1]){
      gene_coverage$coding_transitions[i] <- gene_coverage$coding_transitions[i - 1]
    }else{
      gene_coverage$coding_transitions[i] <- gene_coverage$coding_transitions[i - 1] + 1
    }
  }
  
  gene_coding_transition_table <- as.data.frame(array(data = NA, dim = c(max(gene_coverage$coding_transitions),6)))
  colnames(gene_coding_transition_table) <- c("exonic","start_pos","end_pos","start_idx","end_idx","max_coverage")
  
  for(i in 1:nrow(gene_coding_transition_table)){
    gene_coding_transition_table$exonic[i] <- unique(gene_coverage$coding[which(gene_coverage$coding_transitions == i)])
    gene_coding_transition_table$start_pos[i] <- gene_coverage$pos[min(which(gene_coverage$coding_transitions == i))]
    gene_coding_transition_table$end_pos[i] <- gene_coverage$pos[max(which(gene_coverage$coding_transitions == i))]
    gene_coding_transition_table$start_idx[i] <- gene_coverage$idx[min(which(gene_coverage$coding_transitions == i))]
    gene_coding_transition_table$end_idx[i] <- gene_coverage$idx[max(which(gene_coverage$coding_transitions == i))]
    gene_coding_transition_table$max_coverage[i] <- max(gene_coverage$total[which(gene_coverage$coding_transitions == i)])
  }
  
  gene_coverage$pos_transitions <- NA
  for(i in 1:nrow(gene_coverage)){
    if(i == 1){
      gene_coverage$pos_transitions[i] <- 1
    }else if(gene_coverage$pos[i] == gene_coverage$pos[i - 1] + 1){
      gene_coverage$pos_transitions[i] <- gene_coverage$pos_transitions[i - 1]
    }else{
      gene_coverage$pos_transitions[i] <- gene_coverage$pos_transitions[i - 1] + 1
    }
  }
  
  gene_pos_transition_table <- as.data.frame(array(data = NA, dim = c(max(gene_coverage$pos_transitions),4)))
  colnames(gene_pos_transition_table) <- c("start_pos","end_pos","start_idx","end_idx")
  
  for(i in 1:nrow(gene_pos_transition_table)){
    gene_pos_transition_table$start_pos[i] <- gene_coverage$pos[min(which(gene_coverage$pos_transitions == i))]
    gene_pos_transition_table$end_pos[i] <- gene_coverage$pos[max(which(gene_coverage$pos_transitions == i))]
    gene_pos_transition_table$start_idx[i] <- gene_coverage$idx[min(which(gene_coverage$pos_transitions == i))]
    gene_pos_transition_table$end_idx[i] <- gene_coverage$idx[max(which(gene_coverage$pos_transitions == i))]
  }
  
  gene_muts <- muts[which(muts$chr == immune_genes$chr[which(immune_genes$gene == gene)] & muts$pos >= gene_coverage$pos[1] & muts$pos <= gene_coverage$pos[nrow(gene_coverage)] & (muts$gene == gene | is.na(muts$gene))),]
  
  if(nrow(gene_muts) > 0){
    gene_muts$mut_type <- NA
    
    for(i in 1:nrow(gene_muts)){
      if(!(is.na(gene_muts$gene[i]))){
        if(gene_muts$gene[i] == gene){
          gene_muts$mut_type[i] <- gene_muts$impact[i]
        }
      }else if(gene_muts$type[i] == "snv"){
        gene_muts$mut_type[i] <- "non_coding_snv"
      }else{
        gene_muts$mut_type[i] <- "non_coding_other"
      }
    }
    
    gene_start_loss <- which(gene_muts$mut_type != "Synonymous" & substr(gene_muts$aachange,2,nchar(gene_muts$aachange)-1) == "1")
    if(length(gene_start_loss) >= 1){
      gene_muts$mut_type[gene_start_loss] <- "Start_loss"
    }
    
    mut_types <- names(table(gene_muts$mut_type))
    bin_size <- 50
    
    binned_muts <- as.data.frame(array(data = NA, dim = c(ceiling(nrow(gene_coverage) / bin_size) * length(mut_types),3)))
    colnames(binned_muts) <- c("idx","mut_type","count")
    
    binned_muts$idx <- rep(seq(from = (bin_size)/2, to = (ceiling(nrow(gene_coverage) / bin_size) * bin_size) - (bin_size)/2, by = bin_size), each = length(mut_types))
    binned_muts$mut_type <- rep(mut_types, times = ceiling(nrow(gene_coverage) / bin_size))
    
    for(i in 1:nrow(binned_muts)){
      binned_muts$count[i] <- length(which(gene_muts$mut_type == binned_muts$mut_type[i] & gene_muts$pos >= gene_coverage$pos[binned_muts$idx[i] - ((bin_size)/2 - 1)] & gene_muts$pos <= gene_coverage$pos[min(binned_muts$idx[i] + (bin_size)/2,nrow(gene_coverage))]))
    }
    
    gene_idx <- which(ref_cds_gene_list == gene)
    
    if(RefCDS[[gene_idx]]$strand == "-1"){
      strand = "-"
    }else if(RefCDS[[gene_idx]]$strand == "1"){
      strand = "+"
    }
    
    gene_coverage_smooth <- as.data.frame(array(data = NA, dim = c(nrow(gene_coverage),9)))
    colnames(gene_coverage_smooth) <- c("chr","pos","idx","total","pos_transitions","idx_first","idx_last","bin_size","bin_average")
    gene_coverage_smooth$chr <- gene_coverage$chr
    gene_coverage_smooth$pos <- gene_coverage$pos
    gene_coverage_smooth$idx <- gene_coverage$idx
    gene_coverage_smooth$total <- gene_coverage$total
    gene_coverage_smooth$pos_transitions <- gene_coverage$pos_transitions
    
    for(i in 1:nrow(gene_coverage_smooth)){
      gene_coverage_smooth$idx_first[i] <- max(min(gene_coverage_smooth$idx[which(gene_coverage_smooth$pos_transitions == gene_coverage_smooth$pos_transitions[i])]),gene_coverage_smooth$idx[i] - (bin_size/2))
      gene_coverage_smooth$idx_last[i] <- min(max(gene_coverage_smooth$idx[which(gene_coverage_smooth$pos_transitions == gene_coverage_smooth$pos_transitions[i])]),gene_coverage_smooth$idx[i] + (bin_size/2))
      gene_coverage_smooth$bin_size[i] <- (gene_coverage_smooth$idx_last[i] - gene_coverage_smooth$idx_first[i]) + 1
      gene_coverage_smooth$bin_average[i] <- sum(gene_coverage_smooth$total[gene_coverage_smooth$idx_first[i]:gene_coverage_smooth$idx_last[i]]) / gene_coverage_smooth$bin_size[i]
    }
    
    gene_region_burden <- gene_coding_transition_table
    gene_region_burden$chr <- unique(gene_muts$chr)
    gene_region_burden$gene <- gene
    intron_rows <- which(gene_region_burden$exonic == 0)
    intron_rows <- intron_rows[which(!(intron_rows %in% c(1,nrow(gene_region_burden))))]
    
    gene_region_burden <- gene_region_burden[sort(c(1:nrow(gene_region_burden),intron_rows)),]
    
    for(i in 1:(nrow(gene_region_burden) - 1)){
      if(gene_region_burden$start_pos[i] == gene_region_burden$start_pos[i + 1]){
        if((gene_region_burden$end_idx[i] - gene_region_burden$start_idx[i]) + 1 == 1000){
          gene_region_burden$end_pos[i] = gene_region_burden$start_pos[i] + 499
          gene_region_burden$start_pos[i + 1] = gene_region_burden$end_pos[i + 1] - 499
        }else{
          gene_region_burden$end_pos[i] = floor((gene_region_burden$start_pos[i] + gene_region_burden$end_pos[i])/2)
          gene_region_burden$start_pos[i+1] = gene_region_burden$end_pos[i] + 1
        }
      }
    }
    
    gene_region_burden <- gene_region_burden[,c("gene","chr","start_pos","end_pos","exonic")]
    gene_region_burden$total_cov <- NA
    gene_region_burden$mut_count <- NA
    for(i in 1:nrow(gene_region_burden)){
      gene_region_burden$total_cov[i] <- sum(gene_coverage$total[which(gene_coverage$pos >= gene_region_burden$start_pos[i] & gene_coverage$pos <= gene_region_burden$end_pos[i])])
      gene_region_burden$mut_count[i] <- length(which(gene_muts$pos >= gene_region_burden$start_pos[i] & gene_muts$pos <= gene_region_burden$end_pos[i]))
    }
    gene_region_burden$burden <- gene_region_burden$mut_count / gene_region_burden$total_cov
    
    
    if(strand == "-"){
      gene_region_burden <- gene_region_burden[c(nrow(gene_region_burden):1),]
    }
    
    gene_region_burden$region_name <- NA
    utr_regions <- which(gene_region_burden$exonic == 2)
    coding_exon_regions <- which(gene_region_burden$exonic == 1)
    downstream_intron_regions <- coding_exon_regions + 1
    downstream_intron_regions <- downstream_intron_regions[which(!(downstream_intron_regions %in% utr_regions) & !(downstream_intron_regions == nrow(gene_region_burden)))]
    upstream_intron_regions <- coding_exon_regions - 1
    upstream_intron_regions <- upstream_intron_regions[which(!(upstream_intron_regions %in% utr_regions) & !(upstream_intron_regions == 1))]
    
    gene_region_burden$region_name[1] <- paste0(gene,"_upstream_gene_body")
    gene_region_burden$region_name[utr_regions[which(utr_regions < min(coding_exon_regions))]] <- paste0(gene,"_UTR_5prime")
    gene_region_burden$region_name[utr_regions[which(utr_regions > max(coding_exon_regions))]] <- paste0(gene,"_UTR_3prime")
    gene_region_burden$region_name[coding_exon_regions] <- paste0(gene,"_coding_exon_",1:length(coding_exon_regions))
    for(i in 1:length(downstream_intron_regions)){
      gene_region_burden$region_name[downstream_intron_regions[i]] <- paste0(gene_region_burden$region_name[downstream_intron_regions[i] - 1],"_downstream_intronic")
    }
    for(i in 1:length(upstream_intron_regions)){
      gene_region_burden$region_name[upstream_intron_regions[i]] <- paste0(gene_region_burden$region_name[upstream_intron_regions[i] + 1],"_upstream_intronic")
    }
    gene_region_burden$region_name[nrow(gene_region_burden)] <- paste0(gene,"_downstream_gene_body")
    gene_region_burden$region_name[which(is.na(gene_region_burden$region_name))] <- paste0(gene,"_UTR_intronic")
    
    suppressWarnings({
      mut_histo_plot <- ggplot() +
        geom_rect(data = gene_coding_transition_table[which(gene_coding_transition_table$exonic == 1),], aes(xmin = start_idx, xmax = end_idx, ymin = 0, ymax = ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1)), fill = "firebrick", alpha = 0.2) +
        geom_rect(data = gene_coding_transition_table[which(gene_coding_transition_table$exonic == 0),], aes(xmin = start_idx, xmax = end_idx, ymin = 0, ymax = ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1)), fill = "steelblue", alpha = 0.2) +
        geom_rect(data = gene_coding_transition_table[which(gene_coding_transition_table$exonic == 2),], aes(xmin = start_idx, xmax = end_idx, ymin = 0, ymax = ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1)), fill = "chartreuse4", alpha = 0.2) +
        geom_bar(data = binned_muts, aes(x = idx, y = count, fill = mut_type), stat = "identity", width = bin_size) +
        geom_line(data = gene_coverage_smooth, aes(x = idx, y = bin_average * (ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1) / (ceiling((max(total))/10000) * 10000))), alpha = 0.2) +
        scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding_snv" = "#B2DF8A","non_coding_other" = "darkkhaki","Stop_loss" = "red", "Start_loss" = "firebrick")) +
        geom_vline(xintercept = gene_pos_transition_table$start_idx[c(2:nrow(gene_pos_transition_table))]) +
        theme_bw() +
        scale_x_continuous(expand = c(0,0), breaks = gene_pos_transition_table$start_idx, labels = gene_pos_transition_table$start_pos) +
        labs(x = paste0("Chromosome ",gene_coverage$chr[1]), fill = "Impact", y = "Mutation count", title = paste0(gene," - ",title," samples (",round(sum(combined_per_sample_cov[which(combined_per_sample_cov$gene == gene),which(colnames(combined_per_sample_cov) %in% sample_list)])), " mean coding dx)\n",RefCDS[[gene_idx]]$protein_id," (",strand," strand): ",sum(binned_muts$count)," mutations")) +
        theme(plot.title = element_text(hjust = 0.5, size = 15)) +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.text.x = element_text(hjust = 0)) +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(), expand = c(0,0), limits = c(0,ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1)), sec.axis = sec_axis(~ . * ((ceiling(max(gene_coverage_smooth$total)/10000) * 10000)/(ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1))), name = "Cumulative duplex coverage")) +
        theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))
      print(mut_histo_plot)
      if (runman) { ggsave(plot = mut_histo_plot, filename = plotfilename, width = 20, height = 5)}
    })
  }
}

# Coding mutation distribution plot with comparison to COSMIC mutations (annotated by dNdS with same RefCDS)
coding_mutation_plot = function(gene, mutations, sample_list, title, gene_coverage_path, cosmic_dir_path, interpro_domain_dir_path, immune_genes, plotfilename, runman = T){
  muts <- mutations[which(mutations$sampleID %in% sample_list),]
  muts <- muts[which(!(duplicated(paste(muts$donor, muts$chr, muts$pos, muts$ref, muts$mut, sep = "_")))),]
  
  data(list = "refcds_hg19", package = "dndscv")
  
  gene_num <- which(immune_genes$gene == gene)
  
  ref_cds_gene_list <- NULL
  
  for(i in 1:length(RefCDS)){
    ref_cds_gene_list <- c(ref_cds_gene_list,RefCDS[[i]]$gene_name)
  }
  
  gene_coverage <- read.table(gene_coverage_path, header = T, stringsAsFactors = F)
  gene_coverage <- gene_coverage[,c(1:4,which(colnames(gene_coverage) %in% sample_list)),]
  gene_coverage <- gene_coverage[which(gene_coverage$coding == 1),]
  
  if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
    gene_coverage$idx <- 1:nrow(gene_coverage)
  }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
    gene_coverage$idx <- nrow(gene_coverage):1
  }
  
  gene_coverage <- gene_coverage[order(gene_coverage$idx, decreasing = F),]
  gene_coverage$total <- NA
  for(i in 1:nrow(gene_coverage)){
    gene_coverage$total[i] <- sum(gene_coverage[i,which(substr(colnames(gene_coverage),1,2) %in% c("PD","EL"))])
  }
  
  gene_coverage$pos_transitions <- NA
  for(i in 1:nrow(gene_coverage)){
    if(i == 1){
      gene_coverage$pos_transitions[i] <- 1
    }else if(abs(gene_coverage$pos[i] - gene_coverage$pos[i - 1]) == 1){
      gene_coverage$pos_transitions[i] <- gene_coverage$pos_transitions[i - 1]
    }else{
      gene_coverage$pos_transitions[i] <- gene_coverage$pos_transitions[i - 1] + 1
    }
  }
  
  gene_pos_transition_table <- as.data.frame(array(data = NA, dim = c(max(gene_coverage$pos_transitions),4)))
  colnames(gene_pos_transition_table) <- c("start_pos","end_pos","start_idx","end_idx")
  
  for(i in 1:nrow(gene_pos_transition_table)){
    gene_pos_transition_table$start_pos[i] <- gene_coverage$pos[min(which(gene_coverage$pos_transitions == i))]
    gene_pos_transition_table$end_pos[i] <- gene_coverage$pos[max(which(gene_coverage$pos_transitions == i))]
    gene_pos_transition_table$start_idx[i] <- gene_coverage$idx[min(which(gene_coverage$pos_transitions == i))]
    gene_pos_transition_table$end_idx[i] <- gene_coverage$idx[max(which(gene_coverage$pos_transitions == i))]
  }
  
  gene_coverage_smooth <- as.data.frame(array(data = NA, dim = c(nrow(gene_coverage),9)))
  colnames(gene_coverage_smooth) <- c("chr","pos","idx","total","pos_transitions","idx_first","idx_last","bin_size","bin_average")
  gene_coverage_smooth$chr <- gene_coverage$chr
  gene_coverage_smooth$pos <- gene_coverage$pos
  gene_coverage_smooth$idx <- gene_coverage$idx
  gene_coverage_smooth$total <- gene_coverage$total
  gene_coverage_smooth$pos_transitions <- gene_coverage$pos_transitions
  
  gene_cov_bin_size <- 50
  
  for(i in 1:nrow(gene_coverage_smooth)){
    gene_coverage_smooth$idx_first[i] <- max(min(gene_coverage_smooth$idx[which(gene_coverage_smooth$pos_transitions == gene_coverage_smooth$pos_transitions[i])]),gene_coverage_smooth$idx[i] - (gene_cov_bin_size/2))
    gene_coverage_smooth$idx_last[i] <- min(max(gene_coverage_smooth$idx[which(gene_coverage_smooth$pos_transitions == gene_coverage_smooth$pos_transitions[i])]),gene_coverage_smooth$idx[i] + (gene_cov_bin_size/2))
    gene_coverage_smooth$gene_cov_bin_size[i] <- (gene_coverage_smooth$idx_last[i] - gene_coverage_smooth$idx_first[i]) + 1
    gene_coverage_smooth$bin_average[i] <- sum(gene_coverage_smooth$total[gene_coverage_smooth$idx_first[i]:gene_coverage_smooth$idx_last[i]]) / gene_coverage_smooth$gene_cov_bin_size[i]
  }
  
  gene_muts <- muts[which(muts$gene == gene),]
  
  if(nrow(gene_muts) > 0){
    
    for(i in 1:nrow(gene_muts)){
      if(!(gene_muts$pos[i] %in% gene_coverage$pos)){
        diff_in_pos <- abs(gene_coverage$pos - gene_muts$pos[i])
        gene_muts$pos[i] <- gene_coverage$pos[which(diff_in_pos == min(diff_in_pos))][1]
      }
    }
    
    gene_start_loss <- which(gene_muts$impact != "Synonymous" & substr(gene_muts$aachange,2,nchar(gene_muts$aachange)-1) == "1")
    gene_muts$impact[gene_start_loss] <- "Start_loss"
    
    mut_types <- names(table(gene_muts$impact))
    
    fixed_bin_size <- 3
    variable_bin_size <- max(round(immune_genes$CDS_length[gene_num] / 1500),1) * 3
    
    cosmic_muts <- read.table(paste0(cosmic_dir_path,"/COSMIC_v99_WholeGenomeExome_FullAnnotated_chr",unique(gene_muts$chr),".tsv"), sep = "\t", stringsAsFactors=F, header = T)
    cosmic_gene_muts <- cosmic_muts[which(cosmic_muts$gene == gene),]
    
    if(nrow(cosmic_gene_muts) > 0){
      cosmic_start_loss <- which(cosmic_gene_muts$impact != "Synonymous" & substr(cosmic_gene_muts$aachange,2,nchar(cosmic_gene_muts$aachange)-1) == "1")
      cosmic_gene_muts$impact[cosmic_start_loss] <- "Start_loss"
      
      mut_types <- unique(c(mut_types,names(table(cosmic_gene_muts$impact))))
      
      for(i in 1:nrow(cosmic_gene_muts)){
        if(!(cosmic_gene_muts$pos[i] %in% gene_coverage$pos)){
          diff_in_pos <- abs(gene_coverage$pos - cosmic_gene_muts$pos[i])
          cosmic_gene_muts$pos[i] <- gene_coverage$pos[which(diff_in_pos == min(diff_in_pos))][1]
        }
      }
      
      cosmic_sample_info <- read.table(paste0(cosmic_dir_path,"/COSMIC_v99_WholeGenomeExome_SamplePhenotypeInfo.tsv"), sep = "\t", stringsAsFactors = F, header = T, comment = "", quote = "")
      cosmic_blood_samples <- cosmic_sample_info$COSMIC_SAMPLE_ID[which(cosmic_sample_info$PRIMARY_HISTOLOGY %in% c("lymphoid_neoplasm","haematopoietic_neoplasm"))]
      cosmic_blood_gene_muts <- cosmic_gene_muts[which(cosmic_gene_muts$sampleID %in% cosmic_blood_samples),]
    }
    
    fixed_binned_muts <- as.data.frame(array(data = NA, dim = c(ceiling(nrow(gene_coverage) / fixed_bin_size) * length(mut_types),3)))
    colnames(fixed_binned_muts) <- c("idx","mut_type","count")
    
    fixed_binned_muts$idx <- rep(seq(from = (fixed_bin_size)/2, to = (ceiling(nrow(gene_coverage) / fixed_bin_size) * fixed_bin_size) - (fixed_bin_size)/2, by = fixed_bin_size), each = length(mut_types))
    fixed_binned_muts$mut_type <- rep(mut_types, times = ceiling(nrow(gene_coverage) / fixed_bin_size))
    
    if(nrow(cosmic_gene_muts) > 0){
      cosmic_fixed_binned_muts <- fixed_binned_muts
      
      for(i in 1:nrow(cosmic_fixed_binned_muts)){
        if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
          cosmic_fixed_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_fixed_binned_muts$mut_type[i] & cosmic_gene_muts$pos >= gene_coverage$pos[cosmic_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_gene_muts$pos <= gene_coverage$pos[min(cosmic_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
        }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
          cosmic_fixed_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_fixed_binned_muts$mut_type[i] & cosmic_gene_muts$pos <= gene_coverage$pos[cosmic_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_gene_muts$pos >= gene_coverage$pos[min(cosmic_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
        }
      }
      
      if(nrow(cosmic_blood_gene_muts) > 0){
        cosmic_blood_fixed_binned_muts <- fixed_binned_muts
        
        for(i in 1:nrow(cosmic_blood_fixed_binned_muts)){
          if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
            cosmic_blood_fixed_binned_muts$count[i] <- length(which(cosmic_blood_gene_muts$impact == cosmic_blood_fixed_binned_muts$mut_type[i] & cosmic_blood_gene_muts$pos >= gene_coverage$pos[cosmic_blood_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_blood_gene_muts$pos <= gene_coverage$pos[min(cosmic_blood_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
          }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
            cosmic_blood_fixed_binned_muts$count[i] <- length(which(cosmic_blood_gene_muts$impact == cosmic_blood_fixed_binned_muts$mut_type[i] & cosmic_blood_gene_muts$pos <= gene_coverage$pos[cosmic_blood_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_blood_gene_muts$pos >= gene_coverage$pos[min(cosmic_blood_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
          }
        }
      }
    }
    
    for(i in 1:nrow(fixed_binned_muts)){
      if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
        fixed_binned_muts$count[i] <- length(which(gene_muts$impact == fixed_binned_muts$mut_type[i] & gene_muts$pos >= gene_coverage$pos[fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & gene_muts$pos <= gene_coverage$pos[min(fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
      }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
        fixed_binned_muts$count[i] <- length(which(gene_muts$impact == fixed_binned_muts$mut_type[i] & gene_muts$pos <= gene_coverage$pos[fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & gene_muts$pos >= gene_coverage$pos[min(fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
      }
    }
    
    variable_binned_muts <- as.data.frame(array(data = NA, dim = c(ceiling(nrow(gene_coverage) / variable_bin_size) * length(mut_types),3)))
    colnames(variable_binned_muts) <- c("idx","mut_type","count")
    
    variable_binned_muts$idx <- rep(seq(from = (variable_bin_size)/2, to = (ceiling(nrow(gene_coverage) / variable_bin_size) * variable_bin_size) - (variable_bin_size)/2, by = variable_bin_size), each = length(mut_types))
    variable_binned_muts$mut_type <- rep(mut_types, times = ceiling(nrow(gene_coverage) / variable_bin_size))
    
    if(nrow(cosmic_gene_muts) > 0){
      cosmic_variable_binned_muts <- variable_binned_muts
      
      for(i in 1:nrow(cosmic_variable_binned_muts)){
        if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
          cosmic_variable_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_variable_binned_muts$mut_type[i] & cosmic_gene_muts$pos >= gene_coverage$pos[cosmic_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_gene_muts$pos <= gene_coverage$pos[min(cosmic_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
        }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
          cosmic_variable_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_variable_binned_muts$mut_type[i] & cosmic_gene_muts$pos <= gene_coverage$pos[cosmic_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_gene_muts$pos >= gene_coverage$pos[min(cosmic_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
        }
      }
      
      if(nrow(cosmic_blood_gene_muts) > 0){
        cosmic_blood_variable_binned_muts <- variable_binned_muts
        
        for(i in 1:nrow(cosmic_blood_variable_binned_muts)){
          if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
            cosmic_blood_variable_binned_muts$count[i] <- length(which(cosmic_blood_gene_muts$impact == cosmic_blood_variable_binned_muts$mut_type[i] & cosmic_blood_gene_muts$pos >= gene_coverage$pos[cosmic_blood_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_blood_gene_muts$pos <= gene_coverage$pos[min(cosmic_blood_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
          }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
            cosmic_blood_variable_binned_muts$count[i] <- length(which(cosmic_blood_gene_muts$impact == cosmic_blood_variable_binned_muts$mut_type[i] & cosmic_blood_gene_muts$pos <= gene_coverage$pos[cosmic_blood_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_blood_gene_muts$pos >= gene_coverage$pos[min(cosmic_blood_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
          }
        }
      }
    }
    
    for(i in 1:nrow(variable_binned_muts)){
      if(immune_genes$strand[which(immune_genes$gene == gene)] == 1){
        variable_binned_muts$count[i] <- length(which(gene_muts$impact == variable_binned_muts$mut_type[i] & gene_muts$pos >= gene_coverage$pos[variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & gene_muts$pos <= gene_coverage$pos[min(variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
      }else if(immune_genes$strand[which(immune_genes$gene == gene)] == -1){
        variable_binned_muts$count[i] <- length(which(gene_muts$impact == variable_binned_muts$mut_type[i] & gene_muts$pos <= gene_coverage$pos[variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & gene_muts$pos >= gene_coverage$pos[min(variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
      }
    }
    
    interpro_files <- list.files(interpro_domain_dir_path)
    domain_file <- interpro_files[which(grepl(paste0(gene,"_"),interpro_files) & !(grepl("DifferentIsoform",interpro_files)))]
    
    if(length(domain_file) == 1){
      domains <- read.table(file = paste0(interpro_domain_dir_path,"/",domain_file), sep = "\t", stringsAsFactors = F, header = T, quote = "")
      domains <- domains[which(domains$Source.Database == "interpro" & domains$Type %in% c("domain","repeat")),]
      
      if(nrow(domains) > 0){
        rows_to_split <- which(grepl(",",domains$Matches))
        
        if(length(rows_to_split) >= 1){
          for(i in 1:length(rows_to_split)){
            temp_rows <- as.data.frame(array(data = NA, dim = c(lengths(regmatches(domains$Matches[rows_to_split[i]], gregexpr(",", domains$Matches[rows_to_split[i]]))) + 1,ncol(domains))))
            colnames(temp_rows) <- colnames(domains)
            for(j in 1:nrow(temp_rows)){
              temp_rows[j,which(colnames(temp_rows) != "Matches")] <- domains[rows_to_split[i],which(colnames(domains) != "Matches"),]
            }
            temp_rows$Matches <- unlist(strsplit(domains$Matches[rows_to_split[i]], split = ","))
            
            domains <- rbind(domains,temp_rows)
          }
          
          domains <- domains[-(rows_to_split),]
        }
        
        domains$plot <- NA
        domains$aa_start <- NA
        domains$aa_end <- NA
        
        for(i in 1:nrow(domains)){
          domains$aa_start[i] <- as.numeric(unlist(strsplit(domains$Matches[i], split = "\\.."))[1])
          domains$aa_end[i] <- as.numeric(unlist(strsplit(domains$Matches[i], split = "\\.."))[2])
        }
        
        domains$idx_start <- (3 * (domains$aa_start - 1)) + 1
        domains$idx_end <- (3 * (domains$aa_end))
        
        if(unique(domains$Protein.Length) != ((nrow(gene_coverage) / 3) - 1)){
          domains$plot <- 0
        }else{
          domains$plot <- 1
        }
        
        domains$overlaps <- NA
        domains$encapsulates <- NA
        domains$num_encapsulated <- NA
        domains$high_deg_overlap <- NA
        
        for(i in 1:nrow(domains)){
          overlaps <- NULL
          encapsulates <- NULL
          high_deg_overlap <- NULL
          
          for(j in 1:nrow(domains)){
            if(i != j){
              per_aa_i <- seq(from = domains$aa_start[i], to = domains$aa_end[i], by = 1)
              per_aa_j <- seq(from = domains$aa_start[j], to = domains$aa_end[j], by = 1)
              
              if(length(which(per_aa_j %in% per_aa_i)) > 0){
                overlaps <- c(overlaps,j)
              }
              
              if(length(which(per_aa_j %in% per_aa_i)) == length(per_aa_j)){
                encapsulates <- c(encapsulates,j)
              }
              
              if(length(which(per_aa_j %in% per_aa_i)) > 0.8 * length(per_aa_j) & length(which(per_aa_i %in% per_aa_j)) > 0.8 * length(per_aa_i)){
                high_deg_overlap <- c(high_deg_overlap,j)
              }
            }
            
            domains$overlaps[i] <- paste(overlaps, collapse = ",")
            domains$encapsulates[i] <- paste(encapsulates, collapse = ",")
            domains$num_encapsulated[i] <- length(encapsulates)
            domains$high_deg_overlap[i] <- paste(high_deg_overlap, collapse =  ",")
            
            if(domains$plot[i] != 0){
              domains$plot[high_deg_overlap] <- 0
            }
          }
        }
        
        domains <- domains[order(domains$num_encapsulated, decreasing = T),]
        
        if(length(which(domains$plot ==1)) > 0){
          domains_plot <- ggplot(data = domains[which(domains$plot == 1),]) +
            geom_rect(xmin = min(gene_pos_transition_table$start_idx), xmax = max(gene_pos_transition_table$end_idx), ymin = -0.15, ymax = 0.15) +
            geom_rect(aes(xmin = idx_start, xmax = idx_end, fill = Name), color = "black", ymin = -0.5, ymax = 0.5) +
            scale_fill_brewer(palette = "Set3") +
            scale_x_continuous(expand = c(0,6), limits = c(min(gene_pos_transition_table$start_idx),max(gene_pos_transition_table$end_idx))) +
            scale_y_continuous(expand = c(0,0), limits = c(-0.5,0.5)) +
            theme_bw() +
            labs(fill = "Domain type", y = "Domains") +
            theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5), panel.background = element_rect(fill='transparent'))
        }
      }
    }
    
    fixed_binned_muts_limit <- ceiling(max(aggregate(fixed_binned_muts$count ~ fixed_binned_muts$idx, data=fixed_binned_muts, sum)[,2]) * 1.1)
    mut_fixed_histo_no_cosmic_plot <- ggplot() +
      geom_area(data = gene_coverage_smooth, aes(x = idx, y = (bin_average * fixed_binned_muts_limit / (ceiling((max(total))/10000) * 10000))), color = "grey70", fill = "grey", alpha = 0.25) +
      geom_bar(data = fixed_binned_muts, aes(x = idx, y = count, fill = mut_type), stat = "identity", width = fixed_bin_size) +
      scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,6)) +
      labs(x = "", fill = "Mutation type", y = paste0(title,"\nNanoSeq\nmutation\ncount"), title = paste0(gene," - ",(fixed_bin_size / 3)," aa bins\n",sum(fixed_binned_muts$count)," ",title," mutations")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      theme(legend.key=element_rect(colour="black")) +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_text(angle = 0)) +
      scale_y_continuous(expand = c(0,0), 
                         limits = c(0,fixed_binned_muts_limit),
                         breaks = pretty_breaks(),
                         sec.axis = sec_axis(~ . * ((ceiling(max(gene_coverage_smooth$total)/10000) * 10000)/(ceiling(max(aggregate(fixed_binned_muts$count ~ fixed_binned_muts$idx, data=fixed_binned_muts, sum)[,2]) * 1.1))), name = "Cumulative\nduplex coverage"))
    
    variable_binned_muts_limit <- ceiling(max(aggregate(variable_binned_muts$count ~ variable_binned_muts$idx, data=variable_binned_muts, sum)[,2]) * 1.1)
    
    mut_variable_histo_no_cosmic_plot <- ggplot() +
      geom_area(data = gene_coverage_smooth, aes(x = idx, y = (bin_average * variable_binned_muts_limit / (ceiling((max(total))/10000) * 10000))), color = "grey70", fill = "grey", alpha = 0.25) +
      geom_bar(data = variable_binned_muts, aes(x = idx, y = count, fill = mut_type), stat = "identity", width = variable_bin_size) +
      scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,6)) +
      labs(x = "", fill = "Mutation type", y = paste0(title,"\nNanoSeq\nmutation\ncount"), title = paste0(gene," - ",(variable_bin_size / 3)," aa bins\n",sum(variable_binned_muts$count)," ",title," mutations")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      theme(legend.key=element_rect(colour="black")) +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_text(angle = 0)) +
      scale_y_continuous(expand = c(0,0), 
                         limits = c(0,variable_binned_muts_limit),
                         breaks = pretty_breaks(),
                         sec.axis = sec_axis(~ . * ((ceiling(max(gene_coverage_smooth$total)/10000) * 10000)/(ceiling(max(aggregate(variable_binned_muts$count ~ variable_binned_muts$idx, data=variable_binned_muts, sum)[,2]) * 1.1))), name = "Cumulative\nduplex coverage"))
    
    if(nrow(gene_pos_transition_table) >= 2){
      exome_break_plot <- ggplot() +
        geom_rect(data = gene_pos_transition_table[seq(from = 1, to = nrow(gene_pos_transition_table), by = 2),], aes(xmin = start_idx, xmax = end_idx), fill = "firebrick", color = "black", ymin = 0, ymax = 1, alpha = 0.35) +
        geom_rect(data = gene_pos_transition_table[seq(from = 2, to = nrow(gene_pos_transition_table), by = 2),], aes(xmin = start_idx, xmax = end_idx), fill = "firebrick", color = "black", ymin = 0, ymax = 1, alpha = 0.65) +
        scale_x_continuous(expand = c(0,6)) +
        theme_bw() +
        scale_y_continuous(expand = c(0,0)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
        labs(y = "Exons") +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        theme(panel.border = element_rect(color = "black", fill = NA))
    }else{
      exome_break_plot <- ggplot() +
        geom_rect(data = gene_pos_transition_table[1,], aes(xmin = start_idx, xmax = end_idx), fill = "firebrick", color = "black", ymin = 0, ymax = 1, alpha = 0.35) +
        scale_x_continuous(expand = c(0,6)) +
        theme_bw() +
        scale_y_continuous(expand = c(0,0)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
        labs(y = "Exons") +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        theme(panel.border = element_rect(color = "black", fill = NA))
    }
    
    if(nrow(cosmic_gene_muts) > 0){
      
      cosmic_fixed_binned_muts_limit <- ceiling(max(aggregate(cosmic_fixed_binned_muts$count ~ cosmic_fixed_binned_muts$idx, data=cosmic_fixed_binned_muts, sum)[,2]) * 1.1)
      cosmic_fixed_breaks <- pretty_breaks()(-cosmic_fixed_binned_muts_limit:0)
      
      mut_fixed_histo_cosmic_plot <- ggplot() +
        geom_bar(data = cosmic_fixed_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = fixed_bin_size) +
        scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
        theme_bw() +
        scale_x_continuous(expand = c(0,6)) +
        labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_fixed_binned_muts$count)," COSMIC (all cancer types) mutations")) +
        theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
        theme(axis.text.x = element_blank()) +
        theme(axis.title.x = element_blank()) +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        theme(legend.position = "none") +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        scale_y_continuous(expand = c(0,0), 
                           limits = c(-cosmic_fixed_binned_muts_limit,0),
                           breaks = cosmic_fixed_breaks,
                           labels = -(cosmic_fixed_breaks))
      
      cosmic_variable_binned_muts_limit <- ceiling(max(aggregate(cosmic_variable_binned_muts$count ~ cosmic_variable_binned_muts$idx, data=cosmic_variable_binned_muts, sum)[,2]) * 1.1)
      cosmic_variable_breaks <- pretty_breaks()(-cosmic_variable_binned_muts_limit:0)
      
      mut_variable_histo_cosmic_plot <- ggplot() +
        geom_bar(data = cosmic_variable_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = variable_bin_size) +
        scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
        theme_bw() +
        scale_x_continuous(expand = c(0,6)) +
        labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_variable_binned_muts$count)," COSMIC (all cancer types) mutations")) +
        theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
        theme(axis.text.x = element_blank()) +
        theme(axis.title.x = element_blank()) +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        theme(legend.position = "none") +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        scale_y_continuous(expand = c(0,0), 
                           limits = c(-cosmic_variable_binned_muts_limit,0),
                           breaks = cosmic_variable_breaks,
                           labels = -(cosmic_variable_breaks))
      
      if(nrow(cosmic_blood_gene_muts) > 0){
        cosmic_blood_fixed_binned_muts_limit <- ceiling(max(aggregate(cosmic_blood_fixed_binned_muts$count ~ cosmic_blood_fixed_binned_muts$idx, data=cosmic_blood_fixed_binned_muts, sum)[,2]) * 1.1)
        cosmic_blood_fixed_breaks <- pretty_breaks()(-cosmic_blood_fixed_binned_muts_limit:0)
        
        mut_fixed_histo_cosmic_blood_plot <- ggplot() +
          geom_bar(data = cosmic_blood_fixed_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = fixed_bin_size) +
          scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
          theme_bw() +
          scale_x_continuous(expand = c(0,6)) +
          labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_blood_fixed_binned_muts$count)," COSMIC (lymphoid and haematopoietic neoplasms) mutations")) +
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
          theme(axis.text.x = element_blank()) +
          theme(axis.title.x = element_blank()) +
          theme(axis.ticks.x = element_blank()) +
          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
          theme(legend.position = "none") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
          scale_y_continuous(expand = c(0,0),
                             limits = c(-cosmic_blood_fixed_binned_muts_limit,0),
                             breaks = cosmic_blood_fixed_breaks,
                             labels = -(cosmic_blood_fixed_breaks))
        
        cosmic_blood_variable_binned_muts_limit <- ceiling(max(aggregate(cosmic_blood_variable_binned_muts$count ~ cosmic_blood_variable_binned_muts$idx, data=cosmic_blood_variable_binned_muts, sum)[,2]) * 1.1)
        cosmic_blood_variable_breaks <- pretty_breaks()(-cosmic_blood_variable_binned_muts_limit:0)
        
        mut_variable_histo_cosmic_blood_plot <- ggplot() +
          geom_bar(data = cosmic_blood_variable_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = variable_bin_size) +
          scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
          theme_bw() +
          scale_x_continuous(expand = c(0,6)) +
          labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_blood_variable_binned_muts$count)," COSMIC (lymphoid and haematopoietic neoplasms) mutations")) +
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
          theme(axis.text.x = element_blank()) +
          theme(axis.title.x = element_blank()) +
          theme(axis.ticks.x = element_blank()) +
          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
          theme(legend.position = "none") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
          scale_y_continuous(expand = c(0,0),
                             limits = c(-cosmic_blood_variable_binned_muts_limit,0),
                             breaks = cosmic_blood_variable_breaks,
                             labels = -(cosmic_blood_variable_breaks))
      }
    }
    
    if(!(exists("domains_plot"))){
      fixed_no_cosmic_no_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot + plot_layout(heights = c(20,-1.05,1), guides = "collect")
      if (runman) {ggsave(fixed_no_cosmic_no_domain_plot, file = paste0(plotfilename,"_fixed_no_cosmic_no_domain.pdf"),width = 15, height = 5)}
      
      if(exists("mut_fixed_histo_cosmic_plot")){
        fixed_cosmic_all_cancer_no_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_fixed_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
        if (runman) {ggsave(fixed_cosmic_all_cancer_no_domain_plot, file = paste0(plotfilename,"_fixed_cosmic_all_no_domain_blood.pdf"),width = 15, height = 10)}
      }
      if(exists("mut_fixed_histo_cosmic_blood_plot")){
        fixed_cosmic_blood_cancer_no_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_fixed_histo_cosmic_blood_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
        if (runman) {ggsave(fixed_cosmic_blood_cancer_no_domain_plot, file = paste0(plotfilename,"_fixed_cosmic_blood_no_domain_blood.pdf"),width = 15, height = 10)}
        if(variable_bin_size == fixed_bin_size) {print(fixed_cosmic_blood_cancer_no_domain_plot)}
      }
    }else {
      fixed_no_cosmic_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot + plot_layout(heights = c(20,-1.05,1,1), guides = "collect")
      if (runman) {ggsave(fixed_no_cosmic_domain_plot, file = paste0(plotfilename,"_fixed_no_cosmic_domain_blood.pdf"),width = 15, height = 5)}
      
      if(exists("mut_fixed_histo_cosmic_plot")){
        fixed_cosmic_all_cancer_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_fixed_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
        if (runman) {ggsave(fixed_cosmic_all_cancer_domain_plot, file = paste0(plotfilename,"_fixed_cosmic_all_domain_blood.pdf"),width = 15, height = 10)}
      }
      if(exists("mut_fixed_histo_cosmic_blood_plot")){
        fixed_cosmic_blood_cancer_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_fixed_histo_cosmic_blood_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
        if (runman) {ggsave(fixed_cosmic_blood_cancer_domain_plot, file = paste0(plotfilename,"_fixed_cosmic_blood_domain_blood.pdf"),width = 15, height = 10)}
        if(variable_bin_size == fixed_bin_size) {print(fixed_cosmic_blood_cancer_domain_plot)}
      }
    }
    
    if(variable_bin_size != fixed_bin_size){
      if(!(exists("domains_plot"))){
        variable_no_cosmic_no_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot + plot_layout(heights = c(20,-1.05,1), guides = "collect")
        if (runman) {ggsave(variable_no_cosmic_no_domain_plot, file = paste0(plotfilename,"_variable_no_cosmic_no_domain_blood.pdf"),width = 15, height = 5)}
        
        if(exists("mut_fixed_histo_cosmic_plot")){
          variable_cosmic_all_cancer_no_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_variable_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
          if (runman) {ggsave(variable_cosmic_all_cancer_no_domain_plot, file = paste0(plotfilename,"_variable_cosmic_all_no_domain_blood.pdf"),width = 15, height = 10)}
        }
        if(exists("mut_fixed_histo_cosmic_blood_plot")){
          variable_cosmic_blood_cancer_no_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_variable_histo_cosmic_blood_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
          if (runman) {ggsave(variable_cosmic_blood_cancer_no_domain_plot, file = paste0(plotfilename,"_variable_cosmic_blood_no_domain_blood.pdf"),width = 15, height = 10)}
          print(variable_cosmic_blood_cancer_no_domain_plot)
        }
      }else {
        variable_no_cosmic_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot + plot_layout(heights = c(20,-1.05,1,1), guides = "collect")
        if (runman) {ggsave(variable_no_cosmic_domain_plot, file = paste0(plotfilename,"_variable_no_cosmic_domain_blood.pdf"),width = 15, height = 5)}
        
        if(exists("mut_fixed_histo_cosmic_plot")){
          variable_cosmic_all_cancer_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_variable_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
          if (runman) {ggsave(variable_cosmic_all_cancer_domain_plot, file = paste0(plotfilename,"_variable_cosmic_all_domain_blood.pdf"),width = 15, height = 10)}
        }
        if(exists("mut_fixed_histo_cosmic_blood_plot")){
          variable_cosmic_blood_cancer_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_variable_histo_cosmic_blood_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
          if (runman) {ggsave(variable_cosmic_blood_cancer_domain_plot, file = paste0(plotfilename,"_coding_variable_cosmic_blood_domain_blood.pdf"),width = 15, height = 10)}
          print(variable_cosmic_blood_cancer_domain_plot)
        }
      }
    }
  }
}

#' withingenednds
#' 
#' This function uses Poisson and Negative Binomial regression models at single-site level to study selection across different regions (coding and non-coding) within a gene.
#' 
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param mutations Data frame with all the mutations detected in the study (5-column input table as for dndscv: sampleID, chr, pos, ref, mut).
#' @param gene Name of the gene of interest. This function is currently designed to work on a single gene, but combined analyses of multiple genes could be done using the sites output table generated by this function.
#' @param covtable Table with all sites of interest in the gene. This should be a data frame with one row per site and the following columns: chr, pos, dc (duplex depth). Additional columns will not be used.
#' @param dndsout dndscv output object for the dataset. This is mainly used for the MLEs of the substitution model. Running dndscv on all genes in the dataset is recommended unless the gene of interest is believed to have a different substitution model.
#' @param genomeFile Path to a reference fasta file for the genome assembly.
#' @param regionschr Optional data frame with user-defined regions of interest in the gene. This allows the user to define arbitrary regions within a gene (coding or non-coding) from which to calculate omega (selection or obs/exp) values (e.g. protein domains, splicing regulator regions, core promoters, etc). The table should contain the following columns: chr, start, end, wname (a unique name for the w parameter, e.g. wdomain1, wcorepromoter), impacts (e.g. Missense or Missense|Nonsense will restrict the w calculation with Missense or Missense and Nonsense mutations in the region, respectively), layered (1/0; using "0" removes other w parameters influencing the site, whereas using "1" models selection as relative to other w parameters active at these sites).
#' @param regionsaa Optional data frame with user-defined regions of interest in the gene, using aminoacid coordinates. The table should contain the following columns: gene, aa_start, aa_end, w feature name (e.g. wdomain1), impacts.
#' @param fixtheta Pre-calculated overdispersion (theta) parameter. This should be calculated using sitednds(., method="NB").
#' @param refdb Reference database (path to .rda file or a pre-loaded array object in the right format).
#' @param numcode NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate
#' @param normalisefromsyn Normalise the substitution rates based on the synonymous mutations in the gene. Using TRUE is recommended. Using FALSE uses the expected synonymous mutation rate of the gene from the dndscv negative binomial regression model (dndsout$genemuts). 
#' @param syndrivers Vector of known synonymous driver sites defined by their aminoacid position, to be excluded from the background model (e.g. syndrivers = c("T125T","E224E","Q331Q") for TP53).
#' @param exon_flank_length Exon flank length in bp [default = 10]. Using a value higher than 0 will calculate a separate selection (w) coefficient for synonymous mutations in exon flanks.
#' @param intron_flank_length Intron flank length in bp [default = 10]. Intronic sites occurring within these flanks but not already classified as Essential_Splice will receive a separate w parameter.
#' @param sitefilename Optionally, provide a file name to save the table of all annotated sites in the gene. This table is also always contained in the output object. 
#'
#' @return 'withingenednds' returns a list of objects:
#' @return - sites: Table with the annotation of all sites in the gene (from covtable), including all functional annotations in the "regions" input object as well as default annotations (Missense, Nonsense, Essential_Splice, Start_loss, Stop_loss, etc).
#' @return - par.pois: Poisson regression results (not recommended).
#' @return - par.nb: Negative binomial results fitting a new overdispersion parameter to the data (when fixtheta is not provided).
#' @return - par.nbfix: Negative binomial results using the input fixtheta value as recommended.
#' @return - model.pois: Poisson regression object.
#' @return - model.nb: Negative binomial regression object.
#' @return - model.nbfix: Negative binomial regression object.
#'
#' @export

withingenednds = function(mutations, gene, covtable, dndsout, genomeFile, regionschr = NULL, regionsaa = NULL, fixtheta = NULL, normalisefromsyn = TRUE, syndrivers = NULL, exon_flank_length = 10, intron_flank_length = 10, sitefilename = NULL, refdb = "hg19", numcode = 1) {
  
  ## 1. Environment
  message("Running gene-level selection analyses...")
  
  # [Input] Reference database
  refdb_class = class(refdb)
  if ("character" %in% refdb_class) {
    if (refdb == "hg19") {
      data("refcds_hg19", package="dndscv")
    } else {
      load(refdb)
    }
  } else if("array" %in% refdb_class) {
    # use the user-supplied RefCDS object
    RefCDS = refdb
  } else {
    stop("Expected refdb to be \"hg19\", a file path, or a RefCDS-formatted array object.")
  }
  
  # Annotating all possible coding changes in the positions provided in the input covtable
  library("Rsamtools")
  covtable$ref = as.vector(scanFa(genomeFile, GRanges(covtable$chr, IRanges(covtable$pos, covtable$pos))))
  covtable$ref3 = as.vector(scanFa(genomeFile, GRanges(covtable$chr, IRanges(covtable$pos-1, covtable$pos+1))))
  
  allsubs = data.frame(sampleID="allsubs", chr = rep(covtable$chr, each=3), pos = rep(covtable$pos, each=3), ref = rep(covtable$ref, each=3), mut = NA, ref3 = rep(covtable$ref3, each=3), mut3 = NA, ref_cod = NA, mut_cod = NA, ref3_cod = NA, mut3_cod = NA, stringsAsFactors = F)
  for (j in seq(1,nrow(allsubs),3)) {
    allsubs$mut[j:(j+2)] = setdiff(c("A","C","G","T"),allsubs$ref[j])
  }
  allsubs$mut3 = paste(substr(allsubs$ref3,1,1), allsubs$mut, substr(allsubs$ref3,3,3), sep="")
  aux = dndscv(allsubs, gene_list = gene, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf, outp = 1)$annotmuts # Annotated mutations
  
  # Adding duplex depth, functional impact annotation to all possible changes in the input table
  aux$mstr = paste(aux$chr, aux$pos, aux$mut, sep=":") # mutation string ID
  obssubs = dndsout$annotmuts[which(dndsout$annotmuts$ref %in% c("A","C","G","T") & dndsout$annotmuts$mut %in% c("A","C","G","T") & dndsout$annotmuts$gene == gene), ]
  obssubs$mstr = paste(obssubs$chr, obssubs$pos, obssubs$mut, sep=":")
  mutations$mstr = paste(mutations$chr, mutations$pos, mutations$mut, sep=":")
  
  sites = cbind(data.frame(gene=gene, pid=aux$pid[1], stringsAsFactors = F), allsubs[,-1]) # Initialising the sites table
  pos2dc = setNames(covtable$dc, covtable$pos)
  sites$dc = pos2dc[as.character(sites$pos)]
  sites$mstr = paste(sites$chr, sites$pos, sites$mut, sep=":")
  sites$obs = table(mutations$mstr)[sites$mstr]; sites$obs[is.na(sites$obs)] = 0 # Annotating the number of observed mutations at each site
  
  mstr2imp1 = setNames(aux$ntchange, aux$mstr)
  mstr2imp2 = setNames(aux$aachange, aux$mstr)
  mstr2imp3 = setNames(aux$impact, aux$mstr)
  sites$ntchange = mstr2imp1[sites$mstr]
  sites$aachange = mstr2imp2[sites$mstr]
  sites$impact = mstr2imp3[sites$mstr]
  
  # Adding strand and annotating the coding trinucleotide change
  aux2 = unique(dndsout$annotmuts[,c("gene","strand")])
  gene2strand = setNames(aux2[,2],aux2[,1])
  sites$strand = gene2strand[sites$gene]
  
  sites$ref_cod = sites$ref
  sites$ref_cod[sites$strand==-1] = seqinr::comp(sites$ref[sites$strand==-1], forceToLower = F)
  sites$mut_cod = sites$mut
  sites$mut_cod[sites$strand==-1] = seqinr::comp(sites$mut[sites$strand==-1], forceToLower = F)
  
  revcomp = function(seqvec) {
    as.vector(sapply(seqvec, function(x) paste(seqinr::comp(rev(unlist(strsplit(x,split=""))), forceToLower=F), collapse=""))) # Reverse complement of a vector of sequence motifs
  }
  sites$ref3_cod = sites$ref3
  sites$mut3_cod = sites$mut3
  if (any(sites$strand==-1)) {
    sites$ref3_cod[sites$strand==-1] = revcomp(sites$ref3[sites$strand==-1])
    sites$mut3_cod[sites$strand==-1] = revcomp(sites$mut3[sites$strand==-1])
  }
  
  # Adding the trinucleotide substitution rates from the dndsout model
  
  mle_submodel = setNames(dndsout$mle_submodel[,2], dndsout$mle_submodel[,1])
  mle_submodel = c(mle_submodel, "TTT>TGT"=1) # Adding the missing rate (note that all rates in mle_submodel in dNdScv are relative to TTT>TGT)
  mle_submodel = mle_submodel * mle_submodel["t"] # Absolute rates
  sites$r = mle_submodel[paste(sites$ref3_cod, sites$mut3_cod, sep=">")]
  sites$r = sites$r * sites$dc / mean(sites$dc) # Normalising the expected rates at a site by the duplex depth
  
  # Normalising the global expected rates by the estimated mutation rate of the gene using one of two alternative models:
  # 1. normalisefromsyn = TRUE. We normalise "r" using the synonymous mutations observed in the gene (excluding known syn driver sites)
  # 2. normalisefromsyn = FALSE. We normalise "r" using the negative binomial model from dndscv.
  
  if (normalisefromsyn == TRUE) {
    mobs = sum(sites$obs[which(sites$impact=="Synonymous" & !(sites$aachange %in% syndrivers))])
    mexp = sum(sites$r[which(sites$impact=="Synonymous" & !(sites$aachange %in% syndrivers))])
    sites$rnorm = sites$r * mobs / mexp
  } else {
    message("Option not recommended: Normalising the mutation rate of the gene based on the negative regression model of dNdScv")
    sites$rnorm = sites$r * dndsout$genemuts$exp_syn_cv[dndsout$genemuts$gene_name==gene] / dndsout$genemuts$exp_syn[dndsout$genemuts$gene_name==gene]
  }
  
  # Excluding sites with rate = 0 (e.g. due to lack of duplex coverage)
  ratezero = which(sites$rnorm==0)
  if (length(ratezero)>0) {
    sites = sites[-ratezero, ]
    message(sprintf("%0.0f sites have been excluded due to having a duplex depth and/or a predicted mutation rate of 0", length(ratezero)))
  }
  
  ## 2. Creating an index regression matrix with all functional annotations (each row is a site and each column is a parameter in the selection model)
  
  # Annotating Missense, Nonsense, Essential_Splice and Stop_loss mutations
  mutmat = data.frame(t = rep(1,nrow(sites)), wmis = 0, wnon = 0, wspl = 0, wsyndriv = 0, 
                      wintr = 0, woutcds = 0, wexfl = 0, winfl = 0, wnonlastex = 0, wstoploss = 0, wstartloss = 0, stringsAsFactors = F)
  mutmat$wmis[which(sites$impact=="Missense")] = 1
  mutmat$wnon[which(sites$impact=="Nonsense")] = 1
  mutmat$wspl[which(sites$impact=="Essential_Splice")] = 1
  mutmat$wsyndriv[which(sites$impact=="Synonymous" & sites$aachange %in% syndrivers)] = 1
  mutmat$wstoploss[which(sites$impact=="Stop_loss")] = 1
  
  # Annotating Start_loss mutations (mutations in codon 1)
  
  #start_loss = which(sites$impact != "Synonymous" & substr(sites$aachange,2,nchar(sites$aachange)-1) == "1") # There is no need to check for synonymous changes in ATG
  start_loss = which(substr(sites$aachange,2,nchar(sites$aachange)-1) == "1") # There is no need to check for synonymous changes in ATG
  mutmat$wstartloss[start_loss] = 1
  sites$impact[start_loss] = "Start_loss"
  
  # Annotating introns, exon flanks, and intron flanks
  
  refcdsgene = RefCDS[[which(sapply(RefCDS, function(x) x$gene_name)==gene)]]
  exons = refcdsgene$intervals_cds
  esspl = refcdsgene$intervals_splice
  
  gr_sites = GenomicRanges::GRanges(sites$gene, IRanges::IRanges(sites$pos,sites$pos))
  gr_exons = GenomicRanges::GRanges(gene, IRanges::IRanges(exons[,1],exons[,2]))
  gr_cds = GenomicRanges::GRanges(gene, IRanges::IRanges(min(exons),max(exons)))
  gr_outcds = setdiff(gr_sites, gr_cds)
  ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_outcds, type="any", select="all"))
  mutmat$woutcds[unique(ol[,1])] = 1 # Annotation of the intronic sites
  
  # For genes with >1 exon, we also annotate essential splice sites, and intron and exon flanks.
  if (length(esspl)>0) {
    gr_esspl = GenomicRanges::GRanges(gene, IRanges::IRanges(esspl,esspl))
    gr_introns = setdiff(setdiff(gr_cds, gr_exons), gr_esspl)
    exfl = rbind(cbind(exons[2:nrow(exons),1],exons[2:nrow(exons),1]+exon_flank_length-1), cbind(exons[1:(nrow(exons)-1),2]-exon_flank_length+1,exons[1:(nrow(exons)-1),2]))
    gr_exonflanks = GenomicRanges::GRanges(gene, IRanges::IRanges(exfl[,1],exfl[,2]))
    gr_exonflanks = intersect(gr_exonflanks, gr_exons) # Ensuring that the exon flanks do not extend into introns
    infl = rbind(cbind(exons[2:nrow(exons),1]-intron_flank_length,exons[2:nrow(exons),1]-1), cbind(exons[1:(nrow(exons)-1),2]+1,exons[1:(nrow(exons)-1),2]+intron_flank_length))
    gr_intrflanks = GenomicRanges::GRanges(gene, IRanges::IRanges(infl[,1],infl[,2]))
    gr_intrflanks = setdiff(gr_intrflanks, gr_exons) # Removing any overlaps with exonic sequences
    gr_intrflanks = setdiff(gr_intrflanks, gr_esspl) # Removing any overlaps with essential splice site sequences
    gr_intrflanks = intersect(gr_intrflanks, gr_cds) # Removing UTR sequences
    
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_introns, type="any", select="all"))
    mutmat$wintr[unique(ol[,1])] = 1 # Annotation of the intronic sites
    
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_introns, type="any", select="all"))
    mutmat$wintr[unique(ol[,1])] = 1 # Annotation of the intronic sites
    
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_exonflanks, type="any", select="all"))
    mutmat$wexfl[seq(1,nrow(sites)) %in% ol[,1] & sites$impact=="Synonymous"] = 1 # Annotation of the exonic flank sites (only for synonymous mutations)
    
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_intrflanks, type="any", select="all"))
    mutmat$winfl[unique(ol[,1])] = 1 # Annotation of the exonic flank sites
  }
  
  # Annotation of nonsense mutations in the last coding exon
  
  if (nrow(exons)>1) {
    if (refcdsgene$strand==1) {
      lastexon = exons[nrow(exons),]
    } else {
      lastexon = exons[1,]
    }
    
    gr_lastexon = GenomicRanges::GRanges(gene, IRanges::IRanges(min(lastexon),max(lastexon)))
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_lastexon, type="any", select="all"))
    mutmat$wnonlastex[seq(1,nrow(sites)) %in% ol[,1] & sites$impact=="Nonsense"] = 1
    mutmat$wnon[seq(1,nrow(sites)) %in% ol[,1] & sites$impact=="Nonsense"] = 0 # Removing Nonsense mutations in the last exon from the wnon annotation
  }
  
  ## 3. Other annotations from the user
  
  # Regions defined by chr position
  if (!is.null(regionschr)) {
    
    wnames = unique(regionschr$wname)
    badnames = intersect(wnames,unique(colnames(mutmat)))
    if (length(badnames)>0) { stop(sprintf("The following w names are not allowed in the input regions as they match existing parameters: %s", paste(badnames,collapse = ","))) }
    
    for (j in 1:length(wnames)) {
      
      aux = regionschr[which(regionschr$wname==wnames[j]), ]
      if (length(unique(aux$impacts))!=1) { stop("regionschr: different values found in the impacts column for a given feature, please correct your input object") }
      imps = strsplit(unique(aux$impacts), split=",")[[1]]
      if (any(!(imps %in% setdiff(unique(sites$impact),NA)))) { stop("regionschr: invalid impact terms used in the input object, please correct the impact column") }
      
      gr_wname = GenomicRanges::GRanges(gene, IRanges::IRanges(aux$start,aux$end))
      ol = as.data.frame(GenomicRanges::findOverlaps(gr_sites, gr_wname, type="any", select="all")) #replace gr_sites with gr_sites_old to enable saving of gr_sites output
      mutmat[,wnames[j]] = 0 # Initialising this field in the data frame
      if (length(imps)==0) {
        indx = unique(ol[,1]) # If no impacts are indicated, we consider all sites independent of their impact
      } else {
        indx = intersect(unique(ol[,1]), which(sites$impact %in% imps))
      }
      mutmat[indx,wnames[j]] = 1
      
      # Removing previous w parameters if layered=0
      if (aux$layered[1]==0 | aux$layered[1]==FALSE) {
        mutmat[indx, setdiff(colnames(mutmat),c("t",wnames[j]))] = 0 # Removing previous annotations at these sites
      }
    }
  }
  
  # Regions defined by aminoacid position
  if (!is.null(regionsaa)) {
    
    sites$aux = 1:nrow(sites)
    aas = sites[which(!is.na(sites$aachange) & sites$aachange!="."),]
    aas$aapos = as.numeric(substr(aas$aachange,2,nchar(aas$aachange)-1))
    gr_aas = GenomicRanges::GRanges(gene, IRanges::IRanges(aas$aapos,aas$aapos))
    
    wnames = unique(regionsaa$wname)
    badnames = intersect(wnames,unique(colnames(mutmat)))
    if (length(badnames)>0) { stop(sprintf("The following w names are not allowed in the input regions as they match existing parameters: %s", paste(badnames,collapse = ", "))) }
    
    for (j in 1:length(wnames)) {
      
      aux = regionsaa[which(regionsaa$wname==wnames[j]), ]
      if (length(unique(aux$impacts))!=1) { stop("regionschr: different values found in the impacts column for a given feature, please correct your input object") }
      imps = strsplit(unique(aux$impacts), split=",")[[1]]
      if (any(!(imps %in% setdiff(unique(sites$impact),NA)))) { stop("regionschr: invalid impact terms used in the input object, please correct the impact column") }
      
      gr_wname = GenomicRanges::GRanges(gene, IRanges::IRanges(aux$start,aux$end))
      ol = as.data.frame(GenomicRanges::findOverlaps(gr_aas, gr_wname, type="any", select="all"))
      mutmat[,wnames[j]] = 0 # Initialising this field in the data frame
      if (length(imps)==0) {
        indx = aas$aux[unique(ol[,1])] # If no impacts are indicated, we consider all sites independent of their impact
      } else {
        indx = aas$aux[intersect(unique(ol[,1]), which(aas$impact %in% imps))]
      }
      mutmat[indx,wnames[j]] = 1
      
      # Removing previous w parameters if layered=0
      if (aux$layered[1]==0 | aux$layered[1]==FALSE) {
        mutmat[indx, setdiff(colnames(mutmat),c("t",wnames[j]))] = 0 # Removing previous annotations at these sites
      }
    }
    sites = sites[, setdiff(colnames(sites),"aux")]
  }
  
  ## 4. Poisson regression: fitting the selection model
  
  model = glm(formula = sites$obs ~ offset(log(sites$rnorm)) + . -1, data=mutmat, family=poisson(link=log))
  mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
  pvals = coef(summary(model))[,4]
  model.lrt = drop1(model, test= "Chisq")
  pvals.lrt = setNames(model.lrt[[5]], row.names(model.lrt))
  ci = exp(confint.default(model)) # Wald confidence intervals
  par.pois = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2], pval.wald=pvals[rownames(ci)], pval.lrt=pvals.lrt[rownames(ci)], stringsAsFactors = F) # MLEs and Wald CI95% for the selection parameters
  
  # Adding obs/exp statistics to the regression outputs
  nobs = apply(mutmat, 2, function(x) sum(sites$obs[x==1]))
  nexp = apply(mutmat, 2, function(x) sum(sites$rnorm[x==1]))
  par.pois$obs = nobs[par.pois$name]
  par.pois$exp = nexp[par.pois$name]
  
  ## 5. Negative binomial regression
  
  model.nbfix = model.nb = par.nbfix = par.nb = NULL
  
  if (!is.null(fixtheta)) {
    
    model.nbfix = glm(formula = sites$obs ~ offset(log(sites$rnorm)) + . -1, data=mutmat, family=MASS::negative.binomial(fixtheta))
    mle = exp(coefficients(model.nbfix)) # Maximum-likelihood estimates for the rate params
    pvals = coef(summary(model.nbfix))[,4]
    model.lrt = drop1(model.nbfix, test= "Chisq")
    pvals.lrt = setNames(model.lrt[[5]], row.names(model.lrt))
    ci = exp(confint.default(model)) # Wald confidence intervals
    par.nbfix = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2], pval.wald=pvals[rownames(ci)], pval.lrt=pvals.lrt[rownames(ci)], stringsAsFactors = F) # MLEs and Wald CI95% for the selection parameters
    par.nbfix$obs = nobs[par.nbfix$name]
    par.nbfix$exp = nexp[par.nbfix$name]
    
  } else {
    
    model.nb = MASS::glm.nb(formula = as.vector(sites$obs) ~ offset(log(sites$rnorm)) + . -1, data=mutmat)
    mle = exp(coefficients(model.nb)) # Maximum-likelihood estimates for the rate params
    pvals = coef(summary(model.nb))[,4]
    model.lrt = drop1(model.nb, test= "Chisq")
    pvals.lrt = setNames(model.lrt[[5]], row.names(model.lrt))
    ci = exp(confint.default(model)) # Wald confidence intervals
    par.nb = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2], pval.wald=pvals[rownames(ci)], pval.lrt=pvals.lrt[rownames(ci)], stringsAsFactors = F) # MLEs and Wald CI95% for the selection parameters
    par.nb$obs = nobs[par.nb$name]
    par.nb$exp = nexp[par.nb$name]
  }
  
  ## 6. Outputs
  
  # Annotated sites table
  sites_output <- sites
  mutmat_output <- mutmat
  
  sites2 = cbind(sites, mutmat)
  if (!is.null(sitefilename)) {
    sites2 = apply(sites2,2,as.character)
    write.table(sites2, file=sitefilename, col.names = T, row.names = F, quote = F, sep = "\t")
  }
  
  # Inform the user about the sites used as neutral reference by the model
  
  neutral_sites = (rowSums(mutmat[,setdiff(colnames(mutmat),"t")])==0)
  neutral_nsyn = sum(sites$obs[which(neutral_sites & sites$impact=="Synonymous")]) # Number of synonymous mutations used as background for dN/dS
  neutral_nsynsites = length(which(neutral_sites & sites$impact=="Synonymous"))
  neutral_othersites = length(which(neutral_sites & sites$impact!="Synonymous"))
  nonneutral_nsyn = sum(sites$obs[which(neutral_sites==0 & sites$impact=="Synonymous")]) # Number of synonymous mutations excluded from the neutral backgroup by the "w" annotations
  
  message(sprintf("   Sites used as neutral reference: %0.0f synonymous mutations across %0.0f synonymous sites.", neutral_nsyn, neutral_nsynsites))
  if (neutral_othersites>0) {
    message(sprintf("   %0.0f sites not classified as synonymous coding sites were used in the background model. Please ensure that this was intended.", neutral_othersites))
  }
  if (nonneutral_nsyn>0) {
    message(sprintf("   %0.0f synonymous mutations were not used in the neutral background model. This may be due to excluding known synonymous driver mutations (when using syndrivers), excluding synonymous mutations in exon flanks (when using exon_flank_length>0) and/or due to annotations provided by the user. Please ensure that this is the desired behaviour.", nonneutral_nsyn))
  }
  
  # Output object
  out = list(sites = sites2, dnds.pois = par.pois, dnds.nb = par.nb,  dnds.nbfix = par.nbfix, model.pois = model, model.nb = model.nb, model.nbfix = model.nbfix, sites_output = sites_output, mutmat_output = mutmat_output)
  
}

# Wrapper function containing the withingenednds function, necessary pre- and post-processing, multiple testing correction, and saving of outputs. 
withingenednds_custom <- function(gene, per_gene_cov_dir, mutations, gene_panel, samples, refcds_path, custom_regions = NULL, genomeFile) {
  m = unique(mutations[,c("donor","chr","pos","ref","mut")])
  dndsout = dndscv(mutations = m, gene_list = gene_panel, outmats = TRUE, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, refdb = refcds_path, mingenecovs = 0, onesided = FALSE)
  siteout = sitednds(dndsout = dndsout, gene_list = gene_panel, min_recurr = 1, method = "NB")
  fixtheta = siteout$overdisp[1]
  
  ref_cds_gene_list <- NULL
  
  core_promoter_flank <- 200
  
  load(refcds_path)
  ref_cds_gene_list <- NULL
  for(i in 1:length(RefCDS)){
    ref_cds_gene_list <- c(ref_cds_gene_list,RefCDS[[i]]$gene_name)
  }
  
  gene_idx <- which(ref_cds_gene_list == gene)
  strand <- RefCDS[[gene_idx]]$strand
  
  covtable <- list.files(per_gene_cov_dir, recursive = FALSE, full.names = TRUE, pattern = paste0(gene, ".tsv")) %>%
    lapply(., read_tsv) %>%
    do.call(cbind.data.frame, .) 
  
  sample_cov <- covtable[,colnames(covtable) %in% samples]
  covtable = cbind(covtable[,1:4], dc = rowSums(sample_cov))
  
  regionschr <- as.data.frame(array(data = NA, dim = c(1,6)))
  colnames(regionschr) <- c("chr","start","end","wname","impacts","layered")
  
  if(strand == 1){
    tss_idx <- min(which(covtable$coding == 1 | covtable$coding == 2))
    regionschr$chr <- covtable$chr[tss_idx]
    regionschr$start <- max(covtable$pos[tss_idx] - core_promoter_flank, covtable$pos[1])
    regionschr$end <- min(covtable$pos[tss_idx] + (core_promoter_flank - 1), covtable$pos[min(which(covtable$coding == 1))] - 1)
    regionschr$wname <- "wcorepromoter"
    regionschr$impacts <- ""
    regionschr$layered <- 0
    
    regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],covtable$pos[1],covtable$pos[tss_idx] - (core_promoter_flank + 1),"wupstreamcore","",0))
    
    upstream_utr <- covtable[which(covtable$coding == 2 & covtable$pos < covtable$pos[min(which(covtable$coding == 1))] & covtable$pos > min(covtable$pos[tss_idx] + (core_promoter_flank - 1), covtable$pos[min(which(covtable$coding == 1))])),]
    
    if(nrow(upstream_utr) > 1){
      upstream_utr$region <- NA
      for(j in 1:nrow(upstream_utr)){
        if(j == 1){
          upstream_utr$region[j] <- 1
        }else if(upstream_utr$pos[j] == upstream_utr$pos[j - 1] + 1){
          upstream_utr$region[j] <- upstream_utr$region[j - 1]
        }else{
          upstream_utr$region[j] <- upstream_utr$region[j - 1] + 1
        }
      }
      
      for(j in 1:length(unique(upstream_utr$region))){
        regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],min(upstream_utr$pos[which(upstream_utr$region == j)]),max(upstream_utr$pos[which(upstream_utr$region == j)]),"wupstreamUTR","",0))
      }
    }
    
    downstream_utr <- covtable[which(covtable$coding == 2 & covtable$pos > covtable$pos[max(which(covtable$coding == 1))]),]
    
    if(nrow(downstream_utr) > 1){
      downstream_utr$region <- NA
      for(j in 1:nrow(downstream_utr)){
        if(j == 1){
          downstream_utr$region[j] <- 1
        }else if(downstream_utr$pos[j] == downstream_utr$pos[j - 1] + 1){
          downstream_utr$region[j] <- downstream_utr$region[j - 1]
        }else{
          downstream_utr$region[j] <- downstream_utr$region[j - 1] + 1
        }
      }
      
      for(j in 1:length(unique(downstream_utr$region))){
        regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],min(downstream_utr$pos[which(downstream_utr$region == j)]),max(downstream_utr$pos[which(downstream_utr$region == j)]),"wdownstreamUTR","",0))
      }
    }
    
    regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],max(covtable$pos[which(covtable$coding %in% c(1,2))]) + 1,max(covtable$pos),"wdownstreampolyA","",0))
    
  }else if(strand == -1){
    tss_idx <- max(which(covtable$coding == 1 | covtable$coding == 2))
    regionschr$chr <- covtable$chr[tss_idx]
    regionschr$start <- max(covtable$pos[tss_idx] - (core_promoter_flank - 1), covtable$pos[max(which(covtable$coding == 1))] + 1)
    regionschr$end <- min(covtable$pos[tss_idx] + core_promoter_flank, covtable$pos[nrow(covtable)])
    regionschr$wname <- "wcorepromoter"
    regionschr$impacts <- ""
    regionschr$layered <- 0
    
    regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],covtable$pos[tss_idx] + (core_promoter_flank + 1),covtable$pos[nrow(covtable)],"wupstreamcore","",0))
    
    upstream_utr <- covtable[which(covtable$coding == 2 & covtable$pos > covtable$pos[max(which(covtable$coding == 1))] & covtable$pos < max(covtable$pos[tss_idx] - (core_promoter_flank - 1), covtable$pos[max(which(covtable$coding == 1))])),]
    
    if(nrow(upstream_utr) > 1){
      upstream_utr$region <- NA
      for(j in 1:nrow(upstream_utr)){
        if(j == 1){
          upstream_utr$region[j] <- 1
        }else if(upstream_utr$pos[j] == upstream_utr$pos[j - 1] + 1){
          upstream_utr$region[j] <- upstream_utr$region[j - 1]
        }else{
          upstream_utr$region[j] <- upstream_utr$region[j - 1] + 1
        }
      }
      
      for(j in 1:length(unique(upstream_utr$region))){
        regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],min(upstream_utr$pos[which(upstream_utr$region == j)]),max(upstream_utr$pos[which(upstream_utr$region == j)]),"wupstreamUTR","",0))
      }
    }
    
    downstream_utr <- covtable[which(covtable$coding == 2 & covtable$pos < covtable$pos[min(which(covtable$coding == 1))]),]
    
    if(nrow(downstream_utr) > 1){
      downstream_utr$region <- NA
      for(j in 1:nrow(downstream_utr)){
        if(j == 1){
          downstream_utr$region[j] <- 1
        }else if(downstream_utr$pos[j] == downstream_utr$pos[j - 1] + 1){
          downstream_utr$region[j] <- downstream_utr$region[j - 1]
        }else{
          downstream_utr$region[j] <- downstream_utr$region[j - 1] + 1
        }
      }
      
      for(j in 1:length(unique(downstream_utr$region))){
        regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],min(downstream_utr$pos[which(downstream_utr$region == j)]),max(downstream_utr$pos[which(downstream_utr$region == j)]),"wdownstreamUTR","",0))
      }
    }
    
    regionschr <- rbind(regionschr, c(covtable$chr[tss_idx],min(covtable$pos),min(covtable$pos[which(covtable$coding %in% c(1,2))]) - 1,"wdownstreampolyA","",0))
    
  }
  
  if(!(is.null(custom_regions))){
    regionschr <- rbind(regionschr, custom_regions)
  }
  
  regionschr$start <- as.numeric(regionschr$start)
  regionschr$end <- as.numeric(regionschr$end)
  
  # Save pre-defined regions to test
  if (!dir.exists(paste0("./output/withingenednds"))) {dir.create(paste0("./output/withingenednds"))}
  if (!dir.exists(paste0("./output/withingenednds/",gene))) {dir.create(paste0("./output/withingenednds/",gene))}
  if (!dir.exists(paste0("./output/withingenednds/",gene,"regions_chr/",gene))) {dir.create(paste0("./output/withingenednds/",gene,"/regions_chr/"))}
  
  write.table(x = regionschr, file = paste0("./output/withingenednds/",gene,"/regions_chr/",gene,"_regions_chr.tsv"), sep = "\t", row.names = F, col.names = T)
  
  # Run withingenednds
  genedndsout = withingenednds(mutations=m, gene=gene, syndrivers = NULL, covtable=covtable, dndsout=dndsout, genomeFile=genomeFile, regionschr = regionschr, fixtheta = fixtheta, normalisefromsyn = TRUE, exon_flank_length = 10, intron_flank_length = 10, refdb = refcds_path, numcode = 1)
  
  # Multiple testing correction (Andrew 09/01/2025)
  pvalcols <- which(grepl("pval.wald",colnames(genedndsout$dnds.nbfix))) # pval columns
  qnames <- gsub("pval","qval",colnames(genedndsout$dnds.nbfix)[pvalcols]) # Names for the qval columns
  for(j in 1:length(qnames)) {
    genedndsout$dnds.nbfix[,qnames[j]] = p.adjust(genedndsout$dnds.nbfix[,pvalcols[j]], method="BH")
  }
  
  pvalcols <- which(grepl("pval.lrt",colnames(genedndsout$dnds.nbfix))) # pval columns
  qnames <- gsub("pval","qval",colnames(genedndsout$dnds.nbfix)[pvalcols]) # Names for the qval columns
  for(j in 1:length(qnames)) {
    genedndsout$dnds.nbfix[,qnames[j]] = p.adjust(genedndsout$dnds.nbfix[,pvalcols[j]], method="BH")
  }
  
  nbfix_gene <<- genedndsout$dnds.nbfix
  
  # Save outputs
  if (!dir.exists(paste0("./output/withingenednds/",gene,"pre_annotation_dNdS_output/",gene))) {dir.create(paste0("./output/withingenednds/",gene,"/pre_annotation_dNdS_output/"))}
  
  write.table(x = genedndsout$dnds.pois, file = paste0("./output/withingenednds/",gene,"/pre_annotation_dNdS_output/",gene,"_dnds_pois.tsv"), sep = "\t", row.names = F, col.names = T)
  write.table(x = genedndsout$dnds.nbfix, file = paste0("./output/withingenednds/",gene,"/pre_annotation_dNdS_output/",gene,"_dnds_nbfix.tsv"), sep = "\t", row.names = F, col.names = T)
  saveRDS(object = genedndsout$model.pois, file = paste0("./output/withingenednds/",gene,"/pre_annotation_dNdS_output/",gene,"_model_pois.rds"))
  saveRDS(object = genedndsout$model.nbfix, file = paste0("./output/withingenednds/",gene,"/pre_annotation_dNdS_output/",gene,"_model_nbfix.rds"))
  
  write.table(x = genedndsout$sites, file = paste0("./output/withingenednds/",gene,"/pre_annotation_dNdS_output/",gene,"_withingene_sites.tsv"), sep = "\t", col.names = T, row.names = F)
}

# Withingenednds plotter
nbfix_plotter <- function(nbfix, gene, runman, fig.width, fig.height) {
  pal_sig <- brewer.pal(n = 5, name = "YlOrRd")
  nbfix = nbfix[order(nbfix$name),]
  
  y_limit <- ceiling(max(nbfix$cihigh)/10) * 10
  y_limit_ci <- ceiling((max(nbfix$cihigh) + (y_limit / (2 * nrow(nbfix))))/10) * 10
  
  if (runman) { dev.new(width=6, height=4) }
  plot(as.numeric(nbfix$name), nbfix$mle, pch = 19, xlim = c(min(as.numeric(nbfix$name)) - 0.5,max(as.numeric(nbfix$name)) + 0.5), ylim = c(0,y_limit_ci), xaxt="n", yaxt="n", xlab = NA, ylab = NA)
  axis(side = 1, at=as.numeric(nbfix$name), labels=FALSE, tick = FALSE)
  text(x = as.numeric(nbfix$name), y = par("usr")[3] - 1.5, labels = nbfix$name, xpd = NA, srt = 45, adj = 1)
  axis(side = 2, at=seq(0,y_limit_ci,10), labels=seq(0,y_limit_ci,10), las=2)
  abline(h = 1, col = "red", lty = "dashed")
  arrows(x0 = as.numeric(nbfix$name), x1 = as.numeric(nbfix$name), y0 = nbfix$cilow, y1 = nbfix$cihigh, length=0.05, angle=90, code=3)
  if(length(which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp))) > 0){
    rect(xleft = which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), xright = which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp)) + (1/nrow(nbfix)), ybottom = nbfix$cihigh[which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), ytop = nbfix$cihigh[which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp))] + (y_limit / (2 * nrow(nbfix))), border = "white", col=pal_sig[5])
    text(y = rep(nbfix$cihigh[which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), each = 4) + rep(c(1/3,2/3) * ((y_limit / (2 * nrow(nbfix))) - (y_limit / (4 * nrow(nbfix)))), each = 2), x = rep(which(nbfix$qval.lrt <= 0.0001 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), each = 4) + rep(c(1/3,2/3) * (2/nrow(nbfix)), times = 2), labels = "*")
  }
  if(length(which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp))) > 0){
    rect(xleft = which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), xright = which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp)) + (1/nrow(nbfix)), ybottom = nbfix$cihigh[which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), ytop = nbfix$cihigh[which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp))] + (y_limit / (2 * nrow(nbfix))), border = "white", col=pal_sig[4])
    text(y = rep(nbfix$cihigh[which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), each = 3) + c(1/3,2/3,1/3) * ((y_limit / (2 * nrow(nbfix))) - (y_limit / (4 * nrow(nbfix)))), x = rep(which(nbfix$qval.lrt <= 0.001 & nbfix$qval.lrt > 0.0001 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), each = 3) + c(1/3,1/2,2/3) * (2/nrow(nbfix)), labels = "*")
  }
  if(length(which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp))) > 0){
    rect(xleft = which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), xright = which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp)) + (1/nrow(nbfix)), ybottom = nbfix$cihigh[which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), ytop = nbfix$cihigh[which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp))] + (y_limit / (2 * nrow(nbfix))), border = "white", col=pal_sig[3])
    text(y = rep(nbfix$cihigh[which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))) + 0.5 * ((y_limit / (2 * nrow(nbfix))) - (y_limit / (4 * nrow(nbfix)))), each = 2), x = rep(which(nbfix$qval.lrt <= 0.01 & nbfix$qval.lrt > 0.001 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), each = 2) + rep(c(1/3,2/3) * (2/nrow(nbfix)), times = 2), labels = "*")
  }
  if(length(which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp))) > 0){
    rect(xleft = which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), xright = which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp)) + (1/nrow(nbfix)), ybottom = nbfix$cihigh[which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), ytop = nbfix$cihigh[which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp))] + (y_limit / (2 * nrow(nbfix))), border = "white", col=pal_sig[2])
    text(y = nbfix$cihigh[which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))) + 0.5 * ((y_limit / (2 * nrow(nbfix))) - (y_limit / (4 * nrow(nbfix)))), x = which(nbfix$qval.lrt <= 0.1 & nbfix$qval.lrt > 0.01 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)) + (1/2 * (2/nrow(nbfix))), labels = "*")
  }
  if(length(which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp))) > 0){
    rect(xleft = which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)), xright = which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp)) + (1/nrow(nbfix)), ybottom = nbfix$cihigh[which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))), ytop = nbfix$cihigh[which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp))] + (y_limit / (2 * nrow(nbfix))), border = "white", col=pal_sig[1])
    text(y = nbfix$cihigh[which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp))] + (y_limit / (4 * nrow(nbfix))) + 0.5 * ((y_limit / (2 * nrow(nbfix))) - (y_limit / (4 * nrow(nbfix)))), x = which(nbfix$pval.lrt <= 0.01 & nbfix$qval.lrt > 0.1 & (nbfix$obs > nbfix$exp)) - (1/nrow(nbfix)) + (1/2 * (2/nrow(nbfix))), labels = ".")
  }
  
  title(main = gene, ylab = "dN/dS ratio (95% CIs)", xlab = "")
  if (runman) { dev.copy(pdf,paste0("./output/", gene, "_nbfix_plot.pdf"),width=fig.width,height=fig.height); dev.off() }
}

# Function that spatially maps the mutations called in NanoSeq data (and piled up in LCM-WES data) to the microdissected tissue sample
spatial_mapper <- function(LCM_record_file, muts_of_interest, gene, pileups_dir, donor, outline_file, cuts_file, show_plot = F, restrict_to_lymphocytes = T) {
  # Import the outline annotation
  biopsy_xml <- XML::xmlParse(file = outline_file)
  biopsy_list <- XML::xmlToList(biopsy_xml)
  
  # Import the LCM cuts annotations
  cuts_xml <- XML::xmlParse(file = cuts_file)
  cuts_list <- XML::xmlToList(cuts_xml)
  
  outline_array <- as.data.frame(array(data = NA, dim = c(length(biopsy_list[[1]]$annotation$pointlist),4)))
  colnames(outline_array) <- c("Cut_ID","PD_ID","x","y")
  
  # Import your LCM data that indicates what each cut maps to. Ensure the columns Cut_ID (lymphocyes_1), donor (PD12345), PD_ID (PD12345_lo0001), Proceed (Y or N) are present
  ID_map <- readxl::read_xlsx(path = LCM_record_file, sheet = "LCM_thyroid")
  ID_map <- ID_map %>% dplyr::rename(Cut_ID = Histology_cut, donor = PD_ID, PD_ID = Well_ID, Proceed = Proceed_WES)
  ID_map <- ID_map[which(ID_map$donor == donor),]
  ID_map$Cut_ID <- toupper(ID_map$Cut_ID)
  
  # Subset the matrix generated, to the available pileup files present
  files <- list.files(paste0(pileups_dir,gene,"/"), recursive = T)
  files <- gsub(".*logs/", "", files)
  files <- gsub("\\.bam.*", "", files)
  files <- gsub(".*/", "", files)
  files <- gsub("\\_sorted.*", "", files)
  files <- gsub("\\_count.*", "", files)
  files <- unique(files[!grepl("tsv", files)])
  files <- unique(files[!grepl("NA", files)])
  
  ID_map <- ID_map[which(ID_map$PD_ID %in% files),]
  
  # Ensure all microdissections labelled appropriately (as LYMPHOCYTES) and subset to these
  ID_map$Cut_ID <- gsub("LYMPHOID_AGGREGATE", "LYMPHOCYTES_", ID_map$Cut_ID)
  ID_map$Cut_ID <- gsub("__", "_", ID_map$Cut_ID)
  if(restrict_to_lymphocytes){
    ID_map <- ID_map[grepl("LYMPHOCYTES", ID_map$Cut_ID), , drop = FALSE]
  }
  
  # Extract x and y coordinates for the outline of tissue section
  outline_array$Cut_ID <- biopsy_list[[1]]$title
  
  for(i in 1:nrow(outline_array)){
    outline_array[i,"x"] <- biopsy_list[[1]]$annotation$pointlist[[i]]$x
    outline_array[i,"y"] <- biopsy_list[[1]]$annotation$pointlist[[i]]$y
  }
  
  outline_array$x <- as.numeric(outline_array$x)
  outline_array$y <- as.numeric(outline_array$y)
  
  # Depending on the shape you have drawn with on the LCM annotations file (fixed circles or freehand etc). May need to include more types if you have drawn with squares or others.
  freehand_cut_array <- NULL 
  circle_cut_array <- NULL
  
  # Extract the cuts coordinates
  for(i in 1:length(cuts_list)){
    if(cuts_list[[i]]$annotation$.attrs[1] == "freehand"){
      temp_cut_array <- as.data.frame(array(data = NA, dim = c(length(cuts_list[[i]]$annotation$pointlist),4)))
      if(nrow(temp_cut_array) >= 1){
        colnames(temp_cut_array) <- c("Cut_ID","PD_ID","x","y")
        temp_cut_array$Cut_ID <- toupper(cuts_list[[i]]$title)
        for(j in 1:nrow(temp_cut_array)){
          temp_cut_array[j,"x"] <- cuts_list[[i]]$annotation$pointlist[[j]]$x
          temp_cut_array[j,"y"] <- cuts_list[[i]]$annotation$pointlist[[j]]$y
        }
        if(is.null(freehand_cut_array)){
          freehand_cut_array <- temp_cut_array
        }else{
          freehand_cut_array <- rbind(freehand_cut_array,temp_cut_array)
        }
      }
    }else if(cuts_list[[i]]$annotation$.attrs[1] == "circle"){
      temp_cut_array <- as.data.frame(array(data = NA, dim = c(1,5)))
      colnames(temp_cut_array) <- c("Cut_ID","PD_ID","x","y","radius")
      temp_cut_array$Cut_ID <- toupper(cuts_list[[i]]$title)
      temp_cut_array$x <- cuts_list[[i]]$annotation$x
      temp_cut_array$y <- cuts_list[[i]]$annotation$y
      temp_cut_array$radius <- cuts_list[[i]]$annotation$radius
      if(is.null(circle_cut_array)){
        circle_cut_array <- temp_cut_array
      }else{
        circle_cut_array <- rbind(circle_cut_array,temp_cut_array)
      }
    }
  }
  
  freehand_cut_array$x <- as.numeric(freehand_cut_array$x)
  freehand_cut_array$y <- as.numeric(freehand_cut_array$y)
  freehand_cut_array <- as.data.frame(freehand_cut_array)
  
  circle_cut_array$x <- as.numeric(circle_cut_array$x)
  circle_cut_array$y <- as.numeric(circle_cut_array$y)
  circle_cut_array$radius <- as.numeric(circle_cut_array$radius)
  circle_cut_array <- as.data.frame(circle_cut_array)
  
  if(!(is.null(nrow(freehand_cut_array))) & nrow(freehand_cut_array) > 0){
    freehand_cut_array$Cut_ID <- gsub("LYMPHOID_AGGREGATE", "LYMPHOCYTES_", freehand_cut_array$Cut_ID)
    freehand_cut_array$Cut_ID <- gsub("__", "_", freehand_cut_array$Cut_ID)
  }
  
  if(!(is.null(nrow(circle_cut_array))) & nrow(circle_cut_array) > 0){
    circle_cut_array$Cut_ID <- gsub("LYMPHOID_AGGREGATE", "LYMPHOCYTES_", circle_cut_array$Cut_ID)
    circle_cut_array$Cut_ID <- gsub("__", "_", circle_cut_array$Cut_ID)
  }
  
  
  # Annotate the cuts and their coordinates
  for(i in 1:nrow(ID_map)){
    freehand_cut_array$PD_ID[which(freehand_cut_array$Cut_ID == ID_map$Cut_ID[i])] <- ID_map$PD_ID[i]
  }
  
  for(i in 1:nrow(ID_map)){
    circle_cut_array$PD_ID[which(circle_cut_array$Cut_ID == ID_map$Cut_ID[i])] <- ID_map$PD_ID[i]
  }
  
  # Remove NA
  if(length(freehand_cut_array$PD_ID) > 0){
    freehand_cut_array <- freehand_cut_array[!is.na(freehand_cut_array$PD_ID),]
  }
  if(length(circle_cut_array$PD_ID) > 0){
    circle_cut_array <- circle_cut_array[!is.na(circle_cut_array$PD_ID),]
  }
  
  # Select genes of interest and add their coordinates (matching the pileup data)
  gene_pos <- as.data.frame(array(data = NA, dim = c(6,3)))
  colnames(gene_pos) <- c("gene","start","end")
  
  gene_pos[1,] <- c("TNFRSF14",2486578,2497321)
  gene_pos[2,] <- c("CD274",5450003,5471066)
  gene_pos[3,] <- c("CCR6",167525295,167553184)
  gene_pos[4,] <- c("TET2",106067032,106200973)
  gene_pos[5,] <- c("CBL",119076252,119179359)
  gene_pos[6,] <- c("DNMT3A",25455845,25565459)
  
  gene_pos$start <- as.numeric(gene_pos$start)
  gene_pos$end <- as.numeric(gene_pos$end)
  
  muts_of_interest <- muts_of_interest[,c("sampleID","gene","chr","pos","ref","mut","aachange")][which(substr(muts_of_interest$sampleID,1,7) %in% donor & muts_of_interest$gene == gene),]
  
  if (nrow(muts_of_interest) == 0) {
    stop(paste("Spatial mapping stopped as there are no", gene, "mutations to plot", sep = " "))
  }
  
  muts_of_interest <- muts_of_interest[ , !names(muts_of_interest) %in% c("sampleID")] 
  
  # Identify indels and re-code them for pileup
  if (any((nchar(muts_of_interest$mut) - nchar(muts_of_interest$ref)) >= 1)) {
    muts_of_interest[which((nchar(muts_of_interest$mut) - nchar(muts_of_interest$ref)) >= 1),]$mut <- "INS"
  }
  
  if(any((nchar(muts_of_interest$ref) - nchar(muts_of_interest$mut)) >= 1)) {
    muts_of_interest[which((nchar(muts_of_interest$ref) - nchar(muts_of_interest$mut)) >= 1),]$mut <- "X."
  }
  
  del_rows <- which(muts_of_interest$mut == "X.")
  if(length(del_rows) != 0){
    del_old <- muts_of_interest[del_rows,]
    muts_of_interest <- muts_of_interest[-(del_rows),]
    del_new <- as.data.frame(array(data = NA, dim = c(0,ncol(del_old))))
    colnames(del_new) <- colnames(del_old)
    for(x in 1:length(del_rows)){
      del_old$ref[x] <- substr(del_old$ref[x],2,nchar(del_old$ref[x]))
      del_old$pos[x] <- del_old$pos[x] + 1
      if(nchar(del_old$ref[x]) == 1){
        del_new <- rbind(del_new, del_old[x,])
      }else{
        del_tmp <- do.call("rbind", replicate(
          n = nchar(del_old$ref[x]), del_old[x,], simplify = FALSE))
        for(y in 1:nrow(del_tmp)){
          del_tmp$pos[y] <- del_tmp$pos[y] + y - 1
          del_tmp$ref[y] <- substr(del_tmp$ref[y],y,y)
        }
        del_new <- rbind(del_new, del_tmp)
      }
    }
    muts_of_interest <- rbind(muts_of_interest, del_new)
  }
  
  # Identify DNVs and split them up for pileup
  DNVs <- bind_rows(muts_of_interest[which((nchar(muts_of_interest$mut) - nchar(muts_of_interest$ref)) == 0 & nchar(muts_of_interest$ref) == 2 & nchar(muts_of_interest$mut) == 2),], muts_of_interest[which((nchar(muts_of_interest$mut) - nchar(muts_of_interest$ref)) == 0 & nchar(muts_of_interest$ref) == 2 & nchar(muts_of_interest$mut) == 2),])
  DNVs_odds <- DNVs[seq_len(nrow(DNVs))%%2 == 1,]
  DNVs_odds$ref <- substr(DNVs_odds$ref,1,1)
  DNVs_odds$mut <- substr(DNVs_odds$mut,1,1)
  
  DNVs_evens <- DNVs[seq_len(nrow(DNVs))%%2 == 0,]
  DNVs_evens$ref <- substr(DNVs_evens$ref,2,2)
  DNVs_evens$mut <- substr(DNVs_evens$mut,2,2)
  DNVs_evens$pos <- DNVs_evens$pos + 1
  
  muts_of_interest <- dplyr::anti_join(muts_of_interest,DNVs)
  muts_of_interest <- dplyr::bind_rows(muts_of_interest, bind_rows(DNVs_odds, DNVs_evens))
  
  muts_of_interest$pos <- as.numeric(muts_of_interest$pos)
  
  # Remove duplicates introduced through the splitting of DNVs
  muts_of_interest <- muts_of_interest[which(!(duplicated(paste(muts_of_interest$chr,muts_of_interest$pos,muts_of_interest$ref, muts_of_interest$mut, sep = "_")))),]
  muts_of_interest <- muts_of_interest[order(muts_of_interest$pos, decreasing = F),]
  
  # Generate new matrix and columns
  circle_cols <- ncol(circle_cut_array)
  freehand_cols <- ncol(freehand_cut_array)
  
  if (!(is.null(nrow(circle_cut_array))) & nrow(circle_cut_array) > 0) {
    circle_cut_array[,c((circle_cols + 1):(circle_cols + (5 * nrow(muts_of_interest))))] <- NA
    for(i in 1:nrow(muts_of_interest)){
      colnames(circle_cut_array)[c(((5 * (i - 1)) + 1 + circle_cols):((5 * i) + circle_cols))] <- paste0(rep(paste(muts_of_interest[i,],collapse = "_"), times = 5),"_",c("total_dep","ref_count","mut_count","mut_vaf","mut_vaf_binned"))
    }
  }
  
  if (!(is.null(nrow(freehand_cut_array))) & nrow(freehand_cut_array) > 0) {
    freehand_cut_array[,c((freehand_cols + 1):(freehand_cols + (5 * nrow(muts_of_interest))))] <- NA
    for(i in 1:nrow(muts_of_interest)){
      colnames(freehand_cut_array)[c(((5 * (i - 1)) + 1 + freehand_cols):((5 * i) + freehand_cols))] <- paste0(rep(paste(muts_of_interest[i,],collapse = "_"), times = 5),"_",c("total_dep","ref_count","mut_count","mut_vaf","mut_vaf_binned"))
    }
  }
  
  # Read in pileup data
  if (!(is.null(nrow(freehand_cut_array))) & nrow(freehand_cut_array) > 0) {
    for(i in 1:nrow(ID_map)){
      if(ID_map$Proceed[i] == "Y"){
        for(j in 1:nrow(muts_of_interest)){
          
          temp_pileup <- read.csv(paste0(pileups_dir, gene, "/", ID_map$PD_ID[i],"_count.csv"), sep = ",", header = T)
          temp_pileup <- temp_pileup[(muts_of_interest$pos[j] - gene_pos$start[which(gene_pos$gene == muts_of_interest$gene[j])] + 1),]
          temp_pileup <- temp_pileup[,c(1:6)] + temp_pileup[c(7:12)]
          temp_vaf <- temp_pileup[1,c(which(colnames(temp_pileup) == muts_of_interest$mut[j]))] / sum(temp_pileup[which(!(names(temp_pileup) %in% c("INS","ins")))])
          
          freehand_cut_array[which(freehand_cut_array$PD_ID == ID_map$PD_ID[i]),(freehand_cols + (5 * (j - 1)) + 1)] <- sum(temp_pileup)
          freehand_cut_array[which(freehand_cut_array$PD_ID == ID_map$PD_ID[i]),(freehand_cols + (5 * (j - 1)) + 2)] <- temp_pileup[1,c(which(colnames(temp_pileup) == muts_of_interest$ref[j]))]
          freehand_cut_array[which(freehand_cut_array$PD_ID == ID_map$PD_ID[i]),(freehand_cols + (5 * (j - 1)) + 3)] <- temp_pileup[1,c(which(colnames(temp_pileup) == muts_of_interest$mut[j]))]
          freehand_cut_array[which(freehand_cut_array$PD_ID == ID_map$PD_ID[i]),(freehand_cols + (5 * (j - 1)) + 4)] <- temp_vaf
        }
      }
    }
  }
  
  if (!(is.null(nrow(circle_cut_array))) & nrow(circle_cut_array) > 0) {
    for(i in 1:nrow(ID_map)){
      if(ID_map$Proceed[i] == "Y"){
        for(j in 1:nrow(muts_of_interest)){
          
          temp_pileup <- read.csv(paste0(pileups_dir, gene, "/", ID_map$PD_ID[i],"_count.csv"), sep = ",", header = T)
          temp_pileup <- temp_pileup[(muts_of_interest$pos[j] - gene_pos$start[which(gene_pos$gene == muts_of_interest$gene[j])] + 1),]
          temp_pileup <- temp_pileup[,c(1:6)] + temp_pileup[c(7:12)]
          temp_vaf <- temp_pileup[1,c(which(colnames(temp_pileup) == muts_of_interest$mut[j]))]/ sum(temp_pileup)
          
          circle_cut_array[which(circle_cut_array$PD_ID == ID_map$PD_ID[i]),(circle_cols + (5 * (j - 1)) + 1)] <- sum(temp_pileup)
          circle_cut_array[which(circle_cut_array$PD_ID == ID_map$PD_ID[i]),(circle_cols + (5 * (j - 1)) + 2)] <- temp_pileup[1,c(which(colnames(temp_pileup) == muts_of_interest$ref[j]))]
          circle_cut_array[which(circle_cut_array$PD_ID == ID_map$PD_ID[i]),(circle_cols + (5 * (j - 1)) + 3)] <- temp_pileup[1,c(which(colnames(temp_pileup) == muts_of_interest$mut[j]))]
          circle_cut_array[which(circle_cut_array$PD_ID == ID_map$PD_ID[i]),(circle_cols + (5 * (j - 1)) + 4)] <- temp_vaf
        }
      }
    }
  }
  
  # For binned values generate colour scale
  if (!(is.null(nrow(freehand_cut_array))) & nrow(freehand_cut_array) > 0) {
    for(j in 1:nrow(muts_of_interest)){
      freehand_cut_array[,(freehand_cols + (5 * (j - 1)) + 5)] <- cut(freehand_cut_array[,(freehand_cols + (5 * (j - 1)) + 4)], breaks = c(-1, 0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1), labels = c("0","<0.01","0.01-0.025","0.025-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",">0.5"))
    }
  }
  
  if (!(is.null(nrow(circle_cut_array))) & nrow(circle_cut_array) > 0) {
    for(j in 1:nrow(muts_of_interest)){
      circle_cut_array[,(circle_cols + (5 * (j - 1)) + 5)] <- cut(circle_cut_array[,(circle_cols + (5 * (j - 1)) + 4)], breaks = c(-1, 0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1), labels = c("0","<0.01","0.01-0.025","0.025-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",">0.5"))
    }
  }
  
  color_panel <- brewer.pal(n = 9, name = "YlOrRd")
  
  # Remove NA
  freehand_cut_array[is.na(freehand_cut_array)] <- 0
  circle_cut_array[is.na(circle_cut_array)] <- 0
  
  # Generate plots
  if (length(freehand_cut_array) > 4 & length(circle_cut_array) > 4) {
    for(i in 1:nrow(muts_of_interest)){
      if(restrict_to_lymphocytes){
        spatial_plot <- ggplot() +
          geom_polygon(data = freehand_cut_array[which((grepl(pattern = "LYMPHOID_AGGREGATE", x = freehand_cut_array$Cut_ID) | grepl(pattern = "LYMPHOCYTES", x = freehand_cut_array$Cut_ID) | grepl(pattern = "LYMPHOEPITHELIAL", x = freehand_cut_array$Cut_ID)) & !(freehand_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(freehand_cut_array$PD_ID))),], aes(x = x, y = y, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_circle(data = circle_cut_array[which((grepl(pattern = "LYMPHOCYTES", x = circle_cut_array$Cut_ID) | grepl(pattern = "LYMPHOEPITHELIAL", x = circle_cut_array$Cut_ID)) & !(circle_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(circle_cut_array$PD_ID))),], aes(x0 = x, y0 = y, r = radius, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_polygon(data = outline_array, aes(x = x, y = y, group = Cut_ID), fill = NA, color = "black", linewidth = 1.5) +
          # scale_fill_manual(values=c(">0.3"=color_panel[9], "0.25-0.3"=color_panel[8], "0.2-0.25"=color_panel[7], "0.15-0.2"=color_panel[6], "0.1-0.15"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black") +
          scale_fill_manual(values=c(">0.5"=color_panel[9], "0.4-0.5"=color_panel[8], "0.3-0.4"=color_panel[7], "0.2-0.3"=color_panel[6], "0.1-0.2"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black", drop = FALSE) +
          theme_bw() +
          labs(x = "", y = "", fill = "VAF", title = paste0(muts_of_interest[i,1], " g.", muts_of_interest[i,2], ":", muts_of_interest[i,3], muts_of_interest[i,4], ">", muts_of_interest[i,5], " ", muts_of_interest[i,6])) +
          # theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11.5), legend.text = element_text(size = 11.5), legend.position = "bottom") +
          guides(fill=guide_legend(nrow=2, byrow=TRUE))
        ggsave(plot = spatial_plot, filename = paste0(paste0("./output/spatial_mapping/", donor, "/"), paste(muts_of_interest[i,],collapse = "_"),"_VAF_binned_lymphocytes_only.pdf"), width = (10 * ((1.1 * (max(outline_array$x) - min(outline_array$x))) / (1.1 * (max(outline_array$y) - min(outline_array$y))))), height = 10, create.dir = TRUE)
      }else{
        spatial_plot <- ggplot() +
          geom_polygon(data = freehand_cut_array[which(!(freehand_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(freehand_cut_array$PD_ID))),], aes(x = x, y = y, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_circle(data = circle_cut_array[which(!(circle_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(circle_cut_array$PD_ID))),], aes(x0 = x, y0 = y, r = radius, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_polygon(data = outline_array, aes(x = x, y = y, group = Cut_ID), fill = NA, color = "black", linewidth = 1.5) +
          # scale_fill_manual(values=c(">0.3"=color_panel[9], "0.25-0.3"=color_panel[8], "0.2-0.25"=color_panel[7], "0.15-0.2"=color_panel[6], "0.1-0.15"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black") +
          scale_fill_manual(values=c(">0.5"=color_panel[9], "0.4-0.5"=color_panel[8], "0.3-0.4"=color_panel[7], "0.2-0.3"=color_panel[6], "0.1-0.2"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black", drop = FALSE) +
          theme_bw() +
          labs(x = "", y = "", fill = "VAF", title = paste0(muts_of_interest[i,1], " g.", muts_of_interest[i,2], ":", muts_of_interest[i,3], muts_of_interest[i,4], ">", muts_of_interest[i,5], " ", muts_of_interest[i,6])) +
          # theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11.5), legend.text = element_text(size = 11.5), legend.position = "bottom") +
          guides(fill=guide_legend(nrow=2, byrow=TRUE))
        ggsave(plot = spatial_plot, filename = paste0(paste0("./output/spatial_mapping/", donor, "/"), paste(muts_of_interest[i,],collapse = "_"),"_VAF_binned_lymphocytes_only.pdf"), width = (10 * ((1.1 * (max(outline_array$x) - min(outline_array$x))) / (1.1 * (max(outline_array$y) - min(outline_array$y))))), height = 10, create.dir = TRUE)
        
      }
      if (show_plot) {print(spatial_plot)}
    }
  } else if (length(freehand_cut_array) <= 4 & length(circle_cut_array) > 4) {
    for(i in 1:nrow(muts_of_interest)){
      if(restrict_to_lymphocytes){
        spatial_plot <- ggplot() +
          # geom_polygon(data = freehand_cut_array[which((grepl(pattern = "LYMPHOID_AGGREGATE", x = freehand_cut_array$Cut_ID) | grepl(pattern = "LYMPHOCYTES", x = freehand_cut_array$Cut_ID) | grepl(pattern = "LYMPHOEPITHELIAL", x = freehand_cut_array$Cut_ID)) & !(freehand_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(freehand_cut_array$PD_ID))),], aes(x = x, y = y, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_circle(data = circle_cut_array[which((grepl(pattern = "LYMPHOCYTES", x = circle_cut_array$Cut_ID) | grepl(pattern = "LYMPHOEPITHELIAL", x = circle_cut_array$Cut_ID)) & !(circle_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(circle_cut_array$PD_ID))),], aes(x0 = x, y0 = y, r = radius, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_polygon(data = outline_array, aes(x = x, y = y, group = Cut_ID), fill = NA, color = "black", linewidth = 1.5) +
          # scale_fill_manual(values=c(">0.3"=color_panel[9], "0.25-0.3"=color_panel[8], "0.2-0.25"=color_panel[7], "0.15-0.2"=color_panel[6], "0.1-0.15"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black") +
          scale_fill_manual(values=c(">0.5"=color_panel[9], "0.4-0.5"=color_panel[8], "0.3-0.4"=color_panel[7], "0.2-0.3"=color_panel[6], "0.1-0.2"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black", drop = FALSE) +
          theme_bw() +
          labs(x = "", y = "", fill = "VAF", title = paste0(muts_of_interest[i,1], " g.", muts_of_interest[i,2], ":", muts_of_interest[i,3], muts_of_interest[i,4], ">", muts_of_interest[i,5], " ", muts_of_interest[i,6])) +
          # theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11.5), legend.text = element_text(size = 11.5), legend.position = "bottom") +
          guides(fill=guide_legend(nrow=2,byrow=TRUE))
        ggsave(plot = spatial_plot, filename = paste0(paste0("./output/spatial_mapping/", donor, "/"), paste(muts_of_interest[i,],collapse = "_"),"_VAF_binned_lymphocytes_only.pdf"), width = (10 * ((1.1 * (max(outline_array$x) - min(outline_array$x))) / (1.1 * (max(outline_array$y) - min(outline_array$y))))), height = 10, create.dir = TRUE)
      } else{
        spatial_plot <- ggplot() +
          geom_circle(data = circle_cut_array[which(!(circle_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(circle_cut_array$PD_ID))),], aes(x0 = x, y0 = y, r = radius, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_polygon(data = outline_array, aes(x = x, y = y, group = Cut_ID), fill = NA, color = "black", linewidth = 1.5) +
          scale_fill_manual(values=c(">0.5"=color_panel[9], "0.4-0.5"=color_panel[8], "0.3-0.4"=color_panel[7], "0.2-0.3"=color_panel[6], "0.1-0.2"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black", drop = FALSE) +
          theme_bw() +
          labs(x = "", y = "", fill = "VAF", title = paste0(muts_of_interest[i,1], " g.", muts_of_interest[i,2], ":", muts_of_interest[i,3], muts_of_interest[i,4], ">", muts_of_interest[i,5], " ", muts_of_interest[i,6])) +
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11.5), legend.text = element_text(size = 11.5), legend.position = "bottom") +
          guides(fill=guide_legend(nrow=2,byrow=TRUE))
        ggsave(plot = spatial_plot, filename = paste0(paste0("./output/spatial_mapping/", donor, "/"), paste(muts_of_interest[i,],collapse = "_"),"_VAF_binned_lymphocytes_only.pdf"), width = (10 * ((1.1 * (max(outline_array$x) - min(outline_array$x))) / (1.1 * (max(outline_array$y) - min(outline_array$y))))), height = 10, create.dir = TRUE)
      }
      if (show_plot) {print(spatial_plot)}
    }
  } else if (length(freehand_cut_array) > 4 & length(circle_cut_array) <= 4) {
    for(i in 1:nrow(muts_of_interest)){
      if(restrict_to_lymphocytes){
        spatial_plot <- ggplot() +
          geom_polygon(data = freehand_cut_array[which((grepl(pattern = "LYMPHOID_AGGREGATE", x = freehand_cut_array$Cut_ID) | grepl(pattern = "LYMPHOCYTES", x = freehand_cut_array$Cut_ID) | grepl(pattern = "LYMPHOEPITHELIAL", x = freehand_cut_array$Cut_ID)) & !(freehand_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(freehand_cut_array$PD_ID))),], aes(x = x, y = y, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          # geom_circle(data = circle_cut_array[which((grepl(pattern = "LYMPHOCYTES", x = circle_cut_array$Cut_ID) | grepl(pattern = "LYMPHOEPITHELIAL", x = circle_cut_array$Cut_ID)) & !(circle_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(circle_cut_array$PD_ID))),], aes(x0 = x, y0 = y, r = radius, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5) +
          geom_polygon(data = outline_array, aes(x = x, y = y, group = Cut_ID), fill = NA, color = "black", linewidth = 1.5) +
          # scale_fill_manual(values=c(">0.3"=color_panel[9], "0.25-0.3"=color_panel[8], "0.2-0.25"=color_panel[7], "0.15-0.2"=color_panel[6], "0.1-0.15"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black") +
          scale_fill_manual(values=c(">0.5"=color_panel[9], "0.4-0.5"=color_panel[8], "0.3-0.4"=color_panel[7], "0.2-0.3"=color_panel[6], "0.1-0.2"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black", drop = FALSE) +
          theme_bw() +
          labs(x = "", y = "", fill = "VAF", title = paste0(muts_of_interest[i,1], " g.", muts_of_interest[i,2], ":", muts_of_interest[i,3], muts_of_interest[i,4], ">", muts_of_interest[i,5], " ", muts_of_interest[i,6])) +
          # theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11.5), legend.text = element_text(size = 11.5), legend.position = "bottom") +
          guides(fill=guide_legend(nrow=2,byrow=TRUE))
        ggsave(plot = spatial_plot, filename = paste0(paste0("./output/spatial_mapping/", donor, "/"), paste(muts_of_interest[i,],collapse = "_"),"_VAF_binned_lymphocytes_only.pdf"), width = (10 * ((1.1 * (max(outline_array$x) - min(outline_array$x))) / (1.1 * (max(outline_array$y) - min(outline_array$y))))), height = 10, create.dir = TRUE)
      } else{
        spatial_plot <- ggplot() +
          geom_polygon(data = freehand_cut_array[which(!(freehand_cut_array$PD_ID %in% ID_map$PD_ID[which(ID_map$Proceed == "N")]) & !(is.na(freehand_cut_array$PD_ID))),], aes(x = x, y = y, group = Cut_ID, fill = !!(sym(paste(paste(muts_of_interest[i,],collapse = "_"),"mut_vaf_binned",sep = "_")))), color = "gray80", linewidth = 0.5, show.legend = TRUE) +
          geom_polygon(data = outline_array, aes(x = x, y = y, group = Cut_ID), fill = NA, color = "black", linewidth = 1.5) +
          scale_fill_manual(values=c(">0.5"=color_panel[9], "0.4-0.5"=color_panel[8], "0.3-0.4"=color_panel[7], "0.2-0.3"=color_panel[6], "0.1-0.2"=color_panel[5], "0.05-0.1"=color_panel[4], "0.025-0.05"=color_panel[3],"0.01-0.025"=color_panel[2],"<0.01"=color_panel[1],"0"="gray98"), na.value = "black", drop = FALSE) +
          theme_bw() +
          labs(x = "", y = "", fill = "VAF", title = paste0(muts_of_interest[i,1], " g.", muts_of_interest[i,2], ":", muts_of_interest[i,3], muts_of_interest[i,4], ">", muts_of_interest[i,5], " ", muts_of_interest[i,6])) +
          # theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11.5), legend.text = element_text(size = 11.5), legend.position = "bottom") +
          guides(fill=guide_legend(nrow=2,byrow=TRUE))
        ggsave(plot = spatial_plot, filename = paste0(paste0("./output/spatial_mapping/", donor, "/"), paste(muts_of_interest[i,],collapse = "_"),"_VAF_binned_lymphocytes_only.pdf"), width = (10 * ((1.1 * (max(outline_array$x) - min(outline_array$x))) / (1.1 * (max(outline_array$y) - min(outline_array$y))))), height = 10, create.dir = TRUE)
      }
      if (show_plot) {print(spatial_plot)}
    }
  }
}