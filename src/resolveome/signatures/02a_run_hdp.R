# runsub src/resolveome/signatures/02a_run_hdp.R -R -M 10000

# libraries
library(magrittr)
library(hdp)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

# hdp params
chains <- 20
burnin <- 20000

# dirs
mat_dir <- "out/resolveome/signatures/matrices/"
seq_dir <- "out/resolveome/sequoia/"
out_dir <-
  paste0("out/resolveome/signatures/hdp/branches_over_50/chains=", chains,
         ",burnin=", burnin)
dir.create(out_dir, showWarnings = FALSE)

# palette
mut_colours <- c("dodgerblue", "black", "red", "grey70", "olivedrab3", "plum2")
mut_colours_96 <- rep(mut_colours, each = 16)

# read branch matrix, filter to branches with >50 muts to avoid overfitting
trinuc_mut_mat <-
  file.path(mat_dir, "trinuc_mut_mat_hdp.txt") %>%
  read.table(header = TRUE, row.names = 1, check.names = FALSE)
trinuc_mut_mat <-
  trinuc_mut_mat[apply(trinuc_mut_mat, 1, sum) > 50, ]

# get cell annots
cell_shm <-
  "data/resolveome/manual_inspection/20250902_pta_additional_annotation_H1.tsv" %>%
  readr::read_tsv() %>%
  dplyr::select(well_ID, celltype_SHM)
cell_annots <-
  "data/resolveome/manual_inspection/H1_PD63118_pta_additional_annotation.tsv" %>%
  readr::read_tsv() %>%
  dplyr::left_join(cell_shm) %>%
  dplyr::transmute(
    cell_id = well_ID,
    celltype = case_when(celltype_SHM == "Mature B cell" ~ "Mature B cell",
                         celltype_VDJ_recomb == "B cell" ~ "Non-mature B cell",
                         TRUE ~ celltype_VDJ_recomb),
    celltype = factor(celltype,
                      levels = c("Mature B cell", "Non-mature B cell",
                                 "alpha-beta T cell", "not lymphocyte")),
    biallelic_TNFRSF14, biallelic_CD274)

# create sample-patient key (assumes pdids)
samples_to_patients <-
  tibble::tibble(sample = rownames(trinuc_mut_mat)) %>%
  dplyr::mutate(patient = substr(sample, 1, 7))
samples_per_patient <- table(samples_to_patients$patient)

# create branch-cell key
tree <-
  read.tree(file.path(seq_dir, "Patient_both_tree_with_branch_length.tree"))
tree_df <-
  ggtree::fortify(tree)
branches_to_cells <-
  lapply(1:Ntip(tree), function(tip) {
    nodepath(tree, from = tip, to = Ntip(tree) + 1)
  }) %>%
  setNames(tree$tip.label) %>%
  tibble::enframe(name = "id", value = "node") %>%
  tidyr::unnest(cols = node) %>%
  dplyr::left_join(ggtree::fortify(tree))%>%
  dplyr::transmute(
    id,
    cell_id = gsub("_dna.*", "", id),
    branch = node, n_muts = branch.length)

# setup hdp hierarchy
hdp_in <-
  hdp_init(
    # parent-child relationships: 0->1->2 hierarchy
    ppindex = c(0, rep(1, length(samples_per_patient)),
                rep(2:(length(samples_per_patient) + 1),
                times = samples_per_patient)),
    # concentration parameter indices for each dp
    cpindex = c(1, rep(2, length(samples_per_patient)),
                rep(3:(length(samples_per_patient) + 2),
                times = samples_per_patient)),
    # base distribution (uniform across 96 contexts)
    hh = rep(1, 96),
    # gamma hyperparameters for concentration parameters
    alphaa = rep(1, length(samples_per_patient) + 2),
    alphab = rep(1, length(samples_per_patient) + 2))

# attach mutation count data to sample-level dps
# samples are nodes 2 through (number of patients + 1)
hdp_mut <-
  hdp_setdata(hdp_in,
              dpindex = (length(samples_per_patient) + 2):numdp(hdp_in),
              as.data.frame(trinuc_mut_mat))

# run hdp mcmc sampling
hdp_chains <-
  c(1:chains) %>%
  purrr::map(function(i) {
    print(paste("running chain", i, "of", chains))

    hdp_activated <-
      dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc = 10, seed = i * 200)

    hdp_posterior(hdp_activated,
                  burnin = burnin, n = 100, space = 200, cpiter = 3,
                  seed = i * 1e3)
  })

# save chains
saveRDS(hdp_chains, file = file.path(out_dir, "hdp_chains.rds"))
# hdp_chains <- readRDS(file.path(out_dir, "hdp_chains.rds"))

# create multi-chain object for convergence diagnostics
hdp_multi_chains_1 <- hdp_multi_chain(hdp_chains)

# generate quality control plots to assess chain convergence
pdf(file.path(out_dir, "QC_plots.pdf"))
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# likelihood traces (should stabilize after burnin)
p1 <- lapply(chains(hdp_multi_chains_1), plot_lik, bty = "L", start = 1000)
# number of active clusters over time
p2 <- lapply(chains(hdp_multi_chains_1), plot_numcluster, bty = "L")
# proportion of data assigned to clusters
p3 <- lapply(chains(hdp_multi_chains_1), plot_data_assigned, bty = "L")
dev.off()

# extract consensus signature profiles across chains
# this averages signatures with high posterior support and
# filters out those that appear in few chains (low confidence)
hdp_multi_chains <- hdp_extract_components(hdp_multi_chains_1)

# save consensus
saveRDS(hdp_multi_chains, file = file.path(out_dir, "hdp_multi_chains.rds"))
# hdp_multi_chains <- readRDS(file.path(out_dir, "hdp_multi_chains.rds"))

# plot number of mutations per component
pdf(file.path(out_dir, "de_novo_sigs_data_items_assigned_plot_no_hierarchy.pdf"))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
plot_comp_size(hdp_multi_chains, bty = "L")
dev.off()

# get median number of mutations per component
hdp_multi_chains %>%
  comp_categ_counts() %>%
  sapply(rowSums) %>%
  t() %>%
  tibble::as_tibble(rownames = "component") %>%
  tidyr::pivot_longer(
    cols = -c("component"), names_to = "iteration", values_to = "n_muts") %>%
  dplyr::group_by(component) %>%
  dplyr::summarise(median_n_muts = median(n_muts)) %>%
  dplyr::mutate(total_muts = sum(trinuc_mut_mat),
                prop_muts = median_n_muts / total_muts)

# plot extracted signatures
# component 0 is the "background" component (noise/artefacts)
for (i in 0:hdp_multi_chains@numcomp) {
  pdf(file.path(out_dir, paste0("hdp_component_", i, ".pdf")),
      width = 12, height = 4)
  plot_comp_distn(
    hdp_multi_chains,
    cat_names = sapply(strsplit(colnames(mut_count), "\\."), `[`, 4),
    grouping = as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each = 16)),
    col = mut_colours, comp = i, col_nonsig = "grey80",
    show_group_labels = TRUE)
  dev.off()
}

# plot signature activities across all samples
pdf(file.path(out_dir, "signature_attribution.pdf"), width = 10, height = 8)
plot_dp_comp_exposure(
  hdp_multi_chains, 
  # show only sample-level dps (exclude higher-level hierarchy)
  dpindices = (
    length(hdp_multi_chains@comp_dp_counts) -
      nrow(trinuc_mut_mat) + 1):length(hdp_multi_chains@comp_dp_counts),
  incl_nonsig = TRUE, ylab_exp = "signature exposure", leg.title = "signature",
  col = RColorBrewer::brewer.pal(12, "Set3"))
dev.off()

# extract signature-sample assignment matrix
# this shows which signatures are active in which samples
dp_distn <- comp_dp_distn(hdp_multi_chains)
exposures <-
  t(dp_distn$mean[length(samples_per_patient) + 1 + 1:nrow(trinuc_mut_mat), ,
                  drop = FALSE])
rownames(exposures) <- paste0("N", rownames(exposures))
colnames(exposures) <- rownames(trinuc_mut_mat)
write.table(exposures, file.path(out_dir, "mean_assignment_hdp.txt"))

# plot signature activities per cell
# aggregate branch exposures along each cell's lineage path
cell_exposures <-
  branches_to_cells %>%
  dplyr::mutate(branch = paste0("Patient_", branch)) %>%
  dplyr::filter(branch %in% colnames(exposures)) %>%
  dplyr::left_join(
    exposures %>%
      t() %>%
      tibble::as_tibble(rownames = "branch") %>%
      tidyr::pivot_longer(cols = -c("branch"), names_to = "component", values_to = "exposure"),
    by = "branch", relationship = "many-to-many"
  ) %>%
  dplyr::mutate(n_muts_exp = n_muts * exposure) %>%
  dplyr::group_by(cell_id, component) %>%
  dplyr::summarise(n_muts_exp = sum(n_muts_exp), .groups = "drop") %>%
  dplyr::group_by(cell_id) %>%
  dplyr::mutate(n_muts = sum(n_muts_exp)) %>%
  # add cell annotations
  dplyr::left_join(cell_annots)

# plot exposures per cell
pdf(file.path(out_dir, "signature_attribution_per_cell.pdf"), width = 12, height = 8)
p <-
  cell_exposures %>%
  ggplot(aes(x = reorder(cell_id, -n_muts), y = n_muts_exp, fill = component)) +
  guides(x = guide_axis(angle = -90)) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.x = element_blank())
(p + geom_col() + facet_grid(~ celltype, scales = "free_x", space = "free_x") + ggtitle("By celltype")) /
(p + geom_col(position = "fill") + facet_grid(~ celltype, scales = "free_x", space = "free_x"))
(p + geom_col() + facet_grid(~ biallelic_TNFRSF14, scales = "free_x", space = "free_x") + ggtitle("By biallelic_TNFRSF14")) / 
(p + geom_col(position = "fill") + facet_grid(~ biallelic_TNFRSF14, scales = "free_x", space = "free_x"))
(p + geom_col() + facet_grid(~ biallelic_CD274, scales = "free_x", space = "free_x") + ggtitle("By biallelic_CD274")) / 
(p + geom_col(position = "fill") + facet_grid(~ biallelic_CD274, scales = "free_x", space = "free_x"))
dev.off()

# plot exposures per plate
pdf(file.path(out_dir, "signature_attribution_per_plate.pdf"),
    width = 12, height = 8)
p_dat <-
  cell_exposures %>%
  dplyr::mutate(plate = gsub("_.*", "", cell_id))
(p_dat %>%
  ggplot(aes(x = reorder(cell_id, -n_muts), y = n_muts_exp, fill = component)) +
  guides(x = guide_axis(angle = -90)) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.x = element_blank()) +
  geom_col() +
  facet_grid(~ plate, scales = "free_x", space = "free_x") +
  ggtitle("Signature exposure per cell per plate")) /
(p_dat %>%
  dplyr::filter(component == "N5") %>%
  dplyr::mutate(prop_muts_exp = n_muts_exp / n_muts) %>%
  ggplot(aes(x = plate, y = prop_muts_exp)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(height = 0) +
  theme_classic() +
  ggtitle("% N5 exposure per cell per plate") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format"))
dev.off()

# extract signature profiles (96-context distributions)
hdp_sigs <- as.data.frame(t(comp_categ_distn(hdp_multi_chains)$mean))
colnames(hdp_sigs) <- paste0("N", colnames(hdp_sigs))
write.table(hdp_sigs, file.path(out_dir, "hdp_sigs.txt"))
# hdp_sigs <- read.table(file.path(out_dir, "hdp_sigs.txt"), header = TRUE, row.names = 1)

# get ordering of substitutions
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[",
                   rep(sub_vec, each = 16), "]",
                   rep(c("A", "C", "G", "T"), times = 4))

# load cosmic reference signatures
refs <-
  c("v2", "v3.5") %>%
  purrr::set_names() %>%
  purrr::map(function(v) {
    readr::read_tsv(paste0("../../reference/cosmic/COSMIC_", v, "_SBS_GRCh38.txt")) %>%
      dplyr::mutate(Type = factor(Type, levels = full_vec)) %>%
      dplyr::arrange(Type) %>%
      tibble::column_to_rownames("Type") %>%
      as.matrix()
  })

# use v3.5, replace 0 with small value and normalise
ref <- refs$v3.5
ref[is.na(ref) | ref == 0] <- 0.00001
ref <- t(t(ref) / colSums(ref))

# SBSblood from machado 2022
machado <-
  readr::read_tsv("data/signatures/machado_2022/S8_finalsignaturetable.tsv") %>%
  tidyr::pivot_longer(cols = -c("Signature")) %>%
  dplyr::mutate(Type = factor(paste0(substr(name, 1, 1), "[",
                                     substr(name, 2, 2), ">",
                                     substr(name, 6, 6), "]",
                                     substr(name, 3, 3)),
                              levels = full_vec)) %>%
  tidyr::pivot_wider(names_from = "Signature", values_from = "value") %>%
  tibble::column_to_rownames("Type")

# artefact signature ScF from petljak 2019
petljak <-
  readr::read_tsv("data/signatures/petljak_2019/mmc1.tsv") %>%
  dplyr::mutate(
    Type = factor(paste0(substr(`Mutation Subtype`, 1, 1), "[",
                         `Mutation Type`, "]",
                         substr(`Mutation Subtype`, 3, 3)),
                  levels = full_vec)) %>%
  tibble::column_to_rownames("Type")

# artefact signature ScB from lodato 2018 (MDA)
lodato <-
  "data/signatures/lodato_2018/Lodato2018_SignatureData_Aging.csv" %>%
  readr::read_csv() %>%
  dplyr::mutate(
    Type = factor(paste0(substr(`...1`, 1, 1), "[", `...2`, "]",
                         substr(`...1`, 3, 3)), levels = full_vec)) %>%
  tibble::column_to_rownames("Type")

# universal PTA artefact signature from luquette 2022
load("data/signatures/luquette_2022/snv.artifact.signature.v3.rda")
luquette <-
  snv.artifact.signature.v3 %>%
  tibble::enframe(value = "pta_artefact") %>%
  dplyr::mutate(Type = paste0(substr(name, 1, 1), "[", substr(name, 5, 7), "]",
                substr(name, 3, 3))) %>%
  tibble::column_to_rownames("Type")

# add all additional signatures to ref
ref <-
  cbind(
    ref,
    SBS17 = refs$v2[rownames(ref), "Signature_17"],
    machado_2022_SBSblood = machado[rownames(ref), "SBSblood"],
    machado_2022_SignatureIg = machado[rownames(ref), "Signature.Ig"],
    lodato_2018_ScB = lodato[rownames(ref), "B"],
    petljak_2019_ScF = petljak[rownames(ref), "SBS sc_F"],
    luquette_2022_PTA_artefact = luquette[rownames(ref), "pta_artefact"])

sub_vec <-
  c("C>A","C>G","C>T","T>A","T>C","T>G")
full_vec <-
  paste0(rep(c("A", "C", "G", "T"), each = 4), "[", rep(sub_vec, each = 16),
         "]", rep(c("A", "C", "G", "T"), times = 4))
# ref <- ref[full_vec, ]

# calculate cosine similarity between hdp and reference signatures
# cosine similarity ranges 0-1, with 1 = perfect match
cosine_matrix <- data.frame(matrix(nrow = ncol(hdp_sigs), ncol = ncol(ref)))
rownames(cosine_matrix) <- colnames(hdp_sigs)
colnames(cosine_matrix) <- colnames(ref)
for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n, m] <-
      lsa::cosine(x = hdp_sigs[, rownames(cosine_matrix)[n]],
                  y = ref[, colnames(cosine_matrix)[m]])
  }
}

# visualize cosine similarity matrix as heatmap
pdf(file.path(out_dir, "cosine_similarities.pdf"), height = 5, width = 15)
color.palette <- colorRampPalette(c("white", "orange", "purple"))
t(cosine_matrix[dim(cosine_matrix)[1]:1, ]) %>%
  lattice::levelplot(col.regions = color.palette, aspect = "fill", 
                     scales = list(x = list(rot = 90)))
dev.off()

# define candidate reference signatures and components for deconvolution
# (biologically plausible signatures for this tissue/context)
core_sigs <-
  c("SBS1", "SBS5", "SBS9", "SBS17a", "SBS17b", "machado_2022_SBSblood")
candidate_sigs <- 
  list(
    c("luquette_2022_PTA_artefact", core_sigs),
    c("lodato_2018_ScB", core_sigs),
    c("lodato_2018_ScB", "SBS85", "SBS40a", core_sigs),
    c("petljak_2019_ScF", "SBS85", core_sigs))
names(candidate_sigs) <-
  lapply(candidate_sigs, function(x) {paste(x, collapse = ",")}) %>%
  gsub("SBS", "", .)

# plot cosine similarity for candidate signatures
pdf(file.path(out_dir, "cosine_similarities_candidates.pdf"), height = 5, width = 8)
color.palette <- colorRampPalette(c("white", "orange", "purple"))
cosine_subset <- t(cosine_matrix[, unique(unlist(candidate_sigs))])
lattice::levelplot(
  cosine_subset,
  col.regions = color.palette, aspect = "fill",
  scales = list(x = list(rot = 90)),
  panel = function(x, y, z, ...) {
    lattice::panel.levelplot(x, y, z, ...)
    lattice::panel.text(x, y, round(z, 2), cex = 0.8)
  }
)
dev.off()

# perform expectation maximisation
purrr::map2(names(candidate_sigs), candidate_sigs,
            function(subdir, candidate_sigs_i) {
  
  out_dir_i <- file.path(out_dir, paste0("em_", subdir))
  dir.create(out_dir_i, showWarnings = FALSE)

  # save signatures used
  writeLines(candidate_sigs_i,
             file.path(out_dir_i, "em_candidate_signatures.txt"))

  # setup for EM
  signatures <- t(ref[, candidate_sigs_i])
  sample_list <- colnames(hdp_sigs)
  profiles <- hdp_sigs[, sample_list]

  # initialize deconvolution results matrix
  signature_fraction <- matrix(NA, nrow = nrow(signatures),
                              ncol = length(sample_list))
  rownames(signature_fraction) <- rownames(signatures)
  colnames(signature_fraction) <- sample_list
  maxiter <- 1000

  # run em algorithm for each hdp signature (round 1)
  for (j in seq_along(sample_list)) {
    freqs <- profiles[, j]
    freqs[is.na(freqs)] <- 0

    # em algorithm to estimate the signature contribution
    # initialize with random signature weights
    alpha <- runif(nrow(signatures))
    alpha <- alpha / sum(alpha)

    # em iterations until convergence
    for (iter in 1:maxiter) {
      # e-step: calculate expected contributions
      contr <- t(array(alpha, dim = c(nrow(signatures), 96))) * t(signatures)
      probs <- contr / array(rowSums(contr), dim = dim(contr))
      probs <- probs * freqs

      # m-step: update signature weights
      old_alpha <- alpha
      alpha <- colSums(probs) / sum(probs)

      # check convergence
      if (sum(abs(alpha - old_alpha)) < 1e-5) {
        break
      }
    }

    # saving the signature contributions for the sample
    print(paste0("Round 1: ", j, "/", length(sample_list)))
    signature_fraction[, j] <- alpha
  }

  # filter to signatures contributing >10% to any sample
  sigs_deconv_r2 <- list()
  for (n in seq_along(sample_list)) {
    sigs_deconv_r2[[n]] <-
      rownames(signature_fraction)[signature_fraction[, n] > 0.1]
  }
  names(sigs_deconv_r2) <- colnames(signature_fraction)

  # identify samples requiring further deconvolution (multiple signatures)
  sigs_to_deconv <-
    names(sigs_deconv_r2)[unlist(lapply(sigs_deconv_r2, length)) > 1]

  # second round deconvolution with refined signature sets
  ref_sigs_r2 <- sort(unique(unlist(sigs_deconv_r2)))
  signature_fraction_r2 <- matrix(NA, ncol = length(sigs_to_deconv),
                                  nrow = length(ref_sigs_r2))
  rownames(signature_fraction_r2) <- ref_sigs_r2
  colnames(signature_fraction_r2) <- sigs_to_deconv
  reconst_cosines <- list()

  # repeat the deconvolution with the identified constitutive signatures
  n <- 1; for (s in sigs_to_deconv) {
    # use only signatures identified in round 1 for this sample
    gdsigs <- sigs_deconv_r2[[s]]
    signatures_r2 <- t(ref[, gdsigs])

    freqs <- profiles[, s]
    freqs[is.na(freqs)] <- 0

    # em algorithm (same as above but with reduced signature set)
    alpha <- runif(nrow(signatures_r2))
    alpha <- alpha / sum(alpha)
    names(alpha) <- gdsigs

    for (iter in 1:maxiter) {
      contr <- t(array(alpha, dim = c(nrow(signatures_r2), 96))) * t(signatures_r2)
      probs <- contr / array(rowSums(contr), dim = dim(contr))
      probs <- probs * freqs
      old_alpha <- alpha
      alpha <- colSums(probs) / sum(probs)
      names(alpha) <- gdsigs
      if (sum(abs(alpha - old_alpha)) < 1e-5) {
        break
      }
    }

    # saving the signature contributions for the sample
    signature_fraction_r2[gdsigs, n] <- alpha[gdsigs]
    n <- n + 1

    # reconstruct signature and calculate fit quality
    reconsbs <- rep(0, 96)
    for (g in gdsigs) {
      reconsbs <- reconsbs + (ref[, g] * alpha[g])
    }
    cosine_reconst <- lsa::cosine(x = reconsbs, y = hdp_sigs[, s])
    print(paste0(s, ": cosine = ", round(cosine_reconst, 3)))
    reconst_cosines[[s]] <- as.numeric(cosine_reconst[1, 1])

    # plot reconstitution
    pdf(file.path(out_dir_i, paste0("hdp_", s, "_reconstitution.pdf")), height = 10)
    par(mfrow = c(length(alpha) + 2, 1))
    par(mar = c(1, 2, 4, 1))
    barplot(hdp_sigs[, s], col = mut_colours_96, main = paste0("HDP ", s),
            names.arg = "")
    barplot(reconsbs, col = mut_colours_96,
            main = paste0("Reconstituted ", s, " (cosine = ",
                          round(cosine_reconst, 2), ")"))
    for (g in gdsigs) {
      barplot(ref[, g], col = mut_colours_96,
              main = paste0(g, " (", round(alpha[g] * 100, 1), "%)"))
    }
    dev.off()
  }

  # save deconvolution results
  write.table(signature_fraction,
              file.path(out_dir_i, "em_signature_fraction_r1.txt"),
              quote = FALSE, sep = "\t")
  write.table(signature_fraction_r2,
              file.path(out_dir_i, "em_signature_fraction_r2.txt"),
              quote = FALSE, sep = "\t")
  signature_fraction_r2 %>%
    tibble::as_tibble(rownames = "signature") %>%
    tidyr::pivot_longer(cols = -c("signature"),
                        names_to = "component", values_to = "contribution") %>%
    dplyr::filter(!is.na(contribution)) %>%
    dplyr::arrange(component, -contribution) %>%
    readr::write_tsv(file.path(out_dir_i, "em_signature_contributions.tsv"))

  # plot final deconvolution heatmap
  pdf(file.path(out_dir_i, "em_deconvolution_heatmap.pdf"), height = 5, width = 10)
    color.palette <- colorRampPalette(c("white", "orange", "purple"))
    lattice::levelplot(
      signature_fraction[nrow(signature_fraction):1, ],
      col.regions = color.palette, aspect = "fill",
      scales = list(x = list(rot = 90)),
      main = "EM Signature Deconvolution (Round 1)"
    )
  dev.off()

})