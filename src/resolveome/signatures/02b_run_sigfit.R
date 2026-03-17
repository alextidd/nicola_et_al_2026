# runsub src/resolveome/signatures/02b_run_sigfit.R -R -M 40000

# libraries
library(sigfit)
library(magrittr)
library(ggtree)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# dirs
seq_dir <- "out/resolveome/sequoia/"
out_dir <- "out/resolveome/signatures/sigfit/"
dir.create(out_dir, showWarnings = FALSE)

# get ordering of substitutions
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[",
                   rep(sub_vec, each = 16), "]",
                   rep(c("A", "C", "G", "T"), times = 4))

# load trinuc mut matrix
trinuc_mut_mat <-
  read.table("out/resolveome/signatures/matrices/trinuc_mut_mat_hdp.txt",
             check.names = FALSE)

# load reference signatures     
ref <-
  readr::read_tsv("out/resolveome/signatures/ref_sigs/ref_sigs.tsv") %>%
  dplyr::mutate(Type = factor(Type, levels = full_vec)) %>%
  dplyr::arrange(Type) %>%
  tibble::column_to_rownames("Type") %>%
  as.matrix()

# load tree
tree <- ape::read.tree(file.path(seq_dir, "Patient_both_tree_relabelled.tree"))
tree_df <- tibble::as_tibble(ggtree::fortify(tree))

# load mutations per branch (to get indel burdens)
mut_types_per_branch <-
  file.path(seq_dir, "Patient_both_assigned_to_branches.txt") %>%
  readr::read_tsv() %>%
  dplyr::mutate(node = as.numeric(gsub("Patient_", "", SampleID))) %>%
  dplyr::group_by(node) %>%
  dplyr::summarise(n_indels = sum(nchar(Ref) != 1 | nchar(Alt) != 1),
                   n_snvs = sum(nchar(Ref) == 1 & nchar(Alt) == 1))

# get cell annotations
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
    tip_id = cell_ID,
    celltype_VDJ_recomb,
    celltype = dplyr::case_when(
      celltype_SHM == "Mature B cell" ~ "B",
      celltype_VDJ_recomb == "B cell" ~ "B",
      celltype_VDJ_recomb == "alpha-beta T cell" ~ "alpha-beta T",
      TRUE ~ celltype_VDJ_recomb),
    celltype = factor(celltype, levels = c("B", "alpha-beta T",
                                           "not lymphocyte")),
    celltype_maturity = dplyr::case_when(
      celltype_SHM == "Mature B cell" ~ "mature B",
      celltype == "B" ~ "non-mature\nB",
      celltype == "not lymphocyte" ~ "not\nlymphocyte",
      TRUE ~ celltype),
    celltype_maturity = factor(celltype_maturity,
                      levels = c("mature B", "non-mature\nB",
                                 "alpha-beta T", "not\nlymphocyte")),
    mut = dplyr::case_when(
      tip_id %in% 42:47 ~ "TNFRSF14 WT\n+ other mut\n",
      !biallelic_TNFRSF14 %in% c("WT", "FALSE_Het") &
        !biallelic_CD274 %in% c("WT", "FALSE_Het") ~ "TNFRSF14 +\nCD274 mut\n",
      !biallelic_TNFRSF14 %in% c("WT", "FALSE_Het") ~ "TNFRSF14 mut\n",
      !biallelic_CD274 %in% c("WT", "FALSE_Het") ~ "CD274 mut",
      TRUE ~ ""
    ) %>% factor(levels = c("TNFRSF14 +\nCD274 mut\n", "TNFRSF14 mut\n", "CD274 mut", "TNFRSF14 WT\n+ other mut\n", "")),
    celltype_col = dplyr::case_when(
      celltype_VDJ_recomb == "B cell" ~ "blue",
      celltype_VDJ_recomb == "alpha-beta T cell" ~ "red",
      celltype_VDJ_recomb == "not lymphocyte" ~ "grey")) %>%
  dplyr::arrange(celltype_maturity, mut) %>%
  dplyr::mutate(celltype_mut = paste0(mut, celltype_maturity),
                celltype_mut = forcats::fct_inorder(celltype_mut))

# define the final sigs
final_sigs <-
  c("SBS1", "SBS5", "SBS9", "SBS17a", "SBS17b", "SBS85",
    "machado_2022_SBSblood", "petljak_2019_ScF")
final_ref <- t(as.matrix(ref[, final_sigs]))

# keep branches with >50 muts
hdp_counts <- trinuc_mut_mat[rowSums(trinuc_mut_mat) > 50, ]

# run sigfit on each branch separately
sf_exposures <- list()
for (k in rownames(hdp_counts)) {
  print(k)

  k_file <- file.path(out_dir, paste0("sigfit_exposures_", k, ".rds"))
  if (file.exists(k_file)) {
    sf_exposures[[k]] <- readRDS(k_file)
  } else {

    # fit signatures
    sample_counts <- hdp_counts[k, , drop = FALSE]
    fit <- fit_signatures(counts = sample_counts,
                          signatures = final_ref,
                          iter = 20000, warmup = 10000, seed = 1756,
                          model = "poisson", chains = 4)

    # extract sf_exposures
    sf_exposures[[k]] <-
      retrieve_pars(fit, par = "exposures", hpd_prob = 0.95)

    # drop signatures with <5% contribution and refit
    keep_sigs <- colnames(sf_exposures[[k]]$mean)[sf_exposures[[k]]$mean > 0.05]
    if (length(keep_sigs) > 1 & length(keep_sigs) < ncol(sf_exposures[[k]]$mean)) {
      fit <- fit_signatures(counts = sample_counts,
                            signatures = final_ref[keep_sigs, , drop = FALSE],
                            iter = 20000, warmup = 10000,
                            model = "poisson", chains = 4)
      # extract sf_exposures
      sf_exposures[[k]] <-
        retrieve_pars(fit, par = "exposures", hpd_prob = 0.95)
    }

    # save intermediate results
    saveRDS(sf_exposures[[k]], k_file)
  }
}

# save exposures
saveRDS(sf_exposures, file.path(out_dir, "sigfit_exposures_per_branch.rds"))
# sf_exposures <- readRDS(file.path(out_dir, "sigfit_exposures_per_branch.rds"))

# combine exposures into matrix
sf_exp <-
  sf_exposures %>%
  purrr::map(~ .x$mean) %>%
  dplyr::bind_rows() %>%
  dplyr::select(dplyr::all_of(final_sigs)) %>%
  t()
sf_exp[is.na(sf_exp)] <- 0

# add in indels
sf_exp_w_indels <-
  sf_exp %>%
  t() %>%
  tibble::as_tibble(rownames = "sample") %>%
  dplyr::mutate(node = as.numeric(gsub("Patient_", "", sample))) %>%
  dplyr::left_join(mut_types_per_branch, by = "node") %>%
  tidyr::pivot_longer(cols = dplyr::all_of(final_sigs),
                      names_to = "signature", values_to = "exposure") %>%
  # get number of snvs per signature per branch
  dplyr::mutate(n_snvs_exp = exposure * n_snvs) %>%
  # recalculate exposures as % of all muts (snvs + indels)
  dplyr::mutate(n_muts = n_snvs + n_indels,
                exp_snvs = n_snvs_exp / n_muts,
                indels = n_indels / n_muts) %>%
  dplyr::select(sample, signature, exp_snvs, indels) %>%
  tidyr::pivot_wider(names_from = "signature", values_from = "exp_snvs",
                     id_cols = c("sample", "indels")) %>%
  tibble::column_to_rownames("sample") %>%
  t()

# save exposures with indels
saveRDS(sf_exp_w_indels, file.path(out_dir, "sigfit_exposures_per_branch_with_indels.rds"))

# get exposure colours
sig_cols <- c("black", brewer.pal(n = length(final_sigs), "Set2"))
names(sig_cols) <- rownames(sf_exp_w_indels)

# get celltype colours
celltype_cols <-
  tibble::tibble(tip_id = as.numeric(tree$tip.label)) %>%
  dplyr::left_join(cell_annots) %>%
  dplyr::pull(celltype_col)

# create tree plot with fitted signatures colored along branches
pdf(file.path(out_dir, "tree_with_branch_length_sigfit.pdf"),
    height = 10, width = 10)

# create tree
plot(tree, cex = 0.7, label.offset = 0.003 * max(tree_df$x),
     tip.color = celltype_cols, font = 1)

# for each sample, draw rectangles showing signature proportions
for (sample in colnames(sf_exp_w_indels)) {
  n <- as.numeric(substr(sample, 9, nchar(sample)))
  x_end <- tree_df$x[n]
  x_start <- tree_df$x[tree_df$parent[n]]
  x_intv <- x_end - x_start
  y <- node.height(tree)[n]
  tipnum <- sum(tree_df$isTip)

  # stack signature exposures proportionally along branch length
  for (s in rownames(sf_exp_w_indels)) {
    x_end <- x_start + sf_exp_w_indels[s, sample] * x_intv
    rect(ybottom = y - min(0.015 * tipnum, 0.3),
         ytop = y + min(0.015 * tipnum, 0.3),
         xleft = x_start, xright = x_end, col = sig_cols[s], lwd = 0.25)
    x_start <- x_end
  }
}

# axis and legend
axisPhylo(side = 1, backward = FALSE)
legend("topright", title = "signatures", legend = names(sig_cols),
       fill = sig_cols, bty = "n", cex = 0.7, ncol = 1, xjust = 0.5)
legend("bottomright", title = "celltype",
       legend = c("B", "alpha-beta T", "not lymphocyte"),
       col    = c("blue", "red", "grey"),
       pch    = 19,
       pt.cex = 1,
       cex = 0.7,
       bty    = "n")

dev.off()

# get child nodes
children <-
  tree_df %>%
  dplyr::group_by(parent) %>%
  dplyr::summarise(child = list(node)) %>%
  dplyr::rename(node = parent)

# get exposures per branch for missing branches
per_branch_exp <-
  sf_exp_w_indels %>%
  t() %>%
  tibble::as_tibble(rownames = "sample") %>%
  dplyr::mutate(node = as.numeric(gsub("Patient_", "", sample)),
                over_50 = TRUE) %>%
  dplyr::full_join(tree_df) %>%
  # get child nodes
  dplyr::left_join(children) %>%
  dplyr::mutate(over_50 = ifelse(is.na(over_50), FALSE, over_50)) %>%
  split(.$over_50)

# remove empty columns, unnest child nodes
per_branch_exp$`FALSE` <-
  per_branch_exp$`FALSE` %>%
  dplyr::select(-dplyr::all_of(c("sample", "indels", final_sigs))) %>%
  tidyr::unnest(cols = c("child"))

# if a branch with <50 muts has children with >50 muts, assign the average of
# the children's exposures to the parent branch
per_branch_exp_from_child <-
  per_branch_exp$`FALSE` %>%
  # get children
  dplyr::inner_join(
    per_branch_exp$`TRUE` %>%
      dplyr::select(dplyr::all_of(c("node", "indels", final_sigs))),
    by = c("child" = "node")) %>%
  dplyr::group_by(dplyr::across(-c(child, dplyr::all_of(c("indels", final_sigs))))) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(c("indels", final_sigs)), mean), .groups = "drop")

# if a branch with <50 muts has no children with >50 muts, keep moving up the 
# tree until finding exposures to assign
per_branch_exp_grandchild <-
  per_branch_exp$`FALSE` %>%
  dplyr::anti_join(per_branch_exp_from_child, by = "node") %>%
  # get grandchildren
  dplyr::left_join(
    children %>%
      tidyr::unnest(cols = c("child")) %>%
      dplyr::rename(child = node, grandchild = child)
  ) %>%
  # add grandchildren's exposures
  dplyr::inner_join(
    per_branch_exp$`TRUE` %>%
      dplyr::select(dplyr::all_of(c("node", "indels", final_sigs))),
    by = c("grandchild" = "node")) %>%
  dplyr::group_by(dplyr::across(-c(child, grandchild, dplyr::all_of(c("indels", final_sigs))))) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(c("indels", final_sigs)), mean), .groups = "drop")

# check there are no remaining branches without exposures
per_branch_exp$`FALSE` %>%
  dplyr::anti_join(per_branch_exp_from_child, by = "node") %>%
  dplyr::anti_join(per_branch_exp_grandchild, by = "node")

# combine filled in exposures
per_branch_exp_final <-
  list(per_branch_exp$`TRUE`,
       per_branch_exp_from_child,
       per_branch_exp_grandchild) %>%
  dplyr::bind_rows()

# save exposures per branch with filled in branches
per_branch_exp_final %>%
  saveRDS(file.path(out_dir, "per_branch_exposures_filled.rds"))

# get all nodes leading to each tip
nodes_to_tips <-
  lapply(1:Ntip(tree), function(tip) {
    nodepath(tree, from = tip, to = Ntip(tree) + 1)
  }) %>%
  setNames(tree$tip.label) %>%
  tibble::enframe(name = "tip_id", value = "node") %>%
  tidyr::unnest(cols = node)

# get exposures per cell
exp_per_cell <-
  nodes_to_tips %>%
  dplyr::left_join(per_branch_exp_final, by = c("node")) %>%
  # calculate mutations per exposure
  tidyr::pivot_longer(cols = dplyr::all_of(c("indels", final_sigs)),
                      names_to = "signature", values_to = "exposure") %>%
  dplyr::mutate(n_muts_exp = exposure * branch.length) %>%
  dplyr::group_by(tip_id, signature) %>%
  dplyr::summarise(n_muts_exp = sum(n_muts_exp)) %>%
  # calculate mutations per cell
  dplyr::group_by(tip_id) %>%
  dplyr::mutate(total_n_muts = sum(n_muts_exp)) %>%
  dplyr::ungroup() %>%
  # add cell annots
  dplyr::left_join(cell_annots %>% dplyr::mutate(tip_id = as.character(tip_id))) %>%
  # factor axes
  dplyr::mutate(tip_id = forcats::fct_reorder(tip_id, -total_n_muts),
                signature = factor(signature, levels = rev(c("indels", final_sigs)))) %>%
  # calculate proportions and numbers of mutations
  dplyr::mutate(`n mutations` = n_muts_exp,
                `% mutations` = n_muts_exp / total_n_muts) %>%
  tidyr::pivot_longer(cols = c("n mutations", "% mutations"),
                      names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("n mutations", "% mutations")))

# save exposures per cell
exp_per_cell %>%
  saveRDS(file.path(out_dir, "exposures_per_cell.rds"))

# plot exposures per cell
pdf(file.path(out_dir, "signature_burden_per_cell.pdf"), width = 17, height = 8)
p <-
  exp_per_cell %>%
  ggplot(aes(x = tip_id, y = value, fill = signature)) +
  geom_col(width = 1.05) +
  scale_fill_manual(values = sig_cols) +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.spacing.y = unit(0.6, "lines"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.text.y = element_text(angle = 90)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = c(0, 0))
p +
  ggh4x::facet_nested(metric ~ celltype_mut, scales = "free", space = "free_x") +
  ggtitle("by celltype and TNFRSF14/CD274 status")
dev.off()

# plot celltype spectra
mut_type_cols <-
  c("C>A" = "dodgerblue", "C>G" = "black", "C>T" = "red",
    "T>A" = "grey70", "T>C" = "olivedrab3", "T>G" = "plum2")
p_dat <-
  trinuc_mut_mat %>%
  tibble::as_tibble(rownames = "sample") %>%
  dplyr::mutate(node = as.numeric(gsub("Patient_", "", sample))) %>%
  dplyr::left_join(nodes_to_tips) %>%
  dplyr::group_by(tip_id) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(colnames(trinuc_mut_mat)), sum)) %>%
  dplyr::left_join(cell_annots %>% dplyr::mutate(tip_id = as.character(tip_id))) %>%
  dplyr::mutate(celltype_mat = gsub("\n", " ", celltype_maturity) %>%
    factor(levels = c("mature B", "non-mature B", "alpha-beta T", "not lymphocyte"))) %>%
  dplyr::group_by(celltype_mat) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(colnames(trinuc_mut_mat)), sum)) %>%
  tidyr::pivot_longer(cols = -c("celltype_mat")) %>%
  dplyr::mutate(trinuc = paste0(substr(name, 5, 5), substr(name, 1, 1), substr(name, 7, 7)),
                mut_type = substr(name, 1, 3))

pdf(file.path(out_dir, "celltype_trinuc_spectra.pdf"),
    width = 12, height = 6)
p_dat %>%
  ggplot(aes(x = trinuc, y = value, fill = mut_type)) +
  geom_col() +
  scale_fill_manual(values = mut_type_cols) +
  facet_grid(celltype_mat ~ mut_type, scales = "free", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(family = "mono"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        legend.position = "none", strip.text.y = element_text(angle = 90)) +
  guides(x = guide_axis(angle = -90)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
dev.off()
