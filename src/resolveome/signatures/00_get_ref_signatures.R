# libraries
library(magrittr)
library(ggplot2)
library(tidyverse)

# dirs
out_dir <- "out/resolveome/signatures/ref_sigs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# get ordering of substitutions
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[",
                   rep(sub_vec, each = 16), "]",
                   rep(c("A", "C", "G", "T"), times = 4))

# load cosmic reference signatures
refs <-
  c("v2", "v3.4", "v3.5") %>%
  purrr::set_names() %>%
  purrr::map(function(v) {
    readr::read_tsv(paste0("../../reference/cosmic/COSMIC_", v, "_SBS_GRCh38.txt")) %>%
      dplyr::mutate(Type = factor(Type, levels = full_vec)) %>%
      tibble::column_to_rownames("Type") %>%
      as.matrix()
  })

# use v3.5
ref <- refs$v3.5

# SBSblood and SignatureIg from machado 2022
machado <-
  readr::read_tsv("data/signatures/machado_2022/S8_finalsignaturetable.tsv") %>%
  tidyr::pivot_longer(cols = -c("Signature")) %>%
  dplyr::mutate(Type = factor(paste0(substr(name, 1, 1), "[",
                                     substr(name, 2, 2), ">",
                                     substr(name, 6, 6), "]",
                                     substr(name, 3, 3)),
                              levels = full_vec)) %>%
  tidyr::pivot_wider(names_from = "Signature", values_from = "value") %>%
  dplyr::select(Type,
                machado_2022_SBSblood = SBSblood,
                machado_2022_SignatureIg = Signature.Ig)

# artefact signature ScF from petljak 2019
petljak <-
  readr::read_tsv("data/signatures/petljak_2019/mmc1.tsv") %>%
  dplyr::mutate(
    Type = factor(paste0(substr(`Mutation Subtype`, 1, 1), "[",
                         `Mutation Type`, "]",
                         substr(`Mutation Subtype`, 3, 3)),
                  levels = full_vec)) %>%
  dplyr::select(Type, petljak_2019_ScF = `SBS sc_F`)

# artefact signature ScB from lodato 2018 (MDA)
lodato <-
  "data/signatures/lodato_2018/Lodato2018_SignatureData_Aging.csv" %>%
  readr::read_csv() %>%
  dplyr::mutate(
    Type = factor(paste0(substr(`...1`, 1, 1), "[", `...2`, "]",
                         substr(`...1`, 3, 3)), levels = full_vec)) %>%
  dplyr::select(Type, lodato_2018_ScB = B)

# universal PTA artefact signature from luquette 2022
load("data/signatures/luquette_2022/snv.artifact.signature.v3.rda")
luquette <-
  snv.artifact.signature.v3 %>%
  tibble::enframe(value = "pta_artefact") %>%
  dplyr::mutate(Type = paste0(substr(name, 1, 1), "[", substr(name, 5, 7), "]",
                substr(name, 3, 3))) %>%
  dplyr::select(Type, luquette_2022_PTA_artefact = pta_artefact)

# combine custom refs
custom_ref <-
  ref %>%
  dplyr::select(Type) %>%
  dplyr::left_join(machado) %>%
  dplyr::left_join(lodato) %>%
  dplyr::left_join(petljak) %>%
  dplyr::left_join(luquette)

# add all additional signatures to ref
ref <-
  ref %>%
  dplyr::left_join(custom_ref)

# replace 0 with small value and normalise
ref <- ref %>% tibble::column_to_rownames("Type") %>% as.matrix()
ref[is.na(ref) | ref == 0] <- 0.00001
ref <- t(t(ref) / colSums(ref))
ref <- ref %>% tibble::as_tibble(rownames = "Type")

# save ref
ref %>%
  readr::write_tsv(file.path(out_dir, "ref_sigs.tsv"))

# plot all custom reference signatures
pdf(file.path(out_dir, "custom_ref_sigs.pdf"), width = 12, height = 3)
custom_ref %>%
  tidyr::pivot_longer(cols = -Type,
                      names_to = "signature", values_to = "weight") %>%
  dplyr::mutate(sub = stringr::str_extract(Type, "(?<=\\[).*?(?=\\])")) %>%
  dplyr::arrange(sub, Type) %>%
  dplyr::mutate(Type = forcats::fct_inorder(Type)) %>%
  split(.$signature) %>%
  purrr::map(function(df) {
    df %>%
      ggplot(aes(x = Type, y = weight, fill = sub)) +
      geom_col() +
      theme_classic() +
      guides(x = guide_axis(angle = -90)) +
      scale_fill_manual(values = c(`C>A` = "dodgerblue", `C>G` = "black", `C>T` = "red", 
        `T>A` = "grey70", `T>C` = "olivedrab3", `T>G` = "plum2")) +
      scale_y_continuous(expand = c(0,0)) +
      ggtitle(unique(df$signature)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(family = "mono"),
            legend.position = "none")
  })
dev.off()

# do pairwise cosine similarity between all ref sigs
mat <-
  ref %>%
  select(-Type) %>%
  as.matrix()
cos_sim <- lsa::cosine(mat)
colnames(cos_sim) <- colnames(mat)
rownames(cos_sim) <- colnames(mat)

# plot cosine similarity heatmaps
cos_sim_df <-
  tibble::as_tibble(cos_sim, rownames = "sig1") %>%
  pivot_longer(-sig1, names_to = "sig2", values_to = "cosine_similarity")

pdf(file.path(out_dir, "ref_sigs_cosine_similarity.pdf"),
    width = 13, height = 12)
cos_sim_df %>%
  ggplot(aes(sig1, sig2, fill = cosine_similarity)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank()
  )
dev.off()

pdf(file.path(out_dir, "custom_ref_sigs_cosine_similarity.pdf"),
    width = 6, height = 5)
cos_sim_df %>%
  filter(sig1 %in% colnames(custom_ref), sig2 %in% colnames(custom_ref)) %>%
  ggplot(aes(sig1, sig2, fill = cosine_similarity)) +
  geom_tile() +
  geom_text(aes(label = round(cosine_similarity, 3)), color = "black", size = 3) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank(), legend.position = "none")
dev.off()