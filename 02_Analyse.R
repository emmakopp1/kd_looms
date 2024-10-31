library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(FactoMineR)
library(tidyverse)

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |>
  select(group, lng_label)

kd_lgs_phylo <- read.tree(here("output/trees/kd-lgs_prunedk.trees"))
kd_looms_phylo <- read.tree(here("output/trees/kd-looms_prunedk.trees"))

pca_looms <- ReadAsPhyDat(here("data/nexus/kd-looms_pruned.nex")) |>
  as_tibble() |>
  rename_with(~ str_replace_all(.x, "_", " ")) |>
  rename(any_of(setNames(kd_looms$group, kd_looms$lng_label))) |>
  mutate(across(everything(), ~ as.numeric(.x))) |>
  as.matrix() |>
  t() |>
  PCA(graph = FALSE, scale.unit = FALSE)
pc1_looms <- pca_looms$ind$coord[, 1]

pca_lgs <- ReadAsPhyDat(here("data/nexus/kd-lgs_pruned.nex")) |>
  as_tibble() |>
  mutate(across(everything(), ~ as.numeric(.x))) |>
  as.matrix() |>
  t() |>
  PCA(graph = FALSE, scale.unit = FALSE)
pc1_lgs <- pca_lgs$ind$coord[, 1]

set.seed(973829350)
kd_lgs_on_looms_k <- map_df(1:length(kd_looms_phylo), function(i) {
  ksig <- phylosig(kd_looms_phylo[[i]], pc1_lgs, method = "K", test = TRUE)
  tibble(k = ksig$K, p = ksig$P)
})
kd_lgs_on_looms_k |>
  write_csv(here("output/data/kd-lgs_on_looms_k.csv"))

kd_looms_on_lgs_k <- map_df(1:length(kd_lgs_phylo), function(i) {
  ksig <- phylosig(kd_lgs_phylo[[i]], pc1_looms, method = "K", test = TRUE)
  tibble(k = ksig$K, p = ksig$P)
})
kd_looms_on_lgs_k |>
  write_csv(here("output/data/kd-looms_on_lgs_k.csv"))

kd_lgs_on_looms_k |>
  summarise(across(everything(), mean))
kd_looms_on_lgs_k |>
  summarise(across(everything(), mean))