library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(FactoMineR)
library(tidyverse)
library(ggtree)
library(castor)

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |>
  select(group, lng_label)

kd_lgs_phylo <- read.tree(here("output/trees/kd-lgs_prunedk.trees"))
kd_looms_phylo <- read.tree(here("output/trees/kd-looms_prunedk.trees"))


# Prepare PCA components to compute K -------------------------------------

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


# Compute Blomberg's K ----------------------------------------------------

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
  mutate(data = "lgs on looms", .before = everything()) |>
  bind_rows(mutate(kd_looms_on_lgs_k, data = "looms on lgs", .before = everything())) |>
  group_by(data) |>
  summarise(across(everything(), mean)) |>
  write_csv(here("output/data/kd-k_summary.csv"))


# Correlation test following Brown (2017) ---------------------------------

set.seed(2340706)

## Average
lg_dist_ms <- map(
  kd_lgs_phylo[seq(1, length(kd_lgs_phylo), length.out = 1000)],
  ~ get_all_pairwise_distances(.x,
    only_clades = sort(.x$tip.label)
  ) / max(node.depth.edgelength(.x))
)
lg_dist_m <- Reduce("+", lg_dist_ms) / length(lg_dist_ms)

looms_dist_ms <- map(
  kd_looms_phylo[seq(1, length(kd_looms_phylo), length.out = 1000)],
  ~ get_all_pairwise_distances(.x,
    only_clades = sort(.x$tip.label)
  ) / max(node.depth.edgelength(.x))
)
looms_dist_m <- Reduce("+", looms_dist_ms) / length(looms_dist_ms)

mt_mean <- mantel.test(lg_dist_m, looms_dist_m, nperm = 1e4)

## On consensus trees
kd_lgs_pruned_tips <- ReadAsPhyDat(here("data/nexus/kd-lgs_pruned.nex")) |>
  as_tibble() |>
  colnames()
kd_lgs_cs <- read.tree(here("output/trees/kd-lgs_bcov_relaxed_ht_pos_consensus.tree")) |>
  keep.tip(kd_lgs_pruned_tips)
if (!is.rooted(kd_lgs_cs)) {
  kd_lgs_cs$root.edge.length <- 0
}

kd_looms_cs <- read.tree(here("output/trees/kd-looms_bcov1111_relaxed_ht_consensus.tree")) |>
  fortify() |>
  mutate(label = str_replace_all(label, "_", " ")) |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  as.phylo() |>
  keep.tip(kd_lgs_pruned_tips)
if (!is.rooted(kd_looms_cs)) {
  kd_looms_cs$root.edge.length <- 0
}

lg_cs_dist_m <- get_all_pairwise_distances(
  kd_lgs_cs,
  only_clades = sort(kd_lgs_cs$tip.label)
) / max(node.depth.edgelength(kd_lgs_cs))
looms_cs_dist_m <- get_all_pairwise_distances(
  kd_looms_cs,
  only_clades = sort(kd_looms_cs$tip.label)
) / max(node.depth.edgelength(kd_looms_cs))

mt_cs <- mantel.test(lg_cs_dist_m, looms_cs_dist_m, nperm = 1e4)

mt_mean |>
  as_tibble() |>
  mutate(method = "mean") |>
  bind_rows(mt_cs |>
    as_tibble() |>
    mutate(method = "consensus")) |>
  write_csv(here("output/data/kd-mantel.csv"))
