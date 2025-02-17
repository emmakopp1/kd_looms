library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(FactoMineR)
library(tidyverse)
library(ape)
library(castor)

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

<<<<<<< HEAD

=======
>>>>>>> refs/remotes/origin/main
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
  write_csv("output/data/kd-k_summary.csv")


# correlation test according to brown 2017 
# methodologie : average on all phylogenies
kd_lg_tree <- read.tree(here('output/trees/kd-lgs_pruned.trees'))
M_lg <- length(kd_lg_tree)
kd_lg_tree <- kd_lg_tree[seq(0.2* M_lg,M_lg, length = 1000)]
M_lg <- length(kd_lg_tree)

#kd_lg_tree <- kd_pruned_consensus_lg # obtained from visualise.R
N <- Ntip(kd_lg_tree)[1]
sorted_tips <- sort(kd_lg_tree[[1]]$tip.label)

# pariwise mean distance
lg_distance_matrix <- outer(
  sorted_tips, 
  sorted_tips, 
  Vectorize(function(tip_i, tip_j) {
    distances <- sapply(kd_lg_tree, function(tree) {
      get_pairwise_distances(tree, tip_i, tip_j)/max(node.depth.edgelength(tree))
    })
    
    # Calculer la moyenne des distances sur tous les arbres
    mean(distances)
  })
)

lg_distance_matrix <- round(lg_distance_matrix,3)

# looms 
kd_looms_tree <- read.tree(here('output/trees/kd-looms_prunedk.trees'))
M_looms <- length(kd_looms_tree)
kd_looms_tree <- kd_looms_tree[seq(0.2* M_looms, M_looms, length = 1000)]
M_looms <- length(kd_looms_tree)


looms_distance_matrix <- outer(
  sorted_tips, 
  sorted_tips, 
  Vectorize(function(tip_i, tip_j) {
    distances <- sapply(kd_looms_tree, function(tree) {
      get_pairwise_distances(tree, tip_i, tip_j)/max(node.depth.edgelength(tree))
    })
    
    # mean over all trees
    mean(distances)
  })
)

looms_distance_matrix <- round(looms_distance_matrix,3)

# mantel test
mantel.test(looms_distance_matrix, lg_distance_matrix, nperm = 1e4)


# methodology : on consensus trees 
kd_lg_tree_cs <- read.tree(here("output/trees/kd-lgs_bcov_relaxed_ht_pos_pruned_cs.tree"))
N <- Ntip(kd_lg_tree_cs)[1]
sorted_tips <- sort(kd_lg_tree_cs$tip.label)

lg_distance_matrix_cs <- outer(
  sorted_tips, 
  sorted_tips, 
  Vectorize(function(tip_i, tip_j) {
    get_pairwise_distances(kd_lg_tree_cs, tip_i, tip_j)/max(node.depth.edgelength(kd_lg_tree_cs))
  })
)
lg_distance_matrix_cs <- round(lg_distance_matrix_cs,3)

# looms 
kd_looms_tree_cs <- read.tree(here("output/trees/kd-looms_bcov1111_strict_ht_pruned_cs.tree"))

looms_distance_matrix_cs <- outer(
  sorted_tips, 
  sorted_tips, 
  Vectorize(function(tip_i, tip_j) {
    get_pairwise_distances(kd_looms_tree_cs, tip_i, tip_j)/max(node.depth.edgelength(kd_looms_tree_cs))
  })
)
looms_distance_matrix_cs <- round(looms_distance_matrix_cs,3)

# mantel test
mantel.test(looms_distance_matrix_cs, lg_distance_matrix_cs, nperm = 1e4)













