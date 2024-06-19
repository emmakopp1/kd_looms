library(here)
library(tidyverse)
library(phangorn)
library(phytools)
library(tracerer)

dir.create(here("output/trees"))
ntrees <- 1000


# Consensus trees for looms ---------------------------------------------------------------------------------------

tree1000 <- read.nexus(here("data/by_level/loom1000/kd_loom1000.trees"))
tree1000 <- tree1000[seq(2, length(tree1000), by = round(length(tree1000) / ntrees))]
cs_tree1000 <- consensus(tree1000, p = .5, rooted = TRUE)
cs_tree1000 <- consensus.edges(tree1000, consensus.tree = cs_tree1000, rooted = TRUE)
write.tree(cs_tree1000, here("output/trees/kd_loom1000_consensus.tree"))

tree1111 <- read.nexus(here("data/by_level/loom1111/kd_loom1111.trees"))
tree1111 <- tree1111[seq(2, length(tree1111), by = round(length(tree1111) / ntrees))]
cs_tree1111 <- consensus(tree1111, p = .5, rooted = TRUE)
cs_tree1111 <- consensus.edges(tree1111, consensus.tree = cs_tree1111, rooted = TRUE)
write.tree(cs_tree1111, here("output/trees/kd_loom1111_consensus.tree"))

tree8421 <- read.nexus(here("data/by_level/loom8421/kd_loom8421.trees"))
tree8421 <- tree8421[seq(2, length(tree8421), by = round(length(tree8421) / ntrees))]
cs_tree8421 <- consensus(tree8421, p = .5, rooted = TRUE)
cs_tree8421 <- consensus.edges(tree8421, consensus.tree = cs_tree8421, rooted = TRUE)
write.tree(cs_tree8421, here("output/trees/kd_loom8421_consensus.tree"))


# Consensus tree for languages ------------------------------------------------------------------------------------

lgs_trees <- read.nexus(here("data/languages/kd2.trees"))
lgs_trees <- lgs_trees[seq(2, length(lgs_trees), by = round(length(lgs_trees) / ntrees))]
kd_lgs_cs <- consensus(lgs_trees, p = .5, rooted = TRUE)
kd_lgs_cs <- consensus.edges(lgs_trees, consensus.tree = kd_lgs_cs, rooted = TRUE)
write.tree(kd_lgs_cs, here("output/trees/kd_lgs_consensus.tree"))


# Ages in languages trees  ----------------------------------------------------------------------------------------

getMRCA_age <- function(tree, tips) {
  tips <- if (is.character(tips)) which(tree$tip.label %in% tips) else tips
  mrca <- ifelse(length(tips) > 1, getMRCA(tree, tips), tips)
  root_age <- max(node.depth.edgelength(tree))
  root_age - node.depth.edgelength(tree)[mrca]
}

lgs_trees[[1]]$tip.label
getMRCA_age(lgs_trees, lgs_trees[[1]]$tip.label)
kd_lgs_ages <- lgs_trees |>
  seq_along() |>
  map_df(~ tibble(
    `Kra-Dai` = getMRCA_age(lgs_trees[[.x]], lgs_trees[[1]]$tip.label),
    `Kam-Tai` = getMRCA_age(lgs_trees[[.x]], str_subset(lgs_trees[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")),
    `Tai-Yay` = getMRCA_age(lgs_trees[[.x]], str_subset(lgs_trees[[1]]$tip.label, "^(Tc|Tn|Tsw)"))
  )) |>
  pivot_longer(everything(), names_to = "group", values_to = "age")
write_csv(kd_lgs_ages, here("output/kd_lgs_ages.csv"))


# Mutation rates --------------------------------------------------------------------------------------------------

mutationrate_bylevel <- parse_beast_tracelog_file(here("data/beast/loom_ctmc_variable_trait/loom_ctmc_variable_trait.log")) |>
  as_tibble() |>
  select(mutationRate.s.level1, mutationRate.s.level2, mutationRate.s.level3, mutationRate.s.level4) |>
  rowid_to_column()

mutationrate_bylevel_tb <- mutationrate_bylevel |>
  filter(!(rowid %in% 1:((nrow(mutationrate_bylevel) - 1) * .1))) |>
  pivot_longer(-rowid, names_to = "level", values_to = "rate") |>
  mutate(level = str_remove_all(level, "[^0-9]"))
