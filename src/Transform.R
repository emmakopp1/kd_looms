library(here)
library(tidyverse)
library(phangorn)
library(phytools)
# library(treeio)

dir.create(here("output/trees"))

tree1000 <- read.nexus(here("data/by_level/loom1000/kd_loom1000.trees"))
tree1000 <- tree1000[seq(2, length(tree1000), by = round(length(tree1000) / 1000))]
cs_tree1000 <- consensus(tree1000, p = .5, rooted = TRUE)
cs_tree1000 <- consensus.edges(tree1000, consensus.tree = cs_tree1000, rooted = TRUE)
write.tree(cs_tree1000, here("output/trees/kd_loom1000_consensus.tree"))

tree1111 <- read.nexus(here("data/by_level/loom1111/kd_loom1111.trees"))
tree1111 <- tree1111[seq(2, length(tree1111), by = round(length(tree1111) / 1000))]
cs_tree1111 <- consensus(tree1111, p = .5, rooted = TRUE)
cs_tree1111 <- consensus.edges(tree1111, consensus.tree = cs_tree1111, rooted = TRUE)
write.tree(cs_tree1111, here("output/trees/kd_loom1111_consensus.tree"))

tree8421 <- read.nexus(here("data/by_level/loom8421/kd_loom8421.trees"))
tree8421 <- tree8421[seq(2, length(tree8421), by = round(length(tree8421) / 1000))]
cs_tree8421 <- consensus(tree8421, p = .5, rooted = TRUE)
cs_tree8421 <- consensus.edges(tree8421, consensus.tree = cs_tree8421, rooted = TRUE)
write.tree(cs_tree8421, here("output/trees/kd_loom8421_consensus.tree"))
