library(here)
library(ape)
library(phytools)
library(TreeTools)

remove_burnin <- function(trees, burnin_rate) {
  n <- as.numeric(length(trees))
  return(trees[as.integer(n * burnin_rate):n])
}
myconsensus <- function(trees) {
  consensus.edges(trees, consensus.tree = consensus(trees, p = .5), rooted = T)
}

# Weigth 1111 : Burnin 90%
path1111 <- here("output/by_level/loom1111/kd_loom1111.trees")
trees1111 <- read.nexus(path1111)
trees1111 <- remove_burnin(trees1111, 0.9)
consensus_tree1111 <- myconsensus(trees1111)
plot(consensus_tree1111, main = "1111")
write.tree(consensus_tree1111, here("output/by_level/loom_bcov_1111/loom_bcov_1111_consensus.tree"))

# Read the consensus tree
consensus_tree1111 = read.tree(here("output/by_level/loom_bcov_1111/loom_bcov_1111_consensus.tree"))

# Weigth 8421 : Burnin 90%
path8421 <- here("output/by_level/loom8421/kd_loom8421.trees")
trees8421 <- read.nexus(path8421)
trees8421 <- remove_burnin(trees8421, 0.9)
consensus_tree8421 <- myconsensus(trees8421)
plot(consensus_tree8421, main = "8421")
write.tree(consensus_tree1111, here("output/by_level/loom8421/consensus_tree8421.tree"))

# Weigth 1000 : Burnin 90%
path1000 <- here("output/by_level/loom1000/kd_loom1000.trees")
trees1000 <- read.nexus(path1000)
trees1000 <- remove_burnin(trees1000, 0.9)
consensus_tree1000 <- myconsensus(trees1000)
plot(consensus_tree1000, main = "1000")
write.tree(consensus_tree1000, here("output/by_level/loom1000/consensus_tree1000.tree"))


# Data with 4 sites 
path_4sites <- here("output/by_level/loom_ctmc_variable_trait/loom_ctmc_variable_trait.trees")
trees_4sites <- read.nexus(path_4sites)
trees_4sites <- remove_burnin(trees_4sites, 0.9)
consensus_4sites <- myconsensus(trees_4sites)
plot(consensus_4sites, main = "4 rates sites")
write.tree(consensus_4sites, here("output/by_level/loom_ctmc_variable_trait/consensus_tree_4sites.tree"))


plot.phylo(consensus_4sites, type = "phylogram", edge.width = 2, cex = 0.6)



# Plot mirror
obj <- cophylo(
  consensus_tree1111,
  consensus_4sites,
  rotate = T
)

plot(obj,
  link.type = "curved",
  link.lwd = 2,
  link.lty = "solid",
  link.col = "darkblue", cex = 1
)

par(mfrow=c(1,2))
plot.phylo(consensus_tree1111,consensus_4sites)

