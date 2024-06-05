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
write.tree(consensus_tree1111, here(here("output/by_level/loom1111/consensus_tree1111.tree")))


# Weigth 8421 : Burnin 90%
path8421 <- here("output/by_level/loom8421/kd_loom8421.trees")
trees8421 <- read.nexus(path8421)
trees8421 <- remove_burnin(trees8421, 0.9)
consensus_tree8421 <- myconsensus(trees8421)
plot(consensus_tree8421, main = "8421")
write.tree(consensus_tree1111, here(here("output/by_level/loom8421/consensus_tree8421.tree")))

# Weigth 1000 : Burnin 90%
path1000 <- here("output/by_level/loom1000/kd_loom1000.trees")
trees1000 <- read.nexus(path1000)
trees1000 <- remove_burnin(trees1000, 0.9)
consensus_tree1000 <- myconsensus(trees1000)
plot(consensus_tree1000, main = "1000")
write.tree(consensus_tree1000, here(here("output/by_level/loom1000/consensus_tree1000.tree")))



# Plot mirror
obj <- cophylo(
  consensus_tree1111,
  consensus_tree1000,
  rotate = T
)

plot(obj,
  link.type = "curved",
  link.lwd = 2,
  link.lty = "solid",
  link.col = "darkblue", cex = 0.01
)
