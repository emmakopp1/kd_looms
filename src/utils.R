library(here)
source(here("code/init.R"))

remove_burnin = function(trees,burnin_rate){
  #'remove_burnin
  #'
  #'Remove the burn-in of a sample of trees
  #'@param trees
  #'@param burnin_rate the rate of burn-in we apply
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}
myconsensus = function(trees){
  #' myconcensus
  #' 
  #' Compute the consensus tree of multiples trees
  #' @param trees The trees.
  #' @example 
  #' path = "/Users/kopp/Documents/chr_paper/beast/bantu-ctmc-strict-bd/ctmc-strict-bd.trees"
  #' tree = read.nexus(path)
  #' tree = remove_burnin(tree,0.9)
  #' consensus_tree = myconsensus(tree)
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}

library(ape)

path1 = "/Users/kopp/Documents/kd_looms/output/loom1000/kd_looms_1000.trees"
tree1 = read.nexus(path1)
tree1 = remove_burnin(tree1,0.4)
consensus_tree1 = myconsensus(tree1)
#plot(consensus_tree1)


path2 = "/Users/kopp/Documents/kd_looms/output/loom4100/kd_looms_4100.trees"
tree2 = read.nexus(path2)
tree2 = remove_burnin(tree2,0.9)
consensus_tree2 = myconsensus(tree2)
#plot(consensus_tree2)


# Plot mirror 
obj = cophylo(
  consensus_tree1,
  consensus_tree2,
  rotate = T)

plot(obj,
     link.type="curved",
     link.lwd=2,
     link.lty="solid",
     link.col="darkblue",cex=0.01)


# write files 
write.tree(consensus_tree1,"/Users/kopp/Documents/kd_looms/output/loom1000/consensus.tree")
write.tree(consensus_tree2,"/Users/kopp/Documents/kd_looms/output/loom4100/consensus.tree")
