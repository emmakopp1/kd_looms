library(ape)          # For working with phylogenetic trees
library(phytools)     # For calculating phylogenetic signal
library(geiger)       # For additional phylogenetic tools
library(FactoMineR)   # For performing PCA
library(here)
library(readxl)
library(rwty)
library(TreeTools)
library(readr)
library(stringr)
library(dplyr)
library(purrr)

# Import data ------------------------------------------------------------------

lgs_to_looms = read_excel(
  here("data/Kra-DaiLooms.xlsx"),
  skip = 5,
  sheet=1) |> 
  select(Looms = 2, KD = 3)


# Functions --------------------------------------------------------------------

# Burnin
remove_burnin <- function(trees, burnin_rate) {
  n <- as.numeric(length(trees))
  return(trees[as.integer(n * burnin_rate + 1):n])
}

# Consensus tree 
myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}

# Pruning
prune_trees <- function(trees) {
  common_leaves <- intersect(trees[[2]]$tip.label,trees[[1]]$tip.label)
  return(map(trees, ~ drop.tip(.x, setdiff(.x$tip.label, common_leaves))))
}


# Change tip labels of the looms according into the corresponding language
change_tip_label <- function(tiplabel){
  L = length(tiplabel)
  for (l in 1:L){
    match <- which(str_detect(lgs_to_looms$KD, tiplabel[l]))
    if (length(match) > 0) {
      tiplabel[l] <- lgs_to_looms$Looms[match]
    }}
  return(tiplabel)
  }
  

# ------------- Trees -----------------------
# Languages
phylo_lg = readNexus(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"))
phylo_lg = remove_burnin(phylo_lg,0.9)


# Looms
phylo_looms = readNexus(here("data/kd-looms/kd-looms_bcov1111_strict/kd-looms_bcov1111_strict.trees"))
phylo_looms = remove_burnin(phylo_looms, 0.9)


phylo_lg2 = phylo_lg[[1]]
phylo_looms2 = phylo_looms[[1]]


# Change tips labels
# Languages
for (i in 1:length(phylo_lg)) {
  tree <- phylo_lg[[i]] 
  new_labels <- change_tip_label(tree$tip.label)
  tree$tip.label = new_labels
  phylo_lg2[[i]] <- tree
}

# Looms
for (i in 1:length(phylo_looms)) {
  tree <- phylo_looms[[i]] 
  new_labels <- tree |> TipLabels() |> str_remove_all("_")
  tree$tip.label = new_labels
  phylo_looms2[[i]] <- tree
}

# Common taxa
common_leaves <- intersect(phylo_looms2[[1]]$tip.label, phylo_lg2[[1]]$tip.label)

for (i in 1:length(phylo_lg)){
  phylo_lg2[[i]] = keep.tip(phylo_lg2[[i]],common_leaves)
}

for (i in 1:length(phylo_looms)){
  phylo_looms2[[i]] = keep.tip(phylo_looms2[[i]],common_leaves)
}

# Test of both trees elaged
#par(mfrow=c(1,2))
#plot(phylo_lg2[[5]])
#plot(phylo_looms2[[7]])


# Nexus -----------------------------------------------------------------------
# Looms
data_looms = read.nexus.data(here("data/nexus/kd-looms_1111.nex"))
data_looms = t(data.frame(data_looms))
mode(data_looms) = "integer"
data_looms = as.data.frame(data_looms)
rownames(data_looms) = str_remove(rownames(data_looms),'_')

# Ensure the species in the data match the tree's tip labels
data_looms <- data_looms[match(phylo_looms2[[1]]$tip.label, rownames(data_looms)), ]

# languages
# Compute consensus tree and load it 
#trees_lg <- readNexus(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"))
#trees_lg = remove_burnin(trees_lg, 0.9)
#tree_lg = myconsensus(trees_lg)
#write.tree(tree_lg,here("/Users/kopp/Documents/kd_looms/output/ConsensusTree/lg_bcov.tree"))

data_lg = read.nexus.data(here("data/nexus/kd-lgs.nex"))
data_lg = t(data.frame(data_lg))
mode(data_lg) = "integer"
data_lg = as.data.frame(data_lg)

rownames(data_lg) = rownames(data_lg) |>
  map_chr(~ {
    match <- which(str_detect(lgs_to_looms$KD, fixed(.x)))
    
    if (length(match) > 0) {
      lgs_to_looms$Looms[match]
    } else {
      .x
    }
  })
# Ensure the species in the data match the tree's tip labels
data_lg <- data_lg[match(phylo_lg2[[1]]$tip.label, rownames(data_lg)), ]


# ----------- Perform PCA on binary data ----------
#looms
# Using FactoMineR, scale.unit = FALSE prevents scaling of binary variables
pca_result_looms <- PCA(data_looms, graph = FALSE, scale.unit = FALSE)
pc1_looms <- pca_result_looms$ind$coord[,1]



# languages 
# Perform PCA on binary data
# Using FactoMineR, scale.unit = FALSE prevents scaling of binary variables
pca_result_lg <- PCA(data_lg, graph = FALSE, scale.unit = FALSE)
pc1_lg <- pca_result_lg$ind$coord[,1]



# Blomberg's K for elaged trees

# Data looms / Tree looms
i_loom_loom = matrix(NA, ncol=2, nrow= length(phylo_looms2))

for (i in 1:length(phylo_looms2)){
  res <- phylosig(phylo_looms2[[i]], pc1_looms, method = "K", test = TRUE)
  i_loom_loom[i,] <- c(res$K,res$P)
}

mean(i_loom_loom[,1])
mean(i_loom_loom[,2])

# Data lg / Tree lg
i_lg_lg = matrix(NA, ncol=2, nrow= length(phylo_lg2))
for (i in 1:length(phylo_lg2)){
  res <- phylosig(phylo_lg2[[i]], pc1_lg, method = "K", test = TRUE)
  i_lg_lg[i,] <- c(res$K,res$P)
}

mean(i_lg_lg[,1])
mean(i_lg_lg[,2])


# Data lg / Tree looms
i_lg_loom = matrix(NA, ncol=2, nrow= length(phylo_looms2))
for (i in 1:length(phylo_looms2)){
  res <- phylosig(phylo_looms2[[i]], pc1_lg, method = "K", test = TRUE)
  i_lg_loom[i,] <- c(res$K,res$P)
}

mean(i_lg_loom[,1])
mean(i_lg_loom[,2])

# Data looms / Tree lg
i_loom_lg = matrix(NA, ncol=2, nrow= length(phylo_lg2))

for (i in 1:length(phylo_lg2)){
  res <- phylosig(phylo_lg2[[i]], pc1_looms, method = "K", test = TRUE)
  i_loom_lg[i,] <- c(res$K,res$P)
}

mean(i_loom_lg[,1])
mean(i_loom_lg[,2])


res$sim.K










