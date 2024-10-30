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

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  select(group, lng, lng_group_code) |>
  filter(!is.na(lng_group_code)) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |>
  mutate(group = str_replace_all(group, " ", "_"))
kd_lgs <- read_csv(here("data/kd-lgs/kd-lgs_datapoints.csv")) |>
  select(lng, lng_group_code) |>
  filter(!is.na(lng_group_code)) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng))
kd_lgs_looms <- inner_join(kd_looms, kd_lgs) 
  #mutate(across(where(is.character), ~ str_replace_all(., "Chiangmai", "Chengmai"))) |>
  #mutate(lng = if_else(row_number() == 10, "BeChengmaiMY", lng))


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


# ------------- Trees -----------------------
# Languages
phylo_lg = readNexus(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"))
phylo_lg = remove_burnin(phylo_lg,0.9)


# Looms
phylo_looms = readNexus(here("data/kd-looms/kd-looms_bcov1111_strict/kd-looms_bcov1111_strict.trees"))
phylo_looms = remove_burnin(phylo_looms, 0.9)


phylo_lg2 = phylo_lg[[1]]
phylo_looms2 = phylo_looms[[1]]


# Common taxa
tip_to_keep_lgs = intersect(phylo_lg2$tip.label,kd_lgs_looms$lng_label)
tip_to_keep_looms = intersect(phylo_looms2$tip.label,kd_lgs_looms$group)

keep.tip(phylo_lg[[1]],tip_to_keep_lgs)

for (i in 2:length(phylo_lg)){
  phylo_lg2[[i]] = keep.tip(phylo_lg[[i]],tip_to_keep_lgs)
}

for (i in 2:length(phylo_looms)){
  phylo_looms2[[i]] = keep.tip(phylo_looms[[i]],tip_to_keep_looms)
}

# Test of both trees elaged
#par(mfrow=c(1,2))
#plot(phylo_lg2[[10]])
#plot(phylo_looms2[[7]])


# Nexus -----------------------------------------------------------------------
# Looms
data_looms = read.nexus.data(here("data/nexus/kd-looms_pruned_filtered.nex")) |>
  as.data.frame() |>
  t() 
mode(data_looms) = "integer"
data_looms = as.data.frame(data_looms)

data_looms <- data_looms[match(phylo_looms2[[2]]$tip.label, rownames(data_looms)), ]


# Languages
data_lg = read.nexus.data(here("data/nexus/kd-lngs_pruned_filtered.nex")) |>
  as.data.frame() |>
  t() 
mode(data_lg) = "integer"
data_lg = as.data.frame(data_lg)

# Ensure the species in the data match the tree's tip labels
data_lg <- data_lg[match(phylo_lg2[[2]]$tip.label, rownames(data_lg)), ]


# ----------- Perform PCA on binary data ----------
# Name of looms become languages
rownames(data_looms) = rownames(data_looms) |>
  map_chr(~ {
    match <- which(str_detect(kd_lgs_looms$group, fixed(.x)))
    
    if (length(match) > 0) {
      kd_lgs_looms$lng_label[match]
    } else {
      .x
    }
  })

#looms
# Using FactoMineR, scale.unit = FALSE prevents scaling of binary variables
pca_result_looms <- PCA(data_looms, graph = FALSE, scale.unit = FALSE)
pc1_looms <- pca_result_looms$ind$coord[,1]



# languages 
# Perform PCA on binary data
# Using FactoMineR, scale.unit = FALSE prevents scaling of binary variables
pca_result_lg <- PCA(data_lg, graph = FALSE, scale.unit = FALSE)
pc1_lg <- pca_result_lg$ind$coord[,1]



# Blomberg's K for pruned trees

# Data lg / Tree looms
phylo_looms3 = c()

for (i in 1:length(phylo_looms)){
  match <- which(str_detect(kd_lgs_looms$group,phylo_looms[[i]]$tip.label))
  phylo_looms3[[1]]$edge = phylo_looms[[1]]$edge
  phylo_looms3[[1]]$edge.length = phylo_looms[[1]]$edge.length 
  phylo_looms3[[1]]$Nnode = phylo_looms[[1]]$Nnode
  phylo_looms3[[1]]$root.edge = phylo_looms[[1]]$root.edge
  phylo_looms3[[i]]$tip.label = kd_lgs_looms$lng_label[match]
  }

i_lg_loom = matrix(NA, ncol=2, nrow= length(phylo_looms2))
for (i in 2:length(phylo_looms2)){
  res <- phylosig(phylo_looms2[[i]], pc1_lg, method = "K", test = TRUE)
  i_lg_loom[i,] <- c(res$K,res$P)
}

mean(i_lg_loom[,1])
mean(i_lg_loom[,2])

# Data looms / Tree lg
i_loom_lg = matrix(NA, ncol=2, nrow= length(phylo_lg2))

for (i in 2:length(phylo_lg2)){
  res <- phylosig(phylo_lg2[[i]], pc1_looms, method = "K", test = TRUE)
  i_loom_lg[i,] <- c(res$K,res$P)
}

mean(i_loom_lg[,1], na.rm=T)
mean(i_loom_lg[,2], na.rm=T)

#res$sim.K









