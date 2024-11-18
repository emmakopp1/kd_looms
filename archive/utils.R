rm(list = ls())
library(ggplot2)
library(here)
library(ape)
library(phytools)
library(TreeTools)
library(stats)
library(readxl)
library(purrr)
library(phangorn)
library(readODS)
library(tidyverse)
library(tracerer)

# ---------------------------   Inter group distances ---------------------------

# Import
lgs_to_looms = read_excel(
  here("/Users/kopp/Documents/kd_looms/data/Kra-DaiLooms.xlsx"),
  skip = 5,
  sheet=1) |> 
  select(Looms = 2, KD = 3)
  

# Function
remove_burnin <- function(trees, burnin_rate) {
  n <- as.numeric(length(trees))
  return(trees[as.integer(n * burnin_rate + 1):n])
}

# Change tip labels of the looms according to the corresponding to the closest language in kd-loom.ods file
change_tip_label <- function(tree) {
  tree$tip.label |>
    map_chr(~ {
      match <- which(str_detect(lgs_to_looms$KD, fixed(.x)))

      if (length(match) > 0) {
        lgs_to_looms$Looms[match]
      } else {
        .x
      }
    })
}

# Pruning
prune_trees <- function(trees) {
  common_leaves <- intersect(trees[[2]]$tip.label,trees[[1]]$tip.label)
  return(map(trees, ~ drop.tip(.x, setdiff(.x$tip.label, common_leaves))))
}

# Inter group distance
compute_trees <- function(){
  list(
    here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"),
    here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees")
  ) |>
    map(~ {
      arbres <- readNexus(.x) |>
        remove_burnin(burnin_rate = 0.9)})
}


compute_inter_rf_distance <- function(phylo) {
  phylo |>
    map(~ sample(.x, size = 1)) |>
    map(~ {
      .x[[names(.x)[1]]]
    }) |>
    modify_at(1, ~ {
      .x$tip.label <- change_tip_label(.x)
      .x
    }) |> 
    modify_at(2, ~ {
      .x$tip.label <- str_remove(.x$tip.label,'_')
      .x
    }) |>
    prune_trees() |>
    reduce(RF.dist,normalize = TRUE, rooted = T, check.labels = T)
}

# Répéter l'expérience 100 fois et calculer la moyenne
phylo_lg_loom = compute_trees()
inter_rf_distances <- replicate(1000, {
  phylo_lg_loom |> compute_inter_rf_distance()
  }) 

inter_rf_distances |> median()

# For two trees
trees = list(
  here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"),
  here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees")
) |>
  map(~ {
    arbres <- readNexus(.x) |>
      remove_burnin(burnin_rate = 0.9) |>
      sample(size = 1)
  }) |>
  map(~ {
    .x[[names(.x)[1]]]
  }) |>
  modify_at(1, ~ {
    .x$tip.label <- change_tip_label(.x)
    .x
  }) |> 
  modify_at(2, ~ {
    .x$tip.label <-str_remove(.x$tip.label,'_')
    .x
  }) |>
  prune_trees() 


par(mfrow=c(1,2))
plot(trees[[2]])
plot(trees[[1]])


RF.dist(trees[[1]],trees[[2]])


# ---------------------------   Intra group distances ---------------------------
# For languages 
path_lg <- here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees")

# Lire et prétraiter les arbres une seule fois
trees_lg <- path_lg |>
  readNexus() |>
  remove_burnin(burnin_rate = 0.9)

# Calculer les distances RF intra-groupe pour 100 paires aléatoires
intra_rf_distances_lg <- replicate(1000, {
  trees_lg |> 
  sample(size = 2, replace = FALSE) |> 
  RF.dist(normalize = TRUE, rooted = T, check.labels = T) 
  })


# For looms
path_looms = here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees")


trees_looms <- path_looms |>
  readNexus() |>
  remove_burnin(burnin_rate = 0.9)

# Calculer les distances RF intra-groupe pour 100 paires aléatoires
intra_rf_distances_looms <- replicate(1000, {
  trees_looms |> 
    sample(size = 2, replace = FALSE) |> 
    RF.dist(normalize = TRUE, rooted = T, check.labels = T) 
})

intra_rf_distances_looms |> median()
#intra_rf_distances_looms |> var()

par(mfrow=c(1,3))
hist(intra_rf_distances_looms, breaks = 10, xlim = c(0,1))
hist(intra_rf_distances_lg, breaks = 10, xlim = c(0,1))
hist(inter_rf_distances, breaks = 10, xlim = c(0,1))


# For merged nexus files, delete the full 0 columns

library(ape)
library(tidyverse)
library(here)

# Registration process : delete nulls columns --------------------------------------------------------------------
# Languages 
read.nexus.data(here("data/kd-lgs/kd-lgs_elaged/kd-lg-elaged.nex"))|>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "-"))) |>
  mutate(across(everything(), as.numeric))  |>
  mutate(total = rowSums(across(everything()), na.rm=T)) |>
  filter(total != 0) |> 
  select(-total) |>
  mutate(across(everything(), as.character))  |>
  mutate(across(everything(), ~ replace_na(.x, "?"))) |>
  write.nexus.data(
    here("data/kd-lgs/kd-lgs_elaged/kd-lg-elaged-filtered.nex"),
    format = "STANDARD")


# Looms
read.nexus.data(here("data/kd-looms/kd-looms_elaged/kd-looms_elaged-vf.nex"))|>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "-"))) |>
  mutate(across(everything(), as.numeric))  |>
  mutate(total = rowSums(across(everything()), na.rm=T)) |>
  filter(total != 0) |> 
  select(-total) |>
  mutate(across(everything(), as.character))  |>
  mutate(across(everything(), ~ replace_na(.x, "?"))) |>
  write.nexus.data(
    here("data/kd-looms/kd-looms_elaged/kd-looms_elaged.nex"),
    format = "STANDARD")



# Merged data
read.nexus.data(here("data/kd-merged/merged_vf.nex"))|>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "-"))) |>
  mutate(across(everything(), as.numeric))  |>
  mutate(total = rowSums(across(everything()), na.rm=T)) |>
  filter(total != 0) |> 
  select(-total) |>
  mutate(across(everything(), as.character))  |>
  mutate(across(everything(), ~ replace_na(.x, "?"))) |>
  write.nexus.data(
    here("data/kd-merged/merged.nex"),
    format = "STANDARD")





