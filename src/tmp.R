library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(FactoMineR)
library(tidyverse)

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |> 
  select(group, lng_label)

# pca_looms <- 
ReadAsPhyDat(here("data/nexus/kd-looms_pruned.nex")) |> 
  as_tibble() |>
  rename_with(~ str_replace_all(.x, "_", " ")) |>
  rename(any_of(setNames(kd_looms$group, kd_looms$lng_label))) |>
  mutate(across(everything(), ~ as.numeric(.x))) |>
  as.matrix() |>
  t() |> 
  PCA(graph = FALSE, scale.unit = FALSE)
pc1_looms <- pca_looms$ind$coord[,1]

pca_lgs <- ReadAsPhyDat(here("data/nexus/kd-lgs_pruned.nex")) |> 
  as_tibble() |>
  mutate(across(everything(), ~ as.numeric(.x))) |> 
  as.matrix() |>
  t() |>
  PCA(graph = FALSE, scale.unit = FALSE)
pc1_lgs <- pca_lgs$ind$coord[,1]
