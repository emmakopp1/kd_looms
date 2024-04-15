rm(list=ls())
library(tidyverse)
library(dplyr)
library(readxl)
library(phangorn)
library(TreeTools)
library(here)
setwd("/Users/kopp/Library/CloudStorage/OneDrive-UniversitéParis-Dauphine/these/KD_loom/")

kd_xl <- read_excel("data/Kra-DaiLooms28master-2.xlsx", sheet = "Digit", skip = 1) |>
  filter(Character != "Level") |>
  select(!starts_with("...")) |>
  mutate(Character = str_replace_all(Character, " ", "_"))|>
  column_to_rownames("Character")

simple <- colnames(kd_xl) |>
  str_subset("^FM|HP|RP|RT|RW|MH")
  
complex <- colnames(kd_xl) |>
  str_subset("^PS")

basics <- setdiff(colnames(kd_xl),union(simple,complex))

kd_looms_dt_all <- kd_xl |>
  as.matrix() |>
  MatrixToPhyDat()
write.phyDat(kd_looms_dt_all, "data/kd_looms_dt_all.nex", format = "nexus")

kd_looms_dt_basic <- kd_xl |>
  select(!any_of(c(simple, complex))) |>
  as.matrix() |>
  MatrixToPhyDat()
write.phyDat(kd_looms_dt_basic, "data/kd_looms_dt_basic.nex", format = "nexus")

kd_looms_dt_complex <- kd_xl |>
  select(all_of(complex)) |>
  as.matrix() |>
  MatrixToPhyDat()
write.phyDat(kd_looms_dt_complex, "data/kd_looms_dt_complex.nex", format = "nexus")




