# Import file
library(here)
library(stringr)
library(tidyverse)
library(readxl)
library(phangorn)
library(TreeTools)
library(dplyr)
library(ape)
library(phytools)

source(here("src/functions.R"))

# Init
path <- here("data/Kra-DaiLooms28master-2.xlsx")
weight_level <- c(5, 0, 0, 0)
output_path <- paste0("output/NexusFiles/", "looms_", paste(weight_level, collapse = ""), ".nex")



kd_xl_weighted <- duplication(weights = weight_level,path = path)

# Export
kd_xl_weighted <- kd_xl_weighted |>
  as.matrix() |>
  MatrixToPhyDat()

write.phyDat(kd_xl_weighted, output_path, format = "nexus")
