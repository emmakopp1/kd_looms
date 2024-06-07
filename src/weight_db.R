# Import file
library(here)
source(here("src/functions.R"))
library(stringr)
library(tidyverse)
library(readxl)
library(phangorn)
library(TreeTools)
library(dplyr)
library(ape)
library(phytools)


# Init
path <- here("data/Kra-DaiLooms28master-2.xlsx")
weight_level <- c(30, 0, 0, 0)
output_path <- paste0("output/NexusFiles/", "looms_", paste(weight_level, collapse = ""), ".nex")

kd_xl_weighted <- duplication_with_option(
  weights = weight_level,
  option = "level",
  path = path
)

data_to_nexus(kd_xl_weighted, output_path)
