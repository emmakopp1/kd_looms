# Import file 
library(here)
source(here("src/functions.R"))
library(stringr)

# Init
path = here("data/Kra-DaiLooms28master-2.xlsx")
weight_level = c(20,0,0,0) 
output_path = paste0(
  "output/by_level/loom",
  paste(weight_level,collapse =''),
  "/kd_loom",
  paste(weight_level,collapse =''),
  ".nex")

kd_xl_weighted = duplication_with_option(
  weights = weight_level,
  option = "level",
  path = path)

# Save data
dir.create(paste0(
  "output/by_level/loom",
  paste(weight_level,collapse ='')))

data_to_nexus(
  kd_xl_weighted,here(output_path),
  option="level")



