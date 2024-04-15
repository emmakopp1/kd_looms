# Import file 
library(here)
source(here("src/functions.R"))

# Init
path = here("data/Kra-DaiLooms28master-2.xlsx")
weight_type = c(0,1,1)
#weight_level = c(1,0,0,0) 

kd_xl_weighted = duplication_with_option(
  weights = weight_type,
  option = "type",
  path = path)

# Save data
data_to_nexus(kd_xl_weighted,here("output/by_type/loom011/kd_loom011.nex"),option="type")
# Create dataframes of each levels
#kd_by_level_df(path)
#kd_xl_weighted = duplicate_looms(c(4,1,1,1))

