rm(list=ls())

# Import file 
library(here)
source(here("code/functions.R"))
source(here("code/init.R"))


# Path to data
path = here("data/Kra-DaiLooms28master-2.xlsx")

# Create nexus file of each levels
kd_by_level(path)

# Create dataframes of each levels
kd_by_level_df(path)
kd_xl_weighted = duplicate_looms(c(4,1,1,1))
data_to_nexus(kd_xl_weighted,here("outputs/test.nex"))
