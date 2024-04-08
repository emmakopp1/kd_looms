rm(list=ls())
library(tidyverse)
library(readxl)
library(phangorn)
library(TreeTools)
library(ape)
library(phytools)
library(here)
setwd("/Users/kopp/Library/CloudStorage/OneDrive-UniversitéParis-Dauphine/these/kd_loom/")

# Functions
duplicate_df = function(data, N) {
  #'@description
    #'This function duplicate a dataframe n times

  return(data[,rep(seq_len(ncol(data)), N) ])
}


duplicate_looms = function(weights){
  
  file_names <- c("kd_level1", "kd_level2", "kd_level3","kd_level4")
  
  # Separate initial dataset into basic, simlpe and complex datasets
  kd_xl_weighted <- do.call(
    cbind, 
    lapply(seq_along(weights), 
           function(i) duplicate_df(get(file_names[i]), weights[i]))
  )
  
  return(kd_xl_weighted)
}

# Weight
weights = c(4,1,0,0)
kd_xl_weighted = duplicate_looms(weights)

kd_xl_weighted <- kd_xl_weighted |>
  as.matrix() |>
  MatrixToPhyDat()

write.phyDat(kd_xl_weighted, "loom4100/kd_looms_4100.nex", format = "nexus")
