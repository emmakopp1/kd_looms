rm(list=ls())
library(tidyverse)
library(readxl)
library(phangorn)
library(TreeTools)
library(dplyr)
setwd("/Users/kopp/Library/CloudStorage/OneDrive-UniversitéParis-Dauphine/these/KD_loom/")

kd_xls <- read_excel("data/Kra-DaiLooms28master-2.xlsx", sheet = "Digit", skip = 1) %>%
  select(-starts_with("...")) %>%
  mutate(Character = str_replace_all(Character, " ", "_"))%>%
  column_to_rownames("Character")


kd_by_level = function(kd_xls){
  #'@description
    #'This function create 4 nexus files where each file contains the looms of 
    #'level i.
  
  for (i in 1:4){
    kd_level <- kd_xls %>%
      t() %>%
      as.data.frame() %>%
      filter(Level == i) %>%
      select(-Level) %>%
      t() %>%
      as.data.frame() %>%
      as.matrix() %>%
      MatrixToPhyDat()
    
    write.phyDat(
      kd_level, 
      sprintf("data/data_level/kd_level%i.nex", i), 
      format = "nexus")
  }
  
}

# Create nexus files
kd_by_level(kd_xls)


# Create data frames
kd_by_level_df = function(kd_xls){
  #'@description
    #'This function returns 4 dataframes each one containing the looms of level
    #'i for i in 1:4
  
  kd_to_level = function(i,kd_xls){
    kd_level <- kd_xls %>%
      t() %>%
      as.data.frame() %>%
      filter(Level == i) %>%
      select(-Level) %>%
      t() %>%
      as.data.frame() 
    return(kd_level)
  }
  
  for (i in 1:4){
    var_name = paste0("kd_level", i)
    assign(var_name, kd_to_level(i,kd_xls), envir = .GlobalEnv)
  }
  return(list(
    kd_level1=kd_level1,
    kd_level2=kd_level2,
    kd_level3=kd_level3,
    kd_level4=kd_level4
    ))
}

res = kd_by_level_df(kd_xls)

colnames(kd_level1)
