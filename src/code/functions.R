rm(list=ls())
library(here)
source(here("code/init.R"))



read_xls = function(path){
  #'@description
    #'This function readq the data from the Excel datas
  #'@param path Input path
  #'@returns a dataframe

  kd_xls <- read_excel(path, sheet = "Digit", skip = 1) %>%
    select(-starts_with("...")) %>%
    mutate(Character = str_replace_all(Character, " ", "_"))%>%
    column_to_rownames("Character")
  
  return(kd_xls)
  }
  

kd_by_level = function(path){
  #'@description
    #'This function create 4 nexus files where each file contains a level of looms
  #'@param path path of the data
  
  kd_xls = read_xls(path)
  
  for (i in 1:4){
    kd_level <- kd_xls %>%
      t() %>%
      as.data.frame() %>%
      dplyr::filter(Level == i) %>%
      select(-Level) %>%
      t() %>%
      as.data.frame() %>%
      as.matrix() %>%
      MatrixToPhyDat()
    
    # Write the data in a nexus file
    write.phyDat(
      kd_level, 
      sprintf("data/kd_level%i.nex", i), 
      format = "nexus")
  }}



kd_by_level_df = function(path){
  #'@description
    #'This function returns 4 dataframes each one containing the looms of level
  #'@param path path of the data
  
  # Import the data
  kd_xls = read_xls(path)
  
  # Select the data with level i 
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
  
  # Apply the upper function for all level 1 to 4.
  for (i in 1:4){
    var_name = paste0("kd_level", i)
    assign(var_name, kd_to_level(i,kd_xls), envir = .GlobalEnv)
  }
  
  return()
}


duplicate_df = function(data, N) {
  #'@description
  #'This function duplicate a dataframe n times
  return(data[,rep(seq_len(ncol(data)), N) ])
}


duplicate_looms = function(weights){
  #'@description
    #'This function output a new dataframe with a linear combination of each levels
    #'of the datas. For example if weights = c(2,0,1,0) the final dataset will contain
    #'2*kd_level1 + 1* kd_level3.
    #'
  
  file_names <- c("kd_level1", "kd_level2", "kd_level3","kd_level4")
  
  kd_xl_weighted <- do.call(
    cbind, 
    lapply(seq_along(weights), 
           function(i) duplicate_df(get(file_names[i]), weights[i]))
  )
  
  return(kd_xl_weighted)
}

data_to_nexus = function(df,path_out){
  #'@description
    #'This function write the data on a nexus file
  #'@param df 
  #'@param path_out path to stock the nexus file
  
  df <- df |>
    as.matrix() |>
    MatrixToPhyDat()
  
  write.phyDat(df, path_out, format = "nexus")
}



