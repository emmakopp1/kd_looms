library(here)
source(here("src/init.R"))

duplication_with_option = function(weights,option,path){
  # Read data
  kd_xl = read_xls(path)
  
  # Create data frames 
  if (option=='level'){data = create_kd_by_level_df(path)}
  
  if (option =='type'){
    # Group different types of looms
    simple <- colnames(kd_xl) |>
      str_subset("^FM|HP|RP|RT|RW|MH")
    
    complex <- colnames(kd_xl) |>
      str_subset("^PS")
    
    basics <- setdiff(colnames(kd_xl),union(simple,complex))
    data = list(simple=kd_xl[,simple],basic=kd_xl[,basics],complex=kd_xl[,complex])
  }
  
  # Duplication
  kd_xl_weighted = duplicate_looms(
    weights=weights,
    option=option,
    data=data)
  
  return(kd_xl_weighted)
}

read_xls = function(path){
  #'@description
    #'This function readq the data from the Excel datas.
  #'@param path Input path
  #'@returns a dataframe

  kd_xls <- read_excel(path, sheet = "Digit", skip = 1) %>%
    select(-starts_with("...")) %>%
    mutate(Character = str_replace_all(Character, " ", "_"))%>%
    column_to_rownames("Character")
  
  return(kd_xls)
  }

create_kd_by_level_df = function(path){
  #'@description
    #'This function returns 4 dataframes each one containing the looms of level
  #'@param path path of the data
  
  # Import the data
  kd_xls = read_xls(path)
  
  # Instanciate a dataframe for a given level i
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
  
  # Create four data-frames one for each level
  for (i in 1:4){
    var_name = paste0("kd_level", i)
    assign(var_name, kd_to_level(i,kd_xls), envir = .GlobalEnv)
  }
  data = list(kd_level1 =kd_level1,
              kd_level2 =kd_level2,
              kd_level3 =kd_level3,
              kd_level4 =kd_level4)
  
  return(data)
  }

duplicate_looms = function(weights,option,data){

  if (option == "level"){
    file_names <- c("kd_level1", "kd_level2", "kd_level3","kd_level4")
  }
  
  if (option == "type"){
    file_names = c("simple","basic","complex")
  }
  
  kd_xl_weighted <- do.call(
    cbind, 
    lapply(seq_along(weights), 
           function(i) duplicate_df(data[[file_names[i]]], weights[i]))
  )
  
  return(kd_xl_weighted)
}

duplicate_df = function(data, N) {
  #'@description
  #'This function duplicate a dataframe n times
  return(data[,rep(seq_len(ncol(data)), N) ])
}

data_to_nexus = function(df,path_out,option){
  #'@description
    #'This function write the data on a nexus file
  #'@param df dataframe
  #'@param path_out path to stock the nexus file

  if (option=="type"){
    df = df %>%
      slice(-which(rownames(df) == "Level"))
  }
  df <- df |>
    as.matrix() |>
    MatrixToPhyDat()
  
  
    
  write.phyDat(df, path_out, format = "nexus")
}

kd_by_level_nex = function(path){
  #'@description
  #'This function create 4 nexus files where each file contains a level of looms
  #'@param path path of the data
  
  kd_xls = read_xls(path)
  
  for (i in 1:4){
    # Pre-process data
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
  }
  return()}
