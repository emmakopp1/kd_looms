library(here)

duplication_with_option = function(weights,option,path){
  # Read data
  kd_xl = read_xls(path)
  rownames = row.names(kd_xl)
  kd_xl = as.data.frame(sapply(kd_xl,as.numeric))
  row.names(kd_xl) = rownames
  
  # Create data frames 
  if (option=='level'){data = create_kd_by_level_df(path)}
  
  if (option =='type'){
    # Group different types of looms
    simple <- colnames(kd_xl) |>
      str_subset("^FM|HP|RP|RT|RW|MH")
    
    complex <- colnames(kd_xl) |>
      str_subset("^PS")
    
    basics <- setdiff(colnames(kd_xl),union(simple,complex))
    data <- list(simple=kd_xl[,simple],basic=kd_xl[,basics],complex=kd_xl[,complex])
  }
  
  # Duplication
  kd_xl_weighted = duplicate_looms(
    weights=weights,
    option=option,
    data=data)
  
  return(kd_xl_weighted)
}

read_xls = function(path){

  kd_xls <- read_excel(path, sheet = "Digit", skip = 1) %>%
    select(-starts_with("...")) %>%
    mutate(Character = str_replace_all(Character, " ", "_"))%>%
    column_to_rownames("Character")
  
  return(kd_xls)
  }

create_kd_by_level_df = function(path){

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
    file_names = c("basic","simple","complex")
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


