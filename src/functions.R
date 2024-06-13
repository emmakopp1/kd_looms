library(here)

duplication <- function(weights, path) {
  # Read data
  kd_xl <- read_excel(path, sheet = "Digit", skip = 1) %>%
    select(-starts_with("...")) %>%
    mutate(Character = str_replace_all(Character, " ", "_")) %>%
    column_to_rownames("Character") %>%
    sapply(as.numeric) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(var = "Character") %>%
    column_to_rownames("Character")

  # Create data frames
  data <- create_kd_by_level_df(path)

  # Duplication
  kd_xl_weighted <- duplicate_looms(weights, data)

  return(kd_xl_weighted)
}

create_kd_by_level_df <- function(path) {

  # Instanciate a dataframe for a given level i
  kd_to_level <- function(i, kd_xls) {
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
  for (i in 1:4) {
    var_name <- paste0("kd_level", i)
    assign(var_name, kd_to_level(i, kd_xls), envir = .GlobalEnv)
  }
  data <- list(
    kd_level1 = kd_level1,
    kd_level2 = kd_level2,
    kd_level3 = kd_level3,
    kd_level4 = kd_level4
  )

  return(data)
}


duplicate_looms <- function(weights, data) {
  file_names <- c("kd_level1", "kd_level2", "kd_level3", "kd_level4")
  
  seq_along(weights) %>%
    lapply(function(i) {
      data[[file_names[i]]] %>% duplicate_df(weights[i])
    }) %>%
    do.call(cbind, .)
}

duplicate_df <- function(data, N) {
  return(data[, rep(seq_len(ncol(data)), N)])
}

read_xls <- function(path) {
  kd_xls <- read_excel(path, sheet = "Digit", skip = 1) %>%
    select(-starts_with("...")) %>%
    mutate(Character = str_replace_all(Character, " ", "_")) %>%
    column_to_rownames("Character")

  return(kd_xls)
}
