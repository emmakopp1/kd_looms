library(here)
library(readxl)
library(phangorn)
library(TreeTools)
library(dplyr)
library(ape)
library(tidyverse)
library(phytools)


# Init
path <- here("data/Kra-DaiLooms.xlsx")
weights <- c(2,1,1,1)
output_path <- paste0("output/NexusFiles/", "looms_", paste(weights, collapse = ""), ".nex")

# Function to replicate
replicate_rows <- function(df, times) {
  df[rep(seq_len(nrow(df)), times), ]
}


kd_xl <- read_excel(path) %>%
  mutate(across(-Character, as.character)) %>%
  slice(1)%>%
  pivot_longer(
    cols = -Character,
    names_to = "Trait",
    values_to = "Level"
  ) %>%
  select(-Character) %>%
  right_join(
    read_excel(path) %>% mutate(across(-Character, as.character))%>% slice(-1) %>% pivot_longer(
      cols = -Character,
      names_to = "Trait",
      values_to = "Value"
    ),
    by = "Trait"
  ) %>%
  rename(Language = Character) %>%
  select(Language, Trait, Level, Value) %>%
  group_split(Level) %>%
  map_dfr(~ {
    if (.x$Level[1] == 1) {
      replicate_rows(.x, weights[1])
    } else if (.x$Level[1] == 2) {
      replicate_rows(.x, weights[2])
    } else if (.x$Level[1] == 3) {
      replicate_rows(.x, weights[3])
    } else if (.x$Level[1] == 4) {
      replicate_rows(.x, weights[4])
    } else {
      .x
    }})%>%
  select(Language, Trait, Value)%>%
  arrange(Language,Trait) %>%
  pivot_wider(
    names_from = Trait,
    values_from = Value,
    values_fn = list
  )
  

# Export
#kd_xl_weighted <- kd_xl_weighted |>
#  as.matrix() |>
#  MatrixToPhyDat()

#write.phyDat(kd_xl, output_path, format = "nexus")


