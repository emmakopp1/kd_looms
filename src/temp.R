library(here)
library(readxl)
library(phangorn)
library(TreeTools)
library(dplyr)
library(ape)
library(tidyverse)
library(phytools)

path <- here("data/Kra-DaiLooms.xlsx")
tt=read_excel(path)

kd_xl <- read_excel(path) %>%
  mutate(across(-Character, as.character)) %>%
  pivot_longer(
    cols = -Character,
    names_to = "trait",
    values_to = "value"
    ) %>%
  mutate(Level = ifelse(Character == "Level", value, NA)) %>%
  fill(Level, .direction = "down") %>%
  filter(Character != "Level")






