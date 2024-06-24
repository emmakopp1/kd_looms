library(here)
library(ape)

library(tidyverse)
read_csv(here("data/kd_looms_locs.csv")) |> 
  rename(language = loom) |> 
  mutate(lng_group_code = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |> 
  mutate(language = str_remove(language, "^[A-Z][a-z]+(?=[A-Z])")) |> 
  left_join(read_csv(here("data/kd_looms_languages.csv")) |> 
              mutate(lng_group_code = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |> 
              mutate(language = str_remove(language, "^[A-Z][a-z]+(?=[A-Z])"))
              ) |> 
  mutate(language = ifelse(language == "Liangzhu", NA, language)) |> 
  left_join(tribble(
    ~lng_group_code, ~lng_group,
    "Be", "Ong-Be",
    "Kra", "Kra",
    "Hlai", "Hlai",
    "Ks", "Kam-Sui",
    "Tn", "Northern Tai",
    "Tsw", "Southwestern Tai",
    "Tc", "Central Tai"
  )) |> 
  rename(lng = language, loom_type_code = loom_code) |> 
  select(group, loom_type_code, loom_type, lng, lng_group_code, lng_group, lon, lat) |> 
  write_csv(here("data/kd_looms.csv"))
  
read_csv(here("data/kd_lgs_locs.csv")) |> 
  mutate(lng_group_code = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |> 
  mutate(language = str_remove(language, "^[A-Z][a-z]+(?=[A-Z])")) |> 
  rename(lng = language) |> 
  left_join(tribble(
    ~lng_group_code, ~lng_group,
    "Be", "Ong-Be",
    "Kra", "Kra",
    "Hlai", "Hlai",
    "Ks", "Kam-Sui",
    "Tn", "Northern Tai",
    "Tsw", "Southwestern Tai",
    "Tc", "Central Tai"
  )) |> 
  select(lng_group_code, lng_group, lng, lon, lat) |> 
  write_csv(here("data/kd_lngs.csv"))

