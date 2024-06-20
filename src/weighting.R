library(here)
library(phytools)
library(tidyverse)
library(TreeTools)

weights <- c(2,1,1,1)

looms_m_levels <- read_csv(here("data/kd_looms_matrix_levels.csv")) |>
  mutate(across(everything(), as.character))

chr_levels <- looms_m_levels |>
  filter(Taxon == "Level") |>
  pivot_longer(
    cols = -Taxon,
    names_to = "Character",
    values_to = "Level"
  ) |>
  select(-Taxon)

write_csv(chr_levels, here("output/kd_looms_chr_levels.csv"))

looms_m_weighted <- looms_m_levels |>
  filter(Taxon != "Level") |>
  pivot_longer(
    cols = -Taxon,
    names_to = "Character",
    values_to = "Value"
  ) |>
  left_join(chr_levels) |>
  group_by(Level) |>
  group_map(~ uncount(.x, weights[as.integer(.y$Level)]), .keep = TRUE) |>
  bind_rows() |>
  group_by(Taxon, Character, Level) |>
  mutate(set = row_number(), Character = paste0(Character, "_", set)) |>
  ungroup() |>
  select(-c(Level, set)) |>
  pivot_wider(
    names_from = c(Character),
    values_from = Value
  )

looms_phyDat_weighted <- looms_m_weighted |>
  column_to_rownames("Taxon") |>
  as.matrix() |> 
  MatrixToPhyDat()

# # Init
# path <- here("data/Kra-DaiLooms.xlsx")
# weights <- c(2,1,1,1)
# output_path <- paste0("output/NexusFiles/", "looms_", paste(weights, collapse = ""), ".nex")
#
# # Function to replicate
# replicate_rows <- function(df, times) {
#   df[rep(seq_len(nrow(df)), times), ]
# }
#
#
# kd_xl <- read_excel(path) %>%
#   mutate(across(-Character, as.character)) %>%
#   slice(1)%>%
#   pivot_longer(
#     cols = -Character,
#     names_to = "Trait",
#     values_to = "Level"
#   ) %>%
#   select(-Character) %>%
#   right_join(
#     read_excel(path) %>% mutate(across(-Character, as.character))%>% slice(-1) %>% pivot_longer(
#       cols = -Character,
#       names_to = "Trait",
#       values_to = "Value"
#     ),
#     by = "Trait"
#   ) %>%
#   rename(Language = Character) %>%
#   select(Language, Trait, Level, Value) %>%
#   group_split(Level) %>%
#   map_dfr(~ {
#     if (.x$Level[1] == 1) {
#       replicate_rows(.x, weights[1])
#     } else if (.x$Level[1] == 2) {
#       replicate_rows(.x, weights[2])
#     } else if (.x$Level[1] == 3) {
#       replicate_rows(.x, weights[3])
#     } else if (.x$Level[1] == 4) {
#       replicate_rows(.x, weights[4])
#     } else {
#       .x
#     }})%>%
#   select(Language, Trait, Value)%>%
#   arrange(Language,Trait) %>%
#   pivot_wider(
#     names_from = Trait,
#     values_from = Value,
#     values_fn = list
#   )


# Export
# kd_xl_weighted <- kd_xl_weighted |>
#  as.matrix() |>
#  MatrixToPhyDat()

# write.phyDat(kd_xl, output_path, format = "nexus")
