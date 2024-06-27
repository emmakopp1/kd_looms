library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tidyverse)


# Prepare the nexus files -------------------------------------------------

weights <- c(2, 1, 1, 1)

kd_looms_characters <- read_csv(here("data/kd-looms/kd-looms_characters.csv")) |>
  select(code, level)

kd_looms_matrix <- read_csv(here("data/kd-looms/kd-looms_matrix.csv"))

kd_looms_matrix1111 <-  kd_looms_matrix |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write.phyDat(kd_looms_matrix1111, here("data/kd-looms/kd-looms_bcov1111/kd-looms_1111.nex"), format = "nexus")

kd_looms_matrix1000 <- select(kd_looms_matrix, Taxon, filter(kd_looms_characters, level == 1)$code) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write.phyDat(kd_looms_matrix1000, here("data/kd-looms/kd-looms_bcov1000/kd-looms_1000.nex"), format = "nexus")

kd_looms_matrix_weighted <- kd_looms_matrix |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    cols = -Taxon,
    names_to = "code",
    values_to = "value"
  ) |>
  left_join(kd_looms_characters) |>
  group_by(level) |>
  group_map(~ uncount(.x, weights[.y$level]), .keep = TRUE) |>
  bind_rows() |>
  group_by(Taxon, code, level) |>
  mutate(set = row_number(), code = paste0(code, "_", set)) |>
  ungroup() |>
  select(-c(level, set)) |>
  pivot_wider(
    names_from = c(code),
    values_from = value
  )

kd_looms_matrix_weighted |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()


# Trim raw tree files -----------------------------------------------------

k <- 5000

kd_lgs_bcov <- read.nexus(here("data/raw_trees/kd-lgs_bcov.trees"))
kd_lgs_bcov_trimmed <- kd_lgs_bcov[seq(from = 2, to = length(kd_lgs_bcov), by = (length(kd_lgs_bcov) - 1) / k)]
write.tree(kd_lgs_bcov_trimmed, here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov_trimmed.trees"))

kd_looms_bcov1000 <- read.nexus(here("data/raw_trees/kd-looms_bcov1000.trees"))
kd_looms_bcov1000_trimmed <- kd_looms_bcov1000[seq(from = 2, to = length(kd_looms_bcov1000), by = (length(kd_looms_bcov1000) - 1) / k)]
write.tree(kd_looms_bcov1000_trimmed, here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000_trimmed.trees"))

kd_looms_bcov1111 <- read.nexus(here("data/raw_trees/kd-looms_bcov1111.trees"))
kd_looms_bcov1111_trimmed <- kd_looms_bcov1111[seq(from = 2, to = length(kd_looms_bcov1111), by = (length(kd_looms_bcov1111) - 1) / k)]
write.tree(kd_looms_bcov1111_trimmed, here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111_trimmed.trees"))

kd_looms_bcov8421 <- read.nexus(here("data/raw_trees/kd-looms_bcov8421.trees"))
kd_looms_bcov8421_trimmed <- kd_looms_bcov8421[seq(from = 2, to = length(kd_looms_bcov8421), by = (length(kd_looms_bcov8421) - 1) / k)]
write.tree(kd_looms_bcov8421_trimmed, here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421_trimmed.trees"))

kd_looms_ctmc6 <- read.nexus(here("data/raw_trees/kd-looms_ctmc6.trees"))
kd_looms_ctmc6_trimmed <- kd_looms_ctmc6[seq(from = 2, to = length(kd_looms_ctmc6), by = (length(kd_looms_ctmc6) - 1) / k)]
write.tree(kd_looms_ctmc6_trimmed, here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6_trimmed.trees"))



# read_csv(here("data/kd_looms_locs.csv")) |>
#   rename(language = loom) |>
#   mutate(lng_group_code = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |>
#   mutate(language = str_remove(language, "^[A-Z][a-z]+(?=[A-Z])")) |>
#   left_join(read_csv(here("data/kd_looms_languages.csv")) |>
#               mutate(lng_group_code = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |>
#               mutate(language = str_remove(language, "^[A-Z][a-z]+(?=[A-Z])"))
#               ) |>
#   mutate(language = ifelse(language == "Liangzhu", NA, language)) |>
#   left_join(tribble(
#     ~lng_group_code, ~lng_group,
#     "Be", "Ong-Be",
#     "Kra", "Kra",
#     "Hlai", "Hlai",
#     "Ks", "Kam-Sui",
#     "Tn", "Northern Tai",
#     "Tsw", "Southwestern Tai",
#     "Tc", "Central Tai"
#   )) |>
#   rename(lng = language, loom_type_code = loom_code) |>
#   select(group, loom_type_code, loom_type, lng, lng_group_code, lng_group, lon, lat) |>
#   write_csv(here("data/kd_looms.csv"))
#
# read_csv(here("data/kd_lgs_locs.csv")) |>
#   mutate(lng_group_code = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |>
#   mutate(language = str_remove(language, "^[A-Z][a-z]+(?=[A-Z])")) |>
#   rename(lng = language) |>
#   left_join(tribble(
#     ~lng_group_code, ~lng_group,
#     "Be", "Ong-Be",
#     "Kra", "Kra",
#     "Hlai", "Hlai",
#     "Ks", "Kam-Sui",
#     "Tn", "Northern Tai",
#     "Tsw", "Southwestern Tai",
#     "Tc", "Central Tai"
#   )) |>
#   select(lng_group_code, lng_group, lng, lon, lat) |>
#   write_csv(here("data/kd_lngs.csv"))
