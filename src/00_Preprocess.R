library(here)
library(ape)
library(rwty)
library(tidyverse)

dir.create(here("data/looms"))
dir.create(here("data/languages"))

kd_lgs_bcov <- read.nexus(here("data/raw_trees/kd-lgs_bcov.trees"))
kd_lgs_bcov_trimmed <- kd_lgs_bcov[seq(from = 2, to = length(kd_lgs_bcov), by = (length(kd_lgs_bcov) - 1) / 5000)]
write.tree(kd_lgs_bcov_trimmed, here("data/languages/kd_lgs_bcov_trimmed.trees"))

kd_looms_1000 <- read.nexus(here("data/raw_trees/kd-looms_1000.trees"))
kd_looms_1000_trimmed <- kd_looms_1000[seq(from = 2, to = length(kd_looms_1000), by = (length(kd_looms_1000) - 1) / 5000)]
write.tree(kd_looms_1000_trimmed, here("data/looms/kd_looms_1000_trimmed.trees"))

kd_looms_1111 <- read.nexus(here("data/raw_trees/kd-looms_1111.trees"))
kd_looms_1111_trimmed <- kd_looms_1111[seq(from = 2, to = length(kd_looms_1111), by = (length(kd_looms_1111) - 1) / 5000)]
write.tree(kd_looms_1111_trimmed, here("data/looms/kd_looms_1111_trimmed.trees"))

kd_looms_8421 <- read.nexus(here("data/raw_trees/kd-looms_8421.trees"))
kd_looms_8421_trimmed <- kd_looms_8421[seq(from = 2, to = length(kd_looms_8421), by = (length(kd_looms_8421) - 1) / 5000)]
write.tree(kd_looms_8421_trimmed, here("data/looms/kd_looms_8421_trimmed.trees"))

kd_looms_ctmc6 <- read.nexus(here("data/raw_trees/kd-looms_ctmc6.trees"))
kd_looms_ctmc6_trimmed <- kd_looms_ctmc6[seq(from = 2, to = length(kd_looms_ctmc6), by = (length(kd_looms_ctmc6) - 1) / 5000)]
write.tree(kd_looms_ctmc6_trimmed, here("data/looms/kd-looms_ctmc6_trimmed.trees"))



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
