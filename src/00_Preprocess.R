library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tidyverse)


# Prepare the nexus files -------------------------------------------------

# Function to force nexus files to be in binary format
write_binary_nexus <- function(x, file) {
  write.phyDat(x, file, format = "nexus")
  read_lines(file) |>
    str_replace('symbols="0123456789"', 'symbols="01"') |>
    str_subset("^(?!\\[Data written by)") |>
    write_lines(file)
}

# Languages
kd_lgs_lx <- read_csv(here("data/kd-lgs/kd-lgs_lx.csv")) |>
  select(-concept) |>
  pivot_longer(-c(concept_id, cogid), names_to = "Taxon") |>
  filter(!is.na(value)) |>
  distinct(concept_id, cogid, Taxon) |>
  rowid_to_column() |>
  pivot_wider(
    names_from = Taxon,
    values_from = rowid,
    values_fill = 0,
    values_fn = length
  ) |>
  pivot_longer(cols = !(concept_id | cogid), names_to = "Taxon", values_to = "value") |>
  mutate(value = as.character(value)) |>
  group_by(concept_id, Taxon) |>
  mutate(allzero = sum(value != "0") == 0) |>
  rowwise() |>
  mutate(value = ifelse(allzero, "?", value)) |>
  ungroup() |>
  mutate(id = paste0(concept_id, "_", cogid)) |>
  select(Taxon, id, value) |>
  pivot_wider(names_from = id, values_from = value) |>
  arrange(Taxon)
kd_lgs_matrix <- kd_lgs_lx |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write_binary_nexus(kd_lgs_matrix, here("data/kd-lgs/kd-lgs_bcov/kd-lgs.nex"))

# Looms
kd_looms_characters <- read_csv(here("data/kd-looms/kd-looms_characters.csv")) |>
  select(code, level)
kd_looms_matrix <- read_csv(here("data/kd-looms/kd-looms_matrix.csv")) |>
  mutate(Taxon = str_replace_all(Taxon, " ", "_"))

# Looms, all levels, no weighting
kd_looms_matrix1111 <- kd_looms_matrix |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write_binary_nexus(kd_looms_matrix1111, here("data/kd-looms/kd-looms_bcov1111/kd-looms_1111.nex"))

# Looms, level 1 characters only
kd_looms_matrix1000 <- select(kd_looms_matrix, Taxon, filter(kd_looms_characters, level == 1)$code) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write_binary_nexus(kd_looms_matrix1000, here("data/kd-looms/kd-looms_bcov1000/kd-looms_1000.nex"))

# Looms, weighted characters
weights <- c(2, 1, 1, 1)

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

# Looms, 4 variable rates
kd_looms_matrix_bylevel <- kd_looms_matrix |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    cols = -Taxon,
    names_to = "code",
    values_to = "value"
  ) |>
  left_join(kd_looms_characters) |>
  arrange(level) |>
  select(-level) |>
  pivot_wider(
    names_from = c(code),
    values_from = value
  ) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

kd_looms_partition <- kd_looms_characters |>
  arrange(level) |>
  rowid_to_column() |>
  group_by(level) |>
  summarise(charset = paste0("    charset level", unique(level), " = ", min(rowid), "-", max(rowid), ";")) |>
  pull(charset) |>
  paste0(collapse = "\n")


write_binary_nexus(kd_looms_matrix_bylevel, here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.nex"))
write_file(str_glue("begin assumptions;\n{kd_looms_partition}\nend;\n"), here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.nex"), append = TRUE)


# Looms, 6 variable rates, with instructions for MrBayes
write_binary_nexus(kd_looms_matrix1111, here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6.nex"))
txt_mrbayes <- "begin mrbayes;
	set autoclose=yes nowarn=yes;
	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);
	prset applyto=(all) ratepr=variable;
	lset rates=gamma ngammacat=6;
	mcmcp ngen= 20000000 relburnin=yes burninfrac=0.25  printfreq=10000  samplefreq=1000 nchains=4 savebrlens=yes;
	mcmc;
	sumt relburnin=yes burninfrac=0.25;
	sump;
end;

"
write_file(txt_mrbayes, here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6.nex"), append = TRUE)


# Trim raw tree files -----------------------------------------------------

k <- 5000

kd_lgs_bcov <- read.nexus(here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.trees"))
kd_lgs_bcov_trimmed <- kd_lgs_bcov[seq(from = 2, to = length(kd_lgs_bcov), by = (length(kd_lgs_bcov) - 1) / k)]
write.tree(kd_lgs_bcov_trimmed, here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov_trimmed.trees"))

kd_looms_bcov1000 <- read.nexus(here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000.trees"))
kd_looms_bcov1000_trimmed <- kd_looms_bcov1000[seq(from = 2, to = length(kd_looms_bcov1000), by = (length(kd_looms_bcov1000) - 1) / k)]
write.tree(kd_looms_bcov1000_trimmed, here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000_trimmed.trees"))

kd_looms_bcov1111 <- read.nexus(here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees"))
kd_looms_bcov1111_trimmed <- kd_looms_bcov1111[seq(from = 2, to = length(kd_looms_bcov1111), by = (length(kd_looms_bcov1111) - 1) / k)]
write.tree(kd_looms_bcov1111_trimmed, here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111_trimmed.trees"))

kd_looms_bcov8421 <- read.nexus(here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421.trees"))
kd_looms_bcov8421_trimmed <- kd_looms_bcov8421[seq(from = 2, to = length(kd_looms_bcov8421), by = (length(kd_looms_bcov8421) - 1) / k)]
write.tree(kd_looms_bcov8421_trimmed, here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421_trimmed.trees"))

kd_looms_ctmc4 <- read.nexus(here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.trees"))
kd_looms_ctmc4_trimmed <- kd_looms_ctmc4[seq(from = 2, to = length(kd_looms_ctmc4), by = (length(kd_looms_ctmc4) - 1) / k)]
write.tree(kd_looms_ctmc4_trimmed, here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4_trimmed.trees"))

kd_looms_ctmc6 <- read.nexus(here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6.trees"))
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
