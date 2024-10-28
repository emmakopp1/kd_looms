library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tracerer)
library(tidyverse)
library(fs)

dir.create(here("data/kd-lgs/"))
dir.create(here("data/kd-looms/"))


# Prepare the nexus files ------------------------------------------------------

# Function to force NEXUS files to specify the symbols to be 0 and 1
write_binary_nexus <- function(x, file) {
  write.phyDat(x, file, format = "nexus")
  read_lines(file) |>
    str_replace('symbols="0123456789"', 'symbols="01"') |>
    str_subset("^(?!\\[Data written by)") |>
    write_lines(file)
}

# Languages
kd_lgs_lx <- read_csv(here("data/kd-lgs/kd-lgs_lx.csv")) |>
  pivot_longer(-c(concept, concept_id, cogid), names_to = "Taxon") |>
  filter(!is.na(value)) |>
  distinct(concept, concept_id, cogid, Taxon) |>
  rowid_to_column() |>
  pivot_wider(
    names_from = Taxon,
    values_from = rowid,
    values_fill = 0,
    values_fn = length
  ) |>
  pivot_longer(
    cols = !(concept | concept_id | cogid),
    names_to = "Taxon",
    values_to = "value"
  ) |>
  mutate(value = as.character(value)) |>
  group_by(concept, concept_id, Taxon) |>
  mutate(allzero = sum(value != "0") == 0) |>
  rowwise() |>
  mutate(value = ifelse(allzero, "?", value)) |>
  ungroup() |>
  mutate(id = paste0(concept_id, "_", cogid)) |>
  select(Taxon, concept, concept_id, id, value) |>
  arrange(concept_id)

kd_lgs_matrix <- kd_lgs_lx |>
  select(Taxon, id, value) |>
  pivot_wider(names_from = id, values_from = value) |>
  arrange(Taxon) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_lgs_matrix,
  here("data/nexus/kd-lgs.nex")
)

# Looms
kd_looms_characters <- read_csv(here("data/kd-looms/kd-looms_characters.csv")) |>
  select(code, level, type)
kd_looms_matrix <- read_csv(here("data/kd-looms/kd-looms_matrix.csv")) |>
  mutate(Taxon = str_replace_all(Taxon, " ", "_"))

## Looms, all levels, no weighting
kd_looms_matrix1111 <- kd_looms_matrix |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix1111,
  here("data/nexus/kd-looms_1111.nex")
)

## Looms, level 1 characters only
kd_looms_matrix1000 <- select(
  kd_looms_matrix,
  Taxon,
  filter(kd_looms_characters, level == 1)$code
) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix1000,
  here("data/nexus/kd-looms_1000.nex")
)

## Looms, weighted characters
weights <- c(8, 4, 2, 1)

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
  select(-c(level, set, type)) |>
  pivot_wider(
    names_from = code,
    values_from = value
  )

kd_looms_matrix8421 <- kd_looms_matrix_weighted |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix8421,
  here("data/nexus/kd-looms_8421.nex")
)

## Looms, 4 variable rates
# dir.create(here("data/kd-looms/kd-looms_ctmc4/"))

kd_looms_matrix_bylevel <- kd_looms_matrix |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    cols = -Taxon,
    names_to = "code",
    values_to = "value"
  ) |>
  left_join(kd_looms_characters) |>
  arrange(level) |>
  select(-c(level, type)) |>
  pivot_wider(
    names_from = code,
    values_from = value
  ) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

kd_looms_partition <- kd_looms_characters |>
  arrange(level) |>
  rowid_to_column() |>
  group_by(level) |>
  summarise(charset = paste0(
    "    charset level",
    unique(level),
    " = ",
    min(rowid),
    "-",
    max(rowid), ";"
  )) |>
  pull(charset) |>
  paste0(collapse = "\n")

write_binary_nexus(
  kd_looms_matrix_bylevel,
  here("data/nexus/kd-looms_ctmc4.nex")
)

write_file(str_glue("begin assumptions;\n{kd_looms_partition}\nend;\n"),
  here("data/nexus/kd-looms_ctmc4.nex"),
  append = TRUE
)

## Looms, basic features only

kd_looms_matrix_basic <- select(
  kd_looms_matrix,
  Taxon,
  filter(kd_looms_characters, type == "Loom basics")$code
) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix_basic,
  here("data/nexus/kd-looms_basic.nex")
)

## Looms, pattern features only

kd_looms_matrix_patterns <- select(
  kd_looms_matrix,
  Taxon,
  filter(kd_looms_characters, type != "Loom basics")$code
) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix_patterns,
  here("data/nexus/kd-looms_patterns.nex")
)

# Add the missing "End;" line at the end of the tree files ---------------

tree_files <- dir_ls(here("data"), glob = "*.trees", recurse = TRUE)

# Fonction pour ajouter "End;" si elle n'est pas présente à la fin
add_end_if_missing <- function(file_path) {
  lines <- read_lines(file_path)

  if (last(lines) != "End;") {
    write_lines("End;", file_path, append = TRUE)
    message(glue::glue("Ajout de 'End;' à la fin de {file_path}"))
  } else {
    message(glue::glue("Le fichier {file_path} contient déjà 'End;' à la fin"))
  }
}

# Appliquer la fonction sur chaque fichier
walk(tree_files, add_end_if_missing)


# Check the ESS values ---------------------------------------------------------

burnin <- .1

# Languages
# Binary covarion relaxed clock
kd_lgs_bcov_relaxed_log <- parse_beast_tracelog_file(
  here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-lgs_bcov")
kd_lgs_bcov_relaxed_ess <- kd_lgs_bcov_relaxed_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_lgs_bcov_relaxed_log$Sample) / (nrow(kd_lgs_bcov_relaxed_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_lgs_bcov_relaxed_ess, ESS < 200)

# Binary covarion strict clock
kd_lgs_bcov_strict_log <- parse_beast_tracelog_file(
  here("data/kd-lgs/kd-lgs_bcov_strict/kd-lgs_bcov_strict.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-lgs_bcov")
kd_lgs_bcov_strict_ess <- kd_lgs_bcov_strict_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_lgs_bcov_strict_log$Sample) / (nrow(kd_lgs_bcov_strict_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_lgs_bcov_strict_ess, ESS < 200)

# Binary covarion relaxed clock by concept
kd_lgs_bcov_relaxed_byconcept_log <- parse_beast_tracelog_file(
  here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/kd-lgs_bcov_relaxed_byconcept.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-lgs_bcov")
kd_lgs_bcov_relaxed_byconcept_ess <- kd_lgs_bcov_relaxed_byconcept_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_lgs_bcov_relaxed_byconcept_log$Sample) / (nrow(kd_lgs_bcov_relaxed_byconcept_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_lgs_bcov_relaxed_byconcept_ess, ESS < 200)

# Looms

## Looms 1000
kd_looms_bcov1000_strict_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1000_strict/kd-looms_bcov1000_strict.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov1000")
kd_looms_bcov1000_strict_ess <- kd_looms_bcov1000_strict_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov1000_strict_log$Sample) / (nrow(kd_looms_bcov1000_strict_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov1000_strict_ess, ESS < 200)

## Looms 1111 bcov strict
kd_looms_bcov1111_strict_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1111_strict/kd-looms_bcov1111_strict.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov1111")
kd_looms_bcov1111_strict_ess <- kd_looms_bcov1111_strict_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov1111_strict_log$Sample) / (nrow(kd_looms_bcov1111_strict_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov1111_strict_ess, ESS < 200)

## Looms 1111 bcov relaxed
kd_looms_bcov1111_relaxed_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1111_relaxed/kd-looms_bcov1111_relaxed.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov1111")
kd_looms_bcov1111_relaxed_ess <- kd_looms_bcov1111_relaxed_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov1111_relaxed_log$Sample) / (nrow(kd_looms_bcov1111_relaxed_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov1111_relaxed_ess, ESS < 200)

## Looms 8421
kd_looms_bcov8421_strict_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov8421_strict/kd-looms_bcov8421_strict.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov8421")
kd_looms_bcov84211_strict_ess <- kd_looms_bcov8421_strict_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov8421_strict_log$Sample) / (nrow(kd_looms_bcov8421_strict_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov84211_strict_ess, ESS < 200)

## Looms CTMC strict clock heterogene
kd_looms_ctmc1111_strict_heterogene_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_ctmc1111_strict_heterogene/kd-looms_ctmc1111_strict_heterogene.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd_looms_ctmc4")
kd_looms_ctmc1111_strict_heterogene_ess <- kd_looms_ctmc1111_strict_heterogene_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_ctmc1111_strict_heterogene_log$Sample) / (nrow(kd_looms_ctmc1111_strict_heterogene_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_ctmc1111_strict_heterogene_ess, ESS < 200)

## Looms CTMC strict clock
kd_looms_ctmc1111_strict_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_ctmc1111_strict/kd-looms_ctmc1111_strict.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd_looms_ctmc4")
kd_looms_ctmc1111_strict_ess <- kd_looms_ctmc1111_strict_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_ctmc1111_strict_log$Sample) / (nrow(kd_looms_ctmc1111_strict_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_ctmc1111_strict_ess, ESS < 200)

## Looms, basic features only
kd_looms_bcov_basic_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov_basic/kd-looms_bcov_basic.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov_basic")
kd_looms_bcov_basic_ess <- kd_looms_bcov_basic_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov_basic_log$Sample) / (nrow(kd_looms_bcov_basic_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov_basic_ess, ESS < 200)

## Looms, pattern features only
kd_looms_bcov_patterns_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov_patterns/kd-looms_bcov_patterns.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov_patterns")
kd_looms_bcov_patterns_ess <- kd_looms_bcov_basic_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov_patterns_log$Sample) / (nrow(kd_looms_bcov_patterns_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov_patterns_ess, ESS < 200)

## Binary covarion, strict clock, heterogene
kd_looms_bcov1111_strict_heterogene_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1111_strict_heterogene/kd-looms_bcov1111_strict_heterogene.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov_patterns")
kd_looms_bcov1111_strict_heterogene_ess <- kd_looms_bcov1111_strict_heterogene_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov1111_strict_heterogene_log$Sample) / (nrow(kd_looms_bcov1111_strict_heterogene_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov1111_strict_heterogene_ess, ESS < 200)


# Check the traces -------------------------------------------------------------

bind_rows(
  kd_looms_bcov1111_relaxed_log,
  # kd_looms_bcov1111_strict_heterogene_log,
  kd_looms_bcov8421_strict_log,
  kd_looms_ctmc1111_strict_log,
  # kd_lgs_bcov_relaxed_byconcept_log,
  kd_lgs_bcov_relaxed_log,
  kd_lgs_bcov_strict_log,
  kd_looms_bcov_basic_log,
  kd_looms_bcov_patterns_log,
  kd_looms_bcov1000_strict_log,
  kd_looms_bcov1111_strict_log,
  # kd_looms_ctmc1111_strict_heterogene_log
) |>
  filter(burnin == FALSE) |>
  ggplot(aes(x = Sample, y = posterior)) +
  geom_line(linewidth = .15, color = ggthemes::few_pal("Dark")(2)[1]) +
  facet_wrap(~data, scales = "free") +
  theme_minimal()


# Prune the taxa in the nexus files and keep common taxa only -------------

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  select(group, lng, lng_group_code) |>
  filter(!is.na(lng_group_code)) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |>
  mutate(group = str_replace_all(group, " ", "_"))
kd_lgs <- read_csv(here("data/kd-lgs/kd-lgs_datapoints.csv")) |>
  select(lng, lng_group_code) |>
  filter(!is.na(lng_group_code)) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng))
kd_lgs_looms <- inner_join(kd_looms, kd_lgs)

kd_looms_pruned <- ReadAsPhyDat(here("data/nexus/kd-looms_1111.nex")) |>
  as_tibble() |>
  select(any_of(kd_lgs_looms$group)) |>
  as.matrix() |>
  t() |>
  MatrixToPhyDat()

kd_lgs_pruned <- ReadAsPhyDat(here("data/nexus/kd-lgs.nex")) |>
  as_tibble() |>
  select(any_of(kd_lgs_looms$lng_label)) |>
  as.matrix() |>
  t() |>
  MatrixToPhyDat()

write_binary_nexus(kd_lgs_pruned, here("data/nexus/kd-lgs_pruned.nex"))
write_binary_nexus(kd_looms_pruned, here("data/nexus/kd-looms_pruned.nex"))

# Merged dataset
kd_merged <- kd_looms_pruned |>
  as_tibble() |>
  rename(all_of(setNames(kd_lgs_looms$group, kd_lgs_looms$lng_label))) |>
  bind_rows(as_tibble(kd_lgs_pruned)) |>
  as.matrix() |>
  t() |>
  MatrixToPhyDat()
write_binary_nexus(kd_merged, here("data/nexus/kd-merged.nex"))

# Filtering: delete 0 columns ---------------------------------------------

# Languages
kd_lngs_pruned_filtered <- kd_lgs_pruned |>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "?"))) |>
  mutate(across(everything(), as.numeric)) |>
  mutate(total = rowSums(across(everything()), na.rm = TRUE)) |>
  filter(total != 0) |>
  select(-total) |>
  mutate(across(everything(), as.character)) |>
  mutate(across(everything(), ~ str_replace_na(.x, "?"))) |>
  as.matrix() |>
  t() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_lngs_pruned_filtered,
  here("data/nexus/kd-lngs_pruned_filtered.nex")
)

# Looms
kd_looms_pruned_filtered <- kd_looms_pruned |>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "?"))) |>
  mutate(across(everything(), as.numeric)) |>
  mutate(total = rowSums(across(everything()), na.rm = TRUE)) |>
  filter(total != 0) |>
  select(-total) |>
  mutate(across(everything(), as.character)) |>
  mutate(across(everything(), ~ str_replace_na(.x, "?"))) |>
  as.matrix() |>
  t() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_pruned_filtered,
  here("data/nexus/kd-looms_pruned_filtered.nex")
)

# Merged data
kd_merged_filtered <- kd_merged |>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "?"))) |>
  mutate(across(everything(), as.numeric)) |>
  mutate(total = rowSums(across(everything()), na.rm = TRUE)) |>
  filter(total != 0) |>
  select(-total) |>
  mutate(across(everything(), as.character)) |>
  mutate(across(everything(), ~ str_replace_na(.x, "?"))) |>
  as.matrix() |>
  t() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_merged_filtered,
  here("data/nexus/kd-merged_filtered.nex")
)
