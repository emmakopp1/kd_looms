library(here)
library(phangorn)
library(TreeTools)
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

kd_looms_weighted <- kd_looms_matrix |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    cols = -Taxon,
    names_to = "code",
    values_to = "value"
  ) |>
  left_join(kd_looms_characters) |>
  arrange(level) |>
  group_by(level) |>
  group_map(~ uncount(.x, weights[.y$level]), .keep = TRUE) |>
  bind_rows() |>
  group_by(Taxon, code, level) |>
  mutate(set = row_number(), code = paste0(code, "_", set)) |>
  ungroup()

kd_looms_matrix_weighted <- kd_looms_weighted |>
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

## Looms, weighted characters, partitioned by level
kd_looms_weighted_partition <- kd_looms_weighted |>
  distinct(code, level) |>
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
  kd_looms_matrix8421,
  here("data/nexus/kd-looms_8421_ht.nex")
)

write_file(str_glue("begin assumptions;\n{kd_looms_weighted_partition}\nend;\n"),
  here("data/nexus/kd-looms_8421_ht.nex"),
  append = TRUE
)

## Looms, unweighted characters, partitioned by level
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

## Looms, basic features only, partitioned by level
kd_looms_matrix_basic_ht <- kd_looms_matrix |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    cols = -Taxon,
    names_to = "code",
    values_to = "value"
  ) |>
  left_join(kd_looms_characters) |>
  filter(type == "Loom basics") |>
  arrange(level) |>
  select(-c(level, type)) |>
  pivot_wider(
    names_from = code,
    values_from = value
  ) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix_basic_ht,
  here("data/nexus/kd-looms_basic_ht.nex")
)

kd_looms_basic_partition <- select(
  kd_looms_matrix,
  filter(kd_looms_characters, type == "Loom basics")$code
) |>
  colnames() |>
  enframe(name = NULL, value = "code") |>
  left_join(kd_looms_characters) |>
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

write_file(str_glue("begin assumptions;\n{kd_looms_basic_partition}\nend;\n"),
  here("data/nexus/kd-looms_basic_ht.nex"),
  append = TRUE
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

## Looms, pattern features only, partitioned by level
kd_looms_matrix_patterns_ht <- kd_looms_matrix |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    cols = -Taxon,
    names_to = "code",
    values_to = "value"
  ) |>
  left_join(kd_looms_characters) |>
  filter(type != "Loom basics") |>
  arrange(level) |>
  select(-c(level, type)) |>
  pivot_wider(
    names_from = code,
    values_from = value
  ) |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()

write_binary_nexus(
  kd_looms_matrix_patterns_ht,
  here("data/nexus/kd-looms_patterns_ht.nex")
)

kd_looms_patterns_partition <- select(
  kd_looms_matrix,
  filter(kd_looms_characters, type != "Loom basics")$code
) |>
  colnames() |>
  enframe(name = NULL, value = "code") |>
  left_join(kd_looms_characters) |>
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

write_file(str_glue("begin assumptions;\n{kd_looms_patterns_partition}\nend;\n"),
  here("data/nexus/kd-looms_patterns_ht.nex"),
  append = TRUE
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
