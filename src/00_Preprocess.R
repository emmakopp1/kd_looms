library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tracerer)
library(tidyverse)

dir.create(here("data/kd-lgs/"))
dir.create(here("data/kd-looms/"))

# Prepare the nexus files ------------------------------------------------------

# Function to force nexus files to specify the symbols to be 0 and 1
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
  pivot_longer(
    cols = !(concept_id | cogid),
    names_to = "Taxon",
    values_to = "value"
  ) |>
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
  select(code, level, type)
kd_looms_matrix <- read_csv(here("data/kd-looms/kd-looms_matrix.csv")) |>
  mutate(Taxon = str_replace_all(Taxon, " ", "_"))

# Looms, all levels, no weighting
dir.create(here("data/kd-looms/kd-looms_bcov1111/"))
kd_looms_matrix1111 <- kd_looms_matrix |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write_binary_nexus(
  kd_looms_matrix1111,
  here("data/kd-looms/kd-looms_bcov1111/kd-looms_1111.nex")
)

# Looms, level 1 characters only
dir.create(here("data/kd-looms/kd-looms_bcov1000/"))
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
  here("data/kd-looms/kd-looms_bcov1000/kd-looms_1000.nex")
)

# Looms, weighted characters
weights <- c(8, 4, 2, 1)
dir.create(here("data/kd-looms/kd-looms_bcov8421/"))

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

kd_looms_matrix8421 <- kd_looms_matrix_weighted |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write_binary_nexus(
  kd_looms_matrix8421,
  here("data/kd-looms/kd-looms_bcov8421/kd-looms_8421.nex")
)

# Looms, 4 variable rates
dir.create(here("data/kd-looms/kd-looms_ctmc4/"))

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
  here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.nex")
)
write_file(str_glue("begin assumptions;\n{kd_looms_partition}\nend;\n"),
  here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.nex"),
  append = TRUE
)

# Looms, basic features only
dir.create(here("data/kd-looms/kd-looms_bcov_basic/"))

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
  here("data/kd-looms/kd-looms_bcov_basic/kd-looms_basic.nex")
)

# Looms, pattern features only
dir.create(here("data/kd-looms/kd-looms_bcov_patterns/"))

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
  here("data/kd-looms/kd-looms_bcov_patterns/kd-looms_patterns.nex")
)


# Add the missing "End;" line at the end of the beast tree files ---------------

write_file("End;",
  here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000.trees"),
  append = TRUE
)
write_file("End;",
  here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees"),
  append = TRUE
)
write_file("End;",
  here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421.trees"),
  append = TRUE
)
write_file("End;",
  here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.trees"),
  append = TRUE
)
write_file("End;",
  here("data/kd-looms/kd-looms_bcov_basic/kd-looms_bcov_basic.trees"),
  append = TRUE
)
write_file("End;",
  here("data/kd-looms/kd-looms_bcov_patterns/kd-looms_bcov_patterns.trees"),
  append = TRUE
)
write_file("End;", here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.trees"),
  append = TRUE
)
zip(
  here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.trees.zip"),
  here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.trees")
)

# Check the ESS values ---------------------------------------------------------

burnin <- .1

# Languages

kd_lgs_bcov_log <- parse_beast_tracelog_file(
  here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-lgs_bcov")
kd_lgs_bcov_ess <- kd_lgs_bcov_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_lgs_bcov_log$Sample) / (nrow(kd_lgs_bcov_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_lgs_bcov_ess, ESS < 200)

# Looms 1000

kd_looms_bcov1000_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov1000")
kd_looms_bcov1000_ess <- kd_looms_bcov1000_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov1000_log$Sample) / (nrow(kd_looms_bcov1000_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov1000_ess, ESS < 200)

# Looms 1111

kd_looms_bcov1111_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov1111")
kd_looms_bcov1111_ess <- kd_looms_bcov1111_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov1111_log$Sample) / (nrow(kd_looms_bcov1111_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov1111_ess, ESS < 200)

# Looms 8421

kd_looms_bcov8421_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd-looms_bcov8421")
kd_looms_bcov84211_ess <- kd_looms_bcov8421_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_bcov8421_log$Sample) / (nrow(kd_looms_bcov8421_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_bcov84211_ess, ESS < 200)

# Looms CTMC 4

kd_looms_ctmc4_log <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.log")
) |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "kd_looms_ctmc4")
kd_looms_ctmc4_ess <- kd_looms_ctmc4_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(kd_looms_ctmc4_log$Sample) / (nrow(kd_looms_ctmc4_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_looms_ctmc4_ess, ESS < 200)

# Looms basic

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

# Looms basic

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


# Check the traces -------------------------------------------------------------

bind_rows(
  kd_lgs_bcov_log,
  kd_looms_bcov1000_log,
  kd_looms_bcov1111_log,
  kd_looms_bcov8421_log,
  kd_looms_ctmc4_log,
  kd_looms_bcov_basic_log,
  kd_looms_bcov_patterns_log
) |>
  filter(burnin == FALSE) |>
  ggplot(aes(x = Sample, y = posterior)) +
  geom_line(linewidth = .15, color = ggthemes::few_pal("Dark")(2)[1]) +
  facet_wrap(~data, scales = "free") +
  theme_minimal()
