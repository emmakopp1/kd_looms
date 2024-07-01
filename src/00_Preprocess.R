library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tracerer)
library(tidyverse)


# Prepare the nexus files -------------------------------------------------

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
  select(-c(level, set)) |>
  pivot_wider(
    names_from = c(code),
    values_from = value
  )

kd_looms_matrix8421 <- kd_looms_matrix_weighted |>
  column_to_rownames("Taxon") |>
  as.matrix() |>
  MatrixToPhyDat()
write_binary_nexus(kd_looms_matrix8421, here("data/kd-looms/kd-looms_bcov8421/kd-looms_8421.nex"))

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



# Remove some burn-in, check the ESS values, and trim the raw tree files ------------------------------------------

burnin <- .2
k <- 5000

# Languages

kd_lgs_bcov_log <- parse_beast_tracelog_file(here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.log")) |>
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

kd_lgs_bcov <- read.nexus(here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov.trees"))
kd_lgs_bcov_trimmed <- kd_lgs_bcov[ceiling(length(kd_lgs_bcov) * burnin):length(kd_lgs_bcov)]
kd_lgs_bcov_trimmed <- kd_lgs_bcov_trimmed[seq(from = 2, to = length(kd_lgs_bcov_trimmed), by = round(length(kd_lgs_bcov_trimmed) / k))]
write.tree(kd_lgs_bcov_trimmed, here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov_trimmed.trees"))

# Looms 1000

kd_looms_bcov1000_log <- parse_beast_tracelog_file(here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000.log")) |>
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

kd_looms_bcov1000 <- read.nexus(here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000.trees"))
kd_looms_bcov1000_trimmed <- kd_looms_bcov1000[ceiling(length(kd_looms_bcov1000) * burnin):length(kd_looms_bcov1000)]
kd_looms_bcov1000_trimmed <- kd_looms_bcov1000_trimmed[seq(from = 1, to = length(kd_looms_bcov1000_trimmed), by = length(kd_looms_bcov1000_trimmed) / k)]
write.tree(kd_looms_bcov1000_trimmed, here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000_trimmed.trees"))

# Looms 1111

kd_looms_bcov1111_log <- parse_beast_tracelog_file(here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.log")) |>
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

kd_looms_bcov1111 <- read.nexus(here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees"))
kd_looms_bcov1111_trimmed <- kd_looms_bcov1111[ceiling(length(kd_looms_bcov1111) * burnin):length(kd_looms_bcov1111)]
kd_looms_bcov1111_trimmed <- kd_looms_bcov1111_trimmed[seq(from = 1, to = length(kd_looms_bcov1111_trimmed), by = length(kd_looms_bcov1111_trimmed) / k)]
write.tree(kd_looms_bcov1111_trimmed, here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111_trimmed.trees"))

# Looms 8421

kd_looms_bcov8421_log <- parse_beast_tracelog_file(here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421.log")) |>
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

kd_looms_bcov8421 <- read.nexus(here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421.trees"))
kd_looms_bcov8421_trimmed <- kd_looms_bcov8421[ceiling(length(kd_looms_bcov8421) * burnin):length(kd_looms_bcov8421)]
kd_looms_bcov8421_trimmed <- kd_looms_bcov8421_trimmed[seq(from = 1, to = length(kd_looms_bcov8421_trimmed), by = length(kd_looms_bcov8421_trimmed) / k)]
write.tree(kd_looms_bcov8421_trimmed, here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421_trimmed.trees"))

# Looms CTMC 4

kd_looms_ctmc4_log <- parse_beast_tracelog_file(here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.log")) |>
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

kd_looms_ctmc4 <- read.nexus(here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.trees"))
kd_looms_ctmc4_trimmed <- kd_looms_ctmc4[ceiling(length(kd_looms_ctmc4) * burnin):length(kd_looms_ctmc4)]
kd_looms_ctmc4_trimmed <- kd_looms_ctmc4_trimmed[seq(from = 1, to = length(kd_looms_ctmc4_trimmed), by = length(kd_looms_ctmc4_trimmed) / k)]
write.tree(kd_looms_ctmc4_trimmed, here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4_trimmed.trees"))

# Looms CTMC 6

kd_looms_ctmc6 <- read.nexus(here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6.trees"))
kd_looms_ctmc6_trimmed <- kd_looms_ctmc6[ceiling(length(kd_looms_ctmc6) * burnin):length(kd_looms_ctmc6)]
kd_looms_ctmc6_trimmed <- kd_looms_ctmc6_trimmed[seq(from = 1, to = length(kd_looms_ctmc6_trimmed), by = length(kd_looms_ctmc6_trimmed) / k)]
write.tree(kd_looms_ctmc6_trimmed, here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6_trimmed.trees"))


# Check the traces ------------------------------------------------------------------------------------------------

bind_rows(
  kd_lgs_bcov_log,
  kd_looms_bcov1000_log,
  kd_looms_bcov1111_log,
  kd_looms_bcov8421_log,
  kd_looms_ctmc4_log
) |>
  filter(burnin == FALSE) |>
  ggplot(aes(x = Sample, y = posterior)) +
  geom_line(linewidth = .15, color = ggthemes::few_pal("Dark")(2)[1]) +
  facet_wrap(~data, scales = "free") +
  theme_minimal()
