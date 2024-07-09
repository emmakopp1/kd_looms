library(here)
library(phangorn)
library(phytools)
library(tracerer)
library(HDInterval)
library(tidyverse)

dir.create(here("output/trees"))
dir.create(here("output/data"))

# Consensus trees -------------------------------------------------------------------------------------------------

kd_lgs_bcov <- read.tree(here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov_trimmed.trees"))
kd_lgs_bcov_cs <- consensus(kd_lgs_bcov, p = .5, rooted = TRUE)
kd_lgs_bcov_cs <- consensus.edges(kd_lgs_bcov, consensus.tree = kd_lgs_bcov_cs, rooted = TRUE)
write.tree(kd_lgs_bcov_cs, here("output/trees/kd-lgs_bcov_consensus.tree"))

kd_looms_bcov1000 <- read.nexus(here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000.trees"))
kd_looms_bcov1000_cs <- consensus(kd_looms_bcov1000, p = .5, rooted = TRUE)
kd_looms_bcov1000_cs <- consensus.edges(kd_looms_bcov1000, consensus.tree = kd_looms_bcov1000_cs, rooted = TRUE)
write.tree(kd_looms_bcov1000_cs, here("output/trees/kd-looms_bcov1000_consensus.tree"))

kd_looms_bcov1111 <- read.nexus(here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111.trees"))
kd_looms_bcov1111_cs <- consensus(kd_looms_bcov1111, p = .5, rooted = TRUE)
kd_looms_bcov1111_cs <- consensus.edges(kd_looms_bcov1111, consensus.tree = kd_looms_bcov1111_cs, rooted = TRUE)
write.tree(kd_looms_bcov1111_cs, here("output/trees/kd-looms_bcov1111_consensus.tree"))

kd_looms_bcov8421 <- read.nexus(here("data/kd-looms/kd-looms_bcov8421/kd-looms_bcov8421_trimmed.trees"))
kd_looms_bcov8421_cs <- consensus(kd_looms_bcov8421, p = .5, rooted = TRUE)
kd_looms_bcov8421_cs <- consensus.edges(kd_looms_bcov8421, consensus.tree = kd_looms_bcov8421_cs, rooted = TRUE)
write.tree(kd_looms_bcov8421_cs, here("output/trees/kd-looms_bcov8421_consensus.tree"))

kd_looms_ctmc4 <- read.tree(here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4_trimmed.trees"))
kd_looms_ctmc4_cs <- consensus(kd_looms_ctmc4, p = .5, rooted = TRUE)
kd_looms_ctmc4_cs <- consensus.edges(kd_looms_ctmc4, consensus.tree = kd_looms_ctmc4_cs, rooted = TRUE)
write.tree(kd_looms_ctmc4_cs, here("output/trees/kd-looms_ctmc4_consensus.tree"))

# kd_looms_ctmc6 <- read.nexus(here("data/kd-looms/kd-looms_ctmc6/kd-looms_ctmc6_trimmed.trees"))
# kd_looms_ctmc6_cs <- consensus(kd_looms_ctmc6, p = .5, rooted = TRUE)
# kd_looms_ctmc6_cs <- consensus.edges(kd_looms_ctmc6, consensus.tree = kd_looms_ctmc6_cs, rooted = TRUE)
# write.tree(kd_looms_ctmc6_cs, here("output/trees/kd-looms_ctmc6_consensus.tree"))


# Ages in language trees --------------------------------------------------

getMRCA_age <- function(tree, tips) {
  tips <- if (is.character(tips)) which(tree$tip.label %in% tips) else tips
  mrca <- ifelse(length(tips) > 1, getMRCA(tree, tips), tips)
  root_age <- max(node.depth.edgelength(tree))
  root_age - node.depth.edgelength(tree)[mrca]
}

kd_lgs_ages <- kd_lgs_bcov |>
  seq_along() |>
  map_df(~ tibble(
    `Kra-Dai` = getMRCA_age(kd_lgs_bcov[[.x]], kd_lgs_bcov[[1]]$tip.label),
    `Kam-Tai` = getMRCA_age(kd_lgs_bcov[[.x]], str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")),
    `Tai-Yay` = getMRCA_age(kd_lgs_bcov[[.x]], str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Tc|Tn|Tsw)"))
  )) |>
  pivot_longer(everything(), names_to = "group", values_to = "age")

write_csv(kd_lgs_ages, here("output/data/kd-lgs_ages.csv"))

kd_lgs_ages_summary <- kd_lgs_ages |>
  group_by(group) |>
  summarise(mean = mean(age), median = median(age), sd = sd(age), HPDI_lower = hdi(age)["lower"], HPDI_upper = hdi(age)["upper"]) |>
  left_join(tribble(
    ~group, ~n,
    "Kam-Tai", length(str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")),
    "Tai-Yay", length(str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Tc|Tn|Tsw)")),
    "Kra-Dai", length(kd_lgs_bcov[[1]]$tip.label)
  )) |>
  relocate(n, .after = group) |>
  arrange(-n) |>
  rename(n_lgs = n)

write_csv(kd_lgs_ages_summary, here("output/data/kd-lgs_ages_summary.csv"))


# Mutation rates --------------------------------------------------------------------------------------------------

kd_looms_characters <- read_csv(here("data/kd-looms/kd-looms_characters.csv")) |>
  select(code, level)
burnin <- .2

kd_looms_mu_bylevel <- parse_beast_tracelog_file(here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.log")) |>
  as_tibble() |>
  select(Sample, mutationRate.s.level1, mutationRate.s.level2, mutationRate.s.level3, mutationRate.s.level4)

kd_looms_mu_bylevel_tb <- kd_looms_mu_bylevel |>
  rowid_to_column() |> 
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  filter(burnin == FALSE) |>
  select(-burnin, -rowid) |>
  pivot_longer(-Sample, names_to = "level", values_to = "rate") |>
  mutate(level = str_remove_all(level, "[^0-9]") |> as.numeric())

write_csv(kd_looms_mu_bylevel_tb, here("output/data/kd-looms_mu_bylevel.csv"))

kd_looms_mu_summary <- kd_looms_mu_bylevel_tb |>
  group_by(level) |>
  summarise(mean = mean(rate), median = median(rate), sd = sd(rate), HPDI_lower = hdi(rate)["lower"], HPDI_upper = hdi(rate)["upper"]) |>
  left_join(count(kd_looms_characters, level)) |>
  relocate(n, .after = level) |>
  rename(n_chars = n)

write_csv(kd_looms_mu_summary, here("output/data/kd-looms_mu_summary.csv"))
