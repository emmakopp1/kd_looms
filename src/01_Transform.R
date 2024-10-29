library(here)
library(phangorn)
library(phytools)
library(tracerer)
library(HDInterval)
library(tidyverse)

dir.create(here("output/trees"))
dir.create(here("output/data"))


# Consensus trees --------------------------------------------------------------
burnin <- .1

# Languages binary covarion relaxed uniform rate
tmp <- tempdir()
unzip(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees.zip"),
  junkpaths = TRUE,
  exdir = here("data/kd-lgs/kd-lgs_bcov_relaxed/")
)

kd_lgs_bcov_relaxed_uni <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"))
kd_lgs_bcov_relaxed_uni <- kd_lgs_bcov_relaxed_uni[ceiling(length(kd_lgs_bcov_relaxed_uni) * burnin):length(kd_lgs_bcov_relaxed_uni)]
kd_lgs_bcov_relaxed_uni_cs <- consensus(kd_lgs_bcov_relaxed_uni, p = .5, rooted = TRUE)
kd_lgs_bcov_relaxed_uni_cs <- consensus.edges(kd_lgs_bcov_relaxed_uni,
  consensus.tree = kd_lgs_bcov_relaxed_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_lgs_bcov_relaxed_uni_cs)) {
  kd_lgs_bcov_relaxed_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_lgs_bcov_relaxed_uni_cs,
  here("output/trees/kd-lgs_bcov_relaxed_uni_consensus.tree")
)
unlink(tmp, recursive = TRUE)

# Languages binary covarion relaxed heterogeneous rate
tmp <- tempdir()
unzip(here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/kd-lgs_bcov_relaxed_byconcept.trees.zip"),
  junkpaths = TRUE,
  exdir = here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/")
)
# unzip(here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/kd-lgs_bcov_relaxed_byconcept.log.zip"),
#   junkpaths = TRUE,
#   exdir = here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/")
# )

kd_lgs_bcov_relaxed_ht <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/kd-lgs_bcov_relaxed_byconcept.trees"))
kd_lgs_bcov_relaxed_ht <- kd_lgs_bcov_relaxed_ht[ceiling(length(kd_lgs_bcov_relaxed_ht) * burnin):length(kd_lgs_bcov_relaxed_ht)]
kd_lgs_bcov_relaxed_ht_cs <- consensus(kd_lgs_bcov_relaxed_ht, p = .5, rooted = TRUE)
kd_lgs_bcov_relaxed_ht_cs <- consensus.edges(
  kd_lgs_bcov_relaxed_ht,
  consensus.tree = kd_lgs_bcov_relaxed_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_lgs_bcov_relaxed_ht_cs)) {
  kd_lgs_bcov_relaxed_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_lgs_bcov_relaxed_ht_cs,
  here("output/trees/kd-lgs_bcov_relaxed_ht_consensus.tree")
)

# Language binary covarion strict uniform rate
unzip(here("data/kd-lgs/kd-lgs_bcov_strict/kd-lgs_bcov_strict.trees.zip"),
  junkpaths = TRUE,
  exdir = here("data/kd-lgs/kd-lgs_bcov_strict/")
)
kd_lgs_bcov_strict_uni <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_strict/kd-lgs_bcov_strict.trees"))
kd_lgs_bcov_strict_uni <- kd_lgs_bcov_strict_uni[ceiling(length(kd_lgs_bcov_strict_uni) * burnin):length(kd_lgs_bcov_strict_uni)]
kd_lgs_bcov_strict_uni_cs <- consensus(kd_lgs_bcov_strict_uni, p = .5, rooted = TRUE)
kd_lgs_bcov_strict_uni_cs <- consensus.edges(
  kd_lgs_bcov_strict_uni,
  consensus.tree = kd_lgs_bcov_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_lgs_bcov_strict_uni_cs)) {
  kd_lgs_bcov_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_lgs_bcov_strict_uni_cs,
  here("output/trees/kd-lgs_bcov_strict_uni_consensus.tree")
)

# Looms

## Looms, binary covarion, level 1 characters only, strict uniform rate
kd_looms_bcov1000_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1000_strict/kd-looms_bcov1000_strict.trees")
)
kd_looms_bcov1000_strict_uni <- kd_looms_bcov1000_strict_uni[ceiling(length(kd_looms_bcov1000_strict_uni) * burnin):length(kd_looms_bcov1000_strict_uni)]
kd_looms_bcov1000_strict_uni_cs <- consensus(kd_looms_bcov1000_strict_uni, p = .5, rooted = TRUE)
kd_looms_bcov1000_strict_uni_cs <- consensus.edges(kd_looms_bcov1000_strict_uni,
  consensus.tree = kd_looms_bcov1000_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov1000_strict_uni_cs)) {
  kd_looms_bcov1000_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov1000_strict_uni_cs,
  here("output/trees/kd-looms_bcov1000_strict_uni_consensus.tree")
)

## Looms, binary covarion, all levels, no weighting, strict uniform rate
kd_looms_bcov1111_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_strict/kd-looms_bcov1111_strict.trees")
)
kd_looms_bcov1111_strict_uni <- kd_looms_bcov1111_strict_uni[ceiling(length(kd_looms_bcov1111_strict_uni) * burnin):length(kd_looms_bcov1111_strict_uni)]
kd_looms_bcov1111_strict_uni_cs <- consensus(kd_looms_bcov1111_strict_uni, p = .5, rooted = TRUE)
kd_looms_bcov1111_strict_uni_cs <- consensus.edges(kd_looms_bcov1111_strict_uni,
  consensus.tree = kd_looms_bcov1111_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov1111_strict_uni_cs)) {
  kd_looms_bcov1111_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov1111_strict_uni_cs,
  here("output/trees/kd_looms_bcov1111_strict_uni_consensus.tree")
)

## Looms, binary covarion, all levels, no weighting, strict heterogeneous rate
kd_looms_bcov1111_strict_ht <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_strict_heterogene/kd-looms_bcov1111_strict_heterogene.trees")
)
kd_looms_bcov1111_strict_ht <- kd_looms_bcov1111_strict_ht[ceiling(length(kd_looms_bcov1111_strict_ht) * burnin):length(kd_looms_bcov1111_strict_ht)]
kd_looms_bcov1111_strict_ht_cs <- consensus(kd_looms_bcov1111_strict_ht, p = .5, rooted = TRUE)
kd_looms_bcov1111_strict_ht_cs <- consensus.edges(kd_looms_bcov1111_strict_ht,
  consensus.tree = kd_looms_bcov1111_strict_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov1111_strict_ht_cs)) {
  kd_looms_bcov1111_strict_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov1111_strict_ht_cs,
  here("output/trees/kd-looms_bcov_strict_ht_consensus.tree")
)

## Looms, binary covarion, all levels, no weighting, relaxed uniform rate
kd_looms_bcov1111_relaxed_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_relaxed/kd-looms_bcov1111_relaxed.trees")
)
kd_looms_bcov1111_relaxed_uni <- kd_looms_bcov1111_relaxed_uni[ceiling(length(kd_looms_bcov1111_relaxed_uni) * burnin):length(kd_looms_bcov1111_relaxed_uni)]
kd_looms_bcov1111_relaxed_uni_cs <- consensus(kd_looms_bcov1111_relaxed_uni, p = .5, rooted = TRUE)
kd_looms_bcov1111_relaxed_uni_cs <- consensus.edges(
  kd_looms_bcov1111_relaxed_uni,
  consensus.tree = kd_looms_bcov1111_relaxed_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov1111_relaxed_uni_cs)) {
  kd_looms_bcov1111_relaxed_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov1111_relaxed_uni_cs,
  here("output/trees/kd-looms_bcov1111_relaxed_uni_consensus.tree")
)

## Looms, binary covarion, weighted characters, strict uniform rate
kd_looms_bcov8421_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov8421_strict/kd-looms_bcov8421_strict.trees")
)
kd_looms_bcov8421_strict_uni <- kd_looms_bcov8421_strict_uni[ceiling(length(kd_looms_bcov8421_strict_uni) * burnin):length(kd_looms_bcov8421_strict_uni)]
kd_looms_bcov8421_strict_uni_cs <- consensus(kd_looms_bcov8421_strict_uni, p = .5, rooted = TRUE)
kd_looms_bcov8421_strict_uni_cs <- consensus.edges(kd_looms_bcov8421_strict_uni,
  consensus.tree = kd_looms_bcov8421_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov8421_strict_uni_cs)) {
  kd_looms_bcov8421_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov8421_strict_uni_cs,
  here("output/trees/kd-looms_bcov8421_strict_uni_consensus.tree")
)

## Looms, binary covarion, basic features only, strict uniform rate
kd_looms_bcov_basic_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_basic/kd-looms_bcov_basic.trees")
)
kd_looms_bcov_basic_strict_uni <- kd_looms_bcov_basic_strict_uni[ceiling(length(kd_looms_bcov_basic_strict_uni) * burnin):length(kd_looms_bcov_basic_strict_uni)]
kd_looms_bcov_basic_strict_uni_cs <- consensus(kd_looms_bcov_basic_strict_uni, p = .5, rooted = TRUE)
kd_looms_bcov_basic_strict_uni_cs <- consensus.edges(kd_looms_bcov_basic_strict_uni,
  consensus.tree = kd_looms_bcov_basic_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov_basic_strict_uni_cs)) {
  kd_looms_bcov_basic_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov_basic_strict_uni_cs,
  here("output/trees/kd-looms_bcov_basic_strict_uni_consensus.tree")
)

## Looms, binary covarion, pattern features only, strict uniform rate
kd_looms_bcov_patterns_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_patterns/kd-looms_bcov_patterns.trees")
)
kd_looms_bcov_patterns_strict_uni <- kd_looms_bcov_patterns_strict_uni[ceiling(length(kd_looms_bcov_patterns_strict_uni) * burnin):length(kd_looms_bcov_patterns_strict_uni)]
kd_looms_bcov_patterns_strict_uni_cs <- consensus(kd_looms_bcov_patterns_strict_uni, p = .5, rooted = TRUE)
kd_looms_bcov_patterns_strict_uni_cs <- consensus.edges(kd_looms_bcov_patterns_strict_uni,
  consensus.tree = kd_looms_bcov_patterns_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov_patterns_strict_uni_cs)) {
  kd_looms_bcov_patterns_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov_patterns_strict_uni_cs,
  here("output/trees/kd-looms_bcov_patterns_strict_uni_consensus.tree")
)

## Looms, ctmc strict uniform rate
kd_looms_ctmc1111_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_ctmc1111_strict/kd-looms_ctmc1111_strict.trees")
)
kd_looms_ctmc1111_strict_uni <- kd_looms_ctmc1111_strict_uni[ceiling(length(kd_looms_ctmc1111_strict_uni) * burnin):length(kd_looms_ctmc1111_strict_uni)]
kd_looms_ctmc1111_strict_uni_cs <- consensus(kd_looms_ctmc1111_strict_uni, p = .5, rooted = TRUE)
kd_looms_ctmc1111_strict_uni_cs <- consensus.edges(
  kd_looms_ctmc1111_strict_uni,
  consensus.tree = kd_looms_ctmc1111_strict_uni_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_ctmc1111_strict_uni_cs)) {
  kd_looms_ctmc1111_strict_uni_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_ctmc1111_strict_uni_cs,
  here("output/trees/kd-looms_ctmc1111_strict_uni_consensus.tree")
)

## Looms, ctmc strict heterogeneous rate
kd_looms_ctmc1111_strict_ht <- read.nexus(
  here("data/kd-looms/kd-looms_ctmc1111_strict_heterogene/kd-looms_ctmc1111_strict_heterogene.trees")
)
kd_looms_ctmc1111_strict_ht <- kd_looms_ctmc1111_strict_ht[ceiling(length(kd_looms_ctmc1111_strict_ht) * burnin):length(kd_looms_ctmc1111_strict_ht)]
kd_looms_ctmc1111_strict_ht_cs <- consensus(kd_looms_ctmc1111_strict_ht, p = .5, rooted = TRUE)
kd_looms_ctmc1111_strict_ht_cs <- consensus.edges(
  kd_looms_ctmc1111_strict_ht,
  consensus.tree = kd_looms_ctmc1111_strict_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_ctmc1111_strict_ht_cs)) {
  kd_looms_ctmc1111_strict_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_ctmc1111_strict_ht_cs,
  here("output/trees/kd-looms_ctmc1111_strict_ht_consensus.tree")
)

# Ages in language trees -------------------------------------------------------

getMRCA_age <- function(tree, tips) {
  tips <- if (is.character(tips)) which(tree$tip.label %in% tips) else tips
  mrca <- ifelse(length(tips) > 1, getMRCA(tree, tips), tips)
  root_age <- max(node.depth.edgelength(tree))
  root_age - node.depth.edgelength(tree)[mrca]
}

kd_lgs_ages <- kd_lgs_bcov_relaxed_ht |>
  seq_along() |>
  map_df(~ tibble(
    `Kra-Dai` = getMRCA_age(
      kd_lgs_bcov_relaxed_ht[[.x]],
      kd_lgs_bcov_relaxed_ht[[1]]$tip.label
    ),
    `Kam-Tai` = getMRCA_age(
      kd_lgs_bcov_relaxed_ht[[.x]],
      str_subset(kd_lgs_bcov_relaxed_ht[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")
    ),
    `Tai-Yay` = getMRCA_age(
      kd_lgs_bcov_relaxed_ht[[.x]],
      str_subset(kd_lgs_bcov_relaxed_ht[[1]]$tip.label, "^(Tc|Tn|Tsw)")
    )
  )) |>
  pivot_longer(everything(), names_to = "group", values_to = "age")

write_csv(kd_lgs_ages, here("output/data/kd-lgs_ages.csv"))

kd_lgs_ages_summary <- kd_lgs_ages |>
  group_by(group) |>
  summarise(
    mean = mean(age),
    median = median(age),
    sd = sd(age),
    HPDI_lower = hdi(age)["lower"],
    HPDI_upper = hdi(age)["upper"]
  ) |>
  left_join(tribble(
    ~group, ~n,
    "Kam-Tai", length(str_subset(kd_lgs_bcov_relaxed_ht[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")),
    "Tai-Yay", length(str_subset(kd_lgs_bcov_relaxed_ht[[1]]$tip.label, "^(Tc|Tn|Tsw)")),
    "Kra-Dai", length(kd_lgs_bcov_relaxed_ht[[1]]$tip.label)
  )) |>
  relocate(n, .after = group) |>
  arrange(-n) |>
  rename(n_lgs = n)
write_csv(kd_lgs_ages_summary, here("output/data/kd-lgs_ages_summary.csv"))


# Mutation rates ---------------------------------------------------------------
burnin <- .1

# Languages
kd_lgs_concepts <- read_csv(here("data/kd-lgs/kd-lgs_lx.csv")) |>
  count(concept_id)

kd_lgs_mu_byconcept <- parse_beast_tracelog_file(
  here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/kd-lgs_bcov_relaxed_byconcept.log")
) |>
  as_tibble() |>
  select(Sample, starts_with("mutationRate.s.concept_"))

kd_lgs_mu_byconcept_tb <- kd_lgs_mu_byconcept |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  filter(burnin == FALSE) |>
  select(-burnin, -rowid) |>
  pivot_longer(-Sample, names_to = "concept", values_to = "rate") |>
  mutate(concept = str_remove(concept, "mutationRate.s.concept_"))
# write_csv(kd_lgs_mu_byconcept_tb, here("output/data/kd-lgs_mu_byconcept.csv"))

kd_lgs_mu_summary <- kd_lgs_mu_byconcept_tb |>
  group_by(concept) |>
  summarise(
    mean = mean(rate),
    median = median(rate),
    sd = sd(rate),
    HPDI_lower = hdi(rate)["lower"],
    HPDI_upper = hdi(rate)["upper"]
  ) |>
  mutate(concept_id = str_remove(concept, "_.+$")) |>
  left_join(kd_lgs_concepts) |>
  relocate(n, .after = concept) |>
  rename(n_cogsets = n) |>
  select(-concept_id)
write_csv(kd_lgs_mu_summary, here("output/data/kd-lgs_mu_summary.csv"))

# Looms
kd_looms_characters <- read_csv(here("data/kd-looms/kd-looms_characters.csv")) |>
  select(code, level)

kd_looms_mu_bylevel <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1111_strict_heterogene/kd-looms_bcov1111_strict_heterogene.log")
) |>
  as_tibble() |>
  select(
    Sample,
    mutationRate.s.level1,
    mutationRate.s.level2,
    mutationRate.s.level3,
    mutationRate.s.level4
  )

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
  summarise(
    mean = mean(rate),
    median = median(rate),
    sd = sd(rate),
    HPDI_lower = hdi(rate)["lower"],
    HPDI_upper = hdi(rate)["upper"]
  ) |>
  left_join(count(kd_looms_characters, level)) |>
  relocate(n, .after = level) |>
  rename(n_chars = n)
write_csv(kd_looms_mu_summary, here("output/data/kd-looms_mu_summary.csv"))
