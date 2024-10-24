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

# Languages binary covarion relaxed 
tmp <- tempdir()
unzip(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees.zip"),
  junkpaths = TRUE,
  exdir = here("data/kd-lgs/kd-lgs_bcov_relaxed/")
)

kd_lgs_bcov_relaxed <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed/kd-lgs_bcov_relaxed.trees"))
kd_lgs_bcov_relaxed <- kd_lgs_bcov_relaxed[ceiling(length(kd_lgs_bcov_relaxed) * burnin):length(kd_lgs_bcov_relaxed)]
kd_lgs_bcov_relaxed_cs <- consensus(kd_lgs_bcov_relaxed, p = .5, rooted = TRUE)
kd_lgs_bcov_relaxed_cs <- consensus.edges(kd_lgs_bcov_relaxed,
  consensus.tree = kd_lgs_bcov_relaxed_cs,
  rooted = TRUE
)
kd_lgs_bcov_relaxed_cs$root.edge.length <- 0
write.tree(
  kd_lgs_bcov_relaxed_cs,
  here("output/trees/kd-lgs_bcov_relaxed_consensus.tree")
)
unlink(tmp, recursive = TRUE)

# Language binary covarion relaxed by concept
kd_lgs_bcov_relaxed_heterogene <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed_byconcept/kd-lgs_bcov_byconcept.trees"))
kd_lgs_bcov_relaxed_heterogene <- kd_lgs_bcov_relaxed_heterogene[ceiling(length(kd_lgs_bcov_relaxed_heterogene) * burnin):length(kd_lgs_bcov_relaxed_heterogene)]
kd_lgs_bcov_relaxed_heterogene_cs <- consensus(kd_lgs_bcov_relaxed_heterogene, p = .5, rooted = TRUE)
kd_lgs_bcov_relaxed_heterogene_cs <- consensus.edges(
  kd_lgs_bcov_relaxed_heterogene,
  consensus.tree = kd_lgs_bcov_relaxed_heterogene_cs,
  rooted = TRUE
)

kd_lgs_bcov_relaxed_heterogene_cs$root.edge.length <- 0
write.tree(
  kd_lgs_bcov_relaxed_heterogene_cs,
  here("output/trees/kd-lgs_bcov_relaxed_byconcept_consensus.tree")
)

# Language binary covarion strict 
kd_lgs_bcov_strict <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_strict/kd-lgs_bcov_strict.trees"))
kd_lgs_bcov_strict <- kd_lgs_bcov_strict[ceiling(length(kd_lgs_bcov_strict) * burnin):length(kd_lgs_bcov_strict)]
kd_lgs_bcov_strict_cs <- consensus(kd_lgs_bcov_strict, p = .5, rooted = TRUE)
kd_lgs_bcov_strict_cs <- consensus.edges(
  kd_lgs_bcov_strict,
  consensus.tree = kd_lgs_bcov_strict_cs,
  rooted = TRUE
)

kd_lgs_bcov_strict_cs$root.edge.length <- 0
write.tree(
  kd_lgs_bcov_strict_cs,
  here("output/trees/kd-lgs_bcov_strict_consensus.tree")
)

# Looms
## Looms, level 1 characters only
kd_looms_bcov1000_strict <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1000_strict/kd-looms_bcov1000_strict.trees")
)
kd_looms_bcov1000_strict <- kd_looms_bcov1000_strict[ceiling(length(kd_looms_bcov1000_strict) * burnin):length(kd_looms_bcov1000_strict)]
kd_looms_bcov1000_strict_cs <- consensus(kd_looms_bcov1000_strict, p = .5, rooted = TRUE)
kd_looms_bcov1000_strict_cs <- consensus.edges(kd_looms_bcov1000_strict,
  consensus.tree = kd_looms_bcov1000_strict_cs,
  rooted = TRUE
)
write.tree(
  kd_looms_bcov1000_strict_cs,
  here("output/trees/kd-looms_bcov1000_strict_consensus.tree")
)

## Looms, all levels, no weighting
kd_looms_bcov1111_strict <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_strict/kd-looms_bcov1111_strict.trees")
)
kd_looms_bcov1111_strict <- kd_looms_bcov1111_strict[ceiling(length(kd_looms_bcov1111_strict) * burnin):length(kd_looms_bcov1111_strict)]
kd_looms_bcov1111_strict_cs <- consensus(kd_looms_bcov1111_strict, p = .5, rooted = TRUE)
kd_looms_bcov1111_strict_cs <- consensus.edges(kd_looms_bcov1111_strict,
  consensus.tree = kd_looms_bcov1111_strict_cs,
  rooted = TRUE
)
write.tree(
  kd_looms_bcov1111_strict_cs,
  here("output/trees/kd_looms_bcov1111_strict_consensus.tree")
)

## Looms, weighted characters
kd_looms_bcov8421_strict <- read.nexus(
  here("data/kd-looms/kd-looms_bcov8421_strict/kd-looms_bcov8421_strict.trees")
)
kd_looms_bcov8421_strict <- kd_looms_bcov8421_strict[ceiling(length(kd_looms_bcov8421_strict) * burnin):length(kd_looms_bcov8421_strict)]
kd_looms_bcov8421_strict_cs <- consensus(kd_looms_bcov8421_strict, p = .5, rooted = TRUE)
kd_looms_bcov8421_strict_cs <- consensus.edges(kd_looms_bcov8421_strict,
  consensus.tree = kd_looms_bcov8421_strict_cs,
  rooted = TRUE
)
write.tree(
  kd_looms_bcov8421_strict_cs,
  here("output/trees/kd-looms_bcov8421_strict_consensus.tree")
)

## Looms, binary covarion strict heterogene (bizzare car différent de ctmc4rates)
kd_looms_bcov1111_strict_heterogene <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_strict_heterogene/kd-looms_bcov1111_strict_heterogene.trees")
)
kd_looms_bcov1111_strict_heterogene <- kd_looms_bcov1111_strict_heterogene[ceiling(length(kd_looms_bcov1111_strict_heterogene) * burnin):length(kd_looms_bcov1111_strict_heterogene)]
kd_looms_bcov1111_strict_heterogene_cs <- consensus(kd_looms_bcov1111_strict_heterogene, p = .5, rooted = TRUE)
kd_looms_bcov1111_strict_heterogene_cs <- consensus.edges(kd_looms_bcov1111_strict_heterogene,
  consensus.tree = kd_looms_bcov1111_strict_heterogene_cs,
  rooted = TRUE
)
write.tree(
  kd_looms_bcov1111_strict_heterogene_cs,
  here("output/trees/kd-looms_bcov_strict_heterogene_consensus.tree")
)

## Looms, basic features only
kd_looms_bcov_basic <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_basic/kd-looms_bcov_basic.trees")
)
kd_looms_bcov_basic <- kd_looms_bcov_basic[ceiling(length(kd_looms_bcov_basic) * burnin):length(kd_looms_bcov_basic)]
kd_looms_bcov_basic_cs <- consensus(kd_looms_bcov_basic, p = .5, rooted = TRUE)
kd_looms_bcov_basic_cs <- consensus.edges(kd_looms_bcov_basic,
  consensus.tree = kd_looms_bcov_basic_cs,
  rooted = TRUE
)
write.tree(
  kd_looms_bcov_basic_cs,
  here("output/trees/kd-looms_bcov_basic_consensus.tree")
)

## Looms, pattern features only
kd_looms_bcov_patterns <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_patterns/kd-looms_bcov_patterns.trees")
)
kd_looms_bcov_patterns <- kd_looms_bcov_patterns[ceiling(length(kd_looms_bcov_patterns) * burnin):length(kd_looms_bcov_patterns)]
kd_looms_bcov_patterns_cs <- consensus(kd_looms_bcov_patterns, p = .5, rooted = TRUE)
kd_looms_bcov_patterns_cs <- consensus.edges(kd_looms_bcov_patterns,
  consensus.tree = kd_looms_bcov_patterns_cs,
  rooted = TRUE
)
write.tree(
  kd_looms_bcov_patterns_cs,
  here("output/trees/kd-looms_bcov_patterns_consensus.tree")
)


## Looms, pattern features only
kd_looms_bcov1111_relaxed <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_relaxed/kd-looms_bcov1111_relaxed.trees")
)
kd_looms_bcov1111_relaxed <- kd_looms_bcov1111_relaxed[ceiling(length(kd_looms_bcov1111_relaxed) * burnin):length(kd_looms_bcov1111_relaxed)]
kd_looms_bcov1111_relaxed_cs <- consensus(kd_looms_bcov1111_relaxed, p = .5, rooted = TRUE)
kd_looms_bcov1111_relaxed_cs <- consensus.edges(
  kd_looms_bcov1111_relaxed,
  consensus.tree = kd_looms_bcov1111_relaxed_cs,
  rooted = TRUE
)

write.tree(
  kd_looms_bcov1111_relaxed_cs,
  here("output/trees/kd-looms_bcov1111_relaxed_consensus.tree")
)

## Looms, ctmc strict
kd_looms_ctmc1111_strict <- read.nexus(
  here("data/kd-looms/kd-looms_ctmc1111_strict/kd-looms_ctmc1111_strict.trees")
)
kd_looms_ctmc1111_strict <- kd_looms_ctmc1111_strict[ceiling(length(kd_looms_ctmc1111_strict) * burnin):length(kd_looms_ctmc1111_strict)]
kd_looms_ctmc1111_strict_cs <- consensus(kd_looms_ctmc1111_strict, p = .5, rooted = TRUE)
kd_looms_ctmc1111_strict_cs <- consensus.edges(
  kd_looms_ctmc1111_strict,
  consensus.tree = kd_looms_ctmc1111_strict_cs,
  rooted = TRUE
)

write.tree(
  kd_looms_ctmc1111_strict_cs,
  here("output/trees/kd-looms_ctmc1111_strict_consensus.tree")
)

## Looms, ctmc strict heterogene
kd_looms_ctmc1111_strict_heterogene <- read.nexus(
  here("data/kd-looms/kd-looms_ctmc1111_strict_heterogene/kd-looms_ctmc1111_strict_heterogene.trees")
)
kd_looms_ctmc1111_strict_heterogene <- kd_looms_ctmc1111_strict_heterogene[ceiling(length(kd_looms_ctmc1111_strict_heterogene) * burnin):length(kd_looms_ctmc1111_strict_heterogene)]
kd_looms_ctmc1111_strict_heterogene_cs <- consensus(kd_looms_ctmc1111_strict_heterogene, p = .5, rooted = TRUE)
kd_looms_ctmc1111_strict_heterogene_cs <- consensus.edges(
  kd_looms_ctmc1111_strict_heterogene,
  consensus.tree = kd_looms_ctmc1111_strict_heterogene_cs,
  rooted = TRUE
)

write.tree(
  kd_looms_ctmc1111_strict_heterogene_cs,
  here("output/trees/kd-looms_ctmc1111_strict_heterogene_consensus.tree")
)

# Ages in language trees -------------------------------------------------------

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
    `Kam-Tai` = getMRCA_age(
      kd_lgs_bcov[[.x]],
      str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")
    ),
    `Tai-Yay` = getMRCA_age(
      kd_lgs_bcov[[.x]],
      str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Tc|Tn|Tsw)")
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
    "Kam-Tai", length(str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Ks|Tc|Tn|Tsw)")),
    "Tai-Yay", length(str_subset(kd_lgs_bcov[[1]]$tip.label, "^(Tc|Tn|Tsw)")),
    "Kra-Dai", length(kd_lgs_bcov[[1]]$tip.label)
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
  here("data/kd-looms/kd-looms_ctmc4/kd-looms_ctmc4.log")
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
