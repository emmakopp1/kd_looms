library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tracerer)
library(HDInterval)
library(ggtree)
library(tidyverse)

dir.create(here("output"))
dir.create(here("output/trees"))
dir.create(here("output/data"))


# Consensus trees --------------------------------------------------------------
burnin <- .1

# Languages binary covarion relaxed uniform rate
kd_lgs_bcov_relaxed_uni <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed_uni/kd-lgs_bcov_relaxed_uni.trees"))
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

# Language binary covarion strict uniform rate
kd_lgs_bcov_strict_uni <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_strict_uni/kd-lgs_bcov_strict_uni.trees"))
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

# Language binary covarion strict heterogeneous by part of speech
kd_lgs_bcov_strict_ht_pos <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_strict_ht_pos/kd-lgs_bcov_strict_ht_pos.trees"))
kd_lgs_bcov_strict_ht_pos <- kd_lgs_bcov_strict_ht_pos[ceiling(length(kd_lgs_bcov_strict_ht_pos) * burnin):length(kd_lgs_bcov_strict_ht_pos)]
kd_lgs_bcov_strict_ht_pos_cs <- consensus(kd_lgs_bcov_strict_ht_pos, p = .5, rooted = TRUE)
kd_lgs_bcov_strict_ht_pos_cs <- consensus.edges(
  kd_lgs_bcov_strict_ht_pos,
  consensus.tree = kd_lgs_bcov_strict_ht_pos_cs,
  rooted = TRUE
)
if (!is.rooted(kd_lgs_bcov_strict_ht_pos_cs)) {
  kd_lgs_bcov_strict_ht_pos_cs$root.edge.length <- 0
}
write.tree(
  kd_lgs_bcov_strict_ht_pos_cs,
  here("output/trees/kd-lgs_bcov_strict_ht_pos_consensus.tree")
)

# Language binary covarion relaxed heterogeneous by part of speech
kd_lgs_bcov_relaxed_ht_pos <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed_ht_pos/kd-lgs_bcov_relaxed_ht_pos.trees"))
kd_lgs_bcov_relaxed_ht_pos <- kd_lgs_bcov_relaxed_ht_pos[ceiling(length(kd_lgs_bcov_relaxed_ht_pos) * burnin):length(kd_lgs_bcov_relaxed_ht_pos)]
kd_lgs_bcov_relaxed_ht_pos_cs <- consensus(kd_lgs_bcov_relaxed_ht_pos, p = .5, rooted = TRUE)
kd_lgs_bcov_relaxed_ht_pos_cs <- consensus.edges(
  kd_lgs_bcov_relaxed_ht_pos,
  consensus.tree = kd_lgs_bcov_relaxed_ht_pos_cs,
  rooted = TRUE
)
if (!is.rooted(kd_lgs_bcov_relaxed_ht_pos_cs)) {
  kd_lgs_bcov_relaxed_ht_pos_cs$root.edge.length <- 0
}
write.tree(
  kd_lgs_bcov_relaxed_ht_pos_cs,
  here("output/trees/kd-lgs_bcov_relaxed_ht_pos_consensus.tree")
)

# Looms

## Looms, binary covarion, level 1 characters only, strict uniform rate
kd_looms_bcov1000_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1000_strict_uni/kd-looms_bcov1000_strict_uni.trees")
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
  here("data/kd-looms/kd-looms_bcov1111_strict_uni/kd-looms_bcov1111_strict_uni.trees")
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
  here("output/trees/kd-looms_bcov1111_strict_uni_consensus.tree")
)

## Looms, binary covarion, all levels, no weighting, strict heterogeneous rate
kd_looms_bcov1111_strict_ht <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_strict_ht/kd-looms_bcov1111_strict_ht.trees")
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
  here("output/trees/kd-looms_bcov1111_strict_ht_consensus.tree")
)

## Looms, binary covarion, all levels, no weighting, relaxed uniform rate
kd_looms_bcov1111_relaxed_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_relaxed_uni/kd-looms_bcov1111_relaxed_uni.trees")
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

## Looms, binary covarion, all levels, no weighting, relaxed heterogenous rate
kd_looms_bcov1111_relaxed_ht <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_relaxed_ht/kd_looms_bcov1111_relaxed_ht.trees")
)
kd_looms_bcov1111_relaxed_ht <- kd_looms_bcov1111_relaxed_ht[ceiling(length(kd_looms_bcov1111_relaxed_ht) * burnin):length(kd_looms_bcov1111_relaxed_ht)]
kd_looms_bcov1111_relaxed_ht_cs <- consensus(kd_looms_bcov1111_relaxed_ht, p = .5, rooted = TRUE)
kd_looms_bcov1111_relaxed_ht_cs <- consensus.edges(
  kd_looms_bcov1111_relaxed_ht,
  consensus.tree = kd_looms_bcov1111_relaxed_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov1111_relaxed_ht_cs)) {
  kd_looms_bcov1111_relaxed_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov1111_relaxed_ht_cs,
  here("output/trees/kd-looms_bcov1111_relaxed_ht_consensus.tree")
)

## Looms, binary covarion, weighted characters, strict uniform rate
kd_looms_bcov8421_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov8421_strict_uni/kd-looms_bcov8421_strict_uni.trees")
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

## Looms, binary covarion, weighted characters, strict heterogenous rate
kd_looms_bcov8421_strict_ht <- read.nexus(
  here("data/kd-looms/kd-looms_bcov8421_strict_ht/kd-looms_bcov8421_strict_ht.trees")
)
kd_looms_bcov8421_strict_ht <- kd_looms_bcov8421_strict_ht[ceiling(length(kd_looms_bcov8421_strict_ht) * burnin):length(kd_looms_bcov8421_strict_ht)]
kd_looms_bcov8421_strict_ht_cs <- consensus(kd_looms_bcov8421_strict_ht, p = .5, rooted = TRUE)
kd_looms_bcov8421_strict_ht_cs <- consensus.edges(kd_looms_bcov8421_strict_ht,
  consensus.tree = kd_looms_bcov8421_strict_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov8421_strict_ht_cs)) {
  kd_looms_bcov8421_strict_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov8421_strict_ht_cs,
  here("output/trees/kd-looms_bcov8421_strict_ht_consensus.tree")
)

## Looms, binary covarion, basic features only, strict uniform rate
kd_looms_bcov_basic_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_basic_strict_uni/kd-looms_bcov_basic_strict_uni.trees")
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

## Looms, binary covarion, basic features only, strict heterogeneous rate
kd_looms_bcov_basic_strict_ht <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_basic_strict_ht/kd-looms_bcov_basic_strict_ht.trees")
)
kd_looms_bcov_basic_strict_ht <- kd_looms_bcov_basic_strict_ht[ceiling(length(kd_looms_bcov_basic_strict_ht) * burnin):length(kd_looms_bcov_basic_strict_ht)]
kd_looms_bcov_basic_strict_ht_cs <- consensus(kd_looms_bcov_basic_strict_ht, p = .5, rooted = TRUE)
kd_looms_bcov_basic_strict_ht_cs <- consensus.edges(kd_looms_bcov_basic_strict_ht,
  consensus.tree = kd_looms_bcov_basic_strict_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov_basic_strict_ht_cs)) {
  kd_looms_bcov_basic_strict_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov_basic_strict_ht_cs,
  here("output/trees/kd-looms_bcov_basic_strict_ht_consensus.tree")
)

## Looms, binary covarion, pattern features only, strict uniform rate
kd_looms_bcov_patterns_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_patterns_strict_uni/kd-looms_bcov_patterns_strict_uni.trees")
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

## Looms, binary covarion, pattern features only, strict heterogeneous rate
kd_looms_bcov_patterns_strict_ht <- read.nexus(
  here("data/kd-looms/kd-looms_bcov_patterns_strict_ht/kd-looms_bcov_patterns_strict_ht.trees")
)
kd_looms_bcov_patterns_strict_ht <- kd_looms_bcov_patterns_strict_ht[ceiling(length(kd_looms_bcov_patterns_strict_ht) * burnin):length(kd_looms_bcov_patterns_strict_ht)]
kd_looms_bcov_patterns_strict_ht_cs <- consensus(kd_looms_bcov_patterns_strict_ht, p = .5, rooted = TRUE)
kd_looms_bcov_patterns_strict_ht_cs <- consensus.edges(kd_looms_bcov_patterns_strict_ht,
  consensus.tree = kd_looms_bcov_patterns_strict_ht_cs,
  rooted = TRUE
)
if (!is.rooted(kd_looms_bcov_patterns_strict_ht_cs)) {
  kd_looms_bcov_patterns_strict_ht_cs$root.edge.length <- 0
}
write.tree(
  kd_looms_bcov_patterns_strict_ht_cs,
  here("output/trees/kd-looms_bcov_patterns_strict_ht_consensus.tree")
)

## Looms, ctmc strict uniform rate
kd_looms_ctmc1111_strict_uni <- read.nexus(
  here("data/kd-looms/kd-looms_ctmc1111_strict_uni/kd-looms_ctmc1111_strict_uni.trees")
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
  here("data/kd-looms/kd-looms_ctmc1111_strict_ht/kd-looms_ctmc1111_strict_ht.trees")
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


# Age and probabilities of clades -----------------------------------------

getMRCA_age <- function(tree, tips) {
  tips <- if (is.character(tips)) which(tree$tip.label %in% tips) else tips
  mrca <- ifelse(length(tips) > 1, getMRCA(tree, tips), tips)
  root_age <- max(node.depth.edgelength(tree))
  root_age - node.depth.edgelength(tree)[mrca]
}

kd_lgs_phylo <- list(
  "bcov_relaxed_ht_pos" = kd_lgs_bcov_relaxed_ht_pos,
  "bcov_relaxed_uni" = kd_lgs_bcov_relaxed_uni,
  "bcov_strict_ht_pos" = kd_lgs_bcov_strict_ht_pos,
  "bcov_strict_uni" = kd_lgs_bcov_strict_uni
)
kam_tai <- kd_lgs_bcov_relaxed_ht_pos[[1]]$tip.label |>
  str_subset("^(Ks|Tc|Tn|Tsw)")
tai_yay <- kd_lgs_bcov_relaxed_ht_pos[[1]]$tip.label |>
  str_subset("^(Tc|Tn|Tsw)")
be_kam_tai <- kd_lgs_bcov_relaxed_ht_pos[[1]]$tip.label |>
  str_subset("^(Be|Ks|Tc|Tn|Tsw)")

kd_lgs_clade_ages <- map_df(1:length(kd_lgs_phylo), function(i) {
  map_df(1:length(kd_lgs_phylo[[i]]), ~ tibble(
    model = names(kd_lgs_phylo[i]),
    KraDai_age = getMRCA_age(kd_lgs_phylo[[i]][[.x]], kd_lgs_phylo[[i]][[.x]]$tip.label),
    KraDai_mono = TRUE,
    BeKamTai_age = getMRCA_age(kd_lgs_phylo[[i]][[.x]], be_kam_tai),
    BeKamTai_mono = is.monophyletic(kd_lgs_phylo[[i]][[.x]], be_kam_tai),
    KamTai_age = getMRCA_age(kd_lgs_phylo[[i]][[.x]], kam_tai),
    KamTai_mono = is.monophyletic(kd_lgs_phylo[[i]][[.x]], kam_tai),
    TaiYay_age = getMRCA_age(kd_lgs_phylo[[i]][[.x]], tai_yay),
    TaiYay_mono = is.monophyletic(kd_lgs_phylo[[i]][[.x]], tai_yay)
  )) |>
    rowid_to_column()
}) |>
  pivot_longer(-c(rowid, model)) |>
  separate_wider_delim(name, "_", names = c("group", "type")) |>
  pivot_wider(names_from = type, values_from = value) |>
  mutate(group = str_replace_all(group, "(?<=[a-z])(?=[A-Z])", "-")) |>
  mutate(mono = as.logical(mono))
write_csv(kd_lgs_clade_ages, here("output/data/kd-lgs_clade_ages.csv"))

kd_lgs_clade_ages_summary <- kd_lgs_clade_ages |>
  group_by(model, group) |>
  summarise(
    mean = mean(age),
    median = median(age),
    sd = sd(age),
    HPDI_lower = hdi(age)["lower"],
    HPDI_upper = hdi(age)["upper"],
    monophyletic = mean(mono)
  ) |>
  ungroup() |>
  mutate(n_lgs = case_when(
    group == "Kra-Dai" ~ length(kd_lgs_phylo[[1]][[1]]$tip.label),
    group == "Be-Kam-Tai" ~ length(be_kam_tai),
    group == "Kam-Tai" ~ length(kam_tai),
    group == "Tai-Yay" ~ length(tai_yay)
  ), .after = group) |>
  arrange(model, -n_lgs)
write_csv(kd_lgs_clade_ages_summary, here("output/data/kd-lgs_clade_ages_summary.csv"))


# Mutation rates ---------------------------------------------------------------
burnin <- .1

# Languages
kd_lgs_concepts <- read_csv(here("data/kd-lgs/kd-lgs_lx.csv")) |>
  count(concept_id)
kd_lgs_pos <- read_csv(here("data/kd-lgs/kd-lgs_lx.csv")) |>
  count(pos)

kd_lgs_mu_pos <- parse_beast_tracelog_file(
  here("data/kd-lgs/kd-lgs_bcov_relaxed_ht_pos/kd-lgs_bcov_relaxed_ht_pos.log")
) |>
  as_tibble() |>
  select(Sample, starts_with("mutationRate.s."))
kd_lgs_mu_pos_tb <- kd_lgs_mu_pos |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  filter(burnin == FALSE) |>
  select(-burnin, -rowid) |>
  pivot_longer(-Sample, names_to = "pos", values_to = "rate") |>
  mutate(pos = str_remove(pos, "mutationRate.s."))
write_csv(kd_lgs_mu_pos_tb, here("output/data/kd-lgs_mu_pos.csv"))

kd_lgs_mu_pos_summary <- kd_lgs_mu_pos_tb |>
  group_by(pos) |>
  summarise(
    mean = mean(rate),
    median = median(rate),
    sd = sd(rate),
    HPDI_lower = hdi(rate)["lower"],
    HPDI_upper = hdi(rate)["upper"]
  ) |>
  left_join(kd_lgs_pos) |>
  relocate(n, .after = pos)
write_csv(kd_lgs_mu_pos_summary, here("output/data/kd-lgs_mu_pos_summary.csv"))

# Looms
kd_looms_characters <- read_csv(here("data/kd-looms/kd-looms_characters.csv")) |>
  select(code, level)

kd_looms_mu_bylevel <- parse_beast_tracelog_file(
  here("data/kd-looms/kd-looms_bcov1111_relaxed_ht/kd-looms_bcov1111_relaxed_ht.log")
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


# Model comparison --------------------------------------------------------

## Extract the marginal likelihood values
out_files <- c(
  list.files(
    c(here("data/kd-lgs/model_choice"), here("data/kd-looms/model_choice")),
    "\\.out",
    recursive = TRUE,
    full.names = TRUE
  ),
  here("data/kd-pruned/kd-merged_pruned_bcov_strict_bd_ns/kd-merged_pruned_bcov_strict_bd_ns.out"),
  here("data/kd-pruned/kd-looms_pruned_bcov_strict_bd_ns/kd-looms_pruned_bcov_strict_bd_ns.out"),
  here("data/kd-pruned/kd-lgs_pruned_bcov_strict_bd_ns/kd-lgs_pruned_bcov_strict_bd_ns.out")
)

out_files |>
  map_df(~
    read_lines(.x) |>
      str_subset("Marginal likelihood") |>
      tail(n = 1) |>
      enframe(name = NULL) |>
      mutate(type = ifelse(str_detect(.x, "pruned"), "pruned", "full")) |>
      mutate(data = str_extract(.x, "(?<=(choice|pruned)/kd-)(lgs|looms|merged)")) |>
      mutate(substitution = ifelse(str_detect(.x, "bcov"), "binary covarion", "CTMC")) |>
      mutate(clock = str_extract(.x, "relaxed|strict")) |>
      mutate(rate = ifelse(str_detect(.x, "ht"), "heterogeneous", "uniform")) |>
      mutate(ML = str_extract(value, "(?<=hood: )-[0-9.]+") |> as.numeric()) |>
      mutate(sd = str_extract(value, "(?<=SD=\\()[0-9.]+") |> as.numeric()) |>
      select(-value)) |>
  arrange(type, data, -ML) |>
  write_csv(here("output/data/models_summary.csv"))

## Prune the trees and unify the tip labels to later compute K

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |>
  select(group, lng_label)
kd_lgs_pruned_tips <- ReadAsPhyDat(here("data/nexus/kd-lgs_pruned.nex")) |>
  as_tibble() |>
  colnames()
kd_looms_pruned_tips <- ReadAsPhyDat(here("data/nexus/kd-looms_pruned.nex")) |>
  as_tibble() |>
  colnames()

kd_lgs_phylo <- read.nexus(here("data/kd-lgs/kd-lgs_bcov_relaxed_ht_pos/kd-lgs_bcov_relaxed_ht_pos.trees"))
kd_lgs_phylo <- kd_lgs_phylo[ceiling(length(kd_lgs_phylo) * burnin):length(kd_lgs_phylo)]
kd_lgs_phylo <- keep.tip(kd_lgs_phylo, kd_lgs_pruned_tips)
write.tree(kd_lgs_phylo, here("output/trees/kd-lgs_prunedk.trees"))

kd_looms_phylo <- read.nexus(
  here("data/kd-looms/kd-looms_bcov1111_relaxed_ht/kd_looms_bcov1111_relaxed_ht.trees")
)
kd_looms_phylo <- kd_looms_phylo[ceiling(length(kd_looms_phylo) * burnin):length(kd_looms_phylo)]
kd_looms_phylo <- map(1:length(kd_looms_phylo), ~ kd_looms_phylo[[.x]] |>
  keep.tip(kd_looms_pruned_tips) |>
  fortify() |>
  mutate(label = str_replace_all(label, "_", " ")) |>
  left_join(kd_looms, by = c("label" = "group")) |>
  mutate(label = lng_label) |>
  as.phylo()) |>
  as.multiPhylo()
write.tree(kd_looms_phylo, here("output/trees/kd-looms_prunedk.trees"))
