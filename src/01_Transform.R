library(here)
library(phangorn)
library(phytools)
library(tracerer)
library(HDInterval)
library(tidyverse)


# Consensus trees -------------------------------------------------------------------------------------------------

kd_lgs_bcov <- read.tree(here("data/kd-lgs/kd-lgs_bcov/kd-lgs_bcov_trimmed.trees"))
kd_lgs_bcov_cs <- consensus(kd_lgs_bcov, p = .5, rooted = TRUE)
kd_lgs_bcov_cs <- consensus.edges(kd_lgs_bcov, consensus.tree = kd_lgs_bcov_cs, rooted = TRUE)
write.tree(kd_lgs_bcov_cs, here("output/trees/kd-lgs_bcov_consensus.tree"))

kd_looms_bcov1000 <- read.tree(here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000_trimmed.trees"))
kd_looms_bcov1000_cs <- consensus(kd_looms_bcov1000, p = .5, rooted = TRUE)
kd_looms_bcov1000_cs <- consensus.edges(kd_looms_bcov1000, consensus.tree = kd_looms_bcov1000_cs, rooted = TRUE)
write.tree(kd_looms_bcov1000_cs, here("output/trees/kd-looms_bcov1000_consensus.tree"))

kd_looms_bcov1111 <- read.tree(here("data/kd-looms/kd-looms_bcov1111/kd-looms_bcov1111_trimmed.trees"))
kd_looms_bcov1111_cs <- consensus(kd_looms_bcov1111, p = .5, rooted = TRUE)
kd_looms_bcov1111_cs <- consensus.edges(kd_looms_bcov1111, consensus.tree = kd_looms_bcov1111_cs, rooted = TRUE)
write.tree(kd_looms_bcov1111_cs, here("output/trees/kd-looms_bcov1111_consensus.tree"))

kd_looms_bcov1000 <- read.tree(here("data/kd-looms/kd-looms_bcov1000/kd-looms_bcov1000_trimmed.trees"))
kd_looms_bcov1000_cs <- consensus(kd_looms_bcov1000, p = .5, rooted = TRUE)
kd_looms_bcov1000_cs <- consensus.edges(kd_looms_bcov1000, consensus.tree = kd_looms_bcov1000_cs, rooted = TRUE)
write.tree(kd_looms_bcov1000_cs, here("output/trees/kd-looms_bcov1000_consensus.tree"))
