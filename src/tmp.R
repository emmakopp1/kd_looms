# x <- kd_lgs_bcov_relaxed_ht_cs_tree
# clade_tips <- x$tip.label
#
# mrcas <- kd_lgs |>
#   distinct(lng_group_code, lng_group, color) |>
#   rowwise() |>
#   mutate(mrca = getMRCA(x, str_subset(clade_tips, paste0("^", lng_group_code)))) |>
#   arrange(mrca)
#
# x <- x |>
#   ggtree(ladderize = FALSE, branch.length = "none")
#
# mrcas |>
#   filter(lng_group_code != "Tc") |>
#   # pull(lng_group_code) |>
#   group_by(lng_group_code) |>
#   group_walk(function(y, z) {
#     x <<- collapse(x, y$mrca, "max", fill = y$color)
#   })

# library(here)
# library(phangorn)
# library(tidyverse)
#
# kd_looms_bcov1111_strict_ht_cs_tree <- read.tree(
#   here("output/trees/kd-looms_bcov1111_strict_ht_consensus.tree")
# )
# |>
#   mutate(model = str_replace(model, "bcov", "binary covarion")) |>
#   mutate(model = str_replace(model, "ht", "heterogeneous")) |>
#   mutate(model = str_replace(model, "uni", "uniform")) |>
#   separate_wider_delim(model, "_", names = c("substitution", "clock", "rate"))

