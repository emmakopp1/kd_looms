library(here)
library(utils)
library(phangorn)
library(phytools)
library(TreeTools)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(tidyverse)
library(stringi)
# BiocManager::install("ggtree")
library(ggtree)
nodeid.tbl_tree <- getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- getFromNamespace("parent.tbl_tree", "tidytree")
library(ggthemes)
library(ggridges)
library(ggtext)
library(ggstar)
library(ggforce)
library(ggnewscale)
library(knitr)
library(kableExtra)
library(patchwork)
library(tree)

dir.create(here("output/figures"))
dir.create(here("output/tables"))

wd <- 14
ht <- 25 # 27.6
lwd <- .5 / .pt
base_font <- "Noto Sans Condensed"
base_font1 <- "Noto Sans SemiCondensed"
base_font2 <- "Noto Sans ExtraCondensed"
base_font_size <- 9
ytheme <- theme(
  legend.text = element_text(
    family = base_font,
    size = base_font_size - 1,
    margin = margin(0, 0, 0, 0, unit = "line")
  ),
  legend.title = element_text(
    family = base_font,
    size = base_font_size,
    margin = margin(0, 0, .5, 0, unit = "line")
  ),
  legend.background = element_rect(fill = NA),
  legend.position = c(0, 1),
  legend.justification = c(0, 1),
  legend.key.spacing.y = unit(0, "line"),
  legend.margin = margin(0, 0, 0, 0, unit = "line"),
  legend.box.margin = margin(0, 0, 0, 0, unit = "line"),
  legend.box.spacing = unit(0, "pt")
)
xtheme <- ytheme + theme(aspect.ratio = 1.15)
theme_set(xtheme)

plt <- tableau_color_pal("Classic 10 Medium")(10)
plt2 <- tableau_color_pal("Classic 10")(10)
plt3 <- tableau_color_pal("Classic 10 Light")(10)

plt_tb <- tribble(
  ~loom_type_code, ~lng_group_code, ~order,
  "FBBS", "Hlai", 1,
  NA, "Kra", 10,
  NA, "Be", 8,
  "BFYRH", "Ks", 9,
  NA, "Tn", 6,
  "BFSRH", "Tc", 3,
  "FCBunique", NA, 2,
  "BFcant", NA, 5,
  "FCBcant", NA, 4,
  "FCB", "Tsw", 7,
) |>
  mutate(color_looms = plt[order], color_lgs = plt2[order])

lgs_order <- c(
  "Hlai",
  "Kra",
  "Be",
  "Ks",
  "Tn",
  "Tc",
  "Tsw"
)
looms_order <- c(
  "FBBS",
  "BFSRH",
  "BFYRH",
  "FCBunique",
  "BFcant",
  "FCBcant",
  "FCB"
)

plt_tb_lgs <- plt_tb |>
  filter(!is.na(lng_group_code)) |>
  rename(color = color_lgs) |>
  mutate(lng_group_code = factor(lng_group_code, levels = lgs_order)) |>
  arrange(lng_group_code) |>
  mutate(color = fct_inorder(color)) |>
  select(lng_group_code, color)

plt_tb_looms <- plt_tb |>
  filter(!is.na(loom_type_code)) |>
  rename(color = color_looms) |>
  mutate(loom_type_code = factor(loom_type_code, levels = looms_order)) |>
  arrange(loom_type_code) |>
  mutate(color = fct_inorder(color)) |>
  select(loom_type_code, color)

kd_looms <- read_csv(here("data/kd-looms/kd-looms_datapoints.csv")) |>
  mutate(loom_type_code = factor(loom_type_code,
    levels = levels(plt_tb_looms$loom_type_code)
  )) |>
  mutate(lng_label = paste0(str_replace_na(lng_group_code, ""), lng)) |>
  left_join(plt_tb_looms) |>
  arrange(loom_type_code, lng_group) |>
  mutate(loom_type = fct_inorder(loom_type))

kd_lgs <- read_csv(here("data/kd-lgs/kd-lgs_datapoints.csv")) |>
  mutate(label = paste0(lng_group_code, lng)) |>
  mutate(lng_group_code = factor(lng_group_code,
    levels = levels(plt_tb_lgs$lng_group_code)
  )) |>
  left_join(plt_tb_lgs) |>
  arrange(lng_group_code, lng) |>
  mutate(lng_group = fct_inorder(lng_group))


# Consensus trees --------------------------------------------------------------

# Function to draw all trees with the same style
cs_tree <- function(tr, fontsize = base_font_size) {
  ggtree(tr, ladderize = FALSE, size = lwd) +
    geom_tiplab(aes(fill = color),
      geom = "label",
      label.size = 0,
      label.padding = unit(.1, "lines"),
      family = base_font, size = fontsize / .pt
    ) +
    geom_nodelab(
      family = base_font,
      size = (base_font_size - 1) / .pt,
      hjust = 1.5,
      vjust = -.5
    ) +
    scale_y_reverse() +
    coord_cartesian(clip = "off", expand = FALSE) +
    scale_fill_identity(
      guide = guide_legend(),
      labels = levels(kd_looms$loom_type)
    ) +
    guides(fill = guide_legend(
      title = "Loom type",
      override.aes = aes(label = "     ")
    )) +
    xtheme +
    theme(
      legend.position = c(0, 0),
      legend.justification = c(0, 0)
    )
}


# Languages

## Languages, binary covarion relaxed heterogeneous, by part of speech
kd_lgs_bcov_relaxed_ht_pos_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_relaxed_ht_pos_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_relaxed_ht_pos_cs_tree)) {
  kd_lgs_bcov_relaxed_ht_pos_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_relaxed_ht_pos_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_relaxed_ht_pos_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_relaxed_ht_pos_cs_tree$node.label[1] <- NA
kd_lgs_bcov_relaxed_ht_pos_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_relaxed_ht_pos_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_relaxed_ht_pos_cs_tree_plot <- kd_lgs_bcov_relaxed_ht_pos_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_relaxed_ht_pos_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
ggsave(here("output/figures/kd-lgs_bcov_relaxed_ht_pos_cs_tree.pdf"),
       kd_lgs_bcov_relaxed_ht_pos_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_relaxed_ht_pos_cs_tree.pdf"))

## Languages, binary covarion strict heterogeneous, by part of speech
kd_lgs_bcov_strict_ht_pos_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_strict_ht_pos_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_strict_ht_pos_cs_tree)) {
  kd_lgs_bcov_strict_ht_pos_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_strict_ht_pos_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_strict_ht_pos_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_strict_ht_pos_cs_tree$node.label[1] <- NA
kd_lgs_bcov_strict_ht_pos_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_strict_ht_pos_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_strict_ht_pos_cs_tree_plot <- kd_lgs_bcov_strict_ht_pos_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_strict_ht_pos_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
kd_lgs_bcov_strict_ht_pos_cs_tree_plot <- flip(kd_lgs_bcov_strict_ht_pos_cs_tree_plot, 107, 127)
ggsave(here("output/figures/kd-lgs_bcov_strict_ht_pos_cs_tree.pdf"),
       kd_lgs_bcov_strict_ht_pos_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_strict_ht_pos_cs_tree.pdf"))

## Languages, binary covarion relaxed heterogeneous
kd_lgs_bcov_relaxed_ht_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_relaxed_ht_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_relaxed_ht_cs_tree)) {
  kd_lgs_bcov_relaxed_ht_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_relaxed_ht_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_relaxed_ht_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_relaxed_ht_cs_tree$node.label[1] <- NA
kd_lgs_bcov_relaxed_ht_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_relaxed_ht_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_relaxed_ht_cs_tree_plot <- kd_lgs_bcov_relaxed_ht_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_relaxed_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
kd_lgs_bcov_relaxed_ht_cs_tree_plot <- flip(kd_lgs_bcov_relaxed_ht_cs_tree_plot, 167, 162) |>
  flip(167, 154) |>
  flip(105, 125)
ggsave(here("output/figures/kd-lgs_bcov_relaxed_ht_cs_tree.pdf"),
  kd_lgs_bcov_relaxed_ht_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_relaxed_ht_cs_tree.pdf"))

## Languages, binary covarion relaxed uniform
kd_lgs_bcov_relaxed_uni_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_relaxed_uni_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_relaxed_uni_cs_tree)) {
  kd_lgs_bcov_relaxed_uni_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_relaxed_uni_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_relaxed_uni_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_relaxed_uni_cs_tree$node.label[1] <- NA
kd_lgs_bcov_relaxed_uni_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_relaxed_uni_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_relaxed_uni_cs_tree_plot <- kd_lgs_bcov_relaxed_uni_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_relaxed_uni_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
ggsave(here("output/figures/kd-lgs_bcov_relaxed_uni_cs_tree.pdf"),
  kd_lgs_bcov_relaxed_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_relaxed_uni_cs_tree.pdf"))

## Languages, binary covarion strict uniform
kd_lgs_bcov_strict_uni_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_strict_uni_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_strict_uni_cs_tree)) {
  kd_lgs_bcov_strict_uni_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_strict_uni_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_strict_uni_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_strict_uni_cs_tree$node.label[1] <- NA
kd_lgs_bcov_strict_uni_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_strict_uni_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_strict_uni_cs_tree_plot <- kd_lgs_bcov_strict_uni_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_strict_uni_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
kd_lgs_bcov_strict_uni_cs_tree_plot <- flip(kd_lgs_bcov_strict_uni_cs_tree_plot, 127, 108)
ggsave(here("output/figures/kd-lgs_bcov_strict_uni_cs_tree.pdf"),
  kd_lgs_bcov_strict_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_strict_uni_cs_tree.pdf"))

## Languages, binary covarion strict heterogeneous
kd_lgs_bcov_strict_ht_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_strict_ht_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_strict_ht_cs_tree)) {
  kd_lgs_bcov_strict_ht_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_strict_ht_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_strict_ht_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_strict_ht_cs_tree$node.label[1] <- NA
kd_lgs_bcov_strict_ht_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_strict_ht_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_strict_ht_cs_tree_plot <- kd_lgs_bcov_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
kd_lgs_bcov_strict_ht_cs_tree_plot <- flip(kd_lgs_bcov_strict_ht_cs_tree_plot, 130, 108)
ggsave(here("output/figures/kd-lgs_bcov_strict_ht_cs_tree.pdf"),
       kd_lgs_bcov_strict_ht_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_strict_ht_cs_tree.pdf"))

## Languages, binary covarion strict heterogeneous with 4 rates
kd_lgs_bcov_strict_ht_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_strict_ht_4rates_consensus.tree"))
if (!is.rooted(kd_lgs_bcov_strict_ht_cs_tree)) {
  kd_lgs_bcov_strict_ht_cs_tree$root.edge.length <- 0
}
kd_lgs_bcov_strict_ht_cs_tree$node.label <- round(
  as.numeric(kd_lgs_bcov_strict_ht_cs_tree$node.label),
  2
) * 100
kd_lgs_bcov_strict_ht_cs_tree$node.label[1] <- NA
kd_lgs_bcov_strict_ht_cs_tree$tip.label <- str_replace_all(
  kd_lgs_bcov_strict_ht_cs_tree$tip.label,
  "_",
  " "
)

kd_lgs_bcov_strict_ht_cs_tree_plot <- kd_lgs_bcov_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree(base_font_size - 1) +
  geom_rootedge(
    max(node.depth.edgelength(kd_lgs_bcov_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  scale_fill_identity(
    guide = guide_legend(),
    labels = levels(kd_lgs$lng_group)
  ) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2.25,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
kd_lgs_bcov_strict_ht_cs_tree_plot <- flip(kd_lgs_bcov_strict_ht_cs_tree_plot, 130, 108)
ggsave(here("output/figures/kd-lgs_bcov_strict_ht_cs_tree.pdf"),
       kd_lgs_bcov_strict_ht_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 3, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_bcov_strict_ht_cs_tree.pdf"))

# Looms

## Looms, binary covarion, level 1 characters only, strict uniform
kd_looms_bcov1000_strict_uni_cs_tree <- read.tree(here("output/trees/kd-looms_bcov1000_strict_uni_consensus.tree"))
if (!is.rooted(kd_looms_bcov1000_strict_uni_cs_tree)) {
  kd_looms_bcov1000_strict_uni_cs_tree$root.edge.length <- 0
}
kd_looms_bcov1000_strict_uni_cs_tree$node.label <- round(as.numeric(kd_looms_bcov1000_strict_uni_cs_tree$node.label), 2) * 100
kd_looms_bcov1000_strict_uni_cs_tree$node.label[1] <- NA
kd_looms_bcov1000_strict_uni_cs_tree$tip.label <- str_replace_all(kd_looms_bcov1000_strict_uni_cs_tree$tip.label, "_", " ")

kd_looms_bcov1000_strict_uni_cs_tree_plot <- kd_looms_bcov1000_strict_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov1000_strict_uni_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.5, 0, 0, unit = "line"))
kd_looms_bcov1000_strict_uni_cs_tree_plot <- rotate(kd_looms_bcov1000_strict_uni_cs_tree_plot, 37)
ggsave(here("output/figures/kd-looms_bcov1000_strict_uni_cs_tree.pdf"),
  kd_looms_bcov1000_strict_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov1000_strict_uni_cs_tree.pdf"))

## Looms, binary covarion, all levels, no weighting, relaxed uniform
kd_looms_bcov1111_relaxed_uni_cs_tree <- read.tree(here("output/trees/kd-looms_bcov1111_relaxed_uni_consensus.tree"))
if (!is.rooted(kd_looms_bcov1111_relaxed_uni_cs_tree)) {
  kd_looms_bcov1111_relaxed_uni_cs_tree$root.edge.length <- 0
}
kd_looms_bcov1111_relaxed_uni_cs_tree$node.label <- round(as.numeric(kd_looms_bcov1111_relaxed_uni_cs_tree$node.label), 2) * 100
kd_looms_bcov1111_relaxed_uni_cs_tree$node.label[1] <- NA
kd_looms_bcov1111_relaxed_uni_cs_tree$tip.label <- str_replace_all(kd_looms_bcov1111_relaxed_uni_cs_tree$tip.label, "_", " ")

kd_looms_bcov1111_relaxed_uni_cs_tree_plot <- kd_looms_bcov1111_relaxed_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(
      node.depth.edgelength(kd_looms_bcov1111_relaxed_uni_cs_tree)
    ) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.75, 0, 0, unit = "line"))
kd_looms_bcov1111_relaxed_uni_cs_tree_plot <- flip(kd_looms_bcov1111_relaxed_uni_cs_tree_plot, 32, 43)
kd_looms_bcov1111_relaxed_uni_cs_tree_plot <- flip(kd_looms_bcov1111_relaxed_uni_cs_tree_plot, 50, 48)
ggsave(here("output/figures/kd-looms_bcov1111_relaxed_uni_cs_tree.pdf"),
  kd_looms_bcov1111_relaxed_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov1111_relaxed_uni_cs_tree.pdf"))

## Looms, binary covarion, all levels, no weighting, strict, 4 variable rates
kd_looms_bcov1111_strict_ht_cs_tree <- read.tree(here("output/trees/kd-looms_bcov1111_strict_ht_consensus.tree"))
if (!is.rooted(kd_looms_bcov1111_strict_ht_cs_tree)) {
  kd_looms_bcov1111_strict_ht_cs_tree$root.edge.length <- 0
}
kd_looms_bcov1111_strict_ht_cs_tree$node.label <- round(
  as.numeric(kd_looms_bcov1111_strict_ht_cs_tree$node.label),
  2
) * 100
kd_looms_bcov1111_strict_ht_cs_tree$node.label[1] <- NA
kd_looms_bcov1111_strict_ht_cs_tree$tip.label <- str_replace_all(
  kd_looms_bcov1111_strict_ht_cs_tree$tip.label,
  "_",
  " "
)

kd_looms_bcov1111_strict_ht_cs_tree_plot <- kd_looms_bcov1111_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov1111_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.95, 0, 0, unit = "line"))
kd_looms_bcov1111_strict_ht_cs_tree_plot <- flip(kd_looms_bcov1111_strict_ht_cs_tree_plot, 32, 43)
ggsave(here("output/figures/kd-looms_bcov1111_strict_ht_cs_tree.pdf"),
  kd_looms_bcov1111_strict_ht_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov1111_strict_ht_cs_tree.pdf"))

## Looms, binary covarion, all levels, no weighting, strict uniform
kd_looms_bcov1111_strict_uni_cs_tree <- read.tree(here("output/trees/kd-looms_bcov1111_strict_uni_consensus.tree"))
if (!is.rooted(kd_looms_bcov1111_strict_uni_cs_tree)) {
  kd_looms_bcov1111_strict_uni_cs_tree$root.edge.length <- 0
}
kd_looms_bcov1111_strict_uni_cs_tree$node.label <- round(as.numeric(kd_looms_bcov1111_strict_uni_cs_tree$node.label), 2) * 100
kd_looms_bcov1111_strict_uni_cs_tree$node.label[1] <- NA
kd_looms_bcov1111_strict_uni_cs_tree$tip.label <- str_replace_all(kd_looms_bcov1111_strict_uni_cs_tree$tip.label, "_", " ")

kd_looms_bcov1111_strict_uni_cs_tree_plot <- kd_looms_bcov1111_strict_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(
      node.depth.edgelength(kd_looms_bcov1111_strict_uni_cs_tree)
    ) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.75, 0, 0, unit = "line"))
kd_looms_bcov1111_strict_uni_cs_tree_plot <- flip(kd_looms_bcov1111_strict_uni_cs_tree_plot, 32, 43)
ggsave(here("output/figures/kd-looms_bcov1111_strict_uni_cs_tree.pdf"),
       kd_looms_bcov1111_strict_uni_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov1111_strict_uni_cs_tree.pdf"))

## Looms, CTMC, all levels, no weighting, strict uniform
kd_looms_ctmc1111_stric_uni_cs_tree <- read.tree(here("output/trees/kd-looms_ctmc1111_strict_uni_consensus.tree"))
if (!is.rooted(kd_looms_ctmc1111_stric_uni_cs_tree)) {
  kd_looms_ctmc1111_stric_uni_cs_tree$root.edge.length <- 0
}
kd_looms_ctmc1111_stric_uni_cs_tree$node.label <- round(as.numeric(kd_looms_ctmc1111_stric_uni_cs_tree$node.label), 2) * 100
kd_looms_ctmc1111_stric_uni_cs_tree$node.label[1] <- NA
kd_looms_ctmc1111_stric_uni_cs_tree$tip.label <- str_replace_all(kd_looms_ctmc1111_stric_uni_cs_tree$tip.label, "_", " ")

kd_looms_ctmc1111_stric_uni_cs_tree_plot <- kd_looms_ctmc1111_stric_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(
      node.depth.edgelength(kd_looms_ctmc1111_stric_uni_cs_tree)
    ) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.75, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_ctmc1111_strict_uni_cs_tree.pdf"),
  kd_looms_ctmc1111_stric_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_ctmc1111_strict_uni_cs_tree.pdf"))

## Looms, binary covarion, strict, 4 variable rates
kd_looms_ctmc1111_strict_ht_cs_tree <- read.tree(here("output/trees/kd-looms_ctmc1111_strict_ht_consensus.tree"))
if (!is.rooted(kd_looms_ctmc1111_strict_ht_cs_tree)) {
  kd_looms_ctmc1111_strict_ht_cs_tree$root.edge.length <- 0
}
kd_looms_ctmc1111_strict_ht_cs_tree$node.label <- round(
  as.numeric(kd_looms_ctmc1111_strict_ht_cs_tree$node.label),
  2
) * 100
kd_looms_ctmc1111_strict_ht_cs_tree$node.label[1] <- NA
kd_looms_ctmc1111_strict_ht_cs_tree$tip.label <- str_replace_all(
  kd_looms_ctmc1111_strict_ht_cs_tree$tip.label,
  "_",
  " "
)

kd_looms_ctmc1111_strict_ht_cs_tree_plot <- kd_looms_ctmc1111_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_ctmc1111_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.95, 0, 0, unit = "line"))
kd_looms_ctmc1111_strict_ht_cs_tree_plot <- flip(kd_looms_ctmc1111_strict_ht_cs_tree_plot, 32, 43)
ggsave(here("output/figures/kd-looms_ctmc1111_strict_ht_cs_tree.pdf"),
  kd_looms_ctmc1111_strict_ht_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_ctmc1111_strict_ht_cs_tree.pdf"))

## Looms, binary covarion, weighted characters, strict uniform
kd_looms_bcov8421_strict_uni_cs_tree <- read.tree(here("output/trees/kd-looms_bcov8421_strict_uni_consensus.tree"))
if (!is.rooted(kd_looms_bcov8421_strict_uni_cs_tree)) {
  kd_looms_bcov8421_strict_uni_cs_tree$root.edge.length <- 0
}
kd_looms_bcov8421_strict_uni_cs_tree$node.label <- round(as.numeric(kd_looms_bcov8421_strict_uni_cs_tree$node.label), 2) * 100
kd_looms_bcov8421_strict_uni_cs_tree$node.label[1] <- NA
kd_looms_bcov8421_strict_uni_cs_tree$tip.label <- str_replace_all(kd_looms_bcov8421_strict_uni_cs_tree$tip.label, "_", " ")

kd_looms_bcov8421_strict_uni_cs_tree_plot <- kd_looms_bcov8421_strict_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov8421_strict_uni_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.95, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_bcov8421_strict_uni_cs_tree.pdf"),
  kd_looms_bcov8421_strict_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov8421_strict_uni_cs_tree.pdf"))

## Looms, binary covarion, weighted characters, strict heterogeneous
kd_looms_bcov8421_strict_ht_cs_tree <- read.tree(here("output/trees/kd-looms_bcov8421_strict_ht_consensus.tree"))
if (!is.rooted(kd_looms_bcov8421_strict_ht_cs_tree)) {
  kd_looms_bcov8421_strict_ht_cs_tree$root.edge.length <- 0
}
kd_looms_bcov8421_strict_ht_cs_tree$node.label <- round(as.numeric(kd_looms_bcov8421_strict_ht_cs_tree$node.label), 2) * 100
kd_looms_bcov8421_strict_ht_cs_tree$node.label[1] <- NA
kd_looms_bcov8421_strict_ht_cs_tree$tip.label <- str_replace_all(kd_looms_bcov8421_strict_ht_cs_tree$tip.label, "_", " ")

kd_looms_bcov8421_strict_ht_cs_tree_plot <- kd_looms_bcov8421_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov8421_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.95, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_bcov8421_strict_ht_cs_tree.pdf"),
       kd_looms_bcov8421_strict_ht_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov8421_strict_ht_cs_tree.pdf"))

## Looms, binary covarion, basic features only, strict uni
kd_looms_bcov_basic_strict_uni_cs_tree <- read.tree(here("output/trees/kd-looms_bcov_basic_strict_uni_consensus.tree"))
if (!is.rooted(kd_looms_bcov_basic_strict_uni_cs_tree)) {
  kd_looms_bcov_basic_strict_uni_cs_tree$root.edge.length <- 0
}
kd_looms_bcov_basic_strict_uni_cs_tree$node.label <- round(as.numeric(kd_looms_bcov_basic_strict_uni_cs_tree$node.label), 2) * 100
kd_looms_bcov_basic_strict_uni_cs_tree$node.label[1] <- NA
kd_looms_bcov_basic_strict_uni_cs_tree$tip.label <- str_replace_all(kd_looms_bcov_basic_strict_uni_cs_tree$tip.label, "_", " ")

kd_looms_bcov_basic_strict_uni_cs_tree_plot <- kd_looms_bcov_basic_strict_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov_basic_strict_uni_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.85, 0, 0, unit = "line"))
kd_looms_bcov_basic_strict_uni_cs_tree_plot <- flip(kd_looms_bcov_basic_strict_uni_cs_tree_plot, 32, 41)
ggsave(here("output/figures/kd-looms_bcov_basic_strict_uni_cs_tree.pdf"),
  kd_looms_bcov_basic_strict_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov_basic_strict_uni_cs_tree.pdf"))

## Looms, binary covarion, basic features only, strict heterogeneous
kd_looms_bcov_basic_strict_ht_cs_tree <- read.tree(here("output/trees/kd-looms_bcov_basic_strict_ht_consensus.tree"))
if (!is.rooted(kd_looms_bcov_basic_strict_ht_cs_tree)) {
  kd_looms_bcov_basic_strict_ht_cs_tree$root.edge.length <- 0
}
kd_looms_bcov_basic_strict_ht_cs_tree$node.label <- round(as.numeric(kd_looms_bcov_basic_strict_ht_cs_tree$node.label), 2) * 100
kd_looms_bcov_basic_strict_ht_cs_tree$node.label[1] <- NA
kd_looms_bcov_basic_strict_ht_cs_tree$tip.label <- str_replace_all(kd_looms_bcov_basic_strict_ht_cs_tree$tip.label, "_", " ")

kd_looms_bcov_basic_strict_ht_cs_tree_plot <- kd_looms_bcov_basic_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov_basic_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(plot.margin = margin(0, 3.85, 0, 0, unit = "line"))
kd_looms_bcov_basic_strict_ht_cs_tree_plot <- flip(kd_looms_bcov_basic_strict_ht_cs_tree_plot, 32, 40)
ggsave(here("output/figures/kd-looms_bcov_basic_strict_ht_cs_tree.pdf"),
       kd_looms_bcov_basic_strict_ht_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov_basic_strict_ht_cs_tree.pdf"))

## Looms, binary covarion, pattern features only, strict uniform
kd_looms_bcov_patterns_strict_uni_cs_tree <- read.tree(here("output/trees/kd-looms_bcov_patterns_strict_uni_consensus.tree"))
if (!is.rooted(kd_looms_bcov_patterns_strict_uni_cs_tree)) {
  kd_looms_bcov_patterns_strict_uni_cs_tree$root.edge.length <- 0
}
kd_looms_bcov_patterns_strict_uni_cs_tree$node.label <- round(as.numeric(kd_looms_bcov_patterns_strict_uni_cs_tree$node.label), 2) * 100
kd_looms_bcov_patterns_strict_uni_cs_tree$node.label[1] <- NA
kd_looms_bcov_patterns_strict_uni_cs_tree$tip.label <- str_replace_all(kd_looms_bcov_patterns_strict_uni_cs_tree$tip.label, "_", " ")

kd_looms_bcov_patterns_strict_uni_cs_tree_plot <- kd_looms_bcov_patterns_strict_uni_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov_patterns_strict_uni_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(
    plot.margin = margin(0, 3, 0, 0, unit = "line"),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(-.15, "line"),
    legend.position = c(0, 1),
    legend.justification = c(0, 0.2)
  )
kd_looms_bcov_patterns_strict_uni_cs_tree_plot <- rotate(kd_looms_bcov_patterns_strict_uni_cs_tree_plot, 30)
kd_looms_bcov_patterns_strict_uni_cs_tree_plot <- rotate(kd_looms_bcov_patterns_strict_uni_cs_tree_plot, 31)
ggsave(here("output/figures/kd-looms_bcov_patterns_strict_uni_cs_tree.pdf"),
  kd_looms_bcov_patterns_strict_uni_cs_tree_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov_patterns_strict_uni_cs_tree.pdf"))

## Looms, binary covarion, pattern features only, strict heterogeneous
kd_looms_bcov_patterns_strict_ht_cs_tree <- read.tree(here("output/trees/kd-looms_bcov_patterns_strict_ht_consensus.tree"))
if (!is.rooted(kd_looms_bcov_patterns_strict_ht_cs_tree)) {
  kd_looms_bcov_patterns_strict_ht_cs_tree$root.edge.length <- 0
}
kd_looms_bcov_patterns_strict_ht_cs_tree$node.label <- round(as.numeric(kd_looms_bcov_patterns_strict_ht_cs_tree$node.label), 2) * 100
kd_looms_bcov_patterns_strict_ht_cs_tree$node.label[1] <- NA
kd_looms_bcov_patterns_strict_ht_cs_tree$tip.label <- str_replace_all(kd_looms_bcov_patterns_strict_ht_cs_tree$tip.label, "_", " ")

kd_looms_bcov_patterns_strict_ht_cs_tree_plot <- kd_looms_bcov_patterns_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng_label), label, lng_label)) |>
  cs_tree() +
  geom_rootedge(
    max(node.depth.edgelength(kd_looms_bcov_patterns_strict_ht_cs_tree)) * .025,
    linewidth = lwd
  ) +
  theme(
    plot.margin = margin(0, 3, 0, 0, unit = "line"),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(-.15, "line"),
    legend.position = c(0, 1),
    legend.justification = c(0, 0.2)
  )
kd_looms_bcov_patterns_strict_ht_cs_tree_plot <- flip(kd_looms_bcov_patterns_strict_ht_cs_tree_plot, 41, 42)
kd_looms_bcov_patterns_strict_ht_cs_tree_plot <- rotate(kd_looms_bcov_patterns_strict_ht_cs_tree_plot, 31)
ggsave(here("output/figures/kd-looms_bcov_patterns_strict_ht_cs_tree.pdf"),
       kd_looms_bcov_patterns_strict_ht_cs_tree_plot,
       device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_bcov_patterns_strict_ht_cs_tree.pdf"))


# Age density distribution for languages ---------------------------------------

kd_lgs_ages <- read_csv(here("output/data/kd-lgs_clade_ages.csv"))

kd_lgs_ages_plot <- kd_lgs_ages |>
  filter(model == "bcov_relaxed_ht_pos") |>
  filter(group != "Be-Kam-Tai") |>
  mutate(group = fct(group, levels = c("Kra-Dai", "Kam-Tai", "Tai-Yay"))) |>
  ggplot(aes(x = age * 1000, y = group)) +
  stat_density_ridges(
    aes(fill = .5 - abs(.5 - after_stat(ecdf))),
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE, scale = 1,
    panel_scaling = FALSE,
    color = "gray50", linewidth = lwd
  ) +
  stat_summary(
    geom = "text", fun = "median",
    aes(label = after_stat(round(x))),
    family = base_font,
    size = base_font_size / .pt, vjust = 1.5
  ) +
  scale_fill_distiller(
    palette = "PuBu",
    direction = 1,
    limits = c(0, .5),
    name = "Tail\nprobability"
  ) +
  scale_x_reverse(limits = c(12500, 0)) +
  xlab("Age (years BP)") +
  ylab("Language group") +
  theme_minimal(base_size = base_font_size, base_family = base_font) +
  xtheme +
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "line"),
    aspect.ratio = 0.618,
    legend.position = "right",
    legend.background = element_rect(color = NA)
  )
ggsave(here("output/figures/kd-lgs_ages_plot.pdf"),
  kd_lgs_ages_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_ages_plot.pdf"))

kd_lgs_ages_summary <- read_csv(here("output/data/kd-lgs_clade_ages_summary.csv"))

kd_lgs_ages_summary |>
  filter(model == "bcov_relaxed_ht_pos") |>
  filter(group != "Be-Kam-Tai") |>
  select(-model) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  unite(hdi, HPDI_lower, HPDI_upper, sep = ", ") |>
  mutate(hdi = paste0("[", hdi, "]")) |>
  rename(languages = n_lgs, `95% HPDI` = hdi) |>
  mutate(monophyletic = paste0(monophyletic * 100, "%")) |>
  kbl(
    digits = 2,
    format = "latex", booktabs = TRUE,
    align = c("l", rep("r", 6))
  ) |>
  write_lines(here("output/tables/kd-lgs_ages_summary.tex"))

# Cophylogeny ------------------------------------------------------------------

kd_lgs_cs <- kd_lgs_bcov_relaxed_ht_pos_cs_tree |> 
  # kd_lgs_bcov_relaxed_uni_cs_tree |>
  fortify() |>
  left_join(kd_lgs) |>
  mutate(label = lng) |>
  as.phylo()
if (!is.rooted(kd_lgs_cs)) {
  kd_lgs_cs$root.edge.length <- 0
}

kd_loom_pb <- c(
  "Dai Huayao",
  "Dai Yuxi Yuanjiang",
  "Dai Jinghong",
  "Zhuang Napo",
  "Zhuang Longzhou",
  "Nung An",
  "Tai Phake"
)

kd_looms_cs <- kd_looms_bcov1111_strict_ht_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(label %in% kd_loom_pb, label, lng)) |>
  as.phylo()
if (!is.rooted(kd_looms_cs)) {
  kd_looms_cs$root.edge.length <- 0
}

kd_cophylo <- cophylo(kd_lgs_cs, kd_looms_cs,
  methods = c("pre", "post"), rotate.multi = FALSE
)
kd_cophylo$trees[[2]]$tip.label <- stri_replace_all_fixed(
  str = kd_cophylo$trees[[2]]$tip.label, pattern = kd_looms$group,
  replacement = kd_looms$lng,
  vectorise_all = FALSE
)
kd_cophylo$trees[[2]] <- kd_cophylo$trees[[2]] |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == lng)) |>
  mutate(label = group) |>
  as.phylo()

kd_lng_tree <- ggtree(kd_cophylo$trees[[1]],
  ladderize = FALSE,
  size = lwd,
  branch.length = "none"
) |>
  # rotate(171) |> 
  rotate(getMRCA(kd_lgs_cs, c("Libo", "Qianxi")))
# rotate(getMRCA(kd_lgs_cs, c("KamLP", "KamRJ")))
# rotate(105) |>
  # rotate(176) |>
  # rotate(143) |>
  # rotate(119) |>
  # rotate(132) |>
  # rotate(155) |>
  # flip(105, 113)
# |>
#   flip(167, 154) |>
#   flip(105, 125)
# flip(167, 162) |>
# flip(167, 154) |>
# flip(162, 154)
#   flip(105, 125)
# flip(113, 128) |>
# flip(103,172)
kd_loom_tree <- ggtree(kd_cophylo$trees[[2]],
  ladderize = FALSE,
  size = lwd,
  branch.length = "none"
) |>
  rotate(getMRCA(kd_looms_cs, c("Chiangmai", "Korat")))
  # rotate(43) |>
  # rotate(47)

kd_lng_tree_data <- kd_lng_tree$data |>
  mutate(lng = label) |>
  left_join(select(kd_lgs, -label)) |>
  left_join(select(kd_looms, group, lng))

kd_loom_tree_data <- kd_loom_tree$data |>
  mutate(group = label) |>
  left_join(kd_looms) |>
  arrange(loom_type)

ry <- kd_lng_tree_data$y
kd_loom_tree_data$x <- ((kd_loom_tree_data$x - min(kd_loom_tree_data$x)) / (max(kd_loom_tree_data$x) - min(kd_loom_tree_data$x))) * (max(kd_lng_tree_data$x) - min(kd_lng_tree_data$x)) + min(kd_lng_tree_data$x)
kd_loom_tree_data$x <- max(kd_loom_tree_data$x) - kd_loom_tree_data$x + max(kd_lng_tree_data$x)
kd_loom_tree_data$x <- kd_loom_tree_data$x + (max(c(kd_lng_tree_data$x, kd_loom_tree_data$x)) - min(c(kd_lng_tree_data$x, kd_loom_tree_data$x))) / 100 * 70
kd_loom_tree_data$y <- ((kd_loom_tree_data$y - min(kd_loom_tree_data$y)) / (max(kd_loom_tree_data$y) - min(kd_loom_tree_data$y))) * (max(ry) - min(ry)) + min(ry)

kd_lng_loom_tree <- kd_lng_tree +
  geom_tree(data = kd_loom_tree_data, linewidth = lwd)
kd_lng_loom_tree_data <- bind_rows(kd_lng_tree_data, kd_loom_tree_data) |>
  filter(!is.na(group) & !is.na(lng)) |>
  mutate(pb = group %in% kd_loom_pb)

x1 <- max(filter(kd_lng_tree_data, isTip == TRUE)$x) + diff(range(
  filter(kd_loom_tree_data, isTip == TRUE)$x,
  filter(kd_lng_tree_data, isTip == TRUE)$x
)) / 2.6
x2 <- min(filter(kd_loom_tree_data, isTip == TRUE)$x) - diff(range(
  filter(kd_loom_tree_data, isTip == TRUE)$x,
  filter(kd_lng_tree_data, isTip == TRUE)$x
)) / 3.2

kd_lgs_looms_links <- kd_lng_loom_tree_data |>
  add_count(lng) |>
  filter(n == 2) |>
  group_by(lng) |>
  mutate(x = ifelse(x == min(x), x1 + .25, x2 - .25))

imgs <- paste0(
  "<img src='", here("data/images/"),
  levels(kd_loom_tree_data$loom_type_code),
  ".png' height=60>"
)

kd_cophylo_plot <- kd_lng_loom_tree +
  geom_rootedge(1, linewidth = lwd) +
  geom_segment(
    data = filter(kd_loom_tree_data, parent == node),
    aes(x = x, xend = x + 1, y = y, yend = y),
    linewidth = lwd
  ) +
  geom_line(
    data = kd_lgs_looms_links,
    aes(x, y, group = lng, linetype = pb, color = pb)
  ) +
  scale_color_manual(guide = "none", values = c("grey50", "grey70")) +
  new_scale_colour() +
  geom_tippoint(data = kd_lng_tree_data, aes(color = color, fill = color)) +
  geom_tiplab(
    data = kd_lng_tree_data,
    aes(label = label, color = color, fill = color, x = x1),
    size = (base_font_size - 2) / .pt,
    family = base_font2, hjust = 1,
    show.legend = FALSE
  ) +
  scale_color_identity(
    guide = guide_legend(
      order = 1,
      position = "left",
      override.aes = list(size = 4),
      theme = theme(
        legend.key.spacing.y = unit(.35, "line"),
        legend.margin = margin(0, -1.25, 0, 0, unit = "line")
      )
    ),
    name = "Language group",
    labels = levels(kd_lng_tree_data$lng_group)
  ) +
  xtheme +
  scale_fill_identity(guide = "none") +
  new_scale_colour() +
  new_scale_fill() +
  geom_tippoint(data = kd_loom_tree_data, aes(color = color, fill = color)) +
  geom_tiplab(
    data = filter(kd_loom_tree_data, isTip == TRUE),
    aes(
      label = str_replace_all(lng, "(.+?)(?=[A-Z])", "\\1 ") |> str_wrap(10),
      color = color,
      x = x2
    ),
    hjust = 0,
    size = (base_font_size - 2) / .pt,
    family = base_font2,
    lineheight = 1,
    show.legend = FALSE
  ) +
  scale_color_identity(
    guide = guide_legend(
      order = 2,
      position = "right",
      label.vjust = 0,
      override.aes = list(size = 4),
      theme = theme(
        legend.key.spacing.y = unit(1, "line"),
        legend.text = element_markdown(),
        legend.margin = margin(0, 0, 0, -.5, unit = "line")
      )
    ),
    name = "Loom type",
    labels = paste0(
      str_wrap(levels(kd_loom_tree_data$loom_type), 20),
      "\n",
      imgs
    ) |>
      str_replace_all("\\n", "<br/>")
  ) +
  scale_fill_identity(
    guide = "none"
  ) +
  scale_linetype(guide = "none") +
  theme(
    aspect.ratio = 3,
    legend.box = "horizontal",
    legend.text = element_text(size = base_font_size - 2, family = base_font2),
    legend.key = element_rect(),
    legend.box.margin = margin(0, 0, 0, 0, unit = "line"),
    legend.background = element_blank()
  )
ggsave(here("output/figures/kd_cophylo_plot.pdf"),
  kd_cophylo_plot,
  device = cairo_pdf, width = wd, height = ht * 2, units = "cm"
)
plot_crop(here("output/figures/kd_cophylo_plot.pdf"))

kd_lgs_pruned_tips <- ReadAsPhyDat(here("data/nexus/kd-lgs_pruned.nex")) |>
  as_tibble() |>
  colnames() |>
  str_remove("^[A-Z][a-z]+(?=[A-Z])")

kd_looms_cs <- kd_looms_bcov1111_strict_ht_cs_tree |>
  fortify() |>
  mutate(label = str_replace_all(label, "_", " ")) |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(is.na(lng), label, lng)) |>
  # mutate(label = ifelse(label %in% kd_loom_pb, label, lng)) |>
  as.phylo()
if (!is.rooted(kd_looms_cs)) {
  kd_looms_cs$root.edge.length <- 0
}

kd_cophylo_pruned <- cophylo(
  keep.tip(kd_lgs_cs, kd_lgs_pruned_tips),
  keep.tip(kd_looms_cs, kd_lgs_pruned_tips),
  methods = c("pre", "post"), rotate.multi = FALSE
)

kd_lng_tree_pruned <- ggtree(kd_cophylo_pruned$trees[[1]],
  ladderize = FALSE,
  size = lwd,
  branch.length = "none"
)
kd_lng_tree_pruned <- flip(kd_lng_tree_pruned, 22, 24) |>
  rotate(23) |>
  rotate(30)
kd_loom_tree_pruned <- ggtree(kd_cophylo_pruned$trees[[2]],
  ladderize = FALSE,
  size = lwd,
  branch.length = "none"
)
kd_loom_tree_pruned <- flip(kd_loom_tree_pruned, 8, 27) |>
  rotate(24) |>
  # flip(1, 22) |>
  # rotate(25) |>
  rotate(27) |>
  rotate(28) |>
  rotate(21) |>
  rotate(22)

kd_lng_tree_pruned_data <- kd_lng_tree_pruned$data |>
  mutate(lng = label) |>
  left_join(select(kd_lgs, -label)) |>
  left_join(select(kd_looms, group, lng))

kd_loom_tree_pruned_data <- kd_loom_tree_pruned$data |>
  mutate(lng = label) |>
  left_join(kd_looms) |>
  arrange(loom_type)

kd_loom_tree_pruned_data$x <- ((kd_loom_tree_pruned_data$x - min(kd_loom_tree_pruned_data$x)) / (max(kd_loom_tree_pruned_data$x) - min(kd_loom_tree_pruned_data$x))) * (max(kd_lng_tree_pruned_data$x) - min(kd_lng_tree_pruned_data$x)) + min(kd_lng_tree_pruned_data$x)
kd_loom_tree_pruned_data$x <- max(kd_loom_tree_pruned_data$x) - kd_loom_tree_pruned_data$x + max(kd_lng_tree_pruned_data$x)
kd_loom_tree_pruned_data$x <- kd_loom_tree_pruned_data$x + (max(c(kd_lng_tree_pruned_data$x, kd_loom_tree_pruned_data$x)) - min(c(kd_lng_tree_pruned_data$x, kd_loom_tree_pruned_data$x))) / 100 * 70

kd_lng_loom_tree_pruned <- kd_lng_tree_pruned +
  geom_tree(data = kd_loom_tree_pruned_data, linewidth = lwd)
kd_lng_loom_tree_pruned_data <- bind_rows(kd_lng_tree_pruned_data, kd_loom_tree_pruned_data) |>
  filter(!is.na(group) & !is.na(lng)) |>
  mutate(pb = group %in% kd_loom_pb)

x1 <- max(filter(kd_lng_tree_pruned_data, isTip == TRUE)$x) + diff(range(
  filter(kd_loom_tree_pruned_data, isTip == TRUE)$x,
  filter(kd_lng_tree_pruned_data, isTip == TRUE)$x
)) / 2.8
x2 <- min(filter(kd_loom_tree_pruned_data, isTip == TRUE)$x) - diff(range(
  filter(kd_loom_tree_pruned_data, isTip == TRUE)$x,
  filter(kd_lng_tree_pruned_data, isTip == TRUE)$x
)) / 2.8

kd_lgs_looms_pruned_links <- kd_lng_loom_tree_pruned_data |>
  add_count(lng) |>
  filter(n == 2) |>
  group_by(lng) |>
  mutate(x = ifelse(x == min(x), x1 + .25, x2 - .25))

imgs <- paste0(
  "<img src='", here("data/images/"),
  levels(kd_loom_tree_pruned_data$loom_type_code),
  ".png' height=30>"
)

kd_cophylo_pruned_plot <- kd_lng_loom_tree_pruned +
  geom_rootedge(1, linewidth = lwd) +
  geom_segment(
    data = filter(kd_loom_tree_pruned_data, parent == node),
    aes(x = x, xend = x + 1, y = y, yend = y),
    linewidth = lwd
  ) +
  geom_line(
    data = kd_lgs_looms_pruned_links,
    aes(x, y, group = lng, linetype = pb, color = pb)
  ) +
  scale_color_manual(guide = "none", values = c("grey50", "grey70")) +
  new_scale_colour() +
  geom_tippoint(
    data = kd_lng_tree_pruned_data,
    aes(color = color, fill = color),
    show.legend = TRUE
  ) +
  geom_tiplab(
    data = kd_lng_tree_pruned_data,
    aes(
      label = str_replace_all(lng, "(?<=[a-z])(?=[A-Z])", " ") |> str_wrap(10),
      color = color,
      x = x1
    ),
    size = (base_font_size - 1) / .pt,
    family = base_font2, hjust = 1,
    show.legend = FALSE
  ) +
  scale_color_identity(
    guide = guide_legend(
      order = 1,
      position = "left",
      override.aes = list(size = 4),
      theme = theme(
        legend.key.spacing.y = unit(.35, "line"),
        legend.margin = margin(0, -1.25, 0, 0, unit = "line")
      )
    ),
    name = "Language group",
    drop = FALSE,
    breaks = levels(kd_lng_tree_pruned_data$color),
    labels = levels(kd_lng_tree_pruned_data$lng_group)
  ) +
  xtheme +
  scale_fill_identity(guide = "none") +
  new_scale_colour() +
  new_scale_fill() +
  geom_tippoint(
    data = kd_loom_tree_pruned_data,
    aes(color = color, fill = color)
  ) +
  geom_tiplab(
    data = filter(kd_loom_tree_pruned_data, isTip == TRUE),
    aes(
      label = str_replace_all(lng, "(?<=[a-z])(?=[A-Z])", " ") |> str_wrap(10),
      color = color,
      x = x2
    ),
    hjust = 0,
    size = (base_font_size - 1) / .pt,
    family = base_font2,
    lineheight = 1,
    show.legend = FALSE
  ) +
  scale_color_identity(
    guide = guide_legend(
      order = 2,
      position = "right",
      label.vjust = 0,
      override.aes = list(size = 4),
      theme = theme(
        legend.key.spacing.y = unit(1, "line"),
        legend.text = element_markdown(),
        legend.margin = margin(0, 0, 0, -.5, unit = "line")
      )
    ),
    name = "Loom type",
    labels = paste0(
      str_wrap(levels(droplevels(kd_loom_tree_pruned_data$loom_type)), 20),
      "\n",
      imgs
    ) |>
      str_replace_all("\\n", "<br/>")
  ) +
  scale_fill_identity(
    guide = "none"
  ) +
  scale_linetype(guide = "none") +
  theme(
    aspect.ratio = 1.5,
    legend.box = "horizontal",
    legend.text = element_text(size = base_font_size - 2, family = base_font2),
    legend.key = element_rect(),
    legend.box.margin = margin(0, 0, 0, 0, unit = "line"),
    legend.background = element_blank()
  )
ggsave(here("output/figures/kd_cophylo_pruned_plot.pdf"),
  kd_cophylo_pruned_plot,
  device = cairo_pdf, width = wd, height = ht * 2, units = "cm"
)
plot_crop(here("output/figures/kd_cophylo_pruned_plot.pdf"))


# Mutation rates ---------------------------------------------------------------

kd_looms_mu_bylevel <- read_csv(here("output/data/kd-looms_mu_bylevel.csv"))
kd_looms_mu_plot <- kd_looms_mu_bylevel |>
  ggplot(aes(y = factor(level), x = rate)) +
  stat_density_ridges(
    aes(fill = .5 - abs(.5 - after_stat(ecdf))),
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    scale = 1,
    panel_scaling = FALSE,
    color = "gray50",
    linewidth = lwd
  ) +
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = after_stat(round(x, 2))),
    family = base_font,
    size = base_font_size / .pt,
    vjust = 1.5
  ) +
  ylab("Level") +
  xlab("Mutation rate") +
  xlim(0, 2) +
  scale_fill_distiller(
    palette = "PuBu",
    direction = 1,
    limits = c(0, .5),
    name = "Tail\nprobability"
  ) +
  theme_minimal(base_size = base_font_size, base_family = base_font) +
  xtheme +
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "line"),
    aspect.ratio = 0.618,
    legend.position = "right",
    legend.background = element_rect(color = NA)
  )
ggsave(here("output/figures/kd-looms_mu_plot.pdf"),
  kd_looms_mu_plot,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_mu_plot.pdf"))

kd_looms_mu_summary <- read_csv(here("output/data/kd-looms_mu_summary.csv"))
kd_looms_mu_summary |>
  mutate(across(everything(), ~ round(.x, 2))) |>
  unite(hdi, HPDI_lower, HPDI_upper, sep = ", ") |>
  mutate(hdi = paste0("[", hdi, "]")) |>
  rename(characters = n_chars, `95% HPDI` = hdi) |>
  kbl(
    digits = 2,
    format = "latex", booktabs = TRUE
  ) |>
  write_lines(here("output/tables/kd-looms_mu_summary.tex"))


# Maps -------------------------------------------------------------------------

theme_set(theme_bw() + ytheme + theme(legend.key.spacing.y = unit(.25, "line")))

prj <- "+proj=laea +lon_0=107.11 +lat_0=19.03 +datum=WGS84 +units=m +no_defs"

kd_lgs_pts <- kd_lgs |>
  filter(!is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(prj)

kd_looms_pts <- kd_looms |>
  filter(!is.na(lon)) |>
  arrange(loom_type) |>
  mutate(loom_type = str_replace_all(loom_type, ", ", ",\n")) |>
  mutate(loom_type = fct_inorder(loom_type)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(prj)

kd_bbx_lat <- bind_rows(kd_lgs_pts, kd_looms_pts) |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(50 * 10^3) |>
  st_bbox()

kd_bbx_lon <- bind_rows(kd_lgs_pts, kd_looms_pts) |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(400 * 10^3) |>
  st_bbox()

kd_bbx <- kd_bbx_lon
kd_bbx["ymin"] <- kd_bbx_lat["ymin"]
kd_bbx["ymax"] <- kd_bbx_lat["ymax"]

kd_bbx_sf <- st_transform(st_as_sfc(kd_bbx), prj)

asia <- ne_countries(scale = "medium") |>
  filter(continent %in% c("Asia", "Oceania")) |>
  st_transform(prj) |>
  st_crop(kd_bbx_sf)

country_lbs <- asia |>
  select(label_x, label_y, name_en) |>
  filter(!(name_en %in% c(
    "Hong Kong",
    "Macau",
    "Malaysia",
    "Bhutan",
    "Bangladesh"
  ))) |>
  mutate(name_en = str_remove(name_en, ".+ ")) |>
  mutate(label_y = ifelse(name_en == "China", 28, label_y)) |>
  mutate(label_y = ifelse(name_en == "India", 26.5, label_y)) |>
  mutate(label_y = ifelse(name_en == "Cambodia", 13, label_y)) |>
  mutate(label_x = ifelse(name_en == "India", 92, label_x)) |>
  mutate(label_y = ifelse(name_en == "Philippines", 17, label_y)) |>
  mutate(label_y = ifelse(name_en == "Laos", 19, label_y)) |>
  mutate(label_y = ifelse(name_en == "Thailand", 16, label_y)) |>
  mutate(label_y = ifelse(name_en == "Vietnam", 21.3, label_y)) |>
  mutate(hjust = ifelse(name_en == "India", 0, .5)) |>
  mutate(hjust = ifelse(name_en == "Cambodia", .4, hjust)) |>
  mutate(hjust = ifelse(name_en == "Bangladesh", 0, hjust)) |>
  mutate(hjust = ifelse(name_en == "Taiwan", .6, hjust)) |>
  mutate(hjust = ifelse(name_en == "Philippines", 1.05, hjust)) |>
  st_drop_geometry() |>
  st_as_sf(coords = c("label_x", "label_y"), crs = 4326)

asia_u <- asia |>
  st_union()

bg_map <- asia |>
  ggplot() +
  geom_sf(fill = "white", linetype = "dashed", linewidth = lwd) +
  geom_sf(data = asia_u, fill = NA, linewidth = lwd) +
  geom_sf_text(
    data = country_lbs, aes(label = name_en, hjust = hjust),
    family = base_font,
    size = base_font_size / .pt,
    color = "grey35"
  ) +
  annotation_scale(
    location = "br",
    height = unit(.25, "line"),
    text_family = base_font,
    text_cex = .7
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(3.25, "line"),
    pad_y = unit(1.25, "line"),
    height = unit(1, "line"),
    width = unit(1, "line"),
    style = north_arrow_orienteering(
      text_family = base_font,
      text_size = base_font_size - 2
    )
  ) +
  xlab(NULL) +
  ylab(NULL) +
  coord_sf(crs = prj, expand = FALSE) +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, .5, unit = "line"),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = base_font_size),
    axis.text = element_text(family = base_font, size = (base_font_size - 1)),
    panel.background = element_rect(fill = "grey85")
  )

kd_lgs_map <- bg_map +
  geom_star(
    data = kd_lgs_pts |> df_spatial(),
    aes(x = x, y = y, color = color, starshape = color, fill = color),
    alpha = .85,
    size = 1
  ) +
  scale_starshape_discrete(
    name = "Language group",
    guide = guide_legend(override.aes = list(
      size = 4,
      fill = levels(kd_lgs_pts$color),
      color = levels(kd_lgs_pts$color)
    )),
    labels = levels(kd_lgs_pts$lng_group)
  ) +
  scale_color_identity(name = "Language group") +
  scale_fill_identity(name = "Language group") +
  coord_sf(crs = prj, expand = FALSE)

kd_looms_map <- bg_map +
  geom_star(
    data = kd_looms_pts |> df_spatial(),
    aes(x = x, y = y, color = color, starshape = color, fill = color),
    alpha = .85,
    size = 1
  ) +
  scale_starshape_discrete(
    name = "Loom type",
    guide = guide_legend(override.aes = list(
      size = 4,
      fill = levels(kd_looms_pts$color),
      color = levels(kd_looms_pts$color)
    )),
    labels = levels(kd_looms_pts$loom_type)
  ) +
  scale_color_identity(name = "Loom type") +
  scale_fill_identity(name = "Loom type") +
  coord_sf(crs = prj, expand = FALSE)
kd_looms_map
kd_lgs_map <- set_dim(kd_lgs_map, get_dim(kd_looms_map))

ggsave(here("output/figures/kd-lgs_map.pdf"),
  kd_lgs_map,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-lgs_map.pdf"))

ggsave(here("output/figures/kd-looms_map.pdf"),
  kd_looms_map,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_map.pdf"))

Tsw_imgs <- paste0(
  "<img src='", here("data/images/"),
  unique(filter(kd_looms_pts, lng_group_code == "Tsw")$loom_type_code),
  ".png' height=60>"
)
Tsw_colors <- unique(filter(kd_looms_pts, lng_group_code == "Tsw")$color)
Tsw_shapes <- starshape_pal()(
  length(levels(kd_looms_pts$color))
)[as.numeric(unique(filter(kd_looms_pts, lng_group_code == "Tsw")$color))]
Tsw_labels <- unique(filter(kd_looms_pts, lng_group_code == "Tsw")$loom_type) |>
  str_wrap(20) |>
  paste0("\n", Tsw_imgs) |>
  str_replace_all("\\n", "<br/>")

kd_looms_Tsw_map <- bg_map +
  geom_mark_ellipse(
    data = filter(
      kd_looms_pts,
      lng_group_code == "Tsw" & loom_type_code == "BFYRH"
    ) |>
      df_spatial(),
    aes(x = x, y = y, color = color, fill = color),
    alpha = .2,
    expand = .02,
    linewidth = lwd * 3,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_star(
    data = filter(kd_looms_pts, lng_group_code == "Tsw") |> df_spatial(),
    aes(x = x, y = y, color = color, starshape = color, fill = color),
    alpha = .85,
    size = 2
  ) +
  scale_fill_identity(name = "Loom type") +
  scale_color_identity(name = "Loom type") +
  scale_starshape_manual(
    values = Tsw_shapes,
    name = "Loom type",
    labels = Tsw_labels,
    guide = guide_legend(
      override.aes = list(
        size = 4,
        fill = Tsw_colors,
        color = Tsw_colors
      ),
      theme = theme(
        legend.key.spacing.y = unit(1, "line"),
        legend.text = element_markdown()
      )
    )
  ) +
  coord_sf(crs = prj, expand = FALSE)

kd_looms_Tsw_map <- set_dim(kd_looms_Tsw_map, get_dim(kd_looms_map))
ggsave(here("output/figures/kd-looms_Tsw_map.pdf"),
  kd_looms_Tsw_map,
  device = cairo_pdf, width = wd, height = wd * 2, units = "cm"
)
plot_crop(here("output/figures/kd-looms_Tsw_map.pdf"))


# Model comparison --------------------------------------------------------

models_summary <- read_csv(here("output/data/models_summary.csv"))

models_summary |>
  filter(data == "looms" & type == "full") |>
  select(-data, -type) |>
  arrange(-ML) |>
  mutate(ML = round(ML)) |>
  mutate(` ` = ifelse(ML == max(ML), "\\ding{43}", ""), .before = everything()) |>
  kbl(
    digits = 2,
    format = "latex", booktabs = TRUE, escape = FALSE
  ) |>
  write_lines(here("output/tables/kd-looms_models_summary.tex"))

models_summary |>
  filter(data == "lgs" & type == "full") |>
  select(-data, -type) |>
  arrange(-ML) |>
  mutate(ML = round(ML)) |>
  mutate(` ` = ifelse(ML == max(ML), "\\ding{43}", ""), .before = everything()) |>
  kbl(
    digits = 2,
    format = "latex", booktabs = TRUE, escape = FALSE
  ) |>
  write_lines(here("output/tables/kd-lgs_models_summary.tex"))
