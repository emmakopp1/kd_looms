library(here)
library(tidyverse)
library(phangorn)
library(phytools)
library(treeio)
library(ggtree)
library(knitr)

dir.create(here("output/figures"))
wd <- 14
base_font <- "Noto Sans Condensed"
base_font_size <- 10
# theme_set(
#   theme_minimal(base_family = base_font, base_size = base_font_size) +
#     theme(
#       legend.position = "bottom",
#       aspect.ratio = .618
#     )
# )

tree1000 <- read.nexus(here("data/by_level/loom1000/kd_loom1000.trees"))
cs_tree1000 <- consensus(tree1000, p = .5, rooted = TRUE)
cs_tree1000_edges <- consensus.edges(tree1000, consensus.tree = cs_tree1000, rooted = TRUE)
cs_tree1000_edges$root.edge <- 0
cs_tree1000_edges$node.label <- round(as.numeric(cs_tree1000_edges$node.label), 2) * 100
cs_tree1000_edges$node.label[1] <- NA
cs_tree1000_edges$tip.label <- str_replace_all(cs_tree1000_edges$tip.label, "_", " ")
cs_tree1000_edges$edge.length <- cs_tree1000_edges$edge.length[-length(cs_tree1000_edges$edge.length)]
  
cs_tree1000_edges |> 
  ggtree(ladderize = TRUE) +
  geom_tiplab(family = base_font, size = base_font_size / .pt) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(0, 2.5, 0, 0, unit = "line"), aspect.ratio = 1)
ggsave(here("output/figures/cs_tree1000.pdf"), device = cairo_pdf, width = wd, height = wd, units = "cm")
plot_crop(here("output/figures/cs_tree1000.pdf"))

tree1111 <- read.nexus(here("data/by_level/loom1111/kd_loom1111.trees"))
cs_tree1111 <- consensus(tree1111, p = .5, rooted = TRUE)
cs_tree1111_edges <- consensus.edges(tree1111, consensus.tree = cs_tree1111, rooted = TRUE)
cs_tree1111_edges$root.edge <- 0
cs_tree1111_edges$node.label <- round(as.numeric(cs_tree1111_edges$node.label), 2) * 100
cs_tree1111_edges$node.label[1] <- NA
cs_tree1111_edges$tip.label <- str_replace_all(cs_tree1111_edges$tip.label, "_", " ")
cs_tree1111_edges$edge.length <- cs_tree1111_edges$edge.length[-length(cs_tree1111_edges$edge.length)]

cs_tree1111_edges |> 
  ggtree(ladderize = TRUE) +
  geom_tiplab(family = base_font, size = base_font_size / .pt) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(0, 2.5, 0, 0, unit = "line"), aspect.ratio = 1)
ggsave(here("output/figures/cs_tree1111.pdf"), device = cairo_pdf, width = wd, height = wd, units = "cm")
plot_crop(here("output/figures/cs_tree1111.pdf"))
