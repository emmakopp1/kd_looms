library(here)
library(tidyverse)
library(phangorn)
library(phytools)
library(treeio)
library(ggtree)
library(ggthemes)
library(knitr)

dir.create(here("output/figures"))
wd <- 14
base_font <- "Noto Sans Condensed"
base_font_size <- 10
theme_set(
  theme_minimal(base_family = base_font, base_size = base_font_size) +
    theme(
      legend.text = element_text(family = base_font)
      #       legend.position = "bottom",
      #       aspect.ratio = .618
    )
)

loom_groups_types <- read_csv(here("data/loom_groups_types.csv"))

tree1000 <- read.nexus(here("data/by_level/loom1000/kd_loom1000.trees"))
tree1000 <- tree1000[seq(2, length(tree1000), by = round(length(tree1000) / 1000))]
cs_tree1000 <- consensus(tree1000, p = .5, rooted = TRUE)
cs_tree1000_edges <- consensus.edges(tree1000, consensus.tree = cs_tree1000, rooted = TRUE)
cs_tree1000_edges$root.edge <- 0
cs_tree1000_edges$node.label <- round(as.numeric(cs_tree1000_edges$node.label), 2) * 100
cs_tree1000_edges$node.label[1] <- NA
cs_tree1000_edges$tip.label <- str_replace_all(cs_tree1000_edges$tip.label, "_", " ")
cs_tree1000_edges$edge.length <- cs_tree1000_edges$edge.length[-length(cs_tree1000_edges$edge.length)]


cs_theme <- theme(
  aspect.ratio = 1.15,
  legend.position = "inside",
  legend.justification = c(0, 1),
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
  legend.key.spacing.y = unit(.25, "line"),
  legend.margin = margin(0, 0, 0, 0, unit = "line"),
  legend.box.margin = margin(0, 0, 0, 0, unit = "line"),
  legend.box.spacing = unit(0, "pt")
)

cs_tree1000_plot <- cs_tree1000_edges |>
  fortify() |>
  left_join(loom_groups_types, by = join_by(label == group)) |>
  mutate(type_label = str_replace(type_label, ", ", ",\n")) |>
  ggtree(ladderize = TRUE) +
  geom_tiplab(aes(fill = fct_rev(type_label)), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_few(palette = "Light") +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  theme(plot.margin = margin(0, 2.5, 0, 0, unit = "line")) +
  cs_theme
ggsave(here("output/figures/cs_tree1000.pdf"), cs_tree1000_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree1000.pdf"))

tree1111 <- read.nexus(here("data/by_level/loom1111/kd_loom1111.trees"))
tree1111 <- tree1111[seq(2, length(tree1111), by = round(length(tree1111) / 1000))]
cs_tree1111 <- consensus(tree1111, p = .5, rooted = TRUE)
cs_tree1111_edges <- consensus.edges(tree1111, consensus.tree = cs_tree1111, rooted = TRUE)
cs_tree1111_edges$root.edge <- 0
cs_tree1111_edges$node.label <- round(as.numeric(cs_tree1111_edges$node.label), 2) * 100
cs_tree1111_edges$node.label[1] <- NA
cs_tree1111_edges$tip.label <- str_replace_all(cs_tree1111_edges$tip.label, "_", " ")
cs_tree1111_edges$edge.length <- cs_tree1111_edges$edge.length[-length(cs_tree1111_edges$edge.length)]

cs_tree1111_plot <- cs_tree1111_edges |>
  fortify() |>
  left_join(loom_groups_types, by = join_by(label == group)) |>
  ggtree(ladderize = TRUE) +
  geom_tiplab(aes(fill = fct_rev(type)), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_few(palette = "Light") +
  theme(plot.margin = margin(0, 4, 0, 0, unit = "line"), aspect.ratio = 1.15, legend.position = "none")
ggsave(here("output/figures/cs_tree1111.pdf"), cs_tree1111_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree1111.pdf"))

tree8421 <- read.nexus(here("data/by_level/loom8421/kd_loom8421.trees"))
tree8421 <- tree8421[seq(2, length(tree8421), by = round(length(tree8421) / 1000))]
cs_tree8421 <- consensus(tree8421, p = .5, rooted = TRUE)
cs_tree8421_edges <- consensus.edges(tree8421, consensus.tree = cs_tree8421, rooted = TRUE)
cs_tree8421_edges$root.edge <- 0
cs_tree8421_edges$node.label <- round(as.numeric(cs_tree8421_edges$node.label), 2) * 100
cs_tree8421_edges$node.label[1] <- NA
cs_tree8421_edges$tip.label <- str_replace_all(cs_tree8421_edges$tip.label, "_", " ")
cs_tree8421_edges$edge.length <- cs_tree8421_edges$edge.length[-length(cs_tree8421_edges$edge.length)]

cs_tree8421_plot <- cs_tree8421_edges |>
  fortify() |>
  left_join(loom_groups_types, by = join_by(label == group)) |>
  ggtree(ladderize = TRUE) +
  geom_tiplab(aes(fill = fct_rev(type)), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_few(palette = "Light") +
  theme(plot.margin = margin(0, 5, 0, 0, unit = "line"), aspect.ratio = 1.15, legend.position = "none")
ggsave(here("output/figures/cs_tree8421.pdf"), cs_tree8421_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree8421.pdf"))

# library(patchwork)
# (cs_tree1000_plot + cs_tree1111_plot + cs_tree8421_plot) & theme(plot.margin = margin(1, 5, 0, 0, unit = "line"))
# ggsave(here("output/figures/cs.pdf"), device = cairo_pdf, width = wd*3.5, height = wd*2, units = "cm")
# plot_crop(here("output/figures/cs.pdf"))
