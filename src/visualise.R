library(here)
library(tidyverse)
library(phangorn)
library(ggtree)
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")
library(ggthemes)
library(ggridges)
library(knitr)

dir.create(here("output/figures"))
wd <- 14
lwd <- 1 / .pt
base_font <- "Noto Sans Condensed"
base_font_size <- 9
theme_set(
  theme_minimal(base_family = base_font, base_size = base_font_size) +
    theme(
      legend.text = element_text(family = base_font)
      #       legend.position = "bottom",
      #       aspect.ratio = .618
    )
)

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

loom_groups_types <- read_csv(here("data/loom_groups_types.csv"))

cs_tree1000 <- read.tree(here("output/trees/kd_loom1000_consensus.tree"))
# cs_tree1000_edges$root.edge <- 0
cs_tree1000$node.label <- round(as.numeric(cs_tree1000$node.label), 2) * 100
cs_tree1000$node.label[1] <- NA
cs_tree1000$tip.label <- str_replace_all(cs_tree1000$tip.label, "_", " ")
# cs_tree1000$edge.length <- cs_tree1000$edge.length[-length(cs_tree1000$edge.length)]

cs_tree1000_plot <- cs_tree1000 |>
  fortify() |>
  left_join(loom_groups_types, by = join_by(label == group)) |>
  mutate(type_label = str_replace(type_label, ", ", ",\n")) |>
  ggtree(ladderize = TRUE, size = lwd) +
  geom_tiplab(aes(fill = fct_rev(type_label)), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_few(palette = "Light") +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  theme(plot.margin = margin(0, 2.4, 0, 0, unit = "line")) +
  cs_theme
ggsave(here("output/figures/cs_tree1000.pdf"), cs_tree1000_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree1000.pdf"))

cs_tree1111 <- read.tree(here("output/trees/kd_loom1111_consensus.tree"))
cs_tree1111$node.label <- round(as.numeric(cs_tree1111$node.label), 2) * 100
cs_tree1111$node.label[1] <- NA
cs_tree1111$tip.label <- str_replace_all(cs_tree1111$tip.label, "_", " ")

cs_tree1111_plot <- cs_tree1111 |>
  fortify() |>
  left_join(loom_groups_types, by = join_by(label == group)) |>
  mutate(type_label = str_replace(type_label, ", ", ",\n")) |>
  ggtree(ladderize = FALSE, size = lwd) +
  geom_tiplab(aes(fill = fct_rev(type_label)), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_few(palette = "Light") +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  cs_theme +
  theme(plot.margin = margin(0, 4.98, 0, 0, unit = "line"))
# cs_tree1111_plot + geom_nodelab(aes(label = node))
cs_tree1111_plot <- flip(cs_tree1111_plot, 32, 43)
ggsave(here("output/figures/cs_tree1111.pdf"), cs_tree1111_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree1111.pdf"))

cs_tree8421 <- read.tree(here("output/trees/kd_loom8421_consensus.tree"))
cs_tree8421$node.label <- round(as.numeric(cs_tree8421$node.label), 2) * 100
cs_tree8421$node.label[1] <- NA
cs_tree8421$tip.label <- str_replace_all(cs_tree8421$tip.label, "_", " ")

cs_tree8421_plot <- cs_tree8421 |>
  fortify() |>
  left_join(loom_groups_types, by = join_by(label == group)) |>
  mutate(type_label = str_replace(type_label, ", ", ",\n")) |>
  ggtree(ladderize = TRUE, size = lwd) +
  geom_tiplab(aes(fill = fct_rev(type_label)), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_few(palette = "Light") +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  cs_theme +
  theme(plot.margin = margin(0, 4.9, 0, 0, unit = "line"))
ggsave(here("output/figures/cs_tree8421.pdf"), cs_tree8421_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree8421.pdf"))

# library(patchwork)
# (cs_tree1000_plot + cs_tree1111_plot + cs_tree8421_plot) & theme(plot.margin = margin(1, 5, 0, 0, unit = "line"))
# ggsave(here("output/figures/cs.pdf"), device = cairo_pdf, width = wd*3.5, height = wd*2, units = "cm")
# plot_crop(here("output/figures/cs.pdf"))

kd_lgs_ages <- read_csv(here("output/kd_lgs_ages.csv"))

kd_lgs_ages_plot <- kd_lgs_ages |>
  mutate(group = fct(group, levels = c("Kra-Dai", "Kam-Tai", "Tai-Yay"))) |>
  ggplot(aes(x = age * 1000, y = group, height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white", fill = few_pal("Light")(2)[1]) +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_x_reverse() +
  # scale_x_reverse(limits = c(15000, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  xlab("age (years BP)") +
  ylab("group") +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "line"), aspect.ratio = 0.618)
ggsave(here("output/figures/kd_lgs_ages_plot.pdf"), kd_lgs_ages_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd_lgs_ages_plot.pdf"))
