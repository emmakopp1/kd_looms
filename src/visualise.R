library(here)
library(phangorn)
library(phytools)
library(tidyverse)
library(stringi)
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
xtheme <- theme(
  aspect.ratio = 1.15,
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
  legend.position = c(0, 1),
  legend.justification = c(0, 1),
  legend.key.spacing.y = unit(.25, "line"),
  legend.margin = margin(0, 0, 0, 0, unit = "line"),
  legend.box.margin = margin(0, 0, 0, 0, unit = "line"),
  legend.box.spacing = unit(0, "pt")
  #       legend.position = "bottom",
  #       aspect.ratio = .618
)
theme_set(xtheme)

plt <- tableau_color_pal("Classic 10 Medium")(10)
plt2 <- tableau_color_pal("Classic 10")(10)
plt3 <- tableau_color_pal("Classic 10 Light")(10)
xcol <- few_pal("Light")(2)[1]

kd_looms_languages <- read_csv(here("data/kd_looms_languages.csv")) |>
  mutate(loom_type = str_replace(loom_type, ", ", ",\n")) |>
  mutate(loom_type = fct_rev(loom_type)) |>
  mutate(lng_group = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])")) |>
  arrange(loom_type, lng_group)


# Consensus trees for looms ---------------------------------------------------------------------------------------

cs_tree1000 <- read.tree(here("output/trees/kd_loom1000_consensus.tree"))
# cs_tree1000_edges$root.edge <- 0
cs_tree1000$node.label <- round(as.numeric(cs_tree1000$node.label), 2) * 100
cs_tree1000$node.label[1] <- NA
cs_tree1000$tip.label <- str_replace_all(cs_tree1000$tip.label, "_", " ")
# cs_tree1000$edge.length <- cs_tree1000$edge.length[-length(cs_tree1000$edge.length)]

cs_tree1000_plot <- cs_tree1000 |>
  fortify() |>
  left_join(kd_looms_languages, by = join_by(label == group)) |>
  ggtree(ladderize = TRUE, size = lwd) +
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_manual(values = plt3) +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  theme(plot.margin = margin(0, 2.4, 0, 0, unit = "line")) +
  xtheme
ggsave(here("output/figures/cs_tree1000.pdf"), cs_tree1000_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree1000.pdf"))

cs_tree1111 <- read.tree(here("output/trees/kd_loom1111_consensus.tree"))
cs_tree1111$node.label <- round(as.numeric(cs_tree1111$node.label), 2) * 100
cs_tree1111$node.label[1] <- NA
cs_tree1111$tip.label <- str_replace_all(cs_tree1111$tip.label, "_", " ")

cs_tree1111_plot <- cs_tree1111 |>
  fortify() |>
  left_join(kd_looms_languages, by = join_by(label == group)) |>
  ggtree(ladderize = FALSE, size = lwd) +
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_manual(values = plt3) +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  theme(plot.margin = margin(0, 4.98, 0, 0, unit = "line")) +
  xtheme
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
  left_join(kd_looms_languages, by = join_by(label == group)) |>
  ggtree(ladderize = TRUE, size = lwd) +
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_manual(values = plt3) +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  theme(plot.margin = margin(0, 4.9, 0, 0, unit = "line")) +
  xtheme
ggsave(here("output/figures/cs_tree8421.pdf"), cs_tree8421_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree8421.pdf"))


# Age density distribution for languages --------------------------------------------------------------------------

kd_lgs_ages <- read_csv(here("output/kd_lgs_ages.csv"))

kd_lgs_ages_plot <- kd_lgs_ages |>
  mutate(group = fct(group, levels = c("Kra-Dai", "Kam-Tai", "Tai-Yay"))) |>
  ggplot(aes(x = age * 1000, y = group, height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white", fill = xcol) +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_x_reverse() +
  # scale_x_reverse(limits = c(15000, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  xlab("Age (years BP)") +
  ylab("Language group") +
  theme_minimal(base_size = base_font_size, base_family = base_font) +
  xtheme +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "line"), aspect.ratio = 0.618)
ggsave(here("output/figures/kd_lgs_ages_plot.pdf"), kd_lgs_ages_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd_lgs_ages_plot.pdf"))


# Cophylogeny -----------------------------------------------------------------------------------------------------

kd_lgs_cs <- read.tree(here("output/trees/kd_lgs_consensus.tree"))
kd_lgs_cs$root.edge <- 0

kd_loom_pb <- c("Dai Huayao", "Dai Yuxi Yuanjiang", "Dai Jinghong", "Zhuang Napo", "Zhuang Longzhou", "Nung An", "Tai Phake")

kd_looms_cs <- cs_tree1111 |>
  fortify() |>
  left_join(kd_looms_languages, by = join_by(label == group)) |>
  mutate(label = ifelse(label %in% kd_loom_pb, label, language)) |>
  as.phylo()

kd_cophylo <- cophylo(kd_lgs_cs, kd_looms_cs, methods = c("pre", "post"), rotate.multi = TRUE)
kd_cophylo$trees[[2]]$tip.label <- stringi::stri_replace_all_fixed(
  str = kd_cophylo$trees[[2]]$tip.label, pattern = kd_looms_languages$group,
  replacement = kd_looms_languages$language,
  vectorise_all = FALSE
)
kd_cophylo$trees[[2]] <- kd_cophylo$trees[[2]] |>
  fortify() |>
  left_join(kd_looms_languages, by = join_by(label == language)) |>
  mutate(label = group) |>
  as.phylo()

kd_lng_tree <- ggtree(kd_cophylo$trees[[1]], ladderize = FALSE, size = lwd)
# flip(113, 128) |>
# flip(103,172)
kd_loom_tree <- ggtree(kd_cophylo$trees[[2]], ladderize = FALSE, size = lwd)
# flip(1, 33)
# kd_lng_tree <- ggtree(kd_cophylo$trees[[1]], ladderize = FALSE, size = lwd, branch.length = "none")
# kd_loom_tree <- ggtree(kd_cophylo$trees[[2]], ladderize = FALSE, size = lwd, branch.length = "none")
#   flip(32, 43) |>
#   rotate(35) |>
#   rotate(54) |>
#   flip(50, 44) |>
#   rotate(50) |>
#   rotate(44) |>
#   rotate(45) |>
#   rotate(55)
kd_lng_tree_data <- kd_lng_tree$data |>
  mutate(language = label) |>
  left_join(kd_looms_languages) |>
  mutate(lng_group = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])"))
lnggroup_loom <- tribble(
  ~lng_group, ~lng_group_name, ~loom_code,
  "Be", "Ong-Be", NA,
  "Kra", "Kra", NA,
  "Hlai", "Hlai", "FBBS",
  "Ks", "Kam-Sui", "BFYRH",
  "Tn", "Northern Tai", NA,
  "Tsw", "Southwestern Tai", "2LH",
  "Tc", "Central Tai", "BFSRH"
) |>
  full_join(distinct(kd_lng_tree_data, lng_group)) |>
  filter(!is.na(lng_group)) |>
  full_join(distinct(kd_looms_languages, loom_code, loom_type)) |>
  mutate(loom_type = str_replace_all(loom_type, "\n", " ")) |>
  mutate(loom_type = str_wrap(loom_type, 20)) |>
  mutate(loom_type = fct_rev(loom_type)) |>
  arrange(loom_type, lng_group_name) |>
  mutate(lng_col = plt2) |>
  filter(!is.na(lng_group)) |>
  mutate(lng_col = fct_inorder(lng_col), lng_group = fct_inorder(lng_group))
kd_lng_tree_data <- left_join(kd_lng_tree_data, select(lnggroup_loom, lng_group, lng_col))
kd_looms_languages <- kd_looms_languages |>
  mutate(loom_type = str_replace_all(loom_type, "\n", " ")) |>
  mutate(loom_type = str_wrap(loom_type, 20)) |>
  mutate(loom_type = fct_rev(loom_type))
kd_loom_tree_data <- kd_loom_tree$data |>
  mutate(group = label) |>
  left_join(kd_looms_languages)
# ry <- filter(kd_lng_tree_data, !is.na(group) & !is.na(language))$y
ry <- kd_lng_tree_data$y
kd_loom_tree_data$x <- ((kd_loom_tree_data$x - min(kd_loom_tree_data$x)) / (max(kd_loom_tree_data$x) - min(kd_loom_tree_data$x))) * (max(kd_lng_tree_data$x) - min(kd_lng_tree_data$x)) + min(kd_lng_tree_data$x)
kd_loom_tree_data$x <- max(kd_loom_tree_data$x) - kd_loom_tree_data$x + max(kd_lng_tree_data$x)
kd_loom_tree_data$x <- kd_loom_tree_data$x + (max(c(kd_lng_tree_data$x, kd_loom_tree_data$x)) - min(c(kd_lng_tree_data$x, kd_loom_tree_data$x))) / 100 * 15
kd_loom_tree_data$y <- ((kd_loom_tree_data$y - min(kd_loom_tree_data$y)) / (max(kd_loom_tree_data$y) - min(kd_loom_tree_data$y))) * (max(ry) - min(ry)) + min(ry)

kd_lng_loom_tree <- kd_lng_tree +
  # flip(105, 147) +#|> flip(106,132) +
  geom_tree(data = kd_loom_tree_data, linewidth = lwd)
kd_lng_loom_tree_data <- bind_rows(kd_lng_tree_data, kd_loom_tree_data) %>%
  filter(!is.na(group) & !is.na(language)) |>
  mutate(pb = group %in% kd_loom_pb)

kd_loom_hl <- c(
  getMRCA(kd_cophylo$trees[[2]], c("Dai Huayao", "Dai Yuxi Yuanjiang", "Dai Jinghong")),
  getMRCA(kd_cophylo$trees[[2]], c("Zhuang Napo", "Zhuang Longzhou", "Nung An")),
  which(kd_cophylo$trees[[2]]$tip.label == "Tai Phake")
)

kd_cophylo_plot <- kd_lng_loom_tree +
  # geom_tree(data = kd_lng_tree_data, aes(color = lng_col)) +
  geom_rootedge(.25, linewidth = lwd) +
  geom_segment(data = filter(kd_loom_tree_data, parent == node), aes(x = x, xend = x + .25, y = y, yend = y), linewidth = lwd) +
  geom_line(data = kd_lng_loom_tree_data, aes(x, y, group = language, linetype = pb), color = "grey") +
  # geom_nodelab(aes(label = node), size = 3) +
  # geom_tiplab(size = 3, family = base_font) +
  geom_tippoint(data = kd_lng_tree_data, aes(color = lng_col, fill = lng_col, shape = !is.na(group), size = !is.na(group))) +
  scale_color_identity(guide = guide_legend(order = 1, position = "inside", override.aes = list(size = 4), theme = theme(legend.key.spacing.y = unit(0, "line"))), name = "Language group", labels = lnggroup_loom$lng_group_name) +
  xtheme +
  scale_fill_identity(guide = "none") +
  ggnewscale::new_scale_colour() +
  ggnewscale::new_scale_fill() +
  # geom_nodelab(data = kd_loom_tree_data, aes(label = node), size = 3) +
  # geom_tiplab(data = kd_loom_tree_data, aes(label = node), size = 3, family = base_font) +
  geom_tippoint(data = kd_loom_tree_data, aes(color = loom_type, fill = loom_type, shape = group %in% kd_lng_tree_data$group, size = group %in% kd_lng_tree_data$group)) +
  scale_color_manual(values = plt3, name = "Loom type", guide = guide_legend(order = 2, position = "inside", override.aes = list(size = 4))) +
  scale_fill_manual(values = plt3, guide = "none") +
  scale_shape_manual(values = c(21, 23), guide = "none") +
  scale_size_manual(values = c(.75, 2), guide = "none") +
  scale_linetype(guide = "none") +
  scale_y_reverse() +
  # geom_hilight(data = kd_loom_tree_data, aes(subset = node %in% kd_loom_hl, node = node), fill = NA, color = "grey50", linetype = "dotted", extend = 0, to.bottom = TRUE) +
  theme(legend.spacing.x = unit(14, "line"), aspect.ratio = 1.25, legend.position.inside = c(.5, .96), legend.justification.inside = c(.5, 1), legend.box = "horizontal", legend.background = element_blank())
ggsave(here("output/figures/kd_cophylo_plot.pdf"), device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd_cophylo_plot.pdf"))


# Mutation rates --------------------------------------------------------------------------------------------------

library(ggforce)
library(ggdist)
mutationrate_bylevel_tb <- read_csv(here("output/mutationrate_bylevel.csv"))

mutationrate_bylevel_tb |>
  ggplot(aes(y = factor(level), x = rate)) +
  # geom_density_ridges(fill = NA, color = NA, scale = 1) +
  stat_density_ridges(aes(fill = .5 - abs(.5 - after_stat(ecdf))), geom = "density_ridges_gradient", calc_ecdf = TRUE, scale = 1, panel_scaling = FALSE, color = "gray50", linewidth = lwd) +
  # stat_slab(fill = xcol, limits = c(0,NA), justification = -.1) +
  # geom_boxplot(width = .2, linewidth = lwd, outliers = FALSE, outlier.size = .25, outlier.alpha = .5, outlier.color = "grey50", color = few_pal("Dark")(2)[2], staplewidth = .75) +
  stat_summary(geom = "text", fun = "median", aes(label = round(..x.., 2)), family = base_font, size = base_font_size / .pt, vjust = 2) +
  # stat_summary(geom = "point", fun = "mean", vjust = 2, color = few_pal("Dark")(2)[2], shape = 5) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black") +
  ylab("Level") +
  xlab("Mutation rate") +
  xlim(0, 2) +
  scale_fill_distiller(palette = "YlOrBr", name = "Tail probability") +
  # scale_fill_viridis_c(name = "Tail probability", option = "E", direction = 1) +
  theme_minimal(base_size = base_font_size, base_family = base_font) +
  xtheme +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "line"), aspect.ratio = 0.618, legend.position = "right")
