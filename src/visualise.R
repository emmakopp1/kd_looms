library(here)
library(phangorn)
library(phytools)
library(tidyverse)
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

# plt <- ptol_pal()(12)
plt <- tableau_color_pal('Classic 10 Medium')(10)
plt2 <- tableau_color_pal('Classic 10')(10)

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
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_manual(values = plt) +
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
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_manual(values = plt) +
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
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
  geom_nodelab(family = base_font, size = (base_font_size - 1) / .pt, hjust = 1.5, vjust = -.5) +
  geom_rootedge(.25, linewidth = lwd) +
  coord_cartesian(clip = "off", expand = FALSE) +
  scale_fill_manual(values = plt) +
  guides(fill = guide_legend(title = "Loom type", override.aes = aes(label = "     "))) +
  theme(plot.margin = margin(0, 4.9, 0, 0, unit = "line")) +
  xtheme
ggsave(here("output/figures/cs_tree8421.pdf"), cs_tree8421_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cs_tree8421.pdf"))

# library(patchwork)
# (cs_tree1000_plot + cs_tree1111_plot + cs_tree8421_plot) & theme(plot.margin = margin(1, 5, 0, 0, unit = "line"))
# ggsave(here("output/figures/cs.pdf"), device = cairo_pdf, width = wd*3.5, height = wd*2, units = "cm")
# plot_crop(here("output/figures/cs.pdf"))


# Age density distribution for languages --------------------------------------------------------------------------

kd_lgs_ages <- read_csv(here("output/kd_lgs_ages.csv"))

kd_lgs_ages_plot <- kd_lgs_ages |>
  mutate(group = fct(group, levels = c("Kra-Dai", "Kam-Tai", "Tai-Yay"))) |>
  ggplot(aes(x = age * 1000, y = group, height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white", fill = few_pal("Light")(2)[1]) +
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

kd_looms_cs <- cs_tree1111 |> 
  fortify() |> 
  left_join(kd_looms_languages, by = join_by(label == group)) |> 
  mutate(label = language) |> 
  as.phylo()
  
kd_cophylo <- cophylo(kd_lgs_cs, kd_looms_cs, methods = c("pre","post"), rotate.multi = TRUE)

p1 <- ggtree(kd_cophylo$trees[[1]], ladderize = FALSE)
p2 <- ggtree(kd_cophylo$trees[[2]], ladderize = FALSE)
#   flip(32, 43) |>
#   rotate(35) |>
#   rotate(54) |>
#   flip(50, 44) |>
#   rotate(50) |>
#   rotate(44) |>
#   rotate(45) |>
#   rotate(55)
d1 <- p1$data |>
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
  full_join(distinct(d1, lng_group)) |>
  filter(!is.na(lng_group)) |>
  full_join(distinct(kd_looms_languages, loom_code, loom_type)) |>
  mutate(loom_type = str_replace_all(loom_type, "\n", " ")) |> 
  mutate(loom_type = str_wrap(loom_type, 20)) |> 
  mutate(loom_type = fct_rev(loom_type)) |>
  arrange(loom_type, lng_group_name) |>
  mutate(lng_col = plt2) |>
  filter(!is.na(lng_group)) |> 
  mutate(lng_col = fct_inorder(lng_col), lng_group = fct_inorder(lng_group))
d1 <- left_join(d1, select(lnggroup_loom, lng_group, lng_col))
kd_looms_languages <- kd_looms_languages |> 
  mutate(loom_type = str_replace_all(loom_type, "\n", " ")) |> 
  mutate(loom_type = str_wrap(loom_type, 20)) |> 
  mutate(loom_type = fct_rev(loom_type))
d2 <- p2$data |>
  mutate(language = label) |>
  left_join(kd_looms_languages)
# ry <- filter(d1, !is.na(group) & !is.na(language))$y
ry <- d1$y
d2$x <- ((d2$x - min(d2$x)) / (max(d2$x) - min(d2$x))) * (max(d1$x) - min(d1$x)) + min(d1$x)
d2$x <- max(d2$x) - d2$x + max(d1$x)
d2$x <- d2$x + (max(c(d1$x, d2$x)) - min(c(d1$x, d2$x))) / 100 * 15
d2$y <- ((d2$y - min(d2$y)) / (max(d2$y) - min(d2$y))) * (max(ry) - min(ry)) + min(ry)

pp <- p1 + #|> flip(105, 147) +#|> flip(106,132) +
  geom_tree(data = d2)
dd <- bind_rows(d1, d2) %>%
  filter(!is.na(group) & !is.na(language))

pp + geom_line(aes(x, y, group = language), data = dd, color = "grey") +
  # geom_nodelab(aes(label = node)) +
  geom_tippoint(data = d1, aes(color = lng_col), size = 1.5) +
  scale_color_identity(guide = guide_legend(order = 1, position = "left", override.aes = list(size = 4), theme = theme(legend.key.spacing.y = unit(0, "line"))), name = "Language group", labels = lnggroup_loom$lng_group_name) +
  xtheme +
  ggnewscale::new_scale_colour() +
  # geom_nodelab(data = d2, aes(label = node)) +
  geom_tippoint(data = d2, aes(color = loom_type), size = 1.75) +
  # scale_color_few(palette = "Light", name = "Loom type", guide = guide_legend(order = 2, override.aes = list(size = 4))) +
  scale_color_manual(values = plt, name = "Loom type", guide = guide_legend(order = 2, position = "right", override.aes = list(size = 4))) +
  theme(#legend.position = c(0.5, 1), legend.justification = c(.5, 1), legend.box = "horizontal", 
    legend.spacing.x = unit(14, "line"), aspect.ratio = 2)
ggsave(here("output/figures/cophylogeny.pdf"), device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/cophylogeny.pdf"))

