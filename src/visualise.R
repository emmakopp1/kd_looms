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

kd_looms_languages <- read_csv(here("data/kd_looms_languages.csv")) |>
  mutate(loom_type = str_replace(loom_type, ", ", ",\n")) |>
  mutate(loom_type = fct_rev(loom_type)) |>
  mutate(lng_group = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])"))


# Consensus trees for looms ---------------------------------------------------------------------------------------

cs_tree1000 <- read.tree(here("output/trees/kd_loom1000_consensus.tree"))
# cs_tree1000_edges$root.edge <- 0
cs_tree1000$node.label <- round(as.numeric(cs_tree1000$node.label), 2) * 100
cs_tree1000$node.label[1] <- NA
cs_tree1000$tip.label <- str_replace_all(cs_tree1000$tip.label, "_", " ")
# cs_tree1000$edge.length <- cs_tree1000$edge.length[-length(cs_tree1000$edge.length)]

cs_tree1000_plot <- cs_tree1000 |>
  fortify() |>
  left_join(groups_types, by = join_by(label == group)) |>
  ggtree(ladderize = TRUE, size = lwd) +
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
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
  left_join(groups_types, by = join_by(label == group)) |>
  ggtree(ladderize = FALSE, size = lwd) +
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
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
  left_join(groups_types, by = join_by(label == group)) |>
  ggtree(ladderize = TRUE, size = lwd) +
  geom_tiplab(aes(fill = loom_type), geom = "label", label.size = 0, label.padding = unit(.15, "lines"), family = base_font, size = base_font_size / .pt, alpha = .95) +
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
  xlab("age (years BP)") +
  ylab("group") +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "line"), aspect.ratio = 0.618)
ggsave(here("output/figures/kd_lgs_ages_plot.pdf"), kd_lgs_ages_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd_lgs_ages_plot.pdf"))


# Cophylogeny -----------------------------------------------------------------------------------------------------

kd_lgs_cs <- read.tree(here("output/trees/kd_lgs_consensus.tree"))
kd_lgs_cs$root.edge <- 0

p1 <- ggtree(kd_lgs_cs, ladderize = TRUE)
p2 <- ggtree(cs_tree1111, ladderize = TRUE)
d1 <- p1$data |>
  mutate(language = label) |>
  left_join(kd_looms_languages) |>
  mutate(lng_group = str_extract(language, "^[A-Z][a-z]+(?=[A-Z])"))
lnggroup_loom <- tribble(
  ~lng_group, ~loom_code,
  "Be", NA,
  "Kra", NA,
  "Hlai", "FBBS",
  "Ks", "BFYRH",
  "Tn", NA,
  "Tsw", NA,
  "Tc", "BFSRH"
) |>
  full_join(distinct(d1, lng_group)) |>
  filter(!is.na(lng_group)) |>
  full_join(distinct(kd_looms_languages, loom_code, loom_type)) |>
  arrange(loom_type) |>
  mutate(lng_col = c(few_pal("Dark")(8), "black", "grey30", "grey50")) |>
  filter(!is.na(lng_group))

d2 <- p2$data |>
  mutate(group = label) |>
  left_join(kd_looms_languages)
ry <- filter(d1, !is.na(group) & !is.na(language))$y
d2$x <- ((d2$x - min(d2$x)) / (max(d2$x) - min(d2$x))) * (max(d1$x) - min(d1$x)) + min(d1$x)
d2$x <- max(d2$x) - d2$x + max(d1$x)
d2$x <- d2$x + (max(c(d1$x, d2$x)) - min(c(d1$x, d2$x))) / 100 * 10
d2$y <- ((d2$y - min(d2$y)) / (max(d2$y) - min(d2$y))) * (max(ry) - min(ry)) + min(ry)

pp <- p1 + geom_tree(data = d2)
dd <- bind_rows(d1, d2) %>%
  filter(!is.na(group) & !is.na(language))

pp + geom_line(aes(x, y, group = language), data = dd, color = "grey") +
  geom_tippoint(data = left_join(d1, select(lnggroup_loom, lng_group, lng_col)), aes(color = lng_col), size = 3) +
  scale_color_identity(guide = "legend", name = "Language group", labels = lnggroup_loom$lng_group) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(data = d2, aes(color = loom_type), size = 3) +
  scale_color_few(palette = "Light") +
  theme(legend.position = "bottom")


kd_lgs_cs$tip.label
