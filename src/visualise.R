library(here)
library(phangorn)
library(phytools)
library(sf)
library(rnaturalearth)
library(ggspatial)
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
library(ggtext)
library(knitr)
library(kableExtra)
library(patchwork)

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
  legend.position = c(0, 1),
  legend.justification = c(0, 1),
  legend.key.spacing.y = unit(.25, "line"),
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

lgs_order <- c("Hlai", "Kra", "Be", "Ks", "Tn", "Tc", "Tsw")
looms_order <- c("FBBS", "BFSRH", "BFYRH", "FCBunique", "BFcant", "FCBcant", "FCB")

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
  mutate(loom_type_code = factor(loom_type_code, levels = levels(plt_tb_looms$loom_type_code))) |>
  left_join(plt_tb_looms) |>
  arrange(loom_type_code, lng_group) |>
  mutate(loom_type = fct_inorder(loom_type))

kd_lgs <- read_csv(here("data/kd-lgs/kd-lgs_datapoints.csv")) |>
  mutate(label = paste0(lng_group_code, lng)) |>
  mutate(lng_group_code = factor(lng_group_code, levels = levels(plt_tb_lgs$lng_group_code))) |>
  left_join(plt_tb_lgs) |>
  arrange(lng_group_code, lng) |>
  mutate(lng_group = fct_inorder(lng_group))


# Consensus trees for looms ---------------------------------------------------------------------------------------

cs_tree <- function(tr) {
  ggtree(tr, ladderize = TRUE, size = lwd) +
    geom_tiplab(aes(fill = color),
      geom = "label",
      label.size = 0,
      label.padding = unit(.15, "lines"),
      family = base_font, size = base_font_size / .pt
    ) +
    geom_nodelab(
      family = base_font,
      size = (base_font_size - 1) / .pt,
      hjust = 1.5,
      vjust = -.5
    ) +
    geom_rootedge(.25, linewidth = lwd) +
    scale_y_reverse() +
    coord_cartesian(clip = "off", expand = FALSE) +
    scale_fill_identity(guide = guide_legend(), labels = levels(kd_looms$loom_type)) +
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

kd_lgs_bcov_cs_tree <- read.tree(here("output/trees/kd-lgs_bcov_consensus.tree"))
kd_lgs_bcov_cs_tree$node.label <- round(as.numeric(kd_lgs_bcov_cs_tree$node.label), 2) * 100
kd_lgs_bcov_cs_tree$root.edge.length <- 0
kd_lgs_bcov_cs_tree$node.label[1] <- NA
kd_lgs_bcov_cs_tree$tip.label <- str_replace_all(kd_lgs_bcov_cs_tree$tip.label, "_", " ")

kd_lgs_bcov_cs_tree_plot <- kd_lgs_bcov_cs_tree |>
  fortify() |>
  left_join(kd_lgs, by = join_by(label == label)) |>
  cs_tree() +
  scale_fill_identity(guide = guide_legend(), labels = levels(kd_lgs$lng_group)) +
  guides(fill = guide_legend(
    title = "Language group",
    override.aes = aes(label = "     ")
  )) +
  theme(
    aspect.ratio = 2,
    plot.margin = margin(0, 3.5, 0, 0, unit = "line")
  )
ggsave(here("output/figures/kd-lgs_bcov_cs_tree.pdf"), kd_lgs_bcov_cs_tree_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-lgs_bcov_cs_tree.pdf"))

kd_looms_bcov1000_cs_tree <- read.tree(here("output/trees/kd-looms_bcov1000_consensus.tree"))
kd_looms_bcov1000_cs_tree$node.label <- round(as.numeric(kd_looms_bcov1000_cs_tree$node.label), 2) * 100
kd_looms_bcov1000_cs_tree$node.label[1] <- NA
kd_looms_bcov1000_cs_tree$tip.label <- str_replace_all(kd_looms_bcov1000_cs_tree$tip.label, "_", " ")

kd_looms_bcov1000_cs_tree_plot <- kd_looms_bcov1000_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  cs_tree() +
  theme(plot.margin = margin(0, 2.5, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_bcov1000_cs_tree.pdf"), kd_looms_bcov1000_cs_tree_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-looms_bcov1000_cs_tree.pdf"))

kd_looms_bcov1111_cs_tree <- read.tree(here("output/trees/kd-looms_bcov1111_consensus.tree"))
kd_looms_bcov1111_cs_tree$node.label <- round(as.numeric(kd_looms_bcov1111_cs_tree$node.label), 2) * 100
kd_looms_bcov1111_cs_tree$node.label[1] <- NA
kd_looms_bcov1111_cs_tree$tip.label <- str_replace_all(kd_looms_bcov1111_cs_tree$tip.label, "_", " ")

kd_looms_bcov1111_cs_tree_plot <- kd_looms_bcov1111_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  cs_tree() +
  theme(plot.margin = margin(0, 4.98, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_bcov1111_cs_tree.pdf"), kd_looms_bcov1111_cs_tree_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-looms_bcov1111_cs_tree.pdf"))

kd_looms_bcov8421_cs_tree <- read.tree(here("output/trees/kd-looms_bcov8421_consensus.tree"))
kd_looms_bcov8421_cs_tree$node.label <- round(as.numeric(kd_looms_bcov8421_cs_tree$node.label), 2) * 100
kd_looms_bcov8421_cs_tree$node.label[1] <- NA
kd_looms_bcov8421_cs_tree$tip.label <- str_replace_all(kd_looms_bcov8421_cs_tree$tip.label, "_", " ")

kd_looms_bcov8421_cs_tree_plot <- kd_looms_bcov8421_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  cs_tree() +
  theme(plot.margin = margin(0, 4.9, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_bcov8421_cs_tree.pdf"), kd_looms_bcov8421_cs_tree_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-looms_bcov8421_cs_tree.pdf"))

kd_looms_ctmc4_cs_tree <- read.tree(here("output/trees/kd-looms_ctmc4_consensus.tree"))
kd_looms_ctmc4_cs_tree$node.label <- round(as.numeric(kd_looms_ctmc4_cs_tree$node.label), 2) * 100
kd_looms_ctmc4_cs_tree$node.label[1] <- NA
kd_looms_ctmc4_cs_tree$tip.label <- str_replace_all(kd_looms_ctmc4_cs_tree$tip.label, "_", " ")

kd_looms_ctmc4_cs_tree_plot <- kd_looms_ctmc4_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  cs_tree() +
  theme(plot.margin = margin(0, 4.25, 0, 0, unit = "line"))
ggsave(here("output/figures/kd-looms_ctmc4_cs_tree.pdf"), kd_looms_ctmc4_cs_tree_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-looms_ctmc4_cs_tree.pdf"))


# Age density distribution for languages --------------------------------------------------------------------------

kd_lgs_ages <- read_csv(here("output/data/kd-lgs_ages.csv"))

kd_lgs_ages_plot <- kd_lgs_ages |>
  mutate(group = fct(group, levels = c("Kra-Dai", "Kam-Tai", "Tai-Yay"))) |>
  ggplot(aes(x = age * 1000, y = group)) +
  stat_density_ridges(aes(fill = .5 - abs(.5 - after_stat(ecdf))), geom = "density_ridges_gradient", calc_ecdf = TRUE, scale = 1, panel_scaling = FALSE, color = "gray50", linewidth = lwd) +
  stat_summary(geom = "text", fun = "median", aes(label = round(..x..)), family = base_font, size = base_font_size / .pt, vjust = 1.5) +
  scale_fill_distiller(palette = "PuBu", direction = 1, limits = c(0, .5), name = "Tail\nprobability") +
  xlab("Age (years BP)") +
  ylab("Language group") +
  scale_x_reverse(limits = c(10000, 0)) +
  # scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  xlab("Age (years BP)") +
  ylab("Language group") +
  theme_minimal(base_size = base_font_size, base_family = base_font) +
  xtheme +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "line"), aspect.ratio = 0.618, legend.position = "right")
ggsave(here("output/figures/kd-lgs_ages_plot.pdf"), kd_lgs_ages_plot, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-lgs_ages_plot.pdf"))

kd_lgs_ages_summary <- read_csv(here("output/data/kd-lgs_ages_summary.csv"))

kd_lgs_ages_summary |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  unite(hdi, hdi_lower, hdi_upper, sep = ", ") |>
  mutate(hdi = paste0("[", hdi, "]")) |>
  rename(languages = n_lngs, `95% HPDI` = hdi) |>
  kbl(
    digits = 2,
    format = "latex", booktabs = TRUE
  ) |>
  write_lines(here("output/tables/kd-lgs_ages_summary.tex"))


# Cophylogeny -----------------------------------------------------------------------------------------------------

kd_lgs_cs <- read.tree(here("output/trees/kd-lgs_bcov_consensus.tree")) |>
  fortify() |>
  left_join(kd_lgs) |>
  mutate(label = lng) |>
  as.phylo()
kd_lgs_cs$root.edge <- 0

kd_loom_pb <- c("Dai Huayao", "Dai Yuxi Yuanjiang", "Dai Jinghong", "Zhuang Napo", "Zhuang Longzhou", "Nung An", "Tai Phake")

kd_looms_cs <- kd_looms_bcov1111_cs_tree |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == group)) |>
  mutate(label = ifelse(label %in% kd_loom_pb, label, lng)) |>
  as.phylo()

kd_cophylo <- cophylo(kd_lgs_cs, kd_looms_cs, methods = c("pre", "post"), rotate.multi = TRUE)
kd_cophylo$trees[[2]]$tip.label <- stringi::stri_replace_all_fixed(
  str = kd_cophylo$trees[[2]]$tip.label, pattern = kd_looms$group,
  replacement = kd_looms$lng,
  vectorise_all = FALSE
)
kd_cophylo$trees[[2]] <- kd_cophylo$trees[[2]] |>
  fortify() |>
  left_join(kd_looms, by = join_by(label == lng)) |>
  mutate(label = group) |>
  as.phylo()

kd_lng_tree <- ggtree(kd_cophylo$trees[[1]], ladderize = FALSE, size = lwd, branch.length = "none")
# flip(113, 128) |>
# flip(103,172)
kd_loom_tree <- ggtree(kd_cophylo$trees[[2]], ladderize = FALSE, size = lwd, branch.length = "none")
# flip(1, 33)
#   flip(32, 43) |>
#   rotate(35) |>
#   rotate(54) |>
#   flip(50, 44) |>
#   rotate(50) |>
#   rotate(44) |>
#   rotate(45) |>
#   rotate(55)
# kd_lgs_plt <- tribble(
#   ~lng_group_code, ~loom_type_code,
#   "Be", NA,
#   "Kra", NA,
#   "Hlai", "FBBS",
#   "Ks", "BFYRH",
#   "Tn", NA,
#   "Tsw", "2LH",
#   "Tc", "BFSRH"
# ) |>
#   full_join(distinct(kd_lgs, lng_group_code, lng_group)) |>
#   full_join(distinct(kd_looms, loom_type_code, loom_type)) |>
#   arrange(loom_type, lng_group) |>
#   mutate(lng_col = plt2) |>
#   filter(!is.na(lng_group)) |>
#   mutate(lng_col = fct_inorder(lng_col), lng_group = fct_inorder(lng_group))

kd_lng_tree_data <- kd_lng_tree$data |>
  mutate(lng = label) |>
  left_join(select(kd_lgs, -label)) |>
  left_join(select(kd_looms, group, lng))

kd_loom_tree_data <- kd_loom_tree$data |>
  mutate(group = label) |>
  left_join(kd_looms) |>
  arrange(loom_type)
# mutate(loom_type = str_replace_all(loom_type, "\n", " ")) |>
# mutate(loom_type = str_wrap(loom_type, 20)) |>
# mutate(loom_type = fct_inorder(loom_type))

# ry <- filter(kd_lng_tree_data, !is.na(group) & !is.na(language))$y
ry <- kd_lng_tree_data$y
kd_loom_tree_data$x <- ((kd_loom_tree_data$x - min(kd_loom_tree_data$x)) / (max(kd_loom_tree_data$x) - min(kd_loom_tree_data$x))) * (max(kd_lng_tree_data$x) - min(kd_lng_tree_data$x)) + min(kd_lng_tree_data$x)
kd_loom_tree_data$x <- max(kd_loom_tree_data$x) - kd_loom_tree_data$x + max(kd_lng_tree_data$x)
kd_loom_tree_data$x <- kd_loom_tree_data$x + (max(c(kd_lng_tree_data$x, kd_loom_tree_data$x)) - min(c(kd_lng_tree_data$x, kd_loom_tree_data$x))) / 100 * 70
kd_loom_tree_data$y <- ((kd_loom_tree_data$y - min(kd_loom_tree_data$y)) / (max(kd_loom_tree_data$y) - min(kd_loom_tree_data$y))) * (max(ry) - min(ry)) + min(ry)

kd_lng_loom_tree <- kd_lng_tree +
  # flip(105, 147) +#|> flip(106,132) +
  geom_tree(data = kd_loom_tree_data, linewidth = lwd)
kd_lng_loom_tree_data <- bind_rows(kd_lng_tree_data, kd_loom_tree_data) %>%
  filter(!is.na(group) & !is.na(lng)) |>
  mutate(pb = group %in% kd_loom_pb)

x1 <- max(filter(kd_lng_tree_data, isTip == TRUE)$x) + diff(range(filter(kd_loom_tree_data, isTip == TRUE)$x, filter(kd_lng_tree_data, isTip == TRUE)$x)) / 2.6
x2 <- min(filter(kd_loom_tree_data, isTip == TRUE)$x) - diff(range(filter(kd_loom_tree_data, isTip == TRUE)$x, filter(kd_lng_tree_data, isTip == TRUE)$x)) / 3.2

kd_lgs_looms_links <- kd_lng_loom_tree_data |>
  add_count(lng) |>
  filter(n == 2) |>
  group_by(lng) |>
  mutate(x = ifelse(x == min(x), x1 + .25, x2 - .25))

imgs <- paste0("<img src='", here("data/images/"), levels(kd_loom_tree_data$loom_type_code), ".png' height=60>")

kd_cophylo_plot <- kd_lng_loom_tree +
  geom_rootedge(1, linewidth = lwd) +
  geom_segment(
    data = filter(kd_loom_tree_data, parent == node),
    aes(x = x, xend = x + 1, y = y, yend = y), linewidth = lwd
  ) +
  geom_line(
    data = kd_lgs_looms_links,
    aes(x, y, group = lng, linetype = pb, color = pb)
  ) +
  scale_color_manual(guide = "none", values = c("grey50", "grey70")) +
  ggnewscale::new_scale_colour() +
  # geom_nodelab(aes(label = node), size = 3) +
  geom_tippoint(data = kd_lng_tree_data, aes(color = color, fill = color)) +
  geom_tiplab(data = kd_lng_tree_data, aes(label = label, color = color, fill = color, x = x1), size = (base_font_size - 2) / .pt, family = base_font2, hjust = 1) +
  # geom_tippoint(data = kd_lng_tree_data, aes(
  #   color = color,
  #   fill = color,
  #   shape = !is.na(group),
  #   size = !is.na(group)
  # )) +
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
  ggnewscale::new_scale_colour() +
  ggnewscale::new_scale_fill() +
  # geom_nodelab(data = kd_loom_tree_data, aes(label = node), size = 3) +
  geom_tippoint(data = kd_loom_tree_data, aes(color = color, fill = color)) +
  geom_tiplab(data = filter(kd_loom_tree_data, isTip == TRUE), aes(label = str_replace_all(lng, "(.+?)(?=[A-Z])", "\\1 ") |> str_wrap(10), color = color, x = x2), hjust = 0, size = (base_font_size - 2) / .pt, family = base_font2, lineheight = 1) +
  # geom_tippoint(
  #   data = kd_loom_tree_data,
  #   aes(
  #     color = color,
  #     fill = color,
  #     shape = group %in% kd_lng_tree_data$group,
  #     size = group %in% kd_lng_tree_data$group
  #   )
  # ) +
  scale_color_identity(
    guide = guide_legend(
      order = 2,
      position = "right",
      label.vjust = 0,
      override.aes = list(size = 4),
      theme = theme(
        legend.key.spacing.y = unit(1, "line"),
        # legend.key.height = unit(5, "line"),
        legend.text = element_markdown(),
        legend.margin = margin(0, 0, 0, -.5, unit = "line")
      )
    ),
    name = "Loom type",
    labels = paste0(str_wrap(levels(kd_loom_tree_data$loom_type), 20), "\n", imgs) |> str_replace_all("\\n", "<br/>")
  ) +
  scale_fill_identity(
    guide = "none"
    # name = NULL,
    # guide = guide_legend(
    #   position = "right",
    #   order = 3,
    #   override.aes = list(size = 0, alpha = 0),
    #   theme = theme(
    #     legend.key.spacing.y = unit(1, "line"),
    #     legend.key.height = unit(4, "line"),
    #     legend.key.width = unit(0, "line"),
    #     legend.text = element_markdown(),
    #     legend.margin = margin(2.5, 0, 0, -.5, unit = "line")
    #   )
    # ),
    # labels = imgs
  ) +
  # scale_shape_manual(values = c(21, 23), guide = "none") +
  # scale_size_manual(values = c(.75, 2), guide = "none") +
  scale_linetype(guide = "none") +
  # scale_y_reverse() +
  theme(
    # legend.spacing.x = unit(14, "line"),
    aspect.ratio = 3,
    # legend.position.inside = c(.5, .96),
    # legend.justification.inside = c(.5, 1),
    legend.box = "horizontal",
    # legend.justification.inside = c(.5, 1),
    legend.text = element_text(size = base_font_size - 2, family = base_font2),
    legend.key = element_rect(),
    legend.box.margin = margin(0, 0, 0, 0, unit = "line"),
    legend.background = element_blank()
  )
# kd_cophylo_plot + theme_base()
ggsave(here("output/figures/kd_cophylo_plot.pdf"), kd_cophylo_plot, device = cairo_pdf, width = wd, height = ht * 2, units = "cm")
plot_crop(here("output/figures/kd_cophylo_plot.pdf"))


# Mutation rates --------------------------------------------------------------------------------------------------

kd_looms_mu_bylevel <- read_csv(here("output/data/kd-looms_mu_bylevel.csv"))

kd_looms_mu_plot <- kd_looms_mu_bylevel |>
  ggplot(aes(y = factor(level), x = rate)) +
  stat_density_ridges(aes(fill = .5 - abs(.5 - after_stat(ecdf))), geom = "density_ridges_gradient", calc_ecdf = TRUE, scale = 1, panel_scaling = FALSE, color = "gray50", linewidth = lwd) +
  stat_summary(geom = "text", fun = "median", aes(label = round(..x.., 2)), family = base_font, size = base_font_size / .pt, vjust = 1.5) +
  ylab("Level") +
  xlab("Mutation rate") +
  xlim(0, 2) +
  scale_fill_distiller(palette = "PuBu", direction = 1, limits = c(0, .5), name = "Tail\nprobability") +
  theme_minimal(base_size = base_font_size, base_family = base_font) +
  xtheme +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "line"), aspect.ratio = 0.618, legend.position = "right")
ggsave(here("output/figures/kd_looms_mu_plot.pdf"), kd_looms_mu_plot, device = cairo_pdf, width = wd / 1, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd_looms_mu_plot.pdf"))

kd_looms_mu_summary <- read_csv(here("output/data/kd-looms_mu_summary.csv"))

kd_looms_mu_summary |>
  mutate(across(everything(), ~ round(.x, 2))) |>
  unite(hdi, hdi_lower, hdi_upper, sep = ", ") |>
  mutate(hdi = paste0("[", hdi, "]")) |>
  rename(characters = n_chars, `95% HPDI` = hdi) |>
  kbl(
    digits = 2,
    format = "latex", booktabs = TRUE
  ) |>
  write_lines(here("output/tables/kd-looms_mu_summary.tex"))


# Maps --------------------------------------------------------------------

theme_set(theme_bw() + ytheme)

kd_lgs_pts <- kd_lgs |>
  filter(!is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

kd_looms_pts <- kd_looms |>
  filter(!is.na(lon)) |>
  arrange(loom_type) |>
  # mutate(loom_type = str_replace_all(loom_type, "\\n", " ")) |>
  mutate(loom_type = str_replace_all(loom_type, ", ", ",\n")) |>
  mutate(loom_type = fct_inorder(loom_type)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

kd_bbx_lat <- bind_rows(kd_lgs_pts, kd_looms_pts) |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(30 * 10^3) |>
  st_bbox()

kd_bbx_lon <- bind_rows(kd_lgs_pts, kd_looms_pts) |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(100 * 10^3) |>
  st_bbox()

kd_bbx <- kd_bbx_lon
kd_bbx["ymin"] <- kd_bbx_lat["ymin"]
kd_bbx["ymax"] <- kd_bbx_lat["ymax"]

# prj <- "+proj=cea +lon_0=106 +lat_ts=22 +datum=WGS84 +units=m +no_defs"
# prj <- "+proj=cea +lon_0=103 +lat_ts=17 +datum=WGS84 +units=m +no_defs"
# prj <- "+proj=aea +lon_0=103 +lat_1=9.6666667 +lat_2=24.3333333 +lat_0=17 +datum=WGS84 +units=m +no_defs"
prj <- "+proj=aea +lon_0=102.9741586 +lat_1=9.6945833 +lat_2=25.0240278 +lat_0=17.3593056 +datum=WGS84 +units=m +no_defs"

asia <- ne_countries(scale = "medium") |>
  filter(continent %in% c("Asia", "Oceania")) |>
  st_transform(prj) |>
  st_crop(st_transform(st_as_sfc(kd_bbx), prj))

country_lbs <- asia |>
  select(label_x, label_y, name_en) |>
  filter(!(name_en %in% c("Hong Kong", "Macau", "Malaysia", "Bhutan", "Bangladesh"))) |>
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
    style = north_arrow_orienteering(text_family = base_font, text_size = base_font_size - 2)
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
  geom_sf(data = kd_lgs_pts, aes(color = color), size = 1) +
  scale_color_identity(name = "Language group", guide = guide_legend(override.aes = list(size = 4)), labels = levels(kd_lgs_pts$lng_group)) +
  coord_sf(crs = prj, expand = FALSE)

kd_looms_map <- bg_map +
  geom_sf(data = kd_looms_pts, aes(color = color), size = 1) +
  scale_color_identity(name = "Loom type", guide = guide_legend(override.aes = list(size = 4)), labels = levels(kd_looms_pts$loom_type)) +
  coord_sf(crs = prj, clip = "on", expand = FALSE)

kd_lgs_map <- set_dim(kd_lgs_map, get_dim(kd_looms_map))

ggsave(here("output/figures/kd-lgs_map.pdf"), kd_lgs_map, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-lgs_map.pdf"))

ggsave(here("output/figures/kd-looms_map.pdf"), kd_looms_map, device = cairo_pdf, width = wd, height = wd * 2, units = "cm")
plot_crop(here("output/figures/kd-looms_map.pdf"))
