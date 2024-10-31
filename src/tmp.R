library(TreeTools)
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

keep.tip(kd_lgs_cs, kd_lgs_pruned_tips)
keep.tip(kd_looms_cs, kd_lgs_pruned_tips)

kd_looms_cs$tip.label |> sort()
kd_lgs_pruned_tips |> sort()

kd_cophylo_pruned <- cophylo(keep.tip(kd_lgs_cs, kd_lgs_pruned_tips),
  keep.tip(kd_looms_cs, kd_lgs_pruned_tips),
  methods = c("pre", "post"), rotate.multi = FALSE
)

kd_lng_tree_pruned <- ggtree(kd_cophylo_pruned$trees[[1]],
  ladderize = FALSE,
  size = lwd,
  branch.length = "none"
)

kd_loom_tree_pruned <- ggtree(kd_cophylo_pruned$trees[[2]],
  ladderize = FALSE,
  size = lwd,
  branch.length = "none"
)

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
)) / 2.6
x2 <- min(filter(kd_loom_tree_pruned_data, isTip == TRUE)$x) - diff(range(
  filter(kd_loom_tree_pruned_data, isTip == TRUE)$x,
  filter(kd_lng_tree_pruned_data, isTip == TRUE)$x
)) / 3.2

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
  ggnewscale::new_scale_colour() +
  geom_tippoint(data = kd_lng_tree_pruned_data, aes(color = color, fill = color)) +
  geom_tiplab(
    data = kd_lng_tree_pruned_data,
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
    labels = levels(kd_lng_tree_pruned_data$lng_group)
  ) +
  xtheme +
  scale_fill_identity(guide = "none") +
  ggnewscale::new_scale_colour() +
  ggnewscale::new_scale_fill() +
  geom_tippoint(data = kd_loom_tree_pruned_data, aes(color = color, fill = color)) +
  geom_tiplab(
    data = filter(kd_loom_tree_pruned_data, isTip == TRUE),
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
      str_wrap(levels(kd_loom_tree_pruned_data$loom_type), 20),
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


library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(FactoMineR)
library(ggtree)
library(tidyverse)






pca_looms <- ReadAsPhyDat(here("data/nexus/kd-looms_pruned.nex")) |>
  as_tibble() |>
  rename_with(~ str_replace_all(.x, "_", " ")) |>
  rename(any_of(setNames(kd_looms$group, kd_looms$lng_label))) |>
  mutate(across(everything(), ~ as.numeric(.x))) |>
  as.matrix() |>
  t() |>
  PCA(graph = FALSE, scale.unit = FALSE)
pc1_looms <- pca_looms$ind$coord[, 1]

pca_lgs <- ReadAsPhyDat(here("data/nexus/kd-lgs_pruned.nex")) |>
  as_tibble() |>
  mutate(across(everything(), ~ as.numeric(.x))) |>
  as.matrix() |>
  t() |>
  PCA(graph = FALSE, scale.unit = FALSE)
pc1_lgs <- pca_lgs$ind$coord[, 1]

set.seed(30101005)
kd_lgs_on_looms_k <- map_df(1:length(kd_looms_phylo), function(i) {
  ksig <- phylosig(kd_looms_phylo[[i]], pc1_lgs, method = "K", test = TRUE)
  tibble(k = ksig$K, p = ksig$P)
})
kd_lgs_on_looms_k |>
  write_csv(here("output/data/kd-lgs_on_looms_k.csv"))

kd_looms_on_lgs_k <- map_df(1:length(kd_lgs_phylo), function(i) {
  ksig <- phylosig(kd_lgs_phylo[[i]], pc1_looms, method = "K", test = TRUE)
  tibble(k = ksig$K, p = ksig$P)
})
kd_looms_on_lgs_k |>
  write_csv(here("output/data/kd-looms_on_lgs_k.csv"))

kd_lgs_on_looms_k |>
  summarise(across(everything(), mean))
kd_looms_on_lgs_k |>
  summarise(across(everything(), mean))
