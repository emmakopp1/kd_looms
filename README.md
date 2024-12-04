# Read me

## Description

This is the repository containing the data and code for the article Contrasting modes of cultural evolution: Kra-Dai languages and weaving technologies, by Christopher D. Buckley, Emma Kopp, Thomas Pellard, Robin J. Ryder, and Guillaume Jacques.

The repository has the following structure:

```         
.
├── data # Primary data
│   ├── images # Illusrative drawing of the looms
│   │   ├── BFcant.png
│   │   ├── BFcant.svg
│   │   ├── BFSRH.png
│   │   ├── BFSRH.svg
│   │   ├── BFYRH.png
│   │   ├── BFYRH.svg
│   │   ├── FBBS.png
│   │   ├── FBBS.svg
│   │   ├── FCBcant.png
│   │   ├── FCBcant.svg
│   │   ├── FCB.png
│   │   ├── FCB.svg
│   │   ├── FCBunique.png
│   │   └── FCBunique.svg
│   ├── kd-lgs # Linguistic data
│   │   ├── kd-lgs_bcov_relaxed_ht_pos
│   │   │   ├── kd-lgs_bcov_relaxed_ht_pos.log
│   │   │   ├── kd-lgs_bcov_relaxed_ht_pos.trees
│   │   │   └── kd-lgs_bcov_relaxed_ht_pos.xml
│   │   ├── kd-lgs_bcov_relaxed_uni
│   │   │   ├── kd-lgs_bcov_relaxed_uni.log
│   │   │   ├── kd-lgs_bcov_relaxed_uni.trees
│   │   │   └── kd-lgs_bcov_relaxed_uni.xml
│   │   ├── kd-lgs_bcov_strict_ht_pos
│   │   │   ├── kd-lgs_bcov_strict_ht_pos.log
│   │   │   ├── kd-lgs_bcov_strict_ht_pos.trees
│   │   │   └── kd-lgs_bcov_strict_ht_pos.xml
│   │   ├── kd-lgs_bcov_strict_uni
│   │   │   ├── kd-lgs_bcov_strict_uni.log
│   │   │   ├── kd-lgs_bcov_strict_uni.trees
│   │   │   └── kd-lgs_bcov_strict_uni.xml
│   │   ├── kd-lgs_datapoints.csv # Geographical coordinates of the languages
│   │   └── kd-lgs_lx.csv # Matrix of lexical data
│   ├── kd-looms # Looms data
│   │   ├── kd-looms_bcov1000_strict_uni
│   │   │   ├── kd-looms_bcov1000_strict_uni.log
│   │   │   ├── kd-looms_bcov1000_strict_uni.trees
│   │   │   └── kd-looms_bcov1000_strict_uni.xml
│   │   ├── kd-looms_bcov1111_relaxed_ht
│   │   │   ├── kd-looms_bcov1111_relaxed_ht.log
│   │   │   ├── kd_looms_bcov1111_relaxed_ht.trees
│   │   │   └── kd-looms_bcov1111_relaxed_ht.xml
│   │   ├── kd-looms_bcov1111_relaxed_uni
│   │   │   ├── kd-looms_bcov1111_relaxed_uni.log
│   │   │   ├── kd-looms_bcov1111_relaxed_uni.trees
│   │   │   └── kd-looms_bcov1111_relaxed_uni.xml
│   │   ├── kd-looms_bcov1111_strict_ht
│   │   │   ├── kd-looms_bcov1111_strict_ht.log
│   │   │   ├── kd-looms_bcov1111_strict_ht.out
│   │   │   ├── kd-looms_bcov1111_strict_ht.trees
│   │   │   └── kd-looms_bcov1111_strict_ht.xml
│   │   ├── kd-looms_bcov1111_strict_uni
│   │   │   ├── kd-looms_bcov1111_strict_uni.log
│   │   │   ├── kd-looms_bcov1111_strict_uni.trees
│   │   │   └── kd-looms_bcov1111_strict_uni.xml
│   │   ├── kd-looms_bcov8421_strict_ht
│   │   │   ├── kd-looms_bcov8421_strict_ht.log
│   │   │   ├── kd-looms_bcov8421_strict_ht.trees
│   │   │   └── kd-looms_bcov8421_strict_ht.xml
│   │   ├── kd-looms_bcov8421_strict_uni
│   │   │   ├── kd-looms_bcov8421_strict_uni.log
│   │   │   ├── kd-looms_bcov8421_strict_uni.trees
│   │   │   └── kd-looms_bcov8421_strict_uni.xml
│   │   ├── kd-looms_bcov_basic_strict_ht
│   │   │   ├── kd-looms_bcov_basic_strict_ht.log
│   │   │   ├── kd-looms_bcov_basic_strict_ht.trees
│   │   │   └── kd-looms_bcov_basic_strict_ht.xml
│   │   ├── kd-looms_bcov_basic_strict_uni
│   │   │   ├── kd-looms_bcov_basic_strict_uni.log
│   │   │   ├── kd-looms_bcov_basic_strict_uni.trees
│   │   │   └── kd-looms_bcov_basic_strict_uni.xml
│   │   ├── kd-looms_bcov_patterns_strict_ht
│   │   │   ├── kd-looms_bcov_patterns_strict_ht.log
│   │   │   ├── kd-looms_bcov_patterns_strict_ht.trees
│   │   │   └── kd-looms_bcov_pattern_strict_ht.xml
│   │   ├── kd-looms_bcov_patterns_strict_uni
│   │   │   ├── kd-looms_bcov_patterns_strict_uni.log
│   │   │   ├── kd-looms_bcov_patterns_strict_uni.trees
│   │   │   └── kd-looms_bcov_patterns_strict_uni.xml
│   │   ├── kd-looms_ctmc1111_strict_ht
│   │   │   ├── kd-looms_ctmc1111_strict_ht.log
│   │   │   ├── kd-looms_ctmc1111_strict_ht.trees
│   │   │   └── kd-looms_ctmc1111_strict_ht.xml
│   │   ├── kd-looms_ctmc1111_strict_uni
│   │   │   ├── kd-looms_ctmc1111_strict_uni.log
│   │   │   ├── kd-looms_ctmc1111_strict_uni.trees
│   │   │   └── kd-looms_ctmc1111_strict_uni.xml
│   │   ├── model_choice
│   │   │   ├── kd-looms_bcov1111_relaxed_ht_ns
│   │   │   │   ├── kd-looms_bcov1111_relaxed_ht_ns.log
│   │   │   │   ├── kd_looms_bcov1111_relaxed_ht_ns.trees
│   │   │   │   ├── kd-looms_bcov1111_relaxed_ht_ns.xml
│   │   │   │   └── kd-looms_bcov111_relaxed_ht_ns.out
│   │   │   ├── kd-looms_bcov1111_relaxed_uni_ns
│   │   │   │   ├── kd-looms_bcov1111_relaxed_uni_ns.log
│   │   │   │   ├── kd-looms_bcov1111_relaxed_uni_ns.out
│   │   │   │   ├── kd-looms_bcov1111_relaxed_uni_ns.posterior.log
│   │   │   │   ├── kd-looms_bcov1111_relaxed_uni_ns.posterior.trees
│   │   │   │   ├── kd-looms_bcov1111_relaxed_uni_ns.trees
│   │   │   │   └── kd-looms_bcov1111_relaxed_uni_ns.xml
│   │   │   ├── kd-looms_bcov1111_strict_ht_ns
│   │   │   │   ├── kd-looms_bcov1111_strict_ht_ns.log
│   │   │   │   ├── kd-looms_bcov1111_strict_ht_ns.out
│   │   │   │   ├── kd-looms_bcov1111_strict_ht_ns.posterior.log
│   │   │   │   ├── kd-looms_bcov1111_strict_ht_ns.posterior.trees
│   │   │   │   ├── kd-looms_bcov1111_strict_ht_ns.trees
│   │   │   │   └── kd-looms_bcov1111_strict_ht_ns.xml
│   │   │   ├── kd-looms_bcov1111_strict_uni_ns
│   │   │   │   ├── kd-looms_bcov1111_strict_uni_ns.log
│   │   │   │   ├── kd-looms_bcov1111_strict_uni_ns.out
│   │   │   │   ├── kd-looms_bcov1111_strict_uni_ns.posterior.log
│   │   │   │   ├── kd-looms_bcov1111_strict_uni_ns.posterior.trees
│   │   │   │   ├── kd-looms_bcov1111_strict_uni_ns.trees
│   │   │   │   └── kd-looms_bcov1111_strict_uni_ns.xml
│   │   │   ├── kd-looms_ctmc1111_strict_ht_ns
│   │   │   │   ├── kd-looms_ctmc1111_strict_ht_ns.log
│   │   │   │   ├── kd-looms_ctmc1111_strict_ht_ns.out
│   │   │   │   ├── kd-looms_ctmc1111_strict_ht_ns.posterior.log
│   │   │   │   ├── kd-looms_ctmc1111_strict_ht_ns.posterior.trees
│   │   │   │   ├── kd-looms_ctmc1111_strict_ht_ns.trees
│   │   │   │   └── kd-looms_ctmc1111_strict_ht_ns.xml
│   │   │   └── kd-looms_ctmc1111_strict_uni_ns
│   │   │       ├── kd-looms_ctmc1111_strict_uni_ns.log
│   │   │       ├── kd-looms_ctmc1111_strict_uni_ns.out
│   │   │       ├── kd-looms_ctmc1111_strict_uni_ns.posterior.log
│   │   │       ├── kd-looms_ctmc1111_strict_uni_ns.posterior.trees
│   │   │       ├── kd-looms_ctmc1111_strict_uni_ns.trees
│   │   │       └── kd-looms_ctmc1111_strict_uni_ns.xml
│   │   ├── kd-looms_characters.csv # Information about the looms characters
│   │   ├── kd-looms_data.ods # Raw data and sources for the looms
│   │   ├── kd-looms_datapoints.csv # Geographical coordinates of the looms
│   │   └── kd-looms_matrix.csv # Raw character-state matrix for the looms
│   ├── kd-pruned
│   │   ├── kd-lgs_pruned_bcov_strict_bd
│   │   │   ├── kd-lgs_pruned_bcov_strict_bd.log
│   │   │   ├── kd-lgs_pruned_bcov_strict_bd.trees
│   │   │   └── kd-lgs_pruned_bcov_strict_bd.xml
│   │   ├── kd-lgs_pruned_bcov_strict_bd_ns
│   │   │   ├── kd-lgs_pruned_bcov_strict_bd_ns.log
│   │   │   ├── kd-lgs_pruned_bcov_strict_bd_ns.out
│   │   │   ├── kd-lgs_pruned_bcov_strict_bd_ns.trees
│   │   │   └── kd-lgs_pruned_bcov_strict_bd_ns.xml
│   │   ├── kd-looms_pruned_bcov_strict_bd
│   │   │   ├── kd-looms_pruned_bcov_strict_bd.log
│   │   │   ├── kd-looms_pruned_bcov_strict_bd.trees
│   │   │   └── kd-looms_pruned_bcov_strict_bd.xml
│   │   ├── kd-looms_pruned_bcov_strict_bd_ns
│   │   │   ├── kd-looms_pruned_bcov_strict_bd.log
│   │   │   ├── kd-looms_pruned_bcov_strict_bd_ns.out
│   │   │   ├── kd-looms_pruned_bcov_strict_bd_ns.xml
│   │   │   └── kd-looms_pruned_bcov_strict_bd.trees
│   │   ├── kd-merged_pruned_bcov_strict_bd
│   │   │   ├── kd-merged_pruned_bcov_strict_bd.log
│   │   │   ├── kd-merged_pruned_bcov_strict_bd.trees
│   │   │   └── kd-merged_pruned_bcov_strict_bd.xml
│   │   └── kd-merged_pruned_bcov_strict_bd_ns
│   │       ├── kd-merged_pruned_bcov_strict_bd_ns.log
│   │       ├── kd-merged_pruned_bcov_strict_bd_ns.out
│   │       ├── kd-merged_pruned_bcov_strict_bd_ns.trees
│   │       └── kd-merged_pruned_bcov_strict_bd_ns.xml
│   ├── nexus # Character state matrices in NEXUS format
│   │   ├── kd-lgs.nex
│   │   ├── kd-lgs_pos.nex
│   │   ├── kd-lgs_pruned.nex
│   │   ├── kd-lngs_pruned_filtered.nex
│   │   ├── kd-looms_1000.nex
│   │   ├── kd-looms_1111_ht.nex
│   │   ├── kd-looms_1111.nex
│   │   ├── kd-looms_8421_ht.nex
│   │   ├── kd-looms_8421.nex
│   │   ├── kd-looms_basic_ht.nex
│   │   ├── kd-looms_basic.nex
│   │   ├── kd-looms_patterns_ht.nex
│   │   ├── kd-looms_patterns.nex
│   │   ├── kd-looms_pruned_filtered.nex
│   │   ├── kd-looms_pruned.nex
│   │   ├── kd-merged_filtered.nex
│   │   └── kd-merged.nex
│   ├── kd_looms_languages.csv # Correspondence table of the languages and the looms
│   └── kd_looms_matrix_levels.csv # Character state matrix for looms with level information
├── output # Output produced by running the R code files in src/
│   ├── data # Data summaries
│   │   ├── kd-k_summary.csv
│   │   ├── kd-lgs_clade_ages.csv
│   │   ├── kd-lgs_clade_ages_summary.csv
│   │   ├── kd-lgs_mu_pos.csv
│   │   ├── kd-lgs_mu_pos_summary.csv
│   │   ├── kd-lgs_on_looms_k.csv
│   │   ├── kd-looms_mu_bylevel.csv
│   │   ├── kd-looms_mu_summary.csv
│   │   ├── kd-looms_on_lgs_k.csv
│   │   └── models_summary.csv
│   ├── figures
│   │   ├── kd_cophylo_pruned_plot.pdf
│   │   ├── kd-lgs_ages_plot.pdf
│   │   ├── kd-lgs_bcov_relaxed_ht_pos_cs_tree.pdf
│   │   ├── kd-lgs_bcov_relaxed_uni_cs_tree.pdf
│   │   ├── kd-lgs_bcov_strict_ht_pos_cs_tree.pdf
│   │   ├── kd-lgs_bcov_strict_uni_cs_tree.pdf
│   │   ├── kd-lgs_map.pdf
│   │   ├── kd-lgs_mu_plot.pdf
│   │   ├── kd-looms_bcov1000_strict_uni_cs_tree.pdf
│   │   ├── kd-looms_bcov1111_relaxed_ht_cs_tree.pdf
│   │   ├── kd-looms_bcov1111_relaxed_uni_cs_tree.pdf
│   │   ├── kd-looms_bcov1111_strict_ht_cs_tree.pdf
│   │   ├── kd-looms_bcov1111_strict_uni_cs_tree.pdf
│   │   ├── kd-looms_bcov8421_strict_ht_cs_tree.pdf
│   │   ├── kd-looms_bcov8421_strict_uni_cs_tree.pdf
│   │   ├── kd-looms_bcov_basic_strict_ht_cs_tree.pdf
│   │   ├── kd-looms_bcov_basic_strict_uni_cs_tree.pdf
│   │   ├── kd-looms_bcov_patterns_strict_ht_cs_tree.pdf
│   │   ├── kd-looms_bcov_patterns_strict_uni_cs_tree.pdf
│   │   ├── kd-looms_ctmc1111_strict_ht_cs_tree.pdf
│   │   ├── kd-looms_ctmc1111_strict_uni_cs_tree.pdf
│   │   ├── kd-looms_map.pdf
│   │   ├── kd-looms_mu_plot.pdf
│   │   └── kd-looms_Tsw_map.pdf
│   ├── tables # Summary tables in LaTeX format
│   │   ├── kd-lgs_ages_summary.tex
│   │   ├── kd-lgs_mu_summary.tex
│   │   ├── kd-looms_models_summary.tex
│   │   └── kd-looms_mu_summary.tex
│   └── trees # Consensus trees
│       ├── kd-lgs_bcov_relaxed_ht_pos_consensus.tree
│       ├── kd-lgs_bcov_relaxed_uni_consensus.tree
│       ├── kd-lgs_bcov_strict_ht_pos_consensus.tree
│       ├── kd-lgs_bcov_strict_uni_consensus.tree
│       ├── kd-lgs_prunedk.trees
│       ├── kd-looms_bcov1000_strict_uni_consensus.tree
│       ├── kd-looms_bcov1111_relaxed_ht_consensus.tree
│       ├── kd-looms_bcov1111_relaxed_uni_consensus.tree
│       ├── kd-looms_bcov1111_strict_ht_consensus.tree
│       ├── kd-looms_bcov1111_strict_uni_consensus.tree
│       ├── kd-looms_bcov8421_strict_ht_consensus.tree
│       ├── kd-looms_bcov8421_strict_uni_consensus.tree
│       ├── kd-looms_bcov_basic_strict_ht_consensus.tree
│       ├── kd-looms_bcov_basic_strict_uni_consensus.tree
│       ├── kd-looms_bcov_patterns_strict_ht_consensus.tree
│       ├── kd-looms_bcov_patterns_strict_uni_consensus.tree
│       ├── kd-looms_ctmc1111_strict_ht_consensus.tree
│       ├── kd-looms_ctmc1111_strict_uni_consensus.tree
│       └── kd-looms_prunedk.trees
├── src # R code files
│   ├── 00_Preprocess.R
│   ├── 01_Transform.R
│   ├── 02_Analyse.R
│   └── 03_Visualise.R
├── kd-looms.Rproj # RStudio project file
└── README.md
```

## Running the project

1.  **Preprocess data**
    -   Run `00_Preprocess.R` to create the NEXUS files for languages and looms data, and to adjust and check the BEAST outputs.
2.  **Transform data**
    -   Run `01_Transform.R` to generate consensus trees, compute internal nodes depth, mutation rates, and marginal likelihood when needed.
3.  **Analyse data**
    -   Run `02_Analyse.R` to compute Blomberg's K.
4.  **Visualise data**
    -   Run `03_Visualise.R` to generate the figures and tables used in the article.

## Requirements

The R packages listed below are needed to run the analysis, and they can be installed with the following command in R (<https://www.r-project.org/>):

``` r
install.packages(c(
  "FactoMineR",
  "HDInterval",
  "TreeTools",
  "fs",
  "ggforce",
  "ggnewscale",
  "ggridges",
  "ggspatial",
  "ggstar",
  "ggtext",
  "ggthemes",
  "glue",
  "here",
  "kableExtra",
  "knitr",
  "patchwork",
  "phangorn",
  "phytools",
  "rnaturalearth",
  "sf",
  "stringi",
  "tidyverse",
  "tracerer",
  "utils"
))

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ggtree")
```

The Noto Sans Condensed, Noto Sans SemiCondensed, and Noto Sans ExtraCondensed free fonts are also needed to produce the figures. They are available at <https://notofonts.github.io/latin-greek-cyrillic/>.

Phylogenetic analyses require the BEAST software v2.6.7 (<https://www.beast2.org/>).
The additional packages Babel, CoupledMCMC, and SA for BEAST are also needed.

## License

This dataset and code is licensed under a CC-BY-4.0 license (<https://creativecommons.org/licenses/by/4.0/>).
