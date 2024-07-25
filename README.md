## Structure of the repository

This is the repository containing the data and code for the article Contrasting modes of cultural evolution: Kra-Dai languages and weaving technologies, by Christopher D. Buckley, Emma Kopp, Thomas Pellard, Robin J. Ryder, and Guillaume Jacques.

The repository has the following structure:

```         
.
├── kd-looms.Rproj                                 # RStudio project file
├── README.md
├── data                                          # Primary data
│   ├── Kra-DaiLooms28master-2.xlsx
│   ├── images                                     # Illusrative drawing of the looms
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
│   ├── kd-lgs                                     # Linguistic data
│   │   ├── kd-lgs_datapoints.csv                  # Geographical coordinates of the languages
│   │   ├── kd-lgs_lx.csv                          # Matrix of lexical data
│   │   └── kd-lgs_bcov                            # BEAST input and output files
│   │       ├── kd-lgs_bcov.log
│   │       ├── kd-lgs_bcov.trees.zip
│   │       ├── kd-lgs_bcov.xml
│   │       └── kd-lgs.nex
│   └── kd-looms                                   # Looms data
│       ├── kd-looms_datapoints.csv                # Geographical coordinates of the looms
│       ├── kd-looms_characters.csv                # Information about the looms characters
│       ├── kd-looms_matrix.csv                    # Raw character-state matrix for the looms
│       ├── kd-looms_bcov1000                      # BEAST input and output files
│       │   ├── kd-looms_1000.nex
│       │   ├── kd-looms_bcov1000.log
│       │   ├── kd-looms_bcov1000.trees
│       │   └── kd-looms_bcov1000.xml
│       ├── kd-looms_bcov111 1                     # BEAST input and output files
│       │   ├── kd-looms_1111.nex
│       │   ├── kd-looms_bcov1111.log
│       │   ├── kd-looms_bcov1111.trees
│       │   └── kd-looms_bcov1111.xml
│       ├── kd-looms_bcov8421                      # BEAST input and output files
│       │   ├── kd-looms_8421.nex
│       │   ├── kd-looms_bcov8421.log
│       │   ├── kd-looms_bcov8421.trees
│       │   └── kd-looms_bcov8421.xml
│       ├── kd-looms_bcov_basic                    # BEAST input and output files
│       │   ├── kd-looms_basic.nex
│       │   ├── kd-looms_bcov_basic.log
│       │   ├── kd-looms_bcov_basic.trees
│       │   └── kd-looms_bcov_basic.xml
│       ├── kd-looms_bcov_patterns                 # BEAST input and output files
│       │   ├── kd-looms_bcov_patterns.log
│       │   ├── kd-looms_bcov_patterns.trees
│       │   ├── kd-looms_bcov_patterns.xml
│       │   └── kd-looms_patterns.nex
│       └── kd-looms_ctmc4                         # BEAST input and output files
│           ├── kd-looms_ctmc4.log
│           ├── kd-looms_ctmc4.nex
│           ├── kd-looms_ctmc4.trees
│           └── kd-looms_ctmc4.xml
├── output                                         # Output produced by running the R code files in src/
│   ├── data                                       # Data summaries
│   │   ├── kd-lgs_ages.csv
│   │   ├── kd-lgs_ages_summary.csv
│   │   ├── kd-lgs_mu_summary.csv
│   │   ├── kd-looms_mu_bylevel.csv
│   │   └── kd-looms_mu_summary.csv
│   ├── figures                                    # Figures for the analysis results
│   │   ├── kd_cophylo_plot.pdf
│   │   ├── kd-lgs_ages_plot.pdf
│   │   ├── kd-lgs_bcov_byconcept_cs_tree.pdf
│   │   ├── kd-lgs_bcov_cs_tree.pdf
│   │   ├── kd-lgs_map.pdf
│   │   ├── kd-lgs_mu_plot.pdf
│   │   ├── kd-looms_bcov1000_cs_tree.pdf
│   │   ├── kd-looms_bcov1111_cs_tree.pdf
│   │   ├── kd-looms_bcov8421_cs_tree.pdf
│   │   ├── kd-looms_bcov_basic_cs_tree.pdf
│   │   ├── kd-looms_bcov_patterns_cs_tree.pdf
│   │   ├── kd-looms_ctmc4_cs_tree.pdf
│   │   ├── kd-looms_map.pdf
│   │   ├── kd-looms_mu_plot.pdf
│   │   └── kd-looms_Tsw_map.pdf
│   ├── tables                                     # Summary tables in LaTeX format
│   │   ├── kd-lgs_ages_summary.tex
│   │   └── kd-looms_mu_summary.tex
│   └── trees                                      # Consensus trees
│       ├── kd-lgs_bcov_byconcept_consensus.tree
│       ├── kd-lgs_bcov_consensus.tree
│       ├── kd-looms_bcov1000_consensus.tree
│       ├── kd-looms_bcov1111_consensus.tree
│       ├── kd-looms_bcov8421_consensus.tree
│       ├── kd-looms_bcov_basic_consensus.tree
│       ├── kd-looms_bcov_patterns_consensus.tree
│       └── kd-looms_ctmc4_consensus.tree
└── src                                            # R code files
    ├── 00_Preprocess.R
    ├── 01_Transform.R
    └── 02_Visualise.R
```

The repository `src` contains three files. `00_Preprocess.R` have to objective. It creates the nexus files from the csv files for language and looms data. In particular, for looms data, it weights the according to the differents choosen levels. It also check that traces and ESS outputs after the beast analysis are valids. `01_Transform.R` generates the consensus trees, ages of internal nodes (MRCA's) and mutation rates for analysis where we allow multiples evolution rates. The file `02_Visualize.R` generates the figures in `output`.

## Run the project

1.  **Preprocess data**
    -   Run `00_Preprocess.R` to create the NEXUS files for languages and looms data, and to adjust and check the BEAST outputs.
2.  **Transform data**
    -   Run `01_Transform.R` to generate consensus trees, compute internal nodes depth and mutation rates when needed.
3.  **Visualise data**
    -   Run `02_Visualise.R` to generate the figures and tables.

## Dependencies

The R packages listed below are needed to run the analysis, and they can be installed with the following command in R:

``` r
install.packages(c(
    "phangorn",
    "phytools",
    "TreeTools",
    "tracerer",
    "tidyverse",
    "HDInterval",
    "sf",
    "rnaturalearth",
    "ggspatial",
    "stringi",
    "ggtree",
    "ggthemes",
    "ggridges",
    "ggtext",
    "ggstar",
    "ggforce",
    "knitr",
    "kableExtra",
    "patchwork"))
```

The Noto Sans Condensed, Noto Sans SemiCondensed, and Noto Sans ExtraCondensed free fonts are also needed to produce the figures. They are available at <https://notofonts.github.io/latin-greek-cyrillic/>.
