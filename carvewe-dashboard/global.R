# global.R
# === Core Shiny ===
library(shiny)

# === Data Manipulation ===
library(tidyverse)

# === Visualization ===
library(patchwork)     # Combining plots
library(ggforce)       # ggplot extensions
library(ggtext)        # Rich text in ggplot
library(ggbreak)       # Axis breaks
library(ggridges)      # Ridge plots
library(ggfortify)     # Fortify methods for ggplot
library(gplots)        # Additional plotting (heatmap.2, etc.)
library(ragg)          # Graphics device
library(viridis)       # Specific color palettes

# === Phylogenetic/Tree ===
library(ggtree)        # Tree visualization
library(treeio)        # Tree I/O
library(ggnewscale)    # Multiple color scales in ggplot

# === Analysis ===
library(kohonen)       # SOM (only load once!)
library(vegan)         # Ecological analysis (only load once!)
library(cluster)       # Clustering algorithms
library(mclust)        # Model-based clustering
library(multcompView)  # Multiple comparisons

# === Network/Graph ===
library(igraph)        # Network analysis

# === Utilities ===
library(png)           # PNG graphics

# ============================================
# Set a larger request size for Shiny
# ============================================
options(shiny.maxRequestSize = 100 * 1024^2)
options(verbose = TRUE)

# Source all helpers/scripts
source("carvewe-dashboard/scripts/server_som_clustering.R")
source("carvewe-dashboard/scripts/server_metabolite_analysis.R")

# Source all modules
source("carvewe-dashboard/R/module_upload.R")
source("carvewe-dashboard/R/module_analysis.R")
# source("R/module_results.R")  # When you create it


cat("========================================\n")
cat("GLOBAL.R LOADED SUCCESSFULLY\n")
cat("========================================\n")

