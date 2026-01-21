# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Context

The CarveWe Dashboard (`carvewe-dashboard/`) is a Shiny web application within the larger **CarveWe** repository - a bioinformatics workflow for analyzing genome-scale metabolic models (GEMs) of marine microbes. The parent repository includes:
- CLI tools for genome alignment and metabolic sensitivity prediction
- Publication reproduction scripts (`reproduce_publication/`)
- Processed publication data (`Publication_Data/`)
- Shared utility functions (`reproduce_publication/Codes/CarveWe_functions.R`)

## Overview

The CarveWe Dashboard is an R Shiny web application for analyzing metabolic sensitivity data from genome-scale metabolic models (GEMs). It implements a complete workflow for:
1. Data upload and validation
2. Self-Organizing Map (SOM) neural network training
3. Hierarchical clustering of SOM nodes
4. Metabolite-level visualizations and comparative analyses

The dashboard supports both demo data from the CarveWe publication and user-uploaded sensitivity data.

## Running the Dashboard

### From Repository Root
```bash
./launch_dashboard.sh
```

### From Dashboard Directory
```bash
cd carvewe-dashboard
Rscript -e "shiny::runApp(port=3838, launch.browser=TRUE)"
```

### From R Console
```R
setwd("carvewe-dashboard")
shiny::runApp(port=3838, launch.browser=TRUE)
```

## Installing Dependencies

All required R packages are loaded in `global.R`. To install missing packages:
```R
install.packages(c("shiny", "tidyverse", "kohonen", "vegan", "cluster",
                   "mclust", "patchwork", "ggforce", "ggtext", "ggbreak",
                   "ggridges", "ggfortify", "gplots", "ragg", "viridis",
                   "hues", "ggtree", "treeio", "ggnewscale", "igraph",
                   "multcompView", "png"))
```

## Architecture

### Modular Structure

The app follows Shiny's module pattern with three main components:

1. **Upload Module** (`R/module_upload.R`)
   - Handles demo data loading or user file uploads
   - Validates input data format
   - Generates `feature_data` matrix (genomes × metabolites)

2. **Analysis Module** (`R/module_analysis.R`)
   - Trains SOM on feature_data
   - Performs hierarchical clustering on SOM nodes
   - Assigns genomes to clusters using plurality voting across replicates

3. **Metabolite Module** (`R/module_metabolites.R`)
   - Visualizes metabolite patterns across clusters
   - Generates prototype plots, variance analyses, sensitivity profiles

### File Organization

```
carvewe-dashboard/
├── app.R                    # Main entry point (UI + server)
├── global.R                 # Library loading, script sourcing
├── R/                       # Shiny modules
│   ├── module_upload.R
│   ├── module_analysis.R
│   └── module_metabolites.R
├── scripts/                 # Backend functions
│   ├── server_som_clustering.R      # Core SOM/clustering logic
│   ├── server_metabolite_analysis.R # Metabolite visualization functions
│   ├── preprocess_publication_data.R # Data preprocessing pipeline
│   └── visualization_replicates.R   # Helper plotting functions
└── data/                    # Processed data files
    ├── processed_sensitivities_long.csv
    ├── processed_sensitivities_wide.csv
    └── processed_quality_data.csv
```

### Reactive Values System

All modules share state through the `rv` reactiveValues object:

**Data Storage:**
- `rv$user_data` - Raw uploaded/demo sensitivity data
- `rv$feature_data` - Matrix (genomes × metabolites) for SOM training
- `rv$taxonomy`, `rv$quality`, `rv$environment` - Optional contextual data

**Analysis Results:**
- `rv$som_model`, `rv$som_codes`, `rv$som_dist` - SOM training outputs
- `rv$node_clusters` - Cluster assignments for SOM nodes
- `rv$genome_assignments` - Final genome → cluster mappings with plurality voting
- `rv$gridsize` - SOM grid dimensions

**State Flags:**
- `rv$upload_complete`, `rv$som_complete`, `rv$cluster_complete`

## Key Workflows

### SOM Training Pipeline (scripts/server_som_clustering.R)

1. `generate_feature_data(input_data)` - Converts long-form data to wide matrix
2. `train_som(feature_matrix, grid_size)` - Trains kohonen SOM
3. `perform_hierarchical_clustering(som_codes, k)` - Clusters SOM nodes
4. `assign_genomes_by_plurality(som_model, node_clusters)` - Maps genomes to clusters

### Clustering Methodology

The dashboard uses a two-stage clustering approach:
1. **SOM Training**: Maps high-dimensional metabolite data onto 2D grid
2. **Node Clustering**: Performs hierarchical clustering on SOM prototype vectors
3. **Genome Assignment**: Uses plurality voting when multiple model replicates exist per genome

### Metabolite Analysis Functions (scripts/server_metabolite_analysis.R)

- `prepare_scaled_data()` - Normalizes metabolite values for visualization
- `analyze_prototype_distributions()` - Compares cluster prototypes
- `analyze_replicate_variance()` - Quantifies model ensemble variation
- `analyze_high_sensitivity_metabolites()` - Identifies discriminative features
- `analyze_cluster_metabolite_profiles()` - Ridge plots by cluster
- `analyze_metabolite_bubble_plot()` - Polar coordinate flux visualization

## Data Format Requirements

### Required Input (Sensitivities CSV)
Must contain columns:
- `genome_id` - Unique genome identifier
- `nutrient_class` - Metabolite category name
- `sensitivity_score` - Numeric sensitivity value

Optional column for replicate ensembles:
- `model_id` - Replicate identifier (e.g., "genome1_model1")

### Optional Context Files
- **Taxonomy**: `genome_id, phylum, class, genus, species`
- **Quality**: `genome_id, completeness, contamination`
- **Environment**: `genome_id, temperature_category, depth_category, ...`

## Data Preprocessing

To prepare publication data for the dashboard:

```bash
cd carvewe-dashboard
Rscript scripts/preprocess_publication_data.R
```

This script:
- Loads raw flux data from `Publication_Data/`
- Filters genomes by quality criteria
- Calculates sensitivity scores
- Exports processed CSVs to `data/`

## Key Dependencies

Critical R packages loaded in `global.R`:
- **Shiny ecosystem**: `shiny`
- **Data manipulation**: `tidyverse`
- **SOM/clustering**: `kohonen`, `vegan`, `cluster`, `mclust`
- **Visualization**: `ggplot2`, `patchwork`, `ggforce`, `ggridges`, `ggtree`
- **Network analysis**: `igraph`

## Development Notes

### Adding New Visualizations

1. Add analysis function to `scripts/server_metabolite_analysis.R`
2. Create reactive in `R/module_metabolites.R` using `eventReactive(scale_data(), {...})`
3. Add `plotOutput()` to `metabolitesUI()` and corresponding `renderPlot()` in `metabolitesServer()`

### Modifying SOM Parameters

Default SOM configuration in `train_som()`:
- Grid size: 20×20 (user-configurable via UI)
- Topology: hexagonal, toroidal
- Training iterations: 100
- Learning rate: 0.025 → 0.01

### Debugging Reactive Flow

The code includes extensive `cat()` debug statements. Monitor R console output to trace:
- Data loading progress
- SOM training steps
- Clustering iterations
- Reactive value updates

### Path Dependencies and Working Directory

**CRITICAL**: `global.R` contains hardcoded paths that assume the working directory is the repository root:
```R
source("carvewe-dashboard/scripts/server_som_clustering.R")
source("carvewe-dashboard/R/module_upload.R")
```

When developing:
- Set working directory to repository root (`/Users/Ryan/Documents/Levine_Lab/R_Studio/CarveWe/`)
- Do NOT set working directory to the dashboard subdirectory when sourcing `global.R`
- The `launch_dashboard.sh` script handles this correctly by changing to the dashboard directory before running

**Known Issue**: `global.R` line 51 does not source `R/module_metabolites.R`, but the module is still functional because it's loaded by the Shiny app. If adding new modules, ensure they are sourced in `global.R`.

### Shared Utility Functions

The broader CarveWe repository contains shared utility functions in `reproduce_publication/Codes/CarveWe_functions.R`:
- `bootstrap_sum_pi()` - Bootstrap analysis for cluster statistics
- `clusterMeanDist()` - Distance between cluster means
- `ANOVA.test()` - Tukey HSD testing for significance

These are primarily used in publication reproduction but may be useful for dashboard enhancements.
