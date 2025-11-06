# R/server_metabolite_analysis.R
# Metabolite-level analysis and visualization functions
# Assumes feature_data (metabolite_matrix), genome_assignments, and quality_data available

# ============================================
# THEME AND COLORS
# ============================================

custom_theme <- theme(
  plot.title = element_text(hjust = 0.5),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, color = "black"),
  strip.background = element_blank(),
  strip.text = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  legend.key = element_blank()
)

cluster_colors <- c("#4e7b94", "#547e48", "#a6a541", "#bed196", 
                    "#243a56", "#cd5b34", "#60bb68", "#5279ca")

# ============================================
# DATA PREPARATION
# ============================================

#' Prepare scaled metabolite data
#' @param feature_data Matrix with genomes as rows, metabolites as columns
#' @param genome_assignments Data frame with genome_id and cluster columns
#' @return Pivoted, scaled data frame
prepare_scaled_data <- function(feature_data, genome_assignments) {
  
  all_piv_df <- as_tibble(feature_data) %>%
    mutate(
      clusters = genome_assignments$majority_cluster[match(rownames(feature_data), genome_assignments$genome_id)],
      genomes = rownames(feature_data)
    ) %>%
    pivot_longer(-c(clusters, genomes), names_to = "metab_class", values_to = "scaled_flux")
  
  all_piv_df$clusters <- as.factor(all_piv_df$clusters)
  
  scaled_piv_df <- all_piv_df %>%
    group_by(clusters, metab_class) %>%
    group_by(clusters) %>%
    mutate(
      min = min(scaled_flux),
      max = max(scaled_flux),
      polar_size = (scaled_flux) / (abs(max) + abs(min))
    )
  
  scaled_piv_df$metabolite <- factor(scaled_piv_df$metab_class, 
                                     levels = unique(scaled_piv_df$metab_class))
  
  # Adjust long metabolite names
  scaled_piv_df <- scaled_piv_df %>%
    mutate(metab_class = case_when(
      metab_class == "Phospholipids/Fatty Acids/Triglycerides" ~ 
        "Phospholipids/Fatty\nAcids/Triglycerides",
      metab_class == "Nucleobases/Nucleosides/Nucleotides/Derivatives" ~ 
        "Nucleobases/Nucleosides/\nNucleotides/Derivatives",
      .default = metab_class
    ))
  
  return(scaled_piv_df)
}

#' Prepare SOM prototype data
#' @param som_codes Matrix of SOM prototype vectors
#' @param node_clusters Cluster assignment for each SOM node
#' @param grid_size SOM grid dimension
#' @return Pivoted prototype data frame
prepare_prototype_data <- function(som_codes, node_clusters, grid_size) {
  
  classifiers_df <- as.data.frame(som_codes)
  classifiers_df$nodes <- 1:nrow(som_codes)
  classifiers_df$cluster <- node_clusters
  classifiers_df$column_names <- rep(c(1:grid_size, seq(0.5, grid_size - 0.5, by = 1)), 
                                     times = grid_size)[1:grid_size^2]
  classifiers_df$row_names <- rep(seq(0.5, 0.5 * grid_size, by = 0.5), each = grid_size)
  
  pivot_df <- classifiers_df %>%
    pivot_longer(-c(nodes, cluster, row_names, column_names), 
                 names_to = "metab_class", values_to = "scaled_flux")
  
  pivot_df <- pivot_df %>%
    group_by(nodes) %>%
    mutate(
      min = min(scaled_flux),
      max = max(scaled_flux),
      polar_size = (scaled_flux) / (abs(max) + abs(min))
    )
  
  pivot_df$metabolite <- factor(pivot_df$metab_class, levels = unique(pivot_df$metab_class))
  
  return(list(classifiers = classifiers_df, pivot = pivot_df))
}

# ============================================
# ANALYSIS + PLOTTING FUNCTIONS
# ============================================

#' Analyze and plot SOM prototype distributions per cluster
#' @export
analyze_prototype_distributions <- function(feature_data, genome_assignments) {
  
  scaled_data <- prepare_scaled_data(feature_data, genome_assignments)
  
  plot <- ggplot(data = scaled_data, 
                 aes(x = scaled_flux, y = as.factor(clusters), fill = as.factor(clusters))) +
    geom_density_ridges(scale = 0.9) +
    facet_wrap(~metab_class) +
    coord_cartesian(xlim = c(-0.05, 1.05)) +
    scale_fill_manual(values = cluster_colors) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    theme(legend.position = "none") +
    labs(x = "Growth Sensitivity", y = "Cluster") +
    custom_theme
  
  return(list(data = scaled_data, plot = plot))
}

#' Analyze and plot simple SOM grid
#' @export
analyze_som_grid_simple <- function(som_codes, node_clusters, grid_size) {
  
  proto_data <- prepare_prototype_data(som_codes, node_clusters, grid_size)
  classifiers_df <- proto_data$classifiers
  
  num_clusters <- length(unique(node_clusters))
  
  plot <- ggplot() +
    geom_point(data = classifiers_df, 
               aes(x = column_names, y = row_names, fill = as.factor(cluster)), 
               size = 18, pch = 21, stroke = 1) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    labs(x = "", y = "", fill = "Cluster") +
    scale_fill_manual(values = cluster_colors[1:num_clusters], na.value = "gray80") +
    custom_theme
  
  return(list(data = classifiers_df, plot = plot))
}

#' Analyze replicate variance vs consensus
#' @export
analyze_replicate_variance <- function(feature_data, genome_assignments, quality_data) {
  
  # Assume feature_data rows may have replicates, need to group by genome
  all_metab_df <- as_tibble(feature_data) %>%
    mutate(genome = rownames(feature_data)) %>%
    pivot_longer(-genome, names_to = "metabolite_class", values_to = "growth_req") %>%
    left_join(select(genome_assignments, genome_id, majority_cluster), 
              by = c("genome" = "genome_id")) %>%
    filter(!is.na(majority_cluster))
  
  var_df <- all_metab_df %>%
    group_by(genome, metabolite_class) %>%
    summarise(var = var(growth_req), .groups = "drop") %>%
    group_by(genome) %>%
    summarise(cum_var = sum(var, na.rm = TRUE))
  
  var_df <- left_join(var_df, select(quality_data, genome_id, mean_freq), 
                      by = c("genome" = "genome_id"))
  
  plot <- ggplot() +
    geom_point(data = var_df, aes(x = cum_var, y = mean_freq), size = 2, alpha = 0.75) +
    geom_vline(xintercept = 0.1, linetype = 2) +
    theme_bw() +
    custom_theme +
    labs(x = "Cumulative variance across metabolite classes", y = "Consensus")
  
  return(list(data = var_df, plot = plot))
}

#' Analyze high sensitivity metabolites
#' @export
analyze_high_sensitivity_metabolites <- function(feature_data, genome_assignments) {
  
  all_metab_df <- as_tibble(feature_data) %>%
    mutate(genome = rownames(feature_data)) %>%
    pivot_longer(-genome, names_to = "metabolite_class", values_to = "growth_req") %>%
    left_join(select(genome_assignments, genome_id, majority_cluster), 
              by = c("genome" = "genome_id")) %>%
    filter(!is.na(majority_cluster))
  
  high_sens_df <- all_metab_df %>%
    group_by(genome, metabolite_class) %>%
    summarise(mean = mean(growth_req, na.rm = TRUE), .groups = "drop") %>%
    group_by(metabolite_class) %>%
    summarise(perc = sum(mean >= 0.8) / n())
  
  high_sens_df <- high_sens_df %>%
    mutate(metabolite_class = fct_reorder(metabolite_class, perc, .desc = TRUE))
  
  plot <- ggplot(high_sens_df, aes(x = metabolite_class, y = perc, fill = metabolite_class)) +
    geom_bar(stat = "identity") +
    coord_cartesian(expand = FALSE) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "", y = "Proportion of Genomes", fill = "Metabolite Class") +
    theme_bw() +
    custom_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = c(0.699, 0.75),
      legend.text = element_text(size = 14),
      legend.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(ncol = 2, title.position = "top"))
  
  return(list(data = high_sens_df, plot = plot))
}

#' Analyze cluster metabolite profiles (bar plots)
#' @export
analyze_cluster_metabolite_profiles <- function(feature_data, genome_assignments) {
  
  all_metab_df <- as_tibble(feature_data) %>%
    mutate(genome = rownames(feature_data)) %>%
    pivot_longer(-genome, names_to = "metabolite_class", values_to = "growth_req") %>%
    left_join(select(genome_assignments, genome_id, majority_cluster), 
              by = c("genome" = "genome_id")) %>%
    rename(clusters = majority_cluster) %>%
    filter(!is.na(clusters))
  
  group_metab_df <- all_metab_df %>%
    group_by(clusters, metabolite_class) %>%
    summarise(means = mean(growth_req, na.rm = TRUE), 
              sds = sd(growth_req, na.rm = TRUE), 
              .groups = "drop")
  
  group_metab_df$clusters <- as.factor(group_metab_df$clusters)
  
  # Adjust names
  group_metab_df <- group_metab_df %>%
    mutate(metabolite_class = case_when(
      metabolite_class == "Nucleobases/Nucleosides/Nucleotides/Derivatives" ~ 
        "Nucleobases/Nucleosides/\nNucleotides/Derivatives",
      metabolite_class == "Phospholipids/Fatty Acids/Triglycerides" ~ 
        "Phospholipids/Fatty\nAcids/Triglycerides",
      .default = metabolite_class
    ))
  
  metab_order <- group_metab_df %>%
    group_by(metabolite_class) %>%
    summarise(rank = mean(abs(means), na.rm = TRUE)) %>%
    arrange(desc(rank))
  
  group_metab_df$metabolite_class <- factor(group_metab_df$metabolite_class, 
                                            levels = metab_order$metabolite_class)
  
  num_clusters <- length(unique(group_metab_df$clusters))
  
  plot <- ggplot() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_bar(data = group_metab_df, 
             aes(y = means, x = clusters, fill = clusters), 
             stat = "identity") +
    geom_errorbar(data = group_metab_df, 
                  aes(x = clusters, ymin = means - sds, ymax = means + sds), 
                  width = 0.2) +
    scale_fill_manual(values = cluster_colors[1:num_clusters]) +
    facet_wrap(vars(metabolite_class)) +
    labs(x = "", y = "Growth Sensitivity", fill = "Cluster") +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      strip.text = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      strip.background = element_blank(),
      legend.direction = "horizontal",
      legend.position = "bottom",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25)
    ) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2))
  
  return(list(data = group_metab_df, plot = plot))
}

#' Analyze metabolite bubble plot (subset of key metabolites)
#' @param auto_exclude Boolean, if TRUE excludes metabolite classes where no cluster exceeds 0.1 mean sensitivity
#' @export
analyze_metabolite_bubble_plot <- function(feature_data, genome_assignments, 
                                           auto_exclude = TRUE,
                                           threshold = 0.1) {
  
  profile_result <- analyze_cluster_metabolite_profiles(feature_data, genome_assignments)
  group_metab_df <- profile_result$data
  
  # Dynamically identify classes to exclude
  if (auto_exclude) {
    exclude_classes <- group_metab_df %>%
      group_by(metabolite_class) %>%
      summarise(max_sensitivity = max(means, na.rm = TRUE), .groups = "drop") %>%
      filter(max_sensitivity <= threshold) %>%
      pull(metabolite_class)
    
    cat("Auto-excluding", length(exclude_classes), "metabolite classes with max sensitivity <=", threshold, "\n")
  } else {
    exclude_classes <- character(0)  # Empty vector, exclude nothing
  }
  
  # Scale within each metabolite class
  new_df <- group_metab_df %>%
    group_by(metabolite_class) %>%
    mutate(scaled_flux = ((means + sds) - min(means - sds)) / 
             (max(means + sds) - min(means - sds)))
  
  # Subset
  sub_df <- new_df %>%
    filter(!metabolite_class %in% exclude_classes)
  
  num_clusters <- length(unique(sub_df$clusters))
  
  plot <- ggplot() +
    geom_vline(xintercept = factor(1:num_clusters), lty = 2, linewidth = 0.25) +
    geom_point(data = sub_df, 
               aes(x = clusters, y = metabolite_class, fill = means, size = scaled_flux), 
               pch = 21) +
    scale_fill_viridis_c(option = "plasma") +
    scale_radius(range = c(5, 25), 
                 limits = c(min(sub_df$means), 1), 
                 breaks = c(min(sub_df$means), 0.5, 1), 
                 labels = c("0", "0.5", "1"), 
                 guide = "legend") +
    labs(x = "Cluster", fill = "Growth Sensitivity", size = "Relative Sensitivity") +
    scale_y_discrete(limits = rev) +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(face = "bold", size = 20),
      axis.title.x = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key = element_blank()
    ) +
    guides(
      fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                            barwidth = unit(3.5, "in"), ticks.colour = "black", 
                            frame.colour = "black", order = 1),
      size = guide_legend(title.position = "top", title.hjust = 0.5, order = 2)
    )
  
  return(list(data = sub_df, plot = plot, excluded_classes = exclude_classes))
}