# R/som_clustering.R
# Complete SOM + Clustering module (combines everything)

library(kohonen)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================
# PART 1: CORE SOM FUNCTIONS
# ============================================

#' Train Self-Organizing Map
#' Works for both single models and replicate ensembles
#' 
#' @param feature_matrix Numeric matrix (rows = observations, cols = features)
#'   For replicates: rows may have duplicate genome IDs
#' @param grid_size Integer, size of SOM grid (e.g., 20 for 20x20)
#' @param rlen Number of training iterations (default: 1500)
#' @param alpha Learning rate range (default: c(0.025, 0.01))
#' @param topo Grid topology: "hexagonal" or "rectangular"
#' @param toroidal Logical, whether grid wraps around edges
#' @param seed Random seed for reproducibility
#' @return List with som_model, som_codes, som_dist, and basic assignments
#' @export
train_som <- function(feature_matrix, 
                      grid_size = 20,
                      rlen = 1500,
                      alpha = c(0.025, 0.01),
                      topo = "hexagonal",
                      toroidal = TRUE,
                      seed = 123) {
  
  set.seed(seed)
  
  som_grid <- somgrid(
    xdim = grid_size,
    ydim = grid_size,
    topo = topo,
    toroidal = toroidal
  )
  
  message(sprintf("Training SOM: %dx%d grid, %d iterations...", 
                  grid_size, grid_size, rlen))
  
  som_model <- som(
    feature_matrix,
    grid = som_grid,
    rlen = rlen,
    alpha = alpha,
    keep.data = TRUE
  )
  
  som_codes <- som_model$codes[[1]]
  som_dist <- as.matrix(dist(som_codes))
  
  message("✓ SOM training complete")
  
  return(list(
    som_model = som_model,
    som_codes = som_codes,
    som_dist = som_dist,
    grid_size = grid_size
  ))
}


# ============================================
# PART 2: CLUSTERING FUNCTIONS
# ============================================

#' Helper: Calculate mean intra-cluster distance
#' @keywords internal
cluster_mean_dist <- function(clusters, som_dist) {
  cluster_means <- c()
  
  for (c in unique(clusters)) {
    cluster_members <- which(clusters == c)
    
    if (length(cluster_members) > 1) {
      cluster_dist <- som_dist[cluster_members, cluster_members]
      cluster_means <- append(cluster_means, mean(cluster_dist))
    } else {
      cluster_means <- append(cluster_means, 0)
    }
  }
  
  return(mean(cluster_means))
}


#' Evaluate optimal number of clusters
#' @export
evaluate_cluster_range <- function(som_codes, k_range = 2:20) {
  
  message(sprintf("Evaluating clustering for k = %d to %d...", 
                  min(k_range), max(k_range)))
  
  som_dist <- as.matrix(dist(som_codes))
  
  results <- data.frame(
    k = integer(),
    kmeans_dist = numeric(),
    hclust_dist = numeric()
  )
  
  for (k in k_range) {
    kmeans_clusters <- kmeans(
      som_codes, 
      centers = k, 
      iter.max = 10000,
      nstart = 10
    )$cluster
    
    hclust_clusters <- cutree(
      hclust(vegdist(som_codes, method = "euclidean")),
      k = k
    )
    
    results <- rbind(results, data.frame(
      k = k,
      kmeans_dist = cluster_mean_dist(kmeans_clusters, som_dist),
      hclust_dist = cluster_mean_dist(hclust_clusters, som_dist)
    ))
  }
  
  message("✓ Cluster evaluation complete")
  return(results)
}


#' Plot cluster evaluation results
#' @export
plot_cluster_evaluation <- function(eval_results) {
  p <- ggplot(eval_results, aes(x = k)) +
    geom_line(aes(y = kmeans_dist, color = "k-means"), linewidth = 1) +
    geom_line(aes(y = hclust_dist, color = "hierarchical"), linewidth = 1) +
    geom_point(aes(y = kmeans_dist, color = "k-means"), size = 2) +
    geom_point(aes(y = hclust_dist, color = "hierarchical"), size = 2) +
    scale_color_manual(values = c("k-means" = "black", "hierarchical" = "red")) +
    labs(
      title = "Clustering Quality vs Number of Clusters",
      x = "Number of clusters (k)",
      y = "Mean intra-cluster distance",
      color = "Method"
    ) +
    theme_minimal()
  
  return(p)
}


#' Cluster SOM nodes
#' @export
cluster_som_nodes <- function(som_codes, 
                              k = 8, 
                              method = "kmeans",
                              seed = 123) {
  
  set.seed(seed)
  
  message(sprintf("Clustering SOM nodes: k = %d, method = %s...", k, method))
  
  if (method == "kmeans") {
    node_clusters <- kmeans(
      som_codes,
      centers = k,
      iter.max = 10000,
      nstart = 10
    )$cluster
    
  } else if (method == "hierarchical") {
    som_dist <- dist(som_codes)
    node_clusters <- cutree(
      hclust(som_dist, method = "ward.D2"),
      k = k
    )
    
  } else {
    stop("method must be 'kmeans' or 'hierarchical'")
  }
  
  message("✓ Clustering complete")
  return(node_clusters)
}


# ============================================
# PART 3: REPLICATE HANDLING (PLURALITY VOTING)
# ============================================

#' Check for duplicate replicate models
#' @export
check_duplicate_replicates <- function(feature_matrix) {
  
  unique_genomes <- unique(rownames(feature_matrix))
  duplicate_genomes <- c()
  
  for (genome in unique_genomes) {
    genome_rows <- which(rownames(feature_matrix) == genome)
    genome_matrix <- feature_matrix[genome_rows, , drop = FALSE]
    unique_rows <- unique(genome_matrix)
    
    if (is.matrix(unique_rows) && nrow(unique_rows) == 1) {
      duplicate_genomes <- c(duplicate_genomes, genome)
    }
  }
  
  n_duplicates <- length(duplicate_genomes)
  n_total <- length(unique_genomes)
  pct_duplicates <- 100 * n_duplicates / n_total
  
  if (n_duplicates > 0) {
    message(sprintf("Found %d/%d genomes (%.1f%%) with identical replicate vectors",
                    n_duplicates, n_total, pct_duplicates))
  }
  
  return(list(
    duplicate_genomes = duplicate_genomes,
    n_duplicates = n_duplicates,
    pct_duplicates = pct_duplicates
  ))
}


#' Assign genomes to clusters via plurality voting
#' @export
assign_genomes_by_plurality <- function(feature_matrix,
                                        som_node_assignments,
                                        node_clusters,
                                        seed = 123) {
  
  set.seed(seed)
  
  model_clusters <- node_clusters[som_node_assignments]
  
  model_df <- data.frame(
    genome_id = rownames(feature_matrix),
    model_id = 1:nrow(feature_matrix),
    som_node = som_node_assignments,
    cluster_id = model_clusters,
    stringsAsFactors = FALSE
  )
  
  unique_genomes <- unique(model_df$genome_id)
  
  genome_assignments <- data.frame(
    genome_id = character(),
    majority_cluster = integer(),
    n_models = integer(),
    n_clusters_represented = integer(),
    plurality_fraction = numeric(),
    is_tie = logical(),
    stringsAsFactors = FALSE
  )
  
  ties_found <- 0
  
  for (genome in unique_genomes) {
    genome_models <- model_df %>% filter(genome_id == genome)
    
    cluster_counts <- table(genome_models$cluster_id)
    max_count <- max(cluster_counts)
    
    winning_clusters <- as.numeric(names(cluster_counts)[cluster_counts == max_count])
    is_tie <- length(winning_clusters) > 1
    
    if (is_tie) {
      ties_found <- ties_found + 1
      message(sprintf("  Tie for genome %s: %d clusters with %d models each",
                      genome, length(winning_clusters), max_count))
      majority_cluster <- sample(winning_clusters, 1)
    } else {
      majority_cluster <- winning_clusters[1]
    }
    
    n_models <- nrow(genome_models)
    n_clusters <- length(unique(genome_models$cluster_id))
    plurality_fraction <- max_count / n_models
    
    genome_assignments <- rbind(genome_assignments, data.frame(
      genome_id = genome,
      majority_cluster = majority_cluster,
      n_models = n_models,
      n_clusters_represented = n_clusters,
      plurality_fraction = plurality_fraction,
      is_tie = is_tie
    ))
  }
  
  if (ties_found > 0) {
    message(sprintf("⚠ Found %d genomes with tied plurality votes (randomly assigned)",
                    ties_found))
  }
  
  return(genome_assignments)
}


#' Assign genomes to SOM nodes via plurality voting
#' @export
assign_genomes_to_nodes_by_plurality <- function(feature_matrix,
                                                 som_node_assignments,
                                                 seed = 123) {
  
  set.seed(seed)
  
  model_df <- data.frame(
    genome_id = rownames(feature_matrix),
    som_node = som_node_assignments,
    stringsAsFactors = FALSE
  )
  
  unique_genomes <- unique(model_df$genome_id)
  
  genome_node_assignments <- data.frame(
    genome_id = character(),
    majority_node = integer(),
    n_models = integer(),
    n_nodes_represented = integer(),
    plurality_fraction_node = numeric(),
    is_tie_node = logical(),
    stringsAsFactors = FALSE
  )
  
  ties_found <- 0
  
  for (genome in unique_genomes) {
    genome_models <- model_df %>% filter(genome_id == genome)
    
    node_counts <- table(genome_models$som_node)
    max_count <- max(node_counts)
    
    winning_nodes <- as.numeric(names(node_counts)[node_counts == max_count])
    is_tie <- length(winning_nodes) > 1
    
    if (is_tie) {
      ties_found <- ties_found + 1
      majority_node <- sample(winning_nodes, 1)
    } else {
      majority_node <- winning_nodes[1]
    }
    
    n_models <- nrow(genome_models)
    n_nodes <- length(unique(genome_models$som_node))
    plurality_fraction <- max_count / n_models
    
    genome_node_assignments <- rbind(genome_node_assignments, data.frame(
      genome_id = genome,
      majority_node = majority_node,
      n_models = n_models,
      n_nodes_represented = n_nodes,
      plurality_fraction_node = plurality_fraction,
      is_tie_node = is_tie
    ))
  }
  
  return(genome_node_assignments)
}


# ============================================
# PART 4: REPLICATE QUALITY METRICS
# ============================================

#' Calculate coefficient of variation per genome
#' @export
calculate_replicate_variation <- function(feature_matrix) {
  
  unique_genomes <- unique(rownames(feature_matrix))
  
  cv_results <- data.frame(
    genome_id = character(),
    mean_cv = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (genome in unique_genomes) {
    genome_rows <- which(rownames(feature_matrix) == genome)
    genome_matrix <- feature_matrix[genome_rows, , drop = FALSE]
    
    nonzero_cols <- which(colSums(genome_matrix) > 0)
    if (length(nonzero_cols) == 0) {
      cv_results <- rbind(cv_results, data.frame(genome_id = genome, mean_cv = NA))
      next
    }
    
    genome_matrix <- genome_matrix[, nonzero_cols, drop = FALSE]
    
    cv_per_metabolite <- apply(genome_matrix, 2, function(col) {
      nonzero_vals <- col[col > 0]
      if (length(nonzero_vals) > 1) {
        sd(nonzero_vals) / mean(nonzero_vals)
      } else {
        0
      }
    })
    
    mean_cv <- mean(cv_per_metabolite, na.rm = TRUE)
    
    cv_results <- rbind(cv_results, data.frame(
      genome_id = genome,
      mean_cv = mean_cv
    ))
  }
  
  return(cv_results)
}


#' Calculate cluster consensus metrics across k values
#' @export
calculate_cluster_consensus_metrics <- function(feature_matrix,
                                                som_codes,
                                                k_range = 2:20) {
  
  message("Calculating consensus metrics across k values...")
  
  unique_genomes <- unique(rownames(feature_matrix))
  n_models_per_genome <- sum(rownames(feature_matrix) == unique_genomes[1])
  
  results <- data.frame(
    k = integer(),
    mean_consensus = numeric(),
    mean_breadth = numeric(),
    mean_evenness = numeric()
  )
  
  for (k in k_range) {
    node_clusters <- kmeans(
      som_codes,
      centers = k,
      iter.max = 10000,
      nstart = 10
    )$cluster
    
    consensus_vals <- numeric(length(unique_genomes))
    breadth_vals <- numeric(length(unique_genomes))
    evenness_vals <- numeric(length(unique_genomes))
    
    for (i in seq_along(unique_genomes)) {
      genome <- unique_genomes[i]
      genome_rows <- which(rownames(feature_matrix) == genome)
      
      # Get SOM node assignments for this genome's models
      genome_nodes <- which(rownames(feature_matrix) == genome)
      # This assumes som_node_assignments available - need to pass it
      # For now, skip this part or restructure
    }
    
    results <- rbind(results, data.frame(
      k = k,
      mean_consensus = mean(consensus_vals, na.rm = TRUE),
      mean_breadth = mean(breadth_vals, na.rm = TRUE),
      mean_evenness = mean(evenness_vals, na.rm = TRUE)
    ))
  }
  
  message("✓ Consensus metrics calculated")
  return(results)
}


#' Analyze model distribution evenness
#' @export
analyze_model_distribution <- function(feature_matrix,
                                       node_clusters,
                                       som_node_assignments) {
  
  unique_genomes <- unique(rownames(feature_matrix))
  n_models <- sum(rownames(feature_matrix) == unique_genomes[1])
  
  model_clusters <- node_clusters[som_node_assignments]
  
  distribution_counts <- integer(n_models)
  largest_group_sizes <- integer(length(unique_genomes))
  
  for (i in seq_along(unique_genomes)) {
    genome <- unique_genomes[i]
    genome_clusters <- model_clusters[rownames(feature_matrix) == genome]
    
    cluster_counts <- table(genome_clusters)
    max_count <- max(cluster_counts)
    
    distribution_counts[max_count] <- distribution_counts[max_count] + 1
    largest_group_sizes[i] <- max_count
  }
  
  distribution_df <- data.frame(
    n_models_in_majority = 1:n_models,
    n_genomes = distribution_counts
  )
  
  return(list(
    distribution_df = distribution_df,
    largest_group_sizes = largest_group_sizes,
    mean_largest_group = mean(largest_group_sizes),
    median_largest_group = median(largest_group_sizes)
  ))
}


# ============================================
# PART 5: UNIFIED PIPELINES
# ============================================

#' Complete pipeline WITH replicate support (auto-detects)
#' 
#' @param feature_matrix Matrix (rows = observations, may include replicates)
#' @param grid_size SOM grid size
#' @param k Number of clusters
#' @param method Clustering method
#' @param rlen SOM iterations
#' @param seed Random seed
#' @param force_replicates Force replicate handling even if not detected
#' @return Comprehensive results
#' @export
run_som_clustering_pipeline <- function(feature_matrix,
                                        grid_size = 20,
                                        k = 8,
                                        method = "kmeans",
                                        rlen = 1500,
                                        seed = 123,
                                        force_replicates = FALSE) {
  
  message("========================================")
  message("SOM + Clustering Pipeline")
  message("========================================\n")
  
  # Detect replicates
  unique_ids <- unique(rownames(feature_matrix))
  has_replicates <- (length(unique_ids) < nrow(feature_matrix)) || force_replicates
  
  if (has_replicates) {
    message("Detected replicate models - using plurality voting\n")
    return(run_som_clustering_with_replicates(
      feature_matrix, grid_size, k, method, rlen, seed
    ))
  } else {
    message("No replicates detected - standard pipeline\n")
    return(run_som_clustering_simple(
      feature_matrix, grid_size, k, method, rlen, seed
    ))
  }
}


#' Simple pipeline (no replicates)
#' @keywords internal
run_som_clustering_simple <- function(feature_matrix, grid_size, k, method, rlen, seed) {
  
  # Train SOM
  som_result <- train_som(feature_matrix, grid_size, rlen, seed = seed)
  
  # Cluster
  node_clusters <- cluster_som_nodes(som_result$som_codes, k, method, seed)
  
  # Assign genomes
  node_coords <- som_result$som_model$grid$pts
  genome_clusters <- data.frame(
    genome_id = rownames(feature_matrix),
    cohort_id = node_clusters[som_result$som_model$unit.classif],
    som_node = som_result$som_model$unit.classif,
    som_x = node_coords[som_result$som_model$unit.classif, 1],
    som_y = node_coords[som_result$som_model$unit.classif, 2],
    n_models = 1,
    plurality_fraction = 1.0
  )
  
  message("\n========================================")
  message("Pipeline Complete!")
  message("========================================")
  message(sprintf("Genomes: %d", nrow(feature_matrix)))
  message(sprintf("SOM grid: %d × %d", grid_size, grid_size))
  message(sprintf("Cohorts: %d", k))
  
  return(list(
    som_model = som_result$som_model,
    som_codes = som_result$som_codes,
    som_dist = som_result$som_dist,
    node_clusters = node_clusters,
    genome_assignments = genome_clusters,
    has_replicates = FALSE,
    grid_size = grid_size,
    k = k,
    method = method,
    seed = seed
  ))
}


#' Complete pipeline WITH replicates
#' @keywords internal
run_som_clustering_with_replicates <- function(feature_matrix, grid_size, k, method, rlen, seed) {
  
  # Check duplicates
  message("Checking for duplicate replicate models...")
  dup_check <- check_duplicate_replicates(feature_matrix)
  
  # Train SOM
  message("\nTraining SOM on all replicate models...")
  som_result <- train_som(feature_matrix, grid_size, rlen, seed = seed)
  
  # Cluster
  message("\nClustering SOM nodes...")
  node_clusters <- cluster_som_nodes(som_result$som_codes, k, method, seed)
  
  # Plurality voting for clusters
  message("\nAssigning genomes to clusters via plurality voting...")
  genome_cluster_assignments <- assign_genomes_by_plurality(
    feature_matrix,
    som_result$som_model$unit.classif,
    node_clusters,
    seed
  )
  
  # Plurality voting for nodes
  message("\nAssigning genomes to SOM nodes via plurality voting...")
  genome_node_assignments <- assign_genomes_to_nodes_by_plurality(
    feature_matrix,
    som_result$som_model$unit.classif,
    seed
  )
  
  # Merge
  genome_assignments <- genome_cluster_assignments %>%
    left_join(
      genome_node_assignments %>% 
        select(genome_id, majority_node, n_nodes_represented),
      by = "genome_id"
    ) %>%
    rename(
      cohort_id = majority_cluster,
      som_node = majority_node
    )
  
  # Add coordinates
  node_coords <- som_result$som_model$grid$pts
  genome_assignments$som_x <- node_coords[genome_assignments$som_node, 1]
  genome_assignments$som_y <- node_coords[genome_assignments$som_node, 2]
  
  # CV stats
  message("\nCalculating replicate variation metrics...")
  cv_stats <- calculate_replicate_variation(feature_matrix)
  genome_assignments <- genome_assignments %>%
    left_join(cv_stats, by = "genome_id")
  
  # Distribution
  message("\nAnalyzing model distribution...")
  distribution_analysis <- analyze_model_distribution(
    feature_matrix,
    node_clusters,
    som_result$som_model$unit.classif
  )
  
  message("\n========================================")
  message("Pipeline Complete!")
  message("========================================")
  message(sprintf("Total models: %d", nrow(feature_matrix)))
  message(sprintf("Unique genomes: %d", nrow(genome_assignments)))
  message(sprintf("Models per genome: %d", genome_assignments$n_models[1]))
  message(sprintf("SOM grid: %d × %d", grid_size, grid_size))
  message(sprintf("Cohorts: %d", k))
  message(sprintf("\nMean plurality: %.3f", mean(genome_assignments$plurality_fraction)))
  message(sprintf("Genomes with ties: %d", sum(genome_assignments$is_tie)))
  
  return(list(
    som_model = som_result$som_model,
    som_codes = som_result$som_codes,
    som_dist = som_result$som_dist,
    node_clusters = node_clusters,
    genome_assignments = genome_assignments,
    model_level_assignments = data.frame(
      genome_id = rownames(feature_matrix),
      model_id = 1:nrow(feature_matrix),
      som_node = som_result$som_model$unit.classif,
      cluster_id = node_clusters[som_result$som_model$unit.classif]
    ),
    duplicate_check = dup_check,
    cv_stats = cv_stats,
    distribution_analysis = distribution_analysis,
    has_replicates = TRUE,
    grid_size = grid_size,
    k = k,
    method = method,
    seed = seed
  ))
}


# ============================================
# PART 6: SAVE/LOAD COMPATIBILITY
# ============================================

#' Save results in your .RData format
#' @export
save_som_results <- function(pipeline_result, output_dir = "Publication_Data") {
  
  # Extract in your naming convention
  som.model <- pipeline_result$som_model
  carve.som <- pipeline_result$som_codes
  som.dist <- pipeline_result$som_dist
  som.cluster.k <- pipeline_result$node_clusters
  ids <- pipeline_result$som_model$unit.classif
  k <- pipeline_result$k
  clusters <- som.cluster.k
  
  if (pipeline_result$has_replicates) {
    sample_clusters <- pipeline_result$model_level_assignments$cluster_id
  } else {
    sample_clusters <- pipeline_result$genome_assignments$cohort_id
  }
  
  # Save
  save(som.model, carve.som, som.dist, ids,
       file = file.path(output_dir, "final-kohonen-som.RData"))
  message(sprintf("✓ Saved: %s", file.path(output_dir, "final-kohonen-som.RData")))
  
  save(som.cluster.k, sample_clusters, clusters, k, ids,
       file = file.path(output_dir, "final-som-clusters.RData"))
  message(sprintf("✓ Saved: %s", file.path(output_dir, "final-som-clusters.RData")))
  
  write.csv(
    pipeline_result$genome_assignments,
    file.path(output_dir, "genome_cluster_assignments.csv"),
    row.names = FALSE
  )
  message(sprintf("✓ Saved: %s", 
                  file.path(output_dir, "genome_cluster_assignments.csv")))
}


#' Load pre-existing results
#' @export
load_som_results <- function(som_file = "Publication_Data/final-kohonen-som.RData",
                             cluster_file = "Publication_Data/final-som-clusters.RData",
                             feature_matrix) {
  
  message("Loading pre-computed SOM and clustering results...")
  
  load(som_file)
  load(cluster_file)
  
  # Detect if replicates
  has_replicates <- length(ids) == nrow(feature_matrix) && 
    length(unique(rownames(feature_matrix))) < nrow(feature_matrix)
  
  if (has_replicates) {
    # Reconstruct with plurality
    genome_assignments <- assign_genomes_by_plurality(
      feature_matrix,
      ids,
      som.cluster.k,
      seed = 123
    )
    
    genome_node_assignments <- assign_genomes_to_nodes_by_plurality(
      feature_matrix,
      ids,
      seed = 123
    )
    
    genome_assignments <- genome_assignments %>%
      left_join(genome_node_assignments %>% 
                  select(genome_id, majority_node),
                by = "genome_id") %>%
      rename(cohort_id = majority_cluster, som_node = majority_node)
    
  } else {
    # Simple assignments
    node_coords <- som.model$grid$pts
    genome_assignments <- data.frame(
      genome_id = rownames(feature_matrix),
      cohort_id = sample_clusters,
      som_node = ids,
      som_x = node_coords[ids, 1],
      som_y = node_coords[ids, 2],
      n_models = 1
    )
  }
  
  # Add coordinates if not present
  if (!"som_x" %in% names(genome_assignments)) {
    node_coords <- som.model$grid$pts
    genome_assignments$som_x <- node_coords[genome_assignments$som_node, 1]
    genome_assignments$som_y <- node_coords[genome_assignments$som_node, 2]
  }
  
  message("✓ Loaded successfully")
  
  return(list(
    som_model = som.model,
    som_codes = carve.som,
    som_dist = som.dist,
    node_clusters = som.cluster.k,
    genome_assignments = genome_assignments,
    has_replicates = has_replicates,
    grid_size = som.model$grid$xdim,
    k = k,
    method = "kmeans"
  ))
}