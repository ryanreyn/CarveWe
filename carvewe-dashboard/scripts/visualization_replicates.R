# R/visualization_replicates.R

#' Plot genome plurality distribution
#' 
#' @param genome_assignments Output from assign_genomes_by_plurality()
#' @return ggplot object
#' @export
plot_plurality_distribution <- function(genome_assignments) {
  library(ggplot2)
  
  p <- ggplot(genome_assignments, aes(x = plurality_fraction)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = median(genome_assignments$plurality_fraction),
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Distribution of Plurality Fractions",
      subtitle = sprintf("Median = %.3f", 
                         median(genome_assignments$plurality_fraction)),
      x = "Plurality Fraction (proportion of models in majority cluster)",
      y = "Number of Genomes"
    ) +
    theme_minimal()
  
  return(p)
}


#' Plot model distribution evenness
#' 
#' @param distribution_analysis Output from analyze_model_distribution()
#' @return ggplot object
#' @export
plot_model_distribution_evenness <- function(distribution_analysis) {
  library(ggplot2)
  
  plot_data <- distribution_analysis$distribution_df %>%
    filter(n_genomes > 0)
  
  p <- ggplot(plot_data, aes(x = n_models_in_majority, y = n_genomes)) +
    geom_point(size = 4) +
    geom_line() +
    labs(
      title = "Distribution of Model Agreement",
      x = "Number of Models in Majority Cluster",
      y = "Number of Genomes"
    ) +
    theme_minimal()
  
  return(p)
}


#' Plot coefficient of variation distribution
#' 
#' @param genome_assignments With mean_cv column
#' @return ggplot object
#' @export
plot_cv_distribution <- function(genome_assignments) {
  library(ggplot2)
  
  p <- ggplot(genome_assignments, aes(x = mean_cv)) +
    geom_histogram(bins = 30, fill = "coral", color = "white") +
    geom_vline(xintercept = median(genome_assignments$mean_cv),
               color = "blue", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Replicate Variation per Genome",
      subtitle = "Average coefficient of variation across metabolites",
      x = "Mean Coefficient of Variation",
      y = "Number of Genomes"
    ) +
    theme_minimal()
  
  return(p)
}