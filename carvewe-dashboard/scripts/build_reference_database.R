# scripts/build_reference_database.R

library(DBI)
library(RSQLite)
library(dplyr)
library(readr)

# ============================================
# Configuration
# ============================================

# Paths to your publication data
DATA_DIR <- "../reproduce_publication/Publication_Data"
OUTPUT_DB <- "data/publication_reference.db"

# Create database connection
con <- dbConnect(RSQLite::SQLite(), OUTPUT_DB)

cat("Building reference database from publication data...\n\n")

# ============================================
# 1. Core Data: Genomes
# ============================================

cat("Loading genome list...\n")
# You may need to extract unique genome IDs from your sensitivity data
sensitivities <- read_csv(file.path(DATA_DIR, "sensitivities.csv"))

genomes <- data.frame(
  genome_id = unique(sensitivities$genome_id),
  publication_version = "1.0"
)

dbWriteTable(con, "ref_genomes", genomes, overwrite = TRUE)
cat(sprintf("  ✓ Loaded %d genomes\n", nrow(genomes)))

# ============================================
# 2. Core Data: Sensitivities
# ============================================

cat("Loading sensitivity data...\n")
# Adjust column names to match your actual data
sensitivities_clean <- sensitivities %>%
  select(
    genome_id,
    nutrient,
    growth_replete,
    growth_limited,
    sensitivity_score
    # Add replicate column if you have it
  )

dbWriteTable(con, "ref_sensitivities", sensitivities_clean, overwrite = TRUE)
cat(sprintf("  ✓ Loaded %d sensitivity measurements\n", nrow(sensitivities_clean)))

# ============================================
# 3. Core Data: SOM Coordinates
# ============================================

cat("Loading SOM coordinates...\n")
# If you saved these from your R Markdown
som_coords <- read_csv(file.path(DATA_DIR, "som_coordinates.csv"))

# Adjust columns as needed
som_coords_clean <- som_coords %>%
  select(
    genome_id,
    som_x,
    som_y,
    som_node,
    grid_size = 10  # Or whatever your grid size was
  )

dbWriteTable(con, "ref_som_coordinates", som_coords_clean, overwrite = TRUE)
cat(sprintf("  ✓ Loaded SOM coordinates for %d genomes\n", nrow(som_coords_clean)))

# ============================================
# 4. Core Data: Cohort Assignments
# ============================================

cat("Loading cohort assignments...\n")
# Your 8 metabolic cohorts
cohorts <- read_csv(file.path(DATA_DIR, "cluster_assignments.csv"))

# Define cohort metadata
cohort_metadata <- data.frame(
  cluster_id = 1:8,
  cluster_name = c(
    "Cohort 1: High glucose sensitivity",
    "Cohort 2: Nitrogen specialists",
    "Cohort 3: Generalists",
    "Cohort 4: Iron-limited",
    # ... customize with your actual cohort descriptions
    "Cohort 5: ...",
    "Cohort 6: ...",
    "Cohort 7: ...",
    "Cohort 8: ..."
  ),
  description = c(
    "Detailed description of cohort 1...",
    "Detailed description of cohort 2...",
    # ... add your descriptions
    "...", "...", "...", "...", "...", "..."
  )
)

dbWriteTable(con, "ref_clusters", cohort_metadata, overwrite = TRUE)

# Genome-to-cohort assignments
cohort_assignments <- cohorts %>%
  select(genome_id, cohort_id) %>%
  left_join(cohort_metadata %>% select(cohort_id, cohort_name), 
            by = "cohort_id")

dbWriteTable(con, "ref_cohort_assignments", cohort_assignments, overwrite = TRUE)
cat(sprintf("  ✓ Loaded cohort assignments for %d genomes\n", 
            nrow(cohort_assignments)))

# ============================================
# 5. Compute Cohort Centroids
# ============================================

cat("Computing cohort centroids...\n")

# Prepare feature matrix
feature_matrix <- sensitivities_clean %>%
  select(genome_id, nutrient, sensitivity_score) %>%
  tidyr::pivot_wider(names_from = nutrient, 
                     values_from = sensitivity_score) %>%
  column_to_rownames("genome_id") %>%
  as.matrix()

# Calculate centroid for each cohort
centroids_list <- lapply(1:8, function(cohort_id) {
  genomes_in_cohort <- cohort_assignments %>%
    filter(cohort_id == !!cohort_id) %>%
    pull(genome_id)
  
  cohort_data <- feature_matrix[genomes_in_cohort, , drop = FALSE]
  centroid <- colMeans(cohort_data, na.rm = TRUE)
  
  list(
    cohort_id = cohort_id,
    cohort_name = cohort_metadata$cluster_name[cohort_id],
    centroid = centroid
  )
})

# Store as serialized blobs
centroids_df <- data.frame(
  cohort_id = sapply(centroids_list, function(x) x$cohort_id),
  cohort_name = sapply(centroids_list, function(x) x$cohort_name),
  centroid_blob = I(lapply(centroids_list, function(x) serialize(x$centroid, NULL)))
)

dbWriteTable(con, "ref_cohort_centroids", centroids_df, overwrite = TRUE)
cat("  ✓ Computed centroids for 8 cohorts\n")

# ============================================
# 6. Optional: Taxonomy
# ============================================

if (file.exists(file.path(DATA_DIR, "taxonomy.csv"))) {
  cat("Loading taxonomy data...\n")
  taxonomy <- read_csv(file.path(DATA_DIR, "taxonomy.csv"))
  
  taxonomy_clean <- taxonomy %>%
    select(
      genome_id,
      domain,
      phylum,
      class,
      order_ = order,  # Rename 'order' to avoid SQL keyword
      family,
      genus,
      species,
      taxonomy_source
    )
  
  dbWriteTable(con, "ref_taxonomy", taxonomy_clean, overwrite = TRUE)
  cat(sprintf("  ✓ Loaded taxonomy for %d genomes\n", nrow(taxonomy_clean)))
}

# ============================================
# 7. Optional: Quality Metrics
# ============================================

if (file.exists(file.path(DATA_DIR, "quality_metrics.csv"))) {
  cat("Loading quality metrics...\n")
  quality <- read_csv(file.path(DATA_DIR, "quality_metrics.csv"))
  
  # Add your quality filtering logic
  quality_clean <- quality %>%
    mutate(
      quality_score = completeness - 5 * contamination,  # Example composite
      passes_filter = completeness >= 90 & contamination <= 5,
      filter_reason = case_when(
        completeness < 90 ~ "Low completeness",
        contamination > 5 ~ "High contamination",
        TRUE ~ NA_character_
      )
    )
  
  dbWriteTable(con, "ref_quality", quality_clean, overwrite = TRUE)
  cat(sprintf("  ✓ Loaded quality metrics for %d genomes\n", nrow(quality_clean)))
}

# ============================================
# 8. Optional: Growth Rates
# ============================================

if (file.exists(file.path(DATA_DIR, "growth_rates.csv"))) {
  cat("Loading growth rates...\n")
  growth <- read_csv(file.path(DATA_DIR, "growth_rates.csv"))
  
  dbWriteTable(con, "ref_growth_rates", growth, overwrite = TRUE)
  cat(sprintf("  ✓ Loaded growth rates for %d genomes\n", nrow(growth)))
}

# ============================================
# 9. Optional: Environmental Data
# ============================================

if (file.exists(file.path(DATA_DIR, "environmental_metadata.csv"))) {
  cat("Loading environmental data...\n")
  environment <- read_csv(file.path(DATA_DIR, "environmental_metadata.csv"))
  
  dbWriteTable(con, "ref_environmental", environment, overwrite = TRUE)
  cat(sprintf("  ✓ Loaded environmental data for %d genomes\n", nrow(environment)))
}

# ============================================
# 10. Pre-compute Summary Statistics
# ============================================

cat("Computing summary statistics...\n")

# Cohort-level sensitivity profiles
cohort_profiles <- sensitivities_clean %>%
  inner_join(cohort_assignments, by = "genome_id") %>%
  group_by(cohort_id, nutrient) %>%
  summarise(
    mean_sensitivity = mean(sensitivity_score, na.rm = TRUE),
    median_sensitivity = median(sensitivity_score, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity_score, na.rm = TRUE),
    q25_sensitivity = quantile(sensitivity_score, 0.25, na.rm = TRUE),
    q75_sensitivity = quantile(sensitivity_score, 0.75, na.rm = TRUE),
    n_genomes = n(),
    .groups = "drop"
  )

dbWriteTable(con, "ref_cohort_profiles", cohort_profiles, overwrite = TRUE)
cat("  ✓ Computed cohort sensitivity profiles\n")

# Taxonomy summaries (if available)
if (dbExistsTable(con, "ref_taxonomy")) {
  taxonomy_summary <- dbGetQuery(con, "
    SELECT c.cohort_id, t.phylum, COUNT(*) as n_genomes
    FROM ref_cohort_assignments c
    JOIN ref_taxonomy t ON c.genome_id = t.genome_id
    GROUP BY c.cohort_id, t.phylum
  ")
  
  taxonomy_summary <- taxonomy_summary %>%
    group_by(cohort_id) %>%
    mutate(proportion = n_genomes / sum(n_genomes)) %>%
    ungroup()
  
  dbWriteTable(con, "ref_cohort_taxonomy", taxonomy_summary, overwrite = TRUE)
  cat("  ✓ Computed cohort taxonomy summaries\n")
}

# ============================================
# 11. Create Indexes
# ============================================

cat("Creating database indexes...\n")

dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_sens_genome 
                ON ref_sensitivities(genome_id)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_sens_nutrient 
                ON ref_sensitivities(nutrient)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cohort 
                ON ref_cohort_assignments(cohort_id)")

if (dbExistsTable(con, "ref_taxonomy")) {
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_taxonomy_phylum 
                  ON ref_taxonomy(phylum)")
}

if (dbExistsTable(con, "ref_quality")) {
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_quality_filter 
                  ON ref_quality(passes_filter)")
}

cat("  ✓ Created indexes\n")

# ============================================
# 12. Verify and Report
# ============================================

cat("\n========================================\n")
cat("Database build complete!\n")
cat("========================================\n\n")

# Report statistics
tables <- dbListTables(con)
cat("Tables created:\n")
for (table in tables) {
  n_rows <- dbGetQuery(con, sprintf("SELECT COUNT(*) as n FROM %s", table))$n
  cat(sprintf("  - %s: %d rows\n", table, n_rows))
}

# Get database file size
db_size <- file.size(OUTPUT_DB) / 1024^2  # MB
cat(sprintf("\nDatabase file size: %.2f MB\n", db_size))

# Close connection
dbDisconnect(con)

cat("\nDatabase saved to:", OUTPUT_DB, "\n")
cat("Ready to use with the Shiny dashboard!\n")
