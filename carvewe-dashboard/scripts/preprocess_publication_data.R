# scripts/preprocess_publication_data.R
# Processes raw CarveWe flux data into clean sensitivity measurements
library(dplyr)
library(tidyr)
library(readr)

cat("========================================\n")
cat("CarveWe Publication Data Preprocessing\n")
cat("========================================\n\n")

# ============================================
# CONFIGURATION
# ============================================

OUTPUT_DIR <- "carvewe-dashboard/data"

# ============================================
# Pre-processing quality data
# ============================================
OUTPUT_QUALITY <- file.path(OUTPUT_DIR, "processed_quality_data.csv")

raw_quality <- read_csv("Publication_Data/rxn_info.csv")
quality_data <- raw_quality %>%
  rename(genome_id = ...1) %>%
  select(genome_id, mean_freq)

# ============================================
# Save quality data
# ============================================
write_csv(quality_data, OUTPUT_QUALITY)

# ============================================
# Pre-processing primary feature data
# ============================================
OUTPUT_FILE_LONG <- file.path(OUTPUT_DIR, "processed_sensitivities_long.csv")
OUTPUT_FILE_AVERAGED <- file.path(OUTPUT_DIR, "processed_sensitivities_averaged.csv")
OUTPUT_FILE_WIDE <- file.path(OUTPUT_DIR, "processed_sensitivities_wide.csv")

# ============================================
# 1. LOAD RAW DATA
# ============================================

cat("Step 1: Loading raw data...\n")

# Metabolite classifications
class_file <- read.csv("Publication_Data/classified_metabolites.csv") %>%
  distinct(Metabolite, .keep_all = TRUE) %>%
  filter(!Higher.Level.Classification %in% c("Other", "Inorganic"))

# Model flux data (all 60 models × genomes × conditions)
load(file = "Publication_Data/new-growth50-test-data_all-cats.Rdata")

# Quality-filtered genome lists
load("Publication_Data/hq-bacterial-genomes.RData")
load("Data/low-var-genomes.RData")
load("Data/hq-cyano-genomes.RData")

cat(sprintf("  ✓ Loaded flux data: %d rows\n", nrow(all_dat)))
cat(sprintf("  ✓ Loaded %d metabolite classifications\n", nrow(class_file)))
cat(sprintf("  ✓ Loaded genome filters\n\n"))

# ============================================
# 2. CLEAN METABOLITE CLASS NAMES
# ============================================

cat("Step 2: Standardizing metabolite classifications...\n")

remove_classes <- unique(all_dat$X)[-12]
metab_classes <- unique(all_dat$metab_class)[-c(7, 8, 12)]

match_rows <- match(all_dat$metab_class, metab_classes, nomatch = NA)
all_dat <- all_dat %>%
  mutate(X = case_when(
    X %in% remove_classes ~ metab_classes[match(X, remove_classes)],
    !X %in% remove_classes ~ X
  ))

cat(sprintf("  ✓ Standardized to %d nutrient classes\n", 
            length(unique(all_dat$metab_class))))
cat("  Nutrient classes:\n")
for (nc in sort(unique(all_dat$metab_class))) {
  cat(sprintf("    - %s\n", nc))
}
cat("\n")

# ============================================
# 3. SEPARATE BASELINE VS LIMITED CONDITIONS
# ============================================

cat("Step 3: Computing flux differences from baseline...\n")

# Baseline: original media (nutrient-replete)
original_dat <- all_dat %>% 
  filter(X == "original_media") %>%
  rename(
    genome_id = genome,
    nutrient_class = metab_class,
    model_id = model,
    baseline_flux = averaged_flux
  ) %>%
  select(genome_id, model_id, nutrient_class, baseline_flux)

# Limited: each nutrient class reduced/removed
adjusted_dat <- all_dat %>% 
  filter(X != "original_media") %>%
  rename(
    genome_id = genome,
    nutrient_class = metab_class,
    model_id = model,
    limited_flux = averaged_flux,
    limitation_condition = X
  ) %>%
  select(genome_id, model_id, nutrient_class, limitation_condition, limited_flux)

# Join baseline to limited conditions
combined_dat <- adjusted_dat %>%
  left_join(
    original_dat,
    by = c("genome_id", "model_id", "nutrient_class")
  )

cat(sprintf("  ✓ Baseline measurements: %d\n", nrow(original_dat)))
cat(sprintf("  ✓ Limited condition measurements: %d\n", nrow(adjusted_dat)))
cat(sprintf("  ✓ Matched measurements: %d\n", nrow(combined_dat)))

# Check for missing baseline values
n_missing <- sum(is.na(combined_dat$baseline_flux))
if (n_missing > 0) {
  cat(sprintf("  ⚠ Warning: %d measurements missing baseline flux\n", n_missing))
}
cat("\n")

# ============================================
# 4. QUALITY FILTERING
# ============================================

cat("Step 4: Applying quality filters...\n")

n_start <- nrow(combined_dat)

# 4a. Remove models with no growth (growth flux < 10^-6)
cat("  Filtering no-growth models...\n")
nogrow_models <- combined_dat %>%
  filter(nutrient_class == "Growth" & limited_flux <= 10^-6) %>%
  distinct(limitation_condition, genome_id, model_id) %>%
  mutate(key = paste(limitation_condition, genome_id, model_id))

combined_dat <- combined_dat %>%
  mutate(key = paste(limitation_condition, genome_id, model_id)) %>%
  filter(!key %in% nogrow_models$key) %>%
  select(-key)

cat(sprintf("    Removed %d no-growth models\n", 
            nrow(nogrow_models)))

# 4b. Remove rows with missing baseline flux
combined_dat <- combined_dat %>%
  filter(!is.na(baseline_flux))

cat(sprintf("    Removed %d measurements with missing baseline\n",
            n_start - nrow(combined_dat) - nrow(nogrow_models)))

# 4c. Filter to high-quality bacterial genomes only
combined_dat <- combined_dat %>%
  filter(genome_id %in% bacterial.hq.genomes)

cat(sprintf("    Retained %d high-quality bacterial genomes\n",
            length(unique(combined_dat$genome_id))))

# 4d. Exclude cyanobacteria
n_cyano_removed <- length(unique(combined_dat$genome_id[combined_dat$genome_id %in% cyano_genomes]))
combined_dat <- combined_dat %>%
  filter(!genome_id %in% cyano_genomes)

cat(sprintf("    Excluded %d cyanobacterial genomes\n", n_cyano_removed))

# 4e. Keep only low-variance genomes
n_highvar_removed <- length(unique(combined_dat$genome_id[!combined_dat$genome_id %in% low_var_genomes]))
combined_dat <- combined_dat %>%
  filter(genome_id %in% low_var_genomes)

cat(sprintf("    Excluded %d high-variance genomes\n", n_highvar_removed))

cat(sprintf("\n  ✓ Final dataset: %d genomes, %d measurements\n\n",
            length(unique(combined_dat$genome_id)),
            nrow(combined_dat)))

# ============================================
# 5. COMPUTE GROWTH CORRECTION FACTORS
# ============================================

cat("Step 5: Computing growth-corrected sensitivities...\n")

# Growth correction: standardize by baseline vs limited growth
growth_correction <- combined_dat %>%
  filter(nutrient_class == "Growth") %>%
  mutate(growth_correction_factor = baseline_flux / limited_flux) %>%
  select(limitation_condition, genome_id, model_id, growth_correction_factor)

combined_dat <- combined_dat %>%
  left_join(
    growth_correction,
    by = c("limitation_condition", "genome_id", "model_id")
  )

cat(sprintf("  ✓ Computed growth correction factors\n\n"))

# ============================================
# 6. CALCULATE SENSITIVITY SCORES
# ============================================

cat("Step 6: Calculating sensitivity scores...\n")

# Focus on growth flux as the readout
sensitivity_data <- combined_dat %>%
  filter(nutrient_class == "Growth") %>%
  mutate(
    # Sensitivity score: how much does growth decrease?
    flux_ratio = limited_flux / baseline_flux,
    sensitivity_score = case_when(
      flux_ratio <= 1 ~ 2 * (1 - flux_ratio),  # Growth decreased
      flux_ratio > 1 ~ 0                        # Growth increased (rare)
    )
  ) %>%
  select(
    genome_id,
    model_id,
    nutrient_class = limitation_condition,  # Which nutrient was limited
    baseline_flux,
    limited_flux,
    flux_ratio,
    sensitivity_score,
    growth_correction_factor
  )

cat(sprintf("  ✓ Computed %d sensitivity measurements\n", nrow(sensitivity_data)))
cat(sprintf("  ✓ Across %d genomes\n", length(unique(sensitivity_data$genome_id))))
cat(sprintf("  ✓ For %d nutrient classes\n\n", 
            length(unique(sensitivity_data$nutrient_class))))

# ============================================
# 7. SAVE MULTIPLE FORMATS
# ============================================

cat("Step 7: Saving processed data in multiple formats...\n\n")

# Format 1: LONG FORMAT (per-model, most detailed)
# Best for: databases, flexible analysis
long_format <- sensitivity_data %>%
  select(
    genome_id,
    model_id,
    nutrient_class,
    sensitivity_score,
    baseline_growth = baseline_flux,
    limited_growth = limited_flux
  ) %>%
  arrange(genome_id, model_id, nutrient_class)

write_csv(long_format, OUTPUT_FILE_LONG)
cat(sprintf("✓ Saved long format: %s\n", OUTPUT_FILE_LONG))
cat(sprintf("  Format: genome_id | model_id | nutrient_class | sensitivity_score\n"))
cat(sprintf("  Rows: %d (all %d models per genome)\n\n", nrow(long_format), 
            length(unique(long_format$model_id))))

# Format 2: AVERAGED (across models, recommended for most analyses)
# Best for: SOM clustering, visualization, dashboard
averaged_format <- sensitivity_data %>%
  group_by(genome_id, nutrient_class) %>%
  summarise(
    sensitivity_score = mean(sensitivity_score),
    sensitivity_sd = sd(sensitivity_score),
    baseline_growth_mean = mean(baseline_flux),
    limited_growth_mean = mean(limited_flux),
    n_models = n(),
    .groups = "drop"
  ) %>%
  arrange(genome_id, nutrient_class)

write_csv(averaged_format, OUTPUT_FILE_AVERAGED)
cat(sprintf("✓ Saved averaged format: %s\n", OUTPUT_FILE_AVERAGED))
cat(sprintf("  Format: genome_id | nutrient_class | sensitivity_score | sensitivity_sd\n"))
cat(sprintf("  Rows: %d (averaged across models)\n\n", nrow(averaged_format)))

# Format 3: WIDE FORMAT (for matrix operations)
# Best for: direct input to SOM, quick visual inspection
wide_format <- sensitivity_data %>%
  select(genome_id, nutrient_class, sensitivity_score) %>%
  pivot_wider(
    names_from = nutrient_class,
    values_from = sensitivity_score
  )

write_csv(wide_format, OUTPUT_FILE_WIDE)
cat(sprintf("✓ Saved wide format: %s\n", OUTPUT_FILE_WIDE))
cat(sprintf("  Format: genome_id | [one column per nutrient class]\n"))
cat(sprintf("  Dimensions: %d genomes × %d nutrients\n\n",
            nrow(wide_format),
            ncol(wide_format) - 1))

# ============================================
# 8. SAVE METADATA
# ============================================

cat("Step 8: Saving processing metadata...\n")

metadata <- list(
  processing_date = Sys.time(),
  processing_script = "scripts/preprocess_publication_data.R",
  
  # Input data
  input_files = list(
    flux_data = "Publication_Data/new-growth50-test-data_all-cats.Rdata",
    metabolite_classes = "Publication_Data/classified_metabolites.csv",
    genome_filters = c(
      "Publication_Data/hq-bacterial-genomes.RData",
      "Data/low-var-genomes.RData",
      "Data/hq-cyano-genomes.RData"
    )
  ),
  
  # Output files
  output_files = list(
    long_format = OUTPUT_FILE_LONG,
    averaged_format = OUTPUT_FILE_AVERAGED,
    wide_format = OUTPUT_FILE_WIDE
  ),
  
  # Processing steps
  filtering_criteria = list(
    bacterial_only = TRUE,
    exclude_cyanobacteria = TRUE,
    exclude_high_variance_genomes = TRUE,
    min_growth_flux = 10^-6,
    require_baseline_flux = TRUE
  ),
  
  # Dataset dimensions
  final_dataset = list(
    n_genomes = length(unique(sensitivity_data$genome_id)),
    n_nutrient_classes = length(unique(sensitivity_data$nutrient_class)),
    n_models_per_genome = length(unique(sensitivity_data$model_id)),
    n_measurements_total = nrow(sensitivity_data),
    nutrient_classes = sort(unique(sensitivity_data$nutrient_class))
  ),
  
  # Quality filtering results
  filtering_summary = list(
    genomes_excluded_cyanobacteria = n_cyano_removed,
    genomes_excluded_high_variance = n_highvar_removed,
    models_excluded_no_growth = nrow(nogrow_models)
  ),
  
  # Sensitivity score formula
  sensitivity_calculation = "sensitivity_score = 2 * (1 - limited_growth/baseline_growth) when ratio <= 1, else 0"
)

saveRDS(metadata, file.path(OUTPUT_DIR, "processing_metadata.rds"))
cat(sprintf("  ✓ Saved metadata: %s\n\n", 
            file.path(OUTPUT_DIR, "processing_metadata.rds")))

# ============================================
# 9. SUMMARY REPORT
# ============================================

cat("========================================\n")
cat("PROCESSING COMPLETE!\n")
cat("========================================\n\n")

cat("FINAL DATASET SUMMARY\n")
cat("---------------------\n")
cat(sprintf("Genomes:         %d\n", length(unique(sensitivity_data$genome_id))))
cat(sprintf("Nutrient classes: %d\n", length(unique(sensitivity_data$nutrient_class))))
cat(sprintf("Models/genome:    %d\n", length(unique(sensitivity_data$model_id))))
cat(sprintf("Total measurements: %d\n\n", nrow(sensitivity_data)))

cat("NUTRIENT CLASSES TESTED\n")
cat("-----------------------\n")
nutrient_summary <- averaged_format %>%
  group_by(nutrient_class) %>%
  summarise(
    n_genomes = n(),
    mean_sensitivity = mean(sensitivity_score),
    sd_sensitivity = sd(sensitivity_score),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_sensitivity))

print(nutrient_summary, n = Inf)

cat("\n\nSENSITIVITY SCORE DISTRIBUTION\n")
cat("------------------------------\n")
cat(sprintf("Min:       %.4f\n", min(averaged_format$sensitivity_score)))
cat(sprintf("Q1:        %.4f\n", quantile(averaged_format$sensitivity_score, 0.25)))
cat(sprintf("Median:    %.4f\n", median(averaged_format$sensitivity_score)))
cat(sprintf("Mean:      %.4f\n", mean(averaged_format$sensitivity_score)))
cat(sprintf("Q3:        %.4f\n", quantile(averaged_format$sensitivity_score, 0.75)))
cat(sprintf("Max:       %.4f\n\n", max(averaged_format$sensitivity_score)))

cat("OUTPUT FILES\n")
cat("------------\n")
cat(sprintf("1. %s\n", OUTPUT_FILE_LONG))
cat(sprintf("   → %d rows (per-model detail)\n", nrow(long_format)))
cat(sprintf("\n2. %s\n", OUTPUT_FILE_AVERAGED))
cat(sprintf("   → %d rows (averaged across models)\n", nrow(averaged_format)))
cat(sprintf("   → RECOMMENDED for SOM clustering\n"))
cat(sprintf("\n3. %s\n", OUTPUT_FILE_WIDE))
cat(sprintf("   → %d × %d matrix format\n", nrow(wide_format), ncol(wide_format)))
cat(sprintf("   → Ready for matrix operations\n\n"))

cat("========================================\n")
cat("Ready for analysis!\n")
cat("========================================\n")

# ============================================
# 10. QUICK VALIDATION PLOT
# ============================================

cat("\nGenerating validation plot...\n")

library(ggplot2)

pdf(file.path(OUTPUT_DIR, "sensitivity_distribution.pdf"), width = 10, height = 6)

p1 <- ggplot(averaged_format, aes(x = sensitivity_score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(title = "Distribution of Sensitivity Scores",
       x = "Sensitivity Score", y = "Count") +
  theme_minimal()
print(p1)

p2 <- ggplot(nutrient_summary, 
             aes(x = reorder(nutrient_class, mean_sensitivity), 
                 y = mean_sensitivity)) +
  geom_bar(stat = "identity", fill = "coral") +
  geom_errorbar(aes(ymin = mean_sensitivity - sd_sensitivity,
                    ymax = mean_sensitivity + sd_sensitivity),
                width = 0.2) +
  coord_flip() +
  labs(title = "Mean Sensitivity by Nutrient Class",
       x = "Nutrient Class", y = "Mean Sensitivity Score") +
  theme_minimal()
print(p2)

dev.off()

cat(sprintf("  ✓ Saved validation plots: %s\n", 
            file.path(OUTPUT_DIR, "sensitivity_distribution.pdf")))