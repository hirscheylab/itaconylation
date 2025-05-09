
# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(knitr)

# Import cytokine data with properly formatted factor levels
cytokine_data <- readRDS("src/data/cytokine_data.rds")

# Display basic information about the dataset
cat("Dataset dimensions:", nrow(cytokine_data), "rows and", ncol(cytokine_data), "columns\n\n")

# Examine data structure and variable types
glimpse(cytokine_data)

# Display the first 10 rows for inspection
cytokine_head <- head(cytokine_data, 10)
knitr::kable(cytokine_head, caption = "First 10 rows of cytokine data")

# Check for missing values in each column
missing_values <- colSums(is.na(cytokine_data))
missing_values_df <- data.frame(
  Variable = names(missing_values),
  Missing_Count = missing_values,
  Missing_Percent = round(missing_values / nrow(cytokine_data) * 100, 2)
)

# Display missing values summary
knitr::kable(missing_values_df, caption = "Missing values summary")

# Count number of missing values in the entire dataset
total_missing <- sum(is.na(cytokine_data))
total_cells <- nrow(cytokine_data) * ncol(cytokine_data)
cat("Total missing values:", total_missing, "out of", total_cells, "cells (",
    round(total_missing/total_cells*100, 2), "%)\n\n")

# Summarize sample distribution across experimental conditions
sample_distribution <- cytokine_data %>%
  group_by(Genotype_clean, LPS_treatment_hours) %>%
  summarize(Sample_Count = n(), .groups = "drop")

# Display sample distribution
knitr::kable(sample_distribution, caption = "Sample distribution across experimental conditions")

# Create a summary of the cytokine data
cytokine_summary <- summary(cytokine_data)
print(cytokine_summary)

# Identify the cytokine columns (all columns after protein_mg_ml)
cytokine_cols <- names(cytokine_data)[6:ncol(cytokine_data)]
cat("Cytokines measured in the dataset (n =", length(cytokine_cols), "):\n")
cat(paste(cytokine_cols, collapse = ", "), "\n\n")

# Look at levels of categorical variables
cat("Genotype levels:", paste(levels(cytokine_data$Genotype_clean), collapse = ", "), "\n")
cat("LPS treatment levels:", paste(levels(cytokine_data$LPS_treatment_hours), collapse = ", "), "\n")
cat("Replicate levels:", paste(levels(cytokine_data$Replicate), collapse = ", "), "\n\n")

# Save the processed data for later stages
saveRDS(cytokine_data, "src/tmp_data/cytokine_data_processed.rds")

# Save the summary information
exploration_results <- list(
  data_dimensions = c(rows = nrow(cytokine_data), columns = ncol(cytokine_data)),
  missing_values_summary = missing_values_df,
  sample_distribution = sample_distribution,
  cytokine_columns = cytokine_cols
)
saveRDS(exploration_results, "src/tmp_data/exploration_results.rds")

# Create structure description for the exploration results
exploration_results_structure <- list(
  description = "Data exploration results from cytokine dataset",
  data_dimensions = "Dimensions of the dataset (rows and columns)",
  missing_values_summary = "Summary of missing values by variable",
  sample_distribution = "Distribution of samples across genotypes and LPS treatment conditions",
  cytokine_columns = "Names of all cytokine measurements in the dataset"
)
saveRDS(exploration_results_structure, "src/tmp_data/exploration_results_structure.rds")

# Create data structure description
cytokine_data_structure <- list(
  description = "Cytokine measurements from bone marrow-derived macrophages",
  dimensions = c(rows = nrow(cytokine_data), columns = ncol(cytokine_data)),
  variables = list(
    ID = "Sample identifier",
    Genotype_clean = "Mouse genotype (WT or KO)",
    Replicate = "Biological replicate number",
    LPS_treatment_hours = "Duration of LPS treatment in hours (0, 6, or 24)",
    protein_mg_ml = "Protein concentration in mg/ml",
    cytokines = "Various cytokine measurements in pg/ml"
  ),
  factor_levels = list(
    Genotype_clean = levels(cytokine_data$Genotype_clean),
    LPS_treatment_hours = levels(cytokine_data$LPS_treatment_hours),
    Replicate = levels(cytokine_data$Replicate)
  )
)
saveRDS(cytokine_data_structure, "src/tmp_data/cytokine_data_structure.rds")
