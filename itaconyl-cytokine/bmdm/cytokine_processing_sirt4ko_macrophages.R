
#' @title Import and validate cytokine dataset
#' @description Imports cytokine data from bone marrow-derived macrophages from WT and SIRT4KO mice,
#'              validates the data structure, performs basic QC, and returns a formatted object
#' @param data_dir Directory path containing the raw data files
#' @return List containing the processed cytokine data frame and QC results

process_cytokine_data <- function(data_dir) {
  # Load required libraries
  require(readxl)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  # Construct file path
  file_path <- file.path(data_dir, "cytokines_cleaned.xlsx")
  
  # Validate file existence
  if (!file.exists(file_path)) {
    stop("Cytokine data file not found at: ", file_path)
  }
  
  # Import data
  message("Importing cytokine data...")
  cytokine_data <- readxl::read_excel(file_path)
  
  # Check expected columns
  expected_cols <- c("ID", "Genotype_clean", "Replicate", "LPS_treatment_hours", "protein_mg_ml")
  missing_cols <- setdiff(expected_cols, colnames(cytokine_data))
  if (length(missing_cols) > 0) {
    stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check cytokine columns (all columns except metadata columns)
  cytokine_cols <- setdiff(colnames(cytokine_data), expected_cols)
  if (length(cytokine_cols) == 0) {
    stop("No cytokine measurement columns found in the dataset")
  }
  
  # Ensure proper data types
  cytokine_data <- cytokine_data %>%
    dplyr::mutate(
      ID = as.character(ID),
      Genotype_clean = as.factor(Genotype_clean),
      Replicate = as.factor(Replicate),
      LPS_treatment_hours = as.factor(LPS_treatment_hours),
      protein_mg_ml = as.numeric(protein_mg_ml)
    )
  
  # Check for missing values
  na_counts <- colSums(is.na(cytokine_data))
  cols_with_na <- names(na_counts[na_counts > 0])
  
  # Basic QC - Check for expected experimental design
  design_check <- cytokine_data %>%
    dplyr::group_by(Genotype_clean, LPS_treatment_hours) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  # Normalize cytokine values by protein concentration if needed
  # Here we're not normalizing, just preparing the data structure
  
  # Generate basic QC metrics
  qc_results <- list(
    sample_count = nrow(cytokine_data),
    genotype_counts = table(cytokine_data$Genotype_clean),
    treatment_counts = table(cytokine_data$LPS_treatment_hours),
    missing_values = na_counts[na_counts > 0],
    experimental_design = design_check
  )
  
  # Create basic summary statistics for cytokines
  cytokine_summary <- cytokine_data %>%
    dplyr::group_by(Genotype_clean, LPS_treatment_hours) %>%
    dplyr::summarise(across(all_of(cytokine_cols), 
                            list(mean = ~mean(., na.rm = TRUE), 
                                 sd = ~sd(., na.rm = TRUE))), 
                     .groups = "drop")
  
  # Save the processed data
  saveRDS(cytokine_data, "src/tmp_data_dir/cytokine_data.rds")
  saveRDS(qc_results, "src/tmp_data_dir/cytokine_qc_results.rds")
  saveRDS(cytokine_summary, "src/tmp_data_dir/cytokine_summary.rds")
  
  # Also save as CSV for easier access
  write.csv(cytokine_data, "src/tmp_data_dir/cytokine_data.csv", row.names = FALSE)
  
  # Return processed data and QC results
  return(list(
    data = cytokine_data,
    qc_results = qc_results,
    summary = cytokine_summary
  ))
}

# Execute the function
result <- process_cytokine_data("src/data")

# Create basic QC plots and save them
# Plot distribution of samples across experimental conditions
ggplot2::ggplot(result$data, ggplot2::aes(x = LPS_treatment_hours, fill = Genotype_clean)) +
  ggplot2::geom_bar(position = "dodge") +
  ggplot2::labs(title = "Sample distribution across experimental conditions",
       x = "LPS treatment (hours)",
       y = "Number of samples", 
       fill = "Genotype") +
  ggplot2::theme_minimal()
ggplot2::ggsave("src/tmp_data_dir/sample_distribution.png", width = 8, height = 6)

# Plot protein concentration distribution 
ggplot2::ggplot(result$data, ggplot2::aes(x = Genotype_clean, y = protein_mg_ml, fill = LPS_treatment_hours)) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(title = "Protein concentration by experimental condition",
       x = "Genotype",
       y = "Protein concentration (mg/ml)",
       fill = "LPS treatment (hours)") +
  ggplot2::theme_minimal()
ggplot2::ggsave("src/tmp_data_dir/protein_concentration.png", width = 8, height = 6)

# Output success message
message("Cytokine data successfully processed and saved to src/tmp_data_dir/")
