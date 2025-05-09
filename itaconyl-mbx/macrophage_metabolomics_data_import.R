
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(S4Vectors)
library(ggplot2)
library(stats)

#' Import and Process Metabolomics Data
#' 
#' @param data_dir Directory containing metabolomics data files
#' @param filename Name of the Excel file with metabolomics data
#' @return A SummarizedExperiment object with metabolomics data
#' @details This function reads metabolomics data from an Excel file,
#'          validates its structure, performs basic QC, and returns
#'          a SummarizedExperiment object for further analysis.
import_metabolomics_data <- function(data_dir, filename) {
  # Construct full file path
  file_path <- file.path(data_dir, filename)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the data
  message("Reading metabolomics data...")
  metab_data <- readxl::read_excel(file_path)
  
  # Validate expected columns
  expected_cols <- c("Compound_Name", "sample_id", "abundance", "genotype", "treatment", "replicate")
  missing_cols <- setdiff(expected_cols, colnames(metab_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for missing values
  missing_values <- sum(is.na(metab_data))
  if (missing_values > 0) {
    warning("Dataset contains ", missing_values, " missing values")
  }
  
  # Check data types
  if (!is.numeric(metab_data$abundance)) {
    stop("Abundance values must be numeric")
  }
  
  # Validate factor levels
  genotypes <- unique(metab_data$genotype)
  if (!all(c("WT", "KO") %in% genotypes)) {
    warning("Expected genotypes 'WT' and 'KO' not found. Found: ", 
            paste(genotypes, collapse = ", "))
  }
  
  # Check for all expected treatments
  treatments <- unique(metab_data$treatment)
  expected_treatments <- c("Unstimulated", "LPS_6h", "LPS_24h")
  if (!all(expected_treatments %in% treatments)) {
    warning("Not all expected treatments found. Expected: ", 
            paste(expected_treatments, collapse = ", "), 
            ". Found: ", paste(treatments, collapse = ", "))
  }
  
  # Basic QC: Check number of compounds, samples, and distribution
  n_compounds <- length(unique(metab_data$Compound_Name))
  n_samples <- length(unique(metab_data$sample_id))
  
  message("Dataset contains ", n_compounds, " metabolites across ", n_samples, " samples")
  
  # Reshape data for SummarizedExperiment object
  # Create expression matrix with compounds as rows and samples as columns
  metab_matrix <- metab_data %>%
    dplyr::select(Compound_Name, sample_id, abundance) %>%
    tidyr::pivot_wider(
      names_from = sample_id,
      values_from = abundance,
      id_cols = Compound_Name
    )
  
  # Extract row names (compound names) and remove that column
  compound_names <- metab_matrix$Compound_Name
  metab_matrix <- metab_matrix %>%
    dplyr::select(-Compound_Name) %>%
    as.matrix()
  rownames(metab_matrix) <- compound_names
  
  # Create colData (sample metadata)
  sample_metadata <- metab_data %>%
    dplyr::select(sample_id, genotype, treatment, replicate) %>%
    dplyr::distinct() %>%
    dplyr::arrange(sample_id)
  
  # Ensure matrix columns are in the same order as sample_metadata rows
  if (!all(colnames(metab_matrix) == sample_metadata$sample_id)) {
    metab_matrix <- metab_matrix[, sample_metadata$sample_id, drop = FALSE]
  }
  
  # Create rowData (feature metadata)
  feature_metadata <- data.frame(
    Compound_Name = compound_names,
    row.names = compound_names
  )
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = metab_matrix),
    colData = S4Vectors::DataFrame(sample_metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )
  
  # Perform basic QC visualization and save to output directory
  # Density plot of metabolite abundances
  abundance_density <- metab_data %>%
    ggplot2::ggplot(aes(x = log2(abundance), color = sample_id)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Density plot of log2-transformed metabolite abundances",
      x = "log2(Abundance)",
      y = "Density"
    ) +
    ggplot2::theme(legend.position = "none")
  
  # Save QC plot
  ggplot2::ggsave(
    filename = "src/tmp_data_dir/abundance_density_plot.pdf",
    plot = abundance_density,
    width = 8,
    height = 6
  )
  
  # Basic PCA using stats package
  # Log transform and center data for PCA
  log_matrix <- log2(metab_matrix)
  pca_result <- stats::prcomp(t(log_matrix), center = TRUE, scale. = TRUE)
  
  # Extract the first two principal components
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  pca_data$sample_id <- rownames(pca_data)
  
  # Add metadata
  pca_data <- dplyr::left_join(
    pca_data,
    sample_metadata,
    by = "sample_id"
  )
  
  # Create PCA plot
  pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = genotype, shape = treatment)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "PCA of metabolite profiles",
      x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
    ) +
    ggplot2::theme_minimal()
  
  # Save PCA plot
  ggplot2::ggsave(
    filename = "src/tmp_data_dir/pca_plot.pdf",
    plot = pca_plot,
    width = 8,
    height = 6
  )
  
  message("QC completed. SummarizedExperiment object created with ",
          nrow(se), " metabolites and ", ncol(se), " samples")
  
  return(se)
}

# Execute the import function
data_dir <- "src/data"
filename <- "cleaned_macrophage_metabolite_data.xlsx"

# Import the data
metabolomics_se <- import_metabolomics_data(data_dir, filename)

# Save the SummarizedExperiment object for downstream analysis
saveRDS(metabolomics_se, file = "src/tmp_data_dir/metabolomics_sumexp.rds")

# Also save a simpler format for potential use with other tools
metabolomics_df <- as.data.frame(cbind(
  rowData(metabolomics_se),
  assay(metabolomics_se, "abundance")
))

readr::write_csv(metabolomics_df, file = "src/tmp_data_dir/metabolomics_matrix.csv")

# Save experiment metadata
experiment_metadata <- as.data.frame(colData(metabolomics_se))
readr::write_csv(experiment_metadata, file = "src/tmp_data_dir/experiment_metadata.csv")

# Print summary information about the imported data
message("Metabolomics data import complete")
message("Data saved to src/tmp_data_dir/")
message("Data contains ", nrow(metabolomics_se), " metabolites across ", ncol(metabolomics_se), " samples")
message("Ready for downstream statistical analysis")
