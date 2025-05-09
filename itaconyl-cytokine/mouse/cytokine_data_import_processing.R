
#' Import and Process Cytokine Data
#' 
#' This function imports cytokine data from V-PLEX Mouse Cytokine 19-Plex Kit,
#' performs initial validation and quality checks, and prepares it for downstream analysis.
#'
#' @param data_dir Directory path containing the cytokine data (Excel file)
#' @param filename Name of the Excel file containing cytokine data
#' @return A list containing processed data objects
#' @export
import_cytokine_data <- function(data_dir = "src/data", 
                                filename = "cytokine_data_clean.xlsx") {
  
  # Import required libraries
  require(readxl)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  # Construct the full file path
  file_path <- file.path(data_dir, filename)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }
  
  # Read the Excel file
  message("Reading cytokine data from Excel file...")
  cytokine_data <- readxl::read_excel(file_path)
  
  # Validate data structure
  required_columns <- c("sample_id", "mouse_number", "mouse_id", "genotype", "time")
  missing_columns <- setdiff(required_columns, colnames(cytokine_data))
  
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # Check if cytokine columns exist
  cytokine_columns <- setdiff(colnames(cytokine_data), required_columns)
  
  if (length(cytokine_columns) == 0) {
    stop("No cytokine measurement columns found in the data.")
  }
  
  message(paste("Found", length(cytokine_columns), "cytokine measurements."))
  
  # Basic QC checks
  
  # Check for missing values
  missing_values <- cytokine_data %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    tidyr::pivot_longer(cols = everything(), 
                        names_to = "column", 
                        values_to = "na_count") %>%
    dplyr::filter(na_count > 0)
  
  if (nrow(missing_values) > 0) {
    warning("Missing values found in data: \n", 
            paste(missing_values$column, missing_values$na_count, sep = ": ", collapse = "\n"))
  }
  
  # Check for negative values in cytokine measurements
  negative_values <- cytokine_data %>%
    dplyr::select(all_of(cytokine_columns)) %>%
    dplyr::summarise(across(everything(), ~sum(. < 0, na.rm = TRUE))) %>%
    tidyr::pivot_longer(cols = everything(), 
                        names_to = "cytokine", 
                        values_to = "negative_count") %>%
    dplyr::filter(negative_count > 0)
  
  if (nrow(negative_values) > 0) {
    warning("Negative values found in cytokine measurements: \n", 
            paste(negative_values$cytokine, negative_values$negative_count, sep = ": ", collapse = "\n"))
  }
  
  # Check data distribution and create QC summary
  cytokine_summary <- cytokine_data %>%
    dplyr::select(all_of(cytokine_columns)) %>%
    dplyr::summarise(across(everything(), 
                            list(
                              min = ~min(., na.rm = TRUE),
                              max = ~max(., na.rm = TRUE),
                              mean = ~mean(., na.rm = TRUE),
                              median = ~median(., na.rm = TRUE),
                              sd = ~sd(., na.rm = TRUE),
                              zeros = ~sum(. == 0, na.rm = TRUE)
                            )))
  
  # Convert time to factor with ordered levels
  cytokine_data <- cytokine_data %>%
    dplyr::mutate(time = factor(time, levels = c("0 h", "2 h", "7 h")),
                  genotype = factor(genotype, levels = c("wt", "ko")))
  
  # Create a long format version of the data for easier analysis
  cytokine_data_long <- cytokine_data %>%
    tidyr::pivot_longer(
      cols = all_of(cytokine_columns),
      names_to = "cytokine",
      values_to = "concentration"
    )
  
  # Save processed data
  message("Saving processed data to temporary directory...")
  saveRDS(cytokine_data, file.path("src/tmp_data_dir", "cytokine_data_wide.rds"))
  saveRDS(cytokine_data_long, file.path("src/tmp_data_dir", "cytokine_data_long.rds"))
  saveRDS(cytokine_summary, file.path("src/tmp_data_dir", "cytokine_qc_summary.rds"))
  saveRDS(cytokine_columns, file.path("src/tmp_data_dir", "cytokine_names.rds"))
  
  # Return processed data as a list
  return(list(
    data_wide = cytokine_data,
    data_long = cytokine_data_long,
    qc_summary = cytokine_summary,
    cytokine_names = cytokine_columns
  ))
}

# Execute the function
result <- import_cytokine_data()

# Create a QC plot to visualize data distribution
cytokine_boxplots <- result$data_long %>%
  ggplot2::ggplot(aes(x = cytokine, y = concentration, fill = genotype)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(~time) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
  ggplot2::labs(title = "Distribution of Cytokine Concentrations", 
                y = "Concentration", 
                x = "Cytokine")

# Save QC plot
ggplot2::ggsave("src/tmp_data_dir/cytokine_distribution_qc.pdf", 
                cytokine_boxplots, 
                width = 12, 
                height = 8)

# Print summary message
message("Cytokine data import and initial processing complete. Data ready for analysis.")
