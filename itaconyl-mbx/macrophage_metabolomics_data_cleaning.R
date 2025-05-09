
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

# Read the Excel file
tryCatch({
  metabolite_data <- readxl::read_excel("src/data/Macrophage metabolite data - Hirschey Lab.xlsx")
  
  # Check if data was loaded properly
  if(is.null(metabolite_data) || nrow(metabolite_data) == 0) {
    stop("Failed to load data or data is empty")
  }
  
  # Convert data from wide to long format
  metabolite_data_long <- metabolite_data %>%
    tidyr::pivot_longer(
      cols = -Compound_Name,
      names_to = "sample_id",
      values_to = "abundance"
    )
  
  # Extract metadata from sample_id using regex
  metabolite_data_long <- metabolite_data_long %>%
    dplyr::mutate(
      genotype = stringr::str_extract(sample_id, "WT|KO"),
      treatment = case_when(
        stringr::str_detect(sample_id, "NT") ~ "Unstimulated",
        stringr::str_detect(sample_id, "6hr") ~ "LPS_6hr",
        stringr::str_detect(sample_id, "24hr") ~ "LPS_24hr",
        TRUE ~ NA_character_
      ),
      replicate = stringr::str_extract(sample_id, "[0-9]+$")
    )
  
  # Create a workbook to save the cleaned data
  wb <- openxlsx::createWorkbook()
  
  # Add a worksheet for the long format data
  openxlsx::addWorksheet(wb, "metabolites_long")
  openxlsx::writeData(wb, "metabolites_long", metabolite_data_long)
  
  # Add a worksheet for the original wide format data
  openxlsx::addWorksheet(wb, "metabolites_wide")
  openxlsx::writeData(wb, "metabolites_wide", metabolite_data)
  
  # Save the workbook
  openxlsx::saveWorkbook(wb, "src/tmp_data_clean/cleaned_macrophage_metabolite_data.xlsx", overwrite = TRUE)
  
}, error = function(e) {
  message("Error in data cleaning: ", e$message)
  stop(e$message)
})
