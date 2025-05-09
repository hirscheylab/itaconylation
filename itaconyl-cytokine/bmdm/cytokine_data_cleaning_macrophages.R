
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

# Read the data
tryCatch({
  cytokines_data <- readxl::read_excel("src/data/cytokines.xlsx")
}, error = function(e) {
  stop("Error reading file: ", e$message)
})

# Clean the data
cleaned_data <- cytokines_data %>%
  # Parse the Genotype column
  dplyr::mutate(
    Genotype_clean = stringr::str_extract(Genotype, "wt|ko"),
    Genotype_clean = ifelse(Genotype_clean == "wt", "WT", "KO"),
    Replicate = stringr::str_extract(Genotype, "[0-9]+$"),
    # Rename lps_tx to a more descriptive name and convert to factor
    LPS_treatment_hours = as.factor(lps_tx)
  ) %>%
  # Reorder columns to put metadata first
  dplyr::select(
    ID, Genotype_clean, Replicate, LPS_treatment_hours, protein_mg_ml, 
    dplyr::everything(), 
    -Genotype, -lps_tx
  ) %>%
  # Rename columns to ensure they are R-friendly
  dplyr::rename_with(~ gsub("/", "_", .x)) %>%
  dplyr::rename_with(~ gsub("-", "_", .x))

# Save the cleaned data
tryCatch({
  openxlsx::write.xlsx(cleaned_data, "src/tmp_data_clean/cytokines_cleaned.xlsx", 
                       rowNames = FALSE, overwrite = TRUE)
}, error = function(e) {
  stop("Error writing cleaned data: ", e$message)
})
