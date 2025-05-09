
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)

# Read the data from Excel file
file_path <- "src/data/Itaconylpeptides_2024-11-04_PAG_clean.xlsx"

# Check if the file exists
if (!file.exists(file_path)) {
  stop("Input file not found: ", file_path)
}

# Get sheet names
sheet_names <- readxl::excel_sheets(file_path)

# Read the sheets
proteins_df <- readxl::read_excel(file_path, sheet = "Proteins")
peptides_df <- readxl::read_excel(file_path, sheet = "Peptide Gropus_Itaconyl")

# Clean up protein data
proteins_clean <- proteins_df %>%
  # Ensure condition and treatment are factors with proper levels
  dplyr::mutate(
    condition = factor(condition, levels = c("WT", "KO")),
    treatment = factor(treatment, levels = c("untreated", "LPS"))
  )

# Process peptide data
peptides_clean <- peptides_df %>%
  # Extract number of itaconylation sites from the modifications column
  dplyr::mutate(
    itaconyl_sites = stringr::str_extract(modifications_all_possible_sites, "^\\d+"),
    itaconyl_sites = as.numeric(itaconyl_sites),
    # Extract the modification sites
    modification_sites = stringr::str_extract_all(
      modifications_all_possible_sites, 
      "\\[(.*?)\\]"
    ),
    # Convert to character for storage
    modification_sites = sapply(modification_sites, function(x) paste(x, collapse = "; "))
  ) %>%
  # Join with sample metadata from proteins data
  dplyr::left_join(
    proteins_clean %>% 
      dplyr::select(sample, condition, treatment) %>% 
      dplyr::distinct(),
    by = "sample"
  )

# Reshape protein data for statistical analysis
proteins_long <- proteins_clean %>%
  tidyr::pivot_longer(
    cols = -c(sample, condition, treatment),
    names_to = "protein",
    values_to = "abundance"
  ) %>%
  # Remove any entries with NA abundance
  dplyr::filter(!is.na(abundance))

# Create a workbook to save all clean data
wb <- openxlsx::createWorkbook()

# Add the clean datasets to the workbook
openxlsx::addWorksheet(wb, "Proteins_Clean")
openxlsx::writeData(wb, "Proteins_Clean", proteins_clean)

openxlsx::addWorksheet(wb, "Peptides_Clean")
openxlsx::writeData(wb, "Peptides_Clean", peptides_clean)

openxlsx::addWorksheet(wb, "Proteins_Long")
openxlsx::writeData(wb, "Proteins_Long", proteins_long)

# Create a summary for protein comparison between conditions
protein_summary <- proteins_long %>%
  dplyr::group_by(protein, condition, treatment) %>%
  dplyr::summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(protein, treatment),
    names_from = condition,
    values_from = c(mean_abundance, sd_abundance, n)
  ) %>%
  dplyr::mutate(
    log2FC = log2(mean_abundance_KO / mean_abundance_WT)
  ) %>%
  dplyr::filter(!is.na(log2FC) & is.finite(log2FC))

openxlsx::addWorksheet(wb, "Protein_Summary")
openxlsx::writeData(wb, "Protein_Summary", protein_summary)

# Create a summary for itaconylation sites
itaconyl_summary <- peptides_clean %>%
  dplyr::group_by(itaconyl_sites) %>%
  dplyr::summarise(
    frequency = dplyr::n(),
    .groups = "drop"
  )

openxlsx::addWorksheet(wb, "Itaconyl_Summary")
openxlsx::writeData(wb, "Itaconyl_Summary", itaconyl_summary)

# Save the workbook
openxlsx::saveWorkbook(wb, "src/data/Itaconylpeptides_clean.xlsx", overwrite = TRUE)
