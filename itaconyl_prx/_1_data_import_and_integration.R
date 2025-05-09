
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(123)

#---------------------------------
# Import data files
#---------------------------------

# Import itaconylation summary data
itaconyl_summary <- readxl::read_excel("data/Itaconyl_Summary.xlsx")

# Import peptide data
peptides_clean <- readxl::read_excel("data/Peptides_Clean.xlsx")

# Import protein summary data
protein_summary <- readxl::read_excel("data/Protein_Summary.xlsx")

# Import protein abundance data (wide format)
proteins_clean <- readxl::read_excel("data/Proteins_Clean.xlsx")

# Import protein abundance data (long format)
proteins_long <- readxl::read_excel("data/Proteins_Long.xlsx")

#---------------------------------
# Basic data cleaning and preparation
#---------------------------------

# Clean peptides data
peptides_clean <- readxl::read_excel("data/Peptides_Clean.xlsx")
peptides_clean <- peptides_clean |>
  mutate(condition = case_when(
    str_detect(sample, "wt") ~ "WT",
    str_detect(sample, "ko") ~ "KO",
    str_detect(sample, "control") ~ "control",
    TRUE ~ condition  # Keep existing value for other cases
  )) |>
  mutate(treatment = case_when(
    str_detect(sample, "6h") ~ "6h",
    str_detect(sample, "untreated") ~ "0h",
    str_detect(sample, "24h") ~ "24h",
    str_detect(sample, "48h") ~ "48h",
    str_detect(sample, "pool") ~ "control",
    TRUE ~ condition  # Keep existing value for other cases
  ))

peptides_clean <- peptides_clean %>%
  dplyr::filter(!is.na(abundance),
                condition != "control") %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("WT", "KO")),
    treatment = factor(treatment, levels = c("0h", "6h", "24h", "48h"))
  )

# Clean protein summary data
protein_summary_clean <- protein_summary %>%
  dplyr::filter(!is.na(protein)) %>%
  dplyr::mutate(
    treatment = factor(treatment, levels = c("untreated", "LPS")),
    # Flag significant proteins (|log2FC| > 1 as a preliminary threshold)
    significant = abs(log2FC) > 1
  )

# Clean proteins long data
proteins_long_clean <- proteins_long %>%
  dplyr::filter(!is.na(abundance)) %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("WT", "KO")),
    treatment = factor(treatment, levels = c("untreated", "LPS"))
  )

# Clean proteins wide data
proteins_clean_filtered <- proteins_clean %>%
  dplyr::select(-contains("NA"))  # Remove NA columns if any

#---------------------------------
# Create integrated datasets for further analysis
#---------------------------------

# Prepare peptide data with summarized statistics for each modification site and condition
peptide_stats <- peptides_clean %>%
  dplyr::group_by(gene_symbol, modification_sites, condition, treatment) %>%
  dplyr::summarize(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    n_observations = dplyr::n(),
    .groups = "drop"
  )

# Create protein integrated data with treatment, condition and abundance
protein_integrated <- proteins_long_clean %>%
  dplyr::group_by(protein, condition, treatment) %>%
  dplyr::summarize(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    n_observations = dplyr::n(),
    .groups = "drop"
  )

# Join protein integrated data with fold change information from protein_summary
protein_integrated_with_fc <- protein_integrated %>%
  dplyr::left_join(
    protein_summary_clean %>%
      dplyr::select(protein, treatment, log2FC, significant),
    by = c("protein", "treatment")
  )

# Calculate itaconylation frequency per protein
itaconylation_per_protein <- peptides_clean %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarize(
    total_peptides = dplyr::n(),
    total_itaconyl_sites = sum(itaconyl_sites, na.rm = TRUE),
    avg_itaconyl_sites = mean(itaconyl_sites, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(total_itaconyl_sites))

#---------------------------------
# Display first 10 rows of each processed dataset
#---------------------------------

# Display itaconyl summary
cat("Itaconylation Sites Summary:\n")
print(knitr::kable(itaconyl_summary, digits = 3))

# Display first 10 rows of processed peptides
cat("\nFirst 10 rows of processed peptide data:\n")
print(knitr::kable(head(peptides_clean, 10), digits = 3))

# Display first 10 rows of peptide statistics
cat("\nFirst 10 rows of peptide statistics:\n")
print(knitr::kable(head(peptide_stats, 10), digits = 3))

# Display first 10 rows of protein integrated data with fold change
cat("\nFirst 10 rows of integrated protein data with fold change:\n")
print(knitr::kable(head(protein_integrated_with_fc, 10), digits = 3))

# Display first 10 rows of itaconylation per protein
cat("\nTop 10 proteins by total itaconylation sites:\n")
print(knitr::kable(head(itaconylation_per_protein, 10), digits = 3))

#---------------------------------
# Save processed datasets for later use
#---------------------------------

# Save processed datasets
readr::write_csv(itaconyl_summary, "data/itaconyl_summary.csv")
readr::write_csv(peptides_clean, "data/peptides_clean.csv")
readr::write_csv(peptide_stats, "data/peptide_stats.csv")
readr::write_csv(protein_summary_clean, "data/protein_summary_clean.csv")
readr::write_csv(protein_integrated, "data/protein_integrated.csv")
readr::write_csv(protein_integrated_with_fc, "data/protein_integrated_with_fc.csv")
readr::write_csv(itaconylation_per_protein, "data/itaconylation_per_protein.csv")

# Save protein data
readr::write_csv(proteins_long_clean, "data/proteins_long_clean.csv")

# Create and save an object description for verification
object_descriptions <- data.frame(
  object_name = c(
    "itaconyl_summary", "peptides_clean", "protein_summary_clean",
    "proteins_long_clean", "peptide_stats", "protein_integrated",
    "protein_integrated_with_fc", "itaconylation_per_protein"
  ),
  row_count = c(
    nrow(itaconyl_summary), nrow(peptides_clean), nrow(protein_summary_clean),
    nrow(proteins_long_clean), nrow(peptide_stats), nrow(protein_integrated),
    nrow(protein_integrated_with_fc), nrow(itaconylation_per_protein)
  ),
  column_count = c(
    ncol(itaconyl_summary), ncol(peptides_clean), ncol(protein_summary_clean),
    ncol(proteins_long_clean), ncol(peptide_stats), ncol(protein_integrated),
    ncol(protein_integrated_with_fc), ncol(itaconylation_per_protein)
  ),
  description = c(
    "Summary of itaconylation sites frequency",
    "Cleaned peptide data with modifications",
    "Cleaned protein summary data with significance flags",
    "Cleaned protein abundance data in long format",
    "Peptide statistics grouped by gene, modification, condition, and treatment",
    "Protein integrated data with condition and treatment",
    "Protein data with fold change information",
    "Itaconylation frequency per protein"
  )
)

readr::write_csv(object_descriptions, "data/object_descriptions.csv")
print(knitr::kable(object_descriptions, digits = 3))
