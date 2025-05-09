
# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)
library(gridExtra)
library(knitr)

# Set seed for reproducibility
set.seed(123)

#---------------------------------
# Import data files
#---------------------------------

# Load peptides clean data
peptides_clean <- readr::read_csv("data/peptides_clean.csv", show_col_types = FALSE)

# Load itaconyl summary data
itaconyl_summary <- readr::read_csv("data/itaconyl_summary.csv", show_col_types = FALSE)

# Load affected proteins data for visualization
affected_proteins <- readr::read_csv("data/affected_proteins.csv", show_col_types = FALSE)

# Load top itaconylated proteins for table creation
top_itaconylated_proteins <- readr::read_csv("data/top_itaconylated_proteins.csv", show_col_types = FALSE)

#---------------------------------
# Extract modification site information
#---------------------------------

# Extract amino acid residues from modification sites
modification_site_analysis <- peptides_clean %>%
  dplyr::filter(!is.na(modification_sites)) %>%
  dplyr::mutate(
    # Extract amino acid type (K for lysine, etc.) using regex
    aa_residue = stringr::str_extract_all(
      modification_sites,
      "\\[([A-Z])[0-9]+\\(.*?\\)\\]"
    ),
    # Convert list column to character
    aa_residue = sapply(aa_residue, function(x) {
      paste(stringr::str_extract_all(x, "(?<=\\[)[A-Z](?=[0-9])"), collapse = ",")
    })
  )

# Calculate frequency of different amino acid residues
aa_residue_frequency <- modification_site_analysis %>%
  dplyr::mutate(
    # Split comma-separated residues and unnest to individual rows
    aa_residue = strsplit(aa_residue, ",")
  ) %>%
  tidyr::unnest(aa_residue) %>%
  dplyr::count(aa_residue) %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::rename(residue = aa_residue, frequency = n)

# Calculate frequency of modification sites by condition and treatment
mod_sites_by_condition <- peptides_clean %>%
  dplyr::group_by(condition, treatment) %>%
  dplyr::summarize(
    total_peptides = dplyr::n(),
    total_itaconyl_sites = sum(itaconyl_sites, na.rm = TRUE),
    avg_sites_per_peptide = mean(itaconyl_sites, na.rm = TRUE),
    .groups = "drop"
  )

#---------------------------------
# Analyze distribution of modification patterns
#---------------------------------

# Extract positions from modification sites
position_analysis <- peptides_clean %>%
  dplyr::filter(!is.na(modification_sites)) %>%
  dplyr::mutate(
    # Extract position numbers using regex
    positions = stringr::str_extract_all(
      modification_sites,
      "(?<=\\[[A-Z])[0-9]+(?=\\()"
    ),
    # Convert list column to character
    positions = sapply(positions, paste, collapse = ",")
  )

# Calculate frequency of modification at different positions
position_frequency <- position_analysis %>%
  dplyr::mutate(
    positions = strsplit(positions, ",")
  ) %>%
  tidyr::unnest(positions) %>%
  dplyr::mutate(position = as.numeric(positions)) %>%
  dplyr::count(position) %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::slice_head(n = 20) # Get top 20 positions

#---------------------------------
# Create visualizations
#---------------------------------

# 1. Bar plot of itaconylation sites per peptide frequency
sites_per_peptide_plot <- ggplot2::ggplot(itaconyl_summary,
                                 ggplot2::aes(x = factor(itaconyl_sites),
                                              y = frequency,
                                              fill = factor(itaconyl_sites))) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = frequency), vjust = -0.5) +
  ggplot2::scale_fill_manual(values = c("1" = "#4D6D8E", "2" = "#7AA661", "3" = "#8e6e4d")) +
  ggplot2::labs(
    title = "Distribution of Itaconylation Sites per Peptide",
    x = "Number of Itaconylation Sites",
    y = "Frequency",
    fill = "Sites"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

# 2. Pie chart of itaconylation sites per peptide
pie_data <- itaconyl_summary %>%
  dplyr::mutate(
    percent = frequency / sum(frequency) * 100,
    label = paste(itaconyl_sites, "site(s):", frequency, "(", round(percent, 1), "%)"),
    pos = cumsum(percent) - 0.5 * percent
  )

sites_pie_chart <- ggplot2::ggplot(pie_data, ggplot2::aes(x = "", y = percent, fill = as.factor(itaconyl_sites))) +
  ggplot2::geom_bar(width = 1, stat = "identity") +
  ggplot2::coord_polar("y", start = 0) +
  ggplot2::scale_fill_manual(values = c("1" = "#4D6D8E", "2" = "#7AA661", "3" = "#8e6e4d")) +
  ggplot2::labs(
    title = "Distribution of Itaconylation Sites per Peptide",
    fill = "Number of Sites"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  ggplot2::geom_text(ggplot2::aes(y = pos, label = label), size = 3.5)

# 3. Bar plot of amino acid residue frequency
aa_residue_plot <- ggplot2::ggplot(aa_residue_frequency,
                         ggplot2::aes(x = reorder(residue, -frequency),
                                      y = frequency,
                                      fill = residue)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c("K" = "#4D6D8E", "C" = "#7AA661", "S" = "#8e6e4d",
                                         "T" = "#7b4d8e", "Y" = "#619ca6",
                                         "H" = "#a6617a", "R" = "#3c512f")) +
  ggplot2::labs(
    title = "Frequency of Amino Acid Residues in Itaconylation Sites",
    x = "Amino Acid Residue",
    y = "Frequency",
    fill = "Residue"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

# 4. Bar plot of position frequency (top 20)
position_plot <- ggplot2::ggplot(position_frequency,
                       ggplot2::aes(x = reorder(factor(position), -n),
                                    y = n,
                                    fill = factor(position))) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = rep(c("#4D6D8E", "#7AA661", "#8e6e4d", "#7b4d8e", "#619ca6"), 4)) +
  ggplot2::labs(
    title = "Frequency of Top 20 Modification Positions",
    x = "Position",
    y = "Frequency",
    fill = "Position"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# 5. Stacked bar plot of itaconylation sites by condition and treatment
sites_by_condition_plot <- ggplot2::ggplot(mod_sites_by_condition,
                                 ggplot2::aes(x = interaction(condition, treatment),
                                              y = total_itaconyl_sites,
                                              fill = condition)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = total_itaconyl_sites), vjust = -0.5, size = 3.5) +
  ggplot2::scale_fill_manual(values = c("WT" = "#4D6D8E", "KO" = "#7AA661")) +
  ggplot2::labs(
    title = "Total Itaconylation Sites by Condition and Treatment",
    x = "Condition and Treatment",
    y = "Total Itaconylation Sites",
    fill = "Condition"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# Combine related plots into multi-panel figures
# 1. Distribution plots
distribution_plots <- gridExtra::grid.arrange(
  sites_per_peptide_plot,
  sites_pie_chart,
  ncol = 2
)

# 2. Modification patterns plots
modification_plots <- gridExtra::grid.arrange(
  aa_residue_plot,
  position_plot,
  ncol = 2
)

# Create a table showing top itaconylated proteins
# This is ready from the input data

#---------------------------------
# Save outputs
#---------------------------------

# Save data outputs
readr::write_csv(modification_site_analysis, "data/modification_site_analysis.csv")
readr::write_csv(aa_residue_frequency, "data/aa_residue_frequency.csv")
readr::write_csv(position_frequency, "data/position_frequency.csv")
readr::write_csv(mod_sites_by_condition, "data/mod_sites_by_condition.csv")

# Save plots
ggplot2::ggsave("figures/sites_per_peptide_plot.png", plot = sites_per_peptide_plot, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("figures/sites_pie_chart.png", plot = sites_pie_chart, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("figures/aa_residue_plot.png", plot = aa_residue_plot, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("figures/position_plot.png", plot = position_plot, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("figures/sites_by_condition_plot.png", plot = sites_by_condition_plot, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("figures/distribution_plots.png", plot = distribution_plots, width = 12, height = 6, dpi = 300)
ggplot2::ggsave("figures/modification_plots.png", plot = modification_plots, width = 12, height = 6, dpi = 300)

# Create and save object descriptions for verification
object_descriptions <- data.frame(
  object_name = c(
    "modification_site_analysis", "aa_residue_frequency", "position_frequency",
    "mod_sites_by_condition", "pie_data"
  ),
  row_count = c(
    nrow(modification_site_analysis), nrow(aa_residue_frequency),
    nrow(position_frequency), nrow(mod_sites_by_condition), nrow(pie_data)
  ),
  column_count = c(
    ncol(modification_site_analysis), ncol(aa_residue_frequency),
    ncol(position_frequency), ncol(mod_sites_by_condition), ncol(pie_data)
  ),
  description = c(
    "Peptides with extracted amino acid residue information",
    "Frequency of different amino acid residues in modification sites",
    "Frequency of modification at different positions",
    "Itaconylation sites summarized by condition and treatment",
    "Pie chart data with percentages for visualization"
  )
)

# Save object descriptions
readr::write_csv(object_descriptions, "data/modification_analysis_objects.csv")

#---------------------------------
# Display results
#---------------------------------

# Display object descriptions
cat("Object Descriptions:\n")
print(knitr::kable(object_descriptions, digits = 3))

# Display amino acid residue frequency
cat("\nAmino Acid Residue Frequency:\n")
print(knitr::kable(aa_residue_frequency, digits = 3))

# Display itaconylation sites by condition and treatment
cat("\nItaconylation Sites by Condition and Treatment:\n")
print(knitr::kable(mod_sites_by_condition, digits = 3))

# Display top 10 positions
cat("\nTop 10 Modification Positions:\n")
print(knitr::kable(head(position_frequency, 10), digits = 3))

# Display itaconylation sites summary
cat("\nItaconylation Sites per Peptide Summary:\n")
print(knitr::kable(itaconyl_summary, digits = 3))

# Display top itaconylated proteins
cat("\nTop Itaconylated Proteins:\n")
print(knitr::kable(top_itaconylated_proteins, digits = 3))
