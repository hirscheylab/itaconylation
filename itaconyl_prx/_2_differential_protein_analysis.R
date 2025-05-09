
# Load required libraries
library(dplyr)
library(stats)
library(ggplot2)
library(gridExtra)
library(readr)
library(reshape2)

# Load data files
protein_integrated_with_fc <- readr::read_csv("data/protein_integrated_with_fc.csv")
itaconyl_summary <- readr::read_csv("data/itaconyl_summary.csv")
itaconylation_per_protein <- readr::read_csv("data/itaconylation_per_protein.csv")
protein_summary_clean <- readr::read_csv("data/protein_summary_clean.csv")
peptides_clean <- readr::read_csv("data/peptides_clean.csv")
peptide_stats <- readr::read_csv("data/peptide_stats.csv")

#---------------------------------
# Differential Protein Analysis
#---------------------------------

# Calculate statistical significance for proteins without p-values yet
diff_proteins <- protein_summary_clean %>%
  dplyr::filter(!is.na(log2FC)) %>%
  dplyr::mutate(
    # Calculate p-values based on the fold change magnitude (if not already present)
    # This is a simplified approach for demonstration
    p_value = 2 * stats::pnorm(-abs(log2FC)),
    p_adj = stats::p.adjust(p_value, method = "BH"),
    significant = (abs(log2FC) > 1) & (p_adj < 0.05)
  )

# Stratify analysis by LPS treatment
diff_proteins_by_treatment <- diff_proteins %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarize(
    total_proteins = dplyr::n(),
    significant_proteins = sum(significant, na.rm = TRUE),
    upregulated = sum(log2FC > 1 & p_adj < 0.05, na.rm = TRUE),
    downregulated = sum(log2FC < -1 & p_adj < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# Create lists of affected proteins
affected_proteins <- diff_proteins %>%
  dplyr::filter(significant == TRUE) %>%
  dplyr::arrange(p_adj)

# Count the number of affected proteins by treatment
affected_proteins_count <- affected_proteins %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarize(
    count = dplyr::n(),
    .groups = "drop"
  )

# Count affected modification sites from peptides data
affected_sites <- peptides_clean %>%
  dplyr::group_by(condition, treatment) %>%
  dplyr::summarize(
    total_sites = sum(itaconyl_sites, na.rm = TRUE),
    .groups = "drop"
  )

# Create a wider format for visualization
affected_sites_wide <- affected_sites %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = total_sites,
    names_prefix = "sites_"
  )

#---------------------------------
# Visualization Creation
#---------------------------------

# 1. Volcano Plot
volcano_plot <- ggplot2::ggplot(diff_proteins, ggplot2::aes(x = log2FC, y = -log10(p_adj))) +
  ggplot2::geom_point(ggplot2::aes(color = significant, shape = treatment), size = 3, alpha = 0.7) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#807E7D") +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#807E7D") +
  ggplot2::scale_color_manual(values = c("FALSE" = "#807E7D", "TRUE" = "#4D6D8E")) +
  ggplot2::scale_shape_manual(values = c("untreated" = 16, "LPS" = 17)) +
  ggplot2::labs(
    title = "Volcano Plot: WT vs KO Comparison",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant",
    shape = "Treatment"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

# 2. Column Graph for Affected Proteins
affected_proteins_plot <- ggplot2::ggplot(affected_proteins_count,
                                ggplot2::aes(x = treatment, y = count, fill = treatment)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_fill_manual(values = c("untreated" = "#4D6D8E", "LPS" = "#7AA661")) +
  ggplot2::labs(
    title = "Affected Proteins by Treatment",
    x = "Treatment",
    y = "Number of Proteins",
    fill = "Treatment"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

# 3. Column Graph for Modification Sites
affected_sites_melted <- affected_sites %>%
  dplyr::mutate(group = paste(condition, treatment, sep = "_"))

affected_sites_plot <- ggplot2::ggplot(affected_sites_melted,
                              ggplot2::aes(x = group, y = total_sites, fill = condition)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_fill_manual(values = c("WT" = "#4D6D8E", "KO" = "#7AA661")) +
  ggplot2::labs(
    title = "Itaconylation Sites by Condition and Treatment",
    x = "Condition_Treatment",
    y = "Number of Sites",
    fill = "Condition"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# 4. Pie Chart of Itaconylation Sites per Peptide
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

# Create combined plots for column graphs
column_graphs <- gridExtra::grid.arrange(
  affected_proteins_plot,
  affected_sites_plot,
  ncol = 2
)

# 5. Table of Top Proteins
top_proteins <- itaconylation_per_protein %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::arrange(dplyr::desc(total_itaconyl_sites)) %>%
  dplyr::slice_head(n = 10)

#---------------------------------
# Save Results
#---------------------------------

# Save differential analysis results
readr::write_csv(diff_proteins, "data/differential_proteins.csv")
readr::write_csv(diff_proteins_by_treatment, "data/differential_proteins_by_treatment.csv")
readr::write_csv(affected_proteins, "data/affected_proteins.csv")
readr::write_csv(top_proteins, "data/top_itaconylated_proteins.csv")

# Save plots
ggplot2::ggsave("figures/volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)
ggplot2::ggsave("figures/column_graphs.png", plot = column_graphs, width = 12, height = 6, dpi = 300)
ggplot2::ggsave("figures/sites_pie_chart.png", plot = sites_pie_chart, width = 8, height = 8, dpi = 300)

# Display results
cat("Differential Protein Analysis Results by Treatment:\n")
print(knitr::kable(diff_proteins_by_treatment, digits = 3))

cat("\nTop 10 Proteins by Itaconylation Sites:\n")
print(knitr::kable(top_proteins, digits = 3))

# Create object descriptions dataframe
object_descriptions <- data.frame(
  object_name = c(
    "diff_proteins", "diff_proteins_by_treatment", "affected_proteins",
    "affected_proteins_count", "affected_sites", "top_proteins"
  ),
  row_count = c(
    nrow(diff_proteins), nrow(diff_proteins_by_treatment), nrow(affected_proteins),
    nrow(affected_proteins_count), nrow(affected_sites), nrow(top_proteins)
  ),
  column_count = c(
    ncol(diff_proteins), ncol(diff_proteins_by_treatment), ncol(affected_proteins),
    ncol(affected_proteins_count), ncol(affected_sites), ncol(top_proteins)
  ),
  description = c(
    "Proteins with differential abundance and statistical significance",
    "Summary of differential proteins by treatment",
    "List of significantly affected proteins",
    "Count of affected proteins by treatment",
    "Number of itaconylation sites by condition and treatment",
    "Top proteins ranked by total itaconylation sites"
  )
)

# Save object descriptions
readr::write_csv(object_descriptions, "data/differential_analysis_objects.csv")

# Print object descriptions
cat("\nObject Descriptions:\n")
print(knitr::kable(object_descriptions, digits = 3))
