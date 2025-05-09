
# Load required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readr)
library(knitr)
library(tidyr)

# Set seed for reproducibility
set.seed(123)

#---------------------------------
# Data Import
#---------------------------------

# Import required datasets
differential_proteins <- readr::read_csv("data/differential_proteins.csv", show_col_types = FALSE)
top_itaconylated_proteins <- readr::read_csv("data/top_itaconylated_proteins.csv", show_col_types = FALSE)
itaconylation_per_protein <- readr::read_csv("data/itaconylation_per_protein.csv", show_col_types = FALSE)
affected_proteins <- readr::read_csv("data/affected_proteins.csv", show_col_types = FALSE)
itaconyl_summary <- readr::read_csv("data/itaconyl_summary.csv", show_col_types = FALSE)
mod_sites_by_condition <- readr::read_csv("data/mod_sites_by_condition.csv", show_col_types = FALSE)

#---------------------------------
# Protein Ranking - Fold Change Magnitude
#---------------------------------

# Rank proteins based on absolute fold change magnitude
ranked_by_fc <- differential_proteins %>%
  dplyr::filter(!is.na(log2FC)) %>%
  dplyr::mutate(
    abs_log2FC = abs(log2FC),
    fc_rank = dplyr::row_number(dplyr::desc(abs_log2FC))
  ) %>%
  dplyr::arrange(fc_rank)

# Top 20 proteins by fold change magnitude
top_proteins_by_fc <- ranked_by_fc %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(protein, treatment, log2FC, abs_log2FC, p_value, p_adj, significant, fc_rank)

#---------------------------------
# Protein Ranking - Statistical Significance
#---------------------------------

# Rank proteins based on statistical significance (p_adj)
ranked_by_significance <- differential_proteins %>%
  dplyr::filter(!is.na(p_adj)) %>%
  dplyr::mutate(
    significance_rank = dplyr::row_number(p_adj)
  ) %>%
  dplyr::arrange(significance_rank)

# Top 20 proteins by statistical significance
top_proteins_by_significance <- ranked_by_significance %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(protein, treatment, log2FC, p_value, p_adj, significant, significance_rank)

#---------------------------------
# Protein Ranking - Itaconylation Sites
#---------------------------------

# Rank proteins based on number of itaconylation sites
ranked_by_itaconyl <- itaconylation_per_protein %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::mutate(
    itaconyl_rank = dplyr::row_number(dplyr::desc(total_itaconyl_sites))
  ) %>%
  dplyr::arrange(itaconyl_rank)

# Top 20 proteins by itaconylation sites
top_proteins_by_itaconyl <- ranked_by_itaconyl %>%
  dplyr::slice_head(n = 20)

#---------------------------------
# Create Composite Ranking
#---------------------------------

# Join the different ranking datasets
# First, rename protein column in ranked_by_itaconyl to match others
ranked_by_itaconyl_renamed <- ranked_by_itaconyl %>%
  dplyr::rename(protein = gene_symbol)

# Join FC and significance ranks
composite_ranking <- ranked_by_fc %>%
  dplyr::select(protein, treatment, log2FC, abs_log2FC, p_adj, significant, fc_rank) %>%
  dplyr::left_join(
    ranked_by_significance %>% dplyr::select(protein, treatment, significance_rank),
    by = c("protein", "treatment")
  )

# Add itaconylation rank if available (might have different protein identifier)
# This is a simplified approach - in a real scenario, you'd need a proper ID mapping
composite_ranking <- composite_ranking %>%
  dplyr::left_join(
    ranked_by_itaconyl_renamed %>% dplyr::select(protein, total_itaconyl_sites, itaconyl_rank),
    by = "protein"
  ) %>%
  # Fill NA values for ranks with a high number (indicating low rank)
  dplyr::mutate(
    significance_rank = dplyr::if_else(is.na(significance_rank), max(significance_rank, na.rm = TRUE) + 1, significance_rank),
    itaconyl_rank = dplyr::if_else(is.na(itaconyl_rank), max(itaconyl_rank, na.rm = TRUE) + 1, itaconyl_rank),
    # Create composite score - lower is better
    composite_score = (fc_rank + significance_rank + itaconyl_rank) / 3,
    # Rank based on composite score
    composite_rank = dplyr::row_number(composite_score)
  ) %>%
  dplyr::arrange(composite_rank)

# Top 20 proteins by composite ranking
top_proteins_composite <- composite_ranking %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(protein, treatment, log2FC, abs_log2FC, p_adj, significant, fc_rank, significance_rank, itaconyl_rank, composite_score, composite_rank)

#---------------------------------
# Visualizations - Volcano Plot
#---------------------------------

# Create volcano plot with different shapes for LPS treatment
volcano_plot <- ggplot2::ggplot(differential_proteins, ggplot2::aes(x = log2FC, y = -log10(p_adj))) +
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

# Save volcano plot
ggplot2::ggsave("figures/volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)

#---------------------------------
# Column Graphs - Affected Proteins and Modification Sites
#---------------------------------

# Count affected proteins by treatment
affected_proteins_count <- affected_proteins %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarize(count = dplyr::n(), .groups = "drop")

# Column graph for affected proteins
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

# Prepare data for modification sites column graph
mod_sites_plot_data <- mod_sites_by_condition %>%
  dplyr::filter(!is.na(condition) & !is.na(treatment)) %>%
  dplyr::mutate(group = paste(condition, treatment, sep = "_"))

# Column graph for modification sites
mod_sites_plot <- ggplot2::ggplot(mod_sites_plot_data,
                        ggplot2::aes(x = group, y = total_itaconyl_sites, fill = condition)) +
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

# Combine column graphs
column_graphs <- gridExtra::grid.arrange(
  affected_proteins_plot,
  mod_sites_plot,
  ncol = 2
)

# Save combined column graphs
ggplot2::ggsave("figures/column_graphs.png", plot = column_graphs, width = 12, height = 6, dpi = 300)

#---------------------------------
# Pie Chart - Itaconylation Sites per Peptide
#---------------------------------

# Prepare data for pie chart
pie_data <- itaconyl_summary %>%
  dplyr::mutate(
    percent = frequency / sum(frequency) * 100,
    label = paste(itaconyl_sites, "site(s):", frequency, "(", round(percent, 1), "%)"),
    pos = cumsum(percent) - 0.5 * percent
  )

# Create pie chart
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

# Save pie chart
ggplot2::ggsave("figures/sites_pie_chart.png", plot = sites_pie_chart, width = 8, height = 8, dpi = 300)

#---------------------------------
# Table - Top Proteins
#---------------------------------

# Use top_itaconylated_proteins which is already prepared

#---------------------------------
# Save Results
#---------------------------------

# Save ranking results
readr::write_csv(ranked_by_fc, "data/proteins_ranked_by_fc.csv")
readr::write_csv(ranked_by_significance, "data/proteins_ranked_by_significance.csv")
readr::write_csv(ranked_by_itaconyl, "data/proteins_ranked_by_itaconylation.csv")
readr::write_csv(composite_ranking, "data/proteins_composite_ranking.csv")

# Save top proteins lists
readr::write_csv(top_proteins_by_fc, "data/top_proteins_by_fc.csv")
readr::write_csv(top_proteins_by_significance, "data/top_proteins_by_significance.csv")
readr::write_csv(top_proteins_by_itaconyl, "data/top_proteins_by_itaconylation.csv")
readr::write_csv(top_proteins_composite, "data/top_proteins_composite.csv")

# Create and save object descriptions
ranking_objects <- data.frame(
  object_name = c(
    "ranked_by_fc", "ranked_by_significance", "ranked_by_itaconyl", "composite_ranking",
    "top_proteins_by_fc", "top_proteins_by_significance", "top_proteins_by_itaconyl", "top_proteins_composite"
  ),
  row_count = c(
    nrow(ranked_by_fc), nrow(ranked_by_significance), nrow(ranked_by_itaconyl), nrow(composite_ranking),
    nrow(top_proteins_by_fc), nrow(top_proteins_by_significance), nrow(top_proteins_by_itaconyl), nrow(top_proteins_composite)
  ),
  column_count = c(
    ncol(ranked_by_fc), ncol(ranked_by_significance), ncol(ranked_by_itaconyl), ncol(composite_ranking),
    ncol(top_proteins_by_fc), ncol(top_proteins_by_significance), ncol(top_proteins_by_itaconyl), ncol(top_proteins_composite)
  ),
  description = c(
    "Proteins ranked by fold change magnitude",
    "Proteins ranked by statistical significance",
    "Proteins ranked by number of itaconylation sites",
    "Proteins with composite ranking across multiple metrics",
    "Top 20 proteins by fold change magnitude",
    "Top 20 proteins by statistical significance",
    "Top 20 proteins by number of itaconylation sites",
    "Top 20 proteins by composite ranking"
  )
)

readr::write_csv(ranking_objects, "data/ranking_objects_description.csv")

# Display results for validation
cat("Top 10 Proteins by Fold Change Magnitude:\n")
print(knitr::kable(head(top_proteins_by_fc, 10), digits = 3))

cat("\nTop 10 Proteins by Statistical Significance:\n")
print(knitr::kable(head(top_proteins_by_significance, 10), digits = 3))

cat("\nTop 10 Proteins by Itaconylation Sites:\n")
print(knitr::kable(head(top_proteins_by_itaconyl, 10), digits = 3))

cat("\nTop 10 Proteins by Composite Ranking:\n")
print(knitr::kable(head(top_proteins_composite, 10), digits = 3))

cat("\nRanking Objects Description:\n")
print(knitr::kable(ranking_objects, digits = 3))
