
# Load required libraries
library(limma)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(sva)
library(genefilter)
library(tidyr)
library(pheatmap)
library(knitr)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(123)

# Load input data
cat("Loading data from previous stage or original source...\n")
sumexp <- readRDS("src/data/metabolomics_sumexp.rds")
abundance_matrix <- assay(sumexp, "abundance")
sample_metadata <- as.data.frame(colData(sumexp))

# Create pre-processing statistics (before any transformations/normalizations)
cat("Creating pre-processing statistics...\n")
pre_stats <- data.frame(
  Metric = c("Number of metabolites", "Number of samples",
             "Mean abundance", "Median abundance",
             "Min abundance", "Max abundance",
             "% zeros", "% missing values"),
  Value = c(
    nrow(abundance_matrix),
    ncol(abundance_matrix),
    mean(abundance_matrix, na.rm = TRUE),
    median(abundance_matrix, na.rm = TRUE),
    min(abundance_matrix, na.rm = TRUE),
    max(abundance_matrix, na.rm = TRUE),
    100 * sum(abundance_matrix == 0, na.rm = TRUE) / (nrow(abundance_matrix) * ncol(abundance_matrix)),
    100 * sum(is.na(abundance_matrix)) / (nrow(abundance_matrix) * ncol(abundance_matrix))
  )
)
kable(pre_stats, caption = "Pre-processing statistics", digits = 3)
write.csv(pre_stats, "src/tmp_data/pre_processing_stats.csv", row.names = FALSE)

# Log transformation
cat("Performing log2 transformation of abundance values...\n")
log_abundance_matrix <- log2(abundance_matrix + 1)

# Normalization - using quantile normalization from limma
cat("Performing quantile normalization...\n")
normalized_matrix <- limma::normalizeBetweenArrays(log_abundance_matrix, method = "quantile")

# Visualize effect of normalization
# Convert to long format for plotting
log_abundance_long <- as.data.frame(log_abundance_matrix) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(cols = -feature, names_to = "sample_id", values_to = "log2_abundance") %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

normalized_long <- as.data.frame(normalized_matrix) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(cols = -feature, names_to = "sample_id", values_to = "normalized_abundance") %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

# Boxplots before normalization
p1 <- ggplot2::ggplot(log_abundance_long, ggplot2::aes(x = sample_id, y = log2_abundance, fill = genotype)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_grid(. ~ treatment, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggplot2::labs(title = "Log2 abundance values before normalization",
       x = "Sample", y = "log2(abundance)")
ggplot2::ggsave("src/tmp_figures/before_normalization_boxplot.png", p1, width = 12, height = 7)

# Boxplots after normalization
p2 <- ggplot2::ggplot(normalized_long, ggplot2::aes(x = sample_id, y = normalized_abundance, fill = genotype)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_grid(. ~ treatment, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggplot2::labs(title = "Abundance values after quantile normalization",
       x = "Sample", y = "Normalized abundance")
ggplot2::ggsave("src/tmp_figures/after_normalization_boxplot.png", p2, width = 12, height = 7)

# Density plots before and after normalization
# Before normalization
p3a <- ggplot2::ggplot(log_abundance_long, ggplot2::aes(x = log2_abundance, fill = sample_id)) +
  ggplot2::geom_density(alpha = 0.3) +
  ggplot2::facet_grid(genotype ~ treatment) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(title = "Density distributions before normalization",
       x = "log2(abundance)", y = "Density")
ggplot2::ggsave("src/tmp_figures/before_normalization_density.png", p3a, width = 10, height = 7)

# After normalization
p3b <- ggplot2::ggplot(normalized_long, ggplot2::aes(x = normalized_abundance, fill = sample_id)) +
  ggplot2::geom_density(alpha = 0.3) +
  ggplot2::facet_grid(genotype ~ treatment) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(title = "Density distributions after normalization",
       x = "Normalized abundance", y = "Density")
ggplot2::ggsave("src/tmp_figures/after_normalization_density.png", p3b, width = 10, height = 7)

# Filter out low-quality or low-abundance metabolites
cat("Filtering out low-quality or low-abundance metabolites...\n")

# Define less aggressive filtering criteria: keep metabolites that have expression
# above the lowest quartile in at least 25% of samples
k <- floor(ncol(normalized_matrix) * 0.25)  # 25% of samples
threshold <- quantile(normalized_matrix, 0.25)  # 25th percentile
pass_filter <- genefilter::kOverA(k, threshold)(normalized_matrix)

# Apply filter
filtered_matrix <- normalized_matrix[pass_filter, ]

# Save filtering results
filtering_results <- data.frame(
  Total_metabolites = nrow(normalized_matrix),
  Retained_metabolites = sum(pass_filter),
  Filtered_out_metabolites = sum(!pass_filter),
  Retention_percentage = 100 * sum(pass_filter) / nrow(normalized_matrix)
)
kable(filtering_results, caption = "Filtering results", digits = 2)
write.csv(filtering_results, "src/tmp_data/filtering_results.csv", row.names = FALSE)

# Create post-filtering statistics
post_stats <- data.frame(
  Metric = pre_stats$Metric,
  Before = pre_stats$Value,
  After = c(
    nrow(filtered_matrix),
    ncol(filtered_matrix),
    mean(filtered_matrix, na.rm = TRUE),
    median(filtered_matrix, na.rm = TRUE),
    min(filtered_matrix, na.rm = TRUE),
    max(filtered_matrix, na.rm = TRUE),
    100 * sum(filtered_matrix == 0, na.rm = TRUE) / (nrow(filtered_matrix) * ncol(filtered_matrix)),
    100 * sum(is.na(filtered_matrix)) / (nrow(filtered_matrix) * ncol(filtered_matrix))
  )
)
kable(post_stats, caption = "Statistics before and after filtering", digits = 3)
write.csv(post_stats, "src/tmp_data/post_filter_stats.csv", row.names = FALSE)

# Convert filtered matrix to long format
filtered_long <- as.data.frame(filtered_matrix) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(cols = -feature, names_to = "sample_id", values_to = "filtered_abundance") %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

# Visualize filtered data
p5 <- ggplot2::ggplot(filtered_long, ggplot2::aes(x = sample_id, y = filtered_abundance, fill = genotype)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_grid(. ~ treatment, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggplot2::labs(title = "Abundance values after filtering",
       x = "Sample", y = "Filtered abundance")
ggplot2::ggsave("src/tmp_figures/after_filtering_boxplot.png", p5, width = 12, height = 7)

# Check for batch effects using PCA
cat("Checking for batch effects...\n")
pca_data <- stats::prcomp(t(filtered_matrix), scale = TRUE)
pca_variance <- summary(pca_data)$importance[2,] * 100
pca_results <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  sample_id = colnames(filtered_matrix)
) %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

# PCA plot to check for batch effects
p6 <- ggplot2::ggplot(pca_results, ggplot2::aes(x = PC1, y = PC2, color = genotype, shape = treatment)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "PCA of normalized and filtered data",
       x = paste0("PC1 (", round(pca_variance[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_variance[2], 1), "%)"))
ggplot2::ggsave("src/tmp_figures/pca_after_preprocessing.png", p6, width = 10, height = 8)

# Check correlation between metabolites
cat("Analyzing correlations between metabolites...\n")
# Select top 50 most variable metabolites for visualization if we have more than 50
if (nrow(filtered_matrix) > 50) {
  metabolite_vars <- apply(filtered_matrix, 1, stats::var)
  top_metabolites <- names(sort(metabolite_vars, decreasing = TRUE))[1:50]
  correlation_matrix <- stats::cor(t(filtered_matrix[top_metabolites, ]), method = "spearman")
} else {
  correlation_matrix <- stats::cor(t(filtered_matrix), method = "spearman")
}

# Create a correlation heatmap
color_palette <- grDevices::colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100)
pheatmap::pheatmap(
  correlation_matrix,
  main = "Correlation Heatmap of Metabolites After Preprocessing",
  color = color_palette,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  filename = "src/tmp_figures/metabolite_correlation_heatmap.png",
  width = 12,
  height = 12
)

# Save results for next stages
cat("Saving results for next stages...\n")
saveRDS(filtered_matrix, "src/tmp_data/preprocessed_abundance_matrix.rds")
saveRDS(filtered_long, "src/tmp_data/preprocessed_abundance_long.rds")

# Create a summary of the preprocessing workflow
preprocessing_summary <- list(
  original_dimensions = dim(abundance_matrix),
  normalized_dimensions = dim(normalized_matrix),
  filtered_dimensions = dim(filtered_matrix),
  filtering_criteria = list(k = k, threshold = threshold),
  filtering_results = filtering_results,
  stats_comparison = post_stats
)
saveRDS(preprocessing_summary, "src/tmp_data/preprocessing_summary.rds")

cat("Quality control and preprocessing completed successfully!\n")
