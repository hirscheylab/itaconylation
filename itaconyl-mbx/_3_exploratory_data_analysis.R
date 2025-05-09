
# Load required libraries
library(stats)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(dplyr)
library(tidyr)
library(knitr)
library(here)

# Set seed for reproducibility
set.seed(123)

# Load the preprocessed data from the previous stage
cat("Loading preprocessed data...\n")

# Load the data from source directory
sumexp <- readRDS(here::here("data", "metabolomics_sumexp.rds"))
abundance_matrix <- assay(sumexp, "abundance")
sample_metadata <- as.data.frame(colData(sumexp))

# Log-transform the data
log_abundance_matrix <- log2(abundance_matrix + 1)

# Converting to long format for various analyses
log_abundance_long <- as.data.frame(log_abundance_matrix) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(cols = -feature, names_to = "sample_id", values_to = "log2_abundance") %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

# Check dimensions of the loaded data
cat("Checking dimensions of data...\n")
data_dimensions <- list(
  metabolites = nrow(log_abundance_matrix),
  samples = ncol(log_abundance_matrix)
)
dim_table <- data.frame(
  Dimension = c("Number of metabolites", "Number of samples"),
  Value = c(data_dimensions$metabolites, data_dimensions$samples)
)
knitr::kable(dim_table, caption = "Dimensions of preprocessed data")
write.csv(dim_table, "src/tmp_data/eda_data_dimensions.csv", row.names = FALSE)

# ----------------------------------------------
# Principal Component Analysis (PCA)
# ----------------------------------------------
cat("Performing Principal Component Analysis...\n")

# Transpose matrix for PCA (samples as rows)
# Center and scale the data to handle different scales of metabolites
pca_data <- stats::prcomp(t(log_abundance_matrix), center = TRUE, scale. = TRUE)
pca_variance <- (pca_data$sdev^2 / sum(pca_data$sdev^2)) * 100

# Create a data frame with the first 4 principal components
pca_results <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  PC3 = if(ncol(pca_data$x) >= 3) pca_data$x[, 3] else NA,
  PC4 = if(ncol(pca_data$x) >= 4) pca_data$x[, 4] else NA,
  sample_id = rownames(pca_data$x)
) %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

# Save PCA results
saveRDS(pca_results, "src/tmp_data/pca_results.rds")
write.csv(pca_results, "src/tmp_data/pca_results.csv", row.names = FALSE)

# Create PCA plot colored by genotype, shaped by treatment
p1 <- ggplot2::ggplot(pca_results, ggplot2::aes(x = PC1, y = PC2, color = genotype, shape = treatment)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "PCA of metabolomics data",
    subtitle = paste0("PC1 (", round(pca_variance[1], 1), "%) vs PC2 (", round(pca_variance[2], 1), "%)"),
    x = paste0("PC1 (", round(pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance[2], 1), "%)")
  )
ggplot2::ggsave("src/tmp_figures/pca_plot.png", p1, width = 10, height = 8)

# PCA plot for PC3 vs PC4 (if available)
if(ncol(pca_data$x) >= 4) {
  p2 <- ggplot2::ggplot(pca_results, ggplot2::aes(x = PC3, y = PC4, color = genotype, shape = treatment)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "PCA of metabolomics data",
      subtitle = paste0("PC3 (", round(pca_variance[3], 1), "%) vs PC4 (", round(pca_variance[4], 1), "%)"),
      x = paste0("PC3 (", round(pca_variance[3], 1), "%)"),
      y = paste0("PC4 (", round(pca_variance[4], 1), "%)")
    )
  ggplot2::ggsave("src/tmp_figures/pca_plot_pc3_pc4.png", p2, width = 10, height = 8)
}

# Scree plot of variance explained
pca_scree_data <- data.frame(
  PC = 1:length(pca_variance),
  PercentVariance = pca_variance,
  CumulativeVariance = cumsum(pca_variance)
)

p3 <- ggplot2::ggplot(pca_scree_data, ggplot2::aes(x = PC, y = PercentVariance)) +
  ggplot2::geom_bar(stat = "identity", fill = "#4D6D8E") +
  ggplot2::geom_line(ggplot2::aes(y = CumulativeVariance), color = "#7AA661", size = 1) +
  ggplot2::geom_point(ggplot2::aes(y = CumulativeVariance), color = "#7AA661", size = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Scree plot of PCA variance explained",
    x = "Principal Component",
    y = "Percent variance explained"
  )
ggplot2::ggsave("src/tmp_figures/pca_scree_plot.png", p3, width = 10, height = 6)

# Calculate PC loadings to identify top contributing metabolites
pc_loadings <- data.frame(
  metabolite = rownames(pca_data$rotation),
  PC1 = pca_data$rotation[, 1],
  PC2 = pca_data$rotation[, 2]
)

# Sort by absolute loading on PC1 and PC2
pc1_top <- pc_loadings[order(abs(pc_loadings$PC1), decreasing = TRUE), ][1:20, ]
pc2_top <- pc_loadings[order(abs(pc_loadings$PC2), decreasing = TRUE), ][1:20, ]

# Save top contributors
write.csv(pc1_top, "src/tmp_data/pc1_top_contributors.csv", row.names = FALSE)
write.csv(pc2_top, "src/tmp_data/pc2_top_contributors.csv", row.names = FALSE)

# ----------------------------------------------
# Hierarchical Clustering
# ----------------------------------------------
cat("Performing hierarchical clustering...\n")

# Sample clustering
sample_dist <- stats::dist(t(log_abundance_matrix), method = "euclidean")
sample_hclust <- stats::hclust(sample_dist, method = "complete")

# Create sample annotation for heatmap
sample_anno <- data.frame(
  Genotype = sample_metadata$genotype,
  Treatment = sample_metadata$treatment,
  row.names = sample_metadata$sample_id
)

# Save hierarchical clustering results
saveRDS(sample_hclust, "src/tmp_data/sample_hierarchical_clustering.rds")

# Plot dendrogram for samples
png("src/tmp_figures/sample_dendrogram.png", width = 10, height = 6, units = "in", res = 300)
plot(sample_hclust, main = "Sample Hierarchical Clustering",
     xlab = "", sub = "", cex = 0.8,
     hang = -1, col = "#4D6D8E")
rect.hclust(sample_hclust, k = 6, border = c("#4D6D8E", "#7AA661", "#807E7D", "#8e6e4d", "#7b4d8e", "#619ca6"))
dev.off()

# Metabolite clustering
metabolite_dist <- stats::dist(log_abundance_matrix, method = "euclidean")
metabolite_hclust <- stats::hclust(metabolite_dist, method = "complete")
saveRDS(metabolite_hclust, "src/tmp_data/metabolite_hierarchical_clustering.rds")

# Select top 50 most variable metabolites for visualization
metabolite_vars <- apply(log_abundance_matrix, 1, stats::var)
top50_metabolites <- names(sort(metabolite_vars, decreasing = TRUE))[1:min(50, length(metabolite_vars))]
top50_matrix <- log_abundance_matrix[top50_metabolites, ]

# Create heatmap for top 50 most variable metabolites
pheatmap::pheatmap(
  top50_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = sample_anno,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  color = colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100),
  main = "Heatmap of top 50 most variable metabolites",
  filename = "src/tmp_figures/top50_metabolites_heatmap.png",
  width = 12,
  height = 10
)

# ----------------------------------------------
# Correlation Analysis
# ----------------------------------------------
cat("Generating correlation heatmaps...\n")

# Create correlation matrix for top 30 most variable metabolites (for better visibility)
top30_metabolites <- names(sort(metabolite_vars, decreasing = TRUE))[1:min(30, length(metabolite_vars))]
top30_matrix <- log_abundance_matrix[top30_metabolites, ]
cor_matrix <- stats::cor(t(top30_matrix), method = "spearman")

# Create correlation heatmap
pheatmap::pheatmap(
  cor_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  color = colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100),
  main = "Correlation Heatmap of Top 30 Most Variable Metabolites",
  filename = "src/tmp_figures/metabolite_correlation_heatmap.png",
  width = 12,
  height = 12
)

# Save correlation matrix
saveRDS(cor_matrix, "src/tmp_data/metabolite_correlation_matrix.rds")
write.csv(cor_matrix, "src/tmp_data/metabolite_correlation_matrix.csv")

# ----------------------------------------------
# Distribution Analysis by Experimental Group
# ----------------------------------------------
cat("Examining distributions of metabolites across experimental groups...\n")

# Box plots of metabolite distributions by group
p4 <- ggplot2::ggplot(log_abundance_long, ggplot2::aes(x = treatment, y = log2_abundance, fill = genotype)) +
  ggplot2::geom_boxplot(outlier.size = 1, alpha = 0.8) +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(
    title = "Distribution of metabolite abundances by experimental group",
    x = "Treatment",
    y = "log2(abundance + 1)"
  )
ggplot2::ggsave("src/tmp_figures/metabolite_distribution_by_group.png", p4, width = 10, height = 6)

# Density plots of metabolite distributions by group
p5 <- ggplot2::ggplot(log_abundance_long, ggplot2::aes(x = log2_abundance, fill = genotype)) +
  ggplot2::geom_density(alpha = 0.5) +
  ggplot2::facet_grid(. ~ treatment) +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Density distribution of metabolite abundances",
    x = "log2(abundance + 1)",
    y = "Density"
  )
ggplot2::ggsave("src/tmp_figures/metabolite_density_by_group.png", p5, width = 10, height = 6)

# Calculate group statistics
group_stats <- log_abundance_long %>%
  dplyr::group_by(genotype, treatment) %>%
  dplyr::summarise(
    mean_abundance = mean(log2_abundance, na.rm = TRUE),
    median_abundance = median(log2_abundance, na.rm = TRUE),
    sd_abundance = sd(log2_abundance, na.rm = TRUE),
    min_abundance = min(log2_abundance, na.rm = TRUE),
    max_abundance = max(log2_abundance, na.rm = TRUE),
    q25 = quantile(log2_abundance, 0.25, na.rm = TRUE),
    q75 = quantile(log2_abundance, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

knitr::kable(group_stats, caption = "Summary statistics by experimental group")
write.csv(group_stats, "src/tmp_data/eda_group_statistics.csv", row.names = FALSE)

# Visualize group statistics
p6 <- ggplot2::ggplot(group_stats, ggplot2::aes(x = treatment, y = median_abundance, fill = genotype)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = q25, ymax = q75),
    position = ggplot2::position_dodge(width = 0.9),
    width = 0.25
  ) +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Median abundance by experimental group",
    x = "Treatment",
    y = "Median log2(abundance + 1)"
  )
ggplot2::ggsave("src/tmp_figures/group_median_abundance.png", p6, width = 8, height = 6)

# ----------------------------------------------
# Identify Potential Outliers
# ----------------------------------------------
cat("Identifying potential outlier samples...\n")

# Use robust methods to identify outliers
# Compute median absolute deviation for each sample
sample_mad <- apply(t(log_abundance_matrix), 1, function(x) {
  med <- median(x, na.rm = TRUE)
  mad(x, center = med, na.rm = TRUE)
})

# Calculate Z-scores using MAD
sample_medians <- apply(t(log_abundance_matrix), 1, median, na.rm = TRUE)
mad_z_scores <- abs(sample_medians - median(sample_medians)) / mad(sample_medians)

# Create outlier data frame
outlier_data <- data.frame(
  sample_id = colnames(log_abundance_matrix),
  median_abundance = sample_medians,
  mad = sample_mad,
  mad_z_score = mad_z_scores,
  is_outlier = mad_z_scores > 3.5  # Conservative threshold
) %>%
  dplyr::left_join(sample_metadata, by = "sample_id")

# Save outlier detection results
write.csv(outlier_data, "src/tmp_data/sample_outlier_detection.csv", row.names = FALSE)
saveRDS(outlier_data, "src/tmp_data/sample_outlier_detection.rds")

# Visualize potential outliers using MAD Z-scores
p7 <- ggplot2::ggplot(outlier_data, ggplot2::aes(x = sample_id, y = mad_z_score, fill = genotype)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_hline(yintercept = 3.5, linetype = "dashed", color = "#8e6e4d") +
  ggplot2::facet_grid(. ~ treatment, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggplot2::labs(
    title = "MAD Z-scores by sample",
    subtitle = "Dashed line indicates potential outlier threshold (Z-score > 3.5)",
    x = "Sample ID",
    y = "MAD Z-score"
  )
ggplot2::ggsave("src/tmp_figures/sample_mad_zscore.png", p7, width = 12, height = 7)

# ----------------------------------------------
# Individual Metabolite Analysis
# ----------------------------------------------
cat("Analyzing distribution of top metabolites...\n")

# Extract top 10 most variable metabolites
top10_metabolites <- names(sort(metabolite_vars, decreasing = TRUE))[1:min(10, length(metabolite_vars))]

# Filter data for top 10 metabolites
top10_data <- log_abundance_long %>%
  dplyr::filter(feature %in% top10_metabolites)

# Create boxplots for each of the top 10 metabolites across groups
p8 <- ggplot2::ggplot(top10_data, ggplot2::aes(x = treatment, y = log2_abundance, fill = genotype)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(~feature, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(
    title = "Distribution of top 10 most variable metabolites",
    x = "Treatment",
    y = "log2(abundance + 1)"
  )
ggplot2::ggsave("src/tmp_figures/top10_metabolites_boxplots.png", p8, width = 12, height = 10)

# Save list of top metabolites
write.csv(data.frame(Metabolite = top10_metabolites,
                     Variance = metabolite_vars[top10_metabolites]),
          "src/tmp_data/top10_variable_metabolites.csv", row.names = FALSE)

# Create a summary object to pass to the next stage
eda_summary <- list(
  dimensions = data_dimensions,
  pca_results = pca_results,
  pca_variance = pca_variance,
  top_metabolites = top_metabolites,
  group_statistics = group_stats,
  outlier_detection = outlier_data
)
saveRDS(eda_summary, "src/tmp_data/eda_summary.rds")

cat("Exploratory data analysis completed successfully!\n")
