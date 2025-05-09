
# Load required libraries
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(knitr)
library(pheatmap)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(123)

# Load the SummarizedExperiment object
cat("Loading SummarizedExperiment object...\n")
sumexp <- readRDS(here::here("data", "metabolomics_sumexp.rds"))

# Extract abundance matrix and sample metadata
cat("Extracting abundance matrix and sample metadata...\n")
abundance_matrix <- assay(sumexp, "abundance")
sample_metadata <- as.data.frame(colData(sumexp))
feature_metadata <- as.data.frame(rowData(sumexp))

# Save extracted data
saveRDS(abundance_matrix, "src/tmp_data/abundance_matrix.rds")
saveRDS(sample_metadata, "src/tmp_data/sample_metadata.rds")
if (nrow(feature_metadata) > 0) {
  saveRDS(feature_metadata, "src/tmp_data/feature_metadata.rds")
}

# Check data dimensions
cat("Checking data dimensions...\n")
sumexp_dims <- list(
  samples = ncol(sumexp),
  features = nrow(sumexp),
  sample_vars = colnames(colData(sumexp)),
  feature_vars = if (ncol(rowData(sumexp)) > 0) colnames(rowData(sumexp)) else "None"
)
sumexp_dims_df <- data.frame(
  Dimension = c("Number of samples", "Number of features",
                "Sample metadata variables", "Feature metadata variables"),
  Value = c(
    sumexp_dims$samples,
    sumexp_dims$features,
    paste(sumexp_dims$sample_vars, collapse = ", "),
    paste(sumexp_dims$feature_vars, collapse = ", ")
  )
)
kable(sumexp_dims_df, caption = "SummarizedExperiment Object Dimensions")
write.csv(sumexp_dims_df, "src/tmp_data/sumexp_dimensions.csv", row.names = FALSE)

# Summarize sample groups
cat("Summarizing sample groups...\n")
group_summary <- sample_metadata %>%
  group_by(genotype, treatment) %>%
  summarise(count = n(), .groups = "drop")
kable(group_summary, caption = "Sample counts by experimental group")
write.csv(group_summary, "src/tmp_data/group_summary.csv", row.names = FALSE)

# Assess missing values
cat("Assessing missing values...\n")
missing_per_sample <- colSums(is.na(abundance_matrix))
missing_per_feature <- rowSums(is.na(abundance_matrix))
missing_summary <- data.frame(
  Type = c("Total missing values", "Samples with any missing values",
           "Features with any missing values", "% missing across matrix"),
  Count = c(
    sum(is.na(abundance_matrix)),
    sum(missing_per_sample > 0),
    sum(missing_per_feature > 0),
    round(100 * sum(is.na(abundance_matrix)) / (nrow(abundance_matrix) * ncol(abundance_matrix)), 3)
  )
)
kable(missing_summary, caption = "Summary of missing values")
write.csv(missing_summary, "src/tmp_data/missing_values_summary.csv", row.names = FALSE)

# Visualize missing values by sample
missing_by_sample_df <- data.frame(
  sample_id = colnames(abundance_matrix),
  missing_count = missing_per_sample,
  missing_percent = 100 * missing_per_sample / nrow(abundance_matrix)
) %>%
  left_join(sample_metadata, by = "sample_id")

p1 <- ggplot(missing_by_sample_df, aes(x = sample_id, y = missing_percent, fill = genotype)) +
  geom_bar(stat = "identity") +
  facet_grid(genotype ~ treatment, scales = "free_x") +
  scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Percentage of missing values by sample",
       x = "Sample ID", y = "Missing values (%)")
ggsave("src/tmp_figures/missing_values_by_sample.png", p1, width = 10, height = 7)

# Check for potential outliers using boxplots of log-transformed data
cat("Checking for potential outliers...\n")
# Log-transform abundance values for better visualization
log_abundance <- log2(abundance_matrix + 1)
log_abundance_long <- as.data.frame(log_abundance) %>%
  tibble::rownames_to_column("feature") %>%
  pivot_longer(cols = -feature, names_to = "sample_id", values_to = "log2_abundance") %>%
  left_join(sample_metadata, by = "sample_id")

# Boxplot by sample
p2 <- ggplot(log_abundance_long, aes(x = sample_id, y = log2_abundance, fill = genotype)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  facet_grid(. ~ treatment, scales = "free_x") +
  scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Distribution of log2-transformed abundance values by sample",
       x = "Sample ID", y = "log2(abundance + 1)")
ggsave("src/tmp_figures/log2_abundance_boxplot.png", p2, width = 12, height = 7)

# Density plot to visualize distributions
p3 <- ggplot(log_abundance_long, aes(x = log2_abundance, fill = genotype)) +
  geom_density(alpha = 0.5) +
  facet_grid(genotype ~ treatment) +
  scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  theme_bw() +
  labs(title = "Density of log2-transformed abundance values",
       x = "log2(abundance + 1)", y = "Density")
ggsave("src/tmp_figures/log2_abundance_density.png", p3, width = 10, height = 7)

# Sample correlation heatmap to identify outliers
sample_cor <- cor(abundance_matrix, method = "spearman")
anno_col <- data.frame(
  Genotype = sample_metadata$genotype,
  Treatment = sample_metadata$treatment,
  row.names = sample_metadata$sample_id
)
color_palette <- colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100)
pheatmap::pheatmap(
  sample_cor,
  annotation_col = anno_col,
  main = "Sample Correlation Heatmap",
  color = color_palette,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  filename = "src/tmp_figures/sample_correlation_heatmap.png",
  width = 10,
  height = 10
)

# Summary statistics for each experimental group
cat("Creating summary statistics for each experimental group...\n")
group_stats <- log_abundance_long %>%
  group_by(genotype, treatment) %>%
  summarise(
    mean_log2 = mean(log2_abundance, na.rm = TRUE),
    median_log2 = median(log2_abundance, na.rm = TRUE),
    min_log2 = min(log2_abundance, na.rm = TRUE),
    max_log2 = max(log2_abundance, na.rm = TRUE),
    sd_log2 = sd(log2_abundance, na.rm = TRUE),
    q25_log2 = quantile(log2_abundance, 0.25, na.rm = TRUE),
    q75_log2 = quantile(log2_abundance, 0.75, na.rm = TRUE),
    missing_percent = 100 * sum(is.na(log2_abundance)) / n(),
    .groups = "drop"
  )
kable(group_stats, caption = "Summary statistics by experimental group", digits = 3)
write.csv(group_stats, "src/tmp_data/group_statistics.csv", row.names = FALSE)

# Visualize group statistics
p4 <- ggplot(group_stats, aes(x = treatment, y = median_log2, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = q25_log2, ymax = q75_log2),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
  theme_bw() +
  labs(title = "Median log2-transformed abundance by group",
       x = "Treatment", y = "Median log2(abundance + 1)")
ggsave("src/tmp_figures/group_median_abundance.png", p4, width = 8, height = 6)

# Save processed data for next stages
saveRDS(log_abundance, "src/tmp_data/log_abundance_matrix.rds")
saveRDS(log_abundance_long, "src/tmp_data/log_abundance_long.rds")

# Create a summary object of findings
data_exploration_summary <- list(
  dimensions = sumexp_dims,
  group_counts = group_summary,
  missing_values = missing_summary,
  group_statistics = group_stats
)
saveRDS(data_exploration_summary, "src/tmp_data/data_exploration_summary.rds")

cat("Data import and exploration completed successfully!\n")
