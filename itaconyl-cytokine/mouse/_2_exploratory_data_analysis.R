
# Load required libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(stats)
library(tidyr)
library(readr)

# Function to safely load data with error handling
load_data_safely <- function(file_path, error_message = NULL) {
  if (is.null(error_message)) {
    error_message <- paste("Error loading data from", file_path)
  }
  
  tryCatch({
    if (file.exists(file_path)) {
      readRDS(file_path)
    } else {
      stop(paste("File not found:", file_path))
    }
  }, error = function(e) {
    message(error_message)
    message(e$message)
    NULL
  })
}

# Load data from data directory
cytokine_data_wide <- load_data_safely("src/data/cytokine_data_wide.rds", 
                                      "Error loading wide-format cytokine data")
cytokine_data_long <- load_data_safely("src/data/cytokine_data_long.rds", 
                                      "Error loading long-format cytokine data")
cytokine_names <- load_data_safely("src/data/cytokine_names.rds", 
                                  "Error loading cytokine names")

# Check if data was loaded successfully
if (is.null(cytokine_data_wide) || is.null(cytokine_data_long) || is.null(cytokine_names)) {
  stop("Failed to load required data files")
}

# Print basic data info
cat("Data dimensions:\n")
cat("Wide data:", dim(cytokine_data_wide)[1], "rows,", dim(cytokine_data_wide)[2], "columns\n")
cat("Long data:", dim(cytokine_data_long)[1], "rows,", dim(cytokine_data_long)[2], "columns\n")
cat("Number of cytokines:", length(cytokine_names), "\n")
cat("Cytokine names:", paste(head(cytokine_names, 3), collapse=", "), "...\n")

# 1. Generate summary statistics for each cytokine by genotype and time point
# Create summary statistics
summary_stats <- cytokine_data_long %>%
  group_by(genotype, time, cytokine) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    min = min(concentration, na.rm = TRUE),
    max = max(concentration, na.rm = TRUE),
    .groups = "drop"
  )

# Save summary statistics
saveRDS(summary_stats, "src/tmp_data/eda_summary_stats.rds")

# Print summary stats for first few rows
cat("\nSummary statistics (first 10 rows):\n")
print(head(summary_stats, 10))

# 2. Create boxplots to visualize distributions
# Function to create boxplots for a subset of cytokines
create_cytokine_boxplots <- function(data, cytokines) {
  filtered_data <- data %>%
    filter(cytokine %in% cytokines)
  
  p <- ggplot(filtered_data, aes(x = time, y = concentration, fill = genotype)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ cytokine, scales = "free_y") +
    scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
    theme_bw() +
    labs(title = "Cytokine Concentrations by Genotype and Time",
         x = "Time Point", 
         y = "Concentration (pg/ml)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Get all unique cytokines from the long data format
unique_cytokines <- unique(cytokine_data_long$cytokine)

# Split cytokines into smaller groups for better visualization
group_size <- 5
cytokine_groups <- split(unique_cytokines, ceiling(seq_along(unique_cytokines)/group_size))

# Create and save boxplots for each group
for (i in seq_along(cytokine_groups)) {
  boxplot <- create_cytokine_boxplots(cytokine_data_long, cytokine_groups[[i]])
  
  # Save the boxplot
  filename <- paste0("src/tmp_figures/eda_boxplot_group_", i, ".png")
  ggsave(filename, boxplot, width = 10, height = 8)
  
  # Create log-scale version for better visualization
  log_boxplot <- boxplot + 
    scale_y_log10() +
    labs(title = "Cytokine Concentrations by Genotype and Time (Log Scale)",
         y = "Concentration (pg/ml, log scale)")
  
  # Save the log-scale boxplot
  log_filename <- paste0("src/tmp_figures/eda_boxplot_group_", i, "_log.png")
  ggsave(log_filename, log_boxplot, width = 10, height = 8)
}

# 3. Perform PCA to explore overall patterns in the data
# Extract cytokine columns for PCA
cytokine_columns <- intersect(colnames(cytokine_data_wide), unique_cytokines)

# Prepare data for PCA
pca_data <- cytokine_data_wide[, cytokine_columns]

# Check if we have any missing values
if (any(is.na(pca_data))) {
  # Impute missing values with column means
  for (col in colnames(pca_data)) {
    pca_data[[col]][is.na(pca_data[[col]])] <- mean(pca_data[[col]], na.rm = TRUE)
  }
}

# Log transform data for PCA (adding small constant to handle zeros)
pca_data_log <- log1p(as.matrix(pca_data))

# Perform PCA
pca_result <- prcomp(pca_data_log, scale. = TRUE)

# Extract scores
pca_scores <- as.data.frame(pca_result$x)

# Add metadata
pca_scores$genotype <- cytokine_data_wide$genotype
pca_scores$time <- cytokine_data_wide$time
pca_scores$mouse_id <- cytokine_data_wide$mouse_id

# Calculate variance explained
variance <- pca_result$sdev^2
prop_variance <- variance / sum(variance)
cum_variance <- cumsum(prop_variance)

# Create data frame for plotting
variance_df <- data.frame(
  PC = factor(paste0("PC", 1:length(prop_variance)), 
              levels = paste0("PC", 1:length(prop_variance))),
  Proportion = prop_variance,
  Cumulative = cum_variance
)

# Save PCA results
saveRDS(pca_result, "src/tmp_data/eda_pca_result.rds")
saveRDS(pca_scores, "src/tmp_data/eda_pca_scores.rds")

# Plot variance explained (scree plot)
scree_plot <- ggplot(variance_df[1:min(10, nrow(variance_df)),], aes(x = PC, y = Proportion * 100)) +
  geom_bar(stat = "identity", fill = "#4D6D8E") +
  geom_line(aes(y = Cumulative * 100, group = 1), color = "#7AA661") +
  geom_point(aes(y = Cumulative * 100), color = "#7AA661", size = 3) +
  theme_bw() +
  labs(title = "PCA Variance Explained",
       x = "Principal Component",
       y = "Percentage of Variance Explained") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save scree plot
ggsave("src/tmp_figures/eda_pca_scree_plot.png", scree_plot, width = 8, height = 6)

# Plot first two principal components
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = genotype, shape = time)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
  theme_bw() +
  labs(title = "PCA of Cytokine Profiles",
       x = paste0("PC1 (", round(prop_variance[1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(prop_variance[2] * 100, 1), "% variance)"))

# Save PCA plot
ggsave("src/tmp_figures/eda_pca_plot.png", pca_plot, width = 8, height = 6)

# Create biplot of PC1 vs PC2
loadings <- as.data.frame(pca_result$rotation)
loadings$cytokine <- rownames(loadings)

# Scale loadings for visualization
max_score <- max(abs(c(pca_scores$PC1, pca_scores$PC2)))
scale_factor <- max_score / max(abs(c(loadings$PC1, loadings$PC2))) * 0.8

biplot_data <- loadings %>%
  select(cytokine, PC1, PC2) %>%
  mutate(
    PC1 = PC1 * scale_factor,
    PC2 = PC2 * scale_factor
  )

# Create biplot
biplot <- ggplot() +
  # Add points
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = genotype, shape = time), 
             size = 3, alpha = 0.7) +
  # Add arrows
  geom_segment(data = biplot_data, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#807E7D") +
  # Add labels
  geom_text(data = biplot_data, 
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = cytokine),
            size = 3) +
  scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
  theme_bw() +
  labs(title = "PCA Biplot of Cytokine Data",
       x = paste0("PC1 (", round(prop_variance[1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(prop_variance[2] * 100, 1), "% variance)"))

# Save biplot
ggsave("src/tmp_figures/eda_pca_biplot.png", biplot, width = 10, height = 8)

# 4. Create correlation heatmap between cytokines
# Calculate correlation matrix
cor_matrix <- cor(pca_data_log, method = "spearman", use = "pairwise.complete.obs")

# Save correlation matrix
saveRDS(cor_matrix, "src/tmp_data/eda_correlation_matrix.rds")

# Define colors for the heatmap
heatmap_colors <- colorRampPalette(c("#619ca6", "#FFFFFF", "#a6617a"))(100)

# Create correlation heatmap
pheatmap::pheatmap(
  cor_matrix,
  color = heatmap_colors,
  display_numbers = FALSE,
  fontsize = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Correlation Heatmap of Cytokines",
  filename = "src/tmp_figures/eda_correlation_heatmap.png",
  width = 10,
  height = 8
)

# Create correlation heatmap with values
pheatmap::pheatmap(
  cor_matrix,
  color = heatmap_colors,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 6,
  fontsize = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Correlation Heatmap of Cytokines (with values)",
  filename = "src/tmp_figures/eda_correlation_heatmap_with_values.png",
  width = 10,
  height = 8
)

# Create a heatmap of cytokine concentrations across samples
# Prepare data for the heatmap
sample_labels <- paste(
  cytokine_data_wide$genotype, 
  cytokine_data_wide$time, 
  cytokine_data_wide$mouse_number, 
  sep = "_"
)

# Create heatmap with annotations
# Create annotation data frame
sample_annotation <- data.frame(
  Genotype = cytokine_data_wide$genotype,
  Time = cytokine_data_wide$time,
  row.names = sample_labels
)

# Define annotation colors
annotation_colors <- list(
  Genotype = c(wt = "#4D6D8E", ko = "#7AA661"),
  Time = c("0 h" = "#807E7D", "2 h" = "#8e6e4d", "7 h" = "#7b4d8e")
)

# Generate heatmap of cytokine concentrations
heatmap_data <- as.matrix(pca_data_log)
rownames(heatmap_data) <- sample_labels

pheatmap::pheatmap(
  heatmap_data,
  annotation_row = sample_annotation,
  annotation_colors = annotation_colors,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize = 8,
  main = "Heatmap of Log-Transformed Cytokine Concentrations",
  filename = "src/tmp_figures/eda_concentration_heatmap.png",
  width = 12,
  height = 10
)

# Summary of saved objects
cat("\nSaved objects in exploratory data analysis:\n")
cat("1. eda_summary_stats.rds - Summary statistics for each cytokine by genotype and time\n")
cat("2. eda_pca_result.rds - PCA results object\n")
cat("3. eda_pca_scores.rds - PCA scores for each sample\n")
cat("4. eda_correlation_matrix.rds - Correlation matrix between cytokines\n")

cat("\nSaved figures in exploratory data analysis:\n")
cat("1. eda_boxplot_group_*.png - Boxplots of cytokine distributions by group\n")
cat("2. eda_boxplot_group_*_log.png - Log-scale boxplots of cytokine distributions\n")
cat("3. eda_pca_scree_plot.png - Scree plot of PCA variance explained\n")
cat("4. eda_pca_plot.png - PCA plot (PC1 vs PC2)\n")
cat("5. eda_pca_biplot.png - PCA biplot showing cytokine contributions\n")
cat("6. eda_correlation_heatmap.png - Correlation heatmap between cytokines\n")
cat("7. eda_correlation_heatmap_with_values.png - Correlation heatmap with values\n")
cat("8. eda_concentration_heatmap.png - Heatmap of cytokine concentrations\n")
