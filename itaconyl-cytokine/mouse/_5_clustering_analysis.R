
# Load required libraries
library(stats)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
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

# Load necessary data for clustering analysis
cytokine_data_wide <- load_data_safely("src/data/cytokine_data_wide.rds", 
                                       "Error loading wide-format cytokine data")
cytokine_data_long <- load_data_safely("src/data/cytokine_data_long.rds", 
                                       "Error loading long-format cytokine data")
cytokine_names <- load_data_safely("src/data/cytokine_names.rds", 
                                   "Error loading cytokine names")
time_course_summary <- load_data_safely("src/data/time_course_summary.rds", 
                                        "Error loading time course summary data")

# Check if data was loaded correctly
if (is.null(cytokine_data_wide) || is.null(cytokine_data_long) || is.null(cytokine_names)) {
  stop("Failed to load required data files for clustering analysis")
}

# If time_course_summary is not available, create it from cytokine_data_long
if (is.null(time_course_summary)) {
  time_course_summary <- cytokine_data_long %>%
    dplyr::group_by(genotype, time, cytokine) %>%
    dplyr::summarise(
      mean_concentration = mean(concentration, na.rm = TRUE),
      sd_concentration = sd(concentration, na.rm = TRUE),
      n = dplyr::n(),
      se_concentration = sd_concentration / sqrt(n),
      .groups = "drop"
    )
}

# Display basic information about the data
cat("Data dimensions:\n")
cat("- Number of cytokines:", length(cytokine_names), "\n")
cat("- Time points:", paste(levels(cytokine_data_wide$time), collapse = ", "), "\n\n")

# 1. Prepare data matrix of temporal profiles for WT mice
# ------------------------------------------------------

# Filter for WT mice and reshape data for clustering
wt_profiles <- time_course_summary %>%
  dplyr::filter(genotype == "wt") %>%
  dplyr::select(time, cytokine, mean_concentration) %>%
  tidyr::pivot_wider(
    names_from = time,
    values_from = mean_concentration
  )

# Create a matrix for clustering
wt_matrix <- as.matrix(wt_profiles[, -1])  # Remove cytokine column
rownames(wt_matrix) <- wt_profiles$cytokine

# Check data structure before clustering
cat("Temporal profile matrix for WT mice (first 5 rows):\n")
print(knitr::kable(head(wt_matrix, 5)))

# Log transform to better handle large concentration ranges
wt_matrix_log <- log1p(wt_matrix)  # log(x+1) transformation

# Scale the data (z-score normalization) for clustering
wt_matrix_scaled <- scale(wt_matrix_log)

# Save prepared matrices
saveRDS(wt_matrix, "src/tmp_data/wt_temporal_matrix.rds")
saveRDS(wt_matrix_log, "src/tmp_data/wt_temporal_matrix_log.rds")
saveRDS(wt_matrix_scaled, "src/tmp_data/wt_temporal_matrix_scaled.rds")

# 2. Perform hierarchical clustering of cytokines
# ----------------------------------------------

# Calculate distance matrix
dist_matrix <- stats::dist(wt_matrix_scaled, method = "euclidean")

# Perform hierarchical clustering using different linkage methods
hc_methods <- c("complete", "average", "ward.D2")
hclust_results <- list()

for (method in hc_methods) {
  hclust_results[[method]] <- stats::hclust(dist_matrix, method = method)
}

# Save clustering results
saveRDS(hclust_results, "src/tmp_data/wt_hclust_results.rds")

# Visualize dendrograms for each method
for (method in hc_methods) {
  png(
    filename = paste0("src/tmp_figures/dendrogram_", method, ".png"),
    width = 8, height = 6, units = "in", res = 300
  )
  
  plot(hclust_results[[method]], 
       main = paste("Hierarchical Clustering of Cytokines (", method, ")"),
       xlab = "Cytokines", 
       ylab = "Height",
       sub = "Based on temporal expression patterns in WT mice")
  
  dev.off()
}

# 3. Determine optimal number of clusters
# --------------------------------------

# Define a simplified function to evaluate cluster quality
evaluate_cluster_quality <- function(k, hc_result) {
  # Cut tree to get k clusters
  clusters <- stats::cutree(hc_result, k = k)
  
  # Calculate total within-cluster sum of squares
  within_ss <- 0
  
  # Get cluster means
  cluster_means <- matrix(0, nrow = k, ncol = ncol(wt_matrix_scaled))
  
  for (i in 1:k) {
    if (sum(clusters == i) > 0) {
      cluster_means[i,] <- colMeans(wt_matrix_scaled[clusters == i, , drop = FALSE])
      within_ss <- within_ss + sum(dist(rbind(wt_matrix_scaled[clusters == i, , drop = FALSE], 
                                            cluster_means[i, , drop = FALSE]))^2)
    }
  }
  
  # Calculate between-cluster sum of squares
  between_ss <- 0
  overall_mean <- colMeans(wt_matrix_scaled)
  
  for (i in 1:k) {
    if (sum(clusters == i) > 0) {
      between_ss <- between_ss + sum(clusters == i) * 
        sum((cluster_means[i,] - overall_mean)^2)
    }
  }
  
  # Return ratio of between to within cluster variance (higher is better)
  if (within_ss > 0) {
    return(between_ss / within_ss)
  } else {
    return(0)
  }
}

# Evaluate different numbers of clusters
max_k <- min(8, nrow(wt_matrix_scaled) - 1)  # Maximum number of clusters to consider
cluster_quality <- numeric(max_k - 1)

set.seed(123)  # For reproducibility
for (k in 2:max_k) {
  cluster_quality[k-1] <- evaluate_cluster_quality(k, hclust_results[["ward.D2"]])
}

# Create data frame for plotting
quality_df <- data.frame(
  k = 2:max_k,
  quality = cluster_quality
)

# Plot cluster quality metric
quality_plot <- ggplot2::ggplot(quality_df, ggplot2::aes(x = k, y = quality)) +
  ggplot2::geom_line(color = "#4D6D8E", size = 1) +
  ggplot2::geom_point(color = "#4D6D8E", size = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Cluster Quality for Different Numbers of Clusters",
    x = "Number of Clusters",
    y = "Between/Within Cluster Variance Ratio"
  ) +
  ggplot2::scale_x_continuous(breaks = 2:max_k)

# Save quality plot
ggplot2::ggsave("src/tmp_figures/cluster_quality.png", quality_plot, width = 8, height = 6)

# Determine optimal number of clusters based on quality metric
# Look for "elbow" in the plot, or choose highest reasonable value
quality_increase <- c(0, diff(cluster_quality))
candidate_k <- which(quality_increase == max(quality_increase)) + 1
if (candidate_k < 2 || candidate_k > 5) {
  optimal_k <- 4  # Fallback to a reasonable default
} else {
  optimal_k <- candidate_k
}

cat("Selected number of clusters:", optimal_k, "\n")

# Allow for visual inspection by creating color-coded dendrograms for a few candidate k values
candidate_k_values <- c(3, 4, 5)
for (k in candidate_k_values) {
  png(
    filename = paste0("src/tmp_figures/colored_dendrogram_k", k, ".png"),
    width = 10, height = 7, units = "in", res = 300
  )
  
  # Plot dendrogram with k clusters highlighted
  plot(hclust_results[["ward.D2"]], 
       main = paste("Dendrogram with", k, "Clusters"),
       xlab = "Cytokines",
       ylab = "Height")
  
  # Add colored rectangles for the clusters
  stats::rect.hclust(hclust_results[["ward.D2"]], k = k, 
                     border = c("#4D6D8E", "#7AA661", "#807E7D", "#8e6e4d", "#7b4d8e")[1:k])
  
  dev.off()
}

# Cut the tree to get cluster assignments with the optimal number of clusters
cluster_assignments <- stats::cutree(hclust_results[["ward.D2"]], k = optimal_k)
cluster_df <- data.frame(
  cytokine = names(cluster_assignments),
  cluster = as.factor(cluster_assignments)
)

# Save cluster assignments
saveRDS(cluster_df, "src/tmp_data/cytokine_clusters.rds")

# 4. Create heatmap with clusters
# ------------------------------

# Prepare annotation data
annotation_row <- data.frame(
  Cluster = cluster_df$cluster,
  row.names = cluster_df$cytokine
)

# Define colors for clusters - ensure we have enough colors
cluster_colors <- c("#4D6D8E", "#7AA661", "#807E7D", "#8e6e4d", "#7b4d8e")
cluster_colors <- cluster_colors[1:length(levels(annotation_row$Cluster))]  # Limit to actual number of clusters

# Create named vector for annotation colors
ann_colors <- list(
  Cluster = cluster_colors
)
names(ann_colors$Cluster) <- levels(annotation_row$Cluster)

# Create heatmap of temporal profiles with cluster annotation
pheatmap::pheatmap(
  wt_matrix_log,
  color = colorRampPalette(c("#619ca6", "#FFFFFF", "#a6617a"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_method = "ward.D2",
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  main = "Clustering of Cytokine Temporal Profiles in WT Mice",
  fontsize = 10,
  filename = "src/tmp_figures/cytokine_cluster_heatmap.png",
  width = 10,
  height = 8
)

# 5. Visualize cluster-specific temporal patterns
# ---------------------------------------------

# Join cluster assignments with time course data
wt_time_course_with_clusters <- time_course_summary %>%
  dplyr::filter(genotype == "wt") %>%
  dplyr::left_join(cluster_df, by = "cytokine") %>%
  dplyr::mutate(time_numeric = as.numeric(gsub(" h", "", time)))

# Calculate mean profile for each cluster
cluster_means <- wt_time_course_with_clusters %>%
  dplyr::group_by(cluster, time, time_numeric) %>%
  dplyr::summarise(
    mean_value = mean(mean_concentration, na.rm = TRUE),
    se_value = sd(mean_concentration, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Plot mean profiles for each cluster
cluster_profiles_plot <- ggplot2::ggplot(cluster_means, 
                                        ggplot2::aes(x = time_numeric, y = mean_value, 
                                                    color = cluster, group = cluster)) +
  ggplot2::geom_line(size = 1.2) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - se_value, 
                                     ymax = mean_value + se_value), 
                        width = 0.2) +
  ggplot2::scale_color_manual(values = cluster_colors) +
  ggplot2::scale_x_continuous(breaks = sort(unique(wt_time_course_with_clusters$time_numeric)), 
                             labels = levels(time_course_summary$time)) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Mean Temporal Profiles by Cluster",
    x = "Time after LPS challenge",
    y = "Mean Concentration (pg/ml)",
    color = "Cluster"
  )

# Save cluster profiles plot
ggplot2::ggsave("src/tmp_figures/cluster_mean_profiles.png", cluster_profiles_plot, width = 9, height = 7)

# Create log-scaled version for better visualization
cluster_profiles_plot_log <- cluster_profiles_plot +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Mean Temporal Profiles by Cluster (Log Scale)",
    y = "Mean Concentration (pg/ml, log scale)"
  )

# Save log-scaled cluster profiles plot
ggplot2::ggsave("src/tmp_figures/cluster_mean_profiles_log.png", cluster_profiles_plot_log, width = 9, height = 7)

# 6. Compare cluster patterns between WT and KO
# -------------------------------------------

# Get KO temporal profiles as well
ko_time_course_with_clusters <- time_course_summary %>%
  dplyr::filter(genotype == "ko") %>%
  dplyr::left_join(cluster_df, by = "cytokine") %>%
  dplyr::mutate(time_numeric = as.numeric(gsub(" h", "", time)))

# Combine WT and KO data for comparison
comparison_data <- rbind(
  wt_time_course_with_clusters %>% dplyr::mutate(genotype_label = "WT"),
  ko_time_course_with_clusters %>% dplyr::mutate(genotype_label = "KO")
)

# Save comparison data
saveRDS(comparison_data, "src/tmp_data/wt_ko_cluster_comparison.rds")

# Create comparison plots for each cluster
for (c in levels(cluster_df$cluster)) {
  # Filter data for this cluster
  cluster_data <- comparison_data %>%
    dplyr::filter(cluster == c)
  
  # Create plot
  p <- ggplot2::ggplot(cluster_data, 
                       ggplot2::aes(x = time_numeric, y = mean_concentration, 
                                    color = genotype)) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(genotype, cytokine)), size = 0.7, alpha = 0.7) +
    ggplot2::geom_point(size = 1, alpha = 0.7) +
    ggplot2::facet_wrap(~ cytokine, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661")) +
    ggplot2::scale_x_continuous(breaks = sort(unique(comparison_data$time_numeric)), 
                               labels = levels(time_course_summary$time)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Cluster", c, "Cytokines: WT vs KO Comparison"),
      x = "Time after LPS challenge",
      y = "Concentration (pg/ml)"
    ) +
    ggplot2::theme(legend.position = "bottom")
  
  # Save plot
  ggplot2::ggsave(
    paste0("src/tmp_figures/cluster_", c, "_wt_ko_comparison.png"),
    p, width = 10, height = 8
  )
  
  # Create log-scaled version
  p_log <- p +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = paste("Cluster", c, "Cytokines: WT vs KO Comparison (Log Scale)"),
      y = "Concentration (pg/ml, log scale)"
    )
  
  # Save log-scaled plot
  ggplot2::ggsave(
    paste0("src/tmp_figures/cluster_", c, "_wt_ko_comparison_log.png"),
    p_log, width = 10, height = 8
  )
}

# Create a summary plot comparing WT and KO mean profiles by cluster
cluster_comparison_means <- comparison_data %>%
  dplyr::group_by(genotype, genotype_label, cluster, time, time_numeric) %>%
  dplyr::summarise(
    mean_value = mean(mean_concentration, na.rm = TRUE),
    se_value = sd(mean_concentration, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Plot comparison of mean profiles
cluster_comparison_plot <- ggplot2::ggplot(cluster_comparison_means, 
                                         ggplot2::aes(x = time_numeric, y = mean_value, 
                                                     color = genotype, group = genotype)) +
  ggplot2::geom_line(size = 1.2) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - se_value, 
                                      ymax = mean_value + se_value), 
                         width = 0.2) +
  ggplot2::facet_wrap(~ cluster, scales = "free_y") +
  ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661")) +
  ggplot2::scale_x_continuous(breaks = sort(unique(comparison_data$time_numeric)), 
                              labels = levels(time_course_summary$time)) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Comparison of WT and KO Mean Profiles by Cluster",
    x = "Time after LPS challenge",
    y = "Mean Concentration (pg/ml)",
    color = "Genotype"
  )

# Save comparison plot
ggplot2::ggsave("src/tmp_figures/cluster_wt_ko_mean_comparison.png", cluster_comparison_plot, width = 10, height = 8)

# Create log-scaled version
cluster_comparison_plot_log <- cluster_comparison_plot +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Comparison of WT and KO Mean Profiles by Cluster (Log Scale)",
    y = "Mean Concentration (pg/ml, log scale)"
  )

# Save log-scaled comparison plot
ggplot2::ggsave("src/tmp_figures/cluster_wt_ko_mean_comparison_log.png", cluster_comparison_plot_log, width = 10, height = 8)

# 7. Calculate fold changes from baseline for each cluster
# ---------------------------------------------------

# Calculate fold changes from baseline for each genotype
fold_change_data <- time_course_summary %>%
  dplyr::group_by(genotype, cytokine) %>%
  dplyr::mutate(
    fold_change = mean_concentration / (mean_concentration[time == "0 h"] + 0.01),
    log2_fold_change = log2((mean_concentration + 0.01) / (mean_concentration[time == "0 h"] + 0.01))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(time != "0 h")  # Exclude baseline

# Add cluster information
fold_change_with_clusters <- fold_change_data %>%
  dplyr::left_join(cluster_df, by = "cytokine")

# Save fold change data with clusters
saveRDS(fold_change_with_clusters, "src/tmp_data/fold_change_with_clusters.rds")

# Reshape for heatmap
fold_change_matrix <- fold_change_with_clusters %>%
  dplyr::select(cytokine, genotype, time, log2_fold_change, cluster) %>%
  tidyr::unite("genotype_time", genotype, time, sep = "_") %>%
  tidyr::pivot_wider(
    names_from = genotype_time,
    values_from = log2_fold_change
  )

# Save fold change matrix with clusters
saveRDS(fold_change_matrix, "src/tmp_data/fold_change_matrix_with_clusters.rds")

# Get matrix for heatmap
matrix_cols <- grep("^wt_|^ko_", names(fold_change_matrix), value = TRUE)
heatmap_matrix <- as.matrix(fold_change_matrix[, matrix_cols])
rownames(heatmap_matrix) <- fold_change_matrix$cytokine

# Create annotation for heatmap
heatmap_annotation <- data.frame(
  Cluster = fold_change_matrix$cluster,
  row.names = fold_change_matrix$cytokine
)

# Create heatmap of fold changes with cluster annotation
pheatmap::pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("#619ca6", "#FFFFFF", "#a6617a"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = heatmap_annotation,
  annotation_colors = ann_colors,
  main = "Log2 Fold Changes from Baseline by Cluster",
  fontsize = 10,
  filename = "src/tmp_figures/fold_change_cluster_heatmap.png",
  width = 10,
  height = 8
)

# 8. Create cluster summary statistics
# ----------------------------------

# Calculate summary statistics for each cluster
cluster_summary <- cluster_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    count = dplyr::n(),
    cytokines = paste(cytokine, collapse = ", "),
    .groups = "drop"
  )

# Load response patterns if available
response_patterns <- load_data_safely("src/data/response_patterns.rds",
                                     "Error loading response patterns data")

# If response_patterns is available, add pattern information to cluster summary
if (!is.null(response_patterns)) {
  # Join cluster assignments with response patterns
  cluster_patterns <- cluster_df %>%
    dplyr::left_join(response_patterns, by = "cytokine")
  
  # Summarize patterns by cluster
  pattern_summary <- cluster_patterns %>%
    dplyr::group_by(cluster, overall_pattern) %>%
    dplyr::summarise(
      count = dplyr::n(),
      cytokines = paste(cytokine, collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(cluster, desc(count))
  
  # Calculate dominant pattern for each cluster
  cluster_dominant_patterns <- pattern_summary %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(cluster, dominant_pattern = overall_pattern)
  
  # Add dominant pattern to cluster summary
  cluster_summary <- cluster_summary %>%
    dplyr::left_join(cluster_dominant_patterns, by = "cluster")
  
  # Save pattern summary
  saveRDS(pattern_summary, "src/tmp_data/cluster_pattern_summary.rds")
}

# Save cluster summary
saveRDS(cluster_summary, "src/tmp_data/cluster_summary.rds")

# Display cluster summary
cat("\nCluster Summary:\n")
print(knitr::kable(cluster_summary))

# Summary of files created
cat("\nFiles saved:\n")
cat("1. wt_temporal_matrix.rds - Original concentration matrix for WT mice\n")
cat("2. wt_temporal_matrix_log.rds - Log-transformed concentration matrix\n")
cat("3. wt_temporal_matrix_scaled.rds - Z-score normalized matrix for clustering\n")
cat("4. wt_hclust_results.rds - Hierarchical clustering results with different methods\n")
cat("5. cytokine_clusters.rds - Cluster assignments for each cytokine\n")
cat("6. wt_ko_cluster_comparison.rds - Comparison data between WT and KO by cluster\n")
cat("7. fold_change_with_clusters.rds - Fold changes from baseline with cluster information\n")
cat("8. fold_change_matrix_with_clusters.rds - Fold changes from baseline in matrix format\n")
cat("9. cluster_summary.rds - Summary statistics for each cluster\n")
if (!is.null(response_patterns)) {
  cat("10. cluster_pattern_summary.rds - Response patterns by cluster\n")
}

cat("\nFigures saved:\n")
cat("1. dendrogram_*.png - Hierarchical clustering dendrograms with different linkage methods\n")
cat("2. cluster_quality.png - Cluster quality metric for different numbers of clusters\n")
cat("3. colored_dendrogram_k*.png - Dendrograms with colored branches for different k values\n")
cat("4. cytokine_cluster_heatmap.png - Heatmap of temporal profiles with cluster annotations\n")
cat("5. cluster_mean_profiles.png - Mean profiles for each cluster\n")
cat("6. cluster_mean_profiles_log.png - Log-scaled mean profiles\n")
cat("7. cluster_*_wt_ko_comparison.png - Comparison between WT and KO for each cluster\n")
cat("8. cluster_*_wt_ko_comparison_log.png - Log-scaled comparison between WT and KO\n")
cat("9. cluster_wt_ko_mean_comparison.png - Mean profiles comparison by cluster\n")
cat("10. cluster_wt_ko_mean_comparison_log.png - Log-scaled mean profiles comparison\n")
cat("11. fold_change_cluster_heatmap.png - Heatmap of fold changes with cluster annotations\n")
