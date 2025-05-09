
# Load required libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(dplyr)
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

# Load necessary data files for visualization
cat("Loading data for visualization...\n")

# Load data for volcano plots
diff_analysis_volcano_data <- load_data_safely("src/data/diff_analysis_volcano_data.rds", 
                                               "Error loading volcano plot data")

# Load data for temporal patterns
cytokine_clusters <- load_data_safely("src/data/cytokine_clusters.rds",
                                     "Error loading cytokine cluster data")
time_course_summary <- load_data_safely("src/data/time_course_summary.rds",
                                       "Error loading time course data")
response_patterns <- load_data_safely("src/data/response_patterns.rds",
                                     "Error loading response patterns data")
wt_ko_cluster_comparison <- load_data_safely("src/data/wt_ko_cluster_comparison.rds",
                                            "Error loading cluster comparison data")
fold_change_with_clusters <- load_data_safely("src/data/fold_change_with_clusters.rds",
                                             "Error loading fold change data with clusters")

# Load other necessary data
cytokine_names <- load_data_safely("src/data/cytokine_names.rds",
                                  "Error loading cytokine names")
diff_analysis_top_tables <- load_data_safely("src/data/diff_analysis_top_tables.rds",
                                            "Error loading differential analysis results")

# Check if essential data is available
if (is.null(diff_analysis_volcano_data)) {
  stop("Volcano plot data not available. Unable to create visualizations.")
}

# 1. Create volcano plots for each time point comparison
# ------------------------------------------------------
cat("Creating volcano plots...\n")

create_volcano_plot <- function(data, time_point, p_threshold = 0.05, fc_threshold = 0.5) {
  # Add significance indicators
  plot_data <- data %>%
    dplyr::mutate(
      significant = adj.P.Val < p_threshold & abs(logFC) > fc_threshold,
      label = ifelse(significant | (adj.P.Val < p_threshold * 2 & abs(logFC) > fc_threshold * 1.5), 
                    cytokine, NA)
    )
  
  # Create volcano plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "#807E7D") +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "#807E7D") +
    ggplot2::scale_color_manual(values = c("FALSE" = "#807E7D", "TRUE" = "#7AA661")) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      box.padding = 0.5,
      max.overlaps = 20,
      size = 3.5
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Volcano Plot -", time_point),
      subtitle = paste("KO vs WT Comparison at", time_point),
      x = "Log2 Fold Change (KO/WT)",
      y = "-Log10(Adjusted P-value)"
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}

# Generate individual volcano plots for each time point
volcano_plots <- list()
for (tp in names(diff_analysis_volcano_data)) {
  volcano_plots[[tp]] <- create_volcano_plot(diff_analysis_volcano_data[[tp]], tp)
  
  # Save individual volcano plots
  ggplot2::ggsave(
    filename = paste0("src/tmp_figures/volcano_plot_", gsub(" ", "_", tp), ".png"),
    plot = volcano_plots[[tp]],
    width = 8,
    height = 7,
    dpi = 300
  )
}

# Create combined volcano plot using patchwork
if (length(volcano_plots) > 1) {
  combined_volcanoes <- patchwork::wrap_plots(volcano_plots, ncol = 2)
  
  # Save combined volcano plots
  ggplot2::ggsave(
    filename = "src/tmp_figures/combined_volcano_plots.png",
    plot = combined_volcanoes,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# 2. Create heatmap with hierarchical clustering of temporal patterns
# ------------------------------------------------------------------
cat("Creating heatmaps with hierarchical clustering...\n")

# Check if fold change data with clusters is available
if (!is.null(fold_change_with_clusters)) {
  # Prepare data matrix for heatmap
  fold_change_matrix <- fold_change_with_clusters %>%
    dplyr::select(genotype, time, cytokine, log2_fold_change, cluster) %>%
    dplyr::mutate(genotype_time = paste(genotype, time, sep = "_")) %>%
    tidyr::pivot_wider(
      id_cols = c(cytokine, cluster),
      names_from = genotype_time,
      values_from = log2_fold_change
    )
  
  # Create matrix for heatmap
  fc_cols <- grep("^wt_|^ko_", names(fold_change_matrix), value = TRUE)
  heatmap_matrix <- as.matrix(fold_change_matrix[, fc_cols])
  rownames(heatmap_matrix) <- fold_change_matrix$cytokine
  
  # Create annotation data
  annotation_row <- data.frame(
    Cluster = fold_change_matrix$cluster,
    row.names = fold_change_matrix$cytokine
  )
  
  # Add response pattern annotation if available
  if (!is.null(response_patterns)) {
    # Join response patterns with cytokines
    annotation_row$Response <- factor(response_patterns$overall_pattern[
      match(rownames(annotation_row), response_patterns$cytokine)
    ])
  }
  
  # Define colors for annotations
  cluster_colors <- c("#4D6D8E", "#7AA661", "#807E7D", "#8e6e4d", "#7b4d8e")
  
  # Create named vector for cluster colors
  ann_colors <- list(
    Cluster = setNames(cluster_colors[1:length(unique(annotation_row$Cluster))], 
                     sort(unique(annotation_row$Cluster)))
  )
  
  # Add response pattern colors if available
  if ("Response" %in% colnames(annotation_row)) {
    response_colors <- c(
      "Consistently Enhanced in KO" = "#7AA661",
      "Consistently Reduced in KO" = "#4D6D8E",
      "Early Enhanced in KO" = "#8e6e4d",
      "Late Enhanced in KO" = "#7b4d8e",
      "Early Reduced in KO" = "#619ca6",
      "Late Reduced in KO" = "#a6617a",
      "Similar Response" = "#807E7D"
    )
    
    # Filter to only include observed response patterns
    observed_patterns <- levels(annotation_row$Response)
    if (length(observed_patterns) > 0) {
      ann_colors$Response <- response_colors[observed_patterns]
    }
  }
  
  # Create heatmap with annotations
  pheatmap::pheatmap(
    heatmap_matrix,
    color = colorRampPalette(c("#619ca6", "#FFFFFF", "#a6617a"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    fontsize = 8,
    cellwidth = 15,
    cellheight = 10,
    main = "Log2 Fold Changes from Baseline by Cluster and Genotype",
    filename = "src/tmp_figures/fold_change_heatmap_clustered.png",
    width = 10,
    height = 12
  )
  
  # Save fold change matrix with clusters for future reference
  saveRDS(fold_change_matrix, "src/tmp_data/viz_fold_change_matrix.rds")
}

# 3. Generate time course line plots for visualization by cluster
# --------------------------------------------------------------
cat("Creating time course line plots by cluster...\n")

if (!is.null(wt_ko_cluster_comparison) && !is.null(cytokine_clusters)) {
  # Create a function to generate time course plots for each cluster
  create_cluster_time_course <- function(data, cluster_id) {
    # Filter data for the specified cluster
    cluster_data <- data %>%
      dplyr::filter(cluster == cluster_id)
    
    # Count cytokines in this cluster
    cytokines_in_cluster <- unique(cluster_data$cytokine)
    cluster_title <- paste("Cluster", cluster_id, "-", 
                          length(cytokines_in_cluster), "cytokines")
    
    # Create the plot
    p <- ggplot2::ggplot(cluster_data, 
                        ggplot2::aes(x = time_numeric, y = mean_concentration, 
                                     color = genotype, group = interaction(genotype, cytokine))) +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.8) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::facet_wrap(~ cytokine, scales = "free_y", ncol = 3) +
      ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661"),
                                 labels = c("wt" = "WT", "ko" = "KO")) +
      ggplot2::scale_x_continuous(breaks = c(0, 2, 7),
                                 labels = c("0h", "2h", "7h")) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = cluster_title,
        x = "Time after LPS challenge",
        y = "Concentration (pg/ml)",
        color = "Genotype"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    
    return(p)
  }
  
  # Generate plots for each cluster
  clusters <- sort(unique(cytokine_clusters$cluster))
  cluster_plots <- list()
  
  for (c in clusters) {
    plot <- create_cluster_time_course(wt_ko_cluster_comparison, c)
    cluster_plots[[paste0("cluster_", c)]] <- plot
    
    # Save individual cluster plots
    ggplot2::ggsave(
      filename = paste0("src/tmp_figures/time_course_cluster_", c, ".png"),
      plot = plot,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    # Create log-scale version for better visualization of wide ranges
    log_plot <- plot +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        title = paste("Cluster", c, "- Log Scale"),
        y = "Concentration (pg/ml, log scale)"
      )
    
    # Save log-scale version
    ggplot2::ggsave(
      filename = paste0("src/tmp_figures/time_course_cluster_", c, "_log.png"),
      plot = log_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  # Create a summary plot showing mean profiles by cluster
  cluster_means <- wt_ko_cluster_comparison %>%
    dplyr::group_by(genotype, cluster, time, time_numeric) %>%
    dplyr::summarise(
      mean_concentration = mean(mean_concentration, na.rm = TRUE),
      se_concentration = sd(mean_concentration, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
  
  # Create mean profile plot
  mean_profile_plot <- ggplot2::ggplot(cluster_means, 
                                      ggplot2::aes(x = time_numeric, y = mean_concentration, 
                                                   color = genotype, group = genotype)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_concentration - se_concentration, 
                   ymax = mean_concentration + se_concentration),
      width = 0.2
    ) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661"),
                               labels = c("wt" = "WT", "ko" = "KO")) +
    ggplot2::scale_x_continuous(breaks = c(0, 2, 7),
                               labels = c("0h", "2h", "7h")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Mean Cytokine Response by Cluster",
      x = "Time after LPS challenge",
      y = "Mean Concentration (pg/ml)",
      color = "Genotype"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Save mean profile plot
  ggplot2::ggsave(
    filename = "src/tmp_figures/cluster_mean_profiles.png",
    plot = mean_profile_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Create log-scale version
  log_mean_profile_plot <- mean_profile_plot +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = "Mean Cytokine Response by Cluster (Log Scale)",
      y = "Mean Concentration (pg/ml, log scale)"
    )
  
  # Save log-scale mean profile plot
  ggplot2::ggsave(
    filename = "src/tmp_figures/cluster_mean_profiles_log.png",
    plot = log_mean_profile_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Save cluster means data
  saveRDS(cluster_means, "src/tmp_data/viz_cluster_means.rds")
}

# 4. Create faceted plots comparing WT and KO responses for key cytokines
# ----------------------------------------------------------------------
cat("Creating faceted comparison plots for key cytokines...\n")

if (!is.null(time_course_summary) && !is.null(response_patterns)) {
  # Identify key cytokines with interesting patterns
  interesting_patterns <- c(
    "Consistently Enhanced in KO",
    "Consistently Reduced in KO",
    "Early Enhanced in KO",
    "Late Enhanced in KO"
  )
  
  # Filter for cytokines with interesting patterns
  key_cytokines <- response_patterns %>%
    dplyr::filter(overall_pattern %in% interesting_patterns) %>%
    dplyr::arrange(overall_pattern) %>%
    dplyr::pull(cytokine)
  
  # If we have too many key cytokines, select a subset
  if (length(key_cytokines) > 12) {
    set.seed(123) # For reproducibility
    key_cytokines <- sample(key_cytokines, 12)
  }
  
  # Filter time course data for key cytokines
  key_cytokine_data <- time_course_summary %>%
    dplyr::filter(cytokine %in% key_cytokines) %>%
    dplyr::left_join(
      response_patterns %>% dplyr::select(cytokine, overall_pattern),
      by = "cytokine"
    ) %>%
    dplyr::mutate(
      time_numeric = as.numeric(gsub(" h", "", time)),
      overall_pattern = factor(overall_pattern, levels = interesting_patterns)
    )
  
  # Create faceted plot for key cytokines
  key_cytokines_plot <- ggplot2::ggplot(key_cytokine_data, 
                                       ggplot2::aes(x = time_numeric, y = mean_concentration, 
                                                    color = genotype, group = genotype)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_concentration - se_concentration, 
                   ymax = mean_concentration + se_concentration),
      width = 0.2
    ) +
    ggplot2::facet_wrap(~ cytokine, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661"),
                               labels = c("wt" = "WT", "ko" = "KO")) +
    ggplot2::scale_x_continuous(breaks = c(0, 2, 7),
                               labels = c("0h", "2h", "7h")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Key Cytokines with Different Response Patterns",
      x = "Time after LPS challenge",
      y = "Concentration (pg/ml)",
      color = "Genotype"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "#807E7D"),
      strip.text = ggplot2::element_text(color = "white", face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Save key cytokines plot
  ggplot2::ggsave(
    filename = "src/tmp_figures/key_cytokines_comparison.png",
    plot = key_cytokines_plot,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  # Create log-scale version
  log_key_cytokines_plot <- key_cytokines_plot +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = "Key Cytokines with Different Response Patterns (Log Scale)",
      y = "Concentration (pg/ml, log scale)"
    )
  
  # Save log-scale key cytokines plot
  ggplot2::ggsave(
    filename = "src/tmp_figures/key_cytokines_comparison_log.png",
    plot = log_key_cytokines_plot,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  # Create a plot grouped by response pattern
  pattern_plot <- ggplot2::ggplot(key_cytokine_data, 
                                 ggplot2::aes(x = time_numeric, y = mean_concentration, 
                                             color = genotype, group = interaction(genotype, cytokine))) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(
      overall_pattern ~ cytokine, 
      scales = "free_y"
    ) +
    ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661"),
                               labels = c("wt" = "WT", "ko" = "KO")) +
    ggplot2::scale_x_continuous(breaks = c(0, 2, 7),
                               labels = c("0h", "2h", "7h")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Cytokine Response Patterns: WT vs KO",
      x = "Time after LPS challenge",
      y = "Concentration (pg/ml)",
      color = "Genotype"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background.y = ggplot2::element_rect(fill = "#8e6e4d"),
      strip.text.y = ggplot2::element_text(color = "white", face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Save pattern-grouped plot
  ggplot2::ggsave(
    filename = "src/tmp_figures/response_patterns_cytokines.png",
    plot = pattern_plot,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  # Save key cytokine data for future reference
  saveRDS(key_cytokine_data, "src/tmp_data/viz_key_cytokines.rds")
}

# Create a comprehensive multi-panel visualization
# -----------------------------------------------
cat("Creating comprehensive multi-panel visualization...\n")

# Check if we have all the required plots
if (exists("volcano_plots") && exists("mean_profile_plot") && exists("key_cytokines_plot")) {
  # Select a specific volcano plot (e.g., 2h) if available
  if ("2 h" %in% names(volcano_plots)) {
    # Create multi-panel visualization
    summary_plot <- (volcano_plots[["2 h"]] | mean_profile_plot) / key_cytokines_plot +
      patchwork::plot_annotation(
        title = "Summary of Cytokine Responses in SIRT4KO vs WT Mice",
        subtitle = "Differential expression, cluster patterns, and key cytokine responses",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12)
        )
      )
    
    # Save the combined plot
    ggplot2::ggsave(
      filename = "src/tmp_figures/comprehensive_cytokine_summary.png",
      plot = summary_plot,
      width = 16,
      height = 14,
      dpi = 300
    )
  }
}

cat("Visualization stage completed. All plots saved to src/tmp_figures/\n")
