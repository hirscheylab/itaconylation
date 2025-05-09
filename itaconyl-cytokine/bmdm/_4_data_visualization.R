
# Load required libraries
library(ggplot2)
library(pheatmap)
library(patchwork)
library(stats)
library(dplyr)
library(tidyr)

# -----------------------------------------------
# Import data
# -----------------------------------------------
# Import cytokine data and log2-transformed data
cytokine_data <- readRDS("src/data/cytokine_data.rds")
cytokine_log2 <- readRDS("src/data/cytokine_log2.rds")

# Import statistical analysis results
anova_summary <- readRDS("src/data/anova_summary.rds")
significant_cytokines <- readRDS("src/data/significant_cytokines.rds")

# Identify cytokine columns (all columns after protein_mg_ml)
cytokine_cols <- colnames(cytokine_data)[6:ncol(cytokine_data)]

# Create a long format of the log2-transformed data for easier plotting
cytokine_log2_long <- cytokine_log2 %>%
  tidyr::pivot_longer(
    cols = all_of(cytokine_cols),
    names_to = "Cytokine",
    values_to = "Log2_Value"
  )

# -----------------------------------------------
# 1. Box plots comparing WT vs KO for each treatment condition
# -----------------------------------------------

# Create separate boxplots for each treatment condition to improve readability
lps_levels <- levels(cytokine_data$LPS_treatment_hours)

for (lps in lps_levels) {
  # Skip the 'none' level which doesn't have data
  if (lps != "none") {
    # Filter data for this LPS treatment condition
    condition_data <- cytokine_log2_long %>%
      dplyr::filter(LPS_treatment_hours == lps) %>%
      # Remove NA rows
      dplyr::filter(!is.na(Genotype_clean))
    
    # Create boxplot
    p <- ggplot2::ggplot(condition_data, 
           aes(x = Cytokine, y = Log2_Value, fill = Genotype_clean)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      ggplot2::geom_point(aes(color = Genotype_clean), 
                 position = position_jitterdodge(jitter.width = 0.1), 
                 size = 1.5, alpha = 0.7) +
      ggplot2::scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
      ggplot2::scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)
      ) +
      ggplot2::labs(
        title = paste("Cytokine Levels at", lps, "Hours LPS Treatment"),
        x = "",
        y = "Log2(Cytokine Level + 0.01)",
        fill = "Genotype",
        color = "Genotype"
      )
    
    # Save plot
    ggplot2::ggsave(
      paste0("src/tmp_figures/boxplot_", lps, "_hours.png"),
      p,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
}

# Create a focused plot showing only significant cytokines with interaction effects
sig_cytokines <- significant_cytokines %>%
  dplyr::filter(Interaction_significant_adj) %>%
  dplyr::pull(Cytokine)

if (length(sig_cytokines) > 0) {
  sig_data <- cytokine_log2_long %>%
    dplyr::filter(Cytokine %in% sig_cytokines) %>%
    dplyr::filter(!is.na(Genotype_clean))
  
  p_sig <- ggplot2::ggplot(sig_data, 
         aes(x = LPS_treatment_hours, y = Log2_Value, fill = Genotype_clean)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_point(aes(color = Genotype_clean), 
               position = position_jitterdodge(jitter.width = 0.2), 
               size = 2, alpha = 0.7) +
    ggplot2::facet_wrap(~ Cytokine, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
    ggplot2::scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Significantly Different Cytokines with Interaction Effects",
      x = "LPS Treatment (hours)",
      y = "Log2(Cytokine Level + 0.01)",
      fill = "Genotype",
      color = "Genotype"
    )
  
  ggplot2::ggsave(
    "src/tmp_figures/significant_cytokines_boxplot.png",
    p_sig,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# -----------------------------------------------
# 2. Heatmap of all cytokines across samples
# -----------------------------------------------

# Prepare data for heatmap - create a matrix of log2 values
heatmap_data <- as.matrix(cytokine_log2[, cytokine_cols])
rownames(heatmap_data) <- cytokine_log2$ID

# Create annotation dataframe
annotation_df <- cytokine_log2 %>%
  dplyr::select(ID, Genotype_clean, LPS_treatment_hours) %>%
  dplyr::filter(!is.na(Genotype_clean)) # Remove NA rows

# Make annotation a data frame
annotation_df <- as.data.frame(annotation_df)
rownames(annotation_df) <- annotation_df$ID
annotation_df$ID <- NULL

# Filter heatmap data to match annotation
heatmap_data <- heatmap_data[rownames(annotation_df), ]

# Define colors for annotation
annotation_colors <- list(
  Genotype_clean = c(KO = "#7AA661", WT = "#4D6D8E"),
  LPS_treatment_hours = c(
    "0" = "#807E7D",
    "6" = "#8e6e4d",
    "24" = "#7b4d8e"
  )
)

# Generate heatmap
pheatmap::pheatmap(
  heatmap_data,
  annotation_row = annotation_df,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  scale = "column", # Scale by column (cytokine)
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Heatmap of Cytokine Expression",
  filename = "src/tmp_figures/cytokine_heatmap.png",
  width = 10,
  height = 8
)

# Create a focused heatmap with significant cytokines
if (length(sig_cytokines) > 0) {
  sig_heatmap_data <- heatmap_data[, sig_cytokines, drop = FALSE]
  
  pheatmap::pheatmap(
    sig_heatmap_data,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    scale = "column",
    fontsize_row = 8,
    fontsize_col = 10,
    main = "Heatmap of Significantly Different Cytokines",
    filename = "src/tmp_figures/significant_cytokines_heatmap.png",
    width = 8,
    height = 8
  )
}

# -----------------------------------------------
# 3. PCA to visualize overall sample relationships
# -----------------------------------------------

# Remove missing values for PCA
complete_samples <- cytokine_log2 %>%
  dplyr::filter(!is.na(Genotype_clean))

# Perform PCA on log2-transformed data
pca_data <- complete_samples[, cytokine_cols]
pca_result <- stats::prcomp(pca_data, scale = TRUE)

# Create dataframe for plotting
pca_df <- as.data.frame(pca_result$x) %>%
  dplyr::select(PC1, PC2, PC3) %>%
  dplyr::bind_cols(
    ID = complete_samples$ID,
    Genotype = complete_samples$Genotype_clean,
    LPS = complete_samples$LPS_treatment_hours
  )

# Calculate variance explained
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Create PCA plot
p_pca <- ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2, color = Genotype, shape = LPS)) +
  ggplot2::geom_point(size = 3, alpha = 0.8) +
  ggplot2::scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "PCA of Cytokine Profiles",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    color = "Genotype",
    shape = "LPS Treatment (hours)"
  )

# Save PCA plot
ggplot2::ggsave(
  "src/tmp_figures/pca_plot.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# Create 3D PCA plot data (save coordinates for possible 3D plotting)
pca_3d <- data.frame(
  ID = complete_samples$ID,
  Genotype = complete_samples$Genotype_clean,
  LPS = complete_samples$LPS_treatment_hours,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3]
)

# Save PCA results for future use
saveRDS(pca_result, "src/tmp_data/pca_result.rds")
saveRDS(pca_df, "src/tmp_data/pca_plot_data.rds")
saveRDS(pca_3d, "src/tmp_data/pca_3d_data.rds")

# -----------------------------------------------
# 4. Volcano plots highlighting significantly different cytokines
# -----------------------------------------------

# Extract results to create volcano plots
volcano_data <- anova_summary %>%
  dplyr::mutate(
    Significance = case_when(
      Interaction_significant_adj ~ "Interaction",
      Treatment_significant_adj ~ "Treatment only",
      TRUE ~ "Not significant"
    ),
    NegLogPval = -log10(Treatment_p_adj)
  )

# Create volcano-like plot using -log10(p-value) vs cytokine
p_volcano <- ggplot2::ggplot(volcano_data, 
       aes(x = Cytokine, y = NegLogPval, fill = Significance)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c("Interaction" = "#7b4d8e", 
                              "Treatment only" = "#4D6D8E", 
                              "Not significant" = "#807E7D")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggplot2::labs(
    title = "Statistical Significance of Treatment Effect Across Cytokines",
    x = "Cytokine",
    y = "-log10(adjusted p-value)",
    fill = "Significance Type"
  ) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#000000")

# Save volcano plot
ggplot2::ggsave(
  "src/tmp_figures/volcano_plot.png",
  p_volcano,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------------------------
# 5. Time-course plots showing cytokine changes across LPS treatment durations
# -----------------------------------------------

# Calculate mean and SE for each cytokine by genotype and treatment
time_course_data <- cytokine_log2_long %>%
  dplyr::group_by(Genotype_clean, LPS_treatment_hours, Cytokine) %>%
  dplyr::summarize(
    Mean = mean(Log2_Value, na.rm = TRUE),
    SE = sd(Log2_Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(Genotype_clean)) # Remove NA rows

# Create time-course plots for each cytokine
for (current_cytokine in cytokine_cols) {
  cytokine_data <- time_course_data %>%
    dplyr::filter(Cytokine == current_cytokine)
  
  p <- ggplot2::ggplot(cytokine_data, 
         aes(x = LPS_treatment_hours, y = Mean, color = Genotype_clean, group = Genotype_clean)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    ggplot2::scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Time Course of", current_cytokine, "Expression"),
      x = "LPS Treatment (hours)",
      y = "Log2(Expression + 0.01)",
      color = "Genotype"
    )
  
  # Save time-course plot
  ggplot2::ggsave(
    paste0("src/tmp_figures/time_course_", current_cytokine, ".png"),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# Create a single figure with time courses for significant cytokines (if any)
if (length(sig_cytokines) > 0) {
  sig_time_data <- time_course_data %>%
    dplyr::filter(Cytokine %in% sig_cytokines)
  
  p_sig_time <- ggplot2::ggplot(sig_time_data, 
           aes(x = LPS_treatment_hours, y = Mean, color = Genotype_clean, group = Genotype_clean)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    ggplot2::facet_wrap(~ Cytokine, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Time Course of Significantly Different Cytokines",
      x = "LPS Treatment (hours)",
      y = "Log2(Expression + 0.01)",
      color = "Genotype"
    )
  
  ggplot2::ggsave(
    "src/tmp_figures/significant_cytokines_time_course.png",
    p_sig_time,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# -----------------------------------------------
# Create combined visualizations using patchwork
# -----------------------------------------------

# Create a summary boxplot of all cytokines by treatment
p_summary_box <- ggplot2::ggplot(cytokine_log2_long %>% dplyr::filter(!is.na(Genotype_clean)), 
       aes(x = LPS_treatment_hours, y = Log2_Value, fill = Genotype_clean)) +
  ggplot2::geom_boxplot(alpha = 0.7) +
  ggplot2::facet_wrap(~ Cytokine, scales = "free_y", ncol = 4) +
  ggplot2::scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background = element_rect(fill = "#807E7D"),
    strip.text = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::labs(
    title = "Overview of All Cytokines by Treatment and Genotype",
    x = "LPS Treatment (hours)",
    y = "Log2(Cytokine Level + 0.01)",
    fill = "Genotype"
  )

# Combine PCA and summary boxplot for an overview
combined_plot <- p_pca / p_summary_box + 
  patchwork::plot_layout(heights = c(1, 2))

ggplot2::ggsave(
  "src/tmp_figures/combined_overview.png",
  combined_plot,
  width = 14,
  height = 14,
  dpi = 300
)

# -----------------------------------------------
# Save visualization data for future stages
# -----------------------------------------------

# Save processed data for visualization
visualization_data <- list(
  time_course_data = time_course_data,
  pca_results = pca_result,
  pca_plot_data = pca_df,
  volcano_data = volcano_data
)

saveRDS(visualization_data, "src/tmp_data/visualization_data.rds")

# Create a structure description for the visualization data
visualization_data_structure <- list(
  description = "Processed data for cytokine visualizations",
  time_course_data = "Mean and standard error values for cytokine expression over time, grouped by genotype and treatment",
  pca_results = "Principal component analysis results for the cytokine data",
  pca_plot_data = "Processed PCA results for plotting, including PC coordinates and sample metadata",
  volcano_data = "Data for creating volcano plots showing statistical significance of cytokine differences"
)

saveRDS(visualization_data_structure, "src/tmp_data/visualization_data_structure.rds")

# Create a summary of visualizations created
visualization_summary <- list(
  boxplots = paste0("boxplot_", lps_levels[lps_levels != "none"], "_hours.png"),
  significant_boxplots = "significant_cytokines_boxplot.png",
  heatmaps = c("cytokine_heatmap.png", "significant_cytokines_heatmap.png"),
  pca_plots = "pca_plot.png",
  volcano_plots = "volcano_plot.png",
  time_course_plots = paste0("time_course_", cytokine_cols, ".png"),
  significant_time_course = "significant_cytokines_time_course.png",
  combined_plots = "combined_overview.png"
)

saveRDS(visualization_summary, "src/tmp_data/visualization_summary.rds")
