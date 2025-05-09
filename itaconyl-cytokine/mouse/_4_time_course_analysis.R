
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(pheatmap)

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

# Load necessary data files
cytokine_data_wide <- load_data_safely("src/data/cytokine_data_wide.rds", 
                                       "Error loading wide-format cytokine data")
cytokine_data_long <- load_data_safely("src/data/cytokine_data_long.rds", 
                                       "Error loading long-format cytokine data")
cytokine_names <- load_data_safely("src/data/cytokine_names.rds", 
                                   "Error loading cytokine names")
diff_analysis_top_tables <- load_data_safely("src/data/diff_analysis_top_tables.rds", 
                                             "Error loading differential analysis results")

# Check if data was loaded successfully
if (is.null(cytokine_data_wide) || is.null(cytokine_data_long) || is.null(cytokine_names)) {
  stop("Failed to load required data files")
}

# Print data dimensions
cat("Data dimensions:\n")
cat("- Wide data:", dim(cytokine_data_wide)[1], "rows,", dim(cytokine_data_wide)[2], "columns\n")
cat("- Number of cytokines:", length(cytokine_names), "\n")
cat("- Time points:", paste(levels(cytokine_data_wide$time), collapse=", "), "\n\n")

# 1. Calculate mean cytokine levels over time for each genotype
# -----------------------------------------------------------------
time_course_data <- cytokine_data_long %>%
  dplyr::group_by(genotype, time, cytokine) %>%
  dplyr::summarise(
    mean_concentration = mean(concentration, na.rm = TRUE),
    sd_concentration = sd(concentration, na.rm = TRUE),
    n = dplyr::n(),
    se_concentration = sd_concentration / sqrt(n),
    .groups = "drop"
  )

# Save time course summary data
saveRDS(time_course_data, "src/tmp_data/time_course_summary.rds")

# Display the first 10 rows of time course data
cat("\nFirst 10 rows of time course summary data:\n")
print(knitr::kable(head(time_course_data, 10)))

# 2. Calculate fold changes from baseline for each genotype
# -----------------------------------------------------------------
# Get baseline (0h) mean values for each genotype and cytokine
baseline_values <- time_course_data %>%
  dplyr::filter(time == "0 h") %>%
  dplyr::select(genotype, cytokine, mean_concentration) %>%
  dplyr::rename(baseline_mean = mean_concentration)

# Calculate fold changes from baseline for each genotype
genotype_fold_changes <- time_course_data %>%
  dplyr::left_join(baseline_values, by = c("genotype", "cytokine")) %>%
  dplyr::mutate(
    fold_change = mean_concentration / (baseline_mean + 0.01),  # Add small constant to avoid division by zero
    log2_fold_change = log2((mean_concentration + 0.01) / (baseline_mean + 0.01))
  )

# Save genotype-specific fold changes
saveRDS(genotype_fold_changes, "src/tmp_data/genotype_fold_changes.rds")

# Display the first 10 rows of fold change data
cat("\nFirst 10 rows of genotype fold changes:\n")
print(knitr::kable(head(genotype_fold_changes, 10)))

# 3. Identify cytokines with different baseline levels (time 0h)
# -----------------------------------------------------------------
# If differential analysis results are available, use them
if (!is.null(diff_analysis_top_tables) && "0 h" %in% names(diff_analysis_top_tables)) {
  baseline_diff <- diff_analysis_top_tables[["0 h"]] %>%
    dplyr::arrange(P.Value) %>%
    dplyr::mutate(
      significant = adj.P.Val < 0.1,  # Using relaxed threshold for baseline differences
      difference_type = ifelse(logFC > 0, "Higher in KO", "Higher in WT")
    )
} else {
  # If differential analysis results not available, calculate manually
  baseline_diff <- time_course_data %>%
    dplyr::filter(time == "0 h") %>%
    dplyr::select(genotype, cytokine, mean_concentration) %>%
    tidyr::pivot_wider(
      names_from = genotype,
      values_from = mean_concentration
    ) %>%
    dplyr::mutate(
      logFC = log2((ko + 0.01) / (wt + 0.01)),
      fold_change = ko / (wt + 0.01),
      difference_type = ifelse(logFC > 0, "Higher in KO", "Higher in WT")
    )
}

# Save baseline differences
saveRDS(baseline_diff, "src/tmp_data/baseline_differences.rds")

# Display baseline differences
cat("\nBaseline differences between genotypes:\n")
if ("cytokine" %in% colnames(baseline_diff) && "logFC" %in% colnames(baseline_diff)) {
  baseline_diff_display <- baseline_diff %>%
    dplyr::select(cytokine, logFC, P.Value, adj.P.Val, difference_type) %>%
    dplyr::arrange(desc(abs(logFC)))
  
  print(knitr::kable(head(baseline_diff_display, 10)))
}

# 4. Calculate response magnitudes (log2 fold change from baseline to stimulated state)
# ------------------------------------------------------------------------------------
response_magnitudes <- genotype_fold_changes %>%
  dplyr::filter(time != "0 h") %>%
  dplyr::select(genotype, time, cytokine, log2_fold_change) %>%
  tidyr::pivot_wider(
    names_from = genotype,
    values_from = log2_fold_change
  ) %>%
  dplyr::mutate(
    response_difference = ko - wt,
    difference_type = dplyr::case_when(
      abs(response_difference) < 0.5 ~ "Similar Response",
      response_difference > 0 ~ "Stronger in KO",
      TRUE ~ "Stronger in WT"
    )
  )

# Save response magnitudes
saveRDS(response_magnitudes, "src/tmp_data/response_magnitudes.rds")

# Display first 10 rows of response magnitudes
cat("\nFirst 10 rows of response magnitude differences between genotypes:\n")
print(knitr::kable(head(response_magnitudes, 10)))

# 5. Create time course line plots for each cytokine
# --------------------------------------------------
# Function to create time course plot for a single cytokine
create_time_course_plot <- function(data, cytokine_name) {
  cytokine_data <- data %>%
    dplyr::filter(cytokine == cytokine_name)
  
  ggplot2::ggplot(cytokine_data, 
                 ggplot2::aes(x = time, y = mean_concentration, 
                              color = genotype, group = genotype)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_concentration - se_concentration, 
                                         ymax = mean_concentration + se_concentration), 
                           width = 0.2) +
    ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = cytokine_name,
      x = "Time after LPS challenge",
      y = "Concentration (pg/ml)"
    ) +
    ggplot2::theme(legend.position = "bottom")
}

# Create individual time course plots for each cytokine
time_course_plots <- list()
for (cytokine_name in cytokine_names) {
  time_course_plots[[cytokine_name]] <- create_time_course_plot(time_course_data, cytokine_name)
  
  # Save individual plots - sanitize filename
  safe_name <- gsub("[/\\]", "_", cytokine_name)
  ggplot2::ggsave(
    filename = paste0("src/tmp_figures/time_course_", safe_name, ".png"),
    plot = time_course_plots[[cytokine_name]],
    width = 6,
    height = 5
  )
}

# Create multiple panel plots (4 cytokines per page)
cytokine_groups <- split(cytokine_names, ceiling(seq_along(cytokine_names)/4))

for (i in seq_along(cytokine_groups)) {
  # Create a blank plot
  multi_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  
  # Add up to 4 plots in a 2x2 grid
  for (j in seq_along(cytokine_groups[[i]])) {
    cytokine_name <- cytokine_groups[[i]][j]
    plot_obj <- time_course_plots[[cytokine_name]]
    plot_grob <- ggplot2::ggplotGrob(plot_obj)
    
    # Determine position in grid
    if (j == 1) {
      multi_plot <- multi_plot + 
        ggplot2::annotation_custom(plot_grob, xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1)
    } else if (j == 2) {
      multi_plot <- multi_plot + 
        ggplot2::annotation_custom(plot_grob, xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1)
    } else if (j == 3) {
      multi_plot <- multi_plot + 
        ggplot2::annotation_custom(plot_grob, xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.5)
    } else if (j == 4) {
      multi_plot <- multi_plot + 
        ggplot2::annotation_custom(plot_grob, xmin = 0.5, xmax = 1, ymin = 0, ymax = 0.5)
    }
  }
  
  # Save multi-panel plot
  ggplot2::ggsave(
    filename = paste0("src/tmp_figures/time_course_multi_panel_", i, ".png"),
    plot = multi_plot,
    width = 12,
    height = 10
  )
}

# 6. Create fold change heatmap to visualize temporal patterns
# -----------------------------------------------------------
# Prepare data for fold change heatmap
fold_change_matrix <- genotype_fold_changes %>%
  dplyr::filter(time != "0 h") %>%
  dplyr::select(genotype, time, cytokine, log2_fold_change) %>%
  tidyr::unite("genotype_time", genotype, time, sep = "_") %>%
  tidyr::pivot_wider(
    names_from = genotype_time,
    values_from = log2_fold_change
  )

# Create matrix for heatmap
matrix_data <- as.matrix(fold_change_matrix[, -1])
rownames(matrix_data) <- fold_change_matrix$cytokine

# Set colors for heatmap
heatmap_colors <- colorRampPalette(c("#619ca6", "#FFFFFF", "#a6617a"))(100)

# Create fold change heatmap
pheatmap::pheatmap(
  matrix_data,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize = 10,
  cellwidth = 50,
  cellheight = 12,
  main = "Log2 Fold Changes from Baseline",
  filename = "src/tmp_figures/fold_change_heatmap.png",
  width = 10,
  height = 12
)

# 7. Create response pattern classification based on temporal changes
# ------------------------------------------------------------------
# Create a time course profile for each cytokine and genotype
time_course_profiles <- time_course_data %>%
  dplyr::select(genotype, time, cytokine, mean_concentration) %>%
  tidyr::pivot_wider(
    names_from = time,
    values_from = mean_concentration
  ) %>%
  dplyr::mutate(
    baseline = `0 h`,
    early_response = `2 h` / (`0 h` + 0.01),
    late_response = `7 h` / (`0 h` + 0.01),
    resolution = `7 h` / (`2 h` + 0.01)
  )

# Classify response patterns
response_patterns <- time_course_profiles %>%
  dplyr::group_by(cytokine) %>%
  dplyr::summarise(
    early_response_ratio = early_response[genotype == "ko"] / (early_response[genotype == "wt"] + 0.01),
    late_response_ratio = late_response[genotype == "ko"] / (late_response[genotype == "wt"] + 0.01),
    resolution_ratio = resolution[genotype == "ko"] / (resolution[genotype == "wt"] + 0.01),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    early_response_pattern = dplyr::case_when(
      early_response_ratio > 1.5 ~ "Enhanced in KO",
      early_response_ratio < 0.67 ~ "Reduced in KO",
      TRUE ~ "Similar"
    ),
    late_response_pattern = dplyr::case_when(
      late_response_ratio > 1.5 ~ "Enhanced in KO",
      late_response_ratio < 0.67 ~ "Reduced in KO",
      TRUE ~ "Similar"
    ),
    resolution_pattern = dplyr::case_when(
      resolution_ratio > 1.5 ~ "Slower in KO",
      resolution_ratio < 0.67 ~ "Faster in KO",
      TRUE ~ "Similar"
    ),
    overall_pattern = dplyr::case_when(
      early_response_pattern == "Enhanced in KO" & late_response_pattern == "Enhanced in KO" ~ "Consistently Enhanced in KO",
      early_response_pattern == "Reduced in KO" & late_response_pattern == "Reduced in KO" ~ "Consistently Reduced in KO",
      early_response_pattern == "Enhanced in KO" & late_response_pattern != "Enhanced in KO" ~ "Early Enhanced in KO",
      early_response_pattern != "Enhanced in KO" & late_response_pattern == "Enhanced in KO" ~ "Late Enhanced in KO",
      early_response_pattern == "Reduced in KO" & late_response_pattern != "Reduced in KO" ~ "Early Reduced in KO",
      early_response_pattern != "Reduced in KO" & late_response_pattern == "Reduced in KO" ~ "Late Reduced in KO",
      TRUE ~ "Similar Response"
    )
  )

# Save response patterns
saveRDS(response_patterns, "src/tmp_data/response_patterns.rds")

# Display response pattern summary
cat("\nSummary of cytokine response patterns:\n")
pattern_summary <- response_patterns %>%
  dplyr::group_by(overall_pattern) %>%
  dplyr::summarise(count = n(), .groups = "drop")
print(knitr::kable(pattern_summary))

# 8. Create a heatmap with hierarchical clustering of temporal patterns
# --------------------------------------------------------------------
# Add pattern annotation to heatmap
pattern_annotation <- data.frame(
  Pattern = response_patterns$overall_pattern,
  row.names = response_patterns$cytokine
)

# Set annotation colors
ann_colors <- list(
  Pattern = c(
    "Consistently Enhanced in KO" = "#7AA661",
    "Consistently Reduced in KO" = "#4D6D8E",
    "Early Enhanced in KO" = "#8e6e4d",
    "Late Enhanced in KO" = "#7b4d8e",
    "Early Reduced in KO" = "#619ca6",
    "Late Reduced in KO" = "#a6617a",
    "Similar Response" = "#807E7D"
  )
)

# Create heatmap with annotations
pheatmap::pheatmap(
  matrix_data,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = pattern_annotation,
  annotation_colors = ann_colors,
  fontsize = 10,
  cellwidth = 50,
  cellheight = 12,
  main = "Log2 Fold Changes from Baseline (Clustered by Pattern)",
  filename = "src/tmp_figures/fold_change_heatmap_annotated.png",
  width = 10,
  height = 12
)

# 9. Create log2 fold change comparison plot
# ------------------------------------------
fold_change_comparison <- genotype_fold_changes %>%
  dplyr::filter(time != "0 h") %>%
  dplyr::select(genotype, time, cytokine, log2_fold_change) %>%
  tidyr::pivot_wider(
    names_from = genotype,
    values_from = log2_fold_change
  )

# Add significance indicators if available
if (!is.null(diff_analysis_top_tables)) {
  # Create a safer way to determine significance
  fold_change_comparison$significant <- FALSE
  
  # Update for 2h if available
  if ("2 h" %in% names(diff_analysis_top_tables)) {
    sig_cytokines_2h <- diff_analysis_top_tables[["2 h"]] %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::pull(cytokine)
    
    fold_change_comparison$significant[fold_change_comparison$time == "2 h" & 
                                      fold_change_comparison$cytokine %in% sig_cytokines_2h] <- TRUE
  }
  
  # Update for 7h if available
  if ("7 h" %in% names(diff_analysis_top_tables)) {
    sig_cytokines_7h <- diff_analysis_top_tables[["7 h"]] %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::pull(cytokine)
    
    fold_change_comparison$significant[fold_change_comparison$time == "7 h" & 
                                      fold_change_comparison$cytokine %in% sig_cytokines_7h] <- TRUE
  }
  
  # Add difference column
  fold_change_comparison <- fold_change_comparison %>%
    dplyr::mutate(difference = ko - wt)
}

# Save fold change comparison data
saveRDS(fold_change_comparison, "src/tmp_data/fold_change_comparison.rds")

# Create fold change comparison scatter plot
fold_change_plot <- ggplot2::ggplot(fold_change_comparison, 
                                   ggplot2::aes(x = wt, y = ko, color = time)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#807E7D") +
  ggplot2::geom_point(ggplot2::aes(size = abs(difference))) +
  ggplot2::geom_text(
    data = fold_change_comparison %>% dplyr::filter(significant | abs(difference) > 2),
    ggplot2::aes(label = cytokine),
    vjust = -0.5, hjust = 0.5, size = 3
  ) +
  ggplot2::scale_color_manual(values = c("2 h" = "#8e6e4d", "7 h" = "#7b4d8e")) +
  ggplot2::facet_wrap(~time) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Comparison of Log2 Fold Changes Between Genotypes",
    subtitle = "Points above the line indicate stronger response in KO mice",
    x = "Log2 Fold Change (WT)",
    y = "Log2 Fold Change (KO)"
  )

# Save fold change comparison plot
ggplot2::ggsave(
  filename = "src/tmp_figures/fold_change_comparison.png",
  plot = fold_change_plot,
  width = 10,
  height = 6
)

# 10. Identify cytokines with elevated baseline and different response magnitudes
# ----------------------------------------------------------------------------
# Identify cytokines elevated at baseline in KO
if (!is.null(baseline_diff) && "logFC" %in% colnames(baseline_diff)) {
  elevated_baseline <- baseline_diff %>%
    dplyr::filter(logFC > 0.5) %>%
    dplyr::pull(cytokine)
} else {
  elevated_baseline <- character(0)
}

# Identify cytokines with different response magnitudes
different_response <- response_magnitudes %>%
  dplyr::filter(abs(response_difference) > 1) %>%
  dplyr::pull(cytokine) %>%
  unique()

# Create a categorization of all cytokines
cytokine_categorization <- data.frame(
  cytokine = cytokine_names,
  elevated_at_baseline = cytokine_names %in% elevated_baseline,
  different_response_magnitude = cytokine_names %in% different_response
)

# Add response pattern
cytokine_categorization$response_pattern <- sapply(cytokine_names, function(c) {
  pattern <- response_patterns$overall_pattern[response_patterns$cytokine == c]
  if(length(pattern) == 0) return("Unknown")
  return(pattern)
})

# Save cytokine categorization
saveRDS(cytokine_categorization, "src/tmp_data/cytokine_categorization.rds")

# Display cytokine categorization
cat("\nCytokine categorization summary:\n")
category_summary <- cytokine_categorization %>%
  dplyr::group_by(elevated_at_baseline, different_response_magnitude) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::arrange(desc(elevated_at_baseline), desc(different_response_magnitude))

print(knitr::kable(category_summary))

# 11. Create visualization of cytokines matching the predicted patterns
# --------------------------------------------------------------------
# Create a data frame with mean values for plotting
time_values <- as.numeric(gsub(" h", "", levels(cytokine_data_wide$time)))
expanded_data <- expand.grid(
  time_numeric = time_values,
  genotype = c("wt", "ko"),
  cytokine = cytokine_names
)

# Add time factor for joining
expanded_data$time <- factor(
  paste(expanded_data$time_numeric, "h"),
  levels = levels(cytokine_data_wide$time)
)

# Join with time course data
plotting_data <- expanded_data %>%
  dplyr::left_join(
    time_course_data %>% dplyr::select(genotype, time, cytokine, mean_concentration),
    by = c("genotype", "time", "cytokine")
  )

# Group cytokines by their response patterns
pattern_groups <- list(
  "Elevated at baseline, less change" = cytokine_categorization %>% 
    dplyr::filter(elevated_at_baseline & !different_response_magnitude) %>% 
    dplyr::pull(cytokine),
  "Same at baseline, different response" = cytokine_categorization %>% 
    dplyr::filter(!elevated_at_baseline & different_response_magnitude) %>% 
    dplyr::pull(cytokine),
  "Elevated at baseline, different response" = cytokine_categorization %>% 
    dplyr::filter(elevated_at_baseline & different_response_magnitude) %>% 
    dplyr::pull(cytokine)
)

# Create plots for each pattern group
pattern_plots <- list()
for (pattern_name in names(pattern_groups)) {
  if (length(pattern_groups[[pattern_name]]) == 0) next
  
  pattern_data <- plotting_data %>%
    dplyr::filter(cytokine %in% pattern_groups[[pattern_name]]) %>%
    dplyr::mutate(cytokine = factor(cytokine, levels = pattern_groups[[pattern_name]]))
  
  pattern_plots[[pattern_name]] <- ggplot2::ggplot(
    pattern_data, 
    ggplot2::aes(x = time_numeric, y = mean_concentration, 
                 color = genotype, group = interaction(genotype, cytokine))
  ) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~cytokine, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("wt" = "#4D6D8E", "ko" = "#7AA661")) +
    ggplot2::scale_x_continuous(breaks = time_values, labels = paste(time_values, "h")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = pattern_name,
      x = "Time after LPS challenge (hours)",
      y = "Concentration (pg/ml)"
    )
  
  # Save individual pattern plots
  ggplot2::ggsave(
    filename = paste0("src/tmp_figures/pattern_", gsub(" ", "_", pattern_name), ".png"),
    plot = pattern_plots[[pattern_name]],
    width = 10,
    height = 8
  )
}

# 12. Summary of created objects
# -----------------------------
cat("\nObjects created and saved:\n")
cat("1. time_course_summary.rds - Summary of cytokine concentrations by genotype and time\n")
cat("2. genotype_fold_changes.rds - Fold changes from baseline for each genotype\n")
cat("3. baseline_differences.rds - Differences in cytokine levels at baseline\n")
cat("4. response_magnitudes.rds - Differences in response magnitudes between genotypes\n")
cat("5. response_patterns.rds - Classification of cytokine response patterns\n")
cat("6. fold_change_comparison.rds - Comparison of fold changes between genotypes\n")
cat("7. cytokine_categorization.rds - Categorization of cytokines by response patterns\n")

cat("\nFigures created and saved:\n")
cat("1. time_course_*.png - Individual time course plots for each cytokine\n")
cat("2. time_course_multi_panel_*.png - Multi-panel time course plots\n")
cat("3. fold_change_heatmap.png - Heatmap of log2 fold changes\n")
cat("4. fold_change_heatmap_annotated.png - Heatmap with response pattern annotations\n")
cat("5. fold_change_comparison.png - Scatter plot comparing fold changes between genotypes\n")
cat("6. pattern_*.png - Time course plots for cytokines grouped by response pattern\n")
