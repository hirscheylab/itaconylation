# MOUSE
# Load required libraries
library(limma)
library(dplyr)
library(stats)
library(ggplot2)
library(tidyr)
library(readr)
library(ggrepel)

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
cytokine_data_wide <- load_data_safely(here::here('mouse', 'data', 'cytokine_data_wide.rds'),
                                      "Error loading wide-format cytokine data")
cytokine_data_long <- load_data_safely(here::here('mouse', 'data', 'cytokine_data_long.rds'),
                                      "Error loading long-format cytokine data")
cytokine_names <- load_data_safely(here::here('mouse', 'data', 'cytokine_names.rds'),
                                  "Error loading cytokine names")

# Check if data was loaded successfully
if (is.null(cytokine_data_wide) || is.null(cytokine_data_long) || is.null(cytokine_names)) {
  stop("Failed to load required data files")
}

# Print data dimensions
cat("Data dimensions:\n")
cat("- Wide data:", dim(cytokine_data_wide)[1], "rows,", dim(cytokine_data_wide)[2], "columns\n")
cat("- Number of cytokines:", length(cytokine_names), "\n")
cat("- Cytokines being analyzed:", paste(head(cytokine_names, 5), collapse=", "), "...\n\n")

# Log transform the data prior to differential analysis to stabilize variance
# Add small constant to handle zeros
cytokine_data_wide_log <- cytokine_data_wide %>%
  dplyr::mutate(across(all_of(cytokine_names), ~log2(.+1)))

# Convert cytokine_data_wide to a format suitable for limma
# Create design matrix for each time point separately
time_points <- levels(cytokine_data_wide$time)

# Initialize lists to store results
model_fits <- list()
contrast_matrices <- list()
ebayes_results <- list()
top_tables <- list()
volcano_data <- list()

# Function to create a volcano plot for a specific time point
create_volcano_plot <- function(data, time_point, p_threshold = 0.05, fc_threshold = 0.5) {
  # Add significance indicators
  data <- data %>%
    dplyr::mutate(
      significant = adj.P.Val < p_threshold & abs(logFC) > fc_threshold,
      label = ifelse(significant, cytokine, NA)
    )

  # Create volcano plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "#807E7D") +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "#807E7D") +
    ggplot2::scale_color_manual(values = c("FALSE" = "#807E7D", "TRUE" = "red")) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      box.padding = 0.5,
      max.overlaps = 15,
      size = 3
    ) +
    ggplot2::labs(
      #title = paste("Volcano Plot -", time_point),
      #subtitle = paste("KO vs WT Comparison at", time_point),
      x = "Log2 Fold Change (KO/WT)",
      y = "-Log10(Adjusted P-value)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none"#,
      #plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      #plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}

# Loop through each time point and perform differential expression analysis
for (time_point in time_points) {
  # Filter data for the specific time point
  tp_data <- cytokine_data_wide_log %>%
    dplyr::filter(time == time_point)

  # Skip if insufficient samples
  if (nrow(tp_data) < 4) {
    cat(paste("Skipping", time_point, "- not enough samples for analysis\n"))
    next
  }

  # Create design matrix - compare KO vs WT
  genotype_factor <- factor(tp_data$genotype, levels = c("wt", "ko"))
  design <- stats::model.matrix(~0+genotype_factor)
  colnames(design) <- c("wt", "ko")

  # Create expression matrix
  expr_matrix <- as.matrix(tp_data[, cytokine_names])
  rownames(expr_matrix) <- tp_data$sample_id

  # Fit linear model
  fit <- limma::lmFit(t(expr_matrix), design)

  # Define contrast matrix for KO vs WT comparison
  contrast_matrix <- limma::makeContrasts(ko-wt, levels = design)

  # Calculate contrasts
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)

  # Compute empirical Bayes statistics
  fit2 <- limma::eBayes(fit2, trend = TRUE)

  # Get the top table results
  top_table <- limma::topTable(fit2, number = Inf, sort.by = "P")
  top_table$cytokine <- rownames(top_table)

  # Store results
  model_fits[[time_point]] <- fit
  contrast_matrices[[time_point]] <- contrast_matrix
  ebayes_results[[time_point]] <- fit2
  top_tables[[time_point]] <- top_table

  # Format data for volcano plot
  volcano_data[[time_point]] <- top_table %>%
    dplyr::select(cytokine, logFC, P.Value, adj.P.Val)

  # Print summary of significant cytokines
  sig_count <- sum(top_table$adj.P.Val < 0.05)
  cat(paste("Time point:", time_point, "\n"))
  cat(paste("  - Significantly different cytokines (FDR < 0.05):", sig_count, "\n"))
  if (sig_count > 0) {
    sig_cytokines <- top_table %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::arrange(adj.P.Val)
    cat("  - Top significant cytokines:\n")
    for (i in 1:min(3, nrow(sig_cytokines))) {
      cat(paste("    *", sig_cytokines$cytokine[i],
                "(logFC =", round(sig_cytokines$logFC[i], 3),
                ", adj.P =", format(sig_cytokines$adj.P.Val[i], scientific = TRUE, digits = 3), ")\n"))
    }
  }
  cat("\n")
}

# Save all differential analysis results
saveRDS(top_tables, "mouse/tmp_data/diff_analysis_top_tables.rds")
saveRDS(volcano_data, "mouse/tmp_data/diff_analysis_volcano_data.rds")

# Create volcano plots for each time point
volcano_plots <- lapply(names(volcano_data), function(tp) {
  create_volcano_plot(volcano_data[[tp]], tp)
})

# Save individual volcano plots
for (i in seq_along(volcano_plots)) {
  tp <- names(volcano_data)[i]
  filename <- paste0("mouse/tmp_figures/volcano_plot_", gsub(" ", "_", tp), ".png")
  ggplot2::ggsave(filename, volcano_plots[[i]], width = 3, height = 3, dpi = 300)
}

# Create combined volcano plot for all time points (if we have multiple)
if (length(volcano_plots) > 1) {
  # Create a multi-panel plot manually (avoiding patchwork dependency)
  combined_plot <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotation_custom(
      ggplot2::ggplotGrob(volcano_plots[[1]]),
      xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1
    )

  if (length(volcano_plots) >= 2) {
    combined_plot <- combined_plot +
      ggplot2::annotation_custom(
        ggplot2::ggplotGrob(volcano_plots[[2]]),
        xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1
      )
  }

  if (length(volcano_plots) >= 3) {
    combined_plot <- combined_plot +
      ggplot2::annotation_custom(
        ggplot2::ggplotGrob(volcano_plots[[3]]),
        xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.5
      )
  }

  # Save combined plot
  ggplot2::ggsave("src/tmp_figures/volcano_plots_combined.png", combined_plot,
         width = 12, height = 10, dpi = 300)
}

# Prepare fold change data at each time point for easier downstream analysis
fc_by_time <- list()

for (time_point in time_points) {
  # Skip baseline for fold change calculation
  if (time_point == "0 h") next

  # Filter data for the time point
  time_data <- cytokine_data_wide %>%
    dplyr::filter(time == time_point) %>%
    dplyr::select(sample_id, genotype, all_of(cytokine_names))

  # Calculate average expression for each genotype
  wt_means <- time_data %>%
    dplyr::filter(genotype == "wt") %>%
    dplyr::select(all_of(cytokine_names)) %>%
    dplyr::summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

  ko_means <- time_data %>%
    dplyr::filter(genotype == "ko") %>%
    dplyr::select(all_of(cytokine_names)) %>%
    dplyr::summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

  # Calculate fold changes (KO/WT)
  fc_values <- data.frame(cytokine = cytokine_names)
  fc_values$ko_mean <- sapply(cytokine_names, function(c) ko_means[[c]])
  fc_values$wt_mean <- sapply(cytokine_names, function(c) wt_means[[c]])

  # Calculate fold changes handling zeros
  fc_values$fold_change <- sapply(1:nrow(fc_values), function(i) {
    ko_val <- fc_values$ko_mean[i]
    wt_val <- fc_values$wt_mean[i]

    # Handle zero values with a small constant to avoid division by zero
    if (wt_val == 0) {
      if (ko_val == 0) return(1)  # Both zero means no change
      return(ko_val / 0.01)  # WT is zero, KO has value
    }
    return(ko_val / wt_val)
  })

  # Calculate log2 fold change
  fc_values$log2_fold_change <- log2(fc_values$fold_change)

  # Add p-values from differential analysis
  if (time_point %in% names(top_tables)) {
    fc_values <- fc_values %>%
      dplyr::left_join(
        top_tables[[time_point]] %>% dplyr::select(cytokine, P.Value, adj.P.Val),
        by = "cytokine"
      )
  }

  # Store the fold change data
  fc_by_time[[time_point]] <- fc_values
}

# Save fold change data by time point
saveRDS(fc_by_time, "src/tmp_data/diff_analysis_fc_by_time.rds")

# Create a fold change heatmap table (combining data across time points)
if (length(fc_by_time) > 0) {
  # Combine all fold change data
  fc_combined <- do.call(rbind, lapply(names(fc_by_time), function(tp) {
    fc_by_time[[tp]] %>%
      dplyr::mutate(time = tp)
  }))

  # Create wide format for heatmap (cytokines as rows, time points as columns)
  fc_heatmap_data <- fc_combined %>%
    dplyr::select(cytokine, time, log2_fold_change) %>%
    tidyr::pivot_wider(names_from = time, values_from = log2_fold_change)

  # Save data for heatmap
  saveRDS(fc_heatmap_data, "src/tmp_data/diff_analysis_fc_heatmap_data.rds")
}

# Create a summary table of significant cytokines across time points
sig_summary <- data.frame(
  cytokine = cytokine_names,
  stringsAsFactors = FALSE
)

for (tp in names(top_tables)) {
  if (tp == "0 h") next  # Skip baseline

  tt <- top_tables[[tp]]
  sig_summary[[paste0(tp, "_logFC")]] <- sapply(sig_summary$cytokine, function(c) {
    row <- tt[tt$cytokine == c, ]
    if (nrow(row) > 0) return(row$logFC) else return(NA)
  })

  sig_summary[[paste0(tp, "_adj.P.Val")]] <- sapply(sig_summary$cytokine, function(c) {
    row <- tt[tt$cytokine == c, ]
    if (nrow(row) > 0) return(row$adj.P.Val) else return(NA)
  })

  sig_summary[[paste0(tp, "_significant")]] <- sapply(sig_summary$cytokine, function(c) {
    row <- tt[tt$cytokine == c, ]
    if (nrow(row) > 0) return(row$adj.P.Val < 0.05) else return(FALSE)
  })
}

# Save significant cytokines summary
saveRDS(sig_summary, "mouse/tmp_data/diff_analysis_sig_summary.rds")

# Filter for significant cytokines in any time point
sig_any_tp <- sig_summary %>%
  dplyr::mutate(significant_any = rowSums(dplyr::select(., dplyr::ends_with("_significant")), na.rm = TRUE) > 0) %>%
  dplyr::filter(significant_any)

# Prepare table for display
display_cols <- c("cytokine", grep("_logFC$|_adj.P.Val$", names(sig_summary), value = TRUE))
sig_display <- sig_summary %>%
  dplyr::filter(cytokine %in% sig_any_tp$cytokine) %>%
  dplyr::select(all_of(display_cols))

# Format p-values for display
for (col in grep("adj.P.Val", names(sig_display), value = TRUE)) {
  sig_display[[col]] <- format(sig_display[[col]], scientific = TRUE, digits = 3)
}

# Format fold changes for display
for (col in grep("logFC", names(sig_display), value = TRUE)) {
  sig_display[[col]] <- round(sig_display[[col]], 3)
}

# Display summary of significant cytokines
cat("Summary of significantly different cytokines (KO vs WT):\n")
print(knitr::kable(sig_display))

# Save list of objects created
cat("\nSaved objects:\n")
cat("1. diff_analysis_top_tables.rds - Results of differential expression analysis\n")
cat("2. diff_analysis_volcano_data.rds - Data for volcano plots\n")
cat("3. diff_analysis_fc_by_time.rds - Fold changes by time point\n")
cat("4. diff_analysis_fc_heatmap_data.rds - Fold change data formatted for heatmap\n")
cat("5. diff_analysis_sig_summary.rds - Summary of significant cytokines\n")

cat("\nSaved figures:\n")
cat("1. volcano_plot_*.png - Individual volcano plots for each time point\n")
cat("2. volcano_plots_combined.png - Combined volcano plots for all time points\n")
