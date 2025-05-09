
# Load required libraries
library(limma)
library(SummarizedExperiment)
library(dplyr)
library(stats)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(knitr)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(123)

# Load the SummarizedExperiment object and preprocessed data
cat("Loading data for differential abundance analysis...\n")
sumexp <- readRDS(here::here("data", "metabolomics_sumexp.rds"))
abundance_matrix <- assay(sumexp, "abundance")
sample_metadata <- as.data.frame(colData(sumexp))

# Log-transform the abundance values for better normality
log_abundance <- log2(abundance_matrix + 1)

# Check dimensions of log-transformed data
cat("Checking dimensions of log-transformed data...\n")
cat("Number of metabolites:", nrow(log_abundance), "\n")
cat("Number of samples:", ncol(log_abundance), "\n")

# Create factors for treatment and genotype with the correct order
sample_metadata$treatment_factor <- factor(sample_metadata$treatment,
                                           levels = c("Unstimulated", "LPS_6hr", "LPS_24hr"))
sample_metadata$genotype_factor <- factor(sample_metadata$genotype,
                                          levels = c("WT", "KO"))

# Save a summary of the experimental design
design_summary <- sample_metadata %>%
  dplyr::group_by(genotype, treatment) %>%
  dplyr::summarise(count = n(), .groups = "drop")
knitr::kable(design_summary, caption = "Experimental design summary")
write.csv(design_summary, "tmp_data/experimental_design_summary.csv", row.names = FALSE)

# -----------------------------------------------
# 1. Design the statistical model
# -----------------------------------------------
cat("Designing the statistical model...\n")

# Create a design matrix with clean column names to avoid issues with interaction terms
design <- model.matrix(~ 0 + genotype_factor:treatment_factor, data = sample_metadata)

# Rename columns to be more interpretable
colnames(design) <- c("WT_Unstimulated", "WT_LPS_6hr", "WT_LPS_24hr",
                      "KO_Unstimulated", "KO_LPS_6hr", "KO_LPS_24hr")

# Display the design matrix for verification
design_matrix_display <- as.data.frame(design)
design_matrix_display$sample_id <- sample_metadata$sample_id
knitr::kable(head(design_matrix_display), caption = "Design matrix (first 6 rows)")
write.csv(design_matrix_display, "tmp_data/design_matrix.csv", row.names = FALSE)

# -----------------------------------------------
# 2. Fit linear models for each metabolite
# -----------------------------------------------
cat("Fitting linear models for each metabolite...\n")

# Fit linear models
fit <- limma::lmFit(log_abundance, design)

# Save the initial fit object
saveRDS(fit, "tmp_data/initial_limma_fit.rds")

# -----------------------------------------------
# 3. Make contrasts for differential abundance testing
# -----------------------------------------------
cat("Creating contrasts for differential testing...\n")

# Create contrast matrix using properly named columns
contrast_matrix <- limma::makeContrasts(
  # KO vs WT in each treatment condition
  KO_vs_WT_Unstimulated = KO_Unstimulated - WT_Unstimulated,
  KO_vs_WT_LPS_6hr = KO_LPS_6hr - WT_LPS_6hr,
  KO_vs_WT_LPS_24hr = KO_LPS_24hr - WT_LPS_24hr,

  # Treatment effects within WT
  LPS_6hr_vs_Unstim_WT = WT_LPS_6hr - WT_Unstimulated,
  LPS_24hr_vs_Unstim_WT = WT_LPS_24hr - WT_Unstimulated,
  LPS_24hr_vs_LPS_6hr_WT = WT_LPS_24hr - WT_LPS_6hr,

  # Treatment effects within KO
  LPS_6hr_vs_Unstim_KO = KO_LPS_6hr - KO_Unstimulated,
  LPS_24hr_vs_Unstim_KO = KO_LPS_24hr - KO_Unstimulated,
  LPS_24hr_vs_LPS_6hr_KO = KO_LPS_24hr - KO_LPS_6hr,

  # Interaction effects (difference in treatment response between genotypes)
  Interaction_LPS_6hr_vs_Unstim = (KO_LPS_6hr - KO_Unstimulated) - (WT_LPS_6hr - WT_Unstimulated),
  Interaction_LPS_24hr_vs_Unstim = (KO_LPS_24hr - KO_Unstimulated) - (WT_LPS_24hr - WT_Unstimulated),
  Interaction_LPS_24hr_vs_6hr = (KO_LPS_24hr - KO_LPS_6hr) - (WT_LPS_24hr - WT_LPS_6hr),

  levels = design
)

# Display the contrast matrix (first 10 rows)
if (nrow(contrast_matrix) > 10) {
  knitr::kable(contrast_matrix[1:10, ], caption = "Contrast matrix (first 10 rows)")
} else {
  knitr::kable(contrast_matrix, caption = "Contrast matrix")
}
write.csv(as.data.frame(contrast_matrix), "tmp_data/contrast_matrix.csv", row.names = TRUE)

# Fit the contrasts
fit2 <- limma::contrasts.fit(fit, contrast_matrix)

# -----------------------------------------------
# 4. Empirical Bayes moderation of standard errors
# -----------------------------------------------
cat("Applying empirical Bayes moderation of standard errors...\n")

# Apply empirical Bayes moderation
fit2 <- limma::eBayes(fit2, trend = TRUE)

# Save the fitted and moderated object
saveRDS(fit2, "tmp_data/moderated_limma_fit.rds")

# -----------------------------------------------
# 5. Extract results and apply multiple testing correction
# -----------------------------------------------
cat("Extracting results with multiple testing correction...\n")

# Function to extract and format results for a specific contrast
extract_results <- function(fit, contrast, label) {
  # Extract results using topTable
  results <- limma::topTable(fit, coef = contrast, number = Inf,
                             sort.by = "p", adjust.method = "BH")

  # Add metabolite name and contrast label
  results$metabolite <- rownames(results)
  results$contrast <- label

  # Calculate fold change from log fold change
  results$fold_change <- 2^results$logFC

  # Calculate a significance flag (for convenient filtering)
  results$significant <- results$adj.P.Val < 0.05

  # Reorder columns for better readability
  results <- results %>%
    dplyr::select(metabolite, contrast, logFC, fold_change, AveExpr,
                  t, P.Value, adj.P.Val, significant, B)

  return(results)
}

# Extract results for all contrasts
all_results <- list()
contrast_names <- colnames(contrast_matrix)

for (i in 1:length(contrast_names)) {
  all_results[[i]] <- extract_results(fit2, contrast_names[i], contrast_names[i])
}

# Combine all results into a single data frame
combined_results <- do.call(rbind, all_results)

# Save all results
write.csv(combined_results, "tmp_data/all_differential_abundance_results.csv", row.names = FALSE)
saveRDS(combined_results, "tmp_data/all_differential_abundance_results.rds")

# -----------------------------------------------
# 6. Summarize significant results
# -----------------------------------------------
cat("Summarizing significant results...\n")

# Count significant metabolites for each contrast
sig_summary <- combined_results %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    total_tested = n(),
    significant_padj_05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    significant_padj_01 = sum(adj.P.Val < 0.01, na.rm = TRUE),
    upregulated = sum(adj.P.Val < 0.05 & logFC > 0, na.rm = TRUE),
    downregulated = sum(adj.P.Val < 0.05 & logFC < 0, na.rm = TRUE),
    percent_significant = round(significant_padj_05 / total_tested * 100, 2),
    .groups = "drop"
  )

# Display summary of significant results
knitr::kable(sig_summary, caption = "Summary of significant differential abundance results")
write.csv(sig_summary, "tmp_data/significant_results_summary.csv", row.names = FALSE)

# -----------------------------------------------
# 7. Extract significant results for key comparisons
# -----------------------------------------------
cat("Extracting significant results for key comparisons...\n")

# Extract significant results for KO vs WT comparisons
ko_vs_wt_results <- combined_results %>%
  dplyr::filter(grepl("KO_vs_WT", contrast)) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::arrange(contrast, adj.P.Val)

# Save KO vs WT significant results
write.csv(ko_vs_wt_results, "tmp_data/ko_vs_wt_significant_results.csv", row.names = FALSE)
saveRDS(ko_vs_wt_results, "tmp_data/ko_vs_wt_significant_results.rds")

# Extract significant results for interaction terms
interaction_results <- combined_results %>%
  dplyr::filter(grepl("Interaction", contrast)) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::arrange(contrast, adj.P.Val)

# Save interaction significant results
write.csv(interaction_results, "tmp_data/interaction_significant_results.csv", row.names = FALSE)
saveRDS(interaction_results, "tmp_data/interaction_significant_results.rds")

# Display top significant metabolites for KO vs WT in each condition
for (cond in c("Unstimulated", "LPS_6hr", "LPS_24hr")) {
  cat("\nTop significant metabolites for KO vs WT in", cond, "condition:\n")
  results_subset <- combined_results %>%
    dplyr::filter(contrast == paste0("KO_vs_WT_", cond)) %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(adj.P.Val) %>%
    head(10)

  if (nrow(results_subset) > 0) {
    print(knitr::kable(results_subset[, c("metabolite", "logFC", "fold_change", "adj.P.Val")],
                       digits = c(NA, 3, 3, 6)))
  } else {
    cat("No significant metabolites found.\n")
  }
}

# -----------------------------------------------
# 8. Visualize results
# -----------------------------------------------
cat("Creating visualizations of differential abundance results...\n")

# 1. Volcano plots for each contrast
# Function to create volcano plots
create_volcano_plot <- function(results_df, contrast_name, title_prefix = "Volcano Plot: ") {
  # Filter for the specific contrast
  plot_data <- results_df %>%
    dplyr::filter(contrast == contrast_name)

  # Create significance category for coloring
  plot_data <- plot_data %>%
    dplyr::mutate(significance = dplyr::case_when(
      adj.P.Val < 0.01 & logFC > 0.5 ~ "Highly up",
      adj.P.Val < 0.01 & logFC < -0.5 ~ "Highly down",
      adj.P.Val < 0.05 & logFC > 0.5 ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down",
      TRUE ~ "Not significant"
    ))

  # Set colors for different significance categories
  sig_colors <- c(
    "Highly up" = "red",
    "Up" = "pink",
    "Highly down" = "navy",
    "Down" = "blue",
    "Not significant" = "gray"
  )

  # Define the order you want for the legend
  legend_order <- c("Highly up", "Up", "Not significant", "Down", "Highly down")

  # Create volcano plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = -log10(P.Value), color = significance)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.8) + # Reduced point size
    ggplot2::scale_color_manual(values = sig_colors, breaks = legend_order) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#807E7D", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#807E7D", linewidth = 0.3) +
    ggplot2::labs(
      title = "",
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 7), # Smaller title
      axis.title = ggplot2::element_text(size = 6), # Smaller axis titles
      axis.text = ggplot2::element_text(size = 5),  # Smaller axis text
      legend.title = ggplot2::element_text(size = 6), # Smaller legend title
      legend.text = ggplot2::element_text(size = 5),  # Smaller legend text
      legend.key.size = ggplot2::unit(0.4, "cm"),    # Smaller legend keys
      legend.margin = ggplot2::margin(1, 1, 1, 1),   # Tighter legend margins
      legend.position = "right",
      plot.margin = ggplot2::margin(3, 3, 3, 3)      # Tighter plot margins
    )

  # Find top metabolites (reduce number for smaller plot)
  top_metabolites <- plot_data %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice(9:9) |>  # Reduced from 10 to 5 for smaller plot
    dplyr::mutate(metabolite = stringr::str_replace_all(metabolite, "/", "\n"))

  # Add text labels with ggrepel
  if(nrow(top_metabolites) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_metabolites,
      mapping = ggplot2::aes(x = logFC, y = -log10(P.Value), label = metabolite),
      color = "black",
      size = 1.8,           # Small text size
      segment.size = 0.2,   # Thin leader lines
      max.overlaps = 10,    # Allow some overlap to fit text
      min.segment.length = 0.1, # Shorter leader lines
      box.padding = 0.2,    # Less padding around text
      show.legend = FALSE
    )
  }

  return(p)
}

# Example usage for a multi-panel figure (3x3 layout)
# cond = "LPS_24hr"
# contrast_name <- paste0("KO_vs_WT_", cond)
p <- create_volcano_plot(combined_results, contrast_name, "")  # Removed title prefix to save space
ggplot2::ggsave("tmp_figures/multi_panel_figure.png", p, width = 3, height = 3, units = "in")

cond = "LPS_24hr"
contrast_name <- paste0("KO_vs_WT_", cond)
p <- create_volcano_plot(combined_results, contrast_name, "Volcano Plot: KO vs WT - ")
filename <- paste0("tmp_figures/volcano_", contrast_name, ".png")
ggplot2::ggsave(filename, p, width = 10, height = 8)

# Create volcano plots for KO vs WT comparisons
for (cond in c("Unstimulated", "LPS_6hr", "LPS_24hr")) {
  contrast_name <- paste0("KO_vs_WT_", cond)
  p <- create_volcano_plot(combined_results, contrast_name, "Volcano Plot: KO vs WT - ")

  # Save the plot
  filename <- paste0("tmp_figures/volcano_", contrast_name, ".png")
  ggplot2::ggsave(filename, p, width = 10, height = 8)
}

# Create volcano plots for interaction terms
for (inter in c("Interaction_LPS_6hr_vs_Unstim", "Interaction_LPS_24hr_vs_Unstim", "Interaction_LPS_24hr_vs_6hr")) {
  p <- create_volcano_plot(combined_results, inter, "Volcano Plot: ")

  # Save the plot
  filename <- paste0("tmp_figures/volcano_", inter, ".png")
  ggplot2::ggsave(filename, p, width = 10, height = 8)
}

# 2. Heatmap of significant metabolites across all conditions
# Get the significant metabolites from any contrast
sig_metabolites <- combined_results %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::pull(metabolite) %>%
  unique()

if (length(sig_metabolites) > 0) {
  # If there are too many significant metabolites, limit to top 50 by p-value
  if (length(sig_metabolites) > 50) {
    top_sig_metabolites <- combined_results %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::pull(metabolite) %>%
      unique() %>%
      head(50)
  } else {
    top_sig_metabolites <- sig_metabolites
  }

  # Extract the log-abundance matrix for these metabolites
  sig_abundance <- log_abundance[top_sig_metabolites, ]

  # Create annotation for samples
  sample_annotation <- data.frame(
    Genotype = sample_metadata$genotype,
    Treatment = sample_metadata$treatment,
    row.names = sample_metadata$sample_id
  )

  # Create a heatmap of significant metabolites
  pheatmap::pheatmap(
    sig_abundance,
    annotation_col = sample_annotation,
    main = "Heatmap of Significant Metabolites",
    color = colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    filename = "tmp_figures/heatmap_significant_metabolites.png",
    width = 12,
    height = 10
  )
}

# 3. Plot fold changes for top significant metabolites
# Function to create a fold change plot for top metabolites
create_fc_plot <- function(results_df, contrast_name, n_top = 15) {
  # Filter for the specific contrast and get top metabolites
  plot_data <- results_df %>%
    dplyr::filter(contrast == contrast_name) %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(adj.P.Val) %>%
    head(n_top)

  # If we have significant metabolites, create the plot
  if (nrow(plot_data) > 0) {
    # Reorder metabolites by fold change for better visualization
    plot_data$metabolite <- factor(plot_data$metabolite,
                                   levels = plot_data$metabolite[order(plot_data$logFC)])

    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = metabolite, y = logFC, fill = logFC > 0)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("#8e6e4d", "#4D6D8E"),
                                 labels = c("Down", "Up"),
                                 name = "Regulation") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0("Top Metabolites: ", contrast_name),
        x = "Metabolite",
        y = "log2 Fold Change"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10)
      )

    return(p)
  } else {
    return(NULL)
  }
}

# Create fold change plots for KO vs WT comparisons
for (cond in c("Unstimulated", "LPS_6hr", "LPS_24hr")) {
  contrast_name <- paste0("KO_vs_WT_", cond)
  p <- create_fc_plot(combined_results, contrast_name)

  if (!is.null(p)) {
    # Save the plot
    filename <- paste0("tmp_figures/foldchange_", contrast_name, ".png")
    ggplot2::ggsave(filename, p, width = 10, height = 8)
  }
}

# -----------------------------------------------
# 9. Save additional objects needed for next stages
# -----------------------------------------------
cat("Saving additional objects for next stages...\n")

# Create a differential expression results summary object
da_results_summary <- list(
  design_matrix = design,
  contrast_matrix = contrast_matrix,
  significant_summary = sig_summary,
  ko_vs_wt_results = ko_vs_wt_results,
  interaction_results = interaction_results,
  significant_metabolites = sig_metabolites
)

# Save the summary object
saveRDS(da_results_summary, "tmp_data/differential_abundance_summary.rds")

cat("Differential abundance analysis completed successfully!\n")
