
# Load required libraries
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)
library(kableExtra)
library(SummarizedExperiment)

# Set seed for reproducibility
set.seed(123)

# Load data from previous analysis
cat("Loading data from previous analysis stages...\n")
sumexp <- readRDS("src/data/metabolomics_sumexp.rds")
sample_metadata <- as.data.frame(colData(sumexp))
abundance_matrix <- assay(sumexp, "abundance")
log_abundance_matrix <- log2(abundance_matrix + 1)
all_da_results <- readRDS("src/data/all_differential_abundance_results.rds")
pathway_analysis_results <- readRDS("src/data/pathway_analysis_results.rds")

# Check data dimensions
cat("Number of metabolites:", nrow(abundance_matrix), "\n")
cat("Number of samples:", ncol(abundance_matrix), "\n")

# Define comparisons of interest
genotype_comparisons <- c(
  "KO_vs_WT_Unstimulated", 
  "KO_vs_WT_LPS_6hr", 
  "KO_vs_WT_LPS_24hr"
)

interaction_comparisons <- c(
  "Interaction_LPS_6hr_vs_Unstim",
  "Interaction_LPS_24hr_vs_Unstim",
  "Interaction_LPS_24hr_vs_6hr"
)

# --------------------------------------------------------------
# 1. Generate volcano plots for each comparison
# --------------------------------------------------------------
cat("Creating volcano plots...\n")

# Function to create enhanced volcano plots
create_volcano_plot <- function(results_df, contrast_name, title = NULL, 
                               label_top_n = 10, fc_cutoff = 1, p_cutoff = 0.05) {
  
  # Filter for the specific contrast
  plot_data <- results_df %>%
    dplyr::filter(contrast == contrast_name)
  
  # Create title if not provided
  if (is.null(title)) {
    title <- gsub("_", " ", contrast_name)
  }
  
  # Create significance categories
  plot_data <- plot_data %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        adj.P.Val < 0.01 & abs(logFC) > fc_cutoff ~ "Highly significant",
        adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff ~ "Significant",
        adj.P.Val < p_cutoff ~ "FDR significant only",
        abs(logFC) > fc_cutoff ~ "Large fold-change only",
        TRUE ~ "Not significant"
      )
    )
  
  # Colors for different significance categories
  sig_colors <- c(
    "Highly significant" = "#4D6D8E",     # Blue
    "Significant" = "#7AA661",            # Green
    "FDR significant only" = "#807E7D",   # Gray
    "Large fold-change only" = "#8e6e4d", # Brown
    "Not significant" = "#e0e0e0"         # Light gray
  )
  
  # Find top significant metabolites to label
  top_metabolites <- plot_data %>%
    dplyr::filter(adj.P.Val < p_cutoff) %>%
    dplyr::arrange(adj.P.Val) %>%
    head(label_top_n)
  
  # Create volcano plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "#807E7D", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "#807E7D", alpha = 0.5) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::scale_color_manual(values = sig_colors) +
    ggrepel::geom_text_repel(
      data = top_metabolites,
      ggplot2::aes(label = metabolite),
      size = 3,
      fontface = "italic",
      max.overlaps = 10,
      box.padding = 0.3
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Top ", label_top_n, " significant metabolites labeled"),
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(adjusted p-value)"),
      color = "Significance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
  
  return(p)
}

# Generate volcano plots for genotype comparisons
for (comp in genotype_comparisons) {
  p <- create_volcano_plot(
    all_da_results, comp,
    title = paste("Volcano Plot:", gsub("_", " ", comp)),
    label_top_n = 15
  )
  
  # Save the plot
  filename <- paste0("src/tmp_figures/volcano_", comp, ".png")
  ggplot2::ggsave(filename, p, width = 10, height = 8)
}

# Generate volcano plots for interaction comparisons
for (comp in interaction_comparisons) {
  p <- create_volcano_plot(
    all_da_results, comp,
    title = paste("Volcano Plot:", gsub("_", " ", comp)),
    label_top_n = 15
  )
  
  # Save the plot
  filename <- paste0("src/tmp_figures/volcano_", comp, ".png")
  ggplot2::ggsave(filename, p, width = 10, height = 8)
}

# Save metadata about volcano plots
volcano_plot_info <- data.frame(
  comparison = c(genotype_comparisons, interaction_comparisons),
  plot_type = c(rep("Genotype comparison", length(genotype_comparisons)), 
                rep("Interaction effect", length(interaction_comparisons))),
  file_name = paste0("volcano_", c(genotype_comparisons, interaction_comparisons), ".png")
)
write.csv(volcano_plot_info, "src/tmp_data/volcano_plot_metadata.csv", row.names = FALSE)

# --------------------------------------------------------------
# 2. Create heatmaps of significantly different metabolites
# --------------------------------------------------------------
cat("Creating heatmaps of significantly different metabolites...\n")

# For each genotype comparison
for (comp in genotype_comparisons) {
  # Get significant metabolites for this comparison
  sig_metabolites <- all_da_results %>%
    dplyr::filter(contrast == comp, adj.P.Val < 0.05) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::pull(metabolite)
  
  # If too many significant metabolites, limit to top 50 by p-value
  if (length(sig_metabolites) > 50) {
    sig_metabolites <- all_da_results %>%
      dplyr::filter(contrast == comp, adj.P.Val < 0.05) %>%
      dplyr::arrange(adj.P.Val) %>%
      head(50) %>%
      dplyr::pull(metabolite)
  }
  
  # Create condition subset based on comparison name
  condition <- gsub("KO_vs_WT_", "", comp)
  samples_to_include <- sample_metadata$sample_id[sample_metadata$treatment == condition]
  
  # Check if we have significant metabolites
  if (length(sig_metabolites) > 0) {
    # Extract the log-abundance matrix for these metabolites and samples
    sig_abundance <- log_abundance_matrix[sig_metabolites, samples_to_include]
    
    # Create annotation for samples
    sample_annotation <- data.frame(
      Genotype = sample_metadata$genotype[match(samples_to_include, sample_metadata$sample_id)],
      row.names = samples_to_include
    )
    
    # Create a metabolite annotation with direction
    metabolite_annotation <- all_da_results %>%
      dplyr::filter(contrast == comp, metabolite %in% sig_metabolites) %>%
      dplyr::select(metabolite, logFC) %>%
      dplyr::mutate(Direction = ifelse(logFC > 0, "Up", "Down"))
    
    # Create a proper data frame for metabolite annotation
    metabolite_direction <- data.frame(
      Direction = metabolite_annotation$Direction,
      row.names = metabolite_annotation$metabolite
    )
    
    # Create annotation colors
    genotype_colors <- c("KO" = "#4D6D8E", "WT" = "#7AA661")
    direction_colors <- c("Up" = "#7AA661", "Down" = "#8e6e4d")
    
    ann_colors <- list(
      Genotype = genotype_colors,
      Direction = direction_colors
    )
    
    # Scale data for better visualization
    sig_abundance_scaled <- t(scale(t(sig_abundance)))
    
    # Create and save heatmap
    pheatmap::pheatmap(
      sig_abundance_scaled,
      annotation_col = sample_annotation,
      annotation_row = metabolite_direction,
      annotation_colors = ann_colors,
      main = paste0("Heatmap of Significant Metabolites: ", gsub("_", " ", comp)),
      color = colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      filename = paste0("src/tmp_figures/heatmap_significant_", comp, ".png"),
      width = 12,
      height = 10
    )
    
    # Save the data for this heatmap
    heatmap_data <- list(
      abundance = sig_abundance,
      scaled_abundance = sig_abundance_scaled,
      annotation_col = sample_annotation,
      annotation_row = metabolite_direction
    )
    saveRDS(heatmap_data, paste0("src/tmp_data/heatmap_data_", comp, ".rds"))
  } else {
    cat("No significant metabolites found for", comp, "- skipping heatmap.\n")
  }
}

# Create a combined heatmap with top metabolites from all comparisons
top_metabolites_all <- all_da_results %>%
  dplyr::filter(contrast %in% genotype_comparisons, adj.P.Val < 0.01) %>%
  dplyr::group_by(contrast) %>%
  dplyr::top_n(10, -adj.P.Val) %>%
  dplyr::ungroup() %>%
  dplyr::pull(metabolite) %>%
  unique()

if (length(top_metabolites_all) > 0) {
  # Extract abundance matrix for these metabolites
  combined_abundance <- log_abundance_matrix[top_metabolites_all, ]
  
  # Create annotation for samples
  all_sample_annotation <- data.frame(
    Genotype = sample_metadata$genotype,
    Treatment = sample_metadata$treatment,
    row.names = sample_metadata$sample_id
  )
  
  # Create annotation colors
  genotype_colors <- c("KO" = "#4D6D8E", "WT" = "#7AA661")
  treatment_colors <- c("Unstimulated" = "#807E7D", "LPS_6hr" = "#8e6e4d", "LPS_24hr" = "#7b4d8e")
  
  ann_colors <- list(
    Genotype = genotype_colors,
    Treatment = treatment_colors
  )
  
  # Scale data for better visualization
  combined_abundance_scaled <- t(scale(t(combined_abundance)))
  
  # Create and save combined heatmap
  pheatmap::pheatmap(
    combined_abundance_scaled,
    annotation_col = all_sample_annotation,
    annotation_colors = ann_colors,
    main = "Heatmap of Top Significant Metabolites Across All Comparisons",
    color = colorRampPalette(c("#8e6e4d", "#FFFFFF", "#4D6D8E"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    filename = "src/tmp_figures/heatmap_combined_top_metabolites.png",
    width = 12,
    height = 12
  )
  
  # Save the data for this combined heatmap
  combined_heatmap_data <- list(
    abundance = combined_abundance,
    scaled_abundance = combined_abundance_scaled,
    annotation_col = all_sample_annotation
  )
  saveRDS(combined_heatmap_data, "src/tmp_data/heatmap_data_combined.rds")
}

# --------------------------------------------------------------
# 3. Produce box plots for key metabolites across experimental groups
# --------------------------------------------------------------
cat("Creating boxplots for key metabolites...\n")

# Function to create boxplots for a given metabolite
create_metabolite_boxplot <- function(metabolite_name, data_matrix, metadata, 
                                     title = NULL, y_label = "log2(abundance + 1)") {
  # If the metabolite doesn't exist in the data, return NULL
  if (!metabolite_name %in% rownames(data_matrix)) {
    return(NULL)
  }
  
  # Extract abundance data for the metabolite
  metabolite_data <- data.frame(
    sample_id = colnames(data_matrix),
    abundance = as.numeric(data_matrix[metabolite_name, ])
  ) %>%
    dplyr::left_join(metadata, by = "sample_id")
  
  # Set the order of treatment levels
  metabolite_data$treatment <- factor(metabolite_data$treatment, 
                                    levels = c("Unstimulated", "LPS_6hr", "LPS_24hr"))
  
  # Create title if not provided
  if (is.null(title)) {
    title <- paste("Distribution of", metabolite_name)
  }
  
  # Create boxplot
  p <- ggplot2::ggplot(metabolite_data, 
                      ggplot2::aes(x = treatment, y = abundance, fill = genotype)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, height = 0, size = 1.5, alpha = 0.6, 
                        ggplot2::aes(color = genotype)) +
    ggplot2::scale_fill_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
    ggplot2::scale_color_manual(values = c("KO" = "#4D6D8E", "WT" = "#7AA661")) +
    ggplot2::labs(
      title = title,
      x = "Treatment",
      y = y_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# Identify key metabolites to plot from each comparison
key_metabolites <- list()

# From KO vs WT in unstimulated condition
key_metabolites$unstim <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_Unstimulated", adj.P.Val < 0.01) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(5) %>%
  dplyr::pull(metabolite)

# From KO vs WT in LPS 6hr condition
key_metabolites$lps6 <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_LPS_6hr", adj.P.Val < 0.01) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(5) %>%
  dplyr::pull(metabolite)

# From KO vs WT in LPS 24hr condition
key_metabolites$lps24 <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_LPS_24hr", adj.P.Val < 0.01) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(5) %>%
  dplyr::pull(metabolite)

# From interaction effects
key_metabolites$interaction <- all_da_results %>%
  dplyr::filter(grepl("Interaction", contrast), adj.P.Val < 0.01) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(5) %>%
  dplyr::pull(metabolite)

# Combine all key metabolites and remove duplicates
all_key_metabolites <- unique(c(
  key_metabolites$unstim, 
  key_metabolites$lps6, 
  key_metabolites$lps24, 
  key_metabolites$interaction
))

# Create and save individual boxplots
boxplot_list <- list()
for (metab in all_key_metabolites) {
  plot_title <- paste("Distribution of", metab, "across experimental conditions")
  
  p <- create_metabolite_boxplot(
    metabolite_name = metab,
    data_matrix = log_abundance_matrix,
    metadata = sample_metadata,
    title = plot_title
  )
  
  if (!is.null(p)) {
    boxplot_list[[metab]] <- p
    ggplot2::ggsave(
      paste0("src/tmp_figures/boxplot_", gsub("[^[:alnum:]]", "_", metab), ".png"),
      p,
      width = 8, 
      height = 6
    )
  }
}

# Create composite panels of boxplots for different categories
if (length(key_metabolites$unstim) >= 4) {
  # Unstimulated condition key metabolites
  unstim_plots <- lapply(key_metabolites$unstim[1:4], function(metab) {
    create_metabolite_boxplot(metab, log_abundance_matrix, sample_metadata)
  })
  
  # Check for NULL plots and filter them out
  unstim_plots <- unstim_plots[!sapply(unstim_plots, is.null)]
  
  if (length(unstim_plots) > 0) {
    unstim_panel <- patchwork::wrap_plots(unstim_plots, ncol = 2) +
      patchwork::plot_annotation(
        title = "Key Metabolites: KO vs WT in Unstimulated Condition",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )
    
    ggplot2::ggsave(
      "src/tmp_figures/boxplot_panel_unstimulated.png",
      unstim_panel,
      width = 12, 
      height = 10
    )
  }
}

if (length(key_metabolites$lps24) >= 4) {
  # LPS 24hr condition key metabolites
  lps24_plots <- lapply(key_metabolites$lps24[1:4], function(metab) {
    create_metabolite_boxplot(metab, log_abundance_matrix, sample_metadata)
  })
  
  # Check for NULL plots and filter them out
  lps24_plots <- lps24_plots[!sapply(lps24_plots, is.null)]
  
  if (length(lps24_plots) > 0) {
    lps24_panel <- patchwork::wrap_plots(lps24_plots, ncol = 2) +
      patchwork::plot_annotation(
        title = "Key Metabolites: KO vs WT in LPS 24hr Condition",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )
    
    ggplot2::ggsave(
      "src/tmp_figures/boxplot_panel_lps24hr.png",
      lps24_panel,
      width = 12, 
      height = 10
    )
  }
}

# Save metadata about selected metabolites
box_plot_info <- data.frame(
  metabolite = all_key_metabolites,
  file_name = paste0("boxplot_", gsub("[^[:alnum:]]", "_", all_key_metabolites), ".png"),
  selection_category = ifelse(all_key_metabolites %in% key_metabolites$unstim, "Unstimulated", 
                             ifelse(all_key_metabolites %in% key_metabolites$lps6, "LPS_6hr",
                                   ifelse(all_key_metabolites %in% key_metabolites$lps24, "LPS_24hr", 
                                         "Interaction")))
)
write.csv(box_plot_info, "src/tmp_data/boxplot_metadata.csv", row.names = FALSE)

# --------------------------------------------------------------
# 4. Visualize pathway analysis results
# --------------------------------------------------------------
cat("Visualizing pathway analysis results...\n")

# Function to create pathway enrichment plots
create_pathway_enrichment_plot <- function(enrichment_df, title, show_top_n = 10) {
  # Handle empty or missing data
  if(is.null(enrichment_df) || nrow(enrichment_df) == 0 || 
     !all(c("pathway", "hits", "fold_enrichment", "pvalue") %in% colnames(enrichment_df))) {
    return(NULL)
  }
  
  # Filter to pathways with p < 0.1 for better visualization
  plot_data <- enrichment_df %>%
    dplyr::filter(pvalue < 0.1) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pvalue),
      pathway_label = paste0(pathway, " (", hits, " hits)")
    )
  
  # If no pathways with p < 0.1, use all pathways
  if (nrow(plot_data) == 0) {
    plot_data <- enrichment_df %>%
      dplyr::mutate(
        neg_log10_p = -log10(pvalue),
        pathway_label = paste0(pathway, " (", hits, " hits)")
      )
  }
  
  # For better visualization, limit to top n pathways by p-value
  if (nrow(plot_data) > show_top_n) {
    plot_data <- plot_data %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice(1:show_top_n)
  }
  
  # Reorder pathways for better visualization
  plot_data$pathway <- factor(plot_data$pathway, 
                             levels = plot_data$pathway[order(plot_data$pvalue, decreasing = TRUE)])
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold_enrichment, y = pathway)) +
    ggplot2::geom_bar(
      stat = "identity", 
      ggplot2::aes(fill = neg_log10_p),
      alpha = 0.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0("p = ", sprintf("%.3g", pvalue))),
      hjust = -0.1,
      size = 3.5,
      color = "#807E7D"
    ) +
    ggplot2::scale_fill_gradient(
      low = "#807E7D", 
      high = "#4D6D8E",
      name = expression(-log[10](p-value))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Top ", nrow(plot_data), " enriched pathways"),
      x = "Fold Enrichment",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold")
    )
  
  return(p)
}

# Create pathway visualization for key comparisons
# KO vs WT Unstimulated
ko_unstim_enrichment <- pathway_analysis_results$enrichment_results$ko_unstim
if (!is.null(ko_unstim_enrichment) && nrow(ko_unstim_enrichment) > 0) {
  p_unstim <- create_pathway_enrichment_plot(
    ko_unstim_enrichment,
    title = "Enriched Pathways: KO vs WT (Unstimulated)"
  )
  
  if (!is.null(p_unstim)) {
    ggplot2::ggsave(
      "src/tmp_figures/pathway_enrichment_ko_unstim.png",
      p_unstim,
      width = 10, 
      height = 8
    )
  }
}

# KO vs WT LPS 6hr
ko_lps6_enrichment <- pathway_analysis_results$enrichment_results$ko_lps6
if (!is.null(ko_lps6_enrichment) && nrow(ko_lps6_enrichment) > 0) {
  p_lps6 <- create_pathway_enrichment_plot(
    ko_lps6_enrichment,
    title = "Enriched Pathways: KO vs WT (LPS 6hr)"
  )
  
  if (!is.null(p_lps6)) {
    ggplot2::ggsave(
      "src/tmp_figures/pathway_enrichment_ko_lps6.png",
      p_lps6,
      width = 10, 
      height = 8
    )
  }
}

# KO vs WT LPS 24hr
ko_lps24_enrichment <- pathway_analysis_results$enrichment_results$ko_lps24
if (!is.null(ko_lps24_enrichment) && nrow(ko_lps24_enrichment) > 0) {
  p_lps24 <- create_pathway_enrichment_plot(
    ko_lps24_enrichment,
    title = "Enriched Pathways: KO vs WT (LPS 24hr)"
  )
  
  if (!is.null(p_lps24)) {
    ggplot2::ggsave(
      "src/tmp_figures/pathway_enrichment_ko_lps24.png",
      p_lps24,
      width = 10, 
      height = 8
    )
  }
}

# Create a heatmap of pathway enrichment across conditions
if (!is.null(pathway_analysis_results$summary$combined_pathways)) {
  pathway_summary <- pathway_analysis_results$summary$combined_pathways
  
  if (nrow(pathway_summary) > 0) {
    # Extract and prepare data for heatmap
    heatmap_data <- pathway_summary %>%
      dplyr::select(pathway, pvalue_unstim, pvalue_lps6, pvalue_lps24) %>%
      dplyr::mutate(
        pvalue_unstim = -log10(pvalue_unstim),
        pvalue_lps6 = -log10(pvalue_lps6),
        pvalue_lps24 = -log10(pvalue_lps24)
      )
    
    # Convert to matrix for pheatmap
    rownames(heatmap_data) <- heatmap_data$pathway
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    
    # Limit to top pathways if there are many
    if (nrow(heatmap_matrix) > 15) {
      top_pathways <- pathway_summary %>%
        dplyr::arrange(desc(combined_score)) %>%
        head(15) %>%
        dplyr::pull(pathway)
      
      heatmap_matrix <- heatmap_matrix[top_pathways, ]
    }
    
    # Prepare column annotation
    col_anno <- data.frame(
      Condition = c("Unstimulated", "LPS 6hr", "LPS 24hr"),
      row.names = colnames(heatmap_matrix)
    )
    
    # Create the heatmap
    pheatmap::pheatmap(
      heatmap_matrix,
      annotation_col = col_anno,
      main = "Pathway Enrichment Across Conditions",
      color = colorRampPalette(c("#FFFFFF", "#4D6D8E"))(100),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      number_format = "%.1f",
      fontsize_number = 8,
      fontsize_row = 10,
      fontsize_col = 10,
      filename = "src/tmp_figures/pathway_enrichment_heatmap.png",
      width = 10,
      height = 8
    )
    
    # Save the heatmap data
    saveRDS(list(
      heatmap_data = heatmap_matrix,
      col_anno = col_anno
    ), "src/tmp_data/pathway_enrichment_heatmap_data.rds")
  }
}

# Create bubble plot for pathways
if (!is.null(pathway_analysis_results$summary$combined_pathways)) {
  pathway_summary <- pathway_analysis_results$summary$combined_pathways
  
  if (nrow(pathway_summary) > 0) {
    # Prepare data for bubble plot
    bubble_data <- pathway_summary %>%
      dplyr::select(pathway, hits_unstim, hits_lps6, hits_lps24,
                   pvalue_unstim, pvalue_lps6, pvalue_lps24) 
    
    # Convert to long format
    bubble_data_long <- bubble_data %>%
      tidyr::pivot_longer(
        cols = c(hits_unstim, hits_lps6, hits_lps24),
        names_to = "condition",
        values_to = "hits"
      ) %>%
      dplyr::mutate(
        condition = gsub("hits_", "", condition),
        pvalue = dplyr::case_when(
          condition == "unstim" ~ pvalue_unstim,
          condition == "lps6" ~ pvalue_lps6,
          condition == "lps24" ~ pvalue_lps24
        ),
        neg_log10_p = -log10(pvalue),
        # Replace NAs with 0
        hits = ifelse(is.na(hits), 0, hits),
        condition = factor(condition,
                         levels = c("unstim", "lps6", "lps24"),
                         labels = c("Unstimulated", "LPS 6hr", "LPS 24hr"))
      )
    
    # Limit to top pathways if there are many
    if (length(unique(bubble_data_long$pathway)) > 10) {
      top_pathways <- pathway_summary %>%
        dplyr::arrange(desc(combined_score)) %>%
        head(10) %>%
        dplyr::pull(pathway)
      
      bubble_data_long <- bubble_data_long %>%
        dplyr::filter(pathway %in% top_pathways)
    }
    
    # Create the bubble plot
    p_bubble <- ggplot2::ggplot(
      bubble_data_long,
      ggplot2::aes(x = condition, y = pathway, size = hits, color = neg_log10_p)
    ) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_size_continuous(range = c(1, 10), name = "Number of Hits") +
      ggplot2::scale_color_gradient(
        low = "#807E7D",
        high = "#4D6D8E",
        name = expression(-log[10](p-value))
      ) +
      ggplot2::labs(
        title = "Pathway Enrichment Across Conditions",
        x = "Treatment Condition",
        y = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 10)
      )
    
    ggplot2::ggsave(
      "src/tmp_figures/pathway_enrichment_bubble.png",
      p_bubble,
      width = 12,
      height = 8
    )
    
    # Save the bubble plot data
    saveRDS(bubble_data_long, "src/tmp_data/pathway_enrichment_bubble_data.rds")
  }
}

# --------------------------------------------------------------
# 5. Create summary tables
# --------------------------------------------------------------
cat("Creating summary tables...\n")

# Table 1: Top differentially abundant metabolites in each comparison
top_metabolites_table <- all_da_results %>%
  dplyr::filter(contrast %in% genotype_comparisons, adj.P.Val < 0.05) %>%
  dplyr::group_by(contrast) %>%
  dplyr::top_n(10, -adj.P.Val) %>%
  dplyr::arrange(contrast, adj.P.Val) %>%
  dplyr::select(
    Comparison = contrast,
    Metabolite = metabolite,
    log2FoldChange = logFC,
    FoldChange = fold_change,
    AdjustedPValue = adj.P.Val
  )

# Display the table
knitr::kable(top_metabolites_table, 
             caption = "Top Differentially Abundant Metabolites by Comparison", 
             digits = c(NA, NA, 3, 3, 6))

# Save the table
write.csv(top_metabolites_table, "src/tmp_data/top_metabolites_summary.csv", row.names = FALSE)

# Table 2: Summary of pathway enrichment
pathway_summary_table <- data.frame()

# Process each comparison if the data exists
if (!is.null(pathway_analysis_results$enrichment_results$ko_unstim) && 
    nrow(pathway_analysis_results$enrichment_results$ko_unstim) > 0) {
  
  ko_unstim_pathway <- pathway_analysis_results$enrichment_results$ko_unstim %>%
    dplyr::filter(pvalue < 0.1) %>%
    dplyr::arrange(pvalue) %>%
    head(5) %>%
    dplyr::mutate(Comparison = "KO vs WT Unstimulated") %>%
    dplyr::select(
      Comparison,
      Pathway = pathway,
      Hits = hits,
      FoldEnrichment = fold_enrichment,
      PValue = pvalue
    )
  
  pathway_summary_table <- rbind(pathway_summary_table, ko_unstim_pathway)
}

if (!is.null(pathway_analysis_results$enrichment_results$ko_lps6) && 
    nrow(pathway_analysis_results$enrichment_results$ko_lps6) > 0) {
  
  ko_lps6_pathway <- pathway_analysis_results$enrichment_results$ko_lps6 %>%
    dplyr::filter(pvalue < 0.1) %>%
    dplyr::arrange(pvalue) %>%
    head(5) %>%
    dplyr::mutate(Comparison = "KO vs WT LPS 6hr") %>%
    dplyr::select(
      Comparison,
      Pathway = pathway,
      Hits = hits,
      FoldEnrichment = fold_enrichment,
      PValue = pvalue
    )
  
  pathway_summary_table <- rbind(pathway_summary_table, ko_lps6_pathway)
}

if (!is.null(pathway_analysis_results$enrichment_results$ko_lps24) && 
    nrow(pathway_analysis_results$enrichment_results$ko_lps24) > 0) {
  
  ko_lps24_pathway <- pathway_analysis_results$enrichment_results$ko_lps24 %>%
    dplyr::filter(pvalue < 0.1) %>%
    dplyr::arrange(pvalue) %>%
    head(5) %>%
    dplyr::mutate(Comparison = "KO vs WT LPS 24hr") %>%
    dplyr::select(
      Comparison,
      Pathway = pathway,
      Hits = hits,
      FoldEnrichment = fold_enrichment,
      PValue = pvalue
    )
  
  pathway_summary_table <- rbind(pathway_summary_table, ko_lps24_pathway)
}

# Display the table
if (nrow(pathway_summary_table) > 0) {
  knitr::kable(pathway_summary_table, 
               caption = "Top Enriched Pathways by Comparison", 
               digits = c(NA, NA, 0, 3, 6))
}

# Save the table
write.csv(pathway_summary_table, "src/tmp_data/pathway_enrichment_summary.csv", row.names = FALSE)

# Table 3: Metabolites with consistent differential abundance across conditions
# Get results for all genotype comparisons and combine
unstim_results <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_Unstimulated") %>%
  dplyr::select(metabolite, logFC_unstim = logFC, padj_unstim = adj.P.Val, sig_unstim = significant)

lps6_results <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_LPS_6hr") %>%
  dplyr::select(metabolite, logFC_lps6 = logFC, padj_lps6 = adj.P.Val, sig_lps6 = significant)

lps24_results <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_LPS_24hr") %>%
  dplyr::select(metabolite, logFC_lps24 = logFC, padj_lps24 = adj.P.Val, sig_lps24 = significant)

# Join the results
consistent_metabolites <- unstim_results %>%
  dplyr::full_join(lps6_results, by = "metabolite") %>%
  dplyr::full_join(lps24_results, by = "metabolite")

# Calculate significant counts
consistent_metabolites <- consistent_metabolites %>%
  dplyr::mutate(
    significant_count = (sig_unstim | FALSE) + (sig_lps6 | FALSE) + (sig_lps24 | FALSE)
  )

# Create consistent direction flag
consistent_metabolites <- consistent_metabolites %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    consistent_direction = if(significant_count <= 1) {
      NA_character_
    } else {
      all_same_direction <- TRUE
      directions <- c()
      
      if(!is.na(sig_unstim) && sig_unstim) directions <- c(directions, sign(logFC_unstim))
      if(!is.na(sig_lps6) && sig_lps6) directions <- c(directions, sign(logFC_lps6))
      if(!is.na(sig_lps24) && sig_lps24) directions <- c(directions, sign(logFC_lps24))
      
      if(length(unique(directions)) == 1) "Yes" else "No"
    }
  ) %>%
  dplyr::ungroup()

# Filter for metabolites significant in at least 2 conditions
consistent_metabolites <- consistent_metabolites %>%
  dplyr::filter(significant_count >= 2) %>%
  dplyr::arrange(desc(significant_count), consistent_direction)

# Create a table for metabolites that are consistent across all conditions
consistent_de_table <- consistent_metabolites %>%
  dplyr::filter(significant_count == 3) %>%
  dplyr::select(
    Metabolite = metabolite,
    LFC_Unstimulated = logFC_unstim,
    LFC_LPS_6hr = logFC_lps6,
    LFC_LPS_24hr = logFC_lps24,
    ConsistentDirection = consistent_direction
  )

# Display the table
if (nrow(consistent_de_table) > 0) {
  knitr::kable(consistent_de_table, 
               caption = "Metabolites Consistently Differentially Abundant Across All Conditions", 
               digits = c(NA, 3, 3, 3, NA))
}

# Save the table
write.csv(consistent_de_table, "src/tmp_data/consistent_metabolites.csv", row.names = FALSE)

# Table 4: Treatment-specific metabolites
# Identify metabolites significant in one condition but not others
treatment_specific <- data.frame(
  Metabolite = character(),
  SpecificTo = character(),
  log2FoldChange = numeric(),
  AdjustedPValue = numeric(),
  stringsAsFactors = FALSE
)

# Get all metabolites
all_metabolites <- unique(all_da_results$metabolite)

# For each metabolite, identify if it's specific to a condition
for (metab in all_metabolites) {
  # Get results for each condition
  unstim_specifics <- all_da_results %>%
    dplyr::filter(
      metabolite == metab,
      contrast == "KO_vs_WT_Unstimulated"
    )
  
  lps6_specifics <- all_da_results %>%
    dplyr::filter(
      metabolite == metab,
      contrast == "KO_vs_WT_LPS_6hr"
    )
  
  lps24_specifics <- all_da_results %>%
    dplyr::filter(
      metabolite == metab,
      contrast == "KO_vs_WT_LPS_24hr"
    )
  
  # Check if significant in only one condition
  unstim_specific <- nrow(unstim_specifics) > 0 && unstim_specifics$adj.P.Val < 0.05 &&
    (nrow(lps6_specifics) == 0 || lps6_specifics$adj.P.Val >= 0.05) &&
    (nrow(lps24_specifics) == 0 || lps24_specifics$adj.P.Val >= 0.05)
  
  lps6_specific <- nrow(lps6_specifics) > 0 && lps6_specifics$adj.P.Val < 0.05 &&
    (nrow(unstim_specifics) == 0 || unstim_specifics$adj.P.Val >= 0.05) &&
    (nrow(lps24_specifics) == 0 || lps24_specifics$adj.P.Val >= 0.05)
  
  lps24_specific <- nrow(lps24_specifics) > 0 && lps24_specifics$adj.P.Val < 0.05 &&
    (nrow(unstim_specifics) == 0 || unstim_specifics$adj.P.Val >= 0.05) &&
    (nrow(lps6_specifics) == 0 || lps6_specifics$adj.P.Val >= 0.05)
  
  # Add to dataframe if treatment-specific
  if (unstim_specific) {
    treatment_specific <- rbind(
      treatment_specific,
      data.frame(
        Metabolite = metab,
        SpecificTo = "Unstimulated",
        log2FoldChange = unstim_specifics$logFC,
        AdjustedPValue = unstim_specifics$adj.P.Val,
        stringsAsFactors = FALSE
      )
    )
  } else if (lps6_specific) {
    treatment_specific <- rbind(
      treatment_specific,
      data.frame(
        Metabolite = metab,
        SpecificTo = "LPS_6hr",
        log2FoldChange = lps6_specifics$logFC,
        AdjustedPValue = lps6_specifics$adj.P.Val,
        stringsAsFactors = FALSE
      )
    )
  } else if (lps24_specific) {
    treatment_specific <- rbind(
      treatment_specific,
      data.frame(
        Metabolite = metab,
        SpecificTo = "LPS_24hr",
        log2FoldChange = lps24_specifics$logFC,
        AdjustedPValue = lps24_specifics$adj.P.Val,
        stringsAsFactors = FALSE
      )
    )
  }
}

# Get top 10 treatment-specific metabolites for each condition
if (nrow(treatment_specific) > 0) {
  treatment_specific <- treatment_specific %>%
    dplyr::group_by(SpecificTo) %>%
    dplyr::arrange(AdjustedPValue) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup()
  
  # Display the table
  knitr::kable(treatment_specific, 
               caption = "Top Metabolites Showing Treatment-Specific Differential Abundance", 
               digits = c(NA, NA, 3, 6))
  
  # Save the table
  write.csv(treatment_specific, "src/tmp_data/treatment_specific_metabolites.csv", row.names = FALSE)
}

# Table 5: Key biological findings
if (!is.null(pathway_analysis_results$summary$biological_findings)) {
  biological_findings <- pathway_analysis_results$summary$biological_findings
  
  # Display the table
  knitr::kable(biological_findings, 
               caption = "Summary of Key Biological Findings")
  
  # Create a more detailed summary with additional insights
  detailed_summary <- biological_findings
  
  # Add additional insights if available
  if (exists("consistent_de_table") && nrow(consistent_de_table) > 0) {
    detailed_summary <- rbind(
      detailed_summary,
      data.frame(
        finding = "Consistently differentially abundant metabolites",
        description = paste(head(consistent_de_table$Metabolite, 5), collapse = ", ")
      )
    )
  }
  
  if (exists("treatment_specific") && nrow(treatment_specific) > 0) {
    unstim_specific <- treatment_specific %>% 
      dplyr::filter(SpecificTo == "Unstimulated") %>% 
      dplyr::pull(Metabolite)
    
    lps24_specific <- treatment_specific %>% 
      dplyr::filter(SpecificTo == "LPS_24hr") %>% 
      dplyr::pull(Metabolite)
    
    if (length(unstim_specific) > 0) {
      detailed_summary <- rbind(
        detailed_summary,
        data.frame(
          finding = "Unstimulated-specific metabolites",
          description = paste(head(unstim_specific, 3), collapse = ", ")
        )
      )
    }
    
    if (length(lps24_specific) > 0) {
      detailed_summary <- rbind(
        detailed_summary,
        data.frame(
          finding = "LPS 24hr-specific metabolites",
          description = paste(head(lps24_specific, 3), collapse = ", ")
        )
      )
    }
  }
  
  # Add key pathway insights
  if (nrow(pathway_summary_table) > 0) {
    top_pathway <- pathway_summary_table %>%
      dplyr::arrange(PValue) %>%
      head(1)
    
    detailed_summary <- rbind(
      detailed_summary,
      data.frame(
        finding = paste0("Most significantly enriched pathway: ", top_pathway$Pathway),
        description = paste0(
          "Comparison: ", top_pathway$Comparison,
          ", p-value: ", sprintf("%.2e", top_pathway$PValue),
          ", fold enrichment: ", round(top_pathway$FoldEnrichment, 2)
        )
      )
    )
  }
  
  # Save detailed summary
  write.csv(detailed_summary, "src/tmp_data/detailed_biological_summary.csv", row.names = FALSE)
}

# --------------------------------------------------------------
# 6. Save all visualization metadata
# --------------------------------------------------------------
cat("Saving visualization metadata...\n")

# Save a list of all figures generated
figure_list <- list.files("src/tmp_figures", pattern = "\\.png$", full.names = FALSE)
figure_metadata <- data.frame(
  figure_name = figure_list,
  figure_type = ifelse(
    grepl("volcano", figure_list), "Volcano Plot",
    ifelse(
      grepl("heatmap", figure_list), "Heatmap",
      ifelse(
        grepl("boxplot", figure_list), "Boxplot",
        ifelse(
          grepl("pathway", figure_list), "Pathway Analysis",
          "Other"
        )
      )
    )
  )
)
write.csv(figure_metadata, "src/tmp_data/figure_metadata.csv", row.names = FALSE)

# Save a list of all data files generated
data_list <- list.files("src/tmp_data", pattern = "\\.csv$|\\.rds$", full.names = FALSE)
data_metadata <- data.frame(
  file_name = data_list,
  file_type = ifelse(
    grepl("\\.csv$", data_list), "CSV",
    ifelse(
      grepl("\\.rds$", data_list), "RDS",
      "Other"
    )
  ),
  description = ifelse(
    grepl("metadata", data_list), "Metadata",
    ifelse(
      grepl("summary", data_list), "Summary",
      ifelse(
        grepl("metabolites", data_list), "Metabolite Data",
        ifelse(
          grepl("pathway", data_list), "Pathway Data",
          "Other"
        )
      )
    )
  )
)
write.csv(data_metadata, "src/tmp_data/data_files_metadata.csv", row.names = FALSE)

# Create a complete summary of results
results_summary <- list(
  data_dimensions = list(
    metabolites = nrow(abundance_matrix),
    samples = ncol(abundance_matrix),
    significant_metabolites = sum(all_da_results$significant)
  ),
  comparisons = list(
    genotype_comparisons = genotype_comparisons,
    interaction_comparisons = interaction_comparisons
  ),
  key_metabolites = key_metabolites,
  pathway_summary = if (exists("pathway_summary_table") && nrow(pathway_summary_table) > 0) 
                       pathway_summary_table else NULL,
  consistent_metabolites = if (exists("consistent_de_table") && nrow(consistent_de_table) > 0) 
                             consistent_de_table else NULL,
  treatment_specific = if (exists("treatment_specific") && nrow(treatment_specific) > 0) 
                         treatment_specific else NULL,
  biological_findings = if (exists("detailed_summary")) detailed_summary else NULL
)

# Save the complete summary
saveRDS(results_summary, "src/tmp_data/complete_results_summary.rds")

cat("Visualization and reporting completed successfully!\n")
