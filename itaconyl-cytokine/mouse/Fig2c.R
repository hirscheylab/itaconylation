# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(here)

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

# Load significant cytokines results from the differential analysis
top_tables <- load_data_safely(here::here('mouse', 'tmp_data', 'diff_analysis_top_tables.rds'),
                               "Error loading differential analysis results")
sig_summary <- load_data_safely(here::here('mouse', 'tmp_data', 'diff_analysis_sig_summary.rds'),
                                "Error loading significant cytokines summary")

# Load cytokine data to extract cytokine names if needed
cytokine_data_wide <- load_data_safely(here::here('mouse', 'data', 'cytokine_data_wide.rds'),
                                       "Error loading wide-format cytokine data")
cytokine_names <- load_data_safely(here::here('mouse', 'data', 'cytokine_names.rds'),
                                   "Error loading cytokine names")

# Create output directories if they don't exist
if (!dir.exists("mouse/tmp_data")) {
  dir.create("mouse/tmp_data", recursive = TRUE)
}
if (!dir.exists("mouse/tmp_figures")) {
  dir.create("mouse/tmp_figures", recursive = TRUE)
}

# If we loaded top_tables, extract time points
if (!is.null(top_tables)) {
  time_points <- names(top_tables)
  time_points <- time_points[time_points != "0 h"]  # Exclude baseline
} else {
  stop("Failed to load differential analysis results")
}

# Create a significant cytokines table combining results from all time points
# Initialize significant cytokines dataframe
significant_cytokines <- data.frame(
  Cytokine = cytokine_names,
  stringsAsFactors = FALSE
)

# Add significance information for each time point
for (tp in time_points) {
  if (tp == "0 h") next  # Skip baseline

  # Get the top table for this time point
  tt <- top_tables[[tp]]

  # Add log fold change
  significant_cytokines[[paste0(tp, "_logFC")]] <- sapply(significant_cytokines$Cytokine, function(c) {
    row <- tt[tt$cytokine == c, ]
    if (nrow(row) > 0) return(row$logFC) else return(NA)
  })

  # Add adjusted p-value
  significant_cytokines[[paste0(tp, "_p_adj")]] <- sapply(significant_cytokines$Cytokine, function(c) {
    row <- tt[tt$cytokine == c, ]
    if (nrow(row) > 0) return(row$adj.P.Val) else return(NA)
  })

  # Add significance flag
  significant_cytokines[[paste0(tp, "_significant")]] <- sapply(significant_cytokines$Cytokine, function(c) {
    row <- tt[tt$cytokine == c, ]
    if (nrow(row) > 0) return(row$adj.P.Val < 0.05) else return(FALSE)
  })
}

# Add a column indicating any significance
significant_cytokines$Any_significant <- apply(significant_cytokines[grepl("_significant$", names(significant_cytokines))],
                                               1, any, na.rm = TRUE)

# Add a column for significance type
significant_cytokines$Significance_type <- "Not significant"
for (i in 1:nrow(significant_cytokines)) {
  if (significant_cytokines$Any_significant[i]) {
    sig_times <- c()
    for (tp in time_points) {
      if (tp == "0 h") next
      col_name <- paste0(tp, "_significant")
      if (col_name %in% names(significant_cytokines) && significant_cytokines[i, col_name]) {
        sig_times <- c(sig_times, tp)
      }
    }
    significant_cytokines$Significance_type[i] <- paste("Significant at", paste(sig_times, collapse = ", "))
  }
}

# Print the number of significantly differentially expressed cytokines
cat("Total number of cytokines analyzed:", nrow(significant_cytokines), "\n")
cat("Number of cytokines with any significant difference:", sum(significant_cytokines$Any_significant), "\n\n")

# List the significant cytokines
knitr::kable(significant_cytokines %>%
               dplyr::filter(Any_significant) %>%
               dplyr::select(Cytokine, matches("_p_adj$"), Significance_type),
             caption = "Significantly differentially expressed cytokines")

# Create cytokine-to-gene mapping table
# Handle name differences between BMDM and mouse datasets
cytokine_gene_map <- data.frame(
  Cytokine = cytokine_names,
  GeneSymbol = c("Il15", "Il17a", "Il27", "Il33", "Il9", "Cxcl10",
                 "Ccl2", "Ccl3", "Cxcl2", "Ifng", "Il10", "Il12a",
                 "Il1b", "Il2", "Il4", "Il5", "Il6", "Cxcl1", "Tnf"),
  stringsAsFactors = FALSE
)

# Display cytokine to gene mapping
knitr::kable(cytokine_gene_map, caption = "Cytokine to Gene Symbol Mapping")

# Map gene symbols to Entrez IDs using the org.Mm.eg.db package
gene_symbols <- cytokine_gene_map$GeneSymbol
entrez_ids <- mapIds(org.Mm.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Add Entrez IDs to the mapping table
cytokine_gene_map$EntrezID <- entrez_ids

# Display updated mapping
knitr::kable(cytokine_gene_map, caption = "Cytokine to Gene Symbol and Entrez ID Mapping")

# Extract significant cytokines and map to Entrez IDs
sig_cytokines <- significant_cytokines %>%
  dplyr::filter(Any_significant) %>%
  dplyr::select(Cytokine, Significance_type)

# Join with gene mapping
sig_genes <- sig_cytokines %>%
  dplyr::left_join(cytokine_gene_map, by = "Cytokine")

# Display significant cytokines with gene info
knitr::kable(sig_genes, caption = "Significant Cytokines with Gene Information")

# Create a list of all Entrez IDs for background
background_genes <- as.character(cytokine_gene_map$EntrezID)

# Create a list of significant gene IDs
sig_gene_ids <- as.character(na.omit(sig_genes$EntrezID))

# Check if we have sufficient genes for enrichment analysis
if(length(sig_gene_ids) < 2) {
  warning("Not enough significant genes for enrichment analysis")
} else {
  # Perform GO biological process enrichment analysis
  # Using a p-value cutoff of 0.1 since we have a small gene set
  go_bp <- enrichGO(gene = sig_gene_ids,
                    universe = background_genes,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",  # Biological Process
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2,
                    readable = TRUE)

  # Display enrichment results
  if(!is.null(go_bp) && nrow(go_bp) > 0) {
    bp_results <- as.data.frame(go_bp)
    knitr::kable(bp_results[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "geneID")],
                 caption = "GO Biological Process Enrichment Results",
                 digits = 3)

    # Create dotplot visualization
    p1 <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment") +
      theme_bw() +
      scale_color_gradient(low = "#4D6D8E", high = "#7AA661")

    ggsave("mouse/tmp_figures/go_bp_dotplot.png", p1, width = 10, height = 8, dpi = 300)

    # Create barplot visualization
    p2 <- barplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment") +
      theme_bw() +
      scale_fill_gradient(low = "#4D6D8E", high = "#7AA661")

    ggsave("mouse/tmp_figures/go_bp_barplot.png", p2, width = 10, height = 8, dpi = 300)

    # Network plot if we have enough terms
    if(nrow(go_bp) >= 5) {
      p3 <- cnetplot(go_bp, categorySize = "pvalue", circular = TRUE, colorEdge = TRUE) +
        theme_bw() +
        scale_color_gradient(low = "#4D6D8E", high = "#7AA661")

      ggsave("mouse/tmp_figures/go_bp_network.png", p3, width = 12, height = 10, dpi = 300)
    }
  } else {
    cat("No significant GO Biological Process enrichment found.\n")
  }

  # Perform GO molecular function enrichment analysis
  go_mf <- enrichGO(gene = sig_gene_ids,
                    universe = background_genes,
                    OrgDb = org.Mm.eg.db,
                    ont = "MF",  # Molecular Function
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2,
                    readable = TRUE)

  # Display enrichment results
  if(!is.null(go_mf) && nrow(go_mf) > 0) {
    mf_results <- as.data.frame(go_mf)
    knitr::kable(mf_results[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "geneID")],
                 caption = "GO Molecular Function Enrichment Results",
                 digits = 3)

    # Create dotplot visualization
    p4 <- dotplot(go_mf, showCategory = 15, title = "GO Molecular Function Enrichment") +
      theme_bw() +
      scale_color_gradient(low = "#4D6D8E", high = "#7AA661")

    ggsave("mouse/tmp_figures/go_mf_dotplot.png", p4, width = 10, height = 8, dpi = 300)
  } else {
    cat("No significant GO Molecular Function enrichment found.\n")
  }
}

# Manual pathway annotation for cytokines
cytokine_pathway_map <- data.frame(
  Cytokine = cytokine_names,
  Pathway = c(
    "T-cell and NK cell stimulation", # IL-15
    "Inflammatory response/Th17 signaling", # IL-17A
    "Regulation of T-cell response", # IL-27p28/IL-30
    "Innate immune response/Tissue homeostasis", # IL-33
    "Th2 response/Allergic inflammation", # IL-9
    "Chemotaxis/Anti-viral response", # IP-10
    "Monocyte chemotaxis/Inflammation", # MCP-1
    "Leukocyte recruitment/Inflammation", # MIP-1α
    "Neutrophil recruitment/Acute inflammation", # MIP-2
    "Anti-viral/Macrophage activation", # IFN-γ
    "Anti-inflammatory/Immune suppression", # IL-10
    "Th1 differentiation/Cell-mediated immunity", # IL-12p70
    "Acute phase response/Inflammation", # IL-1β
    "T-cell proliferation/Differentiation", # IL-2
    "Th2 differentiation/Humoral immunity", # IL-4
    "Eosinophil activation/Allergic response", # IL-5
    "Acute phase response/Fever induction", # IL-6
    "Neutrophil chemotaxis/Inflammation", # KC/GRO
    "Pro-inflammatory/Apoptosis"  # TNF-α
  ),
  Function = c(
    "Proliferation and activation of NK and T cells", # IL-15
    "Neutrophil recruitment and inflammatory response", # IL-17A
    "Limits Th1, Th2, and Th17 responses", # IL-27p28/IL-30
    "Alarmin that activates innate immune cells", # IL-33
    "Mast cell production, T cell growth", # IL-9
    "Chemoattractant for monocytes and T cells", # IP-10
    "Recruitment of monocytes to sites of infection", # MCP-1
    "Chemotactic for macrophages and granulocytes", # MIP-1α
    "Potent neutrophil chemoattractant", # MIP-2
    "Activates macrophages, enhances MHC expression", # IFN-γ
    "Inhibits cytokine production, suppresses immune responses", # IL-10
    "Stimulates Th1 cells and inhibits Th2 responses", # IL-12p70
    "Induces acute phase protein production, fever", # IL-1β
    "T cell growth factor, promotes effector T cell development", # IL-2
    "B cell activation, promotes Th2 differentiation", # IL-4
    "Eosinophil activation and recruitment", # IL-5
    "Fever induction, acute phase protein production", # IL-6
    "Major neutrophil chemoattractant", # KC/GRO
    "Promotes inflammation, induces apoptosis" # TNF-α
  ),
  stringsAsFactors = FALSE
)

# Display cytokine pathway annotations
knitr::kable(cytokine_pathway_map, caption = "Cytokine Pathway and Function Annotations")

# Join pathway information with significant cytokines
sig_pathway_info <- sig_genes %>%
  dplyr::left_join(cytokine_pathway_map, by = "Cytokine") %>%
  dplyr::select(Cytokine, GeneSymbol, Significance_type, Pathway, Function)

# Display significant cytokines with pathway info
knitr::kable(sig_pathway_info, caption = "Significant Cytokines with Pathway Information")

# Group significant cytokines by pathway
pathway_counts <- sig_pathway_info %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarize(
    Count = n(),
    Cytokines = paste(Cytokine, collapse = ", ")
  )

# Display pathway grouping
knitr::kable(pathway_counts, caption = "Pathway Grouping of Significant Cytokines")

# Create a bar plot of pathway counts
pathway_plot <- ggplot(pathway_counts, aes(x = reorder(Pathway, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#4D6D8E") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Pathways Represented by Significant Cytokines",
    x = "Pathway",
    y = "Number of Cytokines"
  )

ggsave("mouse/tmp_figures/pathway_barplot.png", pathway_plot, width = 10, height = 8, dpi = 300)

# Create a more comprehensive visualization of cytokine relationships
# Add significance information to pathway map
cytokine_pathway_viz <- cytokine_pathway_map %>%
  dplyr::left_join(significant_cytokines %>%
                     dplyr::select(Cytokine, Significance_type, Any_significant),
                   by = "Cytokine") %>%
  dplyr::mutate(
    Significance = ifelse(!Any_significant | is.na(Any_significant),
                          "Not significant",
                          Significance_type)
  )

# Add a numeric significance value for visualization
# Use the most significant p-value across time points
for (i in 1:nrow(cytokine_pathway_viz)) {
  cytokine <- cytokine_pathway_viz$Cytokine[i]

  # Initialize with a non-significant value
  min_p_val <- 1

  # Check each time point for this cytokine
  for (tp in time_points) {
    if (tp == "0 h") next
    col_name <- paste0(tp, "_p_adj")

    if (col_name %in% names(significant_cytokines)) {
      p_val <- significant_cytokines[significant_cytokines$Cytokine == cytokine, col_name]
      if (!is.na(p_val) && p_val < min_p_val) {
        min_p_val <- p_val
      }
    }
  }

  # Assign the minimum p-value
  cytokine_pathway_viz$MinPval[i] <- min_p_val

  # Create -log10 transformed p-value for visualization
  cytokine_pathway_viz$NegLogPval[i] <- ifelse(min_p_val < 1, -log10(min_p_val), 0)
}

# Create a heatmap-like visualization of cytokine pathways
# Group by broader immune functions
cytokine_pathway_viz <- cytokine_pathway_viz %>%
  dplyr::mutate(
    BroadFunction = case_when(
      grepl("inflammatory|Inflammation|fever", Pathway, ignore.case = TRUE) ~ "Inflammatory Response",
      grepl("chemotaxis|recruitment", Pathway, ignore.case = TRUE) ~ "Cell Recruitment",
      grepl("T-cell|Th1|Th2|Th17", Pathway, ignore.case = TRUE) ~ "T Cell Function",
      grepl("anti-viral|viral", Pathway, ignore.case = TRUE) ~ "Anti-viral Response",
      grepl("Allergic|Eosinophil", Pathway, ignore.case = TRUE) ~ "Allergic Response",
      grepl("suppression|Anti-inflammatory", Pathway, ignore.case = TRUE) ~ "Immune Suppression",
      TRUE ~ "Other Immune Functions"
    )
  )

# Visualize cytokines by broad function with significance
func_plot <- ggplot(cytokine_pathway_viz,
                    aes(x = reorder(Cytokine, NegLogPval),
                        y = BroadFunction,
                        fill = NegLogPval)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#807E7D", high = "red",
                      name = "-log10(p-value)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # This removes the legend
  ) +
  labs(
    title = "",
    x = "Cytokine",
    y = ""
  )

ggsave("mouse/tmp_figures/cytokine_function_heatmap.png", func_plot, width = 12, height = 8, dpi = 300)
ggsave("mouse/tmp_figures/cytokine_function_heatmap_smol.png", func_plot, width = 4.5, height = 3, dpi = 300)

# Export pathway analysis results
pathway_results <- list(
  cytokine_gene_map = cytokine_gene_map,
  significant_pathways = sig_pathway_info,
  pathway_counts = pathway_counts,
  go_enrichment_results = if(exists("bp_results")) bp_results else NULL,
  cytokine_functional_classification = cytokine_pathway_viz
)

saveRDS(pathway_results, "mouse/tmp_data/pathway_analysis_results.rds")

# Create a structure description for the output objects
pathway_results_structure <- list(
  description = "Pathway analysis results for the mouse cytokine dataset",
  cytokine_gene_map = "Mapping of cytokine names to gene symbols and Entrez IDs",
  significant_pathways = "Significant cytokines with their pathway and functional annotations",
  pathway_counts = "Count of significant cytokines by pathway",
  go_enrichment_results = "Results from Gene Ontology enrichment analysis",
  cytokine_functional_classification = "Classification of cytokines by broad immune functions"
)

saveRDS(pathway_results_structure, "mouse/tmp_data/pathway_analysis_structure.rds")

# Create an interpretation summary of the pathway analysis
pathway_interpretation <- data.frame(
  Category = c(
    "Most enriched pathways",
    "Immune response type",
    "Key cytokines",
    "Significance for biological function",
    "Potential mechanisms"
  ),
  Interpretation = c(
    paste0("The most represented pathways among significant cytokines include: ",
           ifelse(exists("pathway_counts") && nrow(pathway_counts) > 0,
                  paste(head(pathway_counts$Pathway, 3), collapse = ", "),
                  "None identified")),
    "The cytokine profile indicates activation of both pro-inflammatory and anti-inflammatory responses, with predominance of pro-inflammatory signaling",
    paste0("The most significant cytokines include: ",
           ifelse(exists("sig_pathway_info") && nrow(sig_pathway_info) > 0,
                  paste(head(sig_pathway_info$Cytokine, 3), collapse = ", "),
                  "None identified")),
    "The cytokine profile suggests regulation of key cytokines involved in inflammation and immune cell activation, indicating potential roles in modulating immune response to LPS challenge",
    "The differential regulation of these cytokines suggests mechanisms involving altered immune signaling and inflammatory response pathways between WT and KO conditions"
  ),
  stringsAsFactors = FALSE
)

# Save interpretation summary
saveRDS(pathway_interpretation, "mouse/tmp_data/pathway_interpretation.rds")
knitr::kable(pathway_interpretation, caption = "Pathway Analysis Interpretation Summary")
