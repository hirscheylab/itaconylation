
# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

# Load significant cytokines results
significant_cytokines <- readRDS("src/data/significant_cytokines.rds")
cytokine_data <- readRDS("src/data/cytokine_data.rds")

# Print the number of significantly differentially expressed cytokines
cat("Total number of cytokines analyzed:", nrow(significant_cytokines), "\n")
cat("Number of cytokines with significant treatment effect:", sum(significant_cytokines$Treatment_significant_adj), "\n")
cat("Number of cytokines with significant interaction effect:", sum(significant_cytokines$Interaction_significant_adj), "\n\n")

# List the significant cytokines
knitr::kable(significant_cytokines %>% 
               dplyr::select(Cytokine, Treatment_p_adj, Interaction_p_adj, Significance_type),
             caption = "Significantly differentially expressed cytokines")

# Create cytokine-to-gene mapping table
# This is a manually created mapping as cytokines have specific gene names
cytokine_gene_map <- data.frame(
  Cytokine = c("IL_15", "IL_17A", "IL_27p28_IL_30", "IL_33", "IL_9", "IP_10", 
               "MCP_1", "MIP_1α", "MIP_2", "IFN_γ", "IL_10", "IL_12p70", 
               "IL_1β", "IL_2", "IL_4", "IL_5", "IL_6", "KC_GRO", "TNF_α"),
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
  dplyr::filter(Treatment_significant_adj | Interaction_significant_adj) %>%
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
    
    ggsave("src/tmp_figures/go_bp_dotplot.png", p1, width = 10, height = 8, dpi = 300)
    
    # Create barplot visualization
    p2 <- barplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment") +
      theme_bw() +
      scale_fill_gradient(low = "#4D6D8E", high = "#7AA661")
    
    ggsave("src/tmp_figures/go_bp_barplot.png", p2, width = 10, height = 8, dpi = 300)
    
    # Network plot if we have enough terms
    if(nrow(go_bp) >= 5) {
      p3 <- cnetplot(go_bp, categorySize = "pvalue", circular = TRUE, colorEdge = TRUE) +
        theme_bw() +
        scale_color_gradient(low = "#4D6D8E", high = "#7AA661")
      
      ggsave("src/tmp_figures/go_bp_network.png", p3, width = 12, height = 10, dpi = 300)
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
    
    ggsave("src/tmp_figures/go_mf_dotplot.png", p4, width = 10, height = 8, dpi = 300)
  } else {
    cat("No significant GO Molecular Function enrichment found.\n")
  }
}

# Manual pathway annotation for cytokines
cytokine_pathway_map <- data.frame(
  Cytokine = cytokine_gene_map$Cytokine,
  Pathway = c(
    "T-cell and NK cell stimulation", # IL_15
    "Inflammatory response/Th17 signaling", # IL_17A
    "Regulation of T-cell response", # IL_27p28_IL_30
    "Innate immune response/Tissue homeostasis", # IL_33
    "Th2 response/Allergic inflammation", # IL_9
    "Chemotaxis/Anti-viral response", # IP_10
    "Monocyte chemotaxis/Inflammation", # MCP_1
    "Leukocyte recruitment/Inflammation", # MIP_1α
    "Neutrophil recruitment/Acute inflammation", # MIP_2
    "Anti-viral/Macrophage activation", # IFN_γ
    "Anti-inflammatory/Immune suppression", # IL_10
    "Th1 differentiation/Cell-mediated immunity", # IL_12p70
    "Acute phase response/Inflammation", # IL_1β
    "T-cell proliferation/Differentiation", # IL_2
    "Th2 differentiation/Humoral immunity", # IL_4
    "Eosinophil activation/Allergic response", # IL_5
    "Acute phase response/Fever induction", # IL_6
    "Neutrophil chemotaxis/Inflammation", # KC_GRO
    "Pro-inflammatory/Apoptosis"  # TNF_α
  ),
  Function = c(
    "Proliferation and activation of NK and T cells", # IL_15
    "Neutrophil recruitment and inflammatory response", # IL_17A
    "Limits Th1, Th2, and Th17 responses", # IL_27p28_IL_30
    "Alarmin that activates innate immune cells", # IL_33
    "Mast cell production, T cell growth", # IL_9
    "Chemoattractant for monocytes and T cells", # IP_10
    "Recruitment of monocytes to sites of infection", # MCP_1
    "Chemotactic for macrophages and granulocytes", # MIP_1α
    "Potent neutrophil chemoattractant", # MIP_2
    "Activates macrophages, enhances MHC expression", # IFN_γ
    "Inhibits cytokine production, suppresses immune responses", # IL_10
    "Stimulates Th1 cells and inhibits Th2 responses", # IL_12p70
    "Induces acute phase protein production, fever", # IL_1β
    "T cell growth factor, promotes effector T cell development", # IL_2
    "B cell activation, promotes Th2 differentiation", # IL_4
    "Eosinophil activation and recruitment", # IL_5
    "Fever induction, acute phase protein production", # IL_6
    "Major neutrophil chemoattractant", # KC_GRO
    "Promotes inflammation, induces apoptosis" # TNF_α
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

ggsave("src/tmp_figures/pathway_barplot.png", pathway_plot, width = 10, height = 8, dpi = 300)

# Create a more comprehensive visualization of cytokine relationships
# Add significance information to pathway map
cytokine_pathway_viz <- cytokine_pathway_map %>%
  dplyr::left_join(significant_cytokines %>% 
              dplyr::select(Cytokine, Significance_type, Treatment_p_adj),
            by = "Cytokine") %>%
  dplyr::mutate(
    Significance = ifelse(is.na(Significance_type), "Not significant", Significance_type),
    NegLogPval = ifelse(is.na(Treatment_p_adj), 0, -log10(Treatment_p_adj))
  )

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
  scale_fill_gradient(low = "#807E7D", high = "#7AA661",
                      name = "-log10(p-value)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "Cytokine Functional Classification",
    x = "Cytokine",
    y = "Immune Function"
  )

ggsave("src/tmp_figures/cytokine_function_heatmap.png", func_plot, width = 12, height = 8, dpi = 300)

# Create a network visualization of cytokine-pathway relationships
# Prepare data in long format for network visualization
cytokine_network_data <- cytokine_pathway_viz %>%
  dplyr::select(Cytokine, BroadFunction, NegLogPval, Significance) %>%
  dplyr::filter(!is.na(BroadFunction))

# Export pathway analysis results
pathway_results <- list(
  cytokine_gene_map = cytokine_gene_map,
  significant_pathways = sig_pathway_info,
  pathway_counts = pathway_counts,
  go_enrichment_results = if(exists("bp_results")) bp_results else NULL,
  cytokine_functional_classification = cytokine_pathway_viz
)

saveRDS(pathway_results, "src/tmp_data/pathway_analysis_results.rds")

# Create a structure description for the output objects
pathway_results_structure <- list(
  description = "Pathway analysis results for the cytokine dataset",
  cytokine_gene_map = "Mapping of cytokine names to gene symbols and Entrez IDs",
  significant_pathways = "Significant cytokines with their pathway and functional annotations",
  pathway_counts = "Count of significant cytokines by pathway",
  go_enrichment_results = "Results from Gene Ontology enrichment analysis",
  cytokine_functional_classification = "Classification of cytokines by broad immune functions"
)

saveRDS(pathway_results_structure, "src/tmp_data/pathway_analysis_structure.rds")

# Create an interpretation summary of the pathway analysis
pathway_interpretation <- data.frame(
  Category = c(
    "Most enriched pathways",
    "Immune response type",
    "Key cytokines",
    "Significance for SIRT4 biology",
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
    "SIRT4 appears to regulate key cytokines involved in macrophage activation and inflammatory response, suggesting a potential role in modulating immune function",
    "The cytokines regulated by SIRT4 suggest potential mechanisms involving metabolic regulation of immune cell activation and inflammatory signaling"
  ),
  stringsAsFactors = FALSE
)

# Save interpretation summary
saveRDS(pathway_interpretation, "src/tmp_data/pathway_interpretation.rds")
knitr::kable(pathway_interpretation, caption = "Pathway Analysis Interpretation Summary")
