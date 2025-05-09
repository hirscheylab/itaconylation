
# Load required libraries
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(knitr)
library(RColorBrewer)
library(stringr)

# Set seed for reproducibility
set.seed(123)

# Load the differential abundance results from previous stage
cat("Loading differential abundance results...\n")
sumexp <- readRDS("src/data/metabolomics_sumexp.rds")
sample_metadata <- as.data.frame(colData(sumexp))
differential_abundance_summary <- readRDS("src/data/differential_abundance_summary.rds")
all_da_results <- readRDS("src/data/all_differential_abundance_results.rds")

# Extract significant metabolites
sig_metabolites <- differential_abundance_summary$significant_metabolites

# Print summary of significant metabolites
cat("Total number of significant metabolites:", length(sig_metabolites), "\n")
if(length(sig_metabolites) > 0) {
  cat("First 10 significant metabolites:", paste(head(sig_metabolites, 10), collapse = ", "), "\n")
}

# Create table of significant metabolites by contrast
sig_by_contrast <- all_da_results %>%
  dplyr::filter(significant == TRUE) %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    significant_count = n(),
    up_regulated = sum(logFC > 0),
    down_regulated = sum(logFC < 0),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(significant_count))

# Print table of significant metabolites by contrast
knitr::kable(sig_by_contrast, caption = "Significant metabolites by contrast")
write.csv(sig_by_contrast, "src/tmp_data/significant_by_contrast.csv", row.names = FALSE)

# -----------------------------------------------------
# 1. Prepare metabolite data for pathway analysis
# -----------------------------------------------------
cat("Preparing metabolite data for pathway analysis...\n")

# Extract significant metabolites by comparison group
ko_unstim_metabolites <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_Unstimulated", significant == TRUE) %>%
  dplyr::arrange(adj.P.Val)

ko_lps6_metabolites <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_LPS_6hr", significant == TRUE) %>%
  dplyr::arrange(adj.P.Val)

ko_lps24_metabolites <- all_da_results %>%
  dplyr::filter(contrast == "KO_vs_WT_LPS_24hr", significant == TRUE) %>%
  dplyr::arrange(adj.P.Val)

int_lps6_unstim <- all_da_results %>%
  dplyr::filter(contrast == "Interaction_LPS_6hr_vs_Unstim", significant == TRUE) %>%
  dplyr::arrange(adj.P.Val)

int_lps24_unstim <- all_da_results %>%
  dplyr::filter(contrast == "Interaction_LPS_24hr_vs_Unstim", significant == TRUE) %>%
  dplyr::arrange(adj.P.Val)

int_lps24_lps6 <- all_da_results %>%
  dplyr::filter(contrast == "Interaction_LPS_24hr_vs_6hr", significant == TRUE) %>%
  dplyr::arrange(adj.P.Val)

# Save the lists of significant metabolites for each comparison
write.csv(ko_unstim_metabolites, "src/tmp_data/ko_vs_wt_unstim_significant.csv", row.names = FALSE)
write.csv(ko_lps6_metabolites, "src/tmp_data/ko_vs_wt_lps6hr_significant.csv", row.names = FALSE)
write.csv(ko_lps24_metabolites, "src/tmp_data/ko_vs_wt_lps24hr_significant.csv", row.names = FALSE)
write.csv(int_lps6_unstim, "src/tmp_data/interaction_lps6hr_vs_unstim_significant.csv", row.names = FALSE)
write.csv(int_lps24_unstim, "src/tmp_data/interaction_lps24hr_vs_unstim_significant.csv", row.names = FALSE)
write.csv(int_lps24_lps6, "src/tmp_data/interaction_lps24hr_vs_lps6hr_significant.csv", row.names = FALSE)

# -----------------------------------------------------
# 2. Categorize metabolites into known pathways manually
# -----------------------------------------------------
cat("Categorizing metabolites into known metabolic pathways...\n")

# Define main metabolic pathways and associated metabolites
# This is a manual curation based on biological knowledge
metabolic_pathways <- list(
  "Glycolysis" = c("glucose", "glucose-6-phosphate", "fructose-6-phosphate", 
                   "fructose-1,6-bisphosphate", "glyceraldehyde 3-phosphate", 
                   "Dihydroxyacetone phosphate", "DHAP", "pyruvate", "lactate", 
                   "1,3-bisphosphoglycerate", "3-phosphoglycerate", "2-phosphoglycerate", 
                   "phosphoenolpyruvate"),
  
  "TCA cycle" = c("citrate", "cis-aconitate", "isocitrate", "2-oxoglutarate", "alpha-ketoglutarate",
                 "succinate", "fumarate", "malate", "oxaloacetate"),
  
  "Pentose Phosphate Pathway" = c("glucose-6-phosphate", "6-phosphogluconolactone", 
                                 "6-phosphogluconate", "ribulose-5-phosphate", 
                                 "ribose-5-phosphate", "xylulose-5-phosphate", 
                                 "sedoheptulose-7-phosphate", "erythrose-4-phosphate",
                                 "Sedoheptulose 7-phosphate"),
  
  "Nucleotide Metabolism" = c("adenosine", "AMP", "ADP", "ATP", "guanosine", "GMP", "GDP", "GTP",
                             "uridine", "UMP", "UDP", "UTP", "cytidine", "CMP", "CDP", "CTP",
                             "thymidine", "TMP", "TDP", "TTP", "deoxyadenosine", "dAMP", "dADP", "dATP",
                             "deoxyguanosine", "dGMP", "dGDP", "dGTP", "deoxycytidine", "dCMP", "dCDP", "dCTP",
                             "deoxythymidine", "dTMP", "dTDP", "dTTP", "Hypoxanthine", "Xanthine", "Urate",
                             "Inosine", "IMP", "CDP", "CMP", "cytidine", "Deoxyuridine", "dGDP/ADP", 
                             "dGMP/AMP", "dUDP", "GMP", "guanosine", "Hypoxanthine", "IMP", 
                             "Thymidine", "UDP", "UMP", "uridine", "UTP", "Xanthosine", "Urate"),
  
  "Amino Acid Metabolism" = c("alanine", "arginine", "asparagine", "aspartate", "cysteine", 
                             "glutamine", "glutamate", "glycine", "histidine", "isoleucine", 
                             "leucine", "lysine", "methionine", "phenylalanine", "proline", 
                             "serine", "threonine", "tryptophan", "tyrosine", "valine",
                             "N-acetyl-aspartate", "citrulline", "ornithine", "putrescine",
                             "spermidine", "spermine", "creatine", "creatinine", "GABA",
                             "taurine", "Isoleucine/Leucine/Norleucine", "Phenylalanine"),
  
  "Lipid Metabolism" = c("palmitate", "palmitoleate", "stearate", "oleate", "linoleate", 
                        "linolenate", "arachidonate", "docosahexaenoate", "ceramide", 
                        "sphingosine", "sphinganine", "cholesterol", "cholesteryl ester",
                        "phosphatidylcholine", "phosphatidylethanolamine", "phosphatidylserine",
                        "phosphatidylinositol", "phosphatidate", "triacylglycerol", "diacylglycerol",
                        "monoacylglycerol", "fatty acid", "glycerol", "glycerol 3-phosphate",
                        "palmitoleate", "arachidonic acid", "oleic acid/elaidate/trans-vaccenate",
                        "clupanodonic acid", "eicosatrienoic acid", "DGLA"),
  
  "Glutathione Metabolism" = c("glutathione", "glutathione disulfide", "cysteinylglycine", "glutamate", 
                              "cysteine", "glycine", "glutamylcysteine", "ophthalmate", "GSSG", "GSH"),
  
  "Vitamin & Cofactor Metabolism" = c("thiamine", "thiamine phosphate", "thiamine pyrophosphate", 
                                    "riboflavin", "FMN", "FAD", "niacin", "NAD+", "NADH", "NADP+", 
                                    "NADPH", "pantothenate", "CoA", "acetyl-CoA", "biotin", 
                                    "pyridoxine", "pyridoxal", "pyridoxamine", "pyridoxal phosphate", 
                                    "folate", "tetrahydrofolate", "5-methyltetrahydrofolate", 
                                    "vitamin B12", "ascorbate", "tocopherol", "phylloquinone", 
                                    "menaquinone", "ubiquinone", "ubiquinol", "choline", 
                                    "NADH", "Pyridoxal", "thiamine", "NAD+", "NADP+"),
  
  "Carnitine Metabolism" = c("carnitine", "acetylcarnitine", "propionylcarnitine", 
                           "butyrylcarnitine", "hexanoylcarnitine", "octanoylcarnitine", 
                           "decanoylcarnitine", "lauroylcarnitine", "myristoylcarnitine", 
                           "palmitoylcarnitine", "stearoylcarnitine", "oleoylcarnitine", 
                           "linoleoylcarnitine", "(R)-carnitine", "hydroxy-isovaleryl carnitine",
                           "Hydroxybutyrylcarnitine", "hydroxyhexadecanoylcarnitine"),
  
  "Itaconate pathway" = c("Itaconate", "cis-aconitate", "citrate", "itaconyl-CoA", "itaconic acid"),
  
  "Purine metabolism" = c("adenosine", "guanosine", "inosine", "xanthosine", "hypoxanthine", 
                        "xanthine", "urate", "adenine", "guanine", "AMP", "GMP", "IMP", "XMP",
                        "Guanine", "Inosine", "Xanthosine", "ADP", "GDP", "IDP", "XDP",
                        "ATP", "GTP", "ITP", "XTP")
)

# Function to map metabolites to pathways
map_metabolites_to_pathways <- function(metabolite_list, pathway_dict) {
  if(length(metabolite_list) == 0) {
    return(list(
      pathway_counts = data.frame(pathway=character(), hits=numeric(), 
                               total_in_pathway=numeric(), coverage=numeric(), 
                               stringsAsFactors=FALSE),
      metabolite_pathway_map = data.frame(metabolite=character(), pathway=character(), 
                                      stringsAsFactors=FALSE),
      pathway_hits = list()
    ))
  }
  
  pathway_hits <- list()
  metabolite_pathway_map <- data.frame(
    metabolite = character(),
    pathway = character(),
    stringsAsFactors = FALSE
  )
  
  # For each pathway, check which metabolites match
  for (pathway_name in names(pathway_dict)) {
    pathway_metabolites <- pathway_dict[[pathway_name]]
    
    # Find matches (using partial matching to handle differences in naming)
    matches <- sapply(metabolite_list, function(metab) {
      any(sapply(pathway_metabolites, function(pm) {
        tolower(metab) == tolower(pm) || 
          grepl(tolower(pm), tolower(metab), fixed = TRUE) ||
          grepl(tolower(metab), tolower(pm), fixed = TRUE)
      }))
    })
    
    matched_metabolites <- metabolite_list[matches]
    
    if (length(matched_metabolites) > 0) {
      pathway_hits[[pathway_name]] <- matched_metabolites
      
      # Add to metabolite-pathway mapping
      for (metab in matched_metabolites) {
        metabolite_pathway_map <- rbind(
          metabolite_pathway_map,
          data.frame(
            metabolite = metab,
            pathway = pathway_name,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  # Count hits per pathway
  pathway_counts <- sapply(pathway_hits, length)
  
  if(length(pathway_counts) > 0) {
    pathway_df <- data.frame(
      pathway = names(pathway_counts),
      hits = as.numeric(pathway_counts),
      stringsAsFactors = FALSE
    )
    
    # Calculate coverage (hits / total pathway metabolites)
    pathway_df$total_in_pathway <- sapply(pathway_df$pathway, function(path) {
      length(pathway_dict[[path]])
    })
    
    pathway_df$coverage <- pathway_df$hits / pathway_df$total_in_pathway
    pathway_df <- pathway_df[order(-pathway_df$hits), ]
  } else {
    pathway_df <- data.frame(
      pathway = character(),
      hits = numeric(),
      total_in_pathway = numeric(),
      coverage = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    pathway_counts = pathway_df,
    metabolite_pathway_map = metabolite_pathway_map,
    pathway_hits = pathway_hits
  ))
}

# Apply the function to each comparison's significant metabolites
ko_unstim_pathways <- map_metabolites_to_pathways(ko_unstim_metabolites$metabolite, metabolic_pathways)
ko_lps6_pathways <- map_metabolites_to_pathways(ko_lps6_metabolites$metabolite, metabolic_pathways)
ko_lps24_pathways <- map_metabolites_to_pathways(ko_lps24_metabolites$metabolite, metabolic_pathways)
int_lps6_unstim_pathways <- map_metabolites_to_pathways(int_lps6_unstim$metabolite, metabolic_pathways)
int_lps24_unstim_pathways <- map_metabolites_to_pathways(int_lps24_unstim$metabolite, metabolic_pathways)
int_lps24_lps6_pathways <- map_metabolites_to_pathways(int_lps24_lps6$metabolite, metabolic_pathways)

# Save pathway mappings
saveRDS(ko_unstim_pathways, "src/tmp_data/ko_unstim_pathway_mapping.rds")
saveRDS(ko_lps6_pathways, "src/tmp_data/ko_lps6_pathway_mapping.rds")
saveRDS(ko_lps24_pathways, "src/tmp_data/ko_lps24_pathway_mapping.rds")
saveRDS(int_lps6_unstim_pathways, "src/tmp_data/int_lps6_unstim_pathway_mapping.rds")
saveRDS(int_lps24_unstim_pathways, "src/tmp_data/int_lps24_unstim_pathway_mapping.rds")
saveRDS(int_lps24_lps6_pathways, "src/tmp_data/int_lps24_lps6_pathway_mapping.rds")

# Display results for first few pathway mappings
if(nrow(ko_unstim_pathways$pathway_counts) > 0) {
  knitr::kable(ko_unstim_pathways$pathway_counts, 
               caption = "Pathway enrichment for KO vs WT in Unstimulated condition")
  write.csv(ko_unstim_pathways$pathway_counts, 
            "src/tmp_data/ko_unstim_pathway_counts.csv", row.names = FALSE)
}

if(nrow(ko_lps24_pathways$pathway_counts) > 0) {
  knitr::kable(ko_lps24_pathways$pathway_counts, 
               caption = "Pathway enrichment for KO vs WT in LPS 24hr condition")
  write.csv(ko_lps24_pathways$pathway_counts, 
            "src/tmp_data/ko_lps24_pathway_counts.csv", row.names = FALSE)
}

# -----------------------------------------------------
# 3. Perform pathway enrichment analysis
# -----------------------------------------------------
cat("Performing pathway enrichment analysis...\n")

# Function to calculate enrichment statistics
# We'll use a hypergeometric test to assess enrichment
calculate_pathway_enrichment <- function(query_metabolites, all_metabolites, pathway_dict) {
  # Handle empty input
  if(length(query_metabolites) == 0) {
    return(data.frame(
      pathway = character(),
      hits = integer(),
      expected = numeric(),
      fold_enrichment = numeric(),
      pvalue = numeric(),
      padj = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Get total number of metabolites
  n_total_metabolites <- length(all_metabolites)
  
  # Initialize results dataframe
  enrichment_results <- data.frame(
    pathway = character(),
    hits = integer(),
    expected = numeric(),
    fold_enrichment = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each pathway, calculate enrichment
  for (pathway_name in names(pathway_dict)) {
    # Get metabolites in this pathway
    pathway_metabolites <- pathway_dict[[pathway_name]]
    
    # Count how many pathway metabolites are in our background
    pathway_in_background <- sum(sapply(all_metabolites, function(metab) {
      any(sapply(pathway_metabolites, function(pm) {
        tolower(metab) == tolower(pm) || 
          grepl(tolower(pm), tolower(metab), fixed = TRUE) ||
          grepl(tolower(metab), tolower(pm), fixed = TRUE)
      }))
    }))
    
    # Skip if no pathway metabolites in background
    if (pathway_in_background == 0) next
    
    # Count hits in query
    hits <- sum(sapply(query_metabolites, function(metab) {
      any(sapply(pathway_metabolites, function(pm) {
        tolower(metab) == tolower(pm) || 
          grepl(tolower(pm), tolower(metab), fixed = TRUE) ||
          grepl(tolower(metab), tolower(pm), fixed = TRUE)
      }))
    }))
    
    # Skip if no hits
    if (hits == 0) next
    
    # Calculate expected hits by chance
    n_query_metabolites <- length(query_metabolites)
    expected <- (n_query_metabolites / n_total_metabolites) * pathway_in_background
    
    # Calculate fold enrichment
    fold_enrichment <- hits / expected
    
    # Calculate p-value using hypergeometric test
    pvalue <- stats::phyper(hits - 1, pathway_in_background, 
                          n_total_metabolites - pathway_in_background, 
                          n_query_metabolites, lower.tail = FALSE)
    
    # Add to results
    enrichment_results <- rbind(
      enrichment_results,
      data.frame(
        pathway = pathway_name,
        hits = hits,
        expected = expected,
        fold_enrichment = fold_enrichment,
        pvalue = pvalue,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Add BH-adjusted p-values
  if (nrow(enrichment_results) > 0) {
    enrichment_results$padj <- stats::p.adjust(enrichment_results$pvalue, method = "BH")
    enrichment_results <- enrichment_results[order(enrichment_results$pvalue), ]
  }
  
  return(enrichment_results)
}

# Get all metabolite names as background
all_metabolites <- unique(all_da_results$metabolite)

# Calculate enrichment for each comparison
ko_unstim_enrichment <- calculate_pathway_enrichment(
  ko_unstim_metabolites$metabolite, all_metabolites, metabolic_pathways)
ko_lps6_enrichment <- calculate_pathway_enrichment(
  ko_lps6_metabolites$metabolite, all_metabolites, metabolic_pathways)
ko_lps24_enrichment <- calculate_pathway_enrichment(
  ko_lps24_metabolites$metabolite, all_metabolites, metabolic_pathways)
int_lps6_unstim_enrichment <- calculate_pathway_enrichment(
  int_lps6_unstim$metabolite, all_metabolites, metabolic_pathways)
int_lps24_unstim_enrichment <- calculate_pathway_enrichment(
  int_lps24_unstim$metabolite, all_metabolites, metabolic_pathways)
int_lps24_lps6_enrichment <- calculate_pathway_enrichment(
  int_lps24_lps6$metabolite, all_metabolites, metabolic_pathways)

# Save enrichment results
write.csv(ko_unstim_enrichment, "src/tmp_data/ko_unstim_enrichment.csv", row.names = FALSE)
write.csv(ko_lps6_enrichment, "src/tmp_data/ko_lps6_enrichment.csv", row.names = FALSE)
write.csv(ko_lps24_enrichment, "src/tmp_data/ko_lps24_enrichment.csv", row.names = FALSE)
write.csv(int_lps6_unstim_enrichment, "src/tmp_data/int_lps6_unstim_enrichment.csv", row.names = FALSE)
write.csv(int_lps24_unstim_enrichment, "src/tmp_data/int_lps24_unstim_enrichment.csv", row.names = FALSE)
write.csv(int_lps24_lps6_enrichment, "src/tmp_data/int_lps24_lps6_enrichment.csv", row.names = FALSE)

# Display results for a few key comparisons
if(nrow(ko_unstim_enrichment) > 0) {
  knitr::kable(ko_unstim_enrichment, 
               caption = "Pathway enrichment for KO vs WT in Unstimulated condition",
               digits = 3)
}

if(nrow(ko_lps24_enrichment) > 0) {
  knitr::kable(ko_lps24_enrichment, 
               caption = "Pathway enrichment for KO vs WT in LPS 24hr condition",
               digits = 3)
}

# -----------------------------------------------------
# 4. Visualize pathway enrichment results
# -----------------------------------------------------
cat("Creating visualizations of pathway enrichment...\n")

# Function to create a dot plot of pathway enrichment
create_pathway_dotplot <- function(enrichment_df, title) {
  # Only include pathways with p < 0.1 for better visualization
  plot_data <- enrichment_df %>%
    dplyr::filter(pvalue < 0.1) %>%
    dplyr::mutate(
      log10pval = -log10(pvalue),
      pathway = factor(pathway, levels = rev(pathway))  # For ordering in plot
    )
  
  if (nrow(plot_data) == 0) {
    # If no pathways with p < 0.1, use all pathways
    plot_data <- enrichment_df %>%
      dplyr::mutate(
        log10pval = -log10(pvalue),
        pathway = factor(pathway, levels = rev(pathway))
      )
  }
  
  # Skip if still no data
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold_enrichment, y = pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = hits, color = log10pval)) +
    ggplot2::scale_color_gradient(low = "#807E7D", high = "#4D6D8E") +
    ggplot2::labs(
      title = title,
      x = "Fold Enrichment",
      y = "Pathway",
      size = "Hits",
      color = "-log10(p-value)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# Create dot plots for each comparison
if(nrow(ko_unstim_enrichment) > 0) {
  p1 <- create_pathway_dotplot(ko_unstim_enrichment, 
                             "Pathway Enrichment: KO vs WT (Unstimulated)")
  if(!is.null(p1)) ggplot2::ggsave("src/tmp_figures/ko_unstim_pathway_dotplot.png", p1, width = 10, height = 8)
}

if(nrow(ko_lps6_enrichment) > 0) {
  p2 <- create_pathway_dotplot(ko_lps6_enrichment, 
                             "Pathway Enrichment: KO vs WT (LPS 6hr)")
  if(!is.null(p2)) ggplot2::ggsave("src/tmp_figures/ko_lps6_pathway_dotplot.png", p2, width = 10, height = 8)
}

if(nrow(ko_lps24_enrichment) > 0) {
  p3 <- create_pathway_dotplot(ko_lps24_enrichment, 
                             "Pathway Enrichment: KO vs WT (LPS 24hr)")
  if(!is.null(p3)) ggplot2::ggsave("src/tmp_figures/ko_lps24_pathway_dotplot.png", p3, width = 10, height = 8)
}

if(nrow(int_lps6_unstim_enrichment) > 0) {
  p4 <- create_pathway_dotplot(int_lps6_unstim_enrichment, 
                             "Pathway Enrichment: Interaction (LPS 6hr vs Unstimulated)")
  if(!is.null(p4)) ggplot2::ggsave("src/tmp_figures/int_lps6_unstim_pathway_dotplot.png", p4, width = 10, height = 8)
}

if(nrow(int_lps24_unstim_enrichment) > 0) {
  p5 <- create_pathway_dotplot(int_lps24_unstim_enrichment, 
                             "Pathway Enrichment: Interaction (LPS 24hr vs Unstimulated)")
  if(!is.null(p5)) ggplot2::ggsave("src/tmp_figures/int_lps24_unstim_pathway_dotplot.png", p5, width = 10, height = 8)
}

if(nrow(int_lps24_lps6_enrichment) > 0) {
  p6 <- create_pathway_dotplot(int_lps24_lps6_enrichment, 
                             "Pathway Enrichment: Interaction (LPS 24hr vs LPS 6hr)")
  if(!is.null(p6)) ggplot2::ggsave("src/tmp_figures/int_lps24_lps6_pathway_dotplot.png", p6, width = 10, height = 8)
}

# Function to create a bar plot of pathway enrichment
create_pathway_barplot <- function(enrichment_df, title) {
  # Only include pathways with p < 0.1 for better visualization
  plot_data <- enrichment_df %>%
    dplyr::filter(pvalue < 0.1) %>%
    dplyr::arrange(pvalue)
  
  if (nrow(plot_data) == 0) {
    # If no pathways with p < 0.1, use top 10 pathways
    plot_data <- enrichment_df %>%
      dplyr::arrange(pvalue) %>%
      head(10)
  }
  
  # Skip if still no data
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # For better visualization, limit to top 15 pathways
  if (nrow(plot_data) > 15) {
    plot_data <- plot_data %>% head(15)
  }
  
  # Create factor for ordering
  plot_data$pathway <- factor(plot_data$pathway, levels = rev(plot_data$pathway))
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold_enrichment, y = pathway)) +
    ggplot2::geom_bar(stat = "identity", fill = "#4D6D8E") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("p=%.3g", pvalue)), 
                      hjust = -0.1, size = 3.5) +
    ggplot2::labs(
      title = title,
      x = "Fold Enrichment",
      y = "Pathway"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# Create bar plots for each comparison
if(nrow(ko_unstim_enrichment) > 0) {
  b1 <- create_pathway_barplot(ko_unstim_enrichment, 
                             "Pathway Enrichment: KO vs WT (Unstimulated)")
  if(!is.null(b1)) ggplot2::ggsave("src/tmp_figures/ko_unstim_pathway_barplot.png", b1, width = 10, height = 8)
}

if(nrow(ko_lps6_enrichment) > 0) {
  b2 <- create_pathway_barplot(ko_lps6_enrichment, 
                             "Pathway Enrichment: KO vs WT (LPS 6hr)")
  if(!is.null(b2)) ggplot2::ggsave("src/tmp_figures/ko_lps6_pathway_barplot.png", b2, width = 10, height = 8)
}

if(nrow(ko_lps24_enrichment) > 0) {
  b3 <- create_pathway_barplot(ko_lps24_enrichment, 
                             "Pathway Enrichment: KO vs WT (LPS 24hr)")
  if(!is.null(b3)) ggplot2::ggsave("src/tmp_figures/ko_lps24_pathway_barplot.png", b3, width = 10, height = 8)
}

# -----------------------------------------------------
# 5. Analyze metabolite direction changes within pathways
# -----------------------------------------------------
cat("Analyzing metabolite direction changes within pathways...\n")

# Function to annotate metabolites with their fold change direction
annotate_pathway_metabolites <- function(enrichment_results, metabolite_data, pathway_dict) {
  # Initialize results list
  pathway_metabolite_changes <- list()
  
  if(nrow(enrichment_results) == 0 || nrow(metabolite_data) == 0) {
    return(pathway_metabolite_changes)
  }
  
  for (i in 1:nrow(enrichment_results)) {
    pathway_name <- enrichment_results$pathway[i]
    pathway_metabolites <- pathway_dict[[pathway_name]]
    
    # Find metabolites in this pathway from our data
    pathway_hits <- sapply(metabolite_data$metabolite, function(metab) {
      any(sapply(pathway_metabolites, function(pm) {
        tolower(metab) == tolower(pm) || 
          grepl(tolower(pm), tolower(metab), fixed = TRUE) ||
          grepl(tolower(metab), tolower(pm), fixed = TRUE)
      }))
    })
    
    # Get the matching metabolites and their info
    pathway_metabolite_data <- metabolite_data[pathway_hits, ]
    
    # Add to results
    pathway_metabolite_changes[[pathway_name]] <- pathway_metabolite_data
  }
  
  return(pathway_metabolite_changes)
}

# Apply function to each comparison
ko_unstim_pathway_changes <- annotate_pathway_metabolites(
  ko_unstim_enrichment, ko_unstim_metabolites, metabolic_pathways)
ko_lps6_pathway_changes <- annotate_pathway_metabolites(
  ko_lps6_enrichment, ko_lps6_metabolites, metabolic_pathways)
ko_lps24_pathway_changes <- annotate_pathway_metabolites(
  ko_lps24_enrichment, ko_lps24_metabolites, metabolic_pathways)

# Save the results
saveRDS(ko_unstim_pathway_changes, "src/tmp_data/ko_unstim_pathway_changes.rds")
saveRDS(ko_lps6_pathway_changes, "src/tmp_data/ko_lps6_pathway_changes.rds")
saveRDS(ko_lps24_pathway_changes, "src/tmp_data/ko_lps24_pathway_changes.rds")

# Create a summary of up/down regulated metabolites per pathway
create_pathway_direction_summary <- function(pathway_changes) {
  if(length(pathway_changes) == 0) {
    return(data.frame(
      pathway = character(),
      total_metabolites = integer(),
      upregulated = integer(),
      downregulated = integer(),
      percent_up = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  summary_df <- data.frame(
    pathway = character(),
    total_metabolites = integer(),
    upregulated = integer(),
    downregulated = integer(),
    percent_up = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pathway_name in names(pathway_changes)) {
    pathway_data <- pathway_changes[[pathway_name]]
    if (nrow(pathway_data) > 0) {
      total <- nrow(pathway_data)
      up <- sum(pathway_data$logFC > 0)
      down <- sum(pathway_data$logFC < 0)
      percent_up <- (up / total) * 100
      
      summary_df <- rbind(
        summary_df,
        data.frame(
          pathway = pathway_name,
          total_metabolites = total,
          upregulated = up,
          downregulated = down,
          percent_up = percent_up,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # Sort by number of total metabolites
  if (nrow(summary_df) > 0) {
    summary_df <- summary_df[order(-summary_df$total_metabolites), ]
  }
  
  return(summary_df)
}

# Create direction summaries
ko_unstim_direction <- create_pathway_direction_summary(ko_unstim_pathway_changes)
ko_lps6_direction <- create_pathway_direction_summary(ko_lps6_pathway_changes)
ko_lps24_direction <- create_pathway_direction_summary(ko_lps24_pathway_changes)

# Save and display the summaries
write.csv(ko_unstim_direction, "src/tmp_data/ko_unstim_direction_summary.csv", row.names = FALSE)
write.csv(ko_lps6_direction, "src/tmp_data/ko_lps6_direction_summary.csv", row.names = FALSE)
write.csv(ko_lps24_direction, "src/tmp_data/ko_lps24_direction_summary.csv", row.names = FALSE)

if(nrow(ko_unstim_direction) > 0) {
  knitr::kable(ko_unstim_direction, 
               caption = "Metabolite direction changes in pathways (KO vs WT, Unstimulated)")
}

if(nrow(ko_lps24_direction) > 0) {
  knitr::kable(ko_lps24_direction, 
               caption = "Metabolite direction changes in pathways (KO vs WT, LPS 24hr)")
}

# -----------------------------------------------------
# 6. Create network visualizations of pathway-metabolite relationships
# -----------------------------------------------------
cat("Creating network visualizations of pathway-metabolite relationships...\n")

# Function to prepare data for a network visualization
prepare_network_data <- function(pathway_changes, enrichment_results, max_pathways = 5, max_metabolites = 30) {
  if(length(pathway_changes) == 0 || nrow(enrichment_results) == 0) {
    return(list(
      nodes = data.frame(
        id = character(),
        label = character(),
        type = character(),
        logFC = numeric(),
        pvalue = numeric(),
        stringsAsFactors = FALSE
      ),
      edges = data.frame(
        from = character(),
        to = character(),
        width = numeric(),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # Select top pathways based on p-value
  top_pathways <- head(enrichment_results[order(enrichment_results$pvalue), ], max_pathways)$pathway
  
  # Initialize node and edge dataframes
  nodes <- data.frame(
    id = character(),
    label = character(),
    type = character(),
    logFC = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  edges <- data.frame(
    from = character(),
    to = character(),
    width = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add pathway nodes
  for (pathway in top_pathways) {
    pathway_index <- which(enrichment_results$pathway == pathway)
    if (length(pathway_index) > 0) {
      pathway_pvalue <- enrichment_results$pvalue[pathway_index]
      nodes <- rbind(
        nodes,
        data.frame(
          id = pathway,
          label = pathway,
          type = "pathway",
          logFC = NA,
          pvalue = pathway_pvalue,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # Add metabolite nodes and edges
  for (pathway in top_pathways) {
    if (!is.null(pathway_changes[[pathway]])) {
      # Get metabolites for this pathway
      pathway_metabolites <- pathway_changes[[pathway]]
      
      # Limit to top metabolites by significance
      if (nrow(pathway_metabolites) > max_metabolites) {
        pathway_metabolites <- pathway_metabolites[order(pathway_metabolites$adj.P.Val), ][1:max_metabolites, ]
      }
      
      # Add nodes
      for (i in 1:nrow(pathway_metabolites)) {
        metab <- pathway_metabolites$metabolite[i]
        # Check if already in nodes
        if (!metab %in% nodes$id) {
          nodes <- rbind(
            nodes,
            data.frame(
              id = metab,
              label = metab,
              type = "metabolite",
              logFC = pathway_metabolites$logFC[i],
              pvalue = pathway_metabolites$P.Value[i],
              stringsAsFactors = FALSE
            )
          )
        }
        
        # Add edge
        edges <- rbind(
          edges,
          data.frame(
            from = pathway,
            to = metab,
            width = 1,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  return(list(nodes = nodes, edges = edges))
}

# Prepare network data for visualization
ko_unstim_network <- prepare_network_data(ko_unstim_pathway_changes, ko_unstim_enrichment)
ko_lps6_network <- prepare_network_data(ko_lps6_pathway_changes, ko_lps6_enrichment)
ko_lps24_network <- prepare_network_data(ko_lps24_pathway_changes, ko_lps24_enrichment)

# Save network data
saveRDS(ko_unstim_network, "src/tmp_data/ko_unstim_network.rds")
saveRDS(ko_lps6_network, "src/tmp_data/ko_lps6_network.rds")
saveRDS(ko_lps24_network, "src/tmp_data/ko_lps24_network.rds")

# Function to create a custom visualization of the network
visualize_network <- function(network_data, title) {
  if(nrow(network_data$nodes) == 0 || nrow(network_data$edges) == 0) {
    return(NULL)
  }
  
  nodes <- network_data$nodes
  edges <- network_data$edges
  
  # Create a long format dataset for plotting
  plot_edges <- edges %>%
    dplyr::left_join(nodes %>% dplyr::select(id, type, logFC), by = c("from" = "id")) %>%
    dplyr::rename(from_type = type, from_logFC = logFC) %>%
    dplyr::left_join(nodes %>% dplyr::select(id, type, logFC), by = c("to" = "id")) %>%
    dplyr::rename(to_type = type, to_logFC = logFC)
  
  # Create a data frame of pathway-metabolite relationships for plotting
  pathway_metab_data <- plot_edges %>%
    dplyr::filter(from_type == "pathway", to_type == "metabolite") %>%
    dplyr::left_join(nodes %>% 
                  dplyr::filter(type == "metabolite") %>% 
                  dplyr::select(id, logFC), 
                by = c("to" = "id"))
  
  # Skip if no data
  if(nrow(pathway_metab_data) == 0) {
    return(NULL)
  }
  
  # Group by pathway and metabolite, summarizing logFC
  grouped_data <- pathway_metab_data %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarise(logFC = dplyr::first(logFC), .groups = "drop")
  
  # Order pathways by number of metabolites
  pathway_counts <- grouped_data %>% 
    dplyr::group_by(from) %>% 
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  pathway_order <- pathway_counts %>% 
    dplyr::arrange(desc(count)) %>% 
    dplyr::pull(from)
  
  grouped_data$from <- factor(grouped_data$from, levels = pathway_order)
  
  # Create the plot
  p <- ggplot2::ggplot(grouped_data, ggplot2::aes(x = from, y = to)) +
    ggplot2::geom_tile(ggplot2::aes(fill = logFC), color = "white") +
    ggplot2::scale_fill_gradient2(low = "#8e6e4d", mid = "white", high = "#4D6D8E", midpoint = 0) +
    ggplot2::labs(
      title = title,
      x = "Pathway",
      y = "Metabolite",
      fill = "log2FC"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# Create network visualizations
n1 <- visualize_network(ko_unstim_network, "Pathway-Metabolite Network (KO vs WT, Unstimulated)")
if(!is.null(n1)) ggplot2::ggsave("src/tmp_figures/ko_unstim_network.png", n1, width = 12, height = 10)

n2 <- visualize_network(ko_lps6_network, "Pathway-Metabolite Network (KO vs WT, LPS 6hr)")
if(!is.null(n2)) ggplot2::ggsave("src/tmp_figures/ko_lps6_network.png", n2, width = 12, height = 10)

n3 <- visualize_network(ko_lps24_network, "Pathway-Metabolite Network (KO vs WT, LPS 24hr)")
if(!is.null(n3)) ggplot2::ggsave("src/tmp_figures/ko_lps24_network.png", n3, width = 12, height = 10)

# -----------------------------------------------------
# 7. Create a comprehensive pathway analysis summary
# -----------------------------------------------------
cat("Creating a comprehensive pathway analysis summary...\n")

# Combine pathway enrichment results across conditions
all_pathways <- unique(c(
  ko_unstim_enrichment$pathway,
  ko_lps6_enrichment$pathway,
  ko_lps24_enrichment$pathway
))

if(length(all_pathways) > 0) {
  pathway_summary <- data.frame(
    pathway = all_pathways,
    stringsAsFactors = FALSE
  )

  # Add enrichment p-values for each condition
  add_pvalues <- function(summary_df, enrichment_df, condition) {
    p_col <- paste0("pvalue_", condition)
    hits_col <- paste0("hits_", condition)
    
    summary_df[[p_col]] <- NA
    summary_df[[hits_col]] <- NA
    
    for (i in 1:nrow(enrichment_df)) {
      pathway <- enrichment_df$pathway[i]
      idx <- which(summary_df$pathway == pathway)
      if (length(idx) > 0) {
        summary_df[[p_col]][idx] <- enrichment_df$pvalue[i]
        summary_df[[hits_col]][idx] <- enrichment_df$hits[i]
      }
    }
    
    return(summary_df)
  }

  pathway_summary <- add_pvalues(pathway_summary, ko_unstim_enrichment, "unstim")
  pathway_summary <- add_pvalues(pathway_summary, ko_lps6_enrichment, "lps6")
  pathway_summary <- add_pvalues(pathway_summary, ko_lps24_enrichment, "lps24")

  # Calculate a combined significance score
  pathway_summary$combined_score <- rowMeans(
    -log10(pathway_summary[, c("pvalue_unstim", "pvalue_lps6", "pvalue_lps24")]), 
    na.rm = TRUE
  )

  # Sort by combined score
  pathway_summary <- pathway_summary[order(-pathway_summary$combined_score), ]

  # Save the summary
  write.csv(pathway_summary, "src/tmp_data/pathway_summary.csv", row.names = FALSE)

  # Display the top pathways
  knitr::kable(head(pathway_summary, 10), 
               caption = "Top pathways across all conditions",
               digits = 3)

  # Create a heatmap of pathway significance across conditions
  create_pathway_heatmap <- function(summary_df, top_n = 15) {
    # Select top pathways
    plot_data <- head(summary_df, min(top_n, nrow(summary_df)))
    
    # Convert to long format for plotting
    plot_data_long <- plot_data %>%
      tidyr::pivot_longer(
        cols = c("pvalue_unstim", "pvalue_lps6", "pvalue_lps24"),
        names_to = "condition",
        values_to = "pvalue"
      ) %>%
      dplyr::mutate(
        condition = stringr::str_replace(condition, "pvalue_", ""),
        condition = factor(condition, levels = c("unstim", "lps6", "lps24")),
        neg_log10_p = -log10(pvalue),
        pathway = factor(pathway, levels = rev(plot_data$pathway))
      )
    
    # Create heatmap
    p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = condition, y = pathway, fill = neg_log10_p)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "#807E7D", high = "#4D6D8E", 
                                 na.value = "white", name = "-log10(p-value)") +
      ggplot2::labs(
        title = "Pathway Significance Across Conditions",
        x = "Condition",
        y = "Pathway"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(size = 12, face = "bold")
      )
    
    return(p)
  }

  # Create the heatmap
  h1 <- create_pathway_heatmap(pathway_summary)
  if(!is.null(h1)) ggplot2::ggsave("src/tmp_figures/pathway_significance_heatmap.png", h1, width = 10, height = 8)
}

# -----------------------------------------------------
# 8. Final summary of biological findings
# -----------------------------------------------------
cat("Creating final summary of biological findings...\n")

# Create a summary of key biological findings
biological_summary <- data.frame(
  finding = c(
    "Top enriched pathway in unstimulated condition",
    "Top enriched pathway in LPS 6hr",
    "Top enriched pathway in LPS 24hr",
    "Consistently enriched pathways across conditions",
    "Pathways showing treatment-specific enrichment",
    "Key upregulated metabolite in KO vs WT (Unstimulated)",
    "Key downregulated metabolite in KO vs WT (Unstimulated)",
    "Key upregulated metabolite in KO vs WT (LPS 24hr)",
    "Key downregulated metabolite in KO vs WT (LPS 24hr)"
  ),
  description = character(9),
  stringsAsFactors = FALSE
)

# Fill in the findings
# Top pathway in unstimulated
if (nrow(ko_unstim_enrichment) > 0) {
  biological_summary$description[1] <- ko_unstim_enrichment$pathway[1]
} else {
  biological_summary$description[1] <- "No significant pathway enrichment"
}

# Top pathway in LPS 6hr
if (nrow(ko_lps6_enrichment) > 0) {
  biological_summary$description[2] <- ko_lps6_enrichment$pathway[1]
} else {
  biological_summary$description[2] <- "No significant pathway enrichment"
}

# Top pathway in LPS 24hr
if (nrow(ko_lps24_enrichment) > 0) {
  biological_summary$description[3] <- ko_lps24_enrichment$pathway[1]
} else {
  biological_summary$description[3] <- "No significant pathway enrichment"
}

# Consistently enriched pathways
if(exists("pathway_summary") && nrow(pathway_summary) > 0) {
  consistent_pathways <- pathway_summary %>%
    dplyr::filter(!is.na(pvalue_unstim) & !is.na(pvalue_lps6) & !is.na(pvalue_lps24)) %>%
    dplyr::filter(pvalue_unstim < 0.1 & pvalue_lps6 < 0.1 & pvalue_lps24 < 0.1) %>%
    dplyr::pull(pathway)
  
  if (length(consistent_pathways) > 0) {
    biological_summary$description[4] <- paste(consistent_pathways, collapse = ", ")
  } else {
    biological_summary$description[4] <- "No consistently enriched pathways"
  }
  
  # Treatment-specific pathways
  unstim_specific <- pathway_summary %>%
    dplyr::filter(!is.na(pvalue_unstim) & pvalue_unstim < 0.05) %>%
    dplyr::filter(is.na(pvalue_lps6) | pvalue_lps6 > 0.1) %>%
    dplyr::filter(is.na(pvalue_lps24) | pvalue_lps24 > 0.1) %>%
    dplyr::pull(pathway)
  
  lps24_specific <- pathway_summary %>%
    dplyr::filter(!is.na(pvalue_lps24) & pvalue_lps24 < 0.05) %>%
    dplyr::filter(is.na(pvalue_unstim) | pvalue_unstim > 0.1) %>%
    dplyr::filter(is.na(pvalue_lps6) | pvalue_lps6 > 0.1) %>%
    dplyr::pull(pathway)
  
  treatment_specific <- c()
  if(length(unstim_specific) > 0) {
    treatment_specific <- c(treatment_specific, paste("Unstimulated:", paste(unstim_specific, collapse = ", ")))
  }
  if(length(lps24_specific) > 0) {
    treatment_specific <- c(treatment_specific, paste("LPS 24hr:", paste(lps24_specific, collapse = ", ")))
  }
  
  if(length(treatment_specific) > 0) {
    biological_summary$description[5] <- paste(treatment_specific, collapse = "; ")
  } else {
    biological_summary$description[5] <- "No treatment-specific pathways identified"
  }
} else {
  biological_summary$description[4] <- "Insufficient data for pathway consistency analysis"
  biological_summary$description[5] <- "Insufficient data for treatment-specific pathway analysis"
}

# Key metabolites in unstimulated condition
if (nrow(ko_unstim_metabolites) > 0) {
  up_metabolites <- ko_unstim_metabolites %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(adj.P.Val)
  
  down_metabolites <- ko_unstim_metabolites %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::arrange(adj.P.Val)
  
  if (nrow(up_metabolites) > 0) {
    up_metabolite <- up_metabolites[1,]
    biological_summary$description[6] <- paste0(
      up_metabolite$metabolite, " (logFC = ", round(up_metabolite$logFC, 2), 
      ", padj = ", sprintf("%.2e", up_metabolite$adj.P.Val), ")"
    )
  } else {
    biological_summary$description[6] <- "No significantly upregulated metabolites"
  }
  
  if (nrow(down_metabolites) > 0) {
    down_metabolite <- down_metabolites[1,]
    biological_summary$description[7] <- paste0(
      down_metabolite$metabolite, " (logFC = ", round(down_metabolite$logFC, 2), 
      ", padj = ", sprintf("%.2e", down_metabolite$adj.P.Val), ")"
    )
  } else {
    biological_summary$description[7] <- "No significantly downregulated metabolites"
  }
} else {
  biological_summary$description[6] <- "No significantly upregulated metabolites"
  biological_summary$description[7] <- "No significantly downregulated metabolites"
}

# Key metabolites in LPS 24hr condition
if (nrow(ko_lps24_metabolites) > 0) {
  up_metabolites <- ko_lps24_metabolites %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(adj.P.Val)
  
  down_metabolites <- ko_lps24_metabolites %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::arrange(adj.P.Val)
  
  if (nrow(up_metabolites) > 0) {
    up_metabolite <- up_metabolites[1,]
    biological_summary$description[8] <- paste0(
      up_metabolite$metabolite, " (logFC = ", round(up_metabolite$logFC, 2), 
      ", padj = ", sprintf("%.2e", up_metabolite$adj.P.Val), ")"
    )
  } else {
    biological_summary$description[8] <- "No significantly upregulated metabolites"
  }
  
  if (nrow(down_metabolites) > 0) {
    down_metabolite <- down_metabolites[1,]
    biological_summary$description[9] <- paste0(
      down_metabolite$metabolite, " (logFC = ", round(down_metabolite$logFC, 2), 
      ", padj = ", sprintf("%.2e", down_metabolite$adj.P.Val), ")"
    )
  } else {
    biological_summary$description[9] <- "No significantly downregulated metabolites"
  }
} else {
  biological_summary$description[8] <- "No significantly upregulated metabolites"
  biological_summary$description[9] <- "No significantly downregulated metabolites"
}

# Save and display the biological summary
write.csv(biological_summary, "src/tmp_data/biological_summary.csv", row.names = FALSE)
knitr::kable(biological_summary, caption = "Summary of key biological findings")

# Save all pathway analysis results in a single list object
pathway_analysis_results <- list(
  pathway_mappings = list(
    ko_unstim = ko_unstim_pathways,
    ko_lps6 = ko_lps6_pathways,
    ko_lps24 = ko_lps24_pathways
  ),
  enrichment_results = list(
    ko_unstim = ko_unstim_enrichment,
    ko_lps6 = ko_lps6_enrichment,
    ko_lps24 = ko_lps24_enrichment,
    int_lps6_unstim = int_lps6_unstim_enrichment,
    int_lps24_unstim = int_lps24_unstim_enrichment,
    int_lps24_lps6 = int_lps24_lps6_enrichment
  ),
  pathway_changes = list(
    ko_unstim = ko_unstim_pathway_changes,
    ko_lps6 = ko_lps6_pathway_changes,
    ko_lps24 = ko_lps24_pathway_changes
  ),
  network_data = list(
    ko_unstim = ko_unstim_network,
    ko_lps6 = ko_lps6_network,
    ko_lps24 = ko_lps24_network
  ),
  summary = list(
    direction_changes = list(
      ko_unstim = ko_unstim_direction,
      ko_lps6 = ko_lps6_direction,
      ko_lps24 = ko_lps24_direction
    ),
    biological_findings = biological_summary
  )
)

# Add combined pathway summary if it exists
if(exists("pathway_summary")) {
  pathway_analysis_results$summary$combined_pathways <- pathway_summary
}

# Save the comprehensive results
saveRDS(pathway_analysis_results, "src/tmp_data/pathway_analysis_results.rds")

cat("Pathway analysis completed successfully!\n")
