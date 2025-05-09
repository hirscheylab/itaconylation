
# Load required libraries
library(stats)
library(limma)
library(dplyr)
library(car)
library(ggplot2)
library(tidyr)
library(knitr)

# Load the cytokine data
cytokine_data <- readRDS("src/data/cytokine_data.rds")

# Identify cytokine columns (all columns after protein_mg_ml)
cytokine_cols <- names(cytokine_data)[6:ncol(cytokine_data)]

# Log2 transform the data (adding small constant to avoid log(0))
small_constant <- 0.01
cytokine_log2 <- cytokine_data %>%
  mutate(across(all_of(cytokine_cols), ~log2(.x + small_constant)))

# Check for missing values and use only complete samples
complete_samples <- cytokine_log2 %>%
  filter(!is.na(Genotype_clean) & !is.na(LPS_treatment_hours))

cat("Using", nrow(complete_samples), "of", nrow(cytokine_log2), "samples with complete metadata\n\n")

# Verify sample distribution
sample_counts <- complete_samples %>%
  group_by(Genotype_clean, LPS_treatment_hours) %>%
  summarize(Count = n(), .groups = "drop")

knitr::kable(sample_counts, caption = "Sample distribution across experimental conditions")

#-------------------------------------------------
# 1. Two-way ANOVA for each cytokine
#-------------------------------------------------

# Initialize lists to store results
anova_results <- list()
anova_summary <- data.frame(
  Cytokine = character(),
  Genotype_p = numeric(),
  Treatment_p = numeric(),
  Interaction_p = numeric(),
  Genotype_significant = logical(),
  Treatment_significant = logical(),
  Interaction_significant = logical(),
  stringsAsFactors = FALSE
)

# Perform two-way ANOVA for each cytokine
for (cytokine in cytokine_cols) {
  # Create formula for the model
  formula <- as.formula(paste(cytokine, "~ Genotype_clean * LPS_treatment_hours"))
  
  # Fit the model
  model <- try(stats::aov(formula, data = complete_samples), silent = TRUE)
  
  # Check if model fitting was successful
  if (!inherits(model, "try-error")) {
    # Store the model
    anova_results[[cytokine]] <- model
    
    # Get ANOVA table
    anova_table <- summary(model)[[1]]
    
    # Extract p-values
    genotype_p <- anova_table["Genotype_clean", "Pr(>F)"]
    treatment_p <- anova_table["LPS_treatment_hours", "Pr(>F)"]
    interaction_p <- anova_table["Genotype_clean:LPS_treatment_hours", "Pr(>F)"]
    
    # Add to summary table
    anova_summary <- rbind(anova_summary, data.frame(
      Cytokine = cytokine,
      Genotype_p = genotype_p,
      Treatment_p = treatment_p,
      Interaction_p = interaction_p,
      Genotype_significant = genotype_p < 0.05,
      Treatment_significant = treatment_p < 0.05,
      Interaction_significant = interaction_p < 0.05,
      stringsAsFactors = FALSE
    ))
  } else {
    warning(paste("ANOVA model fitting failed for cytokine:", cytokine))
  }
}

# Apply multiple testing correction to p-values
anova_summary$Genotype_p_adj <- stats::p.adjust(anova_summary$Genotype_p, method = "BH")
anova_summary$Treatment_p_adj <- stats::p.adjust(anova_summary$Treatment_p, method = "BH")
anova_summary$Interaction_p_adj <- stats::p.adjust(anova_summary$Interaction_p, method = "BH")

# Update significance flags based on adjusted p-values
anova_summary$Genotype_significant_adj <- anova_summary$Genotype_p_adj < 0.05
anova_summary$Treatment_significant_adj <- anova_summary$Treatment_p_adj < 0.05
anova_summary$Interaction_significant_adj <- anova_summary$Interaction_p_adj < 0.05

# Format p-values for display
anova_display <- anova_summary %>%
  mutate(across(ends_with("_p"), ~format.pval(., digits = 3)),
         across(ends_with("_p_adj"), ~format.pval(., digits = 3)))

# Display ANOVA results
knitr::kable(anova_display, 
             caption = "Two-way ANOVA results for each cytokine (Genotype × Treatment)")

# Save the ANOVA results
saveRDS(anova_results, "src/tmp_data/anova_models.rds")
saveRDS(anova_summary, "src/tmp_data/anova_summary.rds")

#-------------------------------------------------
# 2. Limma analysis for robust statistical testing
#-------------------------------------------------

# Prepare data for limma analysis
# Create a matrix of log2-transformed cytokine values
expr_matrix <- as.matrix(complete_samples[, cytokine_cols])
rownames(expr_matrix) <- complete_samples$ID

# Create design matrix
design <- stats::model.matrix(~ Genotype_clean * LPS_treatment_hours, data = complete_samples)

# Fit linear model
fit <- limma::lmFit(t(expr_matrix), design)

# Apply empirical Bayes method to moderate standard errors
fit_eb <- limma::eBayes(fit)

# Extract results for each coefficient
limma_results <- list()
limma_summary <- data.frame()

# Get coefficient names
coef_names <- colnames(design)

# For each coefficient (except intercept)
for (i in 2:length(coef_names)) {
  coef <- coef_names[i]
  
  # Extract results for this coefficient
  top_table <- limma::topTable(fit_eb, coef = i, number = Inf, sort.by = "p")
  limma_results[[coef]] <- top_table
  
  # Add coefficient name
  top_table$Coefficient <- coef
  top_table$Cytokine <- rownames(top_table)
  
  # Append to summary
  limma_summary <- rbind(limma_summary, top_table)
}

# Add significance indicator
limma_summary$Significant <- limma_summary$adj.P.Val < 0.05

# Display the top significant results
top_limma_results <- limma_summary %>%
  filter(Significant) %>%
  arrange(adj.P.Val) %>%
  select(Cytokine, Coefficient, logFC, P.Value, adj.P.Val) %>%
  head(20)

knitr::kable(top_limma_results, 
             caption = "Top significant results from limma analysis",
             digits = 3)

# Save limma results
saveRDS(limma_results, "src/tmp_data/limma_results.rds")
saveRDS(limma_summary, "src/tmp_data/limma_summary.rds")

#-------------------------------------------------
# 3. Post-hoc tests for significant interactions
#-------------------------------------------------

# Identify cytokines with significant interactions after correction
sig_interactions <- anova_summary %>%
  filter(Interaction_significant_adj) %>%
  pull(Cytokine)

# Create function to perform Tukey HSD post-hoc test
perform_tukey <- function(cytokine, data) {
  formula <- as.formula(paste(cytokine, "~ Genotype_clean * LPS_treatment_hours"))
  model <- stats::aov(formula, data = data)
  
  # Perform Tukey HSD test
  tukey_result <- stats::TukeyHSD(model, which = c("Genotype_clean", "LPS_treatment_hours", 
                                           "Genotype_clean:LPS_treatment_hours"))
  return(tukey_result)
}

# Perform post-hoc tests for cytokines with significant interactions
posthoc_results <- list()

if (length(sig_interactions) > 0) {
  for (cytokine in sig_interactions) {
    posthoc_results[[cytokine]] <- try(perform_tukey(cytokine, complete_samples), silent = TRUE)
  }
  
  # Extract and format significant pairwise comparisons
  posthoc_summary <- data.frame()
  
  for (cytokine in names(posthoc_results)) {
    if (!inherits(posthoc_results[[cytokine]], "try-error")) {
      # Extract interaction comparisons
      interactions <- as.data.frame(posthoc_results[[cytokine]]$`Genotype_clean:LPS_treatment_hours`)
      interactions$comparison <- rownames(interactions)
      interactions$cytokine <- cytokine
      
      # Only keep significant comparisons
      sig_comparisons <- interactions %>%
        filter(`p adj` < 0.05) %>%
        select(cytokine, comparison, diff, `p adj`)
      
      posthoc_summary <- rbind(posthoc_summary, sig_comparisons)
    }
  }
  
  if (nrow(posthoc_summary) > 0) {
    # Display significant pairwise comparisons
    knitr::kable(posthoc_summary, 
                caption = "Significant pairwise comparisons from Tukey HSD post-hoc tests",
                col.names = c("Cytokine", "Comparison", "Difference", "Adjusted p-value"),
                digits = 3)
  } else {
    cat("No significant pairwise comparisons found in post-hoc tests.\n")
  }
} else {
  cat("No cytokines with significant interaction terms after multiple testing correction.\n")
}

# Save post-hoc results
saveRDS(posthoc_results, "src/tmp_data/posthoc_results.rds")

#-------------------------------------------------
# 4. Generate tables of significantly different cytokines
#-------------------------------------------------

# A. Cytokines affected by genotype (main effect)
genotype_effect <- anova_summary %>%
  filter(Genotype_significant_adj) %>%
  arrange(Genotype_p_adj) %>%
  select(Cytokine, p_value = Genotype_p, p_adj = Genotype_p_adj)

# B. Cytokines affected by treatment (main effect)
treatment_effect <- anova_summary %>%
  filter(Treatment_significant_adj) %>%
  arrange(Treatment_p_adj) %>%
  select(Cytokine, p_value = Treatment_p, p_adj = Treatment_p_adj)

# C. Cytokines with genotype-treatment interaction
interaction_effect <- anova_summary %>%
  filter(Interaction_significant_adj) %>%
  arrange(Interaction_p_adj) %>%
  select(Cytokine, p_value = Interaction_p, p_adj = Interaction_p_adj)

# Display tables
if (nrow(genotype_effect) > 0) {
  knitr::kable(genotype_effect, 
               caption = "Cytokines significantly affected by genotype (FDR < 0.05)",
               digits = 3)
} else {
  cat("No cytokines significantly affected by genotype after multiple testing correction.\n")
}

if (nrow(treatment_effect) > 0) {
  knitr::kable(treatment_effect, 
               caption = "Cytokines significantly affected by LPS treatment (FDR < 0.05)",
               digits = 3)
} else {
  cat("No cytokines significantly affected by treatment after multiple testing correction.\n")
}

if (nrow(interaction_effect) > 0) {
  knitr::kable(interaction_effect, 
               caption = "Cytokines with significant genotype-treatment interaction (FDR < 0.05)",
               digits = 3)
} else {
  cat("No cytokines with significant genotype-treatment interaction after multiple testing correction.\n")
}

# Create visualizations of significant results
for (cytokine in unique(c(genotype_effect$Cytokine, interaction_effect$Cytokine))) {
  p <- ggplot(complete_samples, aes_string(x = "LPS_treatment_hours", y = cytokine, color = "Genotype_clean", group = "Genotype_clean")) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
    theme_bw() +
    labs(title = paste("Interaction plot for", cytokine),
         x = "LPS treatment (hours)",
         y = paste("log2", cytokine, "expression"),
         color = "Genotype")
  
  # Save the plot
  ggsave(paste0("src/tmp_figures/significant_", cytokine, ".png"), p, width = 8, height = 6, dpi = 300)
}

# Prepare summary of significant findings
all_significant <- anova_summary %>%
  filter(Genotype_significant_adj | Treatment_significant_adj | Interaction_significant_adj) %>%
  mutate(
    Significance_type = case_when(
      Interaction_significant_adj ~ "Interaction",
      Genotype_significant_adj & Treatment_significant_adj ~ "Both main effects",
      Genotype_significant_adj ~ "Genotype only",
      Treatment_significant_adj ~ "Treatment only",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(Significance_type, Genotype_p_adj)

# Save comprehensive results
saveRDS(all_significant, "src/tmp_data/significant_cytokines.rds")

# Create a summary of statistical analysis
stats_summary <- list(
  total_cytokines = length(cytokine_cols),
  genotype_significant = sum(anova_summary$Genotype_significant_adj),
  treatment_significant = sum(anova_summary$Treatment_significant_adj),
  interaction_significant = sum(anova_summary$Interaction_significant_adj),
  any_significant = nrow(all_significant),
  analysis_methods = c("Two-way ANOVA", "limma with empirical Bayes", "Tukey HSD post-hoc tests"),
  multiple_testing_correction = "Benjamini-Hochberg FDR"
)

# Save statistical analysis summary
saveRDS(stats_summary, "src/tmp_data/stats_summary.rds")

# Create structure description for all objects
objects_structure <- list(
  anova_models = "List of ANOVA model objects for each cytokine",
  anova_summary = "Data frame with ANOVA p-values and significance flags for each cytokine",
  limma_results = "List of limma topTable results for each coefficient",
  limma_summary = "Combined data frame of all limma results",
  posthoc_results = "List of Tukey HSD results for cytokines with significant interactions",
  significant_cytokines = "Data frame of all cytokines with any significant effects",
  stats_summary = "Summary statistics of the differential expression analysis"
)

# Save objects structure
saveRDS(objects_structure, "src/tmp_data/objects_structure.rds")

# Print analysis summary
cat("\nAnalysis Summary:\n")
cat("Total cytokines analyzed:", stats_summary$total_cytokines, "\n")
cat("Cytokines with significant genotype effect:", stats_summary$genotype_significant, "\n")
cat("Cytokines with significant treatment effect:", stats_summary$treatment_significant, "\n")
cat("Cytokines with significant interaction:", stats_summary$interaction_significant, "\n")
cat("Cytokines with any significant effect:", stats_summary$any_significant, "\n")
