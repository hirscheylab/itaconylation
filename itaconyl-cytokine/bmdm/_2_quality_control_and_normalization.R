
# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stats)

# Load the cytokine data
cytokine_data <- readRDS("src/data/cytokine_data.rds")

# Load QC results for experimental design verification
cytokine_qc_results <- readRDS("src/data/cytokine_qc_results.rds")

# Identify cytokine columns
cytokine_cols <- colnames(cytokine_data)[6:ncol(cytokine_data)]

# 1. Visualize distribution of cytokine measurements

# Create histograms for each cytokine to visualize distribution
for (cytokine in cytokine_cols) {
  p <- ggplot(cytokine_data, aes_string(x = cytokine)) +
    geom_histogram(fill = "#4D6D8E", color = "#000000", bins = 15) +
    theme_bw() +
    labs(title = paste("Distribution of", cytokine),
         x = paste(cytokine, "(pg/ml)"),
         y = "Frequency")
  
  # Save histogram
  file_name <- paste0("src/tmp_figures/histogram_", cytokine, ".png")
  ggsave(file_name, p, width = 8, height = 6, dpi = 300)
}

# Create boxplots of cytokines by treatment and genotype
cytokine_long <- cytokine_data %>%
  pivot_longer(cols = all_of(cytokine_cols),
               names_to = "Cytokine",
               values_to = "Value")

# Create boxplots for each cytokine by genotype and treatment
p_box <- ggplot(cytokine_long, aes(x = LPS_treatment_hours, y = Value, fill = Genotype_clean)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ Cytokine, scales = "free_y") +
  scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#807E7D"),
        strip.text = element_text(color = "white")) +
  labs(title = "Cytokine levels by genotype and LPS treatment",
       x = "LPS treatment (hours)",
       y = "Concentration (pg/ml)",
       fill = "Genotype")

# Save boxplot
ggsave("src/tmp_figures/cytokine_boxplots.png", p_box, width = 12, height = 10, dpi = 300)

# 2. Identify potential outliers using boxplot method
# Calculate outlier boundaries for each cytokine
outlier_summary <- cytokine_long %>%
  group_by(Cytokine) %>%
  summarize(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower_bound = Q1 - 1.5 * IQR,
    Upper_bound = Q3 + 1.5 * IQR,
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    Potential_outliers = sum(Value < Lower_bound | Value > Upper_bound, na.rm = TRUE),
    .groups = "drop"
  )

# Display outlier summary
knitr::kable(outlier_summary, caption = "Potential outliers based on 1.5×IQR rule", digits = 3)

# Flag potential outliers in the dataset
cytokine_with_outliers <- cytokine_long %>%
  left_join(outlier_summary %>% select(Cytokine, Lower_bound, Upper_bound), by = "Cytokine") %>%
  mutate(is_outlier = Value < Lower_bound | Value > Upper_bound)

# Count outliers by cytokine and treatment/genotype
outlier_counts <- cytokine_with_outliers %>%
  filter(is_outlier) %>%
  group_by(Cytokine, Genotype_clean, LPS_treatment_hours) %>%
  summarize(outlier_count = n(), .groups = "drop")

# 3. Log-transform cytokine values to improve normality
# Create a new dataframe with log2-transformed values
# Add a small constant to avoid log(0)
small_constant <- 0.01

cytokine_log2 <- cytokine_data %>%
  mutate(across(all_of(cytokine_cols), ~log2(.x + small_constant)))

# Visualize log-transformed distributions
cytokine_log2_long <- cytokine_log2 %>%
  pivot_longer(cols = all_of(cytokine_cols),
               names_to = "Cytokine",
               values_to = "Log2_Value")

# Create histograms of log2-transformed values
p_log2_hist <- ggplot(cytokine_log2_long, aes(x = Log2_Value)) +
  geom_histogram(fill = "#7AA661", color = "#000000", bins = 15) +
  facet_wrap(~ Cytokine, scales = "free_x") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#807E7D"),
        strip.text = element_text(color = "white")) +
  labs(title = "Distribution of log2-transformed cytokine values",
       x = "log2(Cytokine + 0.01)",
       y = "Frequency")

# Save log2 histogram plot
ggsave("src/tmp_figures/log2_transformed_histograms.png", p_log2_hist, width = 12, height = 10, dpi = 300)

# Create boxplots of log2-transformed values
p_log2_box <- ggplot(cytokine_log2_long, aes(x = LPS_treatment_hours, y = Log2_Value, fill = Genotype_clean)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ Cytokine, scales = "free_y") +
  scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#807E7D"),
        strip.text = element_text(color = "white")) +
  labs(title = "Log2-transformed cytokine levels by genotype and LPS treatment",
       x = "LPS treatment (hours)",
       y = "log2(Concentration + 0.01)",
       fill = "Genotype")

# Save log2 boxplot
ggsave("src/tmp_figures/log2_cytokine_boxplots.png", p_log2_box, width = 12, height = 10, dpi = 300)

# 4. Scale cytokine values (z-score normalization)
# We'll create z-scores for cytokines grouped by treatment to compare genotypes
cytokine_z_scores <- cytokine_data %>%
  group_by(LPS_treatment_hours) %>%
  mutate(across(all_of(cytokine_cols), ~scale(.x), .names = "{.col}_z")) %>%
  ungroup()

# Extract z-score columns
z_score_cols <- paste0(cytokine_cols, "_z")

# Convert to long format for visualization
cytokine_z_long <- cytokine_z_scores %>%
  pivot_longer(cols = all_of(z_score_cols),
               names_to = "Cytokine_z",
               values_to = "Z_Value") %>%
  mutate(Cytokine = gsub("_z$", "", Cytokine_z))

# Create boxplots of z-scores
p_z_box <- ggplot(cytokine_z_long, aes(x = LPS_treatment_hours, y = Z_Value, fill = Genotype_clean)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ Cytokine, scales = "free_y") +
  scale_fill_manual(values = c("#4D6D8E", "#7AA661")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#807E7D"),
        strip.text = element_text(color = "white")) +
  labs(title = "Z-score normalized cytokine levels by genotype and LPS treatment",
       x = "LPS treatment (hours)",
       y = "Z-score",
       fill = "Genotype")

# Save z-score boxplot
ggsave("src/tmp_figures/z_score_cytokine_boxplots.png", p_z_box, width = 12, height = 10, dpi = 300)

# 5. Verify experimental design balance
experimental_design <- cytokine_qc_results$experimental_design

# Display experimental design
knitr::kable(experimental_design, caption = "Experimental design balance")

# Create summary of sample counts by genotype and treatment
sample_counts <- cytokine_data %>%
  group_by(Genotype_clean, LPS_treatment_hours) %>%
  summarize(Sample_count = n(), .groups = "drop") %>%
  spread(key = Genotype_clean, value = Sample_count, fill = 0)

# Save QC output files
# Save log2-transformed data
saveRDS(cytokine_log2, "src/tmp_data/cytokine_log2.rds")

# Save z-score normalized data
saveRDS(cytokine_z_scores, "src/tmp_data/cytokine_z_scores.rds")

# Save outlier information
saveRDS(list(
  outlier_summary = outlier_summary,
  outlier_counts = outlier_counts,
  cytokine_with_outliers = cytokine_with_outliers
), "src/tmp_data/cytokine_outliers.rds")

# Create QC results object
qc_results <- list(
  data_structure = list(
    raw_dimensions = dim(cytokine_data),
    log2_dimensions = dim(cytokine_log2),
    z_score_dimensions = dim(cytokine_z_scores),
    cytokine_columns = cytokine_cols,
    missing_values = colSums(is.na(cytokine_data))
  ),
  sample_distribution = sample_counts,
  experimental_design = cytokine_qc_results$experimental_design,
  outlier_summary = outlier_summary,
  transformation_notes = list(
    log2_transformation = "Applied log2 transformation with small constant (0.01) to handle zeros",
    z_score_normalization = "Applied z-score normalization within each treatment group"
  )
)

# Save QC results
saveRDS(qc_results, "src/tmp_data/qc_results.rds")

# Create QC results structure description
qc_results_structure <- list(
  description = "Quality control results for cytokine dataset",
  data_structure = "Information about dimensions and missing values",
  sample_distribution = "Sample counts by genotype and treatment",
  experimental_design = "Experimental design balance information",
  outlier_summary = "Summary of potential outliers identified using the 1.5×IQR rule",
  transformation_notes = "Notes about data transformations applied"
)

# Save QC results structure
saveRDS(qc_results_structure, "src/tmp_data/qc_results_structure.rds")
