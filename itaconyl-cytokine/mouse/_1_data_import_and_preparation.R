
# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Create a function for safely loading data with error handling
load_data_safely <- function(file_path, error_message = NULL) {
  if (is.null(error_message)) {
    error_message <- paste("Error loading data from", file_path)
  }
  
  tryCatch({
    if (file.exists(file_path)) {
      readr::read_rds(file_path)
    } else {
      stop(paste("File not found:", file_path))
    }
  }, error = function(e) {
    message(error_message)
    message(e$message)
    NULL
  })
}

# Load data from RDS files
cytokine_data_wide <- load_data_safely("src/data/cytokine_data_wide.rds", 
                                       "Error loading wide-format cytokine data")
cytokine_data_long <- load_data_safely("src/data/cytokine_data_long.rds", 
                                       "Error loading long-format cytokine data")
cytokine_names <- load_data_safely("src/data/cytokine_names.rds", 
                                   "Error loading cytokine names")

# Display the first few rows of the data
cat("Structure of wide-format cytokine data:\n")
print(str(cytokine_data_wide))
cat("\nFirst 10 rows of wide-format cytokine data:\n")
print(knitr::kable(head(cytokine_data_wide, 10)))

cat("\nStructure of long-format cytokine data:\n")
print(str(cytokine_data_long))
cat("\nFirst 10 rows of long-format cytokine data:\n")
print(knitr::kable(head(cytokine_data_long, 10)))

cat("\nCytokine names:\n")
print(cytokine_names)

# Check data completeness
# Summarize missing values in wide format
missing_values_wide <- cytokine_data_wide %>%
  summarise(across(where(is.numeric), ~sum(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to = "Column", values_to = "Missing_Count")

# Count zero values in cytokine measurements
zero_values_wide <- cytokine_data_wide %>%
  summarise(across(all_of(cytokine_names), ~sum(. == 0, na.rm = TRUE))) %>%
  tidyr::pivot_longer(everything(), names_to = "Cytokine", values_to = "Zero_Count")

# Display missing and zero values summaries
cat("\nMissing values summary:\n")
print(knitr::kable(missing_values_wide))

cat("\nZero values summary by cytokine:\n")
print(knitr::kable(zero_values_wide))

# Apply log transformation to handle wide concentration ranges
# Create a copy with log-transformed values (log(x+1) to handle zeros)
cytokine_data_wide_log <- cytokine_data_wide %>%
  mutate(across(all_of(cytokine_names), ~log1p(.)))

# Save the log-transformed data
saveRDS(cytokine_data_wide_log, "src/tmp_data/cytokine_data_wide_log.rds")

# Create long version of log-transformed data for easier plotting
cytokine_data_long_log <- cytokine_data_wide_log %>%
  tidyr::pivot_longer(
    cols = all_of(cytokine_names),
    names_to = "cytokine",
    values_to = "log_concentration"
  )

# Save the long-format log-transformed data
saveRDS(cytokine_data_long_log, "src/tmp_data/cytokine_data_long_log.rds")

# Create visualization of distribution before and after log transformation
# Sample a few cytokines for visualization
selected_cytokines <- sample(cytokine_names, min(5, length(cytokine_names)))

# Create distribution plots for original values
distribution_plot_original <- cytokine_data_long %>%
  filter(cytokine %in% selected_cytokines) %>%
  ggplot(aes(x = concentration, fill = cytokine)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  facet_wrap(~cytokine, scales = "free") +
  scale_fill_manual(values = c("#4D6D8E", "#7AA661", "#807E7D", "#8e6e4d", "#7b4d8e")) +
  theme_bw() +
  labs(title = "Distribution of original cytokine concentrations",
       x = "Concentration (pg/ml)",
       y = "Count") +
  theme(legend.position = "none")

# Create distribution plots for log-transformed values
distribution_plot_log <- cytokine_data_long_log %>%
  filter(cytokine %in% selected_cytokines) %>%
  ggplot(aes(x = log_concentration, fill = cytokine)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  facet_wrap(~cytokine, scales = "free") +
  scale_fill_manual(values = c("#4D6D8E", "#7AA661", "#807E7D", "#8e6e4d", "#7b4d8e")) +
  theme_bw() +
  labs(title = "Distribution of log-transformed cytokine concentrations",
       x = "Log(Concentration + 1)",
       y = "Count") +
  theme(legend.position = "none")

# Save the distribution plots
ggsave("src/tmp_figures/cytokine_distribution_original.png", 
       distribution_plot_original, width = 10, height = 6)
ggsave("src/tmp_figures/cytokine_distribution_log.png", 
       distribution_plot_log, width = 10, height = 6)

# Create derived variables for fold changes from baseline (time 0h)
# First, calculate average baseline values for each genotype and cytokine
baseline_values <- cytokine_data_wide %>%
  filter(time == "0 h") %>%
  group_by(genotype) %>%
  summarise(across(all_of(cytokine_names), ~mean(., na.rm = TRUE))) %>%
  tidyr::pivot_longer(
    cols = all_of(cytokine_names),
    names_to = "cytokine",
    values_to = "baseline_mean"
  )

# Create fold change data
fold_change_data <- cytokine_data_long %>%
  left_join(baseline_values, by = c("genotype", "cytokine")) %>%
  mutate(
    # Add small value to avoid division by zero
    fold_change = concentration / (baseline_mean + 0.01),
    log2_fold_change = log2((concentration + 0.01) / (baseline_mean + 0.01))
  )

# Save the fold change data
saveRDS(fold_change_data, "src/tmp_data/cytokine_fold_change_data.rds")

# Create wide format of fold change data
fold_change_data_wide <- fold_change_data %>%
  select(sample_id, mouse_number, mouse_id, genotype, time, cytokine, fold_change, log2_fold_change) %>%
  tidyr::pivot_wider(
    id_cols = c(sample_id, mouse_number, mouse_id, genotype, time),
    names_from = cytokine,
    values_from = c(fold_change, log2_fold_change),
    names_sep = "_"
  )

# Save the wide-format fold change data
saveRDS(fold_change_data_wide, "src/tmp_data/cytokine_fold_change_data_wide.rds")

# Visualize fold changes for selected cytokines
fold_change_plot <- fold_change_data %>%
  filter(cytokine %in% selected_cytokines & time != "0 h") %>%
  ggplot(aes(x = time, y = log2_fold_change, color = genotype, group = interaction(genotype, cytokine))) +
  geom_boxplot() +
  facet_wrap(~cytokine, scales = "free_y") +
  scale_color_manual(values = c("#4D6D8E", "#7AA661")) +
  theme_bw() +
  labs(title = "Log2 fold changes from baseline by genotype and time",
       x = "Time",
       y = "Log2 Fold Change")

# Save the fold change plot
ggsave("src/tmp_figures/cytokine_fold_change_plot.png", 
       fold_change_plot, width = 10, height = 6)

# Print summary of created objects
cat("\nCreated and saved objects:\n")
cat("1. cytokine_data_wide_log.rds - Log-transformed wide-format data\n")
cat("2. cytokine_data_long_log.rds - Log-transformed long-format data\n") 
cat("3. cytokine_fold_change_data.rds - Fold changes from baseline in long format\n")
cat("4. cytokine_fold_change_data_wide.rds - Fold changes from baseline in wide format\n")
cat("5. cytokine_distribution_original.png - Histogram of original concentration distributions\n")
cat("6. cytokine_distribution_log.png - Histogram of log-transformed concentration distributions\n")
cat("7. cytokine_fold_change_plot.png - Boxplot of log2 fold changes by genotype and time\n")
