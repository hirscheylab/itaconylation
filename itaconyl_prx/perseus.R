library(tidyverse)
library(ggrepel)  # For better label placement

itaconyl_perseus_raw <- readxl::read_excel("data/Itaconylpeptides_2024-11-04_PAG.xlsx", sheet = "Itaconyl_Perseus Output") |>
  relocate(`T: Gene Symbol`, .before = 1)

glimpse(itaconyl_perseus_raw)

# Step 1: Clean column names and pivot to long format
itaconyl_perseus <- itaconyl_perseus_raw %>%
  # Rename the gene symbol column for clarity
  rename(gene_symbol = `T: Gene Symbol`) %>%
  # Select abundance columns and identifier columns
  pivot_longer(
    cols = starts_with("Abundances Normalized"),
    names_to = "sample_info",
    values_to = "abundance"
  ) %>%
  # Extract sample information from column names
  mutate(
    # Remove the "Abundances Normalized " prefix
    sample_info = str_remove(sample_info, "Abundances Normalized "),
    # Extract sample ID (F1, F2, etc.)
    sample_id = str_extract(sample_info, "F\\d+"),
    # Extract sample type (Sample or Control)
    sample_type = str_extract(sample_info, "Sample|Control"),
    # Extract genotype (ko/wt)
    genotype = case_when(
      str_detect(sample_info, "ko") ~ "ko",
      str_detect(sample_info, "wt") ~ "wt",
      str_detect(sample_info, "Pool") ~ "pool",
      TRUE ~ NA_character_
    ),
    # Extract treatment information
    treatment = case_when(
      str_detect(sample_info, "untreated") ~ "untreated",
      str_detect(sample_info, "6h LPS") ~ "6h_LPS",
      str_detect(sample_info, "24h LPS") ~ "24h_LPS",
      str_detect(sample_info, "48h LPS") ~ "48h_LPS",
      str_detect(sample_info, "Pool") ~ "pool",
      TRUE ~ NA_character_
    ),
    # Convert abundance from character to numeric
    abundance = as.numeric(abundance)
  ) %>%
  # Select relevant columns
  select(
    sample_id, sample_type, genotype, treatment,
    gene_symbol, abundance,
    sequence = `T: Sequence`,
    modifications = `T: Modifications`,
    protein_accessions = `T: Protein Accessions`,
    master_protein_accessions = `T: Master Protein Accessions`,
    confidence_byonic = `C: Confidence by Search Engine A6 PMI-Byonic`,
    pep = `N: PEP`,
    q_value = `N: q-Value`,
    psm_count = `N: Number of PSMs`
  )

# Add a new column to extract modification sites from the 'modifications' column
itaconyl_perseus <-
  itaconyl_perseus %>%
  mutate(
    # Extract the modification site (letter and number between "Itaconyl [" and "(100)")
    mod_site = str_extract(modifications, "(?<=Itaconyl \\[)[A-Z][0-9]+(?=\\([0-9.]+\\))"),
    mod_site = ifelse(is.na(mod_site), "Unlocalized", mod_site)
  )

# Preview the tidy data
head(itaconyl_perseus)

# If you want to save the tidy dataset
# write_csv(tidy_data, "tidy_itaconyl_perseus.csv")

# Calculate mean abundances with the same approach as before
gene_means <- itaconyl_perseus %>%
  filter(genotype %in% c("ko", "wt"), !is.na(gene_symbol)) %>%
  group_by(gene_symbol, genotype, treatment) %>%
  summarize(
    mean_abundance = mean(abundance, na.rm = TRUE),
    confidence = first(confidence_byonic),
    sample_count = n(),
    .groups = "drop"
  )

# Calculate fold changes as before
gene_comparisons <- gene_means %>%
  pivot_wider(
    id_cols = c(gene_symbol, treatment, confidence),
    names_from = genotype,
    values_from = c(mean_abundance, sample_count)
  ) %>%
  mutate(
    log2FC = mean_abundance_ko - mean_abundance_wt,  # Data is already log-transformed
    treatment = factor(treatment,
                       levels = c("untreated", "6h_LPS", "24h_LPS", "48h_LPS")),
    # Assign treatment size factors (1 for untreated, 2 for 6h, etc.)
    treatment_size = case_when(
      treatment == "untreated" ~ 1,
      treatment == "6h_LPS" ~ 2,
      treatment == "24h_LPS" ~ 3,
      treatment == "48h_LPS" ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(mean_abundance_ko) & !is.na(mean_abundance_wt))

# Calculate statistics with p-values
volcano_data <- gene_comparisons %>%
  rowwise() %>%
  mutate(
    # Call the safe p-value function (assuming it's defined as in your original code)
    p_result = list(safe_p_value(gene_symbol, treatment, itaconyl_perseus)),
    p_value = p_result$p_value,
    neg_log10_p = -log10(p_value),

    # Define significance categories matching the example image
    significance = case_when(
      p_value < 0.05 & log2FC > 1.2 ~ "Up",
      # p_value < 0.05 & log2FC > 1.2 ~ "Highly up",
      p_value >= 0.05 & log2FC > 1 ~ "Up",
      p_value < 0.05 & log2FC < -1.2 ~ "Highly down",
      p_value >= 0.05 & log2FC < -1 ~ "Down",
      TRUE ~ "Not significant"
    ),

    # Create labels for significant points (similar to example)
    label = ifelse(p_value <= 0.05 | abs(log2FC) > 1.2, gene_symbol, "")
  ) %>%
  ungroup() %>%
  select(-p_result)

# Define a custom color palette to match the example
sig_colors <- c(
  "Highly up" = "red",
  "Up" = "red",
  "Not significant" = "gray",
  "Down" = "blue",
  "Highly down" = "navy"
)

# Define the order for the legend (from highly up to highly down)
legend_order <- c("Highly up", "Up", "Not significant", "Down", "Highly down")

# Create the updated volcano plot
volcano_plot_updated <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_p)) +
  # Add vertical and horizontal threshold lines
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dashed", color = "gray", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", alpha = 0.7) +

  # Add points with coloring by significance and size by treatment
  geom_point(aes(color = significance, size = treatment_size), alpha = 0.8) +

  # Add labels for significant points
  # geom_text_repel(
  #   aes(label = label),
  #   max.overlaps = 15,
  #   box.padding = 0.5,
  #   min.segment.length = 0,
  #   force = 3
  # ) +

  # Custom colors and sizes
  scale_color_manual(values = sig_colors, breaks = legend_order) +
  scale_size_continuous(range = c(1.5, 4),
                        breaks = c(1, 2, 3, 4),
                        labels = c("untreated", "6h LPS", "24h LPS", "48h LPS")) +

  # Adjust theme to match the example
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "lightgray", size = 0.2),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +

  # Add appropriate labels
  labs(
    #title = "Differential Expression Analysis (KO vs WT)",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance",
    size = "Treatment"
  ) +

  # Adjust plot dimensions
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 3))

# Display the plot
print(volcano_plot_updated)

# If you want to save the plot
ggsave("itaconyl_volcano_plot.png", volcano_plot_updated, width = 7, height = 7, dpi = 300, bg = "white")

# Count the total number of unique modification sites per gene
gene_unique_sites <- itaconyl_perseus %>%
  # Select relevant columns
  select(gene_symbol, mod_site) %>%
  # Remove duplicates
  distinct() %>%
  # Group by gene
  group_by(gene_symbol) %>%
  # Count number of unique modification sites
  summarize(num_unique_sites = n()) %>%
  # Arrange by number of sites in descending order
  arrange(desc(num_unique_sites)) %>%
  # Filter out NA or empty gene symbols
  filter(!is.na(gene_symbol), gene_symbol != "")

# View the genes with the most unique modification sites
head(gene_unique_sites, 20)

# Load required packages
library(tidyverse)
library(ggrepel)

# Count proteins by number of modification sites
site_distribution <- gene_unique_sites %>%
  # Create a factor for the number of sites
  mutate(
    site_count_category = case_when(
      num_unique_sites >= 3 ~ "3+",
      TRUE ~ as.character(num_unique_sites)
    ),
    # Convert to factor with specific order
    site_count_category = factor(site_count_category, levels = c("1", "2", "3+"))
  ) %>%
  # Count proteins in each category
  count(site_count_category) %>%
  # Calculate percentages
  mutate(percentage = n / sum(n) * 100)

# Define improved color scheme
site_colors <- c(
  "1" = "lightgray",
  "2" = "pink",
  "3+" = "red"
)

# Create improved pie chart with better legend positioning
improved_pie <- ggplot(site_distribution, aes(x = 0, y = n, fill = site_count_category)) +
  # Use a pie chart
  geom_col(width = 1, position = "stack") +
  coord_polar("y", start = 0) +

  # Use standard geom_text for the large slice
  geom_text(
    data = site_distribution %>% filter(percentage > 90),
    aes(y = sum(site_distribution$n)/2, label = sprintf("%.1f%%", percentage)),
    size = 5
  ) +

  # Better positioning for small slices
  geom_text(
    data = site_distribution %>% filter(percentage < 10 & site_count_category == "2"),
    aes(label = sprintf("%.1f%%", percentage)),
    x = 0, y = sum(site_distribution$n) * 0.98,  # Position at top
    hjust = -0.7, vjust = -1,
    size = 5
  ) +

  # Even further positioning for the smallest slice
  geom_text(
    data = site_distribution %>% filter(percentage < 10 & site_count_category == "3+"),
    aes(label = sprintf("%.1f%%", percentage)),
    x = 0, y = sum(site_distribution$n) * 0.99,  # Position even higher
    hjust = 0.75, vjust = -2.5,
    size = 5
  ) +

  # Use the improved color scheme
  scale_fill_manual(values = site_colors) +

  # Add a clean white background with no border
  theme_void() +
  theme(
    # Place legend in a better position
    legend.position = c(0.95, 0.5),
    legend.justification = c(0, 0.5),
    legend.margin = margin(0, 0, 0, 10),
    legend.box.margin = margin(0, 0, 0, 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.3, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    # Reduce plot margin but ensure there's space for the legend
    plot.margin = margin(5, 40, 5, 5)
  ) +

  # Add appropriate labels
  labs(
    fill = "Number\nof Sites"
  ) +

  # Expand the plot size using the aspect ratio
  theme(aspect.ratio = 1)

# Display the improved pie chart
print(improved_pie)

# Save with better dimensions for the legend placement
ggsave("itaconyl_pie.png", improved_pie, width = 7, height = 7, dpi = 300, bg = "white")
