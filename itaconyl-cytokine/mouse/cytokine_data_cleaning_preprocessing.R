
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

# Read the data
cytokine_data <- tryCatch({
  readxl::read_excel("src/data/20180921_mouse_cytokines_whole_body_R.xlsx")
}, error = function(e) {
  stop("Error reading Excel file: ", e$message)
})

# Inspect column names and clean them
colnames(cytokine_data) <- gsub("[`:]", "", colnames(cytokine_data))

# Clean the data - identify metadata and cytokine columns
meta_cols <- c("Samples", "mouse #", "mouse ID", "genotype", "time")
cytokine_cols <- setdiff(colnames(cytokine_data), meta_cols)

# Rename metadata columns to more R-friendly names
cytokine_data <- cytokine_data %>%
  dplyr::rename_with(~ "sample_id", .cols = "Samples") %>%
  dplyr::rename_with(~ "mouse_number", .cols = "mouse #") %>%
  dplyr::rename_with(~ "mouse_id", .cols = "mouse ID")

# Convert time to ordered factor and genotype to factor
cytokine_data <- cytokine_data %>%
  dplyr::mutate(
    time = factor(time, levels = c("0 h", "2 h", "7 h"), ordered = TRUE),
    genotype = factor(genotype, levels = c("wt", "ko"))
  )

# Save clean data
tryCatch({
  openxlsx::write.xlsx(cytokine_data, file = "src/tmp_data_clean/cytokine_data_clean.xlsx", 
                      rowNames = FALSE, overwrite = TRUE)
}, error = function(e) {
  stop("Error writing Excel file: ", e$message)
})
