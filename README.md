# SIRT4 Controls Macrophage Function and Wound Healing through Control of Protein Itaconylation in Mice

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15372984.svg)](https://doi.org/10.5281/zenodo.15372984)

This repository contains the code and analysis workflows for the paper "SIRT4 Controls Macrophage Function and Wound Healing through Control of Protein Itaconylation in Mice" by Anderson et al.

## Overview

This study identifies a novel enzymatic activity of the mitochondrial sirtuin SIRT4 as a lysine deitaconylase that regulates macrophage inflammatory responses. Using biochemical and proteomics approaches, we demonstrate that SIRT4 efficiently removes itaconyl modifications from target proteins both *in vitro* and *in vivo*. The code in this repository reproduces all analyses and figures from the paper.

## Repository Structure

- **R Analysis Scripts**: Sequential analysis scripts prefixed with numbers:
  - `_1_data_import_exploration.R`: Initial data loading and exploration
  - `_2_quality_control_processing.R`: Quality control steps
  - `_3_exploratory_data_analysis.R`: Exploratory analysis
  - `_4_differential_analysis.R`: Statistical analysis of differences
  - `_5_pathway_analysis.R`: Biological pathway analysis
  - `_6_visualization_reporting.R`: Figure generation for publication

- **Data Processing Scripts**:
  - `macrophage_data_cleaning.R`: Preprocessing of macrophage data
  - `macrophage_data_import.R`: Data import utilities

- **Project Directories**:
  - `itaconyl_prx/`: Itaconyl proteomics analysis
  - `itaconyl-cytokine/`: Cytokine analysis related to itaconylation
  - `itaconyl-mbx/`: Metabolomics analysis related to itaconylation
  - `data/`: Directory for downloaded data (not included in repository)
  - `raw_data/`: Directory for original unprocessed data (not included)
  - `tmp_data/`: Temporary data files generated during analysis
  - `tmp_figures/`: Temporary figures generated during analysis

## Data Availability

All raw data for this project is hosted on Zenodo under DOI: [10.5281/zenodo.15372984](https://doi.org/10.5281/zenodo.15372984).

The dataset includes:
- Mass spectrometry proteomics data for itaconylated proteins
- Cytokine measurements from SIRT4KO and WT mice
- Metabolomics data from macrophages
- Wound healing measurements

## Getting Started

### Prerequisites
- R version 4.0 or higher
- Required R packages (installed automatically by the scripts):
  - tidyverse
  - DESeq2
  - limma
  - ggplot2
  - pheatmap
  - GAUDI (available at github.com/hirscheylab)

### Setup Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/hirscheylab/itaconylation.git
   cd itaconylation
   ```

2. Create a data directory if it doesn't exist:
   ```bash
   mkdir -p data
   ```

3. Download the dataset from Zenodo:
   - Visit [https://doi.org/10.5281/zenodo.15372984](https://doi.org/10.5281/zenodo.15372984)
   - Download the data files and place them in the `data/` directory

4. Run the analysis scripts in numerical order:
   ```bash
   Rscript _1_data_import_exploration.R
   Rscript _2_quality_control_processing.R
   # Continue with remaining scripts...
   ```

### Analysis Workflow

Each script is designed to be run sequentially, with outputs from earlier scripts serving as inputs to later ones:

1. Data import and initial exploration
2. Quality control and preprocessing
3. Exploratory data analysis
4. Differential expression/abundance analysis
5. Pathway analysis and biological interpretation
6. Visualization and figure generation

## Publication Information

- **Title**: SIRT4 Controls Macrophage Function and Wound Healing through Control of Protein Itaconylation in Mice
- **Authors**: Kristin A. Anderson, Beverly deSouza, Pol Castellano-Escuder, Zhihong Lin, Olga R. Ilkayeva, Michael J. Muehlbauer, Christopher B. Newgard, Paul A. Grimsrud, and Matthew D. Hirschey
- **Publication Date**: May 9, 2025
- **DOI**: [10.5281/zenodo.15372984](https://doi.org/10.5281/zenodo.15372984)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions about the code or analysis, please open an issue in this repository or contact:
- Matthew Hirschey - matthew.hirschey@duke.edu
