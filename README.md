# ProjectBioNet

## Overview
This project investigates the co-localization and communication pathways between specific CD8+ T cell and cDC subpopulations in melanoma mouse models. The study aims to understand how cDCs recruit CD8+ T cell progenitors, which then produce effector cells to infiltrate and kill tumors.

## Requirements
To run the code, you need the following software and packages:

- **R** (version 4.0 or higher)
- **RStudio** (recommended for ease of use)
- R packages:
  - `Seurat`
  - `dplyr`
  - `ggplot2`
  - `tibble`
  - `readxl`
  - `tidyverse`
  - `cli`
  - `clusterProfiler`
  - `nichenetr`
  - `org.Mm.eg.db`
  - `BiocManager`

## Installation Instructions
1. **Install R and RStudio**:
   - Download and install R from [CRAN](https://cran.r-project.org/).
   - Download and install RStudio from [RStudio](https://rstudio.com/products/rstudio/download/).

2. **Install Required R Packages**:
   - Open RStudio and run the following commands to install the necessary packages:
     ```r
     install.packages(c("Seurat", "dplyr", "ggplot2", "tibble", "readxl", "tidyverse", "cli"))
     if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
     }
     BiocManager::install("nichenetr")
     BiocManager::install("org.Mm.eg.db")
     BiocManager::install("clusterProfiler")
     ```

## How to Run the Code
1. **Clone the Repository**:
   - Open a terminal and run the following command to clone the repository:
     ```bash
     git clone https://github.com/anisha1706/ProjectBioNet.git
     ```

2. **Open the R Script**:
   - Open RStudio and open the `BioNetProj.R` file from the cloned repository.

3. **Run the Script**:
   - In RStudio, source the script by running:
     ```r
     source("BioNetProj.R")
     ```
   - Ensure that all necessary data files are accessible via URLs provided in the script.

## Input and Output
### Required Input
- **Data Files**: The following data files are required and should be placed accessible via URLs:
  - `GSE231302_tumor_CD8T_cells_metadata_norm.xlsx`
  - `GSM6019670_BD_tumor_dc_filtered_norm.xlsx`

### Produced Output
- **Plots**:
  - Histogram of ligand activities (`p_hist_lig_activity`)
  - Heatmap of ligand activities (`p_ligand_aupr`)
  - Heatmap of ligand-target links (`p_ligand_target`)
  - Heatmap of ligand-receptor links (`p_ligand_receptor`)
  - DotPlot of best upstream ligands (`p_dotplot`)
  - GO Enrichment Analysis for Conditions (`enrichment_results_conditions`)
  - GO Enrichment Analysis for Clusters (`enrichment_results_clusters`)

- **Data Frames**:
  - `ligand_activities`: Data frame containing predicted ligand activities.
  - `degs_conditions_list`: List of data frames containing differentially expressed genes (DEGs) for each condition comparison.
  - `degs_clusters_list`: List of data frames containing DEGs for each cluster comparison.
  - `enrichment_results_conditions`: List of enrichment analysis results for each condition comparison.
  - `enrichment_results_clusters`: List of enrichment analysis results for each cluster comparison.

## Example Output
Example plots can be found in the `figures/` directory. These plots provide visualizations of the ligand activities, ligand-target interactions, and GO enrichment results.

## Session Info
The session information used for this analysis can be found in `sessionInfo.txt`, which includes details about the R environment and package versions used.

## Contact
For any questions, please contact [Anisha Bhandare](mailto:anisha.bhandare@gmail.com).
