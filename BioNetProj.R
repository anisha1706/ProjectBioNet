# Clearing the environment
rm(list = ls())

# Installing necessary packages if not already installed
install.packages(c("Seurat", "dplyr", "ggplot2", "tibble", "readxl"))

# Installing BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Installing necessary Bioconductor packages
BiocManager::install("nichenetr")
BiocManager::install("org.Mm.eg.db")

# Loading the packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(nichenetr)
library(readxl)
library(tidyverse)
library(cli)
library(clusterProfiler)
library(org.Mm.eg.db)

# Updating all packages
update.packages(ask = FALSE)

# Reading data from Excel files
data_CD8T <- read_excel("https://raw.githubusercontent.com/anisha1706/ProjectBioNet/main/GSE231302_tumor_CD8T_cells_metadata_norm.xlsx")
data_cDCs <- read_excel("https://raw.githubusercontent.com/anisha1706/ProjectBioNet/main/GSM6019670_BD_tumor_dc_filtered_norm.xlsx")

# Formatting the CD8+ T cells dataset
colnames(data_CD8T) <- as.character(data_CD8T[1, ])
data_CD8T <- data_CD8T[-1, ]

# Extracting metadata
cd8t_metadata <- data_CD8T[, 1:35]
cdc_metadata <- data_cDCs[, 1:8]

# Transposing the gene expression data to get genes as rows and cells as columns
cd8t_matrix <- as.matrix(t(data_CD8T[, 36:ncol(data_CD8T)]))
rownames(cd8t_matrix) <- colnames(data_CD8T)[36:ncol(data_CD8T)]

cdc_matrix <- as.matrix(t(data_cDCs[, 9:ncol(data_cDCs)]))
rownames(cdc_matrix) <- colnames(data_cDCs)[9:ncol(data_cDCs)]

# Creating Seurat objects
cd8t_seurat <- CreateSeuratObject(counts = cd8t_matrix, project = "CD8T_Cells")
cdc_seurat <- CreateSeuratObject(counts = cdc_matrix, project = "CDC_Cells")

# Adding cell metadata
cd8t_seurat <- AddMetaData(cd8t_seurat, metadata = cd8t_metadata)
cdc_seurat <- AddMetaData(cdc_seurat, metadata = cdc_metadata)

# Setting the global timeout option
options(timeout = 300)

# Defining the organism
organism <- "mouse"

# Function for reading RDS file with error handling
read_rds_url <- function(url) {
  tryCatch({
    readRDS(url(url, "rb"))
  }, error = function(e) {
    message("Error reading RDS file from URL: ", url)
    message("Error: ", e$message)
    NULL
  })
}


# Loading the appropriate networks based on the organism
if (organism == "human") {
  lr_network <- read_rds_url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds")
  ligand_target_matrix <- read_rds_url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds")
  weighted_networks <- read_rds_url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds")
} else if (organism == "mouse") {
  lr_network <- read_rds_url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds")
  ligand_target_matrix <- read_rds_url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds")
  weighted_networks <- read_rds_url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds")
}

# Checking if the files were loaded successfully
if (is.null(lr_network) | is.null(ligand_target_matrix) | is.null(weighted_networks)) {
  stop("Failed to load one or more RDS files. Please check the URLs and try again.")
} else {
  message("Successfully loaded all RDS files.")
}

# Processing loaded data
lr_network <- lr_network %>% distinct(from, to)
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from", "to"))

# Setting cluster identities in Seurat objects
cdc_seurat$seurat_clusters <- as.factor(cdc_seurat$seurat_clusters)  
levels(cdc_seurat$seurat_clusters) <- as.character(0:4)
Idents(cdc_seurat) <- cdc_seurat@meta.data$seurat_clusters

cd8t_seurat$cluster <- as.factor(cd8t_seurat$cluster)  
levels(cd8t_seurat$cluster) <- as.character(1:8)
Idents(cd8t_seurat) <- cd8t_seurat@meta.data$cluster

# Defining sender and receiver cells
sender_cells <- 0:4
receiver_cells <- 1:8

# Getting expressed genes in receiver cells
expressed_genes_receiver <- get_expressed_genes(receiver_cells, cd8t_seurat, pct = 0.05)
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

# Identifying potential ligands
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
list_expressed_genes_sender <- sender_cells %>% unique() %>% lapply(get_expressed_genes, cdc_seurat, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)

# Performing differential expression analysis for conditions
Idents(cd8t_seurat) <- cd8t_seurat$condition
conditions <- c("CD4_Cre_Ptger2dd_Ptger4ff_BRAFV600E", "GzmB_Cre_Ptger2dd_Ptger4ff_BRAFV600E", "Ptger2dd_Ptger4ff_BRAFV600E")

degs_conditions_list <- list()
# Normalizing the data (if not already normalized)
cd8t_seurat <- NormalizeData(cd8t_seurat, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    condition_1 <- conditions[i]
    condition_2 <- conditions[j]
    
    # Ensuring the identifiers exist
    if (!(condition_1 %in% unique(Idents(cd8t_seurat)))) stop(paste("Identifier", condition_1, "not found in the Seurat object."))
    if (!(condition_2 %in% unique(Idents(cd8t_seurat)))) stop(paste("Identifier", condition_2, "not found in the Seurat object."))
    
    # Performing differential expression analysis
    degs <- FindMarkers(cd8t_seurat, ident.1 = condition_1, ident.2 = condition_2, min.pct = 0.25)
    
    # Filtering significant DEGs
    significant_degs <- degs %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
    
    # Storing results
    degs_conditions_list[[paste(condition_1, "vs", condition_2)]] <- significant_degs
  }
}

# Combining DEGs from all condition comparisons
gene_set_conditions <- unique(unlist(lapply(degs_conditions_list, rownames)))

# Performing differential expression analysis for clusters
Idents(cd8t_seurat) <- cd8t_seurat$cluster
clusters <- unique(Idents(cd8t_seurat))

degs_clusters_list <- list()
for (i in 1:(length(clusters) - 1)) {
  for (j in (i + 1):length(clusters)) {
    cluster_1 <- clusters[i]
    cluster_2 <- clusters[j]
    
    # Ensuring the identifiers exist
    if (!(cluster_1 %in% unique(Idents(cd8t_seurat)))) stop(paste("Cluster", cluster_1, "not found in the Seurat object."))
    if (!(cluster_2 %in% unique(Idents(cd8t_seurat)))) stop(paste("Cluster", cluster_2, "not found in the Seurat object."))
    
    # Performing differential expression analysis
    degs <- FindMarkers(cd8t_seurat, ident.1 = cluster_1, ident.2 = cluster_2, min.pct = 0.25)
    
    # Filtering significant DEGs
    significant_degs <- degs %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
    
    # Storing results
    degs_clusters_list[[paste(cluster_1, "vs", cluster_2)]] <- significant_degs
  }
}

# Combining DEGs from all cluster comparisons
gene_set_clusters <- unique(unlist(lapply(degs_clusters_list, rownames)))
gene_set_of_interest <- unique(c(gene_set_conditions, gene_set_clusters))

# Predicting ligand activities
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# Performing ligand activity prediction

ligand_activities <- predict_ligand_activities(
  geneset = gene_set_of_interest,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Arranging ligand activities and adding rank
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

# Checking the result
head(ligand_activities)

# Plotting histogram of ligand activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) +
  geom_histogram(color = 'black', fill = 'darkorange') +
  geom_vline(aes(xintercept = min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color = 'red', linetype = 'dashed', size = 1) +
  labs(x = 'Ligand activity (AUPR)', y = '#ligands') +
       theme_classic()
       
       print(p_hist_lig_activity)

# Selecting best upstream ligands
       
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
print(best_upstream_ligands)



# Visualizing ligand activities

vis_ligand_aupr <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames('test_ligand') %>%
  dplyr::select(aupr_corrected) %>%
  as.matrix()

# Checking the resulting matrix
print(vis_ligand_aupr)

# Plotting heatmap of ligand activities

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     'Prioritized ligands', 'Ligand activity',
                                     legend_title = 'AUPR', color = 'darkorange') +
  theme(axis.text.x.top = element_blank(),  
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm")) +  
  guides(fill = guide_colorbar(label.position = "bottom", 
                               title.position = "top", 
                               title.hjust = 0.5, 
                               label.hjust = 0.5, 
                               barwidth = 10, 
                               barheight = 1))


print(p_ligand_aupr)

# Getting active ligand-target links

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = gene_set_of_interest,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

# Preparing ligand-target visualization

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.4
)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

# Plotting heatmap of ligand-target links

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, 'Prioritized ligands', 'Predicted target genes',
                                       color = 'purple', legend_title = 'Regulatory potential') +
  scale_fill_gradient2(low = 'whitesmoke', high = 'purple')

print(p_ligand_target)

# Getting ligand-receptor links

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig
)

# Preparing ligand-receptor visualization

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = 'both'
)

# Plotting heatmap of ligand-receptor links

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network),
                                         y_name = 'Ligands', x_name = 'Receptors',
                                         color = 'mediumvioletred', legend_title = 'Prior interaction potential')

print(p_ligand_receptor)



# Normalizing the data for CDC Seurat object

cdc_seurat <- NormalizeData(cdc_seurat, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

# Accessing the RNA assay

rna_assay <- cdc_seurat[['RNA']]

# Checking if the ‘data’ layer exists

if ('data' %in% names(rna_assay@layers)) {
  print('The data layer is present in the RNA assay.')
} else {
  
# Manually creating the ‘data’ layer from the normalized counts
  
  rna_assay@layers$data <- rna_assay@layers$counts
  print('The data layer was not found and has been created from the counts layer.')
}

# Updating the Seurat object with the modified assay

cdc_seurat[['RNA']] <- rna_assay

# Ensuring sender cells are clusters that are in seurat_clusters

if (!'seurat_clusters' %in% colnames(cdc_seurat@meta.data)) {
  stop('The seurat_clusters column is not found in the Seurat object’s metadata.')
} else {
  print('The seurat_clusters column is found in the Seurat object’s metadata.')
}

# Checking if the sender cells are correctly referenced

valid_sender_cells <- sender_cells[sender_cells %in% unique(cdc_seurat$seurat_clusters)]
if (length(valid_sender_cells) == 0) {
  stop('None of the sender cells are found in the seurat_clusters column.')
} else {
  print('Valid sender cells found in the seurat_clusters column.')
}

# Subsetting the Seurat object for the sender cells

cdc_seurat_subset <- subset(cdc_seurat, seurat_clusters %in% valid_sender_cells)

# Defining the best upstream ligands 

best_upstream_ligands <- c('Ptprc', 'Cd28', 'Btla', 'Tnf', 'Il12a', 'Cd86', 'Il15', 'Cd80', 'Il3', 'Slamf7', 'Cd274', 'Dlk1', 'Cd48', 'Cd22', 'Psen1')

# Generating the DotPlot

p_dotplot <- DotPlot(cdc_seurat_subset, features = rev(best_upstream_ligands), cols = 'RdYlBu') +
  coord_flip() +
  scale_y_discrete(position = 'right')

print(p_dotplot)

# Setting identities to the ‘condition’ column

Idents(cd8t_seurat) <- cd8t_seurat$condition

# List of conditions

conditions <- c('CD4_Cre_Ptger2dd_Ptger4ff_BRAFV600E', 'GzmB_Cre_Ptger2dd_Ptger4ff_BRAFV600E', 'Ptger2dd_Ptger4ff_BRAFV600E')

# Performing differential expression analysis for each pair of conditions

degs_conditions_list <- list()
for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    condition_1 <- conditions[i]
    condition_2 <- conditions[j]
    
    # Performing differential expression analysis
    degs <- FindMarkers(cd8t_seurat, ident.1 = condition_1, ident.2 = condition_2, min.pct = 0.25)
    
    # Converting row names to a column and filtering significant DEGs
    significant_degs <- rownames_to_column(degs, var = "gene") %>%
      filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
    
    # Storing results
    degs_conditions_list[[paste(condition_1, "vs", condition_2)]] <- significant_degs
  }
}

# Setting identities to the ‘seurat_clusters’ column

Idents(cd8t_seurat) <- cd8t_seurat$cluster

# Getting all unique clusters

clusters <- unique(cd8t_seurat$cluster)

# Initializing a list to store DEG results for clusters

degs_clusters_list <- list()
for (i in 1:(length(clusters) - 1)) {
  for (j in (i + 1):length(clusters)) {
    cluster_1 <- clusters[i]
    cluster_2 <- clusters[j]
    # Performing differential expression analysis
    degs <- FindMarkers(cd8t_seurat, ident.1 = cluster_1, ident.2 = cluster_2, min.pct = 0.25)
    
    # Converting row names to a column and filtering significant DEGs
    significant_degs <- rownames_to_column(degs, var = "gene") %>%
      filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
    
    # Storing results
    degs_clusters_list[[paste(cluster_1, "vs", cluster_2)]] <- significant_degs
  }
}

# Function to perform enrichment analysis

perform_enrichment_analysis <- function(degs_df) {
  
  #Converting gene symbols to Entrez IDs
  gene_conversion <- bitr(degs_df$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
  
  # Merging the conversion results back with significant DEGs
  
  degs_df <- merge(degs_df, gene_conversion, by.x = 'gene', by.y = 'SYMBOL')
  
  # Performing GO enrichment analysis using Entrez IDs
  
  ego <- enrichGO(
    gene = degs_df$ENTREZID,
    keyType = 'ENTREZID',
    OrgDb = org.Mm.eg.db,
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  return(ego)
}



# Initializing list to store enrichment results for conditions

enrichment_results_conditions <- list()

# Performing enrichment analysis for each condition comparison

for (comparison in names(degs_conditions_list)) {
  degs <- degs_conditions_list[[comparison]]
  if (nrow(degs) > 0) {  # Ensuring there are DEGs to analyze
    ego <- perform_enrichment_analysis(degs)
    enrichment_results_conditions[[comparison]] <- ego
  } else {
    enrichment_results_conditions[[comparison]] <- NULL
  }
}

# Initializing list to store enrichment results for clusters

enrichment_results_clusters <- list()

# Performing enrichment analysis for each cluster comparison

for (comparison in names(degs_clusters_list)) {
  degs <- degs_clusters_list[[comparison]]
  if (nrow(degs) > 0) {  # Ensuring there are DEGs to analyze
    ego <- perform_enrichment_analysis(degs)
    enrichment_results_clusters[[comparison]] <- ego
  } else {
    enrichment_results_clusters[[comparison]] <- NULL
  }
}

# Visualizing the results for conditions 

if (!is.null(enrichment_results_conditions[[1]])) {
  dotplot(enrichment_results_conditions[[1]]) + ggtitle('GO Enrichment Analysis for Conditions')
}

# Visualizing the results for clusters

if (!is.null(enrichment_results_clusters[[1]])) {
  dotplot(enrichment_results_clusters[[1]]) + ggtitle('GO Enrichment Analysis for Clusters')
}

# Creating line plot for ligand activities 

make_line_plot(ligand_activities = ligand_activities,
               potential_ligands = potential_ligands_focused) +
  theme(plot.title = element_text(size = 11, hjust = 0.1, margin = margin(0, 0, -5, 0)))
