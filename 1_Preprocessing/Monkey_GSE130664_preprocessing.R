library(Seurat)
library(tidyverse)
library(anndata)

# Set working directory to location of samples quantified using kallisto
# bustools (kb). Each sample's kb output should be contained within its own
# directory, named after the sample. E.g. YF2_OO's output should be in ./YF2_OO.
setwd("/scratch1/zwang474/GSE130664_Monkey/Count_Matrices/")

# Load in transcript to gene .txt file from kb ref for later use (gene names)
gene_info <- read.table("/scratch1/zwang474/Reference_Genomes/Monkey_Index_5.0/Kallisto_Index/transcripts_to_genes.txt", sep = "\t", header = FALSE)
names(gene_info) <- c("Transcript_ID", "Gene_ID", "Gene_Name", "Transcript_Name", "Chromosome", "Starting_Pos", "Ending_Pos", "Orientation")

# Concatenate all samples into one Anndata object
GSE130664_All <- ConcatData()

# Build Seurat object from concatenated data
GSE130664_Seurat <- CreateSeuratObject(counts = t(as.matrix(GSE130664_All$X)), meta.data = GSE130664_All$obs)

# Add useful metadata to Seurat Object (individual name, tissue type, age group, and mitochondrial content)
GSE130664_Seurat <- ExtractMetaData(GSE130664_Seurat)
GSE130664_Seurat[["percent_mt"]] <- PercentageFeatureSet(GSE130664_Seurat, features = ComputeMitoCont(GSE130664_Seurat))

# Visualize QC metrics
VlnPlot(GSE130664_Seurat, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(GSE130664_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "tissue_type")
FeatureScatter(GSE130664_Seurat, feature1 = "nCount_RNA", feature2 = "percent_mt", group.by = "tissue_type")

# Retain only cells with detected genes > 400, detected UMIs > 1000, and mitochondrial content < 30%
GSE130664_Filtered <- subset(GSE130664_Seurat, subset = nFeature_RNA > 400 & percent_mt < 30 & nCount_RNA > 1000)

# Normalization
GSE130664_Filtered <- NormalizeData(GSE130664_Filtered)

# Feature selection
GSE130664_Filtered <- FindVariableFeatures(GSE130664_Filtered, selection.method = "vst", nfeatures = 3000)

# Scaling
all.genes <- rownames(GSE130664_Filtered)
GSE130664_Filtered <- ScaleData(GSE130664_Filtered, features = all.genes)

# PCA
GSE130664_Filtered <- RunPCA(GSE130664_Filtered, features = VariableFeatures(object = GSE130664_Filtered))

# Identify significant PCs
GSE130664_Filtered <- JackStraw(GSE130664_Filtered, dims = 30, num.replicate = 100)
GSE130664_Filtered <- ScoreJackStraw(GSE130664_Filtered, dims = 1:30)
JackStrawPlot(GSE130664_Filtered, dims = 1:30)
ElbowPlot(GSE130664_Filtered, ndims = 50)

# Cluster using first 15 PCs (based on ElbowPlot/JackstrawPlot)
GSE130664_Filtered <- FindNeighbors(GSE130664_Filtered, dims = 1:15)
GSE130664_Filtered <- FindClusters(GSE130664_Filtered, resolution = 0.4)

# Run t-SNE
GSE130664_Filtered <- RunUMAP(GSE130664_Filtered, dims = 1:15)

# Relevant t-SNE Plots
DimPlot(GSE130664_Filtered, reduction = "umap")
DimPlot(GSE130664_Filtered, reduction = "umap", group.by = "tissue_type")

# Identify cluster gene markers
GSE130664_Markers <- FindAllMarkers(GSE130664_Filtered, min.pct = 0.1,  only.pos = TRUE)
GSE130664_Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> T10_Markers
T10_Markers <- FindGeneNames(T10_Markers)
GSE130664_Markers <- FindGeneNames(GSE130664_Markers)

# Assign cluster identities
GSE130664_Filtered <- RenameIdents(GSE130664_Filtered,
                                   "0" = "Stromal",
                                   "1" = "Smooth Muscle",
                                   "2" = "Granulosa",
                                   "3" = "Oocyte",
                                   "4" = "Unknown",
                                   "5" = "T Cell",
                                   "6" = "Endothelial",
                                   "7" = "Macrophage",
                                   "8" = "Unknown",
                                   "9" = "Unknown",
                                   "10" = "Monocyte")

# Labelled t-SNE
DimPlot(GSE130664_Filtered, reduction = "umap", label = TRUE, pt.size = 0.5, label.box = TRUE, repel = TRUE)

################################################################################
# Function definitions

# Concatenate .h5ad files for each sample into one Anndata object
ConcatData <- function() {
  # Get a list of all sample folders in the current working directory
  sample_folders <- list.dirs(".", recursive = FALSE, full.names = FALSE)
  
  # List to store Anndata objects (.h5ad) for each sample
  sample_list <- list()
  
  # Iterate through samples, storing their Anndata objects in sample_list
  for (sample in sample_folders){
    counts_path <- file.path(".", sample, "counts_unfiltered", "adata.h5ad")
    if (counts_path == "./Kallisto_Logs/counts_unfiltered/adata.h5ad") next
    temp_adata <- read_h5ad(counts_path)
    sample_list[[sample]] <- temp_adata
  }
  
  # Concatenate all samples into one Anndata object
  Concatenated_Data <- concat(sample_list, join = "outer", label = "cell", index_unique = "-")
  
  return(Concatenated_Data)
}

# Add useful metadata to Seurat Object (individual name, tissue type, and age group)
ExtractMetaData <- function(SeurObj) {
  # Temporary dataframe for metadata extraction/manipulation
  CellInfo <- SeurObj@meta.data
  CellInfo$individual <- ""
  CellInfo$tissue_type <- ""
  CellInfo$age_group <- ""
  
  for (i in 1:nrow(CellInfo)) {
    # Split the string in the "cell" column by underscore
    cell_parts <- strsplit(as.character(CellInfo[i, "cell"]), "_")[[1]]
    
    # Extract individual name
    individual_name <- cell_parts[1]
    
    # Extract tissue type
    tissue_split <- substr(cell_parts[2], 1, 2)
    
    # Extract age group
    name_split <- strsplit(cell_parts[1], "")[[1]]
    first_character <- name_split[1]
    
    # Assign values to the corresponding columns
    CellInfo[i, "individual"] <- individual_name
    CellInfo[i, "tissue_type"] <- tissue_split
    CellInfo[i, "age_group"] <- first_character
  }
  
  # # Adds processed cell names to the metadata.
  # for (i in seq_along(row_names)) {
  #   current_name <- row_names[i]
  #   processed_name <- barcodes$processed_cell[barcodes$barcode_seq == current_name]
  #   meta$processed_cell[i] <- processed_name
  # }
  # 
  # # Adds annotations to the metadata.
  # for (i in 1:nrow(meta)) {
  #   processed_cell <- meta$processed_cell[i]
  #   if (processed_cell %in% annotations$cell) {
  #     cell_annotation <- annotations$cluster[annotations$cell == processed_cell]
  #     meta$annotation[i] <- cell_annotation
  #   }
  #   else {
  #     meta$annotation[i] <- "N/A"
  #   }
  # }
  
  # Add extracted metadata to Seurat Object
  CellInfoTrim <- subset(CellInfo, select = c("individual", "tissue_type", "age_group"))
  SeurObj <- AddMetaData(SeurObj, CellInfoTrim) 
}

# Return a list of mitochondrial genes for computation
ComputeMitoCont <- function(SeurObj) {
  # Mitochondrial genes for Macaca fascicularis
  mito_genes = c('ENSMFAG00000000001.1',
                 'ENSMFAG00000000002.1',
                 'ENSMFAG00000000003.1',
                 'ENSMFAG00000000004.1',
                 'ENSMFAG00000000005.1',
                 'ND1',
                 'ENSMFAG00000000006.1',
                 'ENSMFAG00000000007.1',
                 'ENSMFAG00000000008.1',
                 'ND2',
                 'ENSMFAG00000000009.1',
                 'ENSMFAG00000000010.1',
                 'ENSMFAG00000000011.1',
                 'ENSMFAG00000000012.1',
                 'ENSMFAG00000000013.1',
                 'COX1',
                 'ENSMFAG00000000014.1',
                 'ENSMFAG00000000015.1',
                 'COX2',
                 'ENSMFAG00000000016.1',
                 'ATP8',
                 'ATP6',
                 'COX3',
                 'ENSMFAG00000000017.1',
                 'ND3',
                 'ENSMFAG00000000018.1',
                 'ND4L',
                 'ND4',
                 'ENSMFAG00000000019.1',
                 'ENSMFAG00000000020.1',
                 'ENSMFAG00000000021.1',
                 'ND5',
                 'ND6',
                 'ENSMFAG00000000022.1',
                 'CYTB',
                 'ENSMFAG00000000023.1',
                 'ENSMFAG00000000024.1')
  
  # Add mitochondrial content (%) to Seurat Object
  mito_genes_1 = intersect(mito_genes, rownames(SeurObj))
  
  return(mito_genes_1)
}

# Finds gene names based on gene IDs returned by cluster marker analysis, requires tr2g.txt from kb ref loaded as gene_info (above)
FindGeneNames <- function(Markers) {
  
  # Iterate through each row in Markers
  for (i in 1:nrow(Markers)) {
    # Get the gene ID of the current row
    current_gene <- Markers$gene[i]
    
    # Find a matching row in gene_info for that gene ID and retrieve its gene name
    gene_name <- gene_info[gene_info$Gene_ID == current_gene, "Gene_Name"]
    
    # If the retrieved gene name is not blank
    if (gene_name[1] != "") {
      # Copy the gene name to "Gene Name" in Markers
      Markers$gene[i] <- gene_name[1]
    } 
    # If no gene name is available, remove kb ref's trailing identifiers
    else {
      Markers$gene[i] <- sub("\\.\\d+$", "", current_gene)
    }
  }
  
  return(Markers)
}

