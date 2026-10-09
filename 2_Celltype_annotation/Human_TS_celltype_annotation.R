options(stringsAsFactors = FALSE)

library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scran)

library(SingleR)

library(scGate)
library(scSorter)

library(dplyr)
library(tidyverse)

# rm(list = ls())

################################################################################
# 10x ovary Human Tabula Sapiens dataset
# Annotate cell type
################################################################################

################################################################################
# 1. Load dataset
################################################################################

load("/Volumes/OIProject_II/1_R/Human/Human_TS/2026-05-02_10x_ovary_Human_TS_integrated_harmony_post_clustering.RData")

obj <- ovary.tabula.sapiens.harmony

DefaultAssay(obj) <- "SCT"

old_names <- rownames(obj[["SCT"]])

ensembl_clean <- sub("\\..*$", "", old_names)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

new_names <- ifelse(is.na(gene_symbols), old_names, gene_symbols)
new_names <- make.unique(new_names)

# Rename RNA assay feature names
rownames(obj[["SCT"]]@counts) <- new_names
rownames(obj[["SCT"]]@data) <- new_names

if (nrow(obj[["SCT"]]@scale.data) > 0) {
  rownames(obj[["SCT"]]@scale.data) <- new_names[
    match(rownames(obj[["SCT"]]@scale.data), old_names)
  ]
}

# Update variable features
VariableFeatures(obj[["SCT"]]) <- new_names[
  match(VariableFeatures(obj[["SCT"]]), old_names)
]

ovary.tabula.sapiens.harmony <- obj

################################################################################
# 2. Gate CD45+ cells (Ptprc expressing cells)
################################################################################

# Define scGate gating model to purify immune cells using PTPRC
my.scGate.model.PTPRC <- gating_model(name = "immune", signature = c("PTPRC"))

ovary.Human.TS.PTPRC <- scGate(data = ovary.tabula.sapiens.harmony, model = my.scGate.model.PTPRC, assay = "RNA", verbose = TRUE)

t.pure <- as.matrix(table(ovary.Human.TS.PTPRC$is.pure, ovary.Human.TS.PTPRC$seurat_clusters))

# Calculate the total cells in each cluster
total <- t.pure[1,] + t.pure[2,]

# Get cluster numbers where Pure cells are more than 80% of total cells in each cluster
high_purity_clusters <- which(t.pure[1,] / total > 0.8)

# Assign cell type to each seurat cluster

ovary.Human.TS.PTPRC <- SetIdent(ovary.Human.TS.PTPRC, value = "seurat_clusters")
ovary.Human.TS <- ovary.Human.TS.PTPRC

ovary.Human.TS@meta.data$celltype.level1 <- rep("nonimmune", dim(ovary.Human.TS@meta.data)[1])
ovary.Human.TS@meta.data$celltype.level1[ovary.Human.TS@meta.data$seurat_clusters == 6] <- "immune"

ovary.Human.TS.immune <- subset(ovary.Human.TS, subset = celltype.level1 == "immune")
ovary.Human.TS.nonimmune <- subset(ovary.Human.TS, subset = celltype.level1 == "nonimmune")

################################################################################
# 4. Annotate immune cells
################################################################################

DefaultAssay(ovary.Human.TS.immune) <- "SCT"

ovary.Human.TS.immune <- FindNeighbors(object = ovary.Human.TS.immune, reduction = "harmony")
ovary.Human.TS.immune <- FindClusters(ovary.Human.TS.immune, resolution = 1)

########## Normalize, SCT, dimreduc data ##########

ovary.Human.TS.immune <- NormalizeData(ovary.Human.TS.immune, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Human.TS.immune <- SCTransform(object = ovary.Human.TS.immune, vars.to.regress = c("nFeature_RNA", "pct_counts_mt"))

save(ovary.Human.TS.immune, file = paste0(Sys.Date(),"_10x_ovary_Human_TS_immune_cells_Seurat_object_clean_postSCT.RData"))

ovary.Human.TS.immune <- RunPCA(ovary.Human.TS.immune, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Human_TS_immune_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.Human.TS.immune, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Human.TS.immune[["pca"]]@stdev / sum(ovary.Human.TS.immune[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.

# Minimum of the two calculation
pcs <- min(co1, co2)

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

###############################################################################
# run dimensionality reduction algorithm
ovary.Human.TS.immune <- RunUMAP(ovary.Human.TS.immune, dims = 1:pcs)
ovary.Human.TS.immune <- FindNeighbors(ovary.Human.TS.immune, dims = 1:pcs)
ovary.Human.TS.immune <- FindClusters(object = ovary.Human.TS.immune)
ovary.Human.TS.immune <- FindClusters(object = ovary.Human.TS.immune, resolution = 0.5)

pdf(paste0(Sys.Date(),"_10x_ovary_Human_TS_immune_cells_UMAP_SeuratClustering_res_3_groupby_library.pdf"), width = 7, height = 5)
DimPlot(ovary.Human.TS.immune, reduction = "umap", group.by = "donor_id", shuffle = TRUE)
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Human_TS_immune_cells_UMAP_SeuratClustering_res_3_with_label.pdf"), width = 7, height = 5)
DimPlot(ovary.Human.TS.immune, reduction = "umap", label = TRUE)
dev.off()

###############################################################################
# 2. Use SingleR - Immgen to annotate immune cells
###############################################################################

DefaultAssay(ovary.Human.TS.immune) <- "RNA"

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Human.TS.immune)

humanAtlas <- HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

# Filter cell types with more than 20 cells in dataset
t.humanAtlas.cellcount <- as.data.frame(table(humanAtlas$label.main))
colnames(t.humanAtlas.cellcount) <- c("celltype", "count")

cellcount.less.than.20 <- t.humanAtlas.cellcount %>% 
  filter(count < 20) %>%
  dplyr::select(celltype)

# Drop cell types with less than 5 cells
humanAtlas <- humanAtlas[, !humanAtlas$label.main %in% cellcount.less.than.20$celltype]

my.singler.immgen = SingleR(test = my.SingleCellExperiment.object,
                            ref  = humanAtlas,
                            assay.type.test = 1,
                            labels = humanAtlas$label.main)

# Add Sample IDs to metadata
my.singler.immgen$meta.data$Library       =    ovary.Human.TS.immune@meta.data$Library                       

save(my.singler.immgen, file = paste(Sys.Date(),"10x_ovary_Human_TS_immune_cells_SingleR_object_Immgen.RData",sep = "_"))

###########################
# Transfer cell annotations to Seurat object

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Human.TS.immune)

# Plot annotation heatmap
annotation_col = data.frame(Labels = factor(my.singler.immgen$labels),
                            Library = as.data.frame(colData(my.SingleCellExperiment.object)[,"donor_id",drop=FALSE]))

# Transfer SingleR annotations to Seurat Object
ovary.Human.TS.immune[["SingleR.HumanAtlas"]] <- my.singler.immgen$labels

table.SingleR.annotation <- table(ovary.Human.TS.immune@meta.data$seurat_clusters, ovary.Human.TS.immune@meta.data$SingleR.HumanAtlas)

###############################################################################
# 3. Transfer final annotations
###############################################################################

ovary.Human.TS.immune <- SetIdent(ovary.Human.TS.immune, value = "seurat_clusters")

ovary.Human.TS.immune@meta.data$celltype.level2 <- rep("Myeloid",dim(ovary.Human.TS.immune@meta.data)[1])

ovary.Human.TS.immune <- SetIdent(ovary.Human.TS.immune, value = "celltype.level2")
levels(ovary.Human.TS.immune) <- c("Myeloid", "DNT", "CD8T", "B")

save(ovary.Human.TS.immune, file = paste(Sys.Date(),"10x_ovary_Human_TS_immune_cells_Seurat_with_final_annotation.RData",sep = "_"))

DefaultAssay(ovary.Human.TS.nonimmune) <- "SCT"

###############################################################################
# 2. Annotate non-immune cells
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions-1
###############################################################################

ovary.Human.TS.nonimmune <- ovary.Human.TS.nonimmune %>% RunUMAP(reduction = "harmony", dims = 1:15) %>% FindNeighbors(reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.5)

ovary.Human.TS.nonimmune <- FindClusters(ovary.Human.TS.nonimmune, resolution = 0.5)

DimPlot(ovary.Human.TS.nonimmune, label = TRUE)

# Find marker genes
hOvary.nonimmune.sct <- PrepSCTFindMarkers(ovary.Human.TS.nonimmune)
hOvary.nonimmune.markers <- FindAllMarkers(hOvary.nonimmune.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

save(hOvary.nonimmune.markers, file = paste(Sys.Date(),"10x_ovary_Human_TS_nonimmune_cells_markers.RData",sep = "_"))

# Plot published marker genes of GSE255690 nonimmune cells
# https://www.frontiersin.org/articles/10.3389/fimmu.2020.559555/full
# A lot of nonimmune cell markers are not expressed in the dataset

pdf(paste(Sys.Date(),"10X_ovary_Human_TS_nonimmune_cells_DotPlot_nonimmune_cell_marker_genes_re-clustered_object.pdf",sep = "_"), height = 15, width = 20)
DotPlot(ovary.Human.TS.nonimmune, features = c("CYP19A1", "FOXL2", "INHA", "CDH2", "AMH", "SERPINE2",   # Granulosa cell marker
                                               "STAR", "CYP17A1",                                       # Theca marker
                                               "DCN", "COL6A3", "LUM", "PDGFRA",                        # Stroma cell marker
                                               "ACTA2", "MYH11", "MCAM", "TAGLN",                       # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                               "FLT1", "VWF", "FLT4", "PROX1",                          # Endothelial cell marker (BEC + LEC)
                                               "CLDN1", "CDH1", "PAX8"))                                # Epithelial cell marker
dev.off()

###############################################################################
# 3. Annotate cell type based on marker gene expression
###############################################################################

ovary.Human.TS.nonimmune@meta.data$celltype.level2 <- rep("NA",dim(ovary.Human.TS.nonimmune@meta.data)[1])
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 0] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 1] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 2] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 3] <- "SMC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 4] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 5] <- "Theca"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 6] <- "BEC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 7] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 8] <- "Granulosa"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 9] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 10] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 11] <- "Epithelial"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 12] <- "Epithelial"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 13] <- "SMC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 14] <- "LEC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 15] <- "Stroma"

Idents(ovary.Human.TS.nonimmune) <- "Markergene.annot"
Idents(ovary.Human.TS.nonimmune) <- "seurat_clusters"


levels(ovary.Human.TS.nonimmune) <- c("Granulosa",
                                      "Theca",
                                      "Stroma",
                                      "SMC",
                                      "BEC",
                                      "LEC",
                                      "Epithelial")

###############################################################################
# 4. Use scSorter to annotate immune cells
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
###############################################################################
# Use curated marker gene list - from GSE255690 Dotplot

my.ovarian.marker.list <- vector(mode = "list", 7)
names(my.ovarian.marker.list) <- c("Granulosa",
                                   "Theca",
                                   "Stroma",
                                   "Mesenchyme",
                                   "BEC",
                                   "LEC",
                                   "Epithelial")

my.ovarian.marker.list$Granulosa <- unique(c("CYP19A1", "FOXL2", "INHA", "CDH2", "AMH", "SERPINE2"))
my.ovarian.marker.list$Theca <- unique(c( "STAR", "CYP17A1"))
my.ovarian.marker.list$Stroma <- unique(c("DCN", "COL6A3", "LUM", "PDGFRA"))
my.ovarian.marker.list$SMC <- unique(c("ACTA2", "MYH11", "MCAM", "TAGLN"))
my.ovarian.marker.list$BEC <- unique(c("FLT1", "VWF"))
my.ovarian.marker.list$LEC <- unique(c("FLT4", "PROX1"))
my.ovarian.marker.list$Epithelial <- unique(c("CLDN1", "CDH1", "PAX8"))

save(my.ovarian.marker.list, file = paste(Sys.Date(),"Ovarian_non-immune_cell_markers_for_scSorter.RData",sep = "_"))

# Filter genes that do not exist in dataset

my.ovarian.marker.list$Granulosa <- intersect(rownames(ovary.Human.TS.nonimmune), my.ovarian.marker.list$Granulosa)
my.ovarian.marker.list$Theca <- intersect(rownames(ovary.Human.TS.nonimmune), my.ovarian.marker.list$Theca)
my.ovarian.marker.list$Stroma <- intersect(rownames(ovary.Human.TS.nonimmune), my.ovarian.marker.list$Stroma)
my.ovarian.marker.list$SMC <- intersect(rownames(ovary.Human.TS.nonimmune), my.ovarian.marker.list$SMC)
my.ovarian.marker.list$BEC <- intersect(rownames(ovary.Human.TS.nonimmune), my.ovarian.marker.list$BEC)
my.ovarian.marker.list$Epithelial <- intersect(rownames(ovary.Human.TS.nonimmune), my.ovarian.marker.list$Epithelial)

### Generate anno table
anno <- data.frame(cbind(rep(names(my.ovarian.marker.list)[1],length(my.ovarian.marker.list[[1]])), my.ovarian.marker.list[[1]], rep(2,length(my.ovarian.marker.list[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.ovarian.marker.list)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.ovarian.marker.list)[i],length(my.ovarian.marker.list[[i]])),
                      "Marker" = my.ovarian.marker.list[[i]],
                      "Weight" = rep(2,length(my.ovarian.marker.list[[i]]))))
}


# Pre-process data for scSorter
topgenes <- head(VariableFeatures(ovary.Human.TS.nonimmune), 2500)
expr <- GetAssayData(ovary.Human.TS.nonimmune)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

# Transfer scSorter annotations to Seurat Object
ovary.Human.TS.nonimmune[["scSorter.labels"]] <- rts$Pred_Type

table.scSorter.annotation <- table(ovary.Human.TS.nonimmune@meta.data$seurat_clusters, ovary.Human.TS.nonimmune@meta.data$scSorter.labels)

ovary.Human.TS.nonimmune <- SetIdent(ovary.Human.TS.nonimmune, value = "scSorter.labels")

###############################################################################
# 5. Import published cell type annotation
###############################################################################

ovary.Human.TS.nonimmune@meta.data$published.annot <- "NA"

published.annotation <- read.table("/Volumes/OIProject_II/1_R/Human_TS/GSE255690_GSE255690_ovary_snRNA-seq_metadata.txt", sep = "\t", header = TRUE, row.names = 1)

# Modify the row names of the published.annotation to match the Seurat format
rownames(published.annotation) <- gsub("\\.([0-9]+)$", "-1_\\1", rownames(published.annotation))

# Extract cell IDs from the Seurat object and published.annotation
seurat_cell_ids <- rownames(ovary.Human.TS.nonimmune@meta.data)
annotation_cell_ids <- rownames(published.annotation)

# Find the overlapping cell IDs
overlap_ids <- intersect(seurat_cell_ids, annotation_cell_ids)

# Modify naming to match current system
# Create a named vector to map old names to new names
name_mapping <- c("SC" = "Stroma",
                  "TC" = "Theca",
                  "GC" = "Granulosa",
                  "IC" = "Immune",
                  "SMC" = "Mesenchyme",
                  "EpiC" = "Epithelial",
                  "LEC" = "LEC",    
                  "BEC" = "BEC"     
)

published.annotation <- published.annotation %>%
  mutate(celltype2 = ifelse(celltype2 %in% names(name_mapping), name_mapping[celltype2], celltype2))

# Transfer annotation to Seurat object
# Create a named vector of annotations from the published.annotation
annotation_vector <- setNames(published.annotation$celltype2, rownames(published.annotation))

# Extract annotations for the overlapping cell IDs
annotations_for_overlap <- annotation_vector[overlap_ids]

# Check the first few entries to ensure the values aren't NA
head(annotations_for_overlap)

# Now, assign these annotations to the Seurat object
ovary.Human.TS.nonimmune$published.annot[overlap_ids] <- annotations_for_overlap

Idents(ovary.Human.TS.nonimmune) <- "published.annot"

levels(ovary.Human.TS.nonimmune) <- c("Granulosa",
                                      "Theca",
                                      "Stroma",
                                      "Mesenchyme",
                                      "BEC",
                                      "LEC",
                                      "Epithelial",
                                      "Immune",
                                      "NA")

DotPlot(ovary.Human.TS.nonimmune, features = c("CYP19A1", "FOXL2", "INHA", "CDH2", "AMH", "SERPINE2",   # Granulosa cell marker
                                               "STAR", "CYP17A1",                                       # Theca marker
                                               "DCN", "COL6A3", "LUM", "PDGFRA",                        # Stroma cell marker
                                               "ACTA2", "MYH11", "MCAM", "TAGLN",                       # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                               "FLT1", "VWF", "FLT4", "PROX1",                          # Endothelial cell marker (BEC + LEC)
                                               "CLDN1", "CDH1", "PAX8"))                                # Epithelial cell marker

###############################################################################
# 6. Transfer final annotations
###############################################################################

ovary.Human.TS.nonimmune@meta.data$celltype.level2 <- rep("NA",dim(ovary.Human.TS.nonimmune@meta.data)[1])
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 0] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 1] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 2] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 3] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 4] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 5] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 6] <- "SMC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 7] <- "SMC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 8] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 9] <- "BEC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 10] <- "BEC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 11] <- "BEC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 12] <- "SMC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 13] <- "Theca"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 14] <- "Stroma"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 15] <- "SMC"
ovary.Human.TS.nonimmune@meta.data$celltype.level2[ovary.Human.TS.nonimmune@meta.data$seurat_clusters %in% 16] <- "Granulosa"

ovary.Human.TS.nonimmune <- SetIdent(ovary.Human.TS.nonimmune, value = "celltype.level2")

levels(ovary.Human.TS.nonimmune) <- c("Granulosa",
                                      "Theca",
                                      "Stroma",
                                      "SMC",
                                      "BEC")

pdf(paste(Sys.Date(),"10X_ovary_Human_TS_nonimmune_cells_Vlnplot_marker_genes_post_annotation.pdf",sep = "_"), height = 15, width = 15)
VlnPlot(ovary.Human.TS.nonimmune, features = c("CYP19A1", "FOXL2", "AMH",                               # Granulosa cell marker
                                               "STAR", "CYP17A1",                                       # Theca marker
                                               "DCN", "COL6A3",                                         # Stroma cell marker
                                               "ACTA2", "MYH11",                                        # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                               "FLT1", "VWF", "FLT4", "PROX1",                          # Endothelial cell marker (BEC + LEC)
                                               "CLDN1", "CDH1"))                                        # Epithelial cell marker
dev.off()

pdf(paste(Sys.Date(),"10X_ovary_Human_TS_nonimmune_cells_UMAP_post_annotation.pdf",sep = "_"), height = 6, width = 8)
DimPlot(ovary.Human.TS.nonimmune, label = TRUE)
dev.off()

save(ovary.Human.TS.nonimmune, file = paste(Sys.Date(),"10x_ovary_Human_TS_nonimmune_cells_Seurat_with_final_annotation.RData",sep = "_"))

###############################################################################
# Combine annotation
###############################################################################

# Transfer annotation
ovary.Human.TS@meta.data$celltype.level2 <- "NA"

my.ovary.Human.TS.annot <- merge(ovary.Human.TS.immune, y = ovary.Human.TS.nonimmune)

# Filter cell IDs
cell_ids_to_keep <- colnames(my.ovary.Human.TS.annot)
ovary.Human.TS.cl <- subset(ovary.Human.TS, cells = cell_ids_to_keep)

ind = match(rownames(ovary.Human.TS.cl@meta.data), rownames(my.ovary.Human.TS.annot@meta.data))
ovary.Human.TS.cl@meta.data[, "celltype.level2"] = my.ovary.Human.TS.annot@meta.data[ind, "celltype.level2"]

save(ovary.Human.TS.cl, file = paste0(Sys.Date(),"_10x_ovary_Human_TS_celltype_annotated_Seurat_object_combined_final.RData"))
save(ovary.Human.TS.nonimmune, file = paste0(Sys.Date(),"_10x_ovary_Human_TS_celltype_annotated_Seurat_object_nonimmune_cells_final.RData"))
save(ovary.Human.TS.immune, file = paste0(Sys.Date(),"_10x_ovary_Human_TS_celltype_annotated_Seurat_object_immune_cells_final.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Human_TS_celltype_annotation_session_Info.txt", sep =""))
sessionInfo()
sink()
