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

################################################################################
# Integrative ovarian aging analysis
# 10x ovary Human GSE255690 dataset
# Annotate cell type
################################################################################

################################################################################
# 1. Load dataset
################################################################################

load("/Volumes/OIProject_II/1_R/2_Integration/Human/Human_GSE255690/2025-06-23_10x_ovary_Human_GSE255690_integrated_harmony_post_clustering.RData")

################################################################################
# 2. Gate PTPRC expressing cells
################################################################################

#Manually define a simple scGate gating model to purify immune cells using a positive marker PTPRC
my.scGate.model.PTPRC <- gating_model(name = "immune", signature = c("PTPRC"))

ovary.Human.GSE255690.PTPRC <- scGate(data = ovary.Human.GSE255690.harmony, model = my.scGate.model.PTPRC, assay = "RNA", verbose = TRUE)

t.pure <- as.matrix(table(ovary.Human.GSE255690.PTPRC$is.pure, ovary.Human.GSE255690.PTPRC$seurat_clusters))

# Calculate the total cells in each cluster
total <- t.pure[1,] + t.pure[2,]

# Get cluster numbers where Pure cells are more than 80% of total cells in each cluster
high_purity_clusters <- which(t.pure[1,] / total > 0.8)

# Print cluster numbers
print(high_purity_clusters)
# 3, 5

# Assign cell type to each seurat cluster
# CD45 positive cells clusters: 3, 5

ovary.Human.GSE255690.PTPRC <- SetIdent(ovary.Human.GSE255690.PTPRC, value = "seurat_clusters")
ovary.Human.GSE255690 <- ovary.Human.GSE255690.PTPRC

ovary.Human.GSE255690@meta.data$celltype.level1 <- rep("nonimmune", dim(ovary.Human.GSE255690@meta.data)[1])
ovary.Human.GSE255690@meta.data$celltype.level1[ovary.Human.GSE255690@meta.data$seurat_clusters == 3] <- "immune"
ovary.Human.GSE255690@meta.data$celltype.level1[ovary.Human.GSE255690@meta.data$seurat_clusters == 5] <- "immune"

ovary.Human.GSE255690.immune <- subset(ovary.Human.GSE255690, subset = celltype.level1 == "immune")
ovary.Human.GSE255690.nonimmune <- subset(ovary.Human.GSE255690, subset = celltype.level1 == "nonimmune")

################################################################################
# 4. Annotate immune cells
################################################################################

DefaultAssay(ovary.Human.GSE255690.immune) <- "SCT"

ovary.Human.GSE255690.immune <- FindNeighbors(object = ovary.Human.GSE255690.immune, reduction = "harmony")
ovary.Human.GSE255690.immune <- FindClusters(ovary.Human.GSE255690.immune, resolution = 1)

########## Normalize, SCT, dimreduc data ##########

ovary.Human.GSE255690.immune <- NormalizeData(ovary.Human.GSE255690.immune, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Human.GSE255690.immune <- SCTransform(object = ovary.Human.GSE255690.immune, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Human.GSE255690.immune, file = paste0(Sys.Date(),"_10x_ovary_GSE255690_immune_cells_Seurat_object_clean_postSCT.RData"))

ovary.Human.GSE255690.immune <- RunPCA(ovary.Human.GSE255690.immune, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_GSE255690_immune_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.Human.GSE255690.immune, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Human.GSE255690.immune[["pca"]]@stdev / sum(ovary.Human.GSE255690.immune[["pca"]]@stdev) * 100

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

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_GSE255690_immune_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Human.GSE255690.immune <- RunUMAP(ovary.Human.GSE255690.immune, dims = 1:pcs)
ovary.Human.GSE255690.immune <- FindNeighbors(ovary.Human.GSE255690.immune, dims = 1:pcs)
ovary.Human.GSE255690.immune <- FindClusters(object = ovary.Human.GSE255690.immune)
ovary.Human.GSE255690.immune <- FindClusters(object = ovary.Human.GSE255690.immune, resolution = 0.5)

###############################################################################
# 2. Use SingleR - Immgen to annotate immune cells
###############################################################################

DefaultAssay(ovary.Human.GSE255690.immune) <- "RNA"

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Human.GSE255690.immune)

# HumanPrimaryCellAtlasData()

humanAtlas <- HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))

# Filter cell types with more than 20 cells in dataset
t.humanAtlas.cellcount <- as.data.frame(table(humanAtlas$label.main))
colnames(t.humanAtlas.cellcount) <- c("celltype", "count")

cellcount.less.than.20 <- t.humanAtlas.cellcount %>% 
  filter(count < 20) %>%
  select(celltype)

# Drop cell types with less than 5 cells
humanAtlas <- humanAtlas[, !humanAtlas$label.main %in% cellcount.less.than.20$celltype]

my.singler.immgen = SingleR(test = my.SingleCellExperiment.object,
                            ref  = humanAtlas,
                            assay.type.test = 1,
                            labels = humanAtlas$label.main)

# Add Sample IDs to metadata
my.singler.immgen$meta.data$Library       =    ovary.Human.GSE255690.immune@meta.data$Library                       # Library IDs

save(my.singler.immgen, file = paste(Sys.Date(),"10x_ovary_Human_GSE255690_immune_cells_SingleR_object_Immgen.RData",sep = "_"))

###########################
# Transfer cell annotations to Seurat object
my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Human.GSE255690.immune)

# Transfer SingleR annotations to Seurat Object
ovary.Human.GSE255690.immune[["SingleR.HumanAtlas"]] <- my.singler.immgen$labels

table.SingleR.annotation <- table(ovary.Human.GSE255690.immune@meta.data$seurat_clusters, ovary.Human.GSE255690.immune@meta.data$SingleR.HumanAtlas)

###############################################################################
# 3. Assess expression of marker genes
###############################################################################

DefaultAssay(ovary.Human.GSE255690.immune) <- "SCT"

pdf(paste(Sys.Date(),"10X_ovary_Human_GSE255690_immune_cells_VlnPlot_T_monocyte_mph_marker_genes.pdf",sep = "_"), height = 10, width = 10)
VlnPlot(ovary.Human.GSE255690.immune, features = c("CD3E", "CD247",             # General T cell marker
                                                   "CD8A", "CD8B", "CD4",       # CD8 and CD4 T cell marker
                                                   "CD163", "MRC1", "CD14",     # Monocyte/Mph markers
                                                   "CD38", "ITGAM"))                     # B markers
dev.off()

# Transfer annotations

ovary.Human.GSE255690.immune <- SetIdent(ovary.Human.GSE255690.immune, value = "seurat_clusters")

ovary.Human.GSE255690.immune@meta.data$celltype.level2 <- rep("NA",dim(ovary.Human.GSE255690.immune@meta.data)[1])
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 0]<- "CD8T"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 1]<- "CD8T"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 2]<- "DNT"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 3]<- "CD8T"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 4]<- "CD8T"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 5]<- "Myeloid"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 6]<- "Myeloid"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 7]<- "Myeloid"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 8]<- "DNT"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 9]<- "Myeloid"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 10]<- "B"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 11]<- "CD8T"
ovary.Human.GSE255690.immune@meta.data$celltype.level2[ovary.Human.GSE255690.immune@meta.data$seurat_clusters %in% 12]<- "Myeloid"

save(ovary.Human.GSE255690.immune, file = paste(Sys.Date(),"10x_ovary_Human_GSE255690_immune_cells_Seurat_with_final_annotation.RData",sep = "_"))

###############################################################################
# 2. Annotate non-immune cells
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions-1
###############################################################################

ovary.Human.GSE255690.nonimmune <- ovary.Human.GSE255690.nonimmune %>% RunUMAP(reduction = "harmony", dims = 1:15) %>% FindNeighbors(reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1) %>% identity()

# Find marker genes
hOvary.nonimmune.sct <- PrepSCTFindMarkers(ovary.Human.GSE255690.nonimmune)
hOvary.nonimmune.markers <- FindAllMarkers(hOvary.nonimmune.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

save(hOvary.nonimmune.markers, file = paste(Sys.Date(),"10x_ovary_Human_GSE255690_nonimmune_cells_markers.RData",sep = "_"))

# Plot published marker genes of GSE255690 nonimmune cells
# https://www.frontiersin.org/articles/10.3389/fimmu.2020.559555/full
# A lot of nonimmune cell markers are not expressed in the dataset

pdf(paste(Sys.Date(),"10X_ovary_Human_GSE255690_nonimmune_cells_DotPlot_nonimmune_cell_marker_genes_re-clustered_object.pdf",sep = "_"), height = 15, width = 20)
DotPlot(ovary.Human.GSE255690.nonimmune, features = c("CYP19A1", "FOXL2", "INHA", "CDH2", "AMH", "SERPINE2",   # Granulosa cell marker
                                                      "STAR", "CYP17A1",                                       # Theca marker
                                                      "DCN", "COL6A3", "LUM", "PDGFRA",                        # Stroma cell marker
                                                      "ACTA2", "MYH11", "MCAM", "TAGLN",                       # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                                      "FLT1", "VWF", "FLT4", "PROX1",                          # Endothelial cell marker (BEC + LEC)
                                                      "CLDN1", "CDH1", "PAX8"))                                # Epithelial cell marker
dev.off()

###############################################################################
# 3. Annotate cell type based on marker gene expression
###############################################################################

ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot <- rep("NA",dim(ovary.Human.GSE255690.nonimmune@meta.data)[1])
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 0] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 1] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 2] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 3] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 4] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 5] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 6] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 7] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 8] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 9] <- "BEC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 10] <- "BEC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 11] <- "BEC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 12] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 13] <- "Theca"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 14] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 15] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$Markergene.annot[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 16] <- "Granulosa"

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

my.ovarian.marker.list$Granulosa <- intersect(rownames(ovary.Human.GSE255690.nonimmune), my.ovarian.marker.list$Granulosa)
my.ovarian.marker.list$Theca <- intersect(rownames(ovary.Human.GSE255690.nonimmune), my.ovarian.marker.list$Theca)
my.ovarian.marker.list$Stroma <- intersect(rownames(ovary.Human.GSE255690.nonimmune), my.ovarian.marker.list$Stroma)
my.ovarian.marker.list$SMC <- intersect(rownames(ovary.Human.GSE255690.nonimmune), my.ovarian.marker.list$SMC)
my.ovarian.marker.list$BEC <- intersect(rownames(ovary.Human.GSE255690.nonimmune), my.ovarian.marker.list$BEC)
my.ovarian.marker.list$Epithelial <- intersect(rownames(ovary.Human.GSE255690.nonimmune), my.ovarian.marker.list$Epithelial)

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
topgenes <- head(VariableFeatures(ovary.Human.GSE255690.nonimmune), 2500)
expr <- GetAssayData(ovary.Human.GSE255690.nonimmune)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

### Viewing Results - The cell type assignment results are stored in the Pred_Type vector.
print(table(rts$Pred_Type))
#    BEC Epithelial  Granulosa        LEC Mesenchyme     Stroma      Theca    Unknown 
#   1305        517       1599        310       1735      16536       1242        566 

# Transfer scSorter annotations to Seurat Object
ovary.Human.GSE255690.nonimmune[["scSorter.labels"]] <- rts$Pred_Type

table.scSorter.annotation <- table(ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters, ovary.Human.GSE255690.nonimmune@meta.data$scSorter.labels)

###############################################################################
# 5. Import published cell type annotation
###############################################################################

ovary.Human.GSE255690.nonimmune@meta.data$published.annot <- "NA"

published.annotation <- read.table("/Volumes/OIProject_II/1_R/Human_GSE255690/GSE255690_GSE255690_ovary_snRNA-seq_metadata.txt", sep = "\t", header = TRUE, row.names = 1)
dim(published.annotation)
# [1] 42568     4

# Modify the row names of the published.annotation to match the Seurat format
rownames(published.annotation) <- gsub("\\.([0-9]+)$", "-1_\\1", rownames(published.annotation))

# Extract cell IDs from the Seurat object and published.annotation
seurat_cell_ids <- rownames(ovary.Human.GSE255690.nonimmune@meta.data)
annotation_cell_ids <- rownames(published.annotation)

# Find the overlapping cell IDs
overlap_ids <- intersect(seurat_cell_ids, annotation_cell_ids)

# Display the number of overlapping cell IDs
length(overlap_ids)
# [1] 31825

# Modify naming to match current system
# Create a named vector to map old names to new names
name_mapping <- c("SC" = "Stroma",
                  "TC" = "Theca",
                  "GC" = "Granulosa",
                  "IC" = "Immune",
                  "SMC" = "Mesenchyme",
                  "EpiC" = "Epithelial",
                  "LEC" = "LEC",      # Since LEC remains the same
                  "BEC" = "BEC"       # Since BEC remains the same
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
ovary.Human.GSE255690.nonimmune$published.annot[overlap_ids] <- annotations_for_overlap

Idents(ovary.Human.GSE255690.nonimmune) <- "published.annot"

levels(ovary.Human.GSE255690.nonimmune) <- c("Granulosa",
                                             "Theca",
                                             "Stroma",
                                             "Mesenchyme",
                                             "BEC",
                                             "LEC",
                                             "Epithelial",
                                             "Immune",
                                             "NA")

DotPlot(ovary.Human.GSE255690.nonimmune, features = c("CYP19A1", "FOXL2", "INHA", "CDH2", "AMH", "SERPINE2",   # Granulosa cell marker
                                                      "STAR", "CYP17A1",                                       # Theca marker
                                                      "DCN", "COL6A3", "LUM", "PDGFRA",                        # Stroma cell marker
                                                      "ACTA2", "MYH11", "MCAM", "TAGLN",                       # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                                      "FLT1", "VWF", "FLT4", "PROX1",                          # Endothelial cell marker (BEC + LEC)
                                                      "CLDN1", "CDH1", "PAX8"))                                # Epithelial cell marker

DimPlot(ovary.Human.GSE255690.nonimmune)

table(ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters, ovary.Human.GSE255690.nonimmune@meta.data$published.annot)

###############################################################################
# 4. Assess cell type annotation by Marker gene expressio, scSorter and published annotation
###############################################################################

# Export annotations as dataframe - Markergene.annotation, scSorter.labels and published.annot

annotations <- ovary.Human.GSE255690.nonimmune@meta.data[, c("Markergene.annot", "scSorter.labels", "published.annot")]
filtered_annotations <- annotations[!is.na(annotations$published.annot), ]

# Compare agreement among three annotations

percentage_agreement <- function(vec1, vec2) {
  sum(vec1 == vec2) / length(vec1) * 100
}

markergene_vs_scSorter <- percentage_agreement(filtered_annotations$Markergene.annot, filtered_annotations$scSorter.labels)          # [1] 82.35615
markergene_vs_publishedannot <- percentage_agreement(filtered_annotations$Markergene.annot, filtered_annotations$published.annot)    # [1] 91.94876
scSorter_vs_publishedannot <- percentage_agreement(filtered_annotations$scSorter.labels, filtered_annotations$published.annot)       # [1] 84.9391

Idents(ovary.Human.GSE255690.nonimmune) <- "Markergene.annot"
DimPlot(ovary.Human.GSE255690.nonimmune)

Idents(ovary.Human.GSE255690.nonimmune) <- "scSorter.labels"
DimPlot(ovary.Human.GSE255690.nonimmune)

Idents(ovary.Human.GSE255690.nonimmune) <- "published.annot"
DimPlot(ovary.Human.GSE255690.nonimmune)

###############################################################################
# 6. Transfer final annotations
###############################################################################

ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2 <- rep("NA",dim(ovary.Human.GSE255690.nonimmune@meta.data)[1])
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 0] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 1] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 2] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 3] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 4] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 5] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 6] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 7] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 8] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 9] <- "BEC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 10] <- "BEC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 11] <- "BEC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 12] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 13] <- "Theca"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 14] <- "Stroma"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 15] <- "SMC"
ovary.Human.GSE255690.nonimmune@meta.data$celltype.level2[ovary.Human.GSE255690.nonimmune@meta.data$seurat_clusters %in% 16] <- "Granulosa"

save(ovary.Human.GSE255690.nonimmune, file = paste(Sys.Date(),"10x_ovary_Human_GSE255690_nonimmune_cells_Seurat_with_final_annotation.RData",sep = "_"))

###############################################################################
# Combine annotation
###############################################################################

# Transfer annotation
ovary.Human.GSE255690@meta.data$celltype.level2 <- "NA"

my.ovary.Human.GSE255690.annot <- merge(ovary.Human.GSE255690.immune, y = ovary.Human.GSE255690.nonimmune)

# Filter cell IDs
cell_ids_to_keep <- colnames(my.ovary.Human.GSE255690.annot)
ovary.Human.GSE255690.cl <- subset(ovary.Human.GSE255690, cells = cell_ids_to_keep)

ind = match(rownames(ovary.Human.GSE255690.cl@meta.data), rownames(my.ovary.Human.GSE255690.annot@meta.data))
ovary.Human.GSE255690.cl@meta.data[, "celltype.level2"] = my.ovary.Human.GSE255690.annot@meta.data[ind, "celltype.level2"]
head(ovary.Human.GSE255690.cl@meta.data)

save(ovary.Human.GSE255690.cl, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_celltype_annotated_Seurat_object_combined_final.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Human_GSE255690_cell_type_annotation_session_info.txt", sep =""))
sessionInfo()
sink()
