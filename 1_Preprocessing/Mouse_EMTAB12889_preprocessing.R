options(stringsAsFactors = FALSE)

library('celda')
library("singleCellTK")
library('Seurat')
library(bitops)
library(sctransform)
library(clustree)
library(scales)
library(dplyr)
library(DoubletFinder)
library(scds)
library(scater)
library(bitops)

# rm(list = ls())

################################################################################
# Integrative ovarian aging analysis
# 10x ovary Mouse E-MTAB-12889
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.Ovary_9M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_1/outs/filtered_feature_bc_matrix/")
counts.Ovary_9M_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_2/outs/filtered_feature_bc_matrix/")
counts.Ovary_9M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_3/outs/filtered_feature_bc_matrix/")
counts.Ovary_9M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_4/outs/filtered_feature_bc_matrix/")
counts.Ovary_12M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_1/outs/filtered_feature_bc_matrix/")
counts.Ovary_12M_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_2/outs/filtered_feature_bc_matrix/")
counts.Ovary_12M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_3/outs/filtered_feature_bc_matrix/")
counts.Ovary_12M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_4/outs/filtered_feature_bc_matrix/")
counts.Ovary_12M_5 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_5/outs/filtered_feature_bc_matrix/")
counts.Ovary_15M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_15M_1/outs/filtered_feature_bc_matrix/")
counts.Ovary_15M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_15M_3/outs/filtered_feature_bc_matrix/")
counts.Ovary_15M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_15M_4/outs/filtered_feature_bc_matrix/")

counts.raw.Ovary_9M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_1/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_9M_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_2/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_9M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_3/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_9M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_9M_4/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_12M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_1/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_12M_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_2/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_12M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_3/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_12M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_4/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_12M_5 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_12M_5/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_15M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_15M_1/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_15M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_15M_3/outs/raw_feature_bc_matrix/")
counts.raw.Ovary_15M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_EMTAB12889/2_Cellranger/Ovary_15M_4/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX

sce.Ovary_9M_1 <- SingleCellExperiment(list(counts = counts.Ovary_9M_1))
sce.Ovary_9M_2 <- SingleCellExperiment(list(counts = counts.Ovary_9M_2))
sce.Ovary_9M_3 <- SingleCellExperiment(list(counts = counts.Ovary_9M_3))
sce.Ovary_9M_4 <- SingleCellExperiment(list(counts = counts.Ovary_9M_4))
sce.Ovary_12M_1 <- SingleCellExperiment(list(counts = counts.Ovary_12M_1))
sce.Ovary_12M_2 <- SingleCellExperiment(list(counts = counts.Ovary_12M_2))
sce.Ovary_12M_3 <- SingleCellExperiment(list(counts = counts.Ovary_12M_3))
sce.Ovary_12M_4 <- SingleCellExperiment(list(counts = counts.Ovary_12M_4))
sce.Ovary_12M_5 <- SingleCellExperiment(list(counts = counts.Ovary_12M_5))
sce.Ovary_15M_1 <- SingleCellExperiment(list(counts = counts.Ovary_15M_1))
sce.Ovary_15M_3 <- SingleCellExperiment(list(counts = counts.Ovary_15M_3))
sce.Ovary_15M_4 <- SingleCellExperiment(list(counts = counts.Ovary_15M_4))

sce.raw.Ovary_9M_1 <- SingleCellExperiment(list(counts = counts.raw.Ovary_9M_1))
sce.raw.Ovary_9M_2 <- SingleCellExperiment(list(counts = counts.raw.Ovary_9M_2))
sce.raw.Ovary_9M_3 <- SingleCellExperiment(list(counts = counts.raw.Ovary_9M_3))
sce.raw.Ovary_9M_4 <- SingleCellExperiment(list(counts = counts.raw.Ovary_9M_4))
sce.raw.Ovary_12M_1 <- SingleCellExperiment(list(counts = counts.raw.Ovary_12M_1))
sce.raw.Ovary_12M_2 <- SingleCellExperiment(list(counts = counts.raw.Ovary_12M_2))
sce.raw.Ovary_12M_3 <- SingleCellExperiment(list(counts = counts.raw.Ovary_12M_3))
sce.raw.Ovary_12M_4 <- SingleCellExperiment(list(counts = counts.raw.Ovary_12M_4))
sce.raw.Ovary_12M_5 <- SingleCellExperiment(list(counts = counts.raw.Ovary_12M_5))
sce.raw.Ovary_15M_1 <- SingleCellExperiment(list(counts = counts.raw.Ovary_15M_1))
sce.raw.Ovary_15M_3 <- SingleCellExperiment(list(counts = counts.raw.Ovary_15M_3))
sce.raw.Ovary_15M_4 <- SingleCellExperiment(list(counts = counts.raw.Ovary_15M_4))

sce.Ovary_9M_1 <- decontX(sce.Ovary_9M_1, background = sce.raw.Ovary_9M_1)
sce.Ovary_9M_2 <- decontX(sce.Ovary_9M_2, background = sce.raw.Ovary_9M_2)
sce.Ovary_9M_3 <- decontX(sce.Ovary_9M_3, background = sce.raw.Ovary_9M_3)
sce.Ovary_9M_4 <- decontX(sce.Ovary_9M_4, background = sce.raw.Ovary_9M_4)
sce.Ovary_12M_1 <- decontX(sce.Ovary_12M_1, background = sce.raw.Ovary_12M_1)
sce.Ovary_12M_2 <- decontX(sce.Ovary_12M_2, background = sce.raw.Ovary_12M_2)
sce.Ovary_12M_3 <- decontX(sce.Ovary_12M_3, background = sce.raw.Ovary_12M_3)
sce.Ovary_12M_4 <- decontX(sce.Ovary_12M_4, background = sce.raw.Ovary_12M_4)
sce.Ovary_12M_5 <- decontX(sce.Ovary_12M_5, background = sce.raw.Ovary_12M_5)
sce.Ovary_15M_1 <- decontX(sce.Ovary_15M_1, background = sce.raw.Ovary_15M_1)
sce.Ovary_15M_3 <- decontX(sce.Ovary_15M_3, background = sce.raw.Ovary_15M_3)
sce.Ovary_15M_4 <- decontX(sce.Ovary_15M_4, background = sce.raw.Ovary_15M_4)

# Save DecontX result
save(sce.Ovary_9M_1, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_9M_1.RData"))
save(sce.Ovary_9M_2, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_9M_2.RData"))
save(sce.Ovary_9M_3, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_9M_3.RData"))
save(sce.Ovary_9M_4, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_9M_4.RData"))
save(sce.Ovary_12M_1, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_12M_1.RData"))
save(sce.Ovary_12M_2, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_12M_2.RData"))
save(sce.Ovary_12M_3, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_12M_3.RData"))
save(sce.Ovary_12M_4, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_12M_4.RData"))
save(sce.Ovary_12M_5, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_12M_5.RData"))
save(sce.Ovary_15M_1, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_15M_1.RData"))
save(sce.Ovary_15M_3, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_15M_3.RData"))
save(sce.Ovary_15M_4, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_DecontX_SCE_object_sce.Ovary_15M_4.RData"))

# Plot UMAP
UMAP.Ovary_9M_1 <- reducedDim(sce.Ovary_9M_1, "decontX_UMAP")
UMAP.Ovary_9M_2 <- reducedDim(sce.Ovary_9M_2, "decontX_UMAP")
UMAP.Ovary_9M_3 <- reducedDim(sce.Ovary_9M_3, "decontX_UMAP")
UMAP.Ovary_9M_4 <- reducedDim(sce.Ovary_9M_4, "decontX_UMAP")
UMAP.Ovary_12M_1 <- reducedDim(sce.Ovary_12M_1, "decontX_UMAP")
UMAP.Ovary_12M_2 <- reducedDim(sce.Ovary_12M_2, "decontX_UMAP")
UMAP.Ovary_12M_3 <- reducedDim(sce.Ovary_12M_3, "decontX_UMAP")
UMAP.Ovary_12M_4 <- reducedDim(sce.Ovary_12M_4, "decontX_UMAP")
UMAP.Ovary_12M_5 <- reducedDim(sce.Ovary_12M_5, "decontX_UMAP")
UMAP.Ovary_15M_1 <- reducedDim(sce.Ovary_15M_1, "decontX_UMAP")
UMAP.Ovary_15M_3 <- reducedDim(sce.Ovary_15M_3, "decontX_UMAP")
UMAP.Ovary_15M_4 <- reducedDim(sce.Ovary_15M_4, "decontX_UMAP")

plotDecontXContamination(sce.Ovary_9M_1)
plotDecontXContamination(sce.Ovary_9M_2)
plotDecontXContamination(sce.Ovary_9M_3)
plotDecontXContamination(sce.Ovary_9M_4)
plotDecontXContamination(sce.Ovary_12M_1)
plotDecontXContamination(sce.Ovary_12M_2)
plotDecontXContamination(sce.Ovary_12M_3)
plotDecontXContamination(sce.Ovary_12M_4)
plotDecontXContamination(sce.Ovary_12M_5)
plotDecontXContamination(sce.Ovary_15M_1)
plotDecontXContamination(sce.Ovary_15M_3)
plotDecontXContamination(sce.Ovary_15M_4)

# Plot example box plots of contamination scores
boxplot(sce.Ovary_9M_1$decontX_contamination)
boxplot(sce.Ovary_9M_2$decontX_contamination)
boxplot(sce.Ovary_9M_3$decontX_contamination)
boxplot(sce.Ovary_9M_4$decontX_contamination)
boxplot(sce.Ovary_12M_1$decontX_contamination)
boxplot(sce.Ovary_12M_2$decontX_contamination)
boxplot(sce.Ovary_12M_3$decontX_contamination)
boxplot(sce.Ovary_12M_4$decontX_contamination)
boxplot(sce.Ovary_12M_5$decontX_contamination)
boxplot(sce.Ovary_15M_1$decontX_contamination)
boxplot(sce.Ovary_15M_3$decontX_contamination)
boxplot(sce.Ovary_15M_4$decontX_contamination)

# Create Seurat objects
seurat.Ovary_9M_1 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_9M_1)), meta.data=as.data.frame(colData(sce.Ovary_9M_1)))
seurat.Ovary_9M_2 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_9M_2)), meta.data=as.data.frame(colData(sce.Ovary_9M_2)))
seurat.Ovary_9M_3 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_9M_3)), meta.data=as.data.frame(colData(sce.Ovary_9M_3)))
seurat.Ovary_9M_4 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_9M_4)), meta.data=as.data.frame(colData(sce.Ovary_9M_4)))
seurat.Ovary_12M_1 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_12M_1)), meta.data=as.data.frame(colData(sce.Ovary_12M_1)))
seurat.Ovary_12M_2 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_12M_2)), meta.data=as.data.frame(colData(sce.Ovary_12M_2)))
seurat.Ovary_12M_3 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_12M_3)), meta.data=as.data.frame(colData(sce.Ovary_12M_3)))
seurat.Ovary_12M_4 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_12M_4)), meta.data=as.data.frame(colData(sce.Ovary_12M_4)))
seurat.Ovary_12M_5 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_12M_5)), meta.data=as.data.frame(colData(sce.Ovary_12M_5)))
seurat.Ovary_15M_1 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_15M_1)), meta.data=as.data.frame(colData(sce.Ovary_15M_1)))
seurat.Ovary_15M_3 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_15M_3)), meta.data=as.data.frame(colData(sce.Ovary_15M_3)))
seurat.Ovary_15M_4 <- CreateSeuratObject(round(decontXcounts(sce.Ovary_15M_4)), meta.data=as.data.frame(colData(sce.Ovary_15M_4)))

ovary.EMTAB12889 <- merge(seurat.Ovary_9M_1, 
                           y =  c(seurat.Ovary_9M_2, seurat.Ovary_9M_3, seurat.Ovary_9M_4, 
                                  seurat.Ovary_12M_1, seurat.Ovary_12M_2, seurat.Ovary_12M_3, seurat.Ovary_12M_4, seurat.Ovary_12M_5,
                                  seurat.Ovary_15M_1, seurat.Ovary_15M_3, seurat.Ovary_15M_4), 
                           add.cell.ids = c("Ovary_9M_1", "Ovary_9M_2", "Ovary_9M_3", "Ovary_9M_4", 
                                            "Ovary_12M_1", "Ovary_12M_2", "Ovary_12M_3", "Ovary_12M_4", "Ovary_12M_5", 
                                            "Ovary_15M_1", "Ovary_15M_3", "Ovary_15M_4"), 
                  project = "10x_ovary_EMTAB12889")

ovary.EMTAB12889
# An object of class Seurat 
# 32285 features across 44176 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

save(ovary.EMTAB12889, file = paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB12889_post-DecontX_Seurat_object.RData"))

# Remove intermediate files
rm(counts.Ovary_9M_1, counts.Ovary_9M_2, counts.Ovary_9M_3, counts.Ovary_9M_4, 
   counts.Ovary_12M_1, counts.Ovary_12M_2, counts.Ovary_12M_3, counts.Ovary_12M_4, counts.Ovary_12M_5, 
   counts.Ovary_15M_1, counts.Ovary_15M_3, counts.Ovary_15M_4, 
   counts.raw.Ovary_9M_1, counts.raw.Ovary_9M_2, counts.raw.Ovary_9M_3, counts.raw.Ovary_9M_4, 
   counts.raw.Ovary_12M_1, counts.raw.Ovary_12M_2, counts.raw.Ovary_12M_3, counts.raw.Ovary_12M_4, counts.raw.Ovary_12M_5, 
   counts.raw.Ovary_15M_1, counts.raw.Ovary_15M_3, counts.raw.Ovary_15M_4,
   sce.Ovary_9M_1, sce.Ovary_9M_2, sce.Ovary_9M_3, sce.Ovary_9M_4, 
   sce.Ovary_12M_1, sce.Ovary_12M_2, sce.Ovary_12M_3, sce.Ovary_12M_4, sce.Ovary_12M_5, 
   sce.Ovary_15M_1, sce.Ovary_15M_3, sce.Ovary_15M_4, 
   sce.raw.Ovary_9M_1, sce.raw.Ovary_9M_2, sce.raw.Ovary_9M_3, sce.raw.Ovary_9M_4, 
   sce.raw.Ovary_12M_1, sce.raw.Ovary_12M_2, sce.raw.Ovary_12M_3, sce.raw.Ovary_12M_4, sce.raw.Ovary_12M_5, 
   sce.raw.Ovary_15M_1, sce.raw.Ovary_15M_3, sce.raw.Ovary_15M_4,
   seurat.Ovary_9M_1, seurat.Ovary_9M_2, seurat.Ovary_9M_3, seurat.Ovary_9M_4, 
   seurat.Ovary_12M_1, seurat.Ovary_12M_2, seurat.Ovary_12M_3, seurat.Ovary_12M_4, seurat.Ovary_12M_5, 
   seurat.Ovary_15M_1, seurat.Ovary_15M_3, seurat.Ovary_15M_4)

################################################################################
# 2. Add metadata to Seurat object
################################################################################

# create Library label
Library <- rep("NA", length(colnames(ovary.EMTAB12889@assays$RNA)))
Library[grep("Ovary_9M_1", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_9M_1"
Library[grep("Ovary_9M_2", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_9M_2"
Library[grep("Ovary_9M_3", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_9M_3"
Library[grep("Ovary_9M_4", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_9M_4"
Library[grep("Ovary_12M_1", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_12M_1"
Library[grep("Ovary_12M_2", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_12M_2"
Library[grep("Ovary_12M_3", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_12M_3"
Library[grep("Ovary_12M_4", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_12M_4"
Library[grep("Ovary_12M_5", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_12M_5"
Library[grep("Ovary_15M_1", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_15M_1"
Library[grep("Ovary_15M_3", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_15M_3"
Library[grep("Ovary_15M_4", colnames(ovary.EMTAB12889@assays$RNA))] <- "Ovary_15M_4"

Library <- data.frame(Library)
rownames(Library) <- colnames(ovary.EMTAB12889@assays$RNA)

# create Age label
Age <- rep("NA", length(colnames(ovary.EMTAB12889@assays$RNA)))
Age[grep("9M", colnames(ovary.EMTAB12889@assays$RNA))] <- "9M"
Age[grep("12M", colnames(ovary.EMTAB12889@assays$RNA))] <- "12M"
Age[grep("15M", colnames(ovary.EMTAB12889@assays$RNA))] <- "15M"

Age <- data.frame(Age)
rownames(Age) <- colnames(ovary.EMTAB12889@assays$RNA)

# update Seurat with metadata
ovary.EMTAB12889 <- AddMetaData(object = ovary.EMTAB12889, metadata = as.vector(Library), col.name = "Library")
ovary.EMTAB12889 <- AddMetaData(object = ovary.EMTAB12889, metadata = as.vector(Age), col.name = "Age")

table(ovary.EMTAB12889@meta.data$Library)
#  Ovary_12M_1 Ovary_12M_2 Ovary_12M_3 Ovary_12M_4 Ovary_12M_5 Ovary_15M_1 Ovary_15M_3 Ovary_15M_4  Ovary_9M_1  Ovary_9M_2  Ovary_9M_3  Ovary_9M_4 
#         2960        4688        5156        2374        4553        4500        4456        4213        2159        1729        4410        2978 

table(ovary.EMTAB12889@meta.data$Age)
#   12M   15M    9M 
# 19731 13169 11276

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500
#                    500 < nCount_RNA < 1e5
#                    percent.mito < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.EMTAB12889 <- SetIdent(ovary.EMTAB12889, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.EMTAB12889@assays$RNA)
num.cells <- Matrix::rowSums(ovary.EMTAB12889@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.EMTAB12889.filt <- subset(ovary.EMTAB12889, features = genes.use)
ovary.EMTAB12889.filt
# An object of class Seurat 
# 19650 features across 44176 samples within 1 assay 
# Active assay: RNA (19650 features, 0 variable features)

############### QC - mitochondrial genes ###############

ovary.EMTAB12889.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.EMTAB12889.filt, pattern = "^mt-")
head(ovary.EMTAB12889.filt@meta.data)
#                           orig.ident nCount_RNA nFeature_RNA decontX_contamination decontX_clusters Library Age percent.mito
# Ovary_9M_1_AAACCCAAGATCGACG-1 SeuratProject       3708         1622           0.136968267                1 Ovary_9M_1  9M    7.2545847
# Ovary_9M_1_AAACCCAAGCGGTAGT-1 SeuratProject        743          495           0.178274486                1 Ovary_9M_1  9M    6.7294751
# Ovary_9M_1_AAACCCACAACGGCCT-1 SeuratProject       5397         2101           0.001004955                1 Ovary_9M_1  9M   14.3598295
# Ovary_9M_1_AAACCCACAAGGTCTT-1 SeuratProject       2034         1044           0.039047662                1 Ovary_9M_1  9M    0.2949853
# Ovary_9M_1_AAACCCACACCTGTCT-1 SeuratProject       1573          867           0.209096912                1 Ovary_9M_1  9M    2.8607756
# Ovary_9M_1_AAACGAAAGCATGGGT-1 SeuratProject      18306         4469           0.007743142                2 Ovary_9M_1  9M    2.0157325

pdf(paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.EMTAB12889.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.EMTAB12889.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.EMTAB12889.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.EMTAB12889.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.EMTAB12889.filt <- subset(ovary.EMTAB12889.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)
ovary.EMTAB12889.filt
# An object of class Seurat 
# 19650 features across 31135 samples within 1 assay 
# Active assay: RNA (19650 features, 0 variable features)

pdf(paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.EMTAB12889.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.EMTAB12889.filt, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.EMTAB12889.filt <- NormalizeData(ovary.EMTAB12889.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.EMTAB12889.filt <- SCTransform(object = ovary.EMTAB12889.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library"))

save(ovary.EMTAB12889.filt, file = paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB128899_Seurat_object_postSCT.RData"))

ovary.EMTAB12889.filt <- RunPCA(ovary.EMTAB12889.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Mouse_EMTAB12889_ElbowPlot.pdf"))
ElbowPlot(ovary.EMTAB12889.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.EMTAB12889.filt[["pca"]]@stdev / sum(ovary.EMTAB12889.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 16

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 16

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_Mouse_EMTAB12889_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.EMTAB12889.filt <- RunUMAP(ovary.EMTAB12889.filt, dims = 1:pcs)
ovary.EMTAB12889.filt <- FindNeighbors(ovary.EMTAB12889.filt, dims = 1:pcs)
ovary.EMTAB12889.filt <- FindClusters(object = ovary.EMTAB12889.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB12889_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.EMTAB12889.filt, reduction = "umap")
dev.off()

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
table(ovary.EMTAB12889.filt$Library)

# Ovary_12M_1 Ovary_12M_2 Ovary_12M_3 Ovary_12M_4 Ovary_12M_5 Ovary_15M_1 Ovary_15M_3 Ovary_15M_4  Ovary_9M_1  Ovary_9M_2  Ovary_9M_3  Ovary_9M_4 
#        2343        3255        4237        1914        2437        3088        3087        3233         956         828        3197        2560 

# Clean memory of intermediates
rm(ovary.EMTAB12889)

################################################################################
# 4. Doublet detection part 1: DoubletFinder
################################################################################

########## 0. Doublet rate calculation ########## 
##################################################
#### 0. Assumed doublet information/to calculate %age for prediction
# Targeted Cell Recovery  # of Cells Loaded	Barcodes Detected	Singlets Multiplets	Multiplet Rate
#  3,000	           4,950           ~3,000        ~2,900	   ~80	    ~2.4%
#  4,000	           6,600           ~3,900	       ~3,800	   ~140     ~3.2%
#  5,000	           8,250           ~4,800        ~4,600	   ~210     ~4.0%
#  6,000            9,900           ~5,700	       ~5,400	   ~300	    ~4.8%
#  7,000            11,550	         ~6,600	       ~6,200	   ~400	    ~5.6%
#  8,000            13,200	         ~7,500	       ~7,000	   ~510	    ~6.4%
#  9,000            14,850	         ~8,400	       ~7,700	   ~640	    ~7.2%
# 10,000            16,500	         ~9,200	       ~8,400	   ~780	    ~8.0%
# 12,000            19,800	         ~10,900	     ~9,800	   ~1,100   ~9.6%
# 14,000            23,100	         ~12,500	     ~11,000	 ~1,500   ~11.2%
# 16,000            26,400	         ~14,000	     ~12,100	 ~1,900   ~12.8%
# 18,000            29,700	         ~15,500	     ~13,100	 ~2,300   ~14.4%
# 20,000            33,000	         ~16,900	     ~14,100	 ~2,800   ~16.0%

pred.10x.dblt <- data.frame( "cell_number" = c(3000,4000,5000,6000,7000,8000,9000, 10000, 12000, 14000, 16000, 18000, 20000),
                             "dblt_rate"   = c(2.4 ,3.2 ,4.0 ,4.8 ,5.6 ,6.4 ,7.2 , 8.0  , 9.6, 11.2, 12.8, 14.4, 16.0))

pred_dblt_lm <- lm(dblt_rate ~ cell_number, data = pred.10x.dblt)

pdf(paste0(Sys.Date(), "_10x_cell_number_vs_doublet_rate.pdf"))
plot(dblt_rate ~ cell_number, data = pred.10x.dblt)
abline(pred_dblt_lm, col = "red", lty = "dashed")
dev.off()

# Use distinct recovered cell count number to estimate doublets in dataset
#predict(pred_dblt_lm, data.frame("cell_number" = 6000))

########################################################
#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
ovary.EMTAB12889.filt.list <- SplitObject(ovary.EMTAB12889.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.EMTAB12889.filt.list, ncol))))/100
pred.dblt.rate
# Ovary_9M_1  Ovary_9M_2  Ovary_9M_3  Ovary_9M_4 Ovary_12M_1 Ovary_12M_2 Ovary_12M_3 Ovary_12M_4 Ovary_12M_5 Ovary_15M_1 Ovary_15M_3 Ovary_15M_4 
#  0.0084128   0.0072864   0.0281336   0.0225280   0.0206184   0.0286440   0.0372856   0.0168432   0.0214456   0.0271744   0.0271656   0.0284504 

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.EMTAB12889.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.EMTAB12889 <- paramSweep_v3(ovary.EMTAB12889.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.EMTAB12889    <- summarizeSweep(sweep.res.list.ovary.EMTAB12889, GT = FALSE)
  bcmvn.ovary.EMTAB12889          <- find.pK(sweep.stats.ovary.EMTAB12889)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.EMTAB12889 <- as.numeric(as.character(bcmvn.ovary.EMTAB12889[as.numeric(bcmvn.ovary.EMTAB12889$pK[bcmvn.ovary.EMTAB12889$BCmetric == max(bcmvn.ovary.EMTAB12889$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.EMTAB12889.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                           ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.EMTAB12889.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.EMTAB12889.filt.list[[i]] <- doubletFinder_v3(ovary.EMTAB12889.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.EMTAB12889, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.EMTAB12889.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.EMTAB12889.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.EMTAB12889.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.EMTAB12889.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.EMTAB12889.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_AC",names(ovary.EMTAB12889.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.EMTAB12889.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.EMTAB12889.DFsinglets <- merge(ovary.EMTAB12889.filt.list[[1]],
                                      y = c(ovary.EMTAB12889.filt.list[[2]], ovary.EMTAB12889.filt.list[[3]], ovary.EMTAB12889.filt.list[[4]],
                                            ovary.EMTAB12889.filt.list[[5]], ovary.EMTAB12889.filt.list[[6]], ovary.EMTAB12889.filt.list[[7]], ovary.EMTAB12889.filt.list[[8]], ovary.EMTAB12889.filt.list[[9]],
                                            ovary.EMTAB12889.filt.list[[10]], ovary.EMTAB12889.filt.list[[11]], ovary.EMTAB12889.filt.list[[12]]),
                                      project = "ovary.EMTAB12889")
ovary.EMTAB12889.DFsinglets
# An object of class Seurat 
# 42968 features across 9030 samples within 2 assays 
# Active assay: SCT (21477 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
ovary.EMTAB12889.DFsinglets@meta.data <- ovary.EMTAB12889.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.EMTAB12889.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.EMTAB12889.filt.list.scds <- lapply(ovary.EMTAB12889.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.EMTAB12889.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.EMTAB12889.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.EMTAB12889.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.EMTAB12889.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.EMTAB12889.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.EMTAB12889.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.EMTAB12889.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.EMTAB12889.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.EMTAB12889.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_AC", names(ovary.EMTAB12889.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.EMTAB12889.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.EMTAB12889.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.EMTAB12889.DFsinglets@meta.data[colnames(ovary.EMTAB12889.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.EMTAB12889.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.EMTAB12889.DFsinglets@meta.data$DoubletFinder, ovary.EMTAB12889.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet      55    1530
# Singlet     750   28800

ovary.EMTAB12889.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.EMTAB12889.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.EMTAB12889.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                             "Doublet", "Singlet")
table(ovary.EMTAB12889.DFsinglets@meta.data$DoubletCall)
# Doublet Singlet 
#    2335   28800 

# re-run dimensionality reduction for plotting purposes
ovary.EMTAB12889.DFsinglets <- SCTransform(object = ovary.EMTAB12889.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library", "decontX_contamination"))
ovary.EMTAB12889.DFsinglets <- RunPCA(ovary.EMTAB12889.DFsinglets, npcs = 50)
ovary.EMTAB12889.DFsinglets <- RunUMAP(ovary.EMTAB12889.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB12889_DoubletCall_UMAP.pdf"))
DimPlot(ovary.EMTAB12889.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.EMTAB12889.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB12889_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.EMTAB12889.DFsinglets   <- subset(ovary.EMTAB12889.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets
ovary.EMTAB12889.DFsinglets
# An object of class Seurat 
# 39297 features across 28800 samples within 2 assays 
# Active assay: SCT (19647 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(ovary.EMTAB12889.DFsinglets@meta.data$Library)
# Ovary_12M_1 Ovary_12M_2 Ovary_12M_3 Ovary_12M_4 Ovary_12M_5 Ovary_15M_1 Ovary_15M_3 Ovary_15M_4  Ovary_9M_1  Ovary_9M_2  Ovary_9M_3  Ovary_9M_4 
#        2189        2989        3825        1803        2279        2852        2853        2972         917         795        2942        2384 

# save filtered/annotated object
save(ovary.EMTAB12889.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB12889_Seurat_object_SINGLETS.RData"))

# ovary.EMTAB12889.DFsinglets <- SetIdent(ovary.EMTAB12889.DFsinglets, value = "seurat_clusters")
# FeaturePlot(ovary.EMTAB12889.DFsinglets, features = c("Ptprc", "Prox1"))
# VlnPlot(ovary.EMTAB12889.DFsinglets, feature = c("nCount_RNA"), pt.size = 0)
# VlnPlot(ovary.EMTAB12889.DFsinglets, feature = c("nFeature_RNA"), pt.size = 0)
# VlnPlot(ovary.EMTAB12889.DFsinglets, feature = c("decontX_contamination"), pt.size = 0)

t.cluster.vs.library <- as.matrix(table(ovary.EMTAB12889.DFsinglets@meta.data$seurat_clusters, ovary.EMTAB12889.DFsinglets@meta.data$Library))

# Filter clusters that are represented by a single library
# First, create a binary matrix
binary_matrix <- t.cluster.vs.library > 0

# Sum across rows to find the number of libraries contributing to each cluster
libraries_per_cluster <- rowSums(binary_matrix)

# Get clusters that have cells from only one library
single_library_clusters <- names(libraries_per_cluster)[libraries_per_cluster == 1]

print(single_library_clusters)     # None

# Filter Seurat object to exclude cluster "34" "35"

ovary.EMTAB12889.cl <- ovary.EMTAB12889.DFsinglets

ovary.EMTAB12889.cl
# An object of class Seurat 
# 39297 features across 28800 samples within 2 assays 
# Active assay: SCT (19647 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# save filtered/annotated object
save(ovary.EMTAB12889.cl, file = paste0(Sys.Date(),"_10x_ovary_Mouse_EMTAB12889_processed_clean_Seurat_object_FINAL.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Mouse_EMTAB12889_preprocessing_session_info.txt", sep =""))
sessionInfo()
sink()
