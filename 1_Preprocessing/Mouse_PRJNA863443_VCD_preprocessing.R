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

################################################################################
# Integrative ovarian aging analysis
# 10x ovary PRJNA863443 VCD data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.CTL_3m_30d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_30d_1/outs/filtered_feature_bc_matrix/")
counts.CTL_3m_30d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_30d_2/outs/filtered_feature_bc_matrix/")
counts.CTL_3m_90d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_90d_1/outs/filtered_feature_bc_matrix/")
counts.CTL_3m_90d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_90d_2/outs/filtered_feature_bc_matrix/")
counts.CTL_10m_30d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_30d_1/outs/filtered_feature_bc_matrix/")
counts.CTL_10m_30d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_30d_2/outs/filtered_feature_bc_matrix/")
counts.CTL_10m_90d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_90d_1/outs/filtered_feature_bc_matrix/")
counts.CTL_10m_90d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_90d_2/outs/filtered_feature_bc_matrix/")

counts.raw.CTL_3m_30d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_30d_1/outs/raw_feature_bc_matrix/")
counts.raw.CTL_3m_30d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_30d_2/outs/raw_feature_bc_matrix/")
counts.raw.CTL_3m_90d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_90d_1/outs/raw_feature_bc_matrix/")
counts.raw.CTL_3m_90d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_3m_90d_2/outs/raw_feature_bc_matrix/")
counts.raw.CTL_10m_30d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_30d_1/outs/raw_feature_bc_matrix/")
counts.raw.CTL_10m_30d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_30d_2/outs/raw_feature_bc_matrix/")
counts.raw.CTL_10m_90d_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_90d_1/outs/raw_feature_bc_matrix/")
counts.raw.CTL_10m_90d_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/VCD/CTL_10m_90d_2/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce.CTL_3m_30d_1 <- SingleCellExperiment(list(counts = counts.CTL_3m_30d_1))
sce.CTL_3m_30d_2 <- SingleCellExperiment(list(counts = counts.CTL_3m_30d_2))
sce.CTL_3m_90d_1 <- SingleCellExperiment(list(counts = counts.CTL_3m_90d_1))
sce.CTL_3m_90d_2 <- SingleCellExperiment(list(counts = counts.CTL_3m_90d_2))
sce.CTL_10m_30d_1 <- SingleCellExperiment(list(counts = counts.CTL_10m_30d_1))
sce.CTL_10m_30d_2 <- SingleCellExperiment(list(counts = counts.CTL_10m_30d_2))
sce.CTL_10m_90d_1 <- SingleCellExperiment(list(counts = counts.CTL_10m_90d_1))
sce.CTL_10m_90d_2 <- SingleCellExperiment(list(counts = counts.CTL_10m_90d_2))

sce.raw.CTL_3m_30d_1 <- SingleCellExperiment(list(counts = counts.raw.CTL_3m_30d_1))
sce.raw.CTL_3m_30d_2 <- SingleCellExperiment(list(counts = counts.raw.CTL_3m_30d_2))
sce.raw.CTL_3m_90d_1 <- SingleCellExperiment(list(counts = counts.raw.CTL_3m_90d_1))
sce.raw.CTL_3m_90d_2 <- SingleCellExperiment(list(counts = counts.raw.CTL_3m_90d_2))
sce.raw.CTL_10m_30d_1 <- SingleCellExperiment(list(counts = counts.raw.CTL_10m_30d_1))
sce.raw.CTL_10m_30d_2 <- SingleCellExperiment(list(counts = counts.raw.CTL_10m_30d_2))
sce.raw.CTL_10m_90d_1 <- SingleCellExperiment(list(counts = counts.raw.CTL_10m_90d_1))
sce.raw.CTL_10m_90d_2 <- SingleCellExperiment(list(counts = counts.raw.CTL_10m_90d_2))

sce.CTL_3m_30d_1 <- decontX(sce.CTL_3m_30d_1, background = sce.raw.CTL_3m_30d_1)
sce.CTL_3m_30d_2 <- decontX(sce.CTL_3m_30d_2, background = sce.raw.CTL_3m_30d_2)
sce.CTL_3m_90d_1 <- decontX(sce.CTL_3m_90d_1, background = sce.raw.CTL_3m_90d_1)
sce.CTL_3m_90d_2 <- decontX(sce.CTL_3m_90d_2, background = sce.raw.CTL_3m_90d_2)
sce.CTL_10m_30d_1 <- decontX(sce.CTL_10m_30d_1, background = sce.raw.CTL_10m_30d_1)
sce.CTL_10m_30d_2 <- decontX(sce.CTL_10m_30d_2, background = sce.raw.CTL_10m_30d_2)
sce.CTL_10m_90d_1 <- decontX(sce.CTL_10m_90d_1, background = sce.raw.CTL_10m_90d_1)
sce.CTL_10m_90d_2 <- decontX(sce.CTL_10m_90d_2, background = sce.raw.CTL_10m_90d_2)

# Save DecontX result
save(sce.CTL_3m_30d_1, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_3m_30d_1.RData"))
save(sce.CTL_3m_30d_2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_3m_30d_2.RData"))
save(sce.CTL_3m_90d_1, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_3m_90d_1.RData"))
save(sce.CTL_3m_90d_2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_3m_90d_2.RData"))
save(sce.CTL_10m_30d_1, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_10m_30d_1.RData"))
save(sce.CTL_10m_30d_2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_10m_30d_2.RData"))
save(sce.CTL_10m_90d_1, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_10m_90d_1.RData"))
save(sce.CTL_10m_90d_2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DecontX_SCE_object_sce.CTL_10m_90d_2.RData"))

# Plot UMAP
UMAP.CTL_3m_30d_1 <- reducedDim(sce.CTL_3m_30d_1, "decontX_UMAP")
UMAP.CTL_3m_30d_2 <- reducedDim(sce.CTL_3m_30d_2, "decontX_UMAP")
UMAP.CTL_3m_90d_1 <- reducedDim(sce.CTL_3m_90d_1, "decontX_UMAP")
UMAP.CTL_3m_90d_2 <- reducedDim(sce.CTL_3m_90d_2, "decontX_UMAP")
UMAP.CTL_10m_30d_1 <- reducedDim(sce.CTL_10m_30d_1, "decontX_UMAP")
UMAP.CTL_10m_30d_2 <- reducedDim(sce.CTL_10m_30d_2, "decontX_UMAP")
UMAP.CTL_10m_90d_1 <- reducedDim(sce.CTL_10m_90d_1, "decontX_UMAP")
UMAP.CTL_10m_90d_2 <- reducedDim(sce.CTL_10m_90d_2, "decontX_UMAP")

plotDecontXContamination(sce.CTL_3m_30d_1)
plotDecontXContamination(sce.CTL_3m_30d_2)
plotDecontXContamination(sce.CTL_3m_90d_1)
plotDecontXContamination(sce.CTL_3m_90d_2)
plotDecontXContamination(sce.CTL_10m_30d_1)
plotDecontXContamination(sce.CTL_10m_30d_2)
plotDecontXContamination(sce.CTL_10m_90d_1)
plotDecontXContamination(sce.CTL_10m_90d_2)

# Create Seurat objects
seurat.CTL_3m_30d_1 <- CreateSeuratObject(round(decontXcounts(sce.CTL_3m_30d_1)), meta.data=as.data.frame(colData(sce.CTL_3m_30d_1)))
seurat.CTL_3m_30d_2 <- CreateSeuratObject(round(decontXcounts(sce.CTL_3m_30d_2)), meta.data=as.data.frame(colData(sce.CTL_3m_30d_2)))
seurat.CTL_3m_90d_1 <- CreateSeuratObject(round(decontXcounts(sce.CTL_3m_90d_1)), meta.data=as.data.frame(colData(sce.CTL_3m_90d_1)))
seurat.CTL_3m_90d_2 <- CreateSeuratObject(round(decontXcounts(sce.CTL_3m_90d_2)), meta.data=as.data.frame(colData(sce.CTL_3m_90d_2)))
seurat.CTL_10m_30d_1 <- CreateSeuratObject(round(decontXcounts(sce.CTL_10m_30d_1)), meta.data=as.data.frame(colData(sce.CTL_10m_30d_1)))
seurat.CTL_10m_30d_2 <- CreateSeuratObject(round(decontXcounts(sce.CTL_10m_30d_2)), meta.data=as.data.frame(colData(sce.CTL_10m_30d_2)))
seurat.CTL_10m_90d_1 <- CreateSeuratObject(round(decontXcounts(sce.CTL_10m_90d_1)), meta.data=as.data.frame(colData(sce.CTL_10m_90d_1)))
seurat.CTL_10m_90d_2 <- CreateSeuratObject(round(decontXcounts(sce.CTL_10m_90d_2)), meta.data=as.data.frame(colData(sce.CTL_10m_90d_2)))

ovary.VCD <- merge(seurat.CTL_3m_30d_1, 
                   y =  c(seurat.CTL_3m_30d_2,
                          seurat.CTL_3m_90d_1,
                          seurat.CTL_3m_90d_2,
                          seurat.CTL_10m_30d_1,
                          seurat.CTL_10m_30d_2,
                          seurat.CTL_10m_90d_1,
                          seurat.CTL_10m_90d_2), 
                   add.cell.ids = c("CTL_3m_30d_1",
                                    "CTL_3m_30d_2",
                                    "CTL_3m_90d_1",
                                    "CTL_3m_90d_2",
                                    "CTL_10m_30d_1",
                                    "CTL_10m_30d_2",
                                    "CTL_10m_90d_1",
                                    "CTL_10m_90d_2"), 
                   project = "10x_ovary_VCD")

ovary.VCD
# An object of class Seurat 
# 32285 features across 122038 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

save(ovary.VCD, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_post-DecontX_Seurat_object.RData"))

# Remove intermediate files
rm(counts.CTL_3m_30d_1, counts.CTL_3m_30d_2, counts.CTL_3m_90d_1, counts.CTL_3m_90d_2,
   counts.CTL_10m_30d_1, counts.CTL_10m_30d_2, counts.CTL_10m_90d_1, counts.CTL_10m_90d_2,
   counts.VCD_3m_30d_1, counts.VCD_3m_30d_2, counts.VCD_3m_90d_1, counts.VCD_3m_90d_2,
   counts.VCD_10m_30d_1, counts.VCD_10m_30d_2, counts.VCD_10m_90d_1, counts.VCD_10m_90d_2,
   counts.raw.CTL_3m_30d_1, counts.raw.CTL_3m_30d_2, counts.raw.CTL_3m_90d_1, counts.raw.CTL_3m_90d_2,
   counts.raw.CTL_10m_30d_1, counts.raw.CTL_10m_30d_2, counts.raw.CTL_10m_90d_1, counts.raw.CTL_10m_90d_2,
   counts.raw.VCD_3m_30d_1, counts.raw.VCD_3m_30d_2, counts.raw.VCD_3m_90d_1, counts.raw.VCD_3m_90d_2,
   counts.raw.VCD_10m_30d_1, counts.raw.VCD_10m_30d_2, counts.raw.VCD_10m_90d_1, counts.raw.VCD_10m_90d_2,
   sce.CTL_3m_30d_1, sce.CTL_3m_30d_2, sce.CTL_3m_90d_1, sce.CTL_3m_90d_2,
   sce.CTL_10m_30d_1, sce.CTL_10m_30d_2, sce.CTL_10m_90d_1, sce.CTL_10m_90d_2,
   sce.VCD_3m_30d_1, sce.VCD_3m_30d_2, sce.VCD_3m_90d_1, sce.VCD_3m_90d_2,
   sce.VCD_10m_30d_1, sce.VCD_10m_30d_2, sce.VCD_10m_90d_1, sce.VCD_10m_90d_2,
   sce.raw.CTL_3m_30d_1, sce.raw.CTL_3m_30d_2, sce.raw.CTL_3m_90d_1, sce.raw.CTL_3m_90d_2,
   sce.raw.CTL_10m_30d_1, sce.raw.CTL_10m_30d_2, sce.raw.CTL_10m_90d_1, sce.raw.CTL_10m_90d_2,
   sce.raw.VCD_3m_30d_1, sce.raw.VCD_3m_30d_2, sce.raw.VCD_3m_90d_1, sce.raw.VCD_3m_90d_2,
   sce.raw.VCD_10m_30d_1, sce.raw.VCD_10m_30d_2, sce.raw.VCD_10m_90d_1, sce.raw.VCD_10m_90d_2,
   seurat.CTL_3m_30d_1, seurat.CTL_3m_30d_2, seurat.CTL_3m_90d_1, seurat.CTL_3m_90d_2,
   seurat.CTL_10m_30d_1, seurat.CTL_10m_30d_2, seurat.CTL_10m_90d_1, seurat.CTL_10m_90d_2,
   seurat.VCD_3m_30d_1, seurat.VCD_3m_30d_2, seurat.VCD_3m_90d_1, seurat.VCD_3m_90d_2,
   seurat.VCD_10m_30d_1, seurat.VCD_10m_30d_2, seurat.VCD_10m_90d_1, seurat.VCD_10m_90d_2)

################################################################################
# 2. Add metadata to Seurat object
################################################################################

# create Library label
Library <- rep("NA", length(colnames(ovary.VCD@assays$RNA)))
Library[grep("CTL_3m_30d_1", colnames(ovary.VCD@assays$RNA))] <- "CTL_3m_30d_1"
Library[grep("CTL_3m_30d_2", colnames(ovary.VCD@assays$RNA))] <- "CTL_3m_30d_2"
Library[grep("CTL_3m_90d_1", colnames(ovary.VCD@assays$RNA))] <- "CTL_3m_90d_1"
Library[grep("CTL_3m_90d_2", colnames(ovary.VCD@assays$RNA))] <- "CTL_3m_90d_2"
Library[grep("CTL_10m_30d_1", colnames(ovary.VCD@assays$RNA))] <- "CTL_10m_30d_1"
Library[grep("CTL_10m_30d_2", colnames(ovary.VCD@assays$RNA))] <- "CTL_10m_30d_2"
Library[grep("CTL_10m_90d_1", colnames(ovary.VCD@assays$RNA))] <- "CTL_10m_90d_1"
Library[grep("CTL_10m_90d_2", colnames(ovary.VCD@assays$RNA))] <- "CTL_10m_90d_2"

Library <- data.frame(Library)
rownames(Library) <- colnames(ovary.VCD@assays$RNA)

# create Age label
Age <- rep("NA", length(colnames(ovary.VCD@assays$RNA)))
Age[grep("3m", colnames(ovary.VCD@assays$RNA))] <- "3m"
Age[grep("10m", colnames(ovary.VCD@assays$RNA))] <- "10m"

Age <- data.frame(Age)
rownames(Age) <- colnames(ovary.VCD@assays$RNA)

# create Duration label
Duration <- rep("NA", length(colnames(ovary.VCD@assays$RNA)))
Duration[grep("30d", colnames(ovary.VCD@assays$RNA))] <- "30d"
Duration[grep("90d", colnames(ovary.VCD@assays$RNA))] <- "90d"

Duration <- data.frame(Duration)
rownames(Duration) <- colnames(ovary.VCD@assays$RNA)

# create Batch label
Batch <- rep("NA", length(colnames(ovary.VCD@assays$RNA)))
Batch[grep("_1", colnames(ovary.VCD@assays$RNA))] <- "VCD_1"
Batch[grep("30d_2", colnames(ovary.VCD@assays$RNA))] <- "VCD_2"
Batch[grep("90d_2", colnames(ovary.VCD@assays$RNA))] <- "VCD_3"

Batch <- data.frame(Batch)
rownames(Batch) <- colnames(ovary.VCD@assays$RNA)

# update Seurat with metadata
ovary.VCD <- AddMetaData(object = ovary.VCD, metadata = as.vector(Library), col.name = "Library")
ovary.VCD <- AddMetaData(object = ovary.VCD, metadata = as.vector(Age), col.name = "Age")
ovary.VCD <- AddMetaData(object = ovary.VCD, metadata = as.vector(Duration), col.name = "Duration")
ovary.VCD <- AddMetaData(object = ovary.VCD, metadata = as.vector(Batch), col.name = "Batch")

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500
#                    500 < nCount_RNA < 1e5
#                    percent.mito < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.VCD <- SetIdent(ovary.VCD, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.VCD@assays$RNA)
num.cells <- Matrix::rowSums(ovary.VCD@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.VCD.filt <- subset(ovary.VCD, features = genes.use)
ovary.VCD.filt
# An object of class Seurat 
# 24489 features across 122038 samples within 1 assay 
# Active assay: RNA (24489 features, 0 variable features)

############### QC - mitochondrial genes ###############

ovary.VCD.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.VCD.filt, pattern = "^mt-")
# head(ovary.VCD.filt@meta.data)

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.VCD.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.VCD.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.VCD.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.VCD.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.VCD.filt <- subset(ovary.VCD.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)
ovary.VCD.filt

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.VCD.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.VCD.filt, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.VCD.filt <- NormalizeData(ovary.VCD.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.VCD.filt <- SCTransform(object = ovary.VCD.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library", "Batch"))

save(ovary.VCD.filt, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_Seurat_object_postSCT.RData"))

ovary.VCD.filt <- RunPCA(ovary.VCD.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_ElbowPlot.pdf"))
ElbowPlot(ovary.VCD.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.VCD.filt[["pca"]]@stdev / sum(ovary.VCD.filt[["pca"]]@stdev) * 100

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
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.VCD.filt <- RunUMAP(ovary.VCD.filt, dims = 1:pcs)
ovary.VCD.filt <- FindNeighbors(ovary.VCD.filt, dims = 1:pcs)
ovary.VCD.filt <- FindClusters(object = ovary.VCD.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.VCD.filt, reduction = "umap")
dev.off()

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
table(ovary.VCD.filt$Library)

# Clean memory of intermediates
rm(ovary.VCD)


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
ovary.VCD.filt.list <- SplitObject(ovary.VCD.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.VCD.filt.list, ncol))))/100
pred.dblt.rate
#   CTL_3m_30d_1  CTL_3m_30d_2  CTL_3m_90d_1  CTL_3m_90d_2 CTL_10m_30d_1 CTL_10m_30d_2 CTL_10m_90d_1 CTL_10m_90d_2  VCD_3m_30d_1  VCD_3m_30d_2  VCD_3m_90d_1  VCD_3m_90d_2 VCD_10m_30d_1 VCD_10m_30d_2 VCD_10m_90d_1 
#      0.0134288     0.0493416     0.0155320     0.0501688     0.0134816     0.0615296     0.0098648     0.0623744     0.0136928     0.0795872     0.0097152     0.0828432     0.0141328     0.0233992     0.0152328 
#  VCD_10m_90d_2 
#      0.0169752 

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.VCD.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.VCD <- paramSweep_v3(ovary.VCD.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.VCD    <- summarizeSweep(sweep.res.list.ovary.VCD, GT = FALSE)
  bcmvn.ovary.VCD          <- find.pK(sweep.stats.ovary.VCD)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.VCD <- as.numeric(as.character(bcmvn.ovary.VCD[as.numeric(bcmvn.ovary.VCD$pK[bcmvn.ovary.VCD$BCmetric == max(bcmvn.ovary.VCD$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.VCD.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.VCD.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.VCD.filt.list[[i]] <- doubletFinder_v3(ovary.VCD.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.VCD, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.VCD.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.VCD.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.VCD.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.VCD.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.VCD.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_VCD",names(ovary.VCD.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.VCD.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.VCD.DFsinglets <- merge(ovary.VCD.filt.list[[1]],
                              y = c(ovary.VCD.filt.list[[2]], ovary.VCD.filt.list[[3]], ovary.VCD.filt.list[[4]],
                                    ovary.VCD.filt.list[[5]], ovary.VCD.filt.list[[6]], ovary.VCD.filt.list[[7]], ovary.VCD.filt.list[[8]]),
                              project = "ovary.VCD")
ovary.VCD.DFsinglets
# An object of class Seurat 
# 48969 features across 60375 samples within 2 assays 
# Active assay: SCT (24484 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
ovary.VCD.DFsinglets@meta.data <- ovary.VCD.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.VCD.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.VCD.filt.list.scds <- lapply(ovary.VCD.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.VCD.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.VCD.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.VCD.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.VCD.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.VCD.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.VCD.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.VCD.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.VCD.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.VCD.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_VCD", names(ovary.VCD.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.VCD.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.VCD.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.VCD.DFsinglets@meta.data[colnames(ovary.VCD.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.VCD.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.VCD.DFsinglets@meta.data$DoubletFinder, ovary.VCD.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet     398    4303
# Singlet    2792   52882

ovary.VCD.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.VCD.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.VCD.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                     "Doublet", "Singlet")
table(ovary.VCD.DFsinglets@meta.data$DoubletCall)
# Doublet Singlet 
#    7493   52882 

# re-run dimensionality reduction for plotting purposes
ovary.VCD.DFsinglets <- SCTransform(object = ovary.VCD.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library", "Batch"))
ovary.VCD.DFsinglets <- RunPCA(ovary.VCD.DFsinglets, npcs = 50)
ovary.VCD.DFsinglets <- RunUMAP(ovary.VCD.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_DoubletCall_UMAP.pdf"))
DimPlot(ovary.VCD.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.VCD.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.VCD.DFsinglets   <- subset(ovary.VCD.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets
ovary.VCD.DFsinglets
# An object of class Seurat 
# 48969 features across 52882 samples within 2 assays 
# Active assay: SCT (24484 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(ovary.VCD.DFsinglets@meta.data$Library)
# CTL_10m_30d_1 CTL_10m_30d_2 CTL_10m_90d_1 CTL_10m_90d_2  CTL_3m_30d_1  CTL_3m_30d_2  CTL_3m_90d_1  CTL_3m_90d_2 VCD_10m_30d_1 VCD_10m_30d_2 VCD_10m_90d_1 VCD_10m_90d_2  VCD_3m_30d_1  VCD_3m_30d_2  VCD_3m_90d_1 
#          1454          6043          1071          6041          1449          4913          1668          5017          1524          2468          1635          1817          1475          7444          1055 
# VCD_3m_90d_2 
#         7808 

# save filtered/annotated object
save(ovary.VCD.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_Seurat_object_SINGLETS.RData"))

ovary.VCD.DFsinglets <- SetIdent(ovary.VCD.DFsinglets, value = "seurat_clusters")
FeaturePlot(ovary.VCD.DFsinglets, features = c("Ptprc", "Prox1"))
VlnPlot(ovary.VCD.DFsinglets, feature = c("nCount_RNA"), pt.size = 0)

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Mouse_PRJNA863443_VCD_preprocessing_session_info.txt", sep =""))
sessionInfo()
sink()