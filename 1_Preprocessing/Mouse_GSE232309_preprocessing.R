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
# Mouse 10x GSE232309 data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.3M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M1/outs/filtered_feature_bc_matrix/")
counts.3M_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M2/outs/filtered_feature_bc_matrix/")
counts.3M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M3/outs/filtered_feature_bc_matrix/")
counts.3M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M4/outs/filtered_feature_bc_matrix/")
counts.9M_5 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M5/outs/filtered_feature_bc_matrix/")
counts.9M_6 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M6/outs/filtered_feature_bc_matrix/")
counts.9M_7 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M7/outs/filtered_feature_bc_matrix/")
counts.9M_8 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M8/outs/filtered_feature_bc_matrix/")

counts.raw.3M_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M1/outs/raw_feature_bc_matrix/")
counts.raw.3M_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M2/outs/raw_feature_bc_matrix/")
counts.raw.3M_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M3/outs/raw_feature_bc_matrix/")
counts.raw.3M_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/3M4/outs/raw_feature_bc_matrix/")
counts.raw.9M_5 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M5/outs/raw_feature_bc_matrix/")
counts.raw.9M_6 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M6/outs/raw_feature_bc_matrix/")
counts.raw.9M_7 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M7/outs/raw_feature_bc_matrix/")
counts.raw.9M_8 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Mouse_GSE232309/9M8/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce.3M_1 <- SingleCellExperiment(list(counts = counts.3M_1))
sce.3M_2 <- SingleCellExperiment(list(counts = counts.3M_2))
sce.3M_3 <- SingleCellExperiment(list(counts = counts.3M_3))
sce.3M_4 <- SingleCellExperiment(list(counts = counts.3M_4))
sce.9M_5 <- SingleCellExperiment(list(counts = counts.9M_5))
sce.9M_6 <- SingleCellExperiment(list(counts = counts.9M_6))
sce.9M_7 <- SingleCellExperiment(list(counts = counts.9M_7))
sce.9M_8 <- SingleCellExperiment(list(counts = counts.9M_8))

sce.raw.3M_1 <- SingleCellExperiment(list(counts = counts.raw.3M_1))
sce.raw.3M_2 <- SingleCellExperiment(list(counts = counts.raw.3M_2))
sce.raw.3M_3 <- SingleCellExperiment(list(counts = counts.raw.3M_3))
sce.raw.3M_4 <- SingleCellExperiment(list(counts = counts.raw.3M_4))
sce.raw.9M_5 <- SingleCellExperiment(list(counts = counts.raw.9M_5))
sce.raw.9M_6 <- SingleCellExperiment(list(counts = counts.raw.9M_6))
sce.raw.9M_7 <- SingleCellExperiment(list(counts = counts.raw.9M_7))
sce.raw.9M_8 <- SingleCellExperiment(list(counts = counts.raw.9M_8))

sce.3M_1 <- decontX(sce.3M_1, background = sce.raw.3M_1)
sce.3M_2 <- decontX(sce.3M_2, background = sce.raw.3M_2)
sce.3M_3 <- decontX(sce.3M_3, background = sce.raw.3M_3)
sce.3M_4 <- decontX(sce.3M_4, background = sce.raw.3M_4)
sce.9M_5 <- decontX(sce.9M_5, background = sce.raw.9M_5)
sce.9M_6 <- decontX(sce.9M_6, background = sce.raw.9M_6)
sce.9M_7 <- decontX(sce.9M_7, background = sce.raw.9M_7)
sce.9M_8 <- decontX(sce.9M_8, background = sce.raw.9M_8)

# Save DecontX result
save(sce.3M_1, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.3M_1.RData"))
save(sce.3M_2, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.3M_2.RData"))
save(sce.3M_3, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.3M_3.RData"))
save(sce.3M_4, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.3M_4.RData"))
save(sce.9M_5, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.9M_5.RData"))
save(sce.9M_6, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.9M_6.RData"))
save(sce.9M_7, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.9M_7.RData"))
save(sce.9M_8, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_DecontX_SCE_object_sce.9M_8.RData"))

# Plot UMAP
UMAP.3M_1 <- reducedDim(sce.3M_1, "decontX_UMAP")
UMAP.3M_2 <- reducedDim(sce.3M_2, "decontX_UMAP")
UMAP.3M_3 <- reducedDim(sce.3M_3, "decontX_UMAP")
UMAP.3M_4 <- reducedDim(sce.3M_4, "decontX_UMAP")
UMAP.9M_5 <- reducedDim(sce.9M_5, "decontX_UMAP")
UMAP.9M_6 <- reducedDim(sce.9M_6, "decontX_UMAP")
UMAP.9M_7 <- reducedDim(sce.9M_7, "decontX_UMAP")
UMAP.9M_8 <- reducedDim(sce.9M_8, "decontX_UMAP")

plotDecontXContamination(sce.3M_1)
plotDecontXContamination(sce.3M_2)
plotDecontXContamination(sce.3M_3)
plotDecontXContamination(sce.3M_4)
plotDecontXContamination(sce.9M_5)
plotDecontXContamination(sce.9M_6)
plotDecontXContamination(sce.9M_7)
plotDecontXContamination(sce.9M_8)

# Create Seurat objects
seurat.3M_1 <- CreateSeuratObject(round(decontXcounts(sce.3M_1)), meta.data=as.data.frame(colData(sce.3M_1)))
seurat.3M_2 <- CreateSeuratObject(round(decontXcounts(sce.3M_2)), meta.data=as.data.frame(colData(sce.3M_2)))
seurat.3M_3 <- CreateSeuratObject(round(decontXcounts(sce.3M_3)), meta.data=as.data.frame(colData(sce.3M_3)))
seurat.3M_4 <- CreateSeuratObject(round(decontXcounts(sce.3M_4)), meta.data=as.data.frame(colData(sce.3M_4)))
seurat.9M_5 <- CreateSeuratObject(round(decontXcounts(sce.9M_5)), meta.data=as.data.frame(colData(sce.9M_5)))
seurat.9M_6 <- CreateSeuratObject(round(decontXcounts(sce.9M_6)), meta.data=as.data.frame(colData(sce.9M_6)))
seurat.9M_7 <- CreateSeuratObject(round(decontXcounts(sce.9M_7)), meta.data=as.data.frame(colData(sce.9M_7)))
seurat.9M_8 <- CreateSeuratObject(round(decontXcounts(sce.9M_8)), meta.data=as.data.frame(colData(sce.9M_8)))


ovary.GSE232309 <- merge(seurat.3M_1, 
                         y =  c(seurat.3M_2,
                                seurat.3M_3,
                                seurat.3M_4,
                                seurat.9M_5,
                                seurat.9M_6,
                                seurat.9M_7,
                                seurat.9M_8), 
                   add.cell.ids = c("3M_1",
                                    "3M_2",
                                    "3M_3",
                                    "3M_4",
                                    "9M_5",
                                    "9M_6",
                                    "9M_7",
                                    "9M_8"), 
                   project = "10x_ovary_GSE232309")

ovary.GSE232309
# An object of class Seurat 
# 32285 features across 17178 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

save(ovary.GSE232309, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_post-DecontX_Seurat_object.RData"))

# Remove intermediate files
rm(counts.3M_1, counts.3M_2, counts.3M_3, counts.3M_4,
   counts.9M_5, counts.9M_6, counts.9M_7, counts.9M_8,
   counts.raw.3M_1, counts.raw.3M_2, counts.raw.3M_3, counts.raw.3M_4,
   counts.raw.9M_5, counts.raw.9M_6, counts.raw.9M_7, counts.raw.9M_8,
   sce.3M_1, sce.3M_2, sce.3M_3, sce.3M_4,
   sce.9M_5, sce.9M_6, sce.9M_7, sce.9M_8,
   sce.raw.3M_1, sce.raw.3M_2, sce.raw.3M_3, sce.raw.3M_4,
   sce.raw.9M_5, sce.raw.9M_6, sce.raw.9M_7, sce.raw.9M_8,
   seurat.3M_1, seurat.3M_2, seurat.3M_3, seurat.3M_4,
   seurat.9M_5, seurat.9M_6, seurat.9M_7, seurat.9M_8)

################################################################################
# 2. Add metadata to Seurat object
################################################################################

# create Library label
Library <- rep("NA", length(colnames(ovary.GSE232309@assays$RNA)))
Library[grep("3M_1", colnames(ovary.GSE232309@assays$RNA))] <- "3M_1"
Library[grep("3M_2", colnames(ovary.GSE232309@assays$RNA))] <- "3M_2"
Library[grep("3M_3", colnames(ovary.GSE232309@assays$RNA))] <- "3M_3"
Library[grep("3M_4", colnames(ovary.GSE232309@assays$RNA))] <- "3M_4"
Library[grep("9M_5", colnames(ovary.GSE232309@assays$RNA))] <- "9M_5"
Library[grep("9M_6", colnames(ovary.GSE232309@assays$RNA))] <- "9M_6"
Library[grep("9M_7", colnames(ovary.GSE232309@assays$RNA))] <- "9M_7"
Library[grep("9M_8", colnames(ovary.GSE232309@assays$RNA))] <- "9M_8"

Library <- data.frame(Library)
rownames(Library) <- colnames(ovary.GSE232309@assays$RNA)

# create Age label
Age <- rep("NA", length(colnames(ovary.GSE232309@assays$RNA)))
Age[grep("3M", colnames(ovary.GSE232309@assays$RNA))] <- "3m"
Age[grep("9M", colnames(ovary.GSE232309@assays$RNA))] <- "9m"

Age <- data.frame(Age)
rownames(Age) <- colnames(ovary.GSE232309@assays$RNA)

# update Seurat with metadata
ovary.GSE232309 <- AddMetaData(object = ovary.GSE232309, metadata = as.vector(Library), col.name = "Library")
ovary.GSE232309 <- AddMetaData(object = ovary.GSE232309, metadata = as.vector(Age), col.name = "Age")

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500
#                    500 < nCount_RNA < 1e5
#                    percent.mito < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.GSE232309 <- SetIdent(ovary.GSE232309, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.GSE232309@assays$RNA)
num.cells <- Matrix::rowSums(ovary.GSE232309@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.GSE232309.filt <- subset(ovary.GSE232309, features = genes.use)
ovary.GSE232309.filt
# An object of class Seurat 
# 20665 features across 17178 samples within 1 assay 
# Active assay: RNA (20665 features, 0 variable features)

############### QC - mitochondrial genes ###############

ovary.GSE232309.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.GSE232309.filt, pattern = "^mt-")
# head(ovary.GSE232309.filt@meta.data)

pdf(paste(Sys.Date(),"10x_ovary_GSE232309_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.GSE232309.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.GSE232309.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.GSE232309.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.GSE232309.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_GSE232309_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.GSE232309.filt <- subset(ovary.GSE232309.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)
ovary.GSE232309.filt
# An object of class Seurat 
# 20665 features across 6201 samples within 1 assay 
# Active assay: RNA (20665 features, 0 variable features)

pdf(paste(Sys.Date(),"10x_ovary_GSE232309_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.GSE232309.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.GSE232309.filt, file = paste(Sys.Date(),"10x_ovary_GSE232309_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.GSE232309.filt <- NormalizeData(ovary.GSE232309.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.GSE232309.filt <- SCTransform(object = ovary.GSE232309.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library"))

save(ovary.GSE232309.filt, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_Seurat_object_postSCT.RData"))

ovary.GSE232309.filt <- RunPCA(ovary.GSE232309.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_GSE232309_ElbowPlot.pdf"))
ElbowPlot(ovary.GSE232309.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.GSE232309.filt[["pca"]]@stdev / sum(ovary.GSE232309.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 18

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 18

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_GSE232309_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.GSE232309.filt <- RunUMAP(ovary.GSE232309.filt, dims = 1:pcs)
ovary.GSE232309.filt <- FindNeighbors(ovary.GSE232309.filt, dims = 1:pcs)
ovary.GSE232309.filt <- FindClusters(object = ovary.GSE232309.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_GSE232309_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.GSE232309.filt, reduction = "umap")
dev.off()

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
table(ovary.GSE232309.filt$Library)
# 3M_1 3M_2 3M_3 3M_4 9M_5 9M_6 9M_7 9M_8 
#  519  425  206  191 1073  866 1769 1152 

# Clean memory of intermediates
rm(ovary.GSE232309)


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
ovary.GSE232309.filt.list <- SplitObject(ovary.GSE232309.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.GSE232309.filt.list, ncol))))/100
pred.dblt.rate
#      3M_1      3M_2      3M_3      3M_4      9M_5      9M_6      9M_7      9M_8 
# 0.0045672 0.0037400 0.0018128 0.0016808 0.0094424 0.0076208 0.0155672 0.0101376 

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.GSE232309.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary <- paramSweep_v3(ovary.GSE232309.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary    <- summarizeSweep(sweep.res.list.ovary, GT = FALSE)
  bcmvn.ovary          <- find.pK(sweep.stats.ovary)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary <- as.numeric(as.character(bcmvn.ovary[as.numeric(bcmvn.ovary$pK[bcmvn.ovary$BCmetric == max(bcmvn.ovary$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.GSE232309.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.GSE232309.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.GSE232309.filt.list[[i]] <- doubletFinder_v3(ovary.GSE232309.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.GSE232309.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.GSE232309.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.GSE232309.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.GSE232309.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.GSE232309.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_GSE232309",names(ovary.GSE232309.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.GSE232309.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.GSE232309.DFsinglets <- merge(ovary.GSE232309.filt.list[[1]],
                                    y = c(ovary.GSE232309.filt.list[[2]], ovary.GSE232309.filt.list[[3]], ovary.GSE232309.filt.list[[4]],
                                          ovary.GSE232309.filt.list[[5]], ovary.GSE232309.filt.list[[6]], ovary.GSE232309.filt.list[[7]], ovary.GSE232309.filt.list[[8]]),
                                    project = "ovary.GSE232309")
ovary.GSE232309.DFsinglets
# An object of class Seurat 
# 41317 features across 6201 samples within 2 assays 
# Active assay: SCT (20652 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
ovary.GSE232309.DFsinglets@meta.data <- ovary.GSE232309.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.GSE232309.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.GSE232309.filt.list.scds <- lapply(ovary.GSE232309.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.GSE232309.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.GSE232309.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.GSE232309.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.GSE232309.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.GSE232309.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.GSE232309.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.GSE232309.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.GSE232309.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.GSE232309.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_GSE232309", names(ovary.GSE232309.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.GSE232309.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.GSE232309.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.GSE232309.DFsinglets@meta.data[colnames(ovary.GSE232309.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.GSE232309.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.GSE232309.DFsinglets@meta.data$DoubletFinder, ovary.GSE232309.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet       6     209
# Singlet      57    5929

ovary.GSE232309.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.GSE232309.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.GSE232309.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                           "Doublet", "Singlet")
table(ovary.GSE232309.DFsinglets@meta.data$DoubletCall)
# Doublet Singlet 
#     272    5929 

# re-run dimensionality reduction for plotting purposes
ovary.GSE232309.DFsinglets <- SCTransform(object = ovary.GSE232309.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library"))
ovary.GSE232309.DFsinglets <- RunPCA(ovary.GSE232309.DFsinglets, npcs = 50)
ovary.GSE232309.DFsinglets <- RunUMAP(ovary.GSE232309.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_GSE232309_DoubletCall_UMAP.pdf"))
DimPlot(ovary.GSE232309.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.GSE232309.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.GSE232309.DFsinglets   <- subset(ovary.GSE232309.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets
ovary.GSE232309.DFsinglets
# An object of class Seurat 
# 41317 features across 5929 samples within 2 assays 
# Active assay: SCT (20652 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(ovary.GSE232309.DFsinglets@meta.data$Library)
# 3M_1 3M_2 3M_3 3M_4 9M_5 9M_6 9M_7 9M_8 
#  502  411  199  185 1026  831 1673 1102 

# save filtered/annotated object
save(ovary.GSE232309.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_GSE232309_Seurat_object_SINGLETS.RData"))

ovary.GSE232309.DFsinglets <- SetIdent(ovary.GSE232309.DFsinglets, value = "seurat_clusters")
FeaturePlot(ovary.GSE232309.DFsinglets, features = c("Ptprc", "Prox1"))
VlnPlot(ovary.GSE232309.DFsinglets, feature = c("nCount_RNA"), pt.size = 0)

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Mouse_GSE232309_preprocessing_session_info.txt", sep =""))
sessionInfo()
sink()
