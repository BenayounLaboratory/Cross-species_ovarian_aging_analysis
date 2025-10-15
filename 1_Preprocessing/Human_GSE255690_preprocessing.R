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
# 10x ovary Human GSE255690 data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.Young_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Young_1/outs/filtered_feature_bc_matrix/")
counts.Young_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Young_3/outs/filtered_feature_bc_matrix/")
counts.Young_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Young_4/outs/filtered_feature_bc_matrix/")
counts.Middle_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Middle_2/outs/filtered_feature_bc_matrix/")
counts.Middle_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Middle_3/outs/filtered_feature_bc_matrix/")
counts.Middle_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Middle_4/outs/filtered_feature_bc_matrix/")
counts.Old_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Old_1/outs/filtered_feature_bc_matrix/")
counts.Old_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Old_2/outs/filtered_feature_bc_matrix/")
counts.Old_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Old_3/outs/filtered_feature_bc_matrix/")

counts.raw.Young_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Young_1/outs/filtered_feature_bc_matrix/")
counts.raw.Young_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Young_3/outs/filtered_feature_bc_matrix/")
counts.raw.Young_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Young_4/outs/filtered_feature_bc_matrix/")
counts.raw.Middle_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Middle_2/outs/filtered_feature_bc_matrix/")
counts.raw.Middle_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Middle_3/outs/filtered_feature_bc_matrix/")
counts.raw.Middle_4 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Middle_4/outs/filtered_feature_bc_matrix/")
counts.raw.Old_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Old_1/outs/filtered_feature_bc_matrix/")
counts.raw.Old_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Old_2/outs/filtered_feature_bc_matrix/")
counts.raw.Old_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE255690/2_cellranger/Old_3/outs/filtered_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce.Young_1 <- SingleCellExperiment(list(counts = counts.Young_1))
sce.Young_3 <- SingleCellExperiment(list(counts = counts.Young_3))
sce.Young_4 <- SingleCellExperiment(list(counts = counts.Young_4))
sce.Middle_2 <- SingleCellExperiment(list(counts = counts.Middle_2))
sce.Middle_3 <- SingleCellExperiment(list(counts = counts.Middle_3))
sce.Middle_4 <- SingleCellExperiment(list(counts = counts.Middle_4))
sce.Old_1 <- SingleCellExperiment(list(counts = counts.Old_1))
sce.Old_2 <- SingleCellExperiment(list(counts = counts.Old_2))
sce.Old_3 <- SingleCellExperiment(list(counts = counts.Old_3))

sce.raw.Young_1 <- SingleCellExperiment(list(counts = counts.raw.Young_1))
sce.raw.Young_3 <- SingleCellExperiment(list(counts = counts.raw.Young_3))
sce.raw.Young_4 <- SingleCellExperiment(list(counts = counts.raw.Young_4))
sce.raw.Middle_2 <- SingleCellExperiment(list(counts = counts.raw.Middle_2))
sce.raw.Middle_3 <- SingleCellExperiment(list(counts = counts.raw.Middle_3))
sce.raw.Middle_4 <- SingleCellExperiment(list(counts = counts.raw.Middle_4))
sce.raw.Old_1 <- SingleCellExperiment(list(counts = counts.raw.Old_1))
sce.raw.Old_2 <- SingleCellExperiment(list(counts = counts.raw.Old_2))
sce.raw.Old_3 <- SingleCellExperiment(list(counts = counts.raw.Old_3))

sce.Young_1 <- decontX(sce.Young_1, background = sce.raw.Young_1)
sce.Young_3 <- decontX(sce.Young_3, background = sce.raw.Young_3)
sce.Young_4 <- decontX(sce.Young_4, background = sce.raw.Young_4)
sce.Middle_2 <- decontX(sce.Middle_2, background = sce.raw.Middle_2)
sce.Middle_3 <- decontX(sce.Middle_3, background = sce.raw.Middle_3)
sce.Middle_4 <- decontX(sce.Middle_4, background = sce.raw.Middle_4)
sce.Old_1 <- decontX(sce.Old_1, background = sce.raw.Old_1)
sce.Old_2 <- decontX(sce.Old_2, background = sce.raw.Old_2)
sce.Old_3 <- decontX(sce.Old_3, background = sce.raw.Old_3)

# Save DecontX result
save(sce.Young_1, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Young_1.Rdata"))
save(sce.Young_3, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Young_3.Rdata"))
save(sce.Young_4, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Young_4.Rdata"))
save(sce.Middle_2, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Middle_2.Rdata"))
save(sce.Middle_3, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Middle_3.Rdata"))
save(sce.Middle_4, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Middle_4.Rdata"))
save(sce.Old_1, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Old_1.Rdata"))
save(sce.Old_2, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Old_2.Rdata"))
save(sce.Old_3, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_DecontX_SCE_object_sce.Old_3.Rdata"))

# Plot UMAP
UMAP.Young_1 <- reducedDim(sce.Young_1, "decontX_UMAP")
UMAP.Young_3 <- reducedDim(sce.Young_3, "decontX_UMAP")
UMAP.Young_4 <- reducedDim(sce.Young_4, "decontX_UMAP")
UMAP.Middle_2 <- reducedDim(sce.Middle_2, "decontX_UMAP")
UMAP.Middle_3 <- reducedDim(sce.Middle_3, "decontX_UMAP")
UMAP.Middle_4 <- reducedDim(sce.Middle_4, "decontX_UMAP")
UMAP.Old_1 <- reducedDim(sce.Old_1, "decontX_UMAP")
UMAP.Old_2 <- reducedDim(sce.Old_2, "decontX_UMAP")
UMAP.Old_3 <- reducedDim(sce.Old_3, "decontX_UMAP")

plotDecontXContamination(sce.Young_1)
plotDecontXContamination(sce.Young_3)
plotDecontXContamination(sce.Young_4)
plotDecontXContamination(sce.Middle_2)
plotDecontXContamination(sce.Middle_3)
plotDecontXContamination(sce.Middle_4)
plotDecontXContamination(sce.Old_1)
plotDecontXContamination(sce.Old_2)
plotDecontXContamination(sce.Old_3)

# Create Seurat objects
seurat.Young_1 <- CreateSeuratObject(round(decontXcounts(sce.Young_1)), meta.data=as.data.frame(colData(sce.Young_1)))
seurat.Young_3 <- CreateSeuratObject(round(decontXcounts(sce.Young_3)), meta.data=as.data.frame(colData(sce.Young_3)))
seurat.Young_4 <- CreateSeuratObject(round(decontXcounts(sce.Young_4)), meta.data=as.data.frame(colData(sce.Young_4)))
seurat.Middle_2 <- CreateSeuratObject(round(decontXcounts(sce.Middle_2)), meta.data=as.data.frame(colData(sce.Middle_2)))
seurat.Middle_3 <- CreateSeuratObject(round(decontXcounts(sce.Middle_3)), meta.data=as.data.frame(colData(sce.Middle_3)))
seurat.Middle_4 <- CreateSeuratObject(round(decontXcounts(sce.Middle_4)), meta.data=as.data.frame(colData(sce.Middle_4)))
seurat.Old_1 <- CreateSeuratObject(round(decontXcounts(sce.Old_1)), meta.data=as.data.frame(colData(sce.Old_1)))
seurat.Old_2 <- CreateSeuratObject(round(decontXcounts(sce.Old_2)), meta.data=as.data.frame(colData(sce.Old_2)))
seurat.Old_3 <- CreateSeuratObject(round(decontXcounts(sce.Old_3)), meta.data=as.data.frame(colData(sce.Old_3)))

seurat.Young_1@meta.data$Library <- "Young_1"
seurat.Young_3@meta.data$Library <- "Young_3"
seurat.Young_4@meta.data$Library <- "Young_4"
seurat.Middle_2@meta.data$Library <- "Middle_2"
seurat.Middle_3@meta.data$Library <- "Middle_3"
seurat.Middle_4@meta.data$Library <- "Middle_4"
seurat.Old_1@meta.data$Library <- "Old_1"
seurat.Old_2@meta.data$Library <- "Old_2"
seurat.Old_3@meta.data$Library <- "Old_3"

seurat.Young_1@meta.data$Age <- "Young"
seurat.Young_3@meta.data$Age <- "Young"
seurat.Young_4@meta.data$Age <- "Young"
seurat.Middle_2@meta.data$Age <- "Middle"
seurat.Middle_3@meta.data$Age <- "Middle"
seurat.Middle_4@meta.data$Age <- "Middle"
seurat.Old_1@meta.data$Age <- "Old"
seurat.Old_2@meta.data$Age <- "Old"
seurat.Old_3@meta.data$Age <- "Old"

ovary.Human.GSE255690 <- merge(seurat.Young_1, 
                               y =  c(seurat.Young_3,
                                      seurat.Young_4,
                                      seurat.Middle_2,
                                      seurat.Middle_3,
                                      seurat.Middle_4,
                                      seurat.Old_1,
                                      seurat.Old_2,
                                      seurat.Old_3), 
                               project = "10x_ovary_Human_GSE255690")

save(ovary.Human.GSE255690, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_post-DecontX_Seurat_object.RData"))

################################################################################
# 2. Seurat QC
# Filter parameters: nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mt < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.Human.GSE255690 <- SetIdent(ovary.Human.GSE255690, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.Human.GSE255690@assays$RNA)
num.cells <- Matrix::rowSums(ovary.Human.GSE255690@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.Human.GSE255690.filt <- subset(ovary.Human.GSE255690, features = genes.use)
ovary.Human.GSE255690.filt
# An object of class Seurat 
# 22439 features across 100640 samples within 1 assay 
# Active assay: RNA (22439 features, 0 variable features)

############### QC - mitochondrial genes ###############

ovary.Human.GSE255690.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.Human.GSE255690.filt, pattern = "^MT-")
head(ovary.Human.GSE255690.filt@meta.data)
#                         orig.ident nCount_RNA nFeature_RNA decontX_contamination decontX_clusters Library Age percent.mito
# AAACCCAAGCAGCCCT-1_1 SeuratProject       4507         1733          4.632868e-05                1 Young_1 Young    10.472598
# AAACCCAAGTGGAAAG-1_1 SeuratProject      17470         3902          6.505248e-05                2 Young_1 Young     4.018317
# AAACCCACACGAGAAC-1_1 SeuratProject       5879         1759          2.898343e-05                3 Young_1 Young     9.508420
# AAACCCACACGGTGAA-1_1 SeuratProject      10184         2262          1.788863e-04                2 Young_1 Young     4.840927
# AAACCCAGTAGACGTG-1_1 SeuratProject       5687         1690          3.837302e-05                2 Young_1 Young     5.310357
# AAACCCAGTCAGGAGT-1_1 SeuratProject       6512         2046          5.793113e-05                2 Young_1 Young     9.966216

pdf(paste(Sys.Date(),"10x_ovary_Human_GSE255690_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Human.GSE255690.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.Human.GSE255690.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.Human.GSE255690.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.Human.GSE255690.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Human_GSE255690_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.Human.GSE255690.filt <- subset(ovary.Human.GSE255690.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)

ovary.Human.GSE255690.filt
#An object of class Seurat 
# 22439 features across 81979 samples within 1 assay 
# Active assay: RNA (22439 features, 0 variable features)

pdf(paste(Sys.Date(),"10x_ovary_Human_GSE255690_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Human.GSE255690.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.Human.GSE255690.filt, file = paste(Sys.Date(),"10x_ovary_Human_GSE255690_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.Human.GSE255690.filt <- NormalizeData(ovary.Human.GSE255690.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Human.GSE255690.filt <- SCTransform(object = ovary.Human.GSE255690.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Human.GSE255690.filt, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_Seurat_object_postSCT.RData"))

ovary.Human.GSE255690.filt <- RunPCA(ovary.Human.GSE255690.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_ElbowPlot.pdf"))
ElbowPlot(ovary.Human.GSE255690.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Human.GSE255690.filt[["pca"]]@stdev / sum(ovary.Human.GSE255690.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 11

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 11

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_Human_GSE255690_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Human.GSE255690.filt <- RunUMAP(ovary.Human.GSE255690.filt, dims = 1:pcs)
ovary.Human.GSE255690.filt <- FindNeighbors(ovary.Human.GSE255690.filt, dims = 1:pcs)
ovary.Human.GSE255690.filt <- FindClusters(object = ovary.Human.GSE255690.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.Human.GSE255690.filt, reduction = "umap")
dev.off()

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
# table(ovary.Human.GSE255690.filt$Library)

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
ovary.Human.GSE255690.filt.list <- SplitObject(ovary.Human.GSE255690.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.Human.GSE255690.filt.list, ncol))))/100
# pred.dblt.rate

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.Human.GSE255690.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.Human.GSE255690 <- paramSweep_v3(ovary.Human.GSE255690.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.Human.GSE255690    <- summarizeSweep(sweep.res.list.ovary.Human.GSE255690, GT = FALSE)
  bcmvn.ovary.Human.GSE255690          <- find.pK(sweep.stats.ovary.Human.GSE255690)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.Human.GSE255690 <- as.numeric(as.character(bcmvn.ovary.Human.GSE255690[as.numeric(bcmvn.ovary.Human.GSE255690$pK[bcmvn.ovary.Human.GSE255690$BCmetric == max(bcmvn.ovary.Human.GSE255690$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.Human.GSE255690.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.Human.GSE255690.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.Human.GSE255690.filt.list[[i]] <- doubletFinder_v3(ovary.Human.GSE255690.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.Human.GSE255690, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.Human.GSE255690.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Human.GSE255690.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.Human.GSE255690.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Human.GSE255690.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.Human.GSE255690.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_Human_GSE255690",names(ovary.Human.GSE255690.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.Human.GSE255690.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.Human.GSE255690.DFsinglets <- merge(ovary.Human.GSE255690.filt.list[[1]],
                                          y = c(ovary.Human.GSE255690.filt.list[[2]], ovary.Human.GSE255690.filt.list[[3]], 
                                                ovary.Human.GSE255690.filt.list[[4]], ovary.Human.GSE255690.filt.list[[5]], ovary.Human.GSE255690.filt.list[[6]], 
                                                ovary.Human.GSE255690.filt.list[[7]], ovary.Human.GSE255690.filt.list[[8]], ovary.Human.GSE255690.filt.list[[9]]),
                                          project = "ovary.Human.GSE255690")
# ovary.Human.GSE255690.DFsinglets

# remove pANN columns that are 10xGenomics library lane specific
ovary.Human.GSE255690.DFsinglets@meta.data <- ovary.Human.GSE255690.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.Human.GSE255690.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.Human.GSE255690.filt.list.scds <- lapply(ovary.Human.GSE255690.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.Human.GSE255690.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.Human.GSE255690.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.Human.GSE255690.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.Human.GSE255690.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.Human.GSE255690.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.Human.GSE255690.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.Human.GSE255690.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.Human.GSE255690.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.Human.GSE255690.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_Human_GSE255690", names(ovary.Human.GSE255690.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.Human.GSE255690.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.Human.GSE255690.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.Human.GSE255690.DFsinglets@meta.data[colnames(ovary.Human.GSE255690.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.Human.GSE255690.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.Human.GSE255690.DFsinglets@meta.data$DoubletFinder, ovary.Human.GSE255690.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet    1592    8092
# Singlet    5919   71333

ovary.Human.GSE255690.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.Human.GSE255690.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.Human.GSE255690.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                                 "Doublet", "Singlet")
table(ovary.Human.GSE255690.DFsinglets@meta.data$DoubletCall)
# Doublet Singlet 
#   14040   67939 

# re-run dimensionality reduction for plotting purposes
ovary.Human.GSE255690.DFsinglets <- SCTransform(object = ovary.Human.GSE255690.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))
ovary.Human.GSE255690.DFsinglets <- RunPCA(ovary.Human.GSE255690.DFsinglets, npcs = 50)
ovary.Human.GSE255690.DFsinglets <- RunUMAP(ovary.Human.GSE255690.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_DoubletCall_UMAP.pdf"))
DimPlot(ovary.Human.GSE255690.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.Human.GSE255690.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.Human.GSE255690.DFsinglets   <- subset(ovary.Human.GSE255690.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets
# ovary.Human.GSE255690.DFsinglets

# save filtered/annotated object
save(ovary.Human.GSE255690.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_Seurat_object_SINGLETS.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Human_GSE255690_processing_session_info.txt", sep =""))
sessionInfo()
sink()
