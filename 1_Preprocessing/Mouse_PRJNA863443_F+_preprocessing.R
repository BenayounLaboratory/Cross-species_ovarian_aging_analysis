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
# 10x ovary PRJNA863443 F+ data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output
counts.Foxl2_wt_old_0 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt/outs/filtered_feature_bc_matrix/")
counts.Foxl2_wt_old_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_old_1/outs/filtered_feature_bc_matrix/")
counts.Foxl2_wt_old_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_old_2/outs/filtered_feature_bc_matrix/")
counts.Foxl2_wt_young_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_young_1/outs/filtered_feature_bc_matrix/")
counts.Foxl2_wt_young_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_young_2/outs/filtered_feature_bc_matrix/")

counts.raw.Foxl2_wt_old_0 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt/outs/raw_feature_bc_matrix/")
counts.raw.Foxl2_wt_old_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_old_1/outs/raw_feature_bc_matrix/")
counts.raw.Foxl2_wt_old_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_old_2/outs/raw_feature_bc_matrix/")
counts.raw.Foxl2_wt_young_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_young_1/outs/raw_feature_bc_matrix/")
counts.raw.Foxl2_wt_young_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/Foxl2/Foxl2_wt_young_2/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce.Foxl2_wt_old_0 <- SingleCellExperiment(list(counts = counts.Foxl2_wt_old_0))
sce.Foxl2_wt_old_1 <- SingleCellExperiment(list(counts = counts.Foxl2_wt_old_1))
sce.Foxl2_wt_old_2 <- SingleCellExperiment(list(counts = counts.Foxl2_wt_old_2))
sce.Foxl2_wt_young_1 <- SingleCellExperiment(list(counts = counts.Foxl2_wt_young_1))
sce.Foxl2_wt_young_2 <- SingleCellExperiment(list(counts = counts.Foxl2_wt_young_2))

sce.raw.Foxl2_wt_old_0 <- SingleCellExperiment(list(counts = counts.raw.Foxl2_wt_old_0))
sce.raw.Foxl2_wt_old_1 <- SingleCellExperiment(list(counts = counts.raw.Foxl2_wt_old_1))
sce.raw.Foxl2_wt_old_2 <- SingleCellExperiment(list(counts = counts.raw.Foxl2_wt_old_2))
sce.raw.Foxl2_wt_young_1 <- SingleCellExperiment(list(counts = counts.raw.Foxl2_wt_young_1))
sce.raw.Foxl2_wt_young_2 <- SingleCellExperiment(list(counts = counts.raw.Foxl2_wt_young_2))

sce.Foxl2_wt_old_0 <- decontX(sce.Foxl2_wt_old_0, background = sce.raw.Foxl2_wt_old_0)
sce.Foxl2_wt_old_1 <- decontX(sce.Foxl2_wt_old_1, background = sce.raw.Foxl2_wt_old_1)
sce.Foxl2_wt_old_2 <- decontX(sce.Foxl2_wt_old_2, background = sce.raw.Foxl2_wt_old_2)
sce.Foxl2_wt_young_1 <- decontX(sce.Foxl2_wt_young_1, background = sce.raw.Foxl2_wt_young_1)
sce.Foxl2_wt_young_2 <- decontX(sce.Foxl2_wt_young_2, background = sce.raw.Foxl2_wt_young_2)

# Save DecontX result
save(sce.Foxl2_wt_old_0, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_DecontX_SCE_object_sce.Foxl2_wt_old_0.RData"))
save(sce.Foxl2_wt_old_1, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_DecontX_SCE_object_sce.Foxl2_wt_old_1.RData"))
save(sce.Foxl2_wt_old_2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_DecontX_SCE_object_sce.Foxl2_wt_old_2.RData"))
save(sce.Foxl2_wt_young_1, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_DecontX_SCE_object_sce.Foxl2_wt_young_1.RData"))
save(sce.Foxl2_wt_young_2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_DecontX_SCE_object_sce.Foxl2_wt_young_2.RData"))

# Plot UMAP
UMAP.Foxl2_wt_old_0 <- reducedDim(sce.Foxl2_wt_old_0, "decontX_UMAP")
UMAP.Foxl2_wt_old_1 <- reducedDim(sce.Foxl2_wt_old_1, "decontX_UMAP")
UMAP.Foxl2_wt_old_2 <- reducedDim(sce.Foxl2_wt_old_2, "decontX_UMAP")
UMAP.Foxl2_wt_young_1 <- reducedDim(sce.Foxl2_wt_young_1, "decontX_UMAP")
UMAP.Foxl2_wt_young_2 <- reducedDim(sce.Foxl2_wt_young_2, "decontX_UMAP")

plotDecontXContamination(sce.Foxl2_wt_old_0)
plotDecontXContamination(sce.Foxl2_wt_old_1)
plotDecontXContamination(sce.Foxl2_wt_old_2)
plotDecontXContamination(sce.Foxl2_wt_young_1)
plotDecontXContamination(sce.Foxl2_wt_young_2)

# Create Seurat objects
seurat.Foxl2_wt_old_0 <- CreateSeuratObject(round(decontXcounts(sce.Foxl2_wt_old_0)), meta.data=as.data.frame(colData(sce.Foxl2_wt_old_0)))
seurat.Foxl2_wt_old_1 <- CreateSeuratObject(round(decontXcounts(sce.Foxl2_wt_old_1)), meta.data=as.data.frame(colData(sce.Foxl2_wt_old_1)))
seurat.Foxl2_wt_old_2 <- CreateSeuratObject(round(decontXcounts(sce.Foxl2_wt_old_2)), meta.data=as.data.frame(colData(sce.Foxl2_wt_old_2)))
seurat.Foxl2_wt_young_1 <- CreateSeuratObject(round(decontXcounts(sce.Foxl2_wt_young_1)), meta.data=as.data.frame(colData(sce.Foxl2_wt_young_1)))
seurat.Foxl2_wt_young_2 <- CreateSeuratObject(round(decontXcounts(sce.Foxl2_wt_young_2)), meta.data=as.data.frame(colData(sce.Foxl2_wt_young_2)))

ovary.Foxl2 <- merge(seurat.Foxl2_wt_young_1, 
                     y =  c(seurat.Foxl2_wt_young_2,
                            seurat.Foxl2_wt_old_0,
                            seurat.Foxl2_wt_old_1,
                            seurat.Foxl2_wt_old_2), 
                     add.cell.ids = c("Foxl2_wt_young_1", "Foxl2_wt_young_2",
                                      "Foxl2_wt_old_0", "Foxl2_wt_old_1", "Foxl2_wt_old_2"),
                     project = "10x_ovary_Foxl2")

save(ovary.Foxl2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_post-DecontX_Seurat_object.RData"))

# Remove intermediate files
rm(counts.Foxl2_wt_old_0, counts.Foxl2_wt_old_1, counts.Foxl2_wt_old_2, 
   counts.Foxl2_wt_young_1, counts.Foxl2_wt_young_2, 
   counts.raw.Foxl2_wt_old_0, counts.raw.Foxl2_wt_old_1, counts.raw.Foxl2_wt_old_2, 
   counts.raw.Foxl2_wt_young_1, counts.raw.Foxl2_wt_young_2, 
   sce.Foxl2_wt_old_0, sce.Foxl2_wt_old_1, sce.Foxl2_wt_old_2, 
   sce.Foxl2_wt_young_1, sce.Foxl2_wt_young_2, 
   sce.raw.Foxl2_wt_old_0, sce.raw.Foxl2_wt_old_1, sce.raw.Foxl2_wt_old_2, 
   sce.raw.Foxl2_wt_young_1, sce.raw.Foxl2_wt_young_2, 
   seurat.Foxl2_wt_old_0, seurat.Foxl2_wt_old_1, seurat.Foxl2_wt_old_2, 
   seurat.Foxl2_wt_young_1, seurat.Foxl2_wt_young_2)

################################################################################
# 2. Add metadata to Seurat object
################################################################################

# create Library label
Library <- rep("NA", length(colnames(ovary.Foxl2@assays$RNA)))
Library[grep("Foxl2_wt_young_1", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_wt_young_1"
Library[grep("Foxl2_wt_young_2", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_wt_young_2"
Library[grep("Foxl2_wt_old_0", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_wt_old_0"
Library[grep("Foxl2_wt_old_1", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_wt_old_1"
Library[grep("Foxl2_wt_old_2", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_wt_old_2"

Library <- data.frame(Library)
rownames(Library) <- colnames(ovary.Foxl2@assays$RNA)

# create Age label
Age <- rep("NA", length(colnames(ovary.Foxl2@assays$RNA)))
Age[grep("young", colnames(ovary.Foxl2@assays$RNA))] <- "Young"
Age[grep("old", colnames(ovary.Foxl2@assays$RNA))] <- "Old"

Age <- data.frame(Age)
rownames(Age) <- colnames(ovary.Foxl2@assays$RNA)

# create Batch label
Batch <- rep("Foxl2_1", length(colnames(ovary.Foxl2@assays$RNA)))
Batch[grep("_1", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_2"
Batch[grep("_2", colnames(ovary.Foxl2@assays$RNA))] <- "Foxl2_3"

Batch <- data.frame(Batch)
rownames(Batch) <- colnames(ovary.Foxl2@assays$RNA)

# update Seurat with metadata
ovary.Foxl2 <- AddMetaData(object = ovary.Foxl2, metadata = as.vector(Library), col.name = "Library")
ovary.Foxl2 <- AddMetaData(object = ovary.Foxl2, metadata = as.vector(Age), col.name = "Age")
ovary.Foxl2 <- AddMetaData(object = ovary.Foxl2, metadata = as.vector(Batch), col.name = "Batch")

table(ovary.Foxl2@meta.data$Library)

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500
#                    500 < nCount_RNA < 1e5
#                    percent.mito < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.Foxl2 <- SetIdent(ovary.Foxl2, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.Foxl2@assays$RNA)
num.cells <- Matrix::rowSums(ovary.Foxl2@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.Foxl2.filt <- subset(ovary.Foxl2, features = genes.use)

############### QC - mitochondrial genes ###############

ovary.Foxl2.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.Foxl2.filt, pattern = "^mt-")
head(ovary.Foxl2.filt@meta.data)
#                             orig.ident nCount_RNA nFeature_RNA decontX_contamination decontX_clusters          Library   Age Genotype   Batch percent.mito
# Foxl2_wt_young_1_AAACCCACAACCGCTG-1 SeuratProject      33243         6397           0.048689830                1 Foxl2_wt_young_1 Young       wt Foxl2_2     1.837981
# Foxl2_wt_young_1_AAACCCACAAGTGATA-1 SeuratProject      19951         4841           0.018818160                2 Foxl2_wt_young_1 Young       wt Foxl2_2     8.876748
# Foxl2_wt_young_1_AAACCCATCAAAGGAT-1 SeuratProject       6522         2975           0.177393723                1 Foxl2_wt_young_1 Young       wt Foxl2_2     1.594603
# Foxl2_wt_young_1_AAACCCATCTCGCGTT-1 SeuratProject      15434         4308           0.084202047                3 Foxl2_wt_young_1 Young       wt Foxl2_2     5.105611
# Foxl2_wt_young_1_AAACGAAAGCAACCAG-1 SeuratProject      21902         5775           0.024149265                4 Foxl2_wt_young_1 Young       wt Foxl2_2     7.638572
# Foxl2_wt_young_1_AAACGAACACCAATTG-1 SeuratProject      44858         6844           0.007175263                2 Foxl2_wt_young_1 Young       wt Foxl2_2     9.240269

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_Foxl2_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Foxl2.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.Foxl2.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.Foxl2.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.Foxl2.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_Foxl2_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.Foxl2.filt <- subset(ovary.Foxl2.filt, subset = nFeature_RNA > 600 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)
# ovary.Foxl2.filt

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_Foxl2_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Foxl2.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.Foxl2.filt, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_Foxl2_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.Foxl2.filt <- NormalizeData(ovary.Foxl2.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Foxl2.filt <- SCTransform(object = ovary.Foxl2.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library"))

save(ovary.Foxl2.filt, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_postSCT.RData"))

ovary.Foxl2.filt <- RunPCA(ovary.Foxl2.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_ElbowPlot.pdf"))
ElbowPlot(ovary.Foxl2.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Foxl2.filt[["pca"]]@stdev / sum(ovary.Foxl2.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 20

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 20

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Foxl2.filt <- RunUMAP(ovary.Foxl2.filt, dims = 1:pcs)
ovary.Foxl2.filt <- FindNeighbors(ovary.Foxl2.filt, dims = 1:pcs)
ovary.Foxl2.filt <- FindClusters(object = ovary.Foxl2.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.Foxl2.filt, reduction = "umap")
dev.off()

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
table(ovary.Foxl2.filt$Library)
#   Foxl2_wt_old_0   Foxl2_wt_old_1   Foxl2_wt_old_2 Foxl2_wt_young_1 Foxl2_wt_young_2 
#             2180             3156             3916             4093             3753 

# Clean memory of intermediates
rm(ovary.Foxl2)

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
ovary.Foxl2.filt.list <- SplitObject(ovary.Foxl2.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.Foxl2.filt.list, ncol))))/100
# pred.dblt.rate

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.Foxl2.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.Foxl2 <- paramSweep_v3(ovary.Foxl2.filt.list[[i]], PCs = 1:pcs, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.Foxl2    <- summarizeSweep(sweep.res.list.ovary.Foxl2, GT = FALSE)
  bcmvn.ovary.Foxl2          <- find.pK(sweep.stats.ovary.Foxl2)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.Foxl2 <- as.numeric(as.character(bcmvn.ovary.Foxl2[as.numeric(bcmvn.ovary.Foxl2$pK[bcmvn.ovary.Foxl2$BCmetric == max(bcmvn.ovary.Foxl2$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.Foxl2.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.Foxl2.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.Foxl2.filt.list[[i]] <- doubletFinder_v3(ovary.Foxl2.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.Foxl2, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.Foxl2.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Foxl2.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.Foxl2.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Foxl2.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.Foxl2.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_Foxl2",names(ovary.Foxl2.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.Foxl2.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.Foxl2.DFsinglets <- merge(ovary.Foxl2.filt.list[[1]],
                                y = c(ovary.Foxl2.filt.list[[2]], ovary.Foxl2.filt.list[[3]], ovary.Foxl2.filt.list[[4]],
                                      ovary.Foxl2.filt.list[[5]]),
                                project = "ovary.Foxl2")
ovary.Foxl2.DFsinglets
# An object of class Seurat 
# 44756 features across 26469 samples within 2 assays 
# Active assay: SCT (22377 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
ovary.Foxl2.DFsinglets@meta.data <- ovary.Foxl2.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.Foxl2.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.Foxl2.filt.list.scds <- lapply(ovary.Foxl2.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.Foxl2.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.Foxl2.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.Foxl2.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.Foxl2.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.Foxl2.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.Foxl2.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.Foxl2.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.Foxl2.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.Foxl2.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_Foxl2", names(ovary.Foxl2.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.Foxl2.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.Foxl2.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.Foxl2.DFsinglets@meta.data[colnames(ovary.Foxl2.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.Foxl2.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.Foxl2.DFsinglets@meta.data$DoubletFinder, ovary.Foxl2.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet      66     898
# Singlet     470   15664

ovary.Foxl2.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.Foxl2.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.Foxl2.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                       "Doublet", "Singlet")
table(ovary.Foxl2.DFsinglets@meta.data$DoubletCall)
# Doublet Singlet 
#    1434   15664 

# re-run dimensionality reduction for plotting purposes
ovary.Foxl2.DFsinglets <- SCTransform(object = ovary.Foxl2.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library", "decontX_contamination", "Batch"))
ovary.Foxl2.DFsinglets <- RunPCA(ovary.Foxl2.DFsinglets, npcs = 50)
ovary.Foxl2.DFsinglets <- RunUMAP(ovary.Foxl2.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_DoubletCall_UMAP.pdf"))
DimPlot(ovary.Foxl2.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.Foxl2.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.Foxl2.DFsinglets   <- subset(ovary.Foxl2.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets

# table(ovary.Foxl2.DFsinglets@meta.data$Library)

# save filtered/annotated object
save(ovary.Foxl2.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_SINGLETS.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Mouse_PRJNA863443_F+_preprocessing_session_info.txt", sep =""))
sessionInfo()
sink()
