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
# 10x ovary Human GSE202601 data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.Human23 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human23/outs/filtered_feature_bc_matrix/")
counts.Human27 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human27/outs/filtered_feature_bc_matrix/")
counts.Human28 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human28/outs/filtered_feature_bc_matrix/")
counts.Human29 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human29/outs/filtered_feature_bc_matrix/")
counts.Human49 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human49/outs/filtered_feature_bc_matrix/")
counts.Human51 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human51/outs/filtered_feature_bc_matrix/")
counts.Human52 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human52/outs/filtered_feature_bc_matrix/")
counts.Human54 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human54/outs/filtered_feature_bc_matrix/")

counts.raw.Human23 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human23/outs/raw_feature_bc_matrix/")
counts.raw.Human27 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human27/outs/raw_feature_bc_matrix/")
counts.raw.Human28 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human28/outs/raw_feature_bc_matrix/")
counts.raw.Human29 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human29/outs/raw_feature_bc_matrix/")
counts.raw.Human49 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human49/outs/raw_feature_bc_matrix/")
counts.raw.Human51 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human51/outs/raw_feature_bc_matrix/")
counts.raw.Human52 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human52/outs/raw_feature_bc_matrix/")
counts.raw.Human54 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Human_GSE202601/2_Cellranger/Human54/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce.Human23 <- SingleCellExperiment(list(counts = counts.Human23))
sce.Human27 <- SingleCellExperiment(list(counts = counts.Human27))
sce.Human28 <- SingleCellExperiment(list(counts = counts.Human28))
sce.Human29 <- SingleCellExperiment(list(counts = counts.Human29))
sce.Human49 <- SingleCellExperiment(list(counts = counts.Human49))
sce.Human51 <- SingleCellExperiment(list(counts = counts.Human51))
sce.Human52 <- SingleCellExperiment(list(counts = counts.Human52))
sce.Human54 <- SingleCellExperiment(list(counts = counts.Human54))

sce.raw.Human23 <- SingleCellExperiment(list(counts = counts.raw.Human23))
sce.raw.Human27 <- SingleCellExperiment(list(counts = counts.raw.Human27))
sce.raw.Human28 <- SingleCellExperiment(list(counts = counts.raw.Human28))
sce.raw.Human29 <- SingleCellExperiment(list(counts = counts.raw.Human29))
sce.raw.Human49 <- SingleCellExperiment(list(counts = counts.raw.Human49))
sce.raw.Human51 <- SingleCellExperiment(list(counts = counts.raw.Human51))
sce.raw.Human52 <- SingleCellExperiment(list(counts = counts.raw.Human52))
sce.raw.Human54 <- SingleCellExperiment(list(counts = counts.raw.Human54))

sce.Human23 <- decontX(sce.Human23, background = sce.raw.Human23)
sce.Human27 <- decontX(sce.Human27, background = sce.raw.Human27)
sce.Human28 <- decontX(sce.Human28, background = sce.raw.Human28)
sce.Human29 <- decontX(sce.Human29, background = sce.raw.Human29)
sce.Human49 <- decontX(sce.Human49, background = sce.raw.Human49)
sce.Human51 <- decontX(sce.Human51, background = sce.raw.Human51)
sce.Human52 <- decontX(sce.Human52, background = sce.raw.Human52)
sce.Human54 <- decontX(sce.Human54, background = sce.raw.Human54)

# Save DecontX result
save(sce.Human23, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human23.RData"))
save(sce.Human27, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human27.RData"))
save(sce.Human28, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human28.RData"))
save(sce.Human29, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human29.RData"))
save(sce.Human49, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human49.RData"))
save(sce.Human51, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human51.RData"))
save(sce.Human52, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human52.RData"))
save(sce.Human54, file = paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_DecontX_SCE_object_sce.Human54.RData"))

# Plot UMAP
UMAP.Human23 <- reducedDim(sce.Human23, "decontX_UMAP")
UMAP.Human27 <- reducedDim(sce.Human27, "decontX_UMAP")
UMAP.Human28 <- reducedDim(sce.Human28, "decontX_UMAP")
UMAP.Human29 <- reducedDim(sce.Human29, "decontX_UMAP")
UMAP.Human49 <- reducedDim(sce.Human49, "decontX_UMAP")
UMAP.Human51 <- reducedDim(sce.Human51, "decontX_UMAP")
UMAP.Human52 <- reducedDim(sce.Human52, "decontX_UMAP")
UMAP.Human54 <- reducedDim(sce.Human54, "decontX_UMAP")

plotDecontXContamination(sce.Human23)
plotDecontXContamination(sce.Human27)
plotDecontXContamination(sce.Human28)
plotDecontXContamination(sce.Human29)
plotDecontXContamination(sce.Human49)
plotDecontXContamination(sce.Human51)
plotDecontXContamination(sce.Human52)
plotDecontXContamination(sce.Human54)

# Create Seurat objects
seurat.Human23 <- CreateSeuratObject(round(decontXcounts(sce.Human23)), meta.data=as.data.frame(colData(sce.Human23)))
seurat.Human27 <- CreateSeuratObject(round(decontXcounts(sce.Human27)), meta.data=as.data.frame(colData(sce.Human27)))
seurat.Human28 <- CreateSeuratObject(round(decontXcounts(sce.Human28)), meta.data=as.data.frame(colData(sce.Human28)))
seurat.Human29 <- CreateSeuratObject(round(decontXcounts(sce.Human29)), meta.data=as.data.frame(colData(sce.Human29)))
seurat.Human49 <- CreateSeuratObject(round(decontXcounts(sce.Human49)), meta.data=as.data.frame(colData(sce.Human49)))
seurat.Human51 <- CreateSeuratObject(round(decontXcounts(sce.Human51)), meta.data=as.data.frame(colData(sce.Human51)))
seurat.Human52 <- CreateSeuratObject(round(decontXcounts(sce.Human52)), meta.data=as.data.frame(colData(sce.Human52)))
seurat.Human54 <- CreateSeuratObject(round(decontXcounts(sce.Human54)), meta.data=as.data.frame(colData(sce.Human54)))

seurat.Human23@meta.data$Library <- "Human23"
seurat.Human27@meta.data$Library <- "Human27"
seurat.Human28@meta.data$Library <- "Human28"
seurat.Human29@meta.data$Library <- "Human29"
seurat.Human49@meta.data$Library <- "Human49"
seurat.Human51@meta.data$Library <- "Human51"
seurat.Human52@meta.data$Library <- "Human52"
seurat.Human54@meta.data$Library <- "Human54"

seurat.Human23@meta.data$Age <- 23
seurat.Human27@meta.data$Age <- 27
seurat.Human28@meta.data$Age <- 28
seurat.Human29@meta.data$Age <- 29
seurat.Human49@meta.data$Age <- 49
seurat.Human51@meta.data$Age <- 51
seurat.Human52@meta.data$Age <- 52
seurat.Human54@meta.data$Age <- 54

ovary.Human.GSE202601 <- merge(seurat.Human23, 
                               y =  c(seurat.Human27,
                                      seurat.Human28,
                                      seurat.Human29,
                                      seurat.Human49,
                                      seurat.Human51,
                                      seurat.Human52,
                                      seurat.Human54), 
                               project = "10x_ovary_Human_GSE202601")

################################################################################
# 2. Add metadata
################################################################################

# create Age label
Age <- rep("NA", length(ovary.Human.GSE202601@meta.data$Library))
Age[grep("Human23", ovary.Human.GSE202601@meta.data$Library)] <- 23
Age[grep("Human27", ovary.Human.GSE202601@meta.data$Library)] <- 27
Age[grep("Human28", ovary.Human.GSE202601@meta.data$Library)] <- 28
Age[grep("Human29", ovary.Human.GSE202601@meta.data$Library)] <- 29
Age[grep("Human49", ovary.Human.GSE202601@meta.data$Library)] <- 49
Age[grep("Human51", ovary.Human.GSE202601@meta.data$Library)] <- 51
Age[grep("Human52", ovary.Human.GSE202601@meta.data$Library)] <- 52
Age[grep("Human54", ovary.Human.GSE202601@meta.data$Library)] <- 54

Age <- data.frame(Age)
rownames(Age) <- colnames(ovary.Human.GSE202601@assays$RNA)

# create Group label
Group <- rep("NA", length(ovary.Human.GSE202601@meta.data$Library))
Group[grep("Human23", ovary.Human.GSE202601@meta.data$Library)] <- "Young"
Group[grep("Human27", ovary.Human.GSE202601@meta.data$Library)] <- "Young"
Group[grep("Human28", ovary.Human.GSE202601@meta.data$Library)] <- "Young"
Group[grep("Human29", ovary.Human.GSE202601@meta.data$Library)] <- "Young"
Group[grep("Human49", ovary.Human.GSE202601@meta.data$Library)] <- "Old"
Group[grep("Human51", ovary.Human.GSE202601@meta.data$Library)] <- "Old"
Group[grep("Human52", ovary.Human.GSE202601@meta.data$Library)] <- "Old"
Group[grep("Human54", ovary.Human.GSE202601@meta.data$Library)] <- "Old"

Group <- data.frame(Group)
rownames(Group) <- colnames(ovary.Human.GSE202601@assays$RNA)

# update Seurat with metadata
ovary.Human.GSE202601 <- AddMetaData(object = ovary.Human.GSE202601, metadata = as.vector(Age), col.name = "Age")
ovary.Human.GSE202601 <- AddMetaData(object = ovary.Human.GSE202601, metadata = as.vector(Group), col.name = "Group")

save(ovary.Human.GSE202601, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_post-DecontX_Seurat_object.RData"))

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mt < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.Human.GSE202601 <- SetIdent(ovary.Human.GSE202601, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.Human.GSE202601@assays$RNA)
num.cells <- Matrix::rowSums(ovary.Human.GSE202601@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.Human.GSE202601.filt <- subset(ovary.Human.GSE202601, features = genes.use)

############### QC - mitochondrial genes ###############

ovary.Human.GSE202601.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.Human.GSE202601.filt, pattern = "^MT-")

pdf(paste(Sys.Date(),"10x_ovary_Human_GSE202601_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Human.GSE202601.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.Human.GSE202601.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.Human.GSE202601.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.Human.GSE202601.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Human_GSE202601_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.Human.GSE202601.filt <- subset(ovary.Human.GSE202601.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)

pdf(paste(Sys.Date(),"10x_ovary_Human_GSE202601_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Human.GSE202601.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.Human.GSE202601.filt, file = paste(Sys.Date(),"10x_ovary_Human_GSE202601_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.Human.GSE202601.filt <- NormalizeData(ovary.Human.GSE202601.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Human.GSE202601.filt <- SCTransform(object = ovary.Human.GSE202601.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Human.GSE202601.filt, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_Seurat_object_postSCT.RData"))

ovary.Human.GSE202601.filt <- RunPCA(ovary.Human.GSE202601.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_ElbowPlot.pdf"))
ElbowPlot(ovary.Human.GSE202601.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Human.GSE202601.filt[["pca"]]@stdev / sum(ovary.Human.GSE202601.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

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
pdf(paste0(Sys.Date(), "_10x_ovary_Human_GSE202601_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Human.GSE202601.filt <- RunUMAP(ovary.Human.GSE202601.filt, dims = 1:pcs)
ovary.Human.GSE202601.filt <- FindNeighbors(ovary.Human.GSE202601.filt, dims = 1:pcs)
ovary.Human.GSE202601.filt <- FindClusters(object = ovary.Human.GSE202601.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.Human.GSE202601.filt, reduction = "umap")
dev.off()

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

########################################################
#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
ovary.Human.GSE202601.filt.list <- SplitObject(ovary.Human.GSE202601.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.Human.GSE202601.filt.list, ncol))))/100

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.Human.GSE202601.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.Human.GSE202601 <- paramSweep_v3(ovary.Human.GSE202601.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.Human.GSE202601    <- summarizeSweep(sweep.res.list.ovary.Human.GSE202601, GT = FALSE)
  bcmvn.ovary.Human.GSE202601          <- find.pK(sweep.stats.ovary.Human.GSE202601)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.Human.GSE202601 <- as.numeric(as.character(bcmvn.ovary.Human.GSE202601[as.numeric(bcmvn.ovary.Human.GSE202601$pK[bcmvn.ovary.Human.GSE202601$BCmetric == max(bcmvn.ovary.Human.GSE202601$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.Human.GSE202601.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.Human.GSE202601.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.Human.GSE202601.filt.list[[i]] <- doubletFinder_v3(ovary.Human.GSE202601.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.Human.GSE202601, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.Human.GSE202601.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Human.GSE202601.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.Human.GSE202601.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Human.GSE202601.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.Human.GSE202601.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_Human_GSE202601",names(ovary.Human.GSE202601.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.Human.GSE202601.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.Human.GSE202601.DFsinglets <- merge(ovary.Human.GSE202601.filt.list[[1]],
                                          y = c(ovary.Human.GSE202601.filt.list[[2]], ovary.Human.GSE202601.filt.list[[3]], ovary.Human.GSE202601.filt.list[[4]],
                                                ovary.Human.GSE202601.filt.list[[5]], ovary.Human.GSE202601.filt.list[[6]], ovary.Human.GSE202601.filt.list[[7]], ovary.Human.GSE202601.filt.list[[8]]),
                                          project = "ovary.Human.GSE202601")

# remove pANN columns that are 10xGenomics library lane specific
ovary.Human.GSE202601.DFsinglets@meta.data <- ovary.Human.GSE202601.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.Human.GSE202601.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.Human.GSE202601.filt.list.scds <- lapply(ovary.Human.GSE202601.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.Human.GSE202601.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.Human.GSE202601.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.Human.GSE202601.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.Human.GSE202601.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.Human.GSE202601.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.Human.GSE202601.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.Human.GSE202601.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.Human.GSE202601.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.Human.GSE202601.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_Human_GSE202601", names(ovary.Human.GSE202601.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.Human.GSE202601.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.Human.GSE202601.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.Human.GSE202601.DFsinglets@meta.data[colnames(ovary.Human.GSE202601.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.Human.GSE202601.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

ovary.Human.GSE202601.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.Human.GSE202601.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.Human.GSE202601.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                                 "Doublet", "Singlet")

# re-run dimensionality reduction for plotting purposes
ovary.Human.GSE202601.DFsinglets <- SCTransform(object = ovary.Human.GSE202601.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))
ovary.Human.GSE202601.DFsinglets <- RunPCA(ovary.Human.GSE202601.DFsinglets, npcs = 50)
ovary.Human.GSE202601.DFsinglets <- RunUMAP(ovary.Human.GSE202601.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_DoubletCall_UMAP.pdf"))
DimPlot(ovary.Human.GSE202601.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.Human.GSE202601.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.Human.GSE202601.DFsinglets   <- subset(ovary.Human.GSE202601.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets

# save filtered/annotated object
save(ovary.Human.GSE202601.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_Seurat_object_SINGLETS.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Human_GSE202601_preprocessing_session_info.txt", sep =""))
sessionInfo()
sink()