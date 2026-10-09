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
# 10x ovary PRJNA863443 AC data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.YF_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/YF_1/outs/filtered_feature_bc_matrix/")
counts.YF_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/YF_2/outs/filtered_feature_bc_matrix/")
counts.YF_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/YF_3/outs/filtered_feature_bc_matrix/")
counts.OF_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/OF_1/outs/filtered_feature_bc_matrix/")
counts.OF_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/OF_2/outs/filtered_feature_bc_matrix/")
counts.OF_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/OF_3/outs/filtered_feature_bc_matrix/")

counts.raw.YF_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/YF_1/outs/raw_feature_bc_matrix/")
counts.raw.YF_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/YF_2/outs/raw_feature_bc_matrix/")
counts.raw.YF_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/YF_3/outs/raw_feature_bc_matrix/")
counts.raw.OF_1 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/OF_1/outs/raw_feature_bc_matrix/")
counts.raw.OF_2 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/OF_2/outs/raw_feature_bc_matrix/")
counts.raw.OF_3 <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Benayoun_lab/AC/OF_3/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX

sce.YF_1 <- SingleCellExperiment(list(counts = counts.YF_1))
sce.YF_2 <- SingleCellExperiment(list(counts = counts.YF_2))
sce.YF_3 <- SingleCellExperiment(list(counts = counts.YF_3))
sce.OF_1 <- SingleCellExperiment(list(counts = counts.OF_1))
sce.OF_2 <- SingleCellExperiment(list(counts = counts.OF_2))
sce.OF_3 <- SingleCellExperiment(list(counts = counts.OF_3))

sce.raw.YF_1 <- SingleCellExperiment(list(counts = counts.raw.YF_1))
sce.raw.YF_2 <- SingleCellExperiment(list(counts = counts.raw.YF_2))
sce.raw.YF_3 <- SingleCellExperiment(list(counts = counts.raw.YF_3))
sce.raw.OF_1 <- SingleCellExperiment(list(counts = counts.raw.OF_1))
sce.raw.OF_2 <- SingleCellExperiment(list(counts = counts.raw.OF_2))
sce.raw.OF_3 <- SingleCellExperiment(list(counts = counts.raw.OF_3))

sce.YF_1 <- decontX(sce.YF_1, background = sce.raw.YF_1)
sce.YF_2 <- decontX(sce.YF_2, background = sce.raw.YF_2)
sce.YF_3 <- decontX(sce.YF_3, background = sce.raw.YF_3)
sce.OF_1 <- decontX(sce.OF_1, background = sce.raw.OF_1)
sce.OF_2 <- decontX(sce.OF_2, background = sce.raw.OF_2)
sce.OF_3 <- decontX(sce.OF_3, background = sce.raw.OF_3)

# Plot UMAP
UMAP.YF_1 <- reducedDim(sce.YF_1, "decontX_UMAP")
UMAP.YF_2 <- reducedDim(sce.YF_2, "decontX_UMAP")
UMAP.YF_3 <- reducedDim(sce.YF_3, "decontX_UMAP")
UMAP.OF_1 <- reducedDim(sce.OF_1, "decontX_UMAP")
UMAP.OF_2 <- reducedDim(sce.OF_2, "decontX_UMAP")
UMAP.OF_3 <- reducedDim(sce.OF_3, "decontX_UMAP")

plotDecontXContamination(sce.YF_1)
plotDecontXContamination(sce.YF_2)
plotDecontXContamination(sce.YF_3)
plotDecontXContamination(sce.OF_1)
plotDecontXContamination(sce.OF_2)
plotDecontXContamination(sce.OF_3)

# Plot example box plots of contamination scores
boxplot(sce.YF_1$decontX_contamination)
boxplot(sce.YF_2$decontX_contamination)
boxplot(sce.YF_3$decontX_contamination)
boxplot(sce.OF_1$decontX_contamination)
boxplot(sce.OF_2$decontX_contamination)
boxplot(sce.OF_3$decontX_contamination)

# Create Seurat objects
seurat.YF_1 <- CreateSeuratObject(round(decontXcounts(sce.YF_1)), meta.data=as.data.frame(colData(sce.YF_1)))
seurat.YF_2 <- CreateSeuratObject(round(decontXcounts(sce.YF_2)), meta.data=as.data.frame(colData(sce.YF_2)))
seurat.YF_3 <- CreateSeuratObject(round(decontXcounts(sce.YF_3)), meta.data=as.data.frame(colData(sce.YF_3)))
seurat.OF_1 <- CreateSeuratObject(round(decontXcounts(sce.OF_1)), meta.data=as.data.frame(colData(sce.OF_1)))
seurat.OF_2 <- CreateSeuratObject(round(decontXcounts(sce.OF_2)), meta.data=as.data.frame(colData(sce.OF_2)))
seurat.OF_3 <- CreateSeuratObject(round(decontXcounts(sce.OF_3)), meta.data=as.data.frame(colData(sce.OF_3)))

ovary.AC <- merge(seurat.YF_1, 
                  y =  c(seurat.YF_2,
                         seurat.YF_3,
                         seurat.OF_1,
                         seurat.OF_2,
                         seurat.OF_3), 
                  add.cell.ids = c("YF_1",
                                   "YF_2",
                                   "YF_3",
                                   "OF_1",
                                   "OF_2",
                                   "OF_3"), 
                  project = "10x_ovary_AC")

################################################################################
# 2. Add metadata to Seurat object
################################################################################

# create Library label
Library <- rep("NA", length(colnames(ovary.AC@assays$RNA)))
Library[grep("YF_1", colnames(ovary.AC@assays$RNA))] <- "YF1"
Library[grep("YF_2", colnames(ovary.AC@assays$RNA))] <- "YF2"
Library[grep("YF_3", colnames(ovary.AC@assays$RNA))] <- "YF3"
Library[grep("OF_1", colnames(ovary.AC@assays$RNA))] <- "OF1"
Library[grep("OF_2", colnames(ovary.AC@assays$RNA))] <- "OF2"
Library[grep("OF_3", colnames(ovary.AC@assays$RNA))] <- "OF3"

Library <- data.frame(Library)
rownames(Library) <- colnames(ovary.AC@assays$RNA)

# create Age label
Age <- rep("NA", length(colnames(ovary.AC@assays$RNA)))
Age[grep("YF", colnames(ovary.AC@assays$RNA))] <- "YF"
Age[grep("OF", colnames(ovary.AC@assays$RNA))] <- "OF"

Age <- data.frame(Age)
rownames(Age) <- colnames(ovary.AC@assays$RNA)

# create Batch label
Batch <- rep("AC_1", length(colnames(ovary.AC@assays$RNA)))

Batch <- data.frame(Batch)
rownames(Batch) <- colnames(ovary.AC@assays$RNA)

# update Seurat with metadata
ovary.AC <- AddMetaData(object = ovary.AC, metadata = as.vector(Library), col.name = "Library")
ovary.AC <- AddMetaData(object = ovary.AC, metadata = as.vector(Age), col.name = "Age")
ovary.AC <- AddMetaData(object = ovary.AC, metadata = as.vector(Batch), col.name = "Batch")

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500
#                    500 < nCount_RNA < 1e5
#                    percent.mito < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.AC <- SetIdent(ovary.AC, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.AC@assays$RNA)
num.cells <- Matrix::rowSums(ovary.AC@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.AC.filt <- subset(ovary.AC, features = genes.use)

############### QC - mitochondrial genes ###############

ovary.AC.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.AC.filt, pattern = "^mt-")

# get clean cells
ovary.AC.filt <- subset(ovary.AC.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)

# Save data
save(ovary.AC.filt, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_AC_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.AC.filt <- NormalizeData(ovary.AC.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.AC.filt <- SCTransform(object = ovary.AC.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library"))

save(ovary.AC.filt, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_Seurat_object_postSCT.RData"))

ovary.AC.filt <- RunPCA(ovary.AC.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_AC_ElbowPlot.pdf"))
ElbowPlot(ovary.AC.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.AC.filt[["pca"]]@stdev / sum(ovary.AC.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# Minimum of the two calculation
pcs <- min(co1, co2)

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

###############################################################################
# run dimensionality reduction algorithm
ovary.AC.filt <- RunUMAP(ovary.AC.filt, dims = 1:pcs)
ovary.AC.filt <- FindNeighbors(ovary.AC.filt, dims = 1:pcs)
ovary.AC.filt <- FindClusters(object = ovary.AC.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.AC.filt, reduction = "umap")
dev.off()

################################################################################
# 4. Doublet detection part 1: DoubletFinder
################################################################################

pred.10x.dblt <- data.frame( "cell_number" = c(3000,4000,5000,6000,7000,8000,9000, 10000, 12000, 14000, 16000, 18000, 20000),
                             "dblt_rate"   = c(2.4 ,3.2 ,4.0 ,4.8 ,5.6 ,6.4 ,7.2 , 8.0  , 9.6, 11.2, 12.8, 14.4, 16.0))

pred_dblt_lm <- lm(dblt_rate ~ cell_number, data = pred.10x.dblt)

# Use distinct recovered cell count number to estimate doublets in dataset
#predict(pred_dblt_lm, data.frame("cell_number" = 6000))

########################################################
#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
ovary.AC.filt.list <- SplitObject(ovary.AC.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.AC.filt.list, ncol))))/100

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.AC.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.AC <- paramSweep_v3(ovary.AC.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.AC    <- summarizeSweep(sweep.res.list.ovary.AC, GT = FALSE)
  bcmvn.ovary.AC          <- find.pK(sweep.stats.ovary.AC)
  
  pk.ovary.AC <- as.numeric(as.character(bcmvn.ovary.AC[as.numeric(bcmvn.ovary.AC$pK[bcmvn.ovary.AC$BCmetric == max(bcmvn.ovary.AC$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.AC.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.AC.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.AC.filt.list[[i]] <- doubletFinder_v3(ovary.AC.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.AC, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.AC.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.AC.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.AC.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.AC.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.AC.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_AC",names(ovary.AC.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.AC.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.AC.DFsinglets <- merge(ovary.AC.filt.list[[1]],
                             y = ovary.AC.filt.list[-1],
                             project = "ovary.AC")

# remove pANN columns that are 10xGenomics library lane specific
ovary.AC.DFsinglets@meta.data <- ovary.AC.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.AC.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.AC.filt.list.scds <- lapply(ovary.AC.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.AC.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.AC.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.AC.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.AC.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.AC.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.AC.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.AC.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.AC.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.AC.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_AC", names(ovary.AC.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.AC.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.AC.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.AC.DFsinglets@meta.data[colnames(ovary.AC.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.AC.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

ovary.AC.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.AC.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.AC.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                    "Doublet", "Singlet")

# re-run dimensionality reduction for plotting purposes
ovary.AC.DFsinglets <- SCTransform(object = ovary.AC.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library", "decontX_contamination"))
ovary.AC.DFsinglets <- RunPCA(ovary.AC.DFsinglets, npcs = 50)
ovary.AC.DFsinglets <- RunUMAP(ovary.AC.DFsinglets, dims = 1:50)

# save annotated object
save(ovary.AC.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.AC.DFsinglets   <- subset(ovary.AC.DFsinglets, subset = DoubletCall %in% "Singlet")  

# save filtered/annotated object
save(ovary.AC.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_Seurat_object_SINGLETS.RData"))

################################################################################
# 7. Assess seurat clusters that are found only from one library
################################################################################

t.cluster.vs.library <- as.matrix(table(ovary.AC.DFsinglets@meta.data$seurat_clusters, ovary.AC.DFsinglets@meta.data$Library))

# Filter clusters that are represented by a single library
# First, create a binary matrix
binary_matrix <- t.cluster.vs.library > 0

# Sum across rows to find the number of libraries contributing to each cluster
libraries_per_cluster <- rowSums(binary_matrix)

# Get clusters that have cells from only one library
single_library_clusters <- names(libraries_per_cluster)[libraries_per_cluster == 1]

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Mouse_PRJNA863443_AC_preprocessing_session_info.txt", sep =""))
sessionInfo()
sink()