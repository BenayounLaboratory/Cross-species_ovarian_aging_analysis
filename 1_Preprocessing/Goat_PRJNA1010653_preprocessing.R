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
# 10x ovary Goat PRJNA1010653 data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output

counts.young <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Goat_PRJNA1010653/young_goat/outs/filtered_feature_bc_matrix/")
counts.aging <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Goat_PRJNA1010653/aging_goat/outs/filtered_feature_bc_matrix/")

counts.raw.young <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Goat_PRJNA1010653/young_goat/outs/raw_feature_bc_matrix/")
counts.raw.aging <- Read10X("/Volumes/OIProject_I/CZI_data/1_Cellranger/Goat_PRJNA1010653/aging_goat/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX

sce.young <- SingleCellExperiment(list(counts = counts.young))
sce.aging <- SingleCellExperiment(list(counts = counts.aging))

sce.raw.young <- SingleCellExperiment(list(counts = counts.raw.young))
sce.raw.aging <- SingleCellExperiment(list(counts = counts.raw.aging))

sce.young <- decontX(sce.young, background = sce.raw.young)
sce.aging <- decontX(sce.aging, background = sce.raw.aging)

# Save DecontX result
save(sce.young, file = paste0(Sys.Date(), "_10x_ovary_Goat_PRJNA1010653_DecontX_SCE_object_sce.young.RData"))
save(sce.aging, file = paste0(Sys.Date(), "_10x_ovary_Goat_PRJNA1010653_DecontX_SCE_object_sce.aging.RData"))

# Plot UMAP
UMAP.young <- reducedDim(sce.young, "decontX_UMAP")
UMAP.aging <- reducedDim(sce.aging, "decontX_UMAP")

plotDecontXContamination(sce.young)
plotDecontXContamination(sce.aging)

# Create Seurat objects
seurat.young <- CreateSeuratObject(round(decontXcounts(sce.young)), meta.data=as.data.frame(colData(sce.young)))
seurat.aging <- CreateSeuratObject(round(decontXcounts(sce.aging)), meta.data=as.data.frame(colData(sce.aging)))

seurat.young@meta.data$Library <- "young"
seurat.aging@meta.data$Library <- "aging"

ovary.Goat.PRJNA1010653 <- merge(seurat.young, 
                                 y =  c(seurat.aging), 
                                 project = "10x_ovary_Goat_PRJNA1010653")

ovary.Goat.PRJNA1010653
# An object of class Seurat 
# 21343 features across 16916 samples within 1 assay 
# Active assay: RNA (21343 features, 0 variable features)

################################################################################
# 2. Add metadata
################################################################################

# create Group label
Group <- rep("NA", length(ovary.Goat.PRJNA1010653@meta.data$Library))
Group[grep("young", ovary.Goat.PRJNA1010653@meta.data$Library)] <- "young"
Group[grep("aging", ovary.Goat.PRJNA1010653@meta.data$Library)] <- "aging"

Group <- data.frame(Group)
rownames(Group) <- colnames(ovary.Goat.PRJNA1010653@assays$RNA)

# update Seurat with metadata
ovary.Goat.PRJNA1010653 <- AddMetaData(object = ovary.Goat.PRJNA1010653, metadata = as.vector(Group), col.name = "Group")

save(ovary.Goat.PRJNA1010653, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_post-DecontX_Seurat_object.RData"))

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mt < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.Goat.PRJNA1010653 <- SetIdent(ovary.Goat.PRJNA1010653, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.Goat.PRJNA1010653@assays$RNA)
num.cells <- Matrix::rowSums(ovary.Goat.PRJNA1010653@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.Goat.PRJNA1010653.filt <- subset(ovary.Goat.PRJNA1010653, features = genes.use)

############### QC - mitochondrial genes ###############

ovary.Goat.PRJNA1010653.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.Goat.PRJNA1010653.filt, pattern = "^MT")
head(ovary.Goat.PRJNA1010653.filt@meta.data)
#                         orig.ident nCount_RNA nFeature_RNA decontX_contamination decontX_clusters Library Group percent.mito
# AAACCCACAGTCAGAG-1_1 SeuratProject       6331         2328             0.5714459                1   young young    0.1105671
# AAACCCACATCATGAC-1_1 SeuratProject      18686         4099             0.1665879                1   young young    0.2033608
# AAACCCACATGGGTCC-1_1 SeuratProject        918          711             0.5158231                1   young young    0.3267974
# AAACCCAGTGTCTTAG-1_1 SeuratProject      24195         4478             0.2745623                1   young young    0.2190535
# AAACCCAGTGTTACAC-1_1 SeuratProject       1534          536             0.9210672                1   young young    0.1303781
# AAACCCAGTTGTCCCT-1_1 SeuratProject        499          391             0.1026314                1   young young    0.0000000

pdf(paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Goat.PRJNA1010653.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.Goat.PRJNA1010653.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.Goat.PRJNA1010653.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.Goat.PRJNA1010653.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.Goat.PRJNA1010653.filt <- subset(ovary.Goat.PRJNA1010653.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)
ovary.Goat.PRJNA1010653.filt

pdf(paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.Goat.PRJNA1010653.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.Goat.PRJNA1010653.filt, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.Goat.PRJNA1010653.filt <- NormalizeData(ovary.Goat.PRJNA1010653.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Goat.PRJNA1010653.filt <- SCTransform(object = ovary.Goat.PRJNA1010653.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Goat.PRJNA1010653.filt, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_Seurat_object_postSCT.RData"))

ovary.Goat.PRJNA1010653.filt <- RunPCA(ovary.Goat.PRJNA1010653.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Goat_PRJNA1010653_ElbowPlot.pdf"))
ElbowPlot(ovary.Goat.PRJNA1010653.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Goat.PRJNA1010653.filt[["pca"]]@stdev / sum(ovary.Goat.PRJNA1010653.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 40

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 17

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 17

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_Goat_PRJNA1010653_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Goat.PRJNA1010653.filt <- RunUMAP(ovary.Goat.PRJNA1010653.filt, dims = 1:pcs)
ovary.Goat.PRJNA1010653.filt <- FindNeighbors(ovary.Goat.PRJNA1010653.filt, dims = 1:pcs)
ovary.Goat.PRJNA1010653.filt <- FindClusters(object = ovary.Goat.PRJNA1010653.filt, resolution = 2)

pdf(paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_UMAP_SeuratClustering_res_2.0.pdf"), width = 7, height = 5)
DimPlot(ovary.Goat.PRJNA1010653.filt, reduction = "umap")
dev.off()

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
# table(ovary.Goat.PRJNA1010653.filt$Library)

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
ovary.Goat.PRJNA1010653.filt.list <- SplitObject(ovary.Goat.PRJNA1010653.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.Goat.PRJNA1010653.filt.list, ncol))))/100
# pred.dblt.rate

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.Goat.PRJNA1010653.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.Goat.PRJNA1010653 <- paramSweep_v3(ovary.Goat.PRJNA1010653.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.Goat.PRJNA1010653    <- summarizeSweep(sweep.res.list.ovary.Goat.PRJNA1010653, GT = FALSE)
  bcmvn.ovary.Goat.PRJNA1010653          <- find.pK(sweep.stats.ovary.Goat.PRJNA1010653)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.Goat.PRJNA1010653 <- as.numeric(as.character(bcmvn.ovary.Goat.PRJNA1010653[as.numeric(bcmvn.ovary.Goat.PRJNA1010653$pK[bcmvn.ovary.Goat.PRJNA1010653$BCmetric == max(bcmvn.ovary.Goat.PRJNA1010653$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.Goat.PRJNA1010653.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.Goat.PRJNA1010653.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.Goat.PRJNA1010653.filt.list[[i]] <- doubletFinder_v3(ovary.Goat.PRJNA1010653.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.Goat.PRJNA1010653, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.Goat.PRJNA1010653.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Goat.PRJNA1010653.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.Goat.PRJNA1010653.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.Goat.PRJNA1010653.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.Goat.PRJNA1010653.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_Goat_PRJNA1010653",names(ovary.Goat.PRJNA1010653.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.Goat.PRJNA1010653.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.Goat.PRJNA1010653.DFsinglets <- merge(ovary.Goat.PRJNA1010653.filt.list[[1]],
                                            y = c(ovary.Goat.PRJNA1010653.filt.list[[2]]),
                                            project = "ovary.Goat.PRJNA1010653")
# ovary.Goat.PRJNA1010653.DFsinglets

# remove pANN columns that are 10xGenomics library lane specific
ovary.Goat.PRJNA1010653.DFsinglets@meta.data <- ovary.Goat.PRJNA1010653.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.Goat.PRJNA1010653.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.Goat.PRJNA1010653.filt.list.scds <- lapply(ovary.Goat.PRJNA1010653.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.Goat.PRJNA1010653.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.Goat.PRJNA1010653.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.Goat.PRJNA1010653.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.Goat.PRJNA1010653.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.Goat.PRJNA1010653.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.Goat.PRJNA1010653.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.Goat.PRJNA1010653.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.Goat.PRJNA1010653.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.Goat.PRJNA1010653.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_Goat_PRJNA1010653", names(ovary.Goat.PRJNA1010653.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.Goat.PRJNA1010653.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.Goat.PRJNA1010653.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.Goat.PRJNA1010653.DFsinglets@meta.data[colnames(ovary.Goat.PRJNA1010653.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.Goat.PRJNA1010653.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.Goat.PRJNA1010653.DFsinglets@meta.data$DoubletFinder, ovary.Goat.PRJNA1010653.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet      26     268
# Singlet     122    5393

ovary.Goat.PRJNA1010653.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.Goat.PRJNA1010653.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.Goat.PRJNA1010653.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                                   "Doublet", "Singlet")
table(ovary.Goat.PRJNA1010653.DFsinglets@meta.data$DoubletCall)
# Doublet Singlet 
#    416    5393 

# re-run dimensionality reduction for plotting purposes
ovary.Goat.PRJNA1010653.DFsinglets <- SCTransform(object = ovary.Goat.PRJNA1010653.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))
ovary.Goat.PRJNA1010653.DFsinglets <- RunPCA(ovary.Goat.PRJNA1010653.DFsinglets, npcs = 50)
ovary.Goat.PRJNA1010653.DFsinglets <- RunUMAP(ovary.Goat.PRJNA1010653.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_DoubletCall_UMAP.pdf"))
DimPlot(ovary.Goat.PRJNA1010653.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.Goat.PRJNA1010653.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.Goat.PRJNA1010653.DFsinglets   <- subset(ovary.Goat.PRJNA1010653.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets
# ovary.Goat.PRJNA1010653.DFsinglets

# table(ovary.Goat.PRJNA1010653.DFsinglets@meta.data$Library)

# save filtered/annotated object
save(ovary.Goat.PRJNA1010653.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_Seurat_object_SINGLETS.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Goat_PRJNA1010653_processing_session_info.txt", sep =""))
sessionInfo()
sink()
