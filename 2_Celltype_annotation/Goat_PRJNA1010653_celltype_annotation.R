options(stringsAsFactors = FALSE)

library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scGate)
library(scProportionTest)
library(scran)
library(SingleR)
library(scSorter)
library(tidyverse)


################################################################################
# Integrative ovarian aging analysis
# Goat PRJNA1010653 dataset
# Annotate cell types
################################################################################

################################################################################
# 1. Load dataset
################################################################################

load("/Volumes/OIProject_II/1_R/1_Pre-processing/Goat/Goat_PRJNA1010653/1_Preprocessing/2025-06-19_10x_ovary_Goat_PRJNA1010653_Seurat_object_SINGLETS.RData")

################################################################################
# 2. Ptprc positive vs. Ptprc negative cells: scGate
################################################################################

# Manually define a simple scGate gating model to purify Ptprc positive cells as immune cells
my.scGate.model.Ptprc <- gating_model(name = "PtprcPos", signature = c("PTPRC"))

#scGate it!
my.ovary.Goat.Ptprc.pos <- scGate(data = ovary.Goat.PRJNA1010653.DFsinglets, model = my.scGate.model.Ptprc)

t.pure <- as.matrix(table(my.ovary.Goat.Ptprc.pos$is.pure, my.ovary.Goat.Ptprc.pos$seurat_clusters))

# Calculate the total cells in each cluster
total <- t.pure[1,] + t.pure[2,]

# Get cluster numbers where Pure cells are more than 80% of total cells in each cluster
high_purity_clusters <- which(t.pure[1,] / total > 0.8)
high_purity_cluster_numbers <- as.numeric(colnames(t.pure)[high_purity_clusters])

# Print cluster numbers
print(high_purity_cluster_numbers)
#[1] [1] 10 14 19 20  4

# Assign cell type to each seurat cluster: Ptprc (CD45) positive cells clusters

my.ovary.Goat.Ptprc.pos <- SetIdent(my.ovary.Goat.Ptprc.pos, value = "seurat_clusters")
ovary.Goat <- my.ovary.Goat.Ptprc.pos

ovary.Goat@meta.data$celltype.level1 <- rep("Ptprc.neg",dim(ovary.Goat@meta.data)[1])
ovary.Goat@meta.data$celltype.level1[ovary.Goat@meta.data$seurat_clusters %in% high_purity_cluster_numbers] <- "Ptprc.pos"

save(ovary.Goat, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_PTPRC_gated_Seurat_object.RData"))

################################################################################
# 2.Ptprc positive sub-celltype annotation: 
#                                           1. Manual marker annotation
#                                           2. SingleR / scType / scSorter
################################################################################

################################
# Subset Ptprc positive cells and re-cluster
################################

ovary.Goat.Ptprc.pos <- subset(ovary.Goat, subset = celltype.level1 == "Ptprc.pos")

DefaultAssay(ovary.Goat.Ptprc.pos) <- "RNA"

########## Normalize, SCT, dimreduc data ##########

ovary.Goat.Ptprc.pos <- NormalizeData(ovary.Goat.Ptprc.pos, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Goat.Ptprc.pos <- SCTransform(object = ovary.Goat.Ptprc.pos, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Goat.Ptprc.pos, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_Ptprc_pos_cells_Seurat_object_clean_postSCT.RData"))

ovary.Goat.Ptprc.pos <- RunPCA(ovary.Goat.Ptprc.pos, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Goat_PRJNA1010653_Ptprc_pos_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.Goat.Ptprc.pos, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Goat.Ptprc.pos[["pca"]]@stdev / sum(ovary.Goat.Ptprc.pos[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 43

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 18

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 18

# Based on these metrics, first 18 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_Goat_PRJNA1010653_Ptprc_pos_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Goat.Ptprc.pos <- RunUMAP(ovary.Goat.Ptprc.pos, dims = 1:pcs)
ovary.Goat.Ptprc.pos <- FindNeighbors(ovary.Goat.Ptprc.pos, dims = 1:pcs)
ovary.Goat.Ptprc.pos <- FindClusters(object = ovary.Goat.Ptprc.pos)
ovary.Goat.Ptprc.pos <- FindClusters(object = ovary.Goat.Ptprc.pos, resolution = 1)

################################
# SingleR - Immgen to annotate Ptprc positive cells
################################

DefaultAssay(ovary.Goat.Ptprc.pos) <- "RNA"

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Goat.Ptprc.pos)

immgen <- ImmGenData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))
immgen
# class: SummarizedExperiment 
# dim: 22134 830 
# metadata(0):
#        assays(1): logcounts
# rownames(22134): Zglp1 Vmn2r65 ... Tiparp Kdm1a
# rowData names(0):
#        colnames(830): GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL
# ... GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
# colData names(3): label.main label.fine label.ont

table(immgen$label.main)
#   B cells      B cells, pro         Basophils                DC Endothelial cells       Eosinophils  Epithelial cells       Fibroblasts               ILC       Macrophages        Mast cells         Microglia 
#        79                 1                 6                88                20                 4                25                21                23                79                20                 3 
# Monocytes       Neutrophils          NK cells               NKT        Stem cells     Stromal cells           T cells               Tgd 
#        33                23                38                22                36                 7               231                71

# Drop cell types with less than 10 cells
immgen <- immgen[,immgen$label.main != 'B cells, pro']
immgen <- immgen[,immgen$label.main != 'Basophils']
immgen <- immgen[,immgen$label.main != 'Eosinophils']
immgen <- immgen[,immgen$label.main != 'Microglia']
immgen <- immgen[,immgen$label.main != 'Stromal cells']

my.singler.immgen = SingleR(test = my.SingleCellExperiment.object,
                            ref  = immgen,
                            assay.type.test = 1,
                            labels = immgen$label.main)

# Add Sample IDs to metadata
my.singler.immgen$meta.data$Library       =    ovary.Goat.Ptprc.pos@meta.data$Library                       # Library IDs

save(my.singler.immgen, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Ptprc_pos_cells_SingleR_object_Immgen.RData",sep = "_"))

###########################
# Transfer cell annotations to Seurat object

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Goat.Ptprc.pos)

# Plot annotation heatmap
annotation_col = data.frame(Labels = factor(my.singler.immgen$labels),
                            Library = as.data.frame(colData(my.SingleCellExperiment.object)[,"Library",drop=FALSE]))

# Transfer SingleR annotations to Seurat Object
ovary.Goat.Ptprc.pos[["SingleR.Immgen"]] <- my.singler.immgen$labels

###########################
# Marker gene expression plots

DefaultAssay(ovary.Goat.Ptprc.pos) <- "SCT"

ovary.Goat.Ptprc.pos <- SetIdent(ovary.Goat.Ptprc.pos, value = "SingleR.Immgen")
levels(ovary.Goat.Ptprc.pos) <- c("Neutrophils", "Macrophages" , "Monocytes", "DC" , "NK cells" ,  "ILC",  "NKT", "T cells",  "Tgd", "B cells", "Endothelial cells", "Fibroblasts", "Mast cells", "Epithelial cells")
DotPlot(ovary.Goat.Ptprc.pos, features = c("Csf3r", "Itgam", "Adgre1", "Itgax", "Klrb1c", "Ccl5", "Itga1", "Gata3", "Il17rb", "Ly6c2", "Cd3e", "Cd8b1", "Cd4", "Cd28", "Tmem176b", "Il7r", "Il2ra", "Cd19"))

table.SingleR.annotation <- table(ovary.Goat.Ptprc.pos@meta.data$seurat_clusters, ovary.Goat.Ptprc.pos@meta.data$SingleR.Immgen)

################################
# scSorter to annotate Ptprc positive cells
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
################################

# Major cell types: Neutrophils, Macrophages, Monocytes, DC, NK cells, ILC, NKT, T cells (CD8 & CD4), B cells
# Use Immgen dataset to curate marker gene list

immgen
# class: SummarizedExperiment 
# dim: 22134 830 ?SingleCellExperiment
# metadata(0):
#   assays(1): logcounts
# rownames(22134): Zglp1 Vmn2r65 ... Tiparp Kdm1a
# rowData names(0):
#   colnames(830): GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL ...
# GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
# colData names(3): label.main label.fine label.ont

se.immgen <- as(immgen, "SingleCellExperiment")
seurat.immgen <- as.Seurat(se.immgen, counts = NULL, data = "logcounts")
seurat.immgen <- SetIdent(seurat.immgen, value = "label.main")

seurat.immgen <- FindVariableFeatures(seurat.immgen, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.immgen)
seurat.immgen <- ScaleData(seurat.immgen, features = all.genes)

Neutrophil.markers <- FindMarkers(seurat.immgen, ident.1 = "Neutrophils", min.pct = 0.25, only.pos = TRUE)
Macrophages.markers <- FindMarkers(seurat.immgen, ident.1 = "Macrophages", min.pct = 0.25, only.pos = TRUE)
Monocytes.markers <- FindMarkers(seurat.immgen, ident.1 = "Monocytes", min.pct = 0.25, only.pos = TRUE)
DC.markers <- FindMarkers(seurat.immgen, ident.1 = "DC", min.pct = 0.25, only.pos = TRUE)
NK.markers <- FindMarkers(seurat.immgen, ident.1 = "NK cells", min.pct = 0.25, only.pos = TRUE)
ILC.markers <- FindMarkers(seurat.immgen, ident.1 = "ILC", min.pct = 0.25, only.pos = TRUE)
NKTcell.markers <- FindMarkers(seurat.immgen, ident.1 = "NKT", min.pct = 0.25, only.pos = TRUE)
Tcell.markers <- FindMarkers(seurat.immgen, ident.1 = "T cells", min.pct = 0.25, only.pos = TRUE)
Bcell.markers <- FindMarkers(seurat.immgen, ident.1 = "B cells", min.pct = 0.25, only.pos = TRUE)

# Curate list of canonical markers for each major cell type
# https://panglaodb.se/markers.html?cell_type=%27NK%20cells%27

my.immune.marker.list <- vector(mode = "list", 9)
names(my.immune.marker.list) <- c("Neutrophils", "Macrophages", "Monocytes", "DC", "NK cells", "ILC", "NKT", "T cells", "B cells")

my.immune.marker.list$Neutrophils <- unique(c(rownames(head(Neutrophil.markers, n = 20)), "Csf3r", "Trem1", "Ptgs2"))     #Ly6g, S100a8 missing from dataset
my.immune.marker.list$Macrophages <- unique(c(rownames(head(Macrophages.markers, n = 20)), "Cd68", "Fcgr1", "Naaa"))
my.immune.marker.list$Monocytes <- unique(c(rownames(head(Monocytes.markers, n = 20)), "Cfp", "Cd14", "Rhoc"))
my.immune.marker.list$DC <- unique(c(rownames(head(DC.markers, n = 20)), "Zbtb46", "Itgax", "Lamp3"))
my.immune.marker.list$'NK cells' <- unique(c(rownames(head(NK.markers, n = 20)), "Ly6c2", "Trdc", "Nkg7", "Klrd1"))
my.immune.marker.list$ILC <- unique(c(rownames(head(ILC.markers, n = 20)), "Il7r", "Gata3"))
my.immune.marker.list$NKT <- unique(c(rownames(head(NKTcell.markers, n = 20)), "Ncam1", "Il2rb", "Cd44"))
my.immune.marker.list$'T cells' <- unique(c(rownames(head(Tcell.markers, n = 20)), "Cd3e", "Trbc2", "Cd3d", "Cd3g"))
my.immune.marker.list$'B cells' <- unique(c(rownames(head(Bcell.markers, n = 20)), "Pxk", "Ms4a1", "Cd75"))

# Filter genes that do not exist in dataset

my.immune.marker.list$Neutrophils <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$Neutrophils)
my.immune.marker.list$Macrophages <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$Macrophages)
my.immune.marker.list$Monocytes <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$Monocytes)
my.immune.marker.list$DC <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$DC)
my.immune.marker.list$'NK cells' <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$'NK cells')
my.immune.marker.list$ILC <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$ILC)
my.immune.marker.list$NKT <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$NKT)
my.immune.marker.list$'T cells' <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$'T cells')
my.immune.marker.list$'B cells' <- intersect(rownames(ovary.Goat.Ptprc.pos), my.immune.marker.list$'B cells')

### Generate anno table
anno<- data.frame(cbind(rep(names(my.immune.marker.list)[1],length(my.immune.marker.list[[1]])), my.immune.marker.list[[1]], rep(2,length(my.immune.marker.list[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.immune.marker.list)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.immune.marker.list)[i],length(my.immune.marker.list[[i]])),
                      "Marker" = my.immune.marker.list[[i]],
                      "Weight" = rep(2,length(my.immune.marker.list[[i]]))))
}

# Manually increase weight for canonical markers for well defined cell types
anno[anno$Marker == "Csf3r",]$Weight <- 6 # Neutrophils
anno[anno$Marker == "Trem1",]$Weight <- 6 # Neutrophils

anno[anno$Marker == "Cd68",]$Weight <- 6 # Macrophages

anno[anno$Marker == "Cfp",]$Weight <- 6 # Monocytes
anno[anno$Marker == "Cd14",]$Weight <- 6 # Monocytes

anno[anno$Marker == "Zbtb46",]$Weight <- 6 # DC

anno[anno$Marker == "Ly6c2",]$Weight <- 6 # NK cells
anno[anno$Marker == "Trdc",]$Weight <- 6 # NK cells
anno[anno$Marker == "Nkg7",]$Weight <- 6 # NK cells

anno[anno$Marker == "Il7r",]$Weight <- 6 # ILC
anno[anno$Marker == "Gata3",]$Weight <- 6 # ILC

anno[anno$Marker == "Cd3e",]$Weight <- 6 # NKT
anno[anno$Marker == "Ly6c2",]$Weight <- 6 # NKT

anno[anno$Marker == "Cd3e",]$Weight <- 6 # T cells
anno[anno$Marker == "Trbc2",]$Weight <- 6 # T cells

anno[anno$Marker == "Cd19",]$Weight <- 6 # B cells
anno[anno$Marker == "Pxk",]$Weight <- 6 # B cells

# Pre-process data for scSorter

DefaultAssay(ovary.Goat.Ptprc.pos) <- "SCT"

topgenes <- head(VariableFeatures(ovary.Goat.Ptprc.pos), 2500)
expr <- GetAssayData(ovary.Goat.Ptprc.pos)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

# Transfer scSorter annotations to Seurat Object
ovary.Goat.Ptprc.pos[["scSorter.Immgen"]] <- rts$Pred_Type

table.scSorter.annotation <- table(ovary.Goat.Ptprc.pos@meta.data$seurat_clusters, ovary.Goat.Ptprc.pos@meta.data$scSorter.Immgen)

################################
# sc-type to annotate Ptprc positive cells
# https://github.com/IanevskiAleksandr/sc-type
################################

# Use same marker list as for scSorter
my.immune.marker.list

cell_markers.tmp = data.frame("cellName" = names(my.immune.marker.list),
                              "geneSymbolmore1" = rep("",length(my.immune.marker.list)),
                              "geneSymbolmore2" = rep("",length(my.immune.marker.list)))

# grab positive markers automatically
for (i in 1:length(my.immune.marker.list)) {
  cell_markers.tmp$geneSymbolmore1[i] <- paste0(my.immune.marker.list[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

# manually edit negative markers
neg.markers <- vector(mode = "list", length = length(my.immune.marker.list))
names(neg.markers) <- names(my.immune.marker.list)

neg.markers[names(neg.markers) != "Neutrophils" ] <- lapply(neg.markers[names(neg.markers) != "Neutrophils" ], function (x){append(x,"Mpx")}) 
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Zap70")}) # Tang, JEM 2017
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Cd4")}) 
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Cd8a")})
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Cd8b1")})
neg.markers[names(neg.markers) != "B cells"     ] <- lapply(neg.markers[names(neg.markers) != "B cells"     ], function (x){append(x,"Cd79a")}) 
neg.markers[names(neg.markers) != "Macrophages" ] <- lapply(neg.markers[names(neg.markers) != "Macrophages" ], function (x){append(x,"Csf1r")})
neg.markers[names(neg.markers) != c("Neutrophils", "Macrophages")] <- lapply(neg.markers[names(neg.markers) != c("Neutrophils", "Macrophages")], function (x){append(x, "Lyz1")}) #lyz cannot be outside of macrophages and neutrophils

# populate negative markers
for (i in 1:length(my.immune.marker.list)) {
  cell_markers.tmp$geneSymbolmore2[i] <- paste0(neg.markers[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

###################################################################################################
######  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   #####
##### Modify code from "gene_sets_prepare" of scType to marker database in R
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

source("/Users/minhookim/Programs/sc-type-1.0/R/sctype_score_.R")

gene_sets_prepare_mOvary <- function(cell_markers.tmp){
  cell_markers = cell_markers.tmp
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

# prepare gene sets
gs_list = gene_sets_prepare_mOvary(cell_markers.tmp)

# Finally, let's assign cell types to each cluster:
# get cell-type by cell matrix
es.max  = sctype_score(scRNAseqData = ovary.Goat.Ptprc.pos[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(ovary.Goat.Ptprc.pos@meta.data$SCT_snn_res.2), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ovary.Goat.Ptprc.pos@meta.data[ovary.Goat.Ptprc.pos@meta.data$SCT_snn_res.2==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ovary.Goat.Ptprc.pos@meta.data$SCT_snn_res.2==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])

ovary.Goat.Ptprc.pos@meta.data$scType_raw = ""

for(i in 1:ncol(es.max)){
  ovary.Goat.Ptprc.pos@meta.data$scType_raw[rownames(ovary.Goat.Ptprc.pos@meta.data) %in% colnames(es.max)[i] ] = as.character(rownames(es.max)[which.max(es.max[,i])])
}

table(ovary.Goat.Ptprc.pos@meta.data$scType_raw)
# B cells          DC         ILC Macrophages   Monocytes Neutrophils    NK cells         NKT     T cells 
#     299         110         363         159          82          57         523         366         813 

table.scType.annotation <- table(ovary.Goat.Ptprc.pos@meta.data$seurat_clusters, ovary.Goat.Ptprc.pos@meta.data$scType_raw)

ovary.Goat.Ptprc.pos <- SetIdent(ovary.Goat.Ptprc.pos, value = "scType_raw")
levels(ovary.Goat.Ptprc.pos) <- c("Neutrophils", "Macrophages" , "Monocytes", "DC" , "NK cells" ,  "ILC",  "NKT", "T cells",  "B cells")
DotPlot(ovary.Goat.Ptprc.pos, features = c("Ptprc", "Csf3r", "Itgam", "Adgre1", "Itgax", "Klrb1c", "Ccl5", "Itga1", "Gata3", "Il17rb", "Ly6c2", "Cd3e", "Cd8b1", "Cd4", "Cd28", "Tmem176b", "Il7r", "Il2ra", "Cd19"))

save(ovary.Goat.Ptprc.pos, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Ptprc_positive_Seurat_object_with_cell_type_annotations_intermediate.RData",sep = "_"))

################################
# Assess cell type annotation by SingleR, scSorter and scType
################################

# SingleR

table.SingleR.annotation <- as.data.frame.matrix(table.SingleR.annotation)

# Manipulate data structure

df_long <- table.SingleR.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.SingleR <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# scSorter

table.scSorter.annotation <- as.data.frame.matrix(table.scSorter.annotation)

# Manipulate data structure

df_long <- table.scSorter.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.scSorter <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# scType

table.scType.annotation <- as.data.frame.matrix(table.scType.annotation)

# Manipulate data structure

df_long <- table.scType.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.scType <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# Combine results

# Check if clusters are in the same order
check_order <- identical(df_max.SingleR$cluster, df_max.scSorter$cluster) && 
  identical(df_max.SingleR$cluster, df_max.scType$cluster)

print(check_order)    # TRUE

annotation.results <- cbind(df_max.SingleR, df_max.scSorter, df_max.scType)
annotation.results <- annotation.results[,c(1,2,3,5,6,8,9)]
colnames(annotation.results) <- c("Cluster", "SingleR_Celltype", "SingleR_Cellcount", "scSorter_Celltype", "scSorter_Cellcount", "scType_Celltype", "scType_Cellcount")

save(annotation.results, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Ptprc_pos_cells_annotation_results.RData",sep = "_"))
write.table(annotation.results, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Ptprc_pos_cells_annotation_results.txt",sep = "_"), quote = FALSE)

################################
# Transfer final annotations
# Refer to ./2024-10-15_10x_ovary_Goat_PRJNA1010653_Ptprc_positive_cells_annotation_results_FINAL.txt
################################

ovary.Goat.Ptprc.pos@meta.data$celltype.level2 <- rep("NA",dim(ovary.Goat.Ptprc.pos@meta.data)[1])
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 0] <- "CD4T"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 1] <- "DNT"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 2] <- "CD8T"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 3] <- "CD4T"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 4] <- "Myeloid"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 5] <- "CD8T"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 6] <- "Myeloid"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 7] <- "CD4T"
ovary.Goat.Ptprc.pos@meta.data$celltype.level2[ovary.Goat.Ptprc.pos@meta.data$seurat_clusters %in% 8] <- "Myeloid"

save(ovary.Goat.Ptprc.pos, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Ptprc_positive_cells_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 3.Ptprc negative sub-celltype annotation: 
#                                           1. Manual marker annotation
#                                           2. SingleR / scType / scSorter
################################################################################

################################
# Subset Ptprc positive cells and re-cluster
################################

ovary.Goat.Ptprc.neg <- subset(ovary.Goat, subset = celltype.level1 == "Ptprc.neg")

DefaultAssay(ovary.Goat.Ptprc.pos) <- "RNA"

########## Normalize, SCT, dimreduc data ##########

ovary.Goat.Ptprc.neg <- NormalizeData(ovary.Goat.Ptprc.neg, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Goat.Ptprc.neg <- SCTransform(object = ovary.Goat.Ptprc.neg, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Goat.Ptprc.neg, file = paste0(Sys.Date(),"_10x_ovary_AC_Ptprc_neg_cells_Seurat_object_clean_postSCT.RData"))

ovary.Goat.Ptprc.neg <- RunPCA(ovary.Goat.Ptprc.neg, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_AC_Ptprc_neg_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.Goat.Ptprc.neg, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Goat.Ptprc.neg[["pca"]]@stdev / sum(ovary.Goat.Ptprc.neg[["pca"]]@stdev) * 100

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
pdf(paste0(Sys.Date(), "_10x_ovary_AC_Ptprc_pos_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Goat.Ptprc.neg <- RunUMAP(ovary.Goat.Ptprc.neg, dims = 1:pcs)
ovary.Goat.Ptprc.neg <- FindNeighbors(ovary.Goat.Ptprc.neg, dims = 1:pcs)
ovary.Goat.Ptprc.neg <- FindClusters(object = ovary.Goat.Ptprc.neg)
ovary.Goat.Ptprc.neg <- FindClusters(object = ovary.Goat.Ptprc.neg, resolution = 1.5)

# Annotate cells

ovary.Goat.Ptprc.neg <- SetIdent(ovary.Goat.Ptprc.neg, value = "seurat_clusters")

VlnPlot(ovary.Goat.Ptprc.neg, features = c("nCount_RNA"))
VlnPlot(ovary.Goat.Ptprc.neg, features = c("nFeature_RNA"))
VlnPlot(ovary.Goat.Ptprc.neg, features = c("percent.mito"))
VlnPlot(ovary.Goat.Ptprc.neg, features = c("decontX_contamination"))

ovary.Goat.Ptprc.neg@meta.data$celltype.level2 <- rep("NA",dim(ovary.Goat.Ptprc.neg@meta.data)[1])
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 0] <- "Theca"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 1] <- "SMC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 2] <- "SMC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 3] <- "Epithelial"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 4] <- "BEC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 5] <- "Stroma"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 6] <- "BEC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 7] <- "Stroma"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 8] <- "Stroma"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 9] <- "Granulosa"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 10] <- "Stroma"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 11] <- "LEC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 12] <- "Theca"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 13] <- "Theca"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 14] <- "SMC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 15] <- "BEC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 16] <- "Granulosa"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 17] <- "BEC"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 18] <- "Epithelial"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 19] <- "Granulosa"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 20] <- "Granulosa"
ovary.Goat.Ptprc.neg@meta.data$celltype.level2[ovary.Goat.Ptprc.neg@meta.data$seurat_clusters %in% 21] <- "SMC"

save(ovary.Goat.Ptprc.neg, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Ptprc_negative_cells_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 4.Combine cell type annotation
################################################################################

# Extract celltype.level2 metadata
celltype.level2.pos <- ovary.Goat.Ptprc.pos@meta.data[, "celltype.level2", drop = FALSE]
celltype.level2.neg <- ovary.Goat.Ptprc.neg@meta.data[, "celltype.level2", drop = FALSE]

# Combine the metadata from both subsetted objects
combined.celltype.level2 <- rbind(celltype.level2.pos, celltype.level2.neg)

# Ensure that the cell IDs match
common_cell_ids <- intersect(rownames(ovary.Goat@meta.data), rownames(combined.celltype.level2))

# Update the celltype.level2 in ovary.Goat for the matching cells
ovary.Goat@meta.data[common_cell_ids, "celltype.level2"] <- combined.celltype.level2[common_cell_ids, "celltype.level2"]

# Verify that the metadata has been successfully updated
table(is.na(ovary.Goat@meta.data$celltype.level2))
# FALSE 
# 7911 

save(ovary.Goat, file = paste(Sys.Date(),"10x_ovary_Goat_PRJNA1010653_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Goat_PRJNA1010653_cell_type_annotation_session_Info.txt", sep =""))
sessionInfo()
sink()
