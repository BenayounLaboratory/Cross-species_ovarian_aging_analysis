options(stringsAsFactors = F)

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
# Mouse EMTAB12889 dataset
# Annotate cell types
################################################################################

################################################################################
# 1. Load data
################################################################################

load("/Volumes/OIProject_II/1_R/Mouse_EMTAB12889/1_Preprocessing/2024-02-26_10x_ovary_Mouse_EMTAB12889_processed_clean_Seurat_object_FINAL.RData")

################################################################################
# 2. Gate CD45+ cells (Ptprc expressing cells)
################################################################################

#Manually define a simple scGate gating model to purify immune cells using a positive marker Ptprc
my.scGate.model.immune <- gating_model(name = "immune", signature = c("Ptprc"))

# scGate it!
# Gate individual batches separately: https://github.com/carmonalab/scGate/issues/12
ovary.Ptprc.pos <- scGate(data = ovary.EMTAB12889.cl, model = my.scGate.model.immune, assay = "RNA", verbose = TRUE)

t.pure <- as.matrix(table(ovary.Ptprc.pos$is.pure, ovary.Ptprc.pos$seurat_clusters))

# Calculate the total cells in each cluster
total <- t.pure[1,] + t.pure[2,]

# Get cluster numbers where Pure cells are more than 80% of total cells in each cluster
high_purity_clusters <- which(t.pure[1,] / total > 0.80)
high_purity_cluster_numbers <- as.numeric(colnames(t.pure)[high_purity_clusters])

# Print cluster numbers
print(high_purity_cluster_numbers)
#  [1]  2 13 17 19 27 29 40 41 45 47 52 54 63

# Assign cell type to each seurat cluster: Ptprc (CD45) positive cells clusters:   [1]  2 13 17 19 27 29 40 41 45 47 52 54 63

ovary.Ptprc.pos <- SetIdent(ovary.Ptprc.pos, value = "seurat_clusters")
ovary.EMTAB12889 <- ovary.Ptprc.pos

ovary.EMTAB12889@meta.data$celltype.level1 <- rep("Ptprc.neg",dim(ovary.EMTAB12889@meta.data)[1])
ovary.EMTAB12889@meta.data$celltype.level1[ovary.EMTAB12889@meta.data$seurat_clusters %in% high_purity_cluster_numbers] <- "Ptprc.pos"

save(ovary.EMTAB12889, file = paste0(Sys.Date(),"_10x_ovary_mouse-EMTAB12889_Ptprc_gated_Seurat_object.RData"))

################################################################################
# 2.Ptprc positive sub-celltype annotation: 
#                                           1. Manual marker annotation
#                                           2. SingleR / scType / scSorter
################################################################################

################################
# Subset Ptprc positive cells and re-cluster
################################

ovary.Ptprc.pos <- subset(ovary.EMTAB12889, subset = celltype.level1 == "Ptprc.pos")

DefaultAssay(ovary.Ptprc.pos) <- "RNA"

########## Normalize, SCT, dimreduc data ##########

ovary.Ptprc.pos <- NormalizeData(ovary.Ptprc.pos, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Ptprc.pos <- SCTransform(object = ovary.Ptprc.pos, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Ptprc.pos, file = paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_pos_cells_Seurat_object_clean_postSCT.RData"))

ovary.Ptprc.pos <- RunPCA(ovary.Ptprc.pos, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_mouse_EMTAB12889Ptprc_pos_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.Ptprc.pos, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Ptprc.pos[["pca"]]@stdev / sum(ovary.Ptprc.pos[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 13

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 13

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_mouse_EMTAB12889Ptprc_pos_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Ptprc.pos <- RunUMAP(ovary.Ptprc.pos, dims = 1:pcs)
ovary.Ptprc.pos <- FindNeighbors(ovary.Ptprc.pos, dims = 1:pcs)
ovary.Ptprc.pos <- FindClusters(object = ovary.Ptprc.pos)
ovary.Ptprc.pos <- FindClusters(object = ovary.Ptprc.pos, resolution = 3)

################################
# SingleR - Immgen to annotate Ptprc positive cells
################################

DefaultAssay(ovary.Ptprc.pos) <- "RNA"

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Ptprc.pos)

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

table(immgen$label.main)
#    B cells                DC Endothelial cells       Fibroblasts               ILC       Macrophages 
#        469               278                 2                 4              1801               313 
# Mast cells         Monocytes       Neutrophils          NK cells               NKT           T cells 
#          9               273               109               242              1333              1086 
#        Tgd 
#         12 

my.singler.immgen = SingleR(test = my.SingleCellExperiment.object,
                            ref  = immgen,
                            assay.type.test = 1,
                            labels = immgen$label.main)

table(my.singler.immgen$labels)
# B cells          DC Fibroblasts         ILC Macrophages  Mast cells   Monocytes Neutrophils    NK cells         NKT     T cells         Tgd 
#.    599         163           1        1713         203          18         512         185         711         837         917          23 

# Add Sample IDs to metadata
my.singler.immgen$meta.data$Library       =    ovary.Ptprc.pos@meta.data$Library                       # Library IDs

save(my.singler.immgen, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_pos_cells_SingleR_object_Immgen.RData",sep = "_"))

###########################
# Transfer cell annotations to Seurat object

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Ptprc.pos)

# Transfer SingleR annotations to Seurat Object
ovary.Ptprc.pos[["SingleR.Immgen"]] <- my.singler.immgen$labels

###########################
# Marker gene expression plots

DefaultAssay(ovary.Ptprc.pos) <- "SCT"

ovary.Ptprc.pos <- SetIdent(ovary.Ptprc.pos, value = "SingleR.Immgen")
levels(ovary.Ptprc.pos) <- c("Neutrophils", "Macrophages" , "Monocytes", "DC" , "NK cells" ,  "ILC",  "NKT", "T cells",  "Tgd", "B cells", "Endothelial cells", "Fibroblasts", "Mast cells", "Epithelial cells", "Stem cells")
DotPlot(ovary.Ptprc.pos, features = c("Csf3r", "Itgam", "Adgre1", "Itgax", "Klrb1c", "Ccl5", "Itga1", "Gata3", "Il17rb", "Ly6c2", "Cd3e", "Cd8b1", "Cd4", "Cd28", "Tmem176b", "Il7r", "Il2ra", "Cd19"))

table.SingleR.annotation <- table(ovary.Ptprc.pos@meta.data$seurat_clusters, ovary.Ptprc.pos@meta.data$SingleR.Immgen)

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

# rownames(head(Neutrophil.markers, n = 20))

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

my.immune.marker.list$Neutrophils <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$Neutrophils)
my.immune.marker.list$Macrophages <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$Macrophages)
my.immune.marker.list$Monocytes <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$Monocytes)
my.immune.marker.list$DC <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$DC)
my.immune.marker.list$'NK cells' <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$'NK cells')
my.immune.marker.list$ILC <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$ILC)
my.immune.marker.list$NKT <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$NKT)
my.immune.marker.list$'T cells' <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$'T cells')
my.immune.marker.list$'B cells' <- intersect(rownames(ovary.Ptprc.pos), my.immune.marker.list$'B cells')

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

DefaultAssay(ovary.Ptprc.pos) <- "SCT"

topgenes <- head(VariableFeatures(ovary.Ptprc.pos), 2500)
expr <- GetAssayData(ovary.Ptprc.pos)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

# Transfer scSorter annotations to Seurat Object
ovary.Ptprc.pos[["scSorter.Immgen"]] <- rts$Pred_Type

table.scSorter.annotation <- table(ovary.Ptprc.pos@meta.data$seurat_clusters, ovary.Ptprc.pos@meta.data$scSorter.Immgen)

################################
# sc-type to annotate Ptprc positive cells
# https://github.com/IanevskiAleksandr/sc-type
################################

# Use same marker list as for scSorter
my.immune.marker.list
names(my.immune.marker.list)
# [1] "Neutrophils" "Macrophages" "Monocytes"   "DC"          "NK cells"    "ILC"         "NKT"         "T cells"     "B cells"   

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
es.max  = sctype_score(scRNAseqData = ovary.Ptprc.pos[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(ovary.Ptprc.pos@meta.data$SCT_snn_res.4), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ovary.Ptprc.pos@meta.data[ovary.Ptprc.pos@meta.data$SCT_snn_res.4==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ovary.Ptprc.pos@meta.data$SCT_snn_res.4==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])

ovary.Ptprc.pos@meta.data$scType_raw = ""

for(i in 1:ncol(es.max)){
  ovary.Ptprc.pos@meta.data$scType_raw[rownames(ovary.Ptprc.pos@meta.data) %in% colnames(es.max)[i] ] = as.character(rownames(es.max)[which.max(es.max[,i])])
}

table.scType.annotation <- table(ovary.Ptprc.pos@meta.data$seurat_clusters, ovary.Ptprc.pos@meta.data$scType_raw)

save(ovary.Ptprc.pos, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_positive_Seurat_object_with_cell_type_annotations_intermediate.RData",sep = "_"))

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

save(annotation.results, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_pos_cells_annotation_results.RData",sep = "_"))
write.table(annotation.results, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_pos_cells_annotation_results.txt",sep = "_"), quote = FALSE)

# Manual marker gene analysis

ovary.Ptprc.pos <- SetIdent(ovary.Ptprc.pos, value = "seurat_clusters")

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_positive_cells_Vlnplot_marker_genes_clusters_1_to_7.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.pos, subset = seurat_clusters %in% c(0:7)), 
        features = c("S100a8",                                        # S100a8 - Neutrophil marker
                     "Adgre1", "C1qa",                                # Monocyte (Cd68, Ly6c2) & Macrophage (Cd63, C1qa) markers
                     "Flt3",                                          # DC marker
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", # ILC
                     "Ncr1", "Gzma", "Gzmb", "Nkg7",                  # NK markers
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))                  # T & B cell markers
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_positive_cells_Vlnplot_marker_genes_clusters_8_to_14.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.pos, subset = seurat_clusters %in% c(8:14)), 
        features = c("S100a8",                                        # S100a8 - Neutrophil marker
                     "Adgre1", "C1qa",                                # Monocyte (Cd68, Ly6c2) & Macrophage (Cd63, C1qa) markers
                     "Flt3",                                          # DC marker
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", # ILC
                     "Ncr1", "Gzma", "Gzmb", "Nkg7",                  # NK markers
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))                  # T & B cell markers
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_positive_cells_Vlnplot_marker_genes_clusters_15_to_17.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.pos, subset = seurat_clusters %in% c(15:17)), 
        features = c("S100a8",                                        # S100a8 - Neutrophil marker
                     "Adgre1", "C1qa",                                 # Monocyte (Cd68, Ly6c2) & Macrophage (Cd63, C1qa) markers
                     "Flt3",                                          # DC marker
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", # ILC
                     "Ncr1", "Gzma", "Gzmb", "Nkg7",                  # NK markers
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))                  # T & B cell markers
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_positive_cells_Vlnplot_marker_genes_clusters_18_to_21.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.pos, subset = seurat_clusters %in% c(18:21)), 
        features = c("S100a8",                                        # S100a8 - Neutrophil marker
                     "Adgre1", "C1qa",                                 # Monocyte (Cd68, Ly6c2) & Macrophage (Cd63, C1qa) markers
                     "Flt3",                                          # DC marker
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", # ILC
                     "Ncr1", "Gzma", "Gzmb", "Nkg7",                  # NK markers
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))                  # T & B cell markers
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_positive_cells_Vlnplot_marker_genes_clusters_22_to_27.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.pos, subset = seurat_clusters %in% c(22:27)), 
        features = c("S100a8",                                        # S100a8 - Neutrophil marker
                     "Adgre1", "C1qa",                                 # Monocyte (Cd68, Ly6c2) & Macrophage (Cd63, C1qa) markers
                     "Flt3",                                          # DC marker
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", # ILC
                     "Ncr1", "Gzma", "Gzmb", "Nkg7",                  # NK markers
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))                  # T & B cell markers
dev.off()

################################
# Transfer final annotations
################################

ovary.Ptprc.pos@meta.data$celltype.level2 <- rep("NA",dim(ovary.Ptprc.pos@meta.data)[1])
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 0] <- "CD8 NKT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 1] <- "NK"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 2] <- "Macrophage"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 3] <- "ILC"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 4] <- "CD8 NKT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 5] <- "ILC"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 6] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 7] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 8] <- "B"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 9] <- "Macrophage"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 10] <- "Macrophage"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 11] <- "NK"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 12] <- "ILC"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 13] <- "Neutrophil"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 14] <- "B"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 15] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 16] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 17] <- "Macrophage"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 18] <- "Macrophage"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 19] <- "B"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 20] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 21] <- "NK"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 22] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 23] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 24] <- "ILC"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 25] <- "DNT"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 26] <- "Macrophage"
ovary.Ptprc.pos@meta.data$celltype.level2[ovary.Ptprc.pos@meta.data$seurat_clusters %in% 27] <- "B"

save(ovary.Ptprc.pos, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_positive_cells_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 3.Ptprc negative sub-celltype annotation: 
#                                           1. Manual marker annotation
#                                           2. SingleR / scType / scSorter
################################################################################

################################
# Subset Ptprc positive cells and re-cluster
################################

ovary.Ptprc.neg <- subset(ovary.EMTAB12889, subset = celltype.level1 == "Ptprc.neg")

DefaultAssay(ovary.Ptprc.pos) <- "RNA"

########## Normalize, SCT, dimreduc data ##########

ovary.Ptprc.neg <- NormalizeData(ovary.Ptprc.neg, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Ptprc.neg <- SCTransform(object = ovary.Ptprc.neg, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.Ptprc.neg, file = paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889_Ptprc_neg_cells_Seurat_object_clean_postSCT.RData"))

ovary.Ptprc.neg <- RunPCA(ovary.Ptprc.neg, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_mouse_EMTAB12889_Ptprc_neg_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.Ptprc.neg, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Ptprc.neg[["pca"]]@stdev / sum(ovary.Ptprc.neg[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 12

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 12

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_mouse_EMTAB12889Ptprc_pos_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Ptprc.neg <- RunUMAP(ovary.Ptprc.neg, dims = 1:pcs)
ovary.Ptprc.neg <- FindNeighbors(ovary.Ptprc.neg, dims = 1:pcs)
ovary.Ptprc.neg <- FindClusters(object = ovary.Ptprc.neg)
ovary.Ptprc.neg <- FindClusters(object = ovary.Ptprc.neg, resolution = 3)

################################
# SingleR - Public datasets to annotate ovarian cells - PMID:36205477, PMID:38200272
################################

###############################
# PMID:36205477
###############################

# Prep reference file

cycling.ovary.seurat <- readRDS("/Volumes/OIProject_II/1_R/SingleR_Ref/ovary_0.rds")
cycling.ovary.seurat
# An object of class Seurat 
# 23657 features across 34712 samples within 1 assay 
# Active assay: RNA (23657 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

cycling.ovary.seurat@meta.data$celltype <- rep('NA')
cycling.ovary.seurat@meta.data$celltype[grep("GC_", cycling.ovary.seurat@meta.data$Level1)] <- "Granulosa"
cycling.ovary.seurat@meta.data$celltype[grep("M_", cycling.ovary.seurat@meta.data$Level1)] <- "Mesenchyme"
cycling.ovary.seurat@meta.data$celltype[grep("Stroma", cycling.ovary.seurat@meta.data$Level1)] <- "Stroma"
cycling.ovary.seurat@meta.data$celltype[grep("Theca", cycling.ovary.seurat@meta.data$Level1)] <- "Theca"
cycling.ovary.seurat@meta.data$celltype[grep("Smooth muscle", cycling.ovary.seurat@meta.data$Level1)] <- "Smooth muscle"
cycling.ovary.seurat@meta.data$celltype[grep("I_", cycling.ovary.seurat@meta.data$Level1)] <- "Immune"
cycling.ovary.seurat@meta.data$celltype[grep("EN_B", cycling.ovary.seurat@meta.data$Level1)] <- "BEC"
cycling.ovary.seurat@meta.data$celltype[grep("EN_L", cycling.ovary.seurat@meta.data$Level1)] <- "LEC"
cycling.ovary.seurat@meta.data$celltype[grep("Epithelium", cycling.ovary.seurat@meta.data$Level1)] <- "Epithelial"
cycling.ovary.seurat@meta.data$celltype[grep("Oocyte", cycling.ovary.seurat@meta.data$Level1)] <- "Oocyte"

cycling.ovary.seurat.2 <- subset(cycling.ovary.seurat, subset = celltype == "Immune", invert = TRUE)

cycle.data.sce <- as.SingleCellExperiment(cycling.ovary.seurat.2)
dim(cycle.data.sce)
# [1] 23657 33063

# Prep test data
DefaultAssay(ovary.Ptprc.neg) <- "RNA"
my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.Ptprc.neg)

# Find genes in common to the data sets and limit analysis to these
commonGenes <- intersect(rownames(cycle.data.sce), rownames(my.SingleCellExperiment.object))
# length(commonGenes)

my.cycling.ovary.ref <- cycle.data.sce[commonGenes,]
# dim(my.cycling.ovary.ref)

my.SingleCellExperiment.object <- my.SingleCellExperiment.object[commonGenes,]
# dim(my.SingleCellExperiment.object)

my.trained.cycling <- trainSingleR(my.cycling.ovary.ref, my.cycling.ovary.ref$celltype)

my.singleR.cycling <- classifySingleR(test = my.SingleCellExperiment.object, trained = my.trained.cycling)

table(my.singleR.cycling$labels)
#  BEC Epithelial  Granulosa        LEC Mesenchyme     Oocyte     Stroma      Theca 
#  461        214       1643         85        251         31       1036        984 

save(my.trained.cycling, file = paste(Sys.Date(),"10X_ovary_mouse_EMTAB12889_Ptprc_negative_cells_SingleR_cycing_data_trained_classifier.RData",sep = "_"))
save(my.singleR.cycling, file = paste(Sys.Date(),"10X_ovary_mouse_EMTAB12889_Ptprc_negative_cells_SingleR_cycing_data_trained_result.RData",sep = "_"))

# Transfer cell type annotation to Seurat
ovary.Ptprc.neg[["SingleR.cycling"]] <- my.singleR.cycling$labels

pdf(paste(Sys.Date(),"10X_ovary_mouse_EMTAB12889_Ptprc_negative_cells_SingleR_UMAP_annotated_cycling.pdf", sep = "_"), height = 7, width = 8)
DimPlot(ovary.Ptprc.neg, reduction = "umap", group.by = "SingleR.cycling")
dev.off()

###########################
# Marker gene expression plots

DefaultAssay(ovary.Ptprc.neg) <- "SCT"

ovary.Ptprc.neg <- SetIdent(ovary.Ptprc.neg, value = "SingleR.cycling")
levels(ovary.Ptprc.neg) <- c("Oocyte",
                             "Granulosa",
                             "Theca",
                             "Stroma",
                             "Mesenchyme",
                             "BEC",
                             "LEC",
                             "Epithelial")
DotPlot(ovary.Ptprc.neg, features = c("Ddx4",                             # Oocyte marker
                                      "Amh", "Fshr","Foxl2", "Cyp19a1",   # Granulosa cell marker
                                      "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                                      "Pdgfra", "Col1a1",                 # Stroma cell marker
                                      "Acta2", "Rgs5", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                      "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                                      "Pax8", "Cdh1"))                    # Epithelial cell marker

DotPlot(ovary.Ptprc.neg, features = c("Amh", "Fshr","Foxl2", "Nr5a2", "Kitl", "Cyp19a1",     # Granulosa cell marker
                                      "Star",                                                # Star is produced in various ovarian cells (granulosa, theca, epithelial, etc)
                                      "Mdk", "Cyp17a1", "Mgp", "Fbln5", "Cyp11a1", "Lhcgr",  # Theca cell marker
                                      "Col1a1", "Lum"))                                      # Stroma cell marker

# table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$SingleR.cycling)

table.SingleR.cycling.annotation <- table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$SingleR.cycling)

###########################
# Marker gene expression plots

DefaultAssay(ovary.Ptprc.neg) <- "SCT"

ovary.Ptprc.neg <- SetIdent(ovary.Ptprc.neg, value = "SingleR.EMTAB12889")
levels(ovary.Ptprc.neg) <- c("Oocyte",
                             "Granulosa",
                             "Theca",
                             "Stroma",
                             "Endothelial",
                             "Epithelial")
DotPlot(ovary.Ptprc.neg, features = c("Ddx4",                             # Oocyte marker
                                      "Amh", "Fshr","Foxl2", "Cyp19a1",   # Granulosa cell marker
                                      "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                                      "Pdgfra", "Col1a1",                 # Stroma cell marker
                                      "Acta2", "Rgs5", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                      "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                                      "Pax8", "Cdh1"))                    # Epithelial cell marker

DotPlot(ovary.Ptprc.neg, features = c("Amh", "Fshr","Foxl2", "Nr5a2", "Kitl", "Cyp19a1",     # Granulosa cell marker
                                      "Star",                                                # Star is produced in various ovarian cells (granulosa, theca, epithelial, etc)
                                      "Mdk", "Cyp17a1", "Mgp", "Fbln5", "Cyp11a1", "Lhcgr",  # Theca cell marker
                                      "Col1a1", "Lum"))                                      # Stroma cell marker

# table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$SingleR.EMTAB12889)

table.SingleR.EMTAB12889.annotation <- table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$SingleR.EMTAB12889)

################################
# scSorter to annotate Ptprc negative cells
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
################################

# Use curated marker gene list
# Marker genes derived from cycling data + manually added from literature

cycling.ovary.seurat.2 <- NormalizeData(cycling.ovary.seurat.2, normalization.method = "LogNormalize", scale.factor = 10000)
cycling.ovary.seurat.2 <- SCTransform(object = cycling.ovary.seurat.2, vars.to.regress = c("nGene", "fraction.mito", "mouse"))
cycling.ovary.seurat.2 <- RunPCA(cycling.ovary.seurat.2, npcs = 50)
cycling.ovary.seurat.2 <- RunUMAP(cycling.ovary.seurat.2, dims = 1:20)
cycling.ovary.seurat.2 <- FindNeighbors(cycling.ovary.seurat.2, dims = 1:20)
cycling.ovary.seurat.2 <- FindClusters(object = cycling.ovary.seurat.2)

cycling.ovary.seurat.2
# An object of class Seurat 
# 46196 features across 33063 samples within 2 assays 
# Active assay: SCT (22539 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

cycling.ovary.seurat.2 <- SetIdent(cycling.ovary.seurat.2, value = "celltype")
cycling.ovary.markers <- FindAllMarkers(cycling.ovary.seurat.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cycling.ovary.markers.top20 <- cycling.ovary.markers %>% 
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

my.ovarian.marker.list <- vector(mode = "list", 9)
names(my.ovarian.marker.list) <- c("Oocyte",
                                   "Granulosa",
                                   "Theca",
                                   "Stroma",
                                   "Mesenchyme",
                                   "BEC",
                                   "LEC",
                                   "Epithelial")

# Curate marker gene list (top 20 marker genes from cycling data + additional markers)
my.ovarian.marker.list$Oocyte <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Oocyte"], "Ddx4", "DazI"))
my.ovarian.marker.list$Granulosa <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Granulosa"], "Amh", "Fshr", "Cyp19a1", "Foxl2", "Bex1", "Cnmd", "Inha", "Gja1", "Amhr2", "KitI", "Cdh2", "Star"))
my.ovarian.marker.list$Theca <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Theca"], "Top2a", "Lhcgr", "Dcn", "Mdk", "Star", "Cyp11a1", "Cyp17a1", "Col1a2", "Col3a1", "Ogn", "Fbln5"))
my.ovarian.marker.list$Stroma <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Stroma"],  "Mdk", "Nr2f2", "Col6a3", "Lum", "Col1a1", "Pdgfra", "Nr2f2"))
my.ovarian.marker.list$Mesenchyme <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Mesenchyme"],  "Acta2", "Mcam", "Rgs5", "Notch3"))
my.ovarian.marker.list$BEC <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "BEC"],  "Flt1", "Vwf"))
my.ovarian.marker.list$LEC <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "LEC"],  "Flt4", "Prox1"))
my.ovarian.marker.list$Epithelial <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Epithelial"],  "Cldn1", "Pax8"))

save(my.ovarian.marker.list, file = paste(Sys.Date(),"Ovarian_non-immune_cell_markers_for_scSorter.RData",sep = "_"))


# Filter genes that do not exist in dataset

my.ovarian.marker.list$Oocyte <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$Oocyte)
my.ovarian.marker.list$Granulosa <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$Granulosa)
my.ovarian.marker.list$Theca <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$Theca)
my.ovarian.marker.list$Stroma <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$Stroma)
my.ovarian.marker.list$Mesenchyme <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$Mesenchyme)
my.ovarian.marker.list$BEC <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$BEC)
my.ovarian.marker.list$LEC <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$LEC)
my.ovarian.marker.list$Epithelial <- intersect(rownames(ovary.Ptprc.neg), my.ovarian.marker.list$Epithelial)

### Generate anno table
anno<- data.frame(cbind(rep(names(my.ovarian.marker.list)[1],length(my.ovarian.marker.list[[1]])), my.ovarian.marker.list[[1]], rep(2,length(my.ovarian.marker.list[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.ovarian.marker.list)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.ovarian.marker.list)[i],length(my.ovarian.marker.list[[i]])),
                      "Marker" = my.ovarian.marker.list[[i]],
                      "Weight" = rep(2,length(my.ovarian.marker.list[[i]]))))
}

# Manually increase weight for canonical markers for well defined cell types
anno[anno$Marker == "Ddx4",]$Weight <- 6
anno[anno$Marker == "Cyp19a1",]$Weight <- 6
anno[anno$Marker == "Amh",]$Weight <- 6
anno[anno$Marker == "Fshr",]$Weight <- 6
anno[anno$Marker == "Foxl2",]$Weight <- 6
anno[anno$Marker == "Cyp11a1",]$Weight <- 6
anno[anno$Marker == "Cyp17a1",]$Weight <- 6
anno[anno$Marker == "Nr2f2",]$Weight <- 6
anno[anno$Marker == "Col1a1",]$Weight <- 6
anno[anno$Marker == "Rgs5",]$Weight <- 6
anno[anno$Marker == "Flt1",]$Weight <- 6
anno[anno$Marker == "Flt4",]$Weight <- 6

# Pre-process data for scSorter
topgenes <- head(VariableFeatures(ovary.Ptprc.neg), 2500)
expr <- GetAssayData(ovary.Ptprc.neg)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

### Viewing Results - The cell type assignment results are stored in the Pred_Type vector.
print(table(rts$Pred_Type))
#   BEC Epithelial  Granulosa        LEC Mesenchyme     Oocyte     Stroma      Theca 
#   505         940         609         168         604         291         377         934         121         156 

# Transfer scSorter annotations to Seurat Object
ovary.Ptprc.neg[["scSorter.labels"]] <- rts$Pred_Type

# table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$scSorter.labels)

table.scSorter.annotation <- table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$scSorter.labels)

ovary.Ptprc.neg <- SetIdent(ovary.Ptprc.neg, value = "scSorter.labels")
levels(ovary.Ptprc.neg) <- c("Oocyte",
                             "Granulosa",
                             "Theca",
                             "Stroma",
                             "Mesenchyme",
                             "Pericyte",
                             "BEC",
                             "LEC",
                             "Epithelial",
                             "Unknown")
DotPlot(ovary.Ptprc.neg, features = c("Ddx4",                             # Oocyte marker
                                      "Amh", "Fshr","Foxl2", "Cyp19a1",   # Granulosa cell marker
                                      "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                                      "Pdgfra", "Col1a1",                 # Stroma cell marker
                                      "Acta2", "Rgs5", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                      "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                                      "Pax8", "Cdh1"))                    # Epithelial cell marker

DotPlot(ovary.Ptprc.neg, features = c("Amh", "Fshr","Foxl2", "Nr5a2", "Kitl", "Cyp19a1",     # Granulosa cell marker
                                      "Star",                                                # Star is produced in various ovarian cells (granulosa, theca, epithelial, etc)
                                      "Mdk", "Cyp17a1", "Mgp", "Fbln5", "Cyp11a1", "Lhcgr",  # Theca cell marker
                                      "Col1a1", "Lum"))                                      # Stroma cell marker

################################
# sc-type to annotate Ptprc negative cells
# https://github.com/IanevskiAleksandr/sc-type
################################

# Use same marker list as for scSorter
my.ovarian.marker.list
names(my.ovarian.marker.list)
# [1] "Oocyte"     "Granulosa"  "Theca"      "Stroma"     "Mesenchyme" "BEC"        "LEC"        "Epithelial"

cell_markers.tmp = data.frame("cellName" = names(my.ovarian.marker.list),
                              "geneSymbolmore1" = rep("",length(my.ovarian.marker.list)),
                              "geneSymbolmore2" = rep("",length(my.ovarian.marker.list)))

# grab positive markers automatically
for (i in 1:length(my.ovarian.marker.list)) {
  cell_markers.tmp$geneSymbolmore1[i] <- paste0(my.ovarian.marker.list[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

# manually edit negative markers
neg.markers <- vector(mode = "list", length = length(my.ovarian.marker.list))
names(neg.markers) <- names(my.ovarian.marker.list)

# populate negative markers
for (i in 1:length(my.ovarian.marker.list)) {
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
es.max  = sctype_score(scRNAseqData = ovary.Ptprc.neg[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(ovary.Ptprc.neg@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ovary.Ptprc.neg@meta.data[ovary.Ptprc.neg@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ovary.Ptprc.neg@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])

ovary.Ptprc.neg@meta.data$scType_raw = ""

for(i in 1:ncol(es.max)){
  ovary.Ptprc.neg@meta.data$scType_raw[rownames(ovary.Ptprc.neg@meta.data) %in% colnames(es.max)[i] ] = as.character(rownames(es.max)[which.max(es.max[,i])])
}

table.scType.annotation <- table(ovary.Ptprc.neg@meta.data$seurat_clusters, ovary.Ptprc.neg@meta.data$scType_raw)

save(ovary.Ptprc.neg, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_negativetive_Seurat_object_with_cell_type_annotations_intermediate.RData",sep = "_"))

################################
# Assess cell type annotation by SingleR, scSorter and scType
################################

# SingleR - cycling

table.SingleR.cycling.annotation <- as.data.frame.matrix(table.SingleR.cycling.annotation)

# Manipulate data structure

df_long <- table.SingleR.cycling.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.SingleR.cycling <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# SingleR - EMTAB12889

table.SingleR.EMTAB12889.annotation <- as.data.frame.matrix(table.SingleR.EMTAB12889.annotation)

# Manipulate data structure

df_long <- table.SingleR.EMTAB12889.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.SingleR.EMTAB12889 <- df_count %>%
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
check_order <- identical(df_max.SingleR.cycling$cluster, df_max.SingleR.EMTAB12889$cluster) && 
  identical(df_max.SingleR.cycling$cluster, df_max.scSorter$cluster) &&
  identical(df_max.SingleR.cycling$cluster, df_max.scType$cluster)

print(check_order)    # TRUE

annotation.results <- cbind(df_max.SingleR.cycling, df_max.SingleR.EMTAB12889, df_max.scSorter, df_max.scType)
annotation.results <- annotation.results[,c(1,2,3,5,6,8,9,11,12)]
colnames(annotation.results) <- c("Cluster", "SingleR_cycling_Celltype", "SingleR_cycling_Cellcount", "SingleR_EMTAB12889_Celltype", "SingleR_EMTAB12889_Cellcount", "scSorter_Celltype", "scSorter_Cellcount", "scType_Celltype", "scType_Cellcount")

save(annotation.results, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_negative_cells_annotation_results.RData",sep = "_"))
write.table(annotation.results, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_negative_cells_annotation_results.txt",sep = "_"), quote = FALSE)

################################
# Manual marker gene analysis
################################

ovary.Ptprc.neg <- SetIdent(ovary.Ptprc.neg, value = "seurat_clusters")

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_negative_cells_Vlnplot_marker_genes_clusters_0_to_7.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.neg, subset = seurat_clusters %in% c(0:7)), 
        features = c("Ddx4",                             # Oocyte marker
                     "Amh", "Foxl2", "Cyp19a1",          # Granulosa cell marker
                     "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                     "Pdgfra", "Col1a1",                 # Stroma cell marker
                     "Rgs5", "Acta2", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                     "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                     "Pax8", "Cdh1"))                    # Epithelial cell marker
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_negative_cells_Vlnplot_marker_genes_clusters_8_to_16.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.neg, subset = seurat_clusters %in% c(8:16)), 
        features = c("Ddx4",                             # Oocyte marker
                     "Amh", "Foxl2", "Cyp19a1",          # Granulosa cell marker
                     "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                     "Pdgfra", "Col1a1",                 # Stroma cell marker
                     "Rgs5", "Acta2", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                     "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                     "Pax8", "Cdh1"))                    # Epithelial cell marker
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_negative_cells_Vlnplot_marker_genes_clusters_17_to24.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.neg, subset = seurat_clusters %in% c(17:24)), 
        features = c("Ddx4",                             # Oocyte marker
                     "Amh", "Foxl2", "Cyp19a1",          # Granulosa cell marker
                     "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                     "Pdgfra", "Col1a1",                 # Stroma cell marker
                     "Rgs5", "Acta2", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                     "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                     "Pax8", "Cdh1"))                    # Epithelial cell marker
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_negative_cells_Vlnplot_marker_genes_clusters_25_to_32.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.neg, subset = seurat_clusters %in% c(25:32)), 
        features = c("Ddx4",                             # Oocyte marker
                     "Amh", "Foxl2", "Cyp19a1",          # Granulosa cell marker
                     "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                     "Pdgfra", "Col1a1",                 # Stroma cell marker
                     "Rgs5", "Acta2", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                     "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                     "Pax8", "Cdh1"))                    # Epithelial cell marker
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_negative_cells_Vlnplot_marker_genes_clusters_33_to_40.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.neg, subset = seurat_clusters %in% c(33:40)), 
        features = c("Ddx4",                             # Oocyte marker
                     "Amh", "Foxl2", "Cyp19a1",          # Granulosa cell marker
                     "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                     "Pdgfra", "Col1a1",                 # Stroma cell marker
                     "Rgs5", "Acta2", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                     "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                     "Pax8", "Cdh1"))                    # Epithelial cell marker
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_mouse_EMTAB12889Ptprc_negative_cells_Vlnplot_marker_genes_clusters_41_to_49.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.Ptprc.neg, subset = seurat_clusters %in% c(41:49)), 
        features = c("Ddx4",                             # Oocyte marker
                     "Amh", "Foxl2", "Cyp19a1",          # Granulosa cell marker
                     "Cyp17a1", "Cyp11a1", "Lhcgr",      # Theca cell marker 
                     "Pdgfra", "Col1a1",                 # Stroma cell marker
                     "Rgs5", "Acta2", "Notch3",          # Mesenchymal cell marker (pericyte + smooth muscle cell)
                     "Vwf", "Flt1", "Flt4", "Prox1",     # Endothelial cell marker (BEC + LEC)
                     "Pax8", "Cdh1"))                    # Epithelial cell marker
dev.off()

################################
# Transfer final annotations
# Refer to ./2024-10-23_10x_ovary_Mouse_EMTAB12889_Ptprc_negative_cells_annotation_results_FINAL.txt
################################

ovary.Ptprc.neg@meta.data$celltype.level2 <- rep("NA",dim(ovary.Ptprc.neg@meta.data)[1])
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 0] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 1] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 2] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 3] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 4] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 5] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 6] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 7] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 8] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 9] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 10] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 11] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 12] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 13] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 14] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 15] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 16] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 17] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 18] <- "SMC"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 19] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 20] <- "BEC"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 21] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 22] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 23] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 24] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 25] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 26] <- "Pericyte"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 27] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 28] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 29] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 30] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 31] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 32] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 33] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 34] <- "SMC"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 35] <- "Pericyte"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 36] <- "Epithelial"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 37] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 38] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 39] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 40] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 41] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 42] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 43] <- "Pericyte"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 44] <- "SMC"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 45] <- "Stroma"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 46] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 47] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 48] <- "Theca"
ovary.Ptprc.neg@meta.data$celltype.level2[ovary.Ptprc.neg@meta.data$seurat_clusters %in% 49] <- "Stroma"

save(ovary.Ptprc.neg, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Ptprc_negative_cells_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 4.Combine cell type annotation
################################################################################

# Extract celltype.level2 metadata
celltype.level2.pos <- ovary.Ptprc.pos@meta.data[, "celltype.level2", drop = FALSE]
celltype.level2.neg <- ovary.Ptprc.neg@meta.data[, "celltype.level2", drop = FALSE]

# Combine the metadata from both subsetted objects
combined.celltype.level2 <- rbind(celltype.level2.pos, celltype.level2.neg)

# Ensure that the cell IDs match
common_cell_ids <- intersect(rownames(ovary.EMTAB12889@meta.data), rownames(combined.celltype.level2))

# Update the celltype.level2 in ovary.EMTAB12889 for the matching cells
ovary.EMTAB12889@meta.data[common_cell_ids, "celltype.level2"] <- combined.celltype.level2[common_cell_ids, "celltype.level2"]

# Verify that the metadata has been successfully updated
table(is.na(ovary.EMTAB12889@meta.data$celltype.level2))
# FALSE 
# 21301

save(ovary.EMTAB12889, file = paste(Sys.Date(),"10x_ovary_Mouse_EMTAB12889_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Mouse_EMTAB12889_cell_type_annotation_session_Info.txt", sep =""))
sessionInfo()
sink()
