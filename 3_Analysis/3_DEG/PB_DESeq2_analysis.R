library('Seurat')
library('muscat')
library('sctransform')
library('DESeq2')
library(sva)
library(limma)
library(dplyr)
library(tibble)
library(SingleCellExperiment)
library('biomaRt')
library(tidyverse)
library('ggplot2')
library(RColorBrewer)

rm(list = ls())

################################################################################
# Integrative ovarian aging analysis
# Perform PB DESeq2 analysis for each dataset
# Humanize DESeq2 objects for non-human datasets
################################################################################

################################################################################
# SECTION 1: Perform PB -> DESeq2
################################################################################

################################################################################
# 1.1.1 Human GSE255690 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-06-26_10x_ovary_Human_GSE255690_celltype_annotated_Seurat_object_combined_final.RData")

DefaultAssay(ovary.Human.GSE255690.cl) <- "RNA"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.Human.GSE255690.cl)
rm(ovary.Human.GSE255690.cl)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age",              
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 6]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
vst.counts <- list()
deseq.res.list <- list()
significant_gene_counts <- list()

age_map_numeric <- function(samples) {
  sapply(samples, function(x) {
    if (grepl("Young", x)) return(23) 
    if (grepl("Middle", x)) return(37.5)
    if (grepl("Old", x)) return(48)
    return(NA)
  })
}

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 6
  counts <- counts[good_genes, ]

  sample_ages <- age_map_numeric(colnames(counts))
  dataDesign <- data.frame(row.names = colnames(counts), age = as.numeric(sample_ages))
  
  # Estimate surrogate variables
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(as.matrix(counts), mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  # Remove unwanted variation
  if (n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = filtered.sva, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  res <- results(dds, name = "age")  
  
  res <- res[!is.na(res$padj), ]
  
  # Store filtered results
  deseq.res.list[[cell_type]] <- res
  
  # Variance-stabilized data
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Human_GSE255690_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Human_GSE255690_PB_DESeq2_object.RData")))

################################################################################
# 1.1.2 Human GSE202601 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-06-26_10x_ovary_Human_GSE202601_celltype_annotated_Seurat_object_combined_final.RData")

DefaultAssay(ovary.Human.GSE202601.cl) <- "RNA"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.Human.GSE202601.cl)
rm(ovary.Human.GSE202601.cl)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age",              
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 6]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list() 

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 6
  counts <- counts[good_genes, ]
  
  # Create DESeq2 dataset
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = as.numeric(ifelse(grepl("Human2", colnames(counts)), 27.5, 51.5))
  )
  
  # Estimate surrogate variables
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(as.matrix(counts), mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  # Remove unwanted variation
  if (n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = filtered.sva, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  res <- results(dds, name = "age")  
  
  res <- res[!is.na(res$padj), ]
  
  # Store filtered results
  deseq.res.list[[cell_type]] <- res
  
  # Variance-stabilized data
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Human_GSE202601_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Human_GSE202601_PB_DESeq2_object.RData")))

################################################################################
# 1.2.1 Mouse Benayoun Lab Aging Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-06-26_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.AC) <- "RNA"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.AC)
rm(ovary.AC)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age",              
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 4]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list()

age_map_numeric <- function(samples) {
  sapply(samples, function(x) {
    if (grepl("YF", x)) return(4)  
    if (grepl("OF", x)) return(20) 
    return(NA)
  })
}

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 3
  counts <- counts[good_genes, ]
  if (nrow(counts) == 0) next
  
  sample_ages <- age_map_numeric(colnames(counts))
  dataDesign <- data.frame(row.names = colnames(counts), age = as.numeric(sample_ages))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  res <- results(dds, name = "age")  
  res <- res[!is.na(res$padj), ]
  
  deseq.res.list[[cell_type]] <- res
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_Aging_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_Aging_PB_DESeq2_object.RData")))

################################################################################
# 1.2.2 Mouse VCD Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.VCD) <- "RNA"

# Define the mapping
age_numeric_map <- c(
  "3m_30d"  = 5,
  "3m_90d"  = 7,
  "10m_30d" = 12,
  "10m_90d" = 14
)

# Extract the Seurat object metadata
meta <- ovary.VCD@meta.data

# Extract the matching key from the Library column
meta$Age_numeric <- sapply(meta$Library, function(x) {
  key <- sub(".*(3m_30d|3m_90d|10m_30d|10m_90d).*", "\\1", x)
  age_numeric_map[[key]]
})

# Add the updated metadata back
ovary.VCD@meta.data <- meta

# Subset only the relevant columns for pseudobulk analysis
sample_metadata <- ovary.VCD@meta.data %>%
  dplyr::select(Library, Age_numeric, Batch) %>%
  distinct()

# Rename columns to match the required format
colnames(sample_metadata) <- c("sample_id", "Age", "Batch")

# Set rownames to match sample IDs
rownames(sample_metadata) <- sample_metadata$sample_id

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.VCD)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age",              
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 6]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
sva.cleaned.counts <- vector(mode = "list", length = length(counts.pb))
names(sva.cleaned.counts) <- names(counts.pb)

vst.counts <- vector("list", length = length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length = length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list() 

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 6
  counts <- counts[good_genes, ]
  if (nrow(counts) == 0) next
  
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = sample_metadata[colnames(counts), "Age"],
    batch = sample_metadata[colnames(counts), "Batch"]
  )
  
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(counts, mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  if (svseq$n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch, covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch, design = mod1)
  }
  
  filtered.sva <- round(2^clean - 0.1)
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = dataDesign, design = ~ age + batch)
  dds <- DESeq(dds)
  res <- results(dds, name = "age")
  res <- res[!is.na(res$padj), ]
  vst_data <- getVarianceStabilizedData(dds)
  
  deseq.res.list[[cell_type]] <- res
  vst.counts[[cell_type]] <- vst_data
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_VCD_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_VCD_PB_DESeq2_object.RData")))

cat("Mouse VCD analysis completed.\n")

################################################################################
# 1.2.3 Mouse Foxl2 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-06-27_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.Foxl2) <- "RNA"

# Define the mapping
age_numeric_map <- c(
  "young"  = 4,
  "old"  = 9
)

# Extract the Seurat object metadata
meta <- ovary.Foxl2@meta.data

# Extract the matching key from the Library column
meta$Age_numeric <- sapply(meta$Library, function(x) {
  key <- sub(".*(young|old).*", "\\1", x)
  age_numeric_map[[key]]
})

# Add the updated metadata back
ovary.Foxl2@meta.data <- meta

# Subset only the relevant columns for pseudobulk analysis
metadata <- ovary.Foxl2@meta.data %>%
  dplyr::select(Library, Age_numeric, Batch) %>%
  distinct()

# Rename columns to match the required format
colnames(metadata) <- c("sample_id", "Age", "Batch")

# Set rownames to match sample IDs
rownames(metadata) <- metadata$sample_id

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.Foxl2)
rm(ovary.Foxl2)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age_numeric",      
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) > 3]

# Drop B cells: Foxl2_wt_old_2 has 0 cells. Analyze independently
celltype.qc <- celltype.qc[celltype.qc != "B"]

counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
sva.cleaned.counts <- vector(mode = "list", length = length(counts.pb))
names(sva.cleaned.counts) <- names(counts.pb)

vst.counts <- vector("list", length = length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length = length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list() 

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 3
  counts <- counts[good_genes, ]
  
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = metadata[colnames(counts), "Age"],
    batch = metadata[colnames(counts), "Batch"]
  )
  
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(counts, mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  if (svseq$n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch, covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch, design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = dataDesign, design = ~ age + batch)
  dds <- DESeq(dds)
  res <- results(dds, name = "age")
  res <- res[!is.na(res$padj), ]
  vst_data <- getVarianceStabilizedData(dds)
  
  deseq.res.list[[cell_type]] <- res
  vst.counts[[cell_type]] <- vst_data
}

# Assess B cells separately
cell_type <- "B"
counts <- pb@assays@data["B"]
counts <- counts$B
counts <- counts[, colnames(counts) != "Foxl2_wt_old_2"]

good_genes <- rowSums(counts > 0) >= 3
counts <- counts[good_genes, ]

dataDesign <- data.frame(
  row.names = colnames(counts),
  age = metadata[colnames(counts), "Age"],
  batch = metadata[colnames(counts), "Batch"]
)

mod1 <- model.matrix(~ age, data = dataDesign)
n.sv <- num.sv(counts, mod1, method = "be")
svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)

if (svseq$n.sv > 0) {
  clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch, covariates = svseq$sv, design = mod1)
} else {
  clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch, design = mod1)
}

filtered.sva <- round(2^clean - 0.1)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = dataDesign, design = ~ age + batch)
dds <- DESeq(dds)
res <- results(dds, name = "age")
res <- res[!is.na(res$padj), ]
vst_data <- getVarianceStabilizedData(dds)

deseq.res.list[[cell_type]] <- res
vst.counts[[cell_type]] <- vst_data

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_Foxl2_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_Foxl2_PB_DESeq2_object.RData")))

################################################################################
# 1.2.4 Mouse EMTAB12889 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-07-23_10x_ovary_Mouse_EMTAB12889_Seurat_object_with_final_annotation.RData")

# Fix naming for consistency
ovary.EMTAB12889@meta.data$celltype.level2[ovary.EMTAB12889@meta.data$celltype.level2 == "CD8 NKT"] <- "CD8NKT"

DefaultAssay(ovary.EMTAB12889) <- "RNA"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.EMTAB12889)
rm(ovary.EMTAB12889)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age",              
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
pb_mds <- pbMDS(pb)

pdf(paste0(Sys.Date(),"_Mouse_EMTAB12889_Muscat_PB_MDS.pdf"))
pb_mds
dev.off()

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 9]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
age_map <- function(samples) {
  as.numeric(sub(".*_(\\d+)M_.*", "\\1", samples))
}

vst.counts <- list()
deseq.res.list <- list()
significant_gene_counts <- list()

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 9
  counts <- counts[good_genes, ]
  if (nrow(counts) == 0) next
  
  sample_ages <- age_map(colnames(counts))
  dataDesign <- data.frame(row.names = colnames(counts), age = as.numeric(sample_ages))
  
  # Estimate surrogate variables
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(as.matrix(counts), mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  # Remove unwanted variation
  if (n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  dds <- DESeqDataSetFromMatrix(countData = filtered.sva, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
  
  res_lm <- results(dds, name = "age")
  deseq.res.list[[cell_type]] <- list(linear_age = res_lm)
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_EMTAB12889_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_EMTAB12889_PB_DESeq2_object.RData")))

################################################################################
# 1.2.5 Mouse GSE232309 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-07-23_10x_ovary_Mouse_GSE232309_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.GSE232309) <- "RNA"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.GSE232309)
rm(ovary.GSE232309)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "Age",              
               sid = "Library",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 6]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list()

age_map_numeric <- function(samples) {
  sapply(samples, function(x) {
    if (grepl("3M", x)) return(3)
    if (grepl("9M", x)) return(9)
    return(NA)
  })
}

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 6
  counts <- counts[good_genes, ]
  if (nrow(counts) == 0) next
  
  sample_ages <- age_map_numeric(colnames(counts))
  dataDesign <- data.frame(row.names = colnames(counts), age = as.numeric(sample_ages))
  
  # Estimate surrogate variables
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(as.matrix(counts), mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  # Remove unwanted variation
  if (n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  dds <- DESeqDataSetFromMatrix(countData = filtered.sva, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  
  res <- results(dds, name = "age") 
  res <- res[!is.na(res$padj), ]
  
  deseq.res.list[[cell_type]] <- res
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_GSE232309_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Mouse_GSE232309_PB_DESeq2_object.RData")))

################################################################################
# 1.3.1 Monkey STRTseq GSE130664 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-07-02_STRTseq_ovary_Monkey_GSE130664_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.Monkey.GSE130664.cl) <- "RNA"

ovary.Monkey.GSE130664.cl@meta.data$celltype.level2[ovary.Monkey.GSE130664.cl@meta.data$celltype.level2 == "Mesenchyme"] <- "SMC"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  heatmaps = file.path("./Heatmaps"),
  mds_plots = file.path("./MDS_Plots"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and aggregate counts
sce <- as.SingleCellExperiment(ovary.Monkey.GSE130664.cl)
rm(ovary.Monkey.GSE130664.cl)

sce <- prepSCE(sce, 
               kid = "celltype.level2",  
               gid = "age_group",              
               sid = "individual",          
               drop = TRUE)              

nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 15, 2, sum) >= 6]
counts.pb <- pb@assays@data[celltype.qc]

# Perform DESeq2 analysis for each cell type
vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list()

for (cell_type in names(counts.pb)) {
  counts <- counts.pb[[cell_type]]
  good_genes <- rowSums(counts > 0) >= 6
  counts <- counts[good_genes, ]
  
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = as.numeric(ifelse(grepl("Y", colnames(counts)), 5, 19))
  )
  
  # Estimate surrogate variables
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(as.matrix(counts), mod1, method = "be")
  svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv)
  
  # Remove unwanted variation
  if (n.sv > 0) {
    clean <- removeBatchEffect(log2(counts + 0.1), covariates = svseq$sv, design = mod1)
  } else {
    clean <- removeBatchEffect(log2(counts + 0.1), design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  # Run DESeq2 with cleaned data
  dds <- DESeqDataSetFromMatrix(countData = filtered.sva, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  res <- results(dds, name = "age")
  
  res <- res[!is.na(res$padj), ]
  
  deseq.res.list[[cell_type]] <- res
  
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Monkey_GSE130664_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Monkey_GSE130664_PB_DESeq2_object.RData")))

################################################################################
# 1.4.1 Goat PRJNA1010653 Analysis
################################################################################

# Load data
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/0_Annotated_seurat_objects/2025-07-01_10x_ovary_Goat_PRJNA1010653_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.Goat) <- "RNA"

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Convert to SingleCellExperiment and generate pseudoreplicates
sce <- as.SingleCellExperiment(ovary.Goat)
rm(ovary.Goat)

sce <- prepSCE(sce,  
               kid = "celltype.level2",  
               gid = "Group",              
               sid = "Library",          
               drop = TRUE)

# Function to generate pseudobulk replicates
generate_pseudobulk_replicates <- function(sce, min_cells = 20) {
  set.seed(123)
  count_list <- list()
  metadata_list <- list()
  celltypes <- unique(sce$cluster_id)
  groups <- unique(sce$group_id)
  for (ct in celltypes) {
    for (grp in groups) {
      sce_sub <- sce[, sce$cluster_id == ct & sce$group_id == grp]
      n_cells <- ncol(sce_sub)
      if (n_cells < min_cells) next
      idx <- sample(seq_len(n_cells))
      half <- floor(n_cells / 2)
      rep1_cells <- idx[1:half]
      rep2_cells <- idx[(half+1):n_cells]
      rep1_counts <- rowSums(assay(sce_sub)[, rep1_cells, drop=FALSE])
      rep2_counts <- rowSums(assay(sce_sub)[, rep2_cells, drop=FALSE])
      sample1_name <- paste0(grp, "_", ct, "_rep1")
      sample2_name <- paste0(grp, "_", ct, "_rep2")
      count_list[[sample1_name]] <- rep1_counts
      count_list[[sample2_name]] <- rep2_counts
      metadata_list[[sample1_name]] <- data.frame(sample = sample1_name, group = grp, celltype = ct, replicate = "rep1")
      metadata_list[[sample2_name]] <- data.frame(sample = sample2_name, group = grp, celltype = ct, replicate = "rep2")
    }
  }
  pseudo_counts <- do.call(cbind, count_list)
  pseudo_coldata <- do.call(rbind, metadata_list)
  rownames(pseudo_coldata) <- pseudo_coldata$sample
  list(counts = pseudo_counts, coldata = pseudo_coldata)
}

pseudo_out <- generate_pseudobulk_replicates(sce, min_cells = 15)
pseudo_counts <- pseudo_out$counts
pseudo_coldata <- pseudo_out$coldata

# Perform DESeq2 analysis for each cell type
celltype.qc <- unique(pseudo_coldata$celltype)

vst.counts <- list()
deseq.res.list <- list()
significant_gene_counts <- list()

for (ct in celltype.qc) {
  ct_samples <- pseudo_coldata[pseudo_coldata$celltype == ct, ]
  counts_ct <- pseudo_counts[, ct_samples$sample]
  good_genes <- rowSums(counts_ct > 0) >= 3
  counts_ct <- counts_ct[good_genes, ]
  if (nrow(counts_ct) == 0) next
  
  dataDesign <- data.frame(
    row.names = rownames(ct_samples),
    age = as.numeric(ifelse(grepl("young", rownames(ct_samples)), 2, 10))
  )
  
  mod1 <- model.matrix(~ age, data = dataDesign)
  n.sv <- num.sv(as.matrix(counts_ct), mod1, method = "be")
  svseq <- svaseq(as.matrix(counts_ct), mod1, n.sv = n.sv)
  clean <- if (n.sv > 0) {
    removeBatchEffect(log2(counts_ct + 0.1), covariates = svseq$sv, design = mod1)
  } else {
    removeBatchEffect(log2(counts_ct + 0.1), design = mod1)
  }
  filtered.sva <- round(2^clean - 0.1)
  
  dds <- DESeqDataSetFromMatrix(countData = filtered.sva, colData = dataDesign, design = ~ age)
  dds <- DESeq(dds)
  res <- results(dds)
  
  res <- res[!is.na(res$padj), ]
  
  deseq.res.list[[ct]] <- res
  
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[ct]] <- vst_data
  
  sig_genes <- rownames(res)[res$padj < 0.05]
  significant_gene_counts[[ct]] <- length(sig_genes)
}

# Save results
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Goat_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Goat_PB_DESeq2_object.RData")))

################################################################################
# SECTION 2: Humanize DESeq2 objects
################################################################################

################################################################################
# 2.1 Humanize DESeq2 Objects
################################################################################

# Import data from all species
# Monkey
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Monkey/2025-08-06/DESeq2_results/2025-08-06_STRTseq_ovary_Monkey_GSE130664_DESeq2_object.RData")
Monkey.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

# Mouse datasets
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Mouse/Mouse_Benayoun_lab_AC/2025-07-23/DESeq2_results/2025-07-23_10x_ovary_Benayoun_lab_AC_PB_DESeq2_object.RData")
Mouse.Aging.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Mouse/Mouse_Benayoun_lab_VCD/2025-07-23/DESeq2_results/2025-07-23_10x_ovary_Benayoun_lab_VCD_controls_only_PB_DESeq2_object.RData")
Mouse.VCD.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Mouse/Mouse_Benayoun_lab_Foxl2/2025-07-23/DESeq2_results/2025-07-23_10x_ovary_Benayoun_lab_Foxl2_wt_only_PB_DESeq2_object.RData")
Mouse.Foxl2.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Mouse/Mouse_EMTAB12889/2025-07-23/DESeq2_results/2025-07-23_PB_DESeq2_object.RData")
Mouse.EMTAB12889.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

# Fix EMTAB12889 structure - extract linear_age results
Mouse.EMTAB12889.deseq.res.list <- lapply(Mouse.EMTAB12889.deseq.res.list, function(x) {
  if (is.list(x) && "linear_age" %in% names(x)) {
    return(x$linear_age)
  } else {
    return(x)
  }
})

load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Mouse/Mouse_GSE232309/2025-07-23/DESeq2_results/2025-07-23_10x_ovary_Mouse_3M_vs_9M_DESeq2_object.RData")
Mouse.GSE232309.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

# Goat
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Goat/2025-08-06/DESeq2_results/2025-08-06_Goat_PB_DESeq2_object.RData")
Goat.deseq.res.list <- deseq.res.list
rm(deseq.res.list)

# Function to humanize DESeq2 results
humanize_deseq_results <- function(deseq_list, ortholog_df, species_name) {
  lapply(names(deseq_list), function(cell_type) {
    res <- deseq_list[[cell_type]]
    
    # Convert DESeqResults object to df
    if (inherits(res, "DESeqResults")) {
      res <- as.data.frame(res)
    }
        
    # Get gene names and check for orthologs
    genes <- rownames(res)
    
    valid_genes <- genes %in% rownames(ortholog_df)
    
    # Filter and rename
    res_filtered <- res[valid_genes, , drop = FALSE]
    rownames(res_filtered) <- ortholog_df[rownames(res_filtered), "hsapiens_homolog_associated_gene_name"]
    
    return(res_filtered)
  })
}

################################################################################
# 2.2 Humanize Monkey genes
################################################################################

# Prep DB
ensembl_monkey <- useEnsembl(biomart = "genes", 
                             dataset = "mfascicularis_gene_ensembl", 
                             mirror = "useast")

# Get unique monkey gene names from your DESeq2 results
monkey.genes <- unique(unlist(lapply(Monkey.deseq.res.list, rownames)))

# Get ortholog mapping
df_genes_monkey <- getBM(
  mart = ensembl_monkey,
  filters = c("external_gene_name"),
  values = monkey.genes,
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_associated_gene_name", 
                 "hsapiens_homolog_orthology_type")
)

# Keep only one-to-one orthologs
df_genes_monkey <- df_genes_monkey %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one") %>%
  filter(hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

rownames(df_genes_monkey) <- df_genes_monkey$external_gene_name

# Convert each DESeq2 result table
hMonkey.deseq.res.list <- humanize_deseq_results(Monkey.deseq.res.list, df_genes_monkey, "Monkey")
names(hMonkey.deseq.res.list) <- names(Monkey.deseq.res.list)

# Save data
save(hMonkey.deseq.res.list, file = paste0(Sys.Date(), "_Monkey_GSE130664_humanized_DESeq2_results_object.RData"))

################################################################################
# 3.3 Humanize Mouse genes
################################################################################

# Prep DB for mouse
ensembl_mouse <- useEnsembl(biomart = "genes", 
                            dataset = "mmusculus_gene_ensembl", 
                            mirror = "useast")

# Get unique mouse gene names from all mouse DESeq2 results
mouse.genes <- unique(c(
  unlist(lapply(Mouse.Aging.deseq.res.list, rownames)),
  unlist(lapply(Mouse.VCD.deseq.res.list, rownames)),
  unlist(lapply(Mouse.Foxl2.deseq.res.list, rownames)),
  unlist(lapply(Mouse.EMTAB12889.deseq.res.list, rownames)),
  unlist(lapply(Mouse.GSE232309.deseq.res.list, rownames))
))

# Get ortholog mapping for mouse
df_genes_mouse <- getBM(
  mart = ensembl_mouse,
  filters = c("external_gene_name"),
  values = mouse.genes,
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_associated_gene_name", 
                 "hsapiens_homolog_orthology_type")
)

# Keep only one-to-one orthologs
df_genes_mouse <- df_genes_mouse %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one") %>%
  filter(hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

rownames(df_genes_mouse) <- df_genes_mouse$external_gene_name

# Convert each mouse DESeq2 result table
hMouse.Aging.deseq.res.list <- humanize_deseq_results(Mouse.Aging.deseq.res.list, df_genes_mouse, "Mouse Aging")
names(hMouse.Aging.deseq.res.list) <- names(Mouse.Aging.deseq.res.list)

hMouse.VCD.deseq.res.list <- humanize_deseq_results(Mouse.VCD.deseq.res.list, df_genes_mouse, "Mouse VCD")
names(hMouse.VCD.deseq.res.list) <- names(Mouse.VCD.deseq.res.list)

hMouse.Foxl2.deseq.res.list <- humanize_deseq_results(Mouse.Foxl2.deseq.res.list, df_genes_mouse, "Mouse Foxl2")
names(hMouse.Foxl2.deseq.res.list) <- names(Mouse.Foxl2.deseq.res.list)

hMouse.EMTAB12889.deseq.res.list <- humanize_deseq_results(Mouse.EMTAB12889.deseq.res.list, df_genes_mouse, "Mouse EMTAB12889")
names(hMouse.EMTAB12889.deseq.res.list) <- names(Mouse.EMTAB12889.deseq.res.list)

hMouse.GSE232309.deseq.res.list <- humanize_deseq_results(Mouse.GSE232309.deseq.res.list, df_genes_mouse, "Mouse GSE232309")
names(hMouse.GSE232309.deseq.res.list) <- names(Mouse.GSE232309.deseq.res.list)

# Save humanized mouse data
save(hMouse.Aging.deseq.res.list, file = paste0(Sys.Date(), "_Mouse_Aging_humanized_DESeq2_results_object.RData"))
save(hMouse.VCD.deseq.res.list, file = paste0(Sys.Date(), "_Mouse_VCD_humanized_DESeq2_results_object.RData"))
save(hMouse.Foxl2.deseq.res.list, file = paste0(Sys.Date(), "_Mouse_Foxl2_humanized_DESeq2_results_object.RData"))
save(hMouse.EMTAB12889.deseq.res.list, file = paste0(Sys.Date(), "_Mouse_EMTAB12889_humanized_DESeq2_results_object.RData"))
save(hMouse.GSE232309.deseq.res.list, file = paste0(Sys.Date(), "_Mouse_GSE232309_humanized_DESeq2_results_object.RData"))

################################################################################
# 2.4 Humanize Goat genes
################################################################################

# Prep DB for goat
ensembl_goat <- useEnsembl(biomart = "genes", 
                           dataset = "chircus_gene_ensembl", 
                           mirror = "useast")

# Get unique goat gene names from your DESeq2 results
goat.genes <- unique(unlist(lapply(Goat.deseq.res.list, rownames)))

# Get ortholog mapping for goat
df_genes_goat <- getBM(
  mart = ensembl_goat,
  filters = c("external_gene_name"),
  values = goat.genes,
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_associated_gene_name", 
                 "hsapiens_homolog_orthology_type")
)

# Keep only one-to-one orthologs
df_genes_goat <- df_genes_goat %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one") %>%
  filter(hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

rownames(df_genes_goat) <- df_genes_goat$external_gene_name

# Convert each goat DESeq2 result table
hGoat.deseq.res.list <- humanize_deseq_results(Goat.deseq.res.list, df_genes_goat, "Goat")
names(hGoat.deseq.res.list) <- names(Goat.deseq.res.list)

# Save humanized goat data
save(hGoat.deseq.res.list, file = paste0(Sys.Date(), "_Goat_humanized_DESeq2_results_object.RData"))

cat("Humanization of DESeq2 objects completed.\n")

################################################################################
# Save session info
sink(file = paste0(Sys.Date(), "_PB_DESeq2_analysis_session_info.txt"))
sessionInfo()
sink()
