options(stringsAsFactors = FALSE)

library(Seurat)
library(Matrix)

################################################################################
# 10x ovary Mouse GSE290742
# Process GSE290742 dataset - 2M vs. 9M
# Raw data only accessible as counts matrix
################################################################################

# Function to keeps gene names in column 1
# Use only numeric count columns before sparse conversion.
df_counts_to_dgC <- function(df, gene_col = 1L) {
  genes <- df[[gene_col]]
  mat <- as.matrix(df[, -gene_col, drop = FALSE])
  storage.mode(mat) <- "numeric"
  out <- as(mat, "dgCMatrix")
  rownames(out) <- genes
  out
}

################################################################################
# 1. Load counts matrix + metadata
################################################################################

# Counts

Young1.counts <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821060_Young1_10X_counts.csv.gz"), as.is = TRUE, check.names = FALSE)
Young2.counts <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821062_Young2_10X_counts.csv.gz"), as.is = TRUE, check.names = FALSE)
Young3.counts <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821064_Young3_10X_counts.csv.gz"), as.is = TRUE, check.names = FALSE)

Aged1.counts <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821066_Aged1_10X_counts.csv.gz"), as.is = TRUE, check.names = FALSE)
Aged2.counts <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821068_Aged2_10X_counts.csv.gz"), as.is = TRUE, check.names = FALSE)
Aged3.counts <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821070_Aged3_10X_counts.csv.gz"), as.is = TRUE, check.names = FALSE)

Young1.counts <- df_counts_to_dgC(Young1.counts)
Young2.counts <- df_counts_to_dgC(Young2.counts)
Young3.counts <- df_counts_to_dgC(Young3.counts)
Aged1.counts  <- df_counts_to_dgC(Aged1.counts)
Aged2.counts  <- df_counts_to_dgC(Aged2.counts)
Aged3.counts  <- df_counts_to_dgC(Aged3.counts)

# Metadata

Young1.metadata <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821060_Young1_10X_Metadata.csv.gz"), as.is = TRUE)
Young2.metadata <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821062_Young2_10X_Metadata.csv.gz"), as.is = TRUE)
Young3.metadata <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821064_Young3_10X_Metadata.csv.gz"), as.is = TRUE)

Aged1.metadata <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821066_Aged1_10X_Metadata.csv.gz"), as.is = TRUE)
Aged2.metadata <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821068_Aged2_10X_Metadata.csv.gz"), as.is = TRUE)
Aged3.metadata <- read.csv(gzfile("/Volumes/OIProject_II/1_R/1_Pre-processing/Mouse/Mouse_GSE290742/1_Preprocessing/GSE290742_RAW/GSM8821070_Aged3_10X_Metadata.csv.gz"), as.is = TRUE)

rownames(Young1.metadata) <- Young1.metadata$X
rownames(Young2.metadata) <- Young2.metadata$X
rownames(Young3.metadata) <- Young3.metadata$X

rownames(Aged1.metadata) <- Aged1.metadata$X
rownames(Aged2.metadata) <- Aged2.metadata$X
rownames(Aged3.metadata) <- Aged3.metadata$X

################################################################################
# 2. Create Seurat objects
################################################################################

Young1 <- CreateSeuratObject(counts = Young1.counts, meta.data = Young1.metadata)
Young2 <- CreateSeuratObject(counts = Young2.counts, meta.data = Young2.metadata)
Young3 <- CreateSeuratObject(counts = Young3.counts, meta.data = Young3.metadata)

Aged1 <- CreateSeuratObject(counts = Aged1.counts, meta.data = Aged1.metadata)
Aged2 <- CreateSeuratObject(counts = Aged2.counts, meta.data = Aged2.metadata)
Aged3 <- CreateSeuratObject(counts = Aged3.counts, meta.data = Aged3.metadata)

samples <- c("Young1", "Young2", "Young3",
             "Aged1", "Aged2", "Aged3")

# Keep only QC-qualified cells
Young1.cl <- subset(Young1, subset = QC_Filter == "Kept")
Young2.cl <- subset(Young2, subset = QC_Filter == "Kept")
Young3.cl <- subset(Young3, subset = QC_Filter == "Kept")
Aged1.cl  <- subset(Aged1, subset = QC_Filter == "Kept")
Aged2.cl  <- subset(Aged2, subset = QC_Filter == "Kept")
Aged3.cl  <- subset(Aged3, subset = QC_Filter == "Kept")

# Merge objects
ovary.GSE290742 <- merge(Young1.cl, y = c(Young2.cl, Young3.cl,
                                          Aged1.cl, Aged2.cl, Aged3.cl),
                         project = "ovary_GSE290742")

save(ovary.GSE290742, file = paste0(Sys.Date(), "_CZI_Mouse_GSE290742_seurat_object_raw_counts.RData"))

################################################################################
sink(file = paste0(Sys.Date(), "_Mouse_GSE290742_preprocessing_session_info.txt"))
sessionInfo()
sink()
