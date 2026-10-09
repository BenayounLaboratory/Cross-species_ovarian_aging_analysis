options(stringsAsFactors = F)

library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library('scProportionTest')

################################################################################
# Integrative ovarian aging analysis
# Plot cell type proportions data
################################################################################

###############################################################################
# Define objects & parameters
###############################################################################

my.freq.level1 <- list()
my.freq.level2 <- list()

# Celltype order
custom_order_level1 <- c("nonimmune", "immune")

custom_order_level2 <- c(
  "Granulosa","Theca","Stroma","SMC","BEC","LEC",
  "Epithelial","Myeloid","DC","ILC","NK","NKT","CD8NKT", 
  "CD8T","CD4T","DNT","DPT","B"
)

# Define functions
compute_freq_table <- function(seurat_obj,
                               level_col,
                               lib_col = "Library",
                               desired_levels,
                               dataset_name = "dataset") {
  
  md <- seurat_obj@meta.data
  
  level_vec <- factor(md[[level_col]], levels = desired_levels)
  lib_vec   <- factor(md[[lib_col]], levels = sort(unique(md[[lib_col]])))
  
  # Counts table
  counts_tbl <- table(level_vec, lib_vec)
  
  # Proportions
  freq_tbl <- prop.table(counts_tbl, margin = 2)
  
  # Convert to matrices 
  counts_mat <- as.matrix(counts_tbl)
  freq_mat   <- as.matrix(freq_tbl)
  
  # For level2, mark entire rows as NA 
  absent_rows <- rowSums(counts_mat, na.rm = TRUE) == 0
  if (any(absent_rows)) {
    freq_mat[absent_rows, ] <- NA_real_
  }
  
  return(list(counts = counts_mat, freq = freq_mat))
}

###############################################################################
# Process human data
###############################################################################

#################################
# Load data
#################################


load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Human/Human_GSE255690/2025-06-23/2025-06-26_10x_ovary_Human_GSE255690_celltype_annotated_Seurat_object_combined_final.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Human/Human_GSE202601/2025-06-26/2025-06-26_10x_ovary_Human_GSE202601_celltype_annotated_Seurat_object_combined_final.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Human/Human_TS/2026-05-04_10x_ovary_Human_TS_celltype_annotated_Seurat_object_combined_final.RData")

#################################
# Generate frequency tables
#################################


# Level 1
res_L1_GSE255690 <- compute_freq_table(
  seurat_obj     = ovary.Human.GSE255690.cl,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Human_GSE255690"
)
my.freq.level1[["Human_GSE255690"]] <- res_L1_GSE255690$freq

res_L1_GSE202601 <- compute_freq_table(
  seurat_obj     = ovary.Human.GSE202601.cl,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Human_GSE202601"
)
my.freq.level1[["Human_GSE202601"]] <- res_L1_GSE202601$freq

res_L1_TS <- compute_freq_table(
  seurat_obj     = ovary.Human.TS.cl,
  level_col      = "celltype.level1",
  lib_col        = "development_stage",
  desired_levels = custom_order_level1,
  dataset_name   = "Human_TS"
)
my.freq.level1[["Human_TS"]] <- res_L1_TS$freq

# Level 2
res_L2_GSE255690 <- compute_freq_table(
  seurat_obj     = ovary.Human.GSE255690.cl,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "GSE255690"
)
my.freq.level2[["GSE255690"]] <- res_L2_GSE255690$freq

res_L2_GSE202601 <- compute_freq_table(
  seurat_obj     = ovary.Human.GSE202601.cl,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "GSE202601"
)
my.freq.level2[["GSE202601"]] <- res_L2_GSE202601$freq

res_L2_TS <- compute_freq_table(
  seurat_obj     = ovary.Human.TS.cl,
  level_col      = "celltype.level2",
  lib_col        = "development_stage",
  desired_levels = custom_order_level2,
  dataset_name   = "Human_TS"
)
my.freq.level2[["Human_TS"]] <- res_L2_TS$freq

#################################
# Assess changes in proportions using scProportionTest
#################################


# Re-annotate age groups

unique(ovary.Human.GSE255690.cl@meta.data$Age)          # [1] "Young"  "Middle" "Old"   
unique(ovary.Human.GSE202601.cl@meta.data$Age)          # [1] "23" "27" "28" "29" "49" "51" "52" "54"
unique(ovary.Human.TS.cl@meta.data$development_stage)   # [1] "55-year-old stage" "26-year-old stage" "56-year-old stage"

ovary.Human.GSE255690.cl@meta.data$Age_group <- "Peripostmenopause"
ovary.Human.GSE255690.cl@meta.data$Age_group[ovary.Human.GSE255690.cl@meta.data$Age == "Young"] <- "Premenopause"

ovary.Human.GSE202601.cl@meta.data$Age_group <- "Peripostmenopause"
ovary.Human.GSE202601.cl@meta.data$Age_group[ovary.Human.GSE202601.cl@meta.data$Age %in% c(23, 27, 28, 29)] <- "Premenopause"

ovary.Human.TS.cl@meta.data$Age_group <- "Peripostmenopause"
ovary.Human.TS.cl@meta.data$Age_group[ovary.Human.TS.cl@meta.data$development_stage == "26-year-old stage"] <- "Premenopause"

# Create prop test objects
ovary.Human.GSE255690.prop_test <- sc_utils(ovary.Human.GSE255690.cl)
ovary.Human.GSE202601.prop_test <- sc_utils(ovary.Human.GSE202601.cl)
ovary.Human.TS.prop_test        <- sc_utils(ovary.Human.TS.cl)

# GSE255690 
ovary.prop_test.Human.GSE255690.level1 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level1",
                                                                sample_1 = "Premenopause",
                                                                sample_2 = "Peripostmenopause",
                                                                sample_identity = "Age_group")

ovary.prop_test.Human.GSE255690.level2 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                           cluster_identity = "celltype.level2",
                                                           sample_1 = "Premenopause",
                                                           sample_2 = "Peripostmenopause",
                                                           sample_identity = "Age_group")

ovary.prop_test.Human.GSE202601.level1 <- permutation_test(ovary.Human.GSE202601.prop_test,
                                                           cluster_identity = "celltype.level1",
                                                           sample_1 = "Premenopause",
                                                           sample_2 = "Peripostmenopause",
                                                           sample_identity = "Age_group")

ovary.prop_test.Human.GSE202601.level2 <- permutation_test(ovary.Human.GSE202601.prop_test,
                                                           cluster_identity = "celltype.level2",
                                                           sample_1 = "Premenopause",
                                                           sample_2 = "Peripostmenopause",
                                                           sample_identity = "Age_group")

ovary.prop_test.Human.TS.level1 <- permutation_test(ovary.Human.TS.prop_test,
                                                            cluster_identity = "celltype.level1",
                                                            sample_1 = "Premenopause",
                                                            sample_2 = "Peripostmenopause",
                                                            sample_identity = "Age_group")

ovary.prop_test.Human.TS.level2 <- permutation_test(ovary.Human.TS.prop_test,
                                                           cluster_identity = "celltype.level2",
                                                           sample_1 = "Premenopause",
                                                           sample_2 = "Peripostmenopause",
                                                           sample_identity = "Age_group")

###############################################################################
# Process mouse data
###############################################################################

#################################
# Load data
#################################

load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-11-05_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2026-04-28_10x_ovary_Foxl2_wt_only_Seurat_object_celltypes_annotated.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_VCD/2025-06-26/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_EMTAB12889/2025-07-22/2025-07-23_10x_ovary_Mouse_EMTAB12889_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_GSE232309/2025-06-30/2025-06-30_10x_ovary_Mouse_GSE232309_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2026-04-22_10x_ovary_Mouse_EMTAB11491_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2026-05-01_10x_ovary_Mouse_GSE290742_Seurat_object_with_final_annotation.RData")

ovary.AC@meta.data$celltype.level2[ovary.AC@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.AC@meta.data$celltype.level2[ovary.AC@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

ovary.Foxl2.wt@meta.data$celltype.level2[ovary.Foxl2.wt@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.Foxl2.wt@meta.data$celltype.level2[ovary.Foxl2.wt@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

ovary.VCD@meta.data$celltype.level2[ovary.VCD@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.VCD@meta.data$celltype.level2[ovary.VCD@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

ovary.EMTAB12889@meta.data$celltype.level2[ovary.EMTAB12889@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.EMTAB12889@meta.data$celltype.level2[ovary.EMTAB12889@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

ovary.EMTAB11491@meta.data$celltype.level2[ovary.EMTAB11491@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.EMTAB11491@meta.data$celltype.level2[ovary.EMTAB11491@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

ovary.GSE232309@meta.data$celltype.level2[ovary.GSE232309@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.GSE232309@meta.data$celltype.level2[ovary.GSE232309@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

ovary.GSE290742@meta.data$celltype.level2[ovary.GSE290742@meta.data$celltype.level2 %in% c("Neutrophil", "Macrophage", "Monocyte")] <- "Myeloid"
ovary.GSE290742@meta.data$celltype.level2[ovary.GSE290742@meta.data$celltype.level2 == "Pericyte"] <- "SMC"

# Age ranges for mouse data
# Benayoun lab aging: 4m vs. 20m
# Benayoun lab Foxl2: 4m vs. 9m
# Benayoun lab VCD: 5m vs 7m vs. 12m vs. 14m
# EMTAB12889: 9m vs. 12m vs. 15m
# GSE232309: 3m vs. 9m

# Group ages together:
# 3-7m: Pre-estropausal
# 9-15m: Peri-estropausal
# 16m~: Post-estropausal

# Rename metadata for consistency
# Benayoun lab aging
ovary.AC$Group <- "PreE"
ovary.AC$Group[ovary.AC$Age == "OF"] <- "PostE"

ovary.AC$celltype.level1[ovary.AC$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.AC$celltype.level1[ovary.AC$celltype.level1 == "Ptprc.pos"] <- "immune"

# Benayoun lab Foxl2
ovary.Foxl2.wt$Group <- "PreE"
ovary.Foxl2.wt$Group[ovary.Foxl2.wt$Age == "Old"] <- "PeriE"
ovary.Foxl2.wt$Group[ovary.Foxl2.wt$Age == "Supold"] <- "PostE"

ovary.Foxl2.wt$celltype.level1[ovary.Foxl2.wt$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.Foxl2.wt$celltype.level1[ovary.Foxl2.wt$celltype.level1 == "Ptprc.pos"] <- "immune"

# Benayoun lab VCD
ovary.VCD$Group <- "PreE"
ovary.VCD$Group[grepl("10m_30d", ovary.VCD$Library)] <- "PeriE"
ovary.VCD$Group[grepl("10m_90d", ovary.VCD$Library)] <- "PeriE"

ovary.VCD$celltype.level1[ovary.VCD$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.VCD$celltype.level1[ovary.VCD$celltype.level1 == "Ptprc.pos"] <- "immune"

# EMTAB12889
ovary.EMTAB12889$Group <- "PeriE"
ovary.EMTAB12889$Group[ovary.EMTAB12889$Age == "15M"] <- "PostE"

ovary.EMTAB12889$celltype.level1[ovary.EMTAB12889$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.EMTAB12889$celltype.level1[ovary.EMTAB12889$celltype.level1 == "Ptprc.pos"] <- "immune"

#EMTAB11491
ovary.EMTAB11491$Group <- "PeriE"
ovary.EMTAB11491$Group[ovary.EMTAB11491$Age == "18M"] <- "PostE"

ovary.EMTAB11491$celltype.level1[ovary.EMTAB11491$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.EMTAB11491$celltype.level1[ovary.EMTAB11491$celltype.level1 == "Ptprc.pos"] <- "immune"

# GSE232309
ovary.GSE232309$Group <- "PreE"
ovary.GSE232309$Group[ovary.GSE232309$Age == "9m"] <- "PeriE"

ovary.GSE232309$celltype.level1[ovary.GSE232309$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.GSE232309$celltype.level1[ovary.GSE232309$celltype.level1 == "Ptprc.pos"] <- "immune"

# GSE290742
ovary.GSE290742$Group <- "PreE"
ovary.GSE290742$Group[ovary.GSE290742$Age == "9mo"] <- "PeriE"

ovary.GSE290742$celltype.level1[ovary.GSE290742$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.GSE290742$celltype.level1[ovary.GSE290742$celltype.level1 == "Ptprc.pos"] <- "immune"

#################################
# Generate frequency tables
#################################

# Level 1
res_L1_Benayoun_Aging <- compute_freq_table(
  seurat_obj     = ovary.AC,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_Benayoun_Aging"
)
my.freq.level1[["Mouse_Benayoun_Aging"]] <- res_L1_Benayoun_Aging$freq

res_L1_Benayoun_Foxl2 <- compute_freq_table(
  seurat_obj     = ovary.Foxl2.wt,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_Benayoun_Foxl2"
)
my.freq.level1[["Mouse_Benayoun_Foxl2"]] <- res_L1_Benayoun_Foxl2$freq

res_L1_Benayoun_VCD <- compute_freq_table(
  seurat_obj     = ovary.VCD,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_Benayoun_VCD"
)
my.freq.level1[["Mouse_Benayoun_VCD"]] <- res_L1_Benayoun_VCD$freq

res_L1_EMTAB12889 <- compute_freq_table(
  seurat_obj     = ovary.EMTAB12889,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_EMTAB12889"
)
my.freq.level1[["Mouse_EMTAB12889"]] <- res_L1_EMTAB12889$freq

res_L1_EMTAB11491 <- compute_freq_table(
  seurat_obj     = ovary.EMTAB11491,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_EMTAB11491"
)
my.freq.level1[["Mouse_EMTAB11491"]] <- res_L1_EMTAB11491$freq

res_L1_GSE232309 <- compute_freq_table(
  seurat_obj     = ovary.GSE232309,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_GSE232309"
)
my.freq.level1[["Mouse_GSE232309"]] <- res_L1_GSE232309$freq

res_L1_GSE290742 <- compute_freq_table(
  seurat_obj     = ovary.GSE290742,
  level_col      = "celltype.level1",
  lib_col        = "SampleID",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_GSE290742"
)
my.freq.level1[["Mouse_GSE290742"]] <- res_L1_GSE290742$freq

# Level 2
res_L2_Benayoun_Aging <- compute_freq_table(
  seurat_obj     = ovary.AC,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_Benayoun_Aging"
)
my.freq.level2[["Mouse_Benayoun_Aging"]] <- res_L2_Benayoun_Aging$freq

res_L2_Benayoun_Foxl2 <- compute_freq_table(
  seurat_obj     = ovary.Foxl2.wt,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_Benayoun_Foxl2"
)
my.freq.level2[["Mouse_Benayoun_Foxl2"]] <- res_L2_Benayoun_Foxl2$freq

res_L2_Benayoun_VCD <- compute_freq_table(
  seurat_obj     = ovary.VCD,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_Benayoun_VCD"
)
my.freq.level2[["Mouse_Benayoun_VCD"]] <- res_L2_Benayoun_VCD$freq

res_L2_EMTAB12889 <- compute_freq_table(
  seurat_obj     = ovary.EMTAB12889,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_EMTAB12889"
)
my.freq.level2[["Mouse_EMTAB12889"]] <- res_L2_EMTAB12889$freq

res_L2_EMTAB11491 <- compute_freq_table(
  seurat_obj     = ovary.EMTAB11491,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_EMTAB11491"
)
my.freq.level2[["Mouse_EMTAB11491"]] <- res_L2_EMTAB11491$freq

res_L2_GSE232309 <- compute_freq_table(
  seurat_obj     = ovary.GSE232309,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_GSE232309"
)
my.freq.level2[["Mouse_GSE232309"]] <- res_L2_GSE232309$freq

res_L2_GSE290742 <- compute_freq_table(
  seurat_obj     = ovary.GSE290742,
  level_col      = "celltype.level2",
  lib_col        = "SampleID",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_GSE290742"
)
my.freq.level2[["Mouse_GSE290742"]] <- res_L2_GSE290742$freq


#################################
# Assess changes in proportions using scProportionTest
#################################

# Create prop test objects
ovary.Mouse.Benayoun.Aging.prop_test <- sc_utils(ovary.AC)
ovary.Mouse.Benayoun.Foxl2.prop_test <- sc_utils(ovary.Foxl2.wt)
ovary.Mouse.Benayoun.VCD.prop_test   <- sc_utils(ovary.VCD)
ovary.Mouse.EMTAB12889.prop_test    <- sc_utils(ovary.EMTAB12889)
ovary.Mouse.EMTAB11491.prop_test    <- sc_utils(ovary.EMTAB11491)
ovary.Mouse.GSE232309.prop_test      <- sc_utils(ovary.GSE232309)
ovary.Mouse.GSE290742.prop_test      <- sc_utils(ovary.GSE290742)

# Benayoun Aging

unique(ovary.AC@meta.data$Group)         # [1] "PreE"      "PostE"

ovary.prop_test.Mouse.Benayoun.Aging.PreEvsPostpostE.level1 <- permutation_test(ovary.Mouse.Benayoun.Aging.prop_test,
                                                                                cluster_identity = "celltype.level1",
                                                                                sample_1 = "PreE",
                                                                                sample_2 = "PostE",
                                                                                sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Aging.PreEvsPostpostE.level2 <- permutation_test(ovary.Mouse.Benayoun.Aging.prop_test,
                                                                                cluster_identity = "celltype.level2",
                                                                                sample_1 = "PreE",
                                                                                sample_2 = "PostE",
                                                                                sample_identity = "Group")
# Benayoun Foxl2

unique(ovary.Foxl2.wt@meta.data$Group)         # [1] "PreE"  "PeriE" "PostE"

ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPeriE.level1 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level1",
                                                                            sample_1 = "PreE",
                                                                            sample_2 = "PeriE",
                                                                            sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPostE.level1 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level1",
                                                                            sample_1 = "PreE",
                                                                            sample_2 = "PostE",
                                                                            sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Foxl2.PeriEvsPostE.level1 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level1",
                                                                            sample_1 = "PeriE",
                                                                            sample_2 = "PostE",
                                                                            sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPeriE.level2 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level2",
                                                                            sample_1 = "PreE",
                                                                            sample_2 = "PeriE",
                                                                            sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPostE.level2 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level2",
                                                                            sample_1 = "PreE",
                                                                            sample_2 = "PostE",
                                                                            sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Foxl2.PeriEvsPostE.level2 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                             cluster_identity = "celltype.level2",
                                                                             sample_1 = "PeriE",
                                                                             sample_2 = "PostE",
                                                                             sample_identity = "Group")

# Benayoun VCD

unique(ovary.VCD@meta.data$Group)         # [1] "PreE"  "PeriE"

ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPeriE.level1 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                          cluster_identity = "celltype.level1",
                                                                          sample_1 = "PreE",
                                                                          sample_2 = "PeriE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPeriE.level2 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                          cluster_identity = "celltype.level2",
                                                                          sample_1 = "PreE",
                                                                          sample_2 = "PeriE",
                                                                          sample_identity = "Group")

# EMTAB12889

unique(ovary.EMTAB12889@meta.data$Group)         # [1] "PeriE" "PostE"

ovary.prop_test.Mouse.EMTAB12889.PeriEvsPostE.level1 <- permutation_test(ovary.Mouse.EMTAB12889.prop_test,
                                                                          cluster_identity = "celltype.level1",
                                                                          sample_1 = "PeriE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.EMTAB12889.PeriEvsPostE.level2 <- permutation_test(ovary.Mouse.EMTAB12889.prop_test,
                                                                          cluster_identity = "celltype.level2",
                                                                          sample_1 = "PeriE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")
# EMTAB11491

unique(ovary.EMTAB11491@meta.data$Group)         # [1] "PeriE" "PostE"

ovary.prop_test.Mouse.EMTAB11491.PeriEvsPostE.level1 <- permutation_test(ovary.Mouse.EMTAB11491.prop_test,
                                                                          cluster_identity = "celltype.level1",
                                                                          sample_1 = "PeriE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.EMTAB11491.PeriEvsPostE.level2 <- permutation_test(ovary.Mouse.EMTAB11491.prop_test,
                                                                          cluster_identity = "celltype.level2",
                                                                          sample_1 = "PeriE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")

# GSE232309

unique(ovary.GSE232309@meta.data$Group)         # [1] "PreE"  "PeriE"

ovary.prop_test.Mouse.GSE232309.PreEvsPeriE.level1 <- permutation_test(ovary.Mouse.GSE232309.prop_test,
                                                                       cluster_identity = "celltype.level1",
                                                                       sample_1 = "PreE",
                                                                       sample_2 = "PeriE",
                                                                       sample_identity = "Group")

ovary.prop_test.Mouse.GSE232309.PreEvsPeriE.level2 <- permutation_test(ovary.Mouse.GSE232309.prop_test,
                                                                       cluster_identity = "celltype.level2",
                                                                       sample_1 = "PreE",
                                                                       sample_2 = "PeriE",
                                                                       sample_identity = "Group")

# GSE290742

unique(ovary.GSE290742@meta.data$Group)         # [1] "PreE"  "PeriE"

ovary.prop_test.Mouse.GSE290742.PreEvsPeriE.level1 <- permutation_test(ovary.Mouse.GSE290742.prop_test,
                                                                       cluster_identity = "celltype.level1",
                                                                       sample_1 = "PreE",
                                                                       sample_2 = "PeriE",
                                                                       sample_identity = "Group")

ovary.prop_test.Mouse.GSE290742.PreEvsPeriE.level2 <- permutation_test(ovary.Mouse.GSE290742.prop_test,
                                                                       cluster_identity = "celltype.level2",
                                                                       sample_1 = "PreE",
                                                                       sample_2 = "PeriE",
                                                                       sample_identity = "Group")

###############################################################################
# Process goat data
###############################################################################

#################################
# Load data
#################################

load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Goat/Goat_PRJNA1010653/2025-07-01/2025-07-01_10x_ovary_Goat_PRJNA1010653_Seurat_object_with_final_annotation.RData")

unique(ovary.Goat@meta.data$Group)       # [1] "young" "aging"

# Modify metadata for consistency

ovary.Goat$celltype.level1[ovary.Goat$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.Goat$celltype.level1[ovary.Goat$celltype.level1 == "Ptprc.pos"] <- "immune"

#################################
# Generate frequency tables
#################################

# Level 1
res_L1_Goat <- compute_freq_table(
  seurat_obj     = ovary.Goat,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_Goat"
)
my.freq.level1[["Goat"]] <- res_L1_Goat$freq

# Level 2
res_L2_Goat <- compute_freq_table(
  seurat_obj     = ovary.Goat,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_Goat"
)
my.freq.level2[["Goat"]] <- res_L2_Goat$freq

#################################
# Assess changes in proportions using scProportionTest
#################################

# Create prop test objects
ovary.Goat.prop_test <- sc_utils(ovary.Goat)

unique(ovary.Goat@meta.data$Group)         # [1] "young" "aging"

ovary.prop_test.Goat.YoungVsAging.level1 <- permutation_test(ovary.Goat.prop_test,
                                                             cluster_identity = "celltype.level1",
                                                             sample_1 = "young",
                                                             sample_2 = "aging",
                                                             sample_identity = "Group")

ovary.prop_test.Goat.YoungVsAging.level2 <- permutation_test(ovary.Goat.prop_test,
                                                             cluster_identity = "celltype.level2",
                                                             sample_1 = "young",
                                                             sample_2 = "aging",
                                                             sample_identity = "Group")

################################################################################
sink(file = paste0(Sys.Date(), "_Celltype_proportion_analysis_session_info.txt"))
sessionInfo()
sink()
