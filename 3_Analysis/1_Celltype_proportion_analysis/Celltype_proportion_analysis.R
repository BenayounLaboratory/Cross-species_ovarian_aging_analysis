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
  "Epithelial","Neutrophil", "Myeloid","DC","ILC","NK","NKT","CD8NKT", "CD8 NKT",
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

# Age ranges for human data
# GSE255690: 23 yrs vs. 38 yrs vs. 48 yrs
# GSE202601: 28 yrs vs. 52 yrs

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

#################################
# Assess changes in proportions using scProportionTest
#################################

# Create prop test objects
ovary.Human.GSE255690.prop_test <- sc_utils(ovary.Human.GSE255690.cl)
ovary.Human.GSE202601.prop_test <- sc_utils(ovary.Human.GSE202601.cl)

# GSE255690 
ovary.prop_test.Human.GSE255690.YvsM.level1 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level1",
                                                                sample_1 = "Young",
                                                                sample_2 = "Middle",
                                                                sample_identity = "Age")

ovary.prop_test.Human.GSE255690.YvsO.level1 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level1",
                                                                sample_1 = "Young",
                                                                sample_2 = "Old",
                                                                sample_identity = "Age")

ovary.prop_test.Human.GSE255690.MvsO.level1 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level1",
                                                                sample_1 = "Middle",
                                                                sample_2 = "Old",
                                                                sample_identity = "Age")

ovary.prop_test.Human.GSE255690.YvsM.level2 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level2",
                                                                sample_1 = "Young",
                                                                sample_2 = "Middle",
                                                                sample_identity = "Age")

ovary.prop_test.Human.GSE255690.YvsO.level2 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level2",
                                                                sample_1 = "Young",
                                                                sample_2 = "Old",
                                                                sample_identity = "Age")

ovary.prop_test.Human.GSE255690.MvsO.level2 <- permutation_test(ovary.Human.GSE255690.prop_test,
                                                                cluster_identity = "celltype.level2",
                                                                sample_1 = "Middle",
                                                                sample_2 = "Old",
                                                                sample_identity = "Age")

ovary.prop_test.Human.GSE202601.level1 <- permutation_test(ovary.Human.GSE202601.prop_test,
                                                           cluster_identity = "celltype.level1",
                                                           sample_1 = "Young",
                                                           sample_2 = "Old",
                                                           sample_identity = "Group")

ovary.prop_test.Human.GSE202601.level2 <- permutation_test(ovary.Human.GSE202601.prop_test,
                                                           cluster_identity = "celltype.level2",
                                                           sample_1 = "Young",
                                                           sample_2 = "Old",
                                                           sample_identity = "Group")

# Clear memory
rm(ovary.Human.GSE255690.cl, ovary.Human.GSE202601.cl)

###############################################################################
# Process mouse data
###############################################################################

#################################
# Load data
#################################

load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_Aging/2025-06-26/2025-06-26_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_Foxl2/2025-06-27/2025-06-27_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_VCD/2025-06-26/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_EMTAB12889/2025-07-22/2025-07-23_10x_ovary_Mouse_EMTAB128899_Seurat_object_with_final_annotation.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_GSE232309/2025-06-30/2025-06-30_10x_ovary_Mouse_GSE232309_Seurat_object_with_final_annotation.RData")

# Age ranges for mouse data
# Benayoun lab aging: 4m vs. 20m
# Benayoun lab Foxl2: 4m vs. 9m
# Benayoun lab VCD: 5m vs 7m vs. 12m vs. 14m
# EMTAB128899: 9m vs. 12m vs. 15m
# GSE232309: 3m vs. 9m

# Group ages together:
# 3-7m: Pre-estropausal
# 9-12m: Peri-estropausal
# 14-15m: Post-estropausal
# 20m: Post-post-estropausal

# Rename metadata for consistency
# Benayoun lab aging
ovary.AC$Group <- "PreE"
ovary.AC$Group[ovary.AC$Age == "OF"] <- "PostpostE"

ovary.AC$celltype.level1[ovary.AC$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.AC$celltype.level1[ovary.AC$celltype.level1 == "Ptprc.pos"] <- "immune"

# Benayoun lab Foxl2
ovary.Foxl2$Group <- "PreE"
ovary.Foxl2$Group[ovary.Foxl2$Age == "Old"] <- "PeriE"

ovary.Foxl2$celltype.level1[ovary.Foxl2$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.Foxl2$celltype.level1[ovary.Foxl2$celltype.level1 == "Ptprc.pos"] <- "immune"

# Benayoun lab VCD
ovary.VCD$Group <- "PreE"
ovary.VCD$Group[grepl("10m_30d", ovary.VCD$Library)] <- "PeriE"
ovary.VCD$Group[grepl("10m_90d", ovary.VCD$Library)] <- "PostE"

ovary.VCD$celltype.level1[ovary.VCD$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.VCD$celltype.level1[ovary.VCD$celltype.level1 == "Ptprc.pos"] <- "immune"

# EMTAB128899
ovary.EMTAB128899$Group <- "PeriE"
ovary.EMTAB128899$Group[ovary.EMTAB128899$Age == "15M"] <- "PostE"

ovary.EMTAB128899$celltype.level1[ovary.EMTAB128899$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.EMTAB128899$celltype.level1[ovary.EMTAB128899$celltype.level1 == "Ptprc.pos"] <- "immune"

# EMTAB128899
ovary.GSE232309$Group <- "PreE"
ovary.GSE232309$Group[ovary.GSE232309$Age == "9m"] <- "PeriE"

ovary.GSE232309$celltype.level1[ovary.GSE232309$celltype.level1 == "Ptprc.neg"] <- "nonimmune"
ovary.GSE232309$celltype.level1[ovary.GSE232309$celltype.level1 == "Ptprc.pos"] <- "immune"

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
  seurat_obj     = ovary.Foxl2,
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

res_L1_EMTAB128899 <- compute_freq_table(
  seurat_obj     = ovary.EMTAB128899,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_EMTAB128899"
)
my.freq.level1[["Mouse_EMTAB128899"]] <- res_L1_EMTAB128899$freq

res_L1_GSE232309 <- compute_freq_table(
  seurat_obj     = ovary.GSE232309,
  level_col      = "celltype.level1",
  lib_col        = "Library",
  desired_levels = custom_order_level1,
  dataset_name   = "Mouse_GSE232309"
)
my.freq.level1[["Mouse_GSE232309"]] <- res_L1_GSE232309$freq

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
  seurat_obj     = ovary.Foxl2,
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

res_L2_EMTAB128899 <- compute_freq_table(
  seurat_obj     = ovary.EMTAB128899,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_EMTAB128899"
)
my.freq.level2[["Mouse_EMTAB128899"]] <- res_L2_EMTAB128899$freq

res_L2_GSE232309 <- compute_freq_table(
  seurat_obj     = ovary.GSE232309,
  level_col      = "celltype.level2",
  lib_col        = "Library",
  desired_levels = custom_order_level2,
  dataset_name   = "Mouse_GSE232309"
)
my.freq.level2[["Mouse_GSE232309"]] <- res_L2_GSE232309$freq

#################################
# Assess changes in proportions using scProportionTest
#################################

# Create prop test objects
ovary.Mouse.Benayoun.Aging.prop_test <- sc_utils(ovary.AC)
ovary.Mouse.Benayoun.Foxl2.prop_test <- sc_utils(ovary.Foxl2)
ovary.Mouse.Benayoun.VCD.prop_test   <- sc_utils(ovary.VCD)
ovary.Mouse.EMTAB128899.prop_test    <- sc_utils(ovary.EMTAB128899)
ovary.Mouse.GSE232309.prop_test      <- sc_utils(ovary.GSE232309)

# Benayoun Aging

unique(ovary.AC@meta.data$Group)         # [1] "PreE"      "PostpostE"

ovary.prop_test.Mouse.Benayoun.Aging.PreEvsPostpostE.level1 <- permutation_test(ovary.Mouse.Benayoun.Aging.prop_test,
                                                                                cluster_identity = "celltype.level1",
                                                                                sample_1 = "PreE",
                                                                                sample_2 = "PostpostE",
                                                                                sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Aging.PreEvsPostpostE.level2 <- permutation_test(ovary.Mouse.Benayoun.Aging.prop_test,
                                                                                cluster_identity = "celltype.level2",
                                                                                sample_1 = "PreE",
                                                                                sample_2 = "PostpostE",
                                                                                sample_identity = "Group")
# Benayoun Foxl2

unique(ovary.Foxl2@meta.data$Group)         # [1] "PreE"  "PeriE"

ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPeriE.level1 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level1",
                                                                            sample_1 = "PreE",
                                                                            sample_2 = "PeriE",
                                                                            sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPeriE.level2 <- permutation_test(ovary.Mouse.Benayoun.Foxl2.prop_test,
                                                                            cluster_identity = "celltype.level2",
                                                                            sample_1 = "PreE",
                                                                            sample_2 = "PeriE",
                                                                            sample_identity = "Group")
# Benayoun VCD

unique(ovary.VCD@meta.data$Group)         # [1] "PreE"  "PeriE" "PostE"

ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPeriE.level1 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                          cluster_identity = "celltype.level1",
                                                                          sample_1 = "PreE",
                                                                          sample_2 = "PeriE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPostE.level1 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                          cluster_identity = "celltype.level1",
                                                                          sample_1 = "PreE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.VCD.PeriEvsPostE.level1 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                           cluster_identity = "celltype.level1",
                                                                           sample_1 = "PeriE",
                                                                           sample_2 = "PostE",
                                                                           sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPeriE.level2 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                          cluster_identity = "celltype.level2",
                                                                          sample_1 = "PreE",
                                                                          sample_2 = "PeriE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPostE.level2 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                          cluster_identity = "celltype.level2",
                                                                          sample_1 = "PreE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.Benayoun.VCD.PeriEvsPostE.level2 <- permutation_test(ovary.Mouse.Benayoun.VCD.prop_test,
                                                                           cluster_identity = "celltype.level2",
                                                                           sample_1 = "PeriE",
                                                                           sample_2 = "PostE",
                                                                           sample_identity = "Group")

# EMTAB128899

unique(ovary.EMTAB128899@meta.data$Group)         # [1] "PeriE" "PostE"

ovary.prop_test.Mouse.EMTAB128899.PeriEvsPostE.level1 <- permutation_test(ovary.Mouse.EMTAB128899.prop_test,
                                                                          cluster_identity = "celltype.level1",
                                                                          sample_1 = "PeriE",
                                                                          sample_2 = "PostE",
                                                                          sample_identity = "Group")

ovary.prop_test.Mouse.EMTAB128899.PeriEvsPostE.level2 <- permutation_test(ovary.Mouse.EMTAB128899.prop_test,
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

# Clear memory
rm(ovary.AC, ovary.Foxl2, ovary.VCD, ovary.EMTAB128899, ovary.GSE232309)

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

# Clear memory
rm(ovary.Goat)

###############################################################################
# Generate plots
###############################################################################

#################################
# Human
#################################

generate_permutation_plot_human <- function(..., groups, title, custom_order, max_range = 6) {
  plot_list <- list(...)
  
  # Standardize column names 
  for (i in seq_along(plot_list)) {
    colnames(plot_list[[i]]@results$permutation) <- c("clusters", "Sample_1", "Sample_2", "obs_log2FD", "pval", "FDR", "boot_mean_log2FD", "boot_CI_2.5", "boot_CI_97.5")
  }
  
  # Combine all plot data and assign group/comparison labels
  plot_data <- do.call(rbind, lapply(seq_along(plot_list), function(i) {
    df <- plot_list[[i]]@results$permutation
    df$groups <- rep(groups[i], nrow(df))       # for shape
    df$comparison <- rep(groups[i], nrow(df))   # for tracking
    return(df)
  }))
  
  # Set factor levels
  plot_data$groups <- factor(plot_data$groups, levels = rev(groups))
  plot_data$comparison <- factor(plot_data$comparison, levels = groups)
  plot_data$clusters <- factor(plot_data$clusters, levels = custom_order)
  
  # Replace Inf values with max/min finite values
  max_value <- max(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  min_value <- min(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  
  plot_data$obs_log2FD <- ifelse(
    plot_data$obs_log2FD == Inf, max_value,
    ifelse(plot_data$obs_log2FD == -Inf, min_value, plot_data$obs_log2FD)
  )
  
  # Annotate significance and color group based on direction and comparison
  plot_data <- plot_data %>%
    mutate(
      significance = case_when(
        FDR < 0.05 & boot_mean_log2FD > 0 ~ "Increased",
        FDR < 0.05 & boot_mean_log2FD < 0 ~ "Decreased",
        TRUE ~ "n.s."
      ),
      color_group = case_when(
        significance == "Increased" ~ sub(".*_vs_", "", comparison),  
        significance == "Decreased" ~ sub("_vs_.*", "", comparison), 
        TRUE ~ "n.s."
      )
    )
  
  # Define color palette based on age enrichment
  color_map <- c(
    "Young" = "lightblue1",
    "Middle" = "lightblue3",
    "Old" = "darkblue",
    "n.s." = "gray"
  )
  
  # Define shape mapping per comparison
  shape_mapping <- c(
    "Young_vs_Middle" = 16, 
    "Middle_vs_Old"   = 17, 
    "Young_vs_Old"    = 15
  )
  
  # Plot
  ggplot(plot_data, aes(
    x = clusters, y = obs_log2FD, group = interaction(groups, clusters)
  )) +
    theme_bw() +
    geom_pointrange(
      aes(
        ymin = boot_CI_2.5,
        ymax = boot_CI_97.5,
        color = color_group,
        shape = groups
      ),
      position = position_dodge2(width = 0.8, padding = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_mapping) +
    coord_flip() +
    ylim(-max_range, max_range) +
    ggtitle(title) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    )
}

# Generate permutation plots with custom orders
# Level 1
plot_Level1_GSE255690 <- generate_permutation_plot_human(
  ovary.prop_test.Human.GSE255690.YvsM.level1,
  ovary.prop_test.Human.GSE255690.MvsO.level1,
  ovary.prop_test.Human.GSE255690.YvsO.level1,
  groups = c("Young_vs_Middle", "Middle_vs_Old", "Young_vs_Old"),
  title = "Human_GSE255690: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

plot_Level1_GSE202601 <- generate_permutation_plot_human(
  ovary.prop_test.Human.GSE202601.level1,
  groups = c("Young_vs_Old"),
  title = "Human_GSE202601: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

# Level 2
plot_Level2_GSE255690 <- generate_permutation_plot_human(
  ovary.prop_test.Human.GSE255690.YvsM.level2,
  ovary.prop_test.Human.GSE255690.MvsO.level2,
  ovary.prop_test.Human.GSE255690.YvsO.level2,
  groups = c("Young_vs_Middle", "Middle_vs_Old", "Young_vs_Old"),
  title = "Human_GSE255690: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

plot_Level2_GSE202601 <- generate_permutation_plot_human(
  ovary.prop_test.Human.GSE202601.level2,
  groups = c("Young_vs_Old"),
  title = "Human_GSE202601: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

#################################
# Mouse
#################################

generate_permutation_plot_mouse <- function(..., groups, title, custom_order, max_range = 6) {
  plot_list <- list(...)
  
  # Standardize column names
  for (i in seq_along(plot_list)) {
    colnames(plot_list[[i]]@results$permutation) <- c("clusters", "Sample_1", "Sample_2", "obs_log2FD", "pval", "FDR", "boot_mean_log2FD", "boot_CI_2.5", "boot_CI_97.5")
  }
  
  # Combine all plot data and assign group/comparison labels
  plot_data <- do.call(rbind, lapply(seq_along(plot_list), function(i) {
    df <- plot_list[[i]]@results$permutation
    df$groups <- rep(groups[i], nrow(df))       # for shape
    df$comparison <- rep(groups[i], nrow(df))   # for tracking
    return(df)
  }))
  
  # Set factor levels
  plot_data$groups <- factor(plot_data$groups, levels = rev(groups))
  plot_data$comparison <- factor(plot_data$comparison, levels = groups)
  plot_data$clusters <- factor(plot_data$clusters, levels = custom_order)
  
  # Replace Inf values with max/min finite values
  max_value <- max(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  min_value <- min(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  
  plot_data$obs_log2FD <- ifelse(
    plot_data$obs_log2FD == Inf, max_value,
    ifelse(plot_data$obs_log2FD == -Inf, min_value, plot_data$obs_log2FD)
  )
  
  # Annotate significance and color group based on direction and comparison
  plot_data <- plot_data %>%
    mutate(
      significance = case_when(
        FDR < 0.05 & boot_mean_log2FD > 0 ~ "Increased",
        FDR < 0.05 & boot_mean_log2FD < 0 ~ "Decreased",
        TRUE ~ "n.s."
      ),
      color_group = case_when(
        significance == "Increased" ~ sub(".*_vs_", "", comparison),   # Group 2 (older)
        significance == "Decreased" ~ sub("_vs_.*", "", comparison),   # Group 1 (younger)
        TRUE ~ "n.s."
      )
    )
  
  # Define color palette based on age enrichment
  color_map <- c(
    "PreE" = "deeppink1",
    "PeriE" = "deeppink2",
    "PostE" = "deeppink3",
    "PostpostE" = "deeppink4",
    "n.s." = "gray"
  )
  
  # Define shape mapping per comparison
  shape_mapping <- c(
    "PreE_vs_PeriE"     = 16, 
    "PreE_vs_PostE"     = 17, 
    "PreE_vs_PostpostE" = 15,
    "PeriE_vs_PostE"    = 18
  )
  
  # Plot
  ggplot(plot_data, aes(
    x = clusters, y = obs_log2FD, group = interaction(groups, clusters)
  )) +
    theme_bw() +
    geom_pointrange(
      aes(
        ymin = boot_CI_2.5,
        ymax = boot_CI_97.5,
        color = color_group,
        shape = groups
      ),
      position = position_dodge2(width = 0.8, padding = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_mapping) +
    coord_flip() +
    ylim(-max_range, max_range) +
    ggtitle(title) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    )
}

# Generate permutation plots with custom orders
# Level 1
plot_Level1_Benayoun_Aging <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.Benayoun.Aging.PreEvsPostpostE.level1,
  groups = c("PreE_vs_PostpostE"),
  title = "Mouse_Benayoun_Aging: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

plot_Level1_Benayoun_Foxl2 <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPeriE.level1,
  groups = c("PreE_vs_PeriE"),
  title = "Mouse_Benayoun_Foxl2: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

plot_Level1_Benayoun_VCD <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPeriE.level1,
  ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPostE.level1,
  ovary.prop_test.Mouse.Benayoun.VCD.PeriEvsPostE.level1,
  groups = c("PreE_vs_PeriE", "PreE_vs_PostE", "PeriE_vs_PostE"),
  title = "Mouse_Benayoun_VCD: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

plot_Level1_EMTAB128899 <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.EMTAB128899.PeriEvsPostE.level1,
  groups = c("PeriE_vs_PostE"),
  title = "Mouse_EMTAB128899: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

plot_Level1_GSE232309 <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.GSE232309.PreEvsPeriE.level1,
  groups = c("PreE_vs_PeriE"),
  title = "Mouse_GSE232309: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

# Level 2
plot_Level2_Benayoun_Aging <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.Benayoun.Aging.PreEvsPostpostE.level2,
  groups = c("PreE_vs_PostpostE"),
  title = "Mouse_Benayoun_Aging: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

plot_Level2_Benayoun_Foxl2 <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.Benayoun.Foxl2.PreEvsPeriE.level2,
  groups = c("PreE_vs_PeriE"),
  title = "Mouse_Benayoun_Foxl2: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

plot_Level2_Benayoun_VCD <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPeriE.level2,
  ovary.prop_test.Mouse.Benayoun.VCD.PreEvsPostE.level2,
  ovary.prop_test.Mouse.Benayoun.VCD.PeriEvsPostE.level2,
  groups = c("PreE_vs_PeriE", "PreE_vs_PostE", "PeriE_vs_PostE"),
  title = "Mouse_Benayoun_VCD: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

plot_Level2_EMTAB128899 <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.EMTAB128899.PeriEvsPostE.level2,
  groups = c("PeriE_vs_PostE"),
  title = "Mouse_EMTAB128899: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

plot_Level2_GSE232309 <- generate_permutation_plot_mouse(
  ovary.prop_test.Mouse.GSE232309.PreEvsPeriE.level2,
  groups = c("PreE_vs_PeriE"),
  title = "Mouse_GSE232309: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 10
)

generate_permutation_plot_goat <- function(..., groups, title, custom_order, max_range = 6) {
  plot_list <- list(...)
  
  # Standardize column names
  for (i in seq_along(plot_list)) {
    colnames(plot_list[[i]]@results$permutation) <- c("clusters", "Sample_1", "Sample_2", "obs_log2FD", "pval", "FDR", "boot_mean_log2FD", "boot_CI_2.5", "boot_CI_97.5")
  }
  
  # Combine all plot data and assign group/comparison labels
  plot_data <- do.call(rbind, lapply(seq_along(plot_list), function(i) {
    df <- plot_list[[i]]@results$permutation
    df$groups <- rep(groups[i], nrow(df))       # for shape
    df$comparison <- rep(groups[i], nrow(df))   # for tracking
    return(df)
  }))
  
  # Set factor levels
  plot_data$groups <- factor(plot_data$groups, levels = rev(groups))
  plot_data$comparison <- factor(plot_data$comparison, levels = groups)
  plot_data$clusters <- factor(plot_data$clusters, levels = custom_order)
  
  # Replace Inf values with max/min finite values
  max_value <- max(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  min_value <- min(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  
  plot_data$obs_log2FD <- ifelse(
    plot_data$obs_log2FD == Inf, max_value,
    ifelse(plot_data$obs_log2FD == -Inf, min_value, plot_data$obs_log2FD)
  )
  
  # Annotate significance and color group based on direction and comparison
  plot_data <- plot_data %>%
    mutate(
      significance = case_when(
        FDR < 0.05 & boot_mean_log2FD > 0 ~ "Increased",
        FDR < 0.05 & boot_mean_log2FD < 0 ~ "Decreased",
        TRUE ~ "n.s."
      ),
      color_group = case_when(
        significance == "Increased" ~ sub(".*_vs_", "", comparison),  
        significance == "Decreased" ~ sub("_vs_.*", "", comparison),  
        TRUE ~ "n.s."
      )
    )
  
  # Define color palette based on age enrichment
  color_map <- c(
    "young" = "gold",
    "aging" = "orange",
    "n.s." = "gray"
  )
  
  # Define shape mapping per comparison
  shape_mapping <- c(
    "young_vs_aging" = 16
  )
  
  # Plot
  ggplot(plot_data, aes(
    x = clusters, y = obs_log2FD, group = interaction(groups, clusters)
  )) +
    theme_bw() +
    geom_pointrange(
      aes(
        ymin = boot_CI_2.5,
        ymax = boot_CI_97.5,
        color = color_group,
        shape = groups
      ),
      position = position_dodge2(width = 0.8, padding = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_mapping) +
    coord_flip() +
    ylim(-max_range, max_range) +
    ggtitle(title) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    )
}

# Generate permutation plots with custom orders
# Level 1
plot_Level1_Goat <- generate_permutation_plot_goat(
  ovary.prop_test.Goat.YoungVsAging.level1,
  groups = c("young_vs_aging"),
  title = "Goat: Level 1 Cell Proportion Changes",
  custom_order = custom_order_level1,
  max_range = 2
)

plot_Level2_Goat <- generate_permutation_plot_goat(
  ovary.prop_test.Goat.YoungVsAging.level2,
  groups = c("young_vs_aging"),
  title = "Goat: Level 2 Cell Proportion Changes",
  custom_order = custom_order_level2,
  max_range = 7
)

################################################################################
# Save plot
################################################################################

pdf(paste(Sys.Date(),"CZI_scProportionTest_data_plots_combined.pdf",sep = "_"), height = 20, width = 50)
gridExtra::grid.arrange(plot_Level1_GSE255690,
                        plot_Level1_GSE202601,
                        plot_Level1_Benayoun_Aging,
                        plot_Level1_Benayoun_Foxl2,
                        plot_Level1_Benayoun_VCD,
                        plot_Level1_EMTAB128899,
                        plot_Level1_GSE232309,
                        plot_Level1_Goat,
                        plot_Level2_GSE255690,
                        plot_Level2_GSE202601,
                        plot_Level2_Benayoun_Aging, 
                        plot_Level2_Benayoun_Foxl2,
                        plot_Level2_Benayoun_VCD,
                        plot_Level2_EMTAB128899,
                        plot_Level2_GSE232309,      
                        plot_Level2_Goat,           
                        nrow = 2, ncol = 8)
dev.off()


pdf(paste(Sys.Date(),"CZI_scProportionTest_data_plots_LEVEL1_combined.pdf",sep = "_"), height = 3, width = 50)
gridExtra::grid.arrange(plot_Level1_GSE255690,
                        plot_Level1_GSE202601,
                        plot_Level1_Benayoun_Aging,
                        plot_Level1_Benayoun_Foxl2,
                        plot_Level1_Benayoun_VCD,
                        plot_Level1_EMTAB128899,
                        plot_Level1_GSE232309,
                        plot_Level1_Goat,
                        nrow = 1, ncol = 8)
dev.off()

pdf(paste(Sys.Date(),"CZI_scProportionTest_data_plots_LEVEL2_combined.pdf",sep = "_"), height = 10, width = 50)
gridExtra::grid.arrange(plot_Level2_GSE255690,
                        plot_Level2_GSE202601,
                        plot_Level2_Benayoun_Aging, 
                        plot_Level2_Benayoun_Foxl2,
                        plot_Level2_Benayoun_VCD,
                        plot_Level2_EMTAB128899,
                        plot_Level2_GSE232309,      
                        plot_Level2_Goat,           
                        nrow = 1, ncol = 8)
dev.off()

################################################################################
sink(file = paste0(Sys.Date(), "_Celltype_proportion_analysis_session_info.txt"))
sessionInfo()
sink()
