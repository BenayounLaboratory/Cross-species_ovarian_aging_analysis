options(stringsAsFactors = F)

library(Seurat)
library(SeuratWrappers)
library('Augur')
library(viridis)

###############################################################################
# Integrative ovarian aging analysis
# Perform Augur
###############################################################################

################################################################################
# SECTION 1: Human datasets analysis
################################################################################

# Load Human datasets
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Human/Human_GSE202601/2025-06-26/2025-06-26_10x_ovary_Human_GSE202601_celltype_annotated_Seurat_object_combined_final.RData")
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Human/Human_GSE255690/2025-06-23/2025-06-26_10x_ovary_Human_GSE255690_celltype_annotated_Seurat_object_combined_final.RData")

# Run Augur for Human GSE202601
augur.ovary.Human.GSE202601 <- calculate_auc(as.matrix(ovary.Human.GSE202601.cl@assays$SCT@data),
                                           ovary.Human.GSE202601.cl@meta.data, 
                                           cell_type_col = "celltype.level2", 
                                           label_col = "Group",
                                           n_threads = 3)

# Save Human GSE202601 results
save(augur.ovary.Human.GSE202601, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE202601_AUGUR_object.RData"))

# Run Augur for Human GSE255690
augur.ovary.Age <- calculate_auc(as.matrix(ovary.Human.GSE255690.cl@assays$SCT@data),
                               ovary.Human.GSE255690.cl@meta.data, 
                               cell_type_col = "celltype.level2", 
                               label_col = "Age",
                               n_threads = 3)

# Save Human GSE255690 results
save(augur.ovary.Age, file = paste0(Sys.Date(),"_10x_ovary_Human_GSE255690_AUGUR_object_Age.RData"))

################################################################################
# SECTION 2: Mouse datasets analysis
################################################################################

# Benayoun aging dataset
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_Aging/2025-06-26/2025-06-26_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

# Benayoun Foxl2 dataset
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_Foxl2/2025-06-27/2025-06-27_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

# Add Group metadata for Foxl2 dataset
Group <- rep("NA", length(colnames(ovary.Foxl2@assays$RNA)))
Group[grep("wt_young", colnames(ovary.Foxl2@assays$RNA))]  <- "wt_young"
Group[grep("wt_old", colnames(ovary.Foxl2@assays$RNA))]    <- "wt_midage"
Group[grep("het_young", colnames(ovary.Foxl2@assays$RNA))] <- "het_young"
Group[grep("het_old", colnames(ovary.Foxl2@assays$RNA))]   <- "het_midage"

Group <- data.frame(Group)
rownames(Group) <- colnames(ovary.Foxl2@assays$RNA)

# Update Seurat object with metadata
ovary.Foxl2 <- AddMetaData(object = ovary.Foxl2, metadata = as.vector(Group), col.name = "Group")

# Benayoun VCD dataset
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_Benayoun_lab_VCD/2025-06-26/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

# Add Group metadata for VCD dataset
ovary.VCD@meta.data$Group <- paste0(ovary.VCD@meta.data$Age, "_", ovary.VCD@meta.data$Duration)

# GSE232309 dataset
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_GSE232309/2025-06-30/2025-06-30_10x_ovary_Mouse_GSE232309_Seurat_object_with_final_annotation.RData")

# EMTAB12889 dataset
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Mouse/Mouse_EMTAB12889/2025-06-30/2025-06-30_10x_ovary_Mouse_EMTAB12889_Seurat_object_with_final_annotation.RData")

# Run Augur analyses for each Mouse dataset
augur.ovary.AC <- calculate_auc(as.matrix(ovary.AC@assays$SCT@data),
                              ovary.AC@meta.data, 
                              cell_type_col = "celltype.level2", 
                              label_col = "Age",
                              n_threads = 3)

save(augur.ovary.AC, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Aging_AUGUR_object.RData"))

augur.ovary.Foxl2 <- calculate_auc(as.matrix(ovary.Foxl2@assays$SCT@data),
                                  ovary.Foxl2@meta.data, 
                                  cell_type_col = "celltype.level2", 
                                  label_col = "Group",
                                  n_threads = 3)

save(augur.ovary.Foxl2, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_wt_Augur_object.RData"))

augur.ovary.VCD.Age.groups <- calculate_auc(as.matrix(ovary.VCD@assays$SCT@data),
                                           ovary.VCD@meta.data, 
                                           cell_type_col = "celltype.level2", 
                                           label_col = "Age",
                                           n_threads = 3)

save(augur.ovary.VCD.Age.groups, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_CTL_only_4-6m_vs_11-13m_AUGUR_object.RData"))

augur.ovary.VCD.Age <- calculate_auc(as.matrix(ovary.VCD@assays$SCT@data),
                                     ovary.VCD@meta.data, 
                                     cell_type_col = "celltype.level2", 
                                     label_col = "Group",
                                     n_threads = 3)

save(augur.ovary.VCD.Age, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_CTL_only_by_age_AUGUR_object.RData"))

augur.ovary.GSE232309 <- calculate_auc(as.matrix(ovary.GSE232309@assays$SCT@data),
                                       ovary.GSE232309@meta.data, 
                                       cell_type_col = "celltype.level2", 
                                       label_col = "Age",
                                       n_threads = 3)

save(augur.ovary.GSE232309, file = paste0(Sys.Date(),"_10x_ovary_Mouse_GSE232309_all_cells_AUGUR_object.RData"))

augur.all.age <- calculate_auc(as.matrix(ovary.EMTAB12889@assays$SCT@data),
                              ovary.EMTAB12889@meta.data,
                              cell_type_col = "celltype.level2",
                              label_col = "Age",
                              n_threads = 4)

save(augur.all.age, file = paste0(Sys.Date(), "_Mouse_EMTAB12889_AUGUR_Multiclass_Age.RData"))

################################################################################
# SECTION 3: Monkey dataset analysis
################################################################################

# Load Monkey dataset
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Monkey/Monkey_GSE130664/2025-07-02/2025-07-02_STRTseq_ovary_Monkey_GSE130664_Seurat_object_with_final_annotation.RData")

# Run Augur for Monkey dataset
augur.ovary.monkey <- calculate_auc(as.matrix(ovary.Monkey.GSE130664.cl@assays$SCT@data),
                                   ovary.Monkey.GSE130664.cl@meta.data, 
                                   cell_type_col = "celltype.level2", 
                                   label_col = "age_group",
                                   n_threads = 3)

# Save Monkey results
save(augur.ovary.monkey, file = paste0(Sys.Date(),"_STRTseq_ovary_Monkey_GSE130664_all_cells_AUGUR_object.RData"))

################################################################################
# SECTION 4: Goat dataset analysis
################################################################################

# Load Goat dataset
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-07-01_10x_ovary_Goat_PRJNA1010653_Seurat_object_with_final_annotation.RData")

# Run Augur for Goat dataset
augur.ovary.Goat <- calculate_auc(as.matrix(ovary.Goat@assays$SCT@data),
                                 ovary.Goat@meta.data, 
                                 cell_type_col = "celltype.level2", 
                                 label_col = "Library",
                                 n_threads = 3)

# Save Goat results
save(augur.ovary.Goat, file = paste0(Sys.Date(),"_10x_ovary_Goat_PRJNA1010653_AUGUR_object.RData"))

################################################################################
# Save session info
sink(file = paste0(Sys.Date(),"_Augur_analysis_session_info.txt"))
sessionInfo()
sink()