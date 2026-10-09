# Load libraries
library(Seurat)
library(future)

plan(sequential)
options(future.globals.maxSize = 80 * 1024^3) 

#################################
# Integrative ovarian aging analysis
# Integrate Seurat objects - Mouse data
#################################

#################################
# 1. Load data
#################################

load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2025-11-05_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2026-04-28_10x_ovary_Foxl2_wt_only_Seurat_object_celltypes_annotated.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2026-05-01_10x_ovary_Mouse_GSE290742_Seurat_object_with_final_annotation.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2026-04-22_10x_ovary_Mouse_EMTAB11491_Seurat_object_with_final_annotation.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2025-06-30_10x_ovary_Mouse_GSE232309_Seurat_object_with_final_annotation.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2025-06-30_10x_ovary_Mouse_EMTAB12889_Seurat_object_with_final_annotation.RData")

#################################
# 2. Update old Seurat objects
#################################

ovary.AC          <- UpdateSeuratObject(ovary.AC)
ovary.VCD         <- UpdateSeuratObject(ovary.VCD)
ovary.Foxl2       <- UpdateSeuratObject(ovary.Foxl2.wt)
ovary.EMTAB11491  <- UpdateSeuratObject(ovary.EMTAB11491)
ovary.EMTAB12889  <- UpdateSeuratObject(ovary.EMTAB12889)
ovary.GSE232309   <- UpdateSeuratObject(ovary.GSE232309)
ovary.GSE290742   <- UpdateSeuratObject(ovary.GSE290742)

#################################
# 3. Make sure each object has a dataset label
#################################

ovary.AC$Dataset           <- "AC"
ovary.VCD$Dataset          <- "VCD"
ovary.Foxl2$Dataset        <- "Foxl2"
ovary.EMTAB11491$Dataset   <- "EMTAB11491"
ovary.EMTAB12889$Dataset   <- "EMTAB12889"
ovary.GSE232309$Dataset    <- "GSE232309"
ovary.GSE290742$Dataset    <- "GSE290742"

obj.list <- list(
  AC = ovary.AC,
  VCD = ovary.VCD,
  Foxl2 = ovary.Foxl2,
  EMTAB11491 = ovary.EMTAB11491,
  EMTAB12889 = ovary.EMTAB12889,
  GSE232309 = ovary.GSE232309,
  GSE290742 = ovary.GSE290742
)

obj.list <- lapply(obj.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(
    x,
    method = "glmGamPoi",
    return.only.var.genes = TRUE,
    conserve.memory = TRUE,
    verbose = FALSE
  )
  x <- RunPCA(x, npcs = 50, verbose = FALSE)
  x
})

# Save intermediate file
save(obj.list, file = "CZI_mouse_combined_object_list_post_SCTransform.RData")

# 4. Select integration features
features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 2000
)

# rm(ovary.AC, ovary.Foxl2, ovary.Foxl2.wt, ovary.VCD, ovary.EMTAB11491, ovary.EMTAB12889, ovary.GSE232309, ovary.GSE290742)

# 5. Prepare for SCT integration
obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features,
  verbose = FALSE
)

save(obj.list, file = "MeMo_object_list_post_prepsctintegration.RData")

# 6. Run PCA using the same integration features
obj.list <- lapply(obj.list, function(x) {
  x <- RunPCA(
    x,
    features = features,
    npcs = 50,
    verbose = FALSE
  )
  return(x)
})

save(obj.list, file = "CZI_object_list_post_PCA.RData")

# 7. Find RPCA anchors
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30,
  k.anchor = 5,
  verbose = FALSE
)

ovary.integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 50,
  preserve.order = TRUE,
  verbose = TRUE
)

# 8. Integrate data
ovary.integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30,
  verbose = FALSE
)

# 9. Downstream analysis
DefaultAssay(ovary.integrated) <- "integrated"

ovary.integrated <- RunPCA(ovary.integrated, npcs = 50, verbose = FALSE)
ovary.integrated <- RunUMAP(ovary.integrated, dims = 1:30)
ovary.integrated <- FindNeighbors(ovary.integrated, dims = 1:30)
ovary.integrated <- FindClusters(ovary.integrated, resolution = 0.5)

save(ovary.integrated, file = "CZI_integrated_object.RData")


#################################
sink(file = "Integration_mouse_datasets_session_info.txt")
sessionInf()
sink()

