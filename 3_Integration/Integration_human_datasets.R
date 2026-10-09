# Load libraries
library(Seurat)
library(future)

plan(sequential)
options(future.globals.maxSize = 80 * 1024^3)  # 80 GB

#################################
# Integrative ovarian aging analysis
# Integrate Seurat objects - Human data
#################################

#################################
# 1. Load data
#################################

load("./2026-06-15_CZI_Human_Seurat_objects_combined.RData")

#################################
# 2. Update old Seurat objects
#################################

ovary.GSE202601         <- UpdateSeuratObject(ovary.Human.GSE202601.cl)
ovary.GSE255690         <- UpdateSeuratObject(ovary.Human.GSE255690.cl)
ovary.tabula.sapiens    <- UpdateSeuratObject(ovary.Human.TS.cl)

#################################
# 3. Add dataset label
#################################

ovary.GSE202601$Dataset          <- "GSE202601"
ovary.GSE255690$Dataset          <- "GSE255690"
ovary.tabula.sapiens$Dataset     <- "TS"

obj.list <- list(
  GSE202601 = ovary.GSE202601,
  GSE255690 = ovary.GSE255690,
  TS = ovary.tabula.sapiens
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

# 4. Select integration features
features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 2000
)

# 5. Prepare for SCT integration
obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features,
  verbose = FALSE
)

save(obj.list, file = "CZI_human_object_list_post_prepsctintegration.RData")

# 6. Run PCA
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

save(anchors, file = "CZI_anchors_by_reference.RData")

ovary.integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 50,
  preserve.order = TRUE,
  verbose = TRUE
)

# 9. Downstream analysis
DefaultAssay(ovary.integrated) <- "integrated"

ovary.integrated <- RunPCA(ovary.integrated, npcs = 50, verbose = FALSE)
ovary.integrated <- RunUMAP(ovary.integrated, dims = 1:30)
ovary.integrated <- FindNeighbors(ovary.integrated, dims = 1:30)
ovary.integrated <- FindClusters(ovary.integrated, resolution = 0.5)

save(ovary.integrated, file = "CZI_human_integrated_object.RData")

#################################
sink(file = "Integration_human_datasets_session_info.txt")
sessionInf()
sink()
