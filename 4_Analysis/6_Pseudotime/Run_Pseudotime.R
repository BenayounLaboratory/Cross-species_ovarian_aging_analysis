library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)

#############################################
# Integrative ovarian aging analysis
# Pseudotime analysis using monocle 3
# Integrated human and mouse dataset
#############################################

#############################################
# 1. Load data
#############################################

# Human data
load("../2026-06-15_CZI_human_integrated_granulosa_cells_post_integration.Rdata")
human.granulosa <- granulosa

# Mouse data
load("../2026-06-15_CZI_mouse_integrated_granulosa_cells_post_integration.Rdata")
mouse.granulosa <- granulosa

mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "YF"] <- 4
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "OF"] <- 20
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Library %in% c("CTL_3m_30d_1", "CTL_3m_30d_2")] <- 4
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Library %in% c("CTL_3m_90d_1", "CTL_3m_90d_2")] <- 6
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Library %in% c("CTL_10m_30d_1", "CTL_10m_30d_2")] <- 11
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Library %in% c("CTL_10m_90d_1", "CTL_10m_90d_2")] <- 13
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "Young"] <- 4
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "Old"] <- 9
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "Supold"] <- 17
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "4M"] <- 4
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "18M"] <- 18
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "3m"] <- 3
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "9m"] <- 9
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "9mo"] <- 9
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "2m"] <- 2
mouse.granulosa@meta.data$Age[mouse.granulosa@meta.data$Age == "2mo"] <- 2

cat("Human granulosa cells:", ncol(human.granulosa), "\n")     
cat("Mouse granulosa cells:", ncol(mouse.granulosa), "\n")    

save(human.granulosa, mouse.granulosa,
     file = "2026-06-24_CZI_human_mouse_granulosa_cells_subsetted_seurat_objects.RData")

# Re-cluster mouse granulosa — full reprocess

mouse.granulosa <- SCTransform(mouse.granulosa, vst.flavor = "v2", verbose = FALSE)
mouse.granulosa <- RunPCA(mouse.granulosa, verbose = FALSE)

ElbowPlot(mouse.granulosa, ndims = 50)

mouse.granulosa <- RunHarmony(mouse.granulosa,
                              group.by.vars  = "Library",
                              reduction      = "pca",
                              reduction.save = "harmony",
                              verbose        = TRUE)
mouse.granulosa <- RunUMAP(mouse.granulosa, reduction = "harmony", dims = 1:30)
mouse.granulosa <- FindNeighbors(mouse.granulosa, reduction = "harmony", dims = 1:30)
mouse.granulosa <- FindClusters(mouse.granulosa, resolution = 0.5)

DimPlot(mouse.granulosa, group.by = "seurat_clusters", label = TRUE)
DimPlot(mouse.granulosa, group.by = "Library")
DimPlot(mouse.granulosa, group.by = "Age",
        cols = colorRampPalette(c("blue", "yellow", "red"))(
          length(unique(mouse.granulosa@meta.data$Age))
        ))


# Convert to numeric
mouse.granulosa@meta.data$Age <- as.numeric(mouse.granulosa@meta.data$Age)

# Now check cell counts per age, properly sorted
mouse.granulosa@meta.data %>%
  count(Age) %>%
  arrange(Age)

#############################################
# 2. Helper function
#############################################

# Wraps the repeated CDS-building, graph-learning, and pseudotime steps
# so the logic stays in one place and both species run identically.

run_pseudotime <- function(seurat_obj,
                           species,          # "human" or "mouse"
                           youngest_ages,    # numeric vector of root age values
                           age_unit,         # "years" or "months" (for plot labels)
                           n_dims    = 30,
                           cores     = 4) {
  
  prefix <- species  
  
  # ── 1. Build CDS ────────────────────────────────────────────────────────────
  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  cell_meta  <- seurat_obj@meta.data
  gene_meta  <- data.frame(gene_short_name = rownames(counts_mat),
                           row.names       = rownames(counts_mat))
  
  cds <- new_cell_data_set(
    expression_data = counts_mat,
    cell_metadata   = cell_meta,
    gene_metadata   = gene_meta
  )
  
  cds <- estimate_size_factors(cds)
  
  # ── 2. Transfer Seurat UMAP ─────────────────────────────────────────────────
  reducedDims(cds)[["UMAP"]] <- Embeddings(seurat_obj, "umap")
  
  # ── 3. Cluster and learn graph ──────────────────────────────────────────────
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)
  
  # ── 4. Sanity check plots ───────────────────────────────────────────────────
  pdf(paste0("01_", prefix, "_trajectory.pdf"), width = 10, height = 8)
  
  plot_cells(cds,
             color_cells_by          = "Age",
             label_groups_by_cluster = FALSE,
             label_leaves            = FALSE,
             label_branch_points     = FALSE,
             trajectory_graph_color  = "grey40") +
    ggtitle(paste(species, "- Age")) |> print()
  
  plot_cells(cds,
             color_cells_by          = "Dataset",
             label_groups_by_cluster = FALSE,
             label_leaves            = FALSE,
             label_branch_points     = FALSE,
             trajectory_graph_color  = "grey40") +
    ggtitle(paste(species, "- Dataset")) |> print()
  
  dev.off()
  
  saveRDS(cds, paste0("01_", prefix, "_cds_learned.rds"))
  
  # ── 5. UMAP by age ───────────────────────────────────────────────────────
  age_levels <- sort(unique(as.numeric(colData(cds)$Age)))
  age_cols   <- colorRampPalette(c("blue", "yellow", "red"))(length(age_levels))
  
  p1 <- DimPlot(seurat_obj,
                group.by = "Age",
                cols     = age_cols) +
    ggtitle(paste(species, "- UMAP by age")) +
    theme(legend.position = "right")
  
  p2 <- plot_cells(cds,
                   color_cells_by          = "Age",
                   label_groups_by_cluster = FALSE,
                   label_leaves            = TRUE,
                   label_branch_points     = TRUE,
                   trajectory_graph_color  = "grey40") +
    ggtitle(paste(species, "- Trajectory colored by age"))
  
  print(p1)
  print(p2)
  
  pdf(paste0("02_", prefix, "_umap_by_age.pdf"), width = 10, height = 8)
  print(p1)
  print(p2)
  dev.off()
  
  # ── 6. Root selection ───────────────────────────────────────────────────────
  
  cds_shiny <- order_cells(cds)
  
  # ── 7. DEG analysis ─────────────────────────────────────────────────────────
  cds <- cds_shiny
  colData(cds)$Dataset <- factor(colData(cds)$Dataset)
  
  graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph",
                               cores = cores)
  
  write.csv(graph_test_res[order(graph_test_res$q_value), ],
            paste0("04_", prefix, "_graph_test_results.csv"))
  
  sig_genes <- rownames(subset(graph_test_res, q_value < 0.05 & morans_I > 0.1))
  
  gene_fits <- fit_models(
    cds[sig_genes, ],
    model_formula_str = "~pseudotime + Dataset",
    cores             = cores
  )
  
  fit_coefs <- coefficient_table(gene_fits)
  
  pseudotime_degs <- fit_coefs %>%
    filter(term == "pseudotime") %>%
    select(gene_short_name, estimate, std_err, test_val,
           normalized_effect, p_value, q_value) %>%
    arrange(q_value)
   
  write.csv(pseudotime_degs,
            paste0("04_", prefix, "_pseudotime_degs.csv"), row.names = FALSE)
  
  fit_coefs_flat <- fit_coefs %>% select(where(~ !is.list(.x)))
  write.csv(fit_coefs_flat,
            paste0("04_", prefix, "_fit_models_all_terms.csv"), row.names = FALSE)
  
  # ── 7. Visualize top genes ──────────────────────────────────────────────────
  top_genes <- pseudotime_degs %>%
    filter(q_value < 0.05) %>%
    slice_min(q_value, n = 6) %>%
    pull(gene_short_name)
  
  pdf(paste0("04_", prefix, "_top_pseudotime_genes.pdf"), width = 12, height = 8)
  
  if (length(top_genes) > 0) {
    plot_cells(cds,
               genes                 = top_genes,
               label_cell_groups     = FALSE,
               show_trajectory_graph = FALSE) |> print()
  }
  
  as_tibble(colData(cds)) %>%
    filter(is.finite(pseudotime(cds))) %>%
    mutate(Age = factor(Age, levels = sort(unique(as.numeric(Age))))) %>%
    ggplot(aes(x = Age,
               y = pseudotime(cds)[colnames(cds)],
               fill = Age)) +
    geom_violin(scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.3) +
    labs(x = paste0("Age (", age_unit, ")"), y = "Pseudotime",
         title = paste(species, "- Pseudotime by age")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) |>
    print()
  
  dev.off()
  
}

#############################################
# 3. Run pseudotime — Human
#############################################

human_results <- run_pseudotime(
  seurat_obj    = human.granulosa,
  species       = "human",
  youngest_ages = c(23), 
  age_unit      = "years"
)

#############################################
# 4. Run pseudotime — Mouse
#############################################

mouse_results <- run_pseudotime(
  seurat_obj    = mouse.granulosa,
  species       = "mouse",
  youngest_ages = c(3),      
  age_unit      = "months"
)

#############################################
# 5. Cross-species DEG overlap
#############################################

human_degs <- human_results$pseudotime_degs %>%
  filter(q_value < 0.05, abs(normalized_effect) > 0.1)

mouse_degs <- mouse_results$pseudotime_degs %>%
  filter(q_value < 0.05, abs(normalized_effect) > 0.1)

# Mouse symbols are title-case; human are uppercase — harmonize for comparison
conserved <- intersect(toupper(human_degs$gene_short_name),
                       toupper(mouse_degs$gene_short_name))

write.csv(data.frame(gene = conserved),
          "05_conserved_pseudotime_degs.csv", row.names = FALSE)

#############################################
sink(file = "Run_Pseudotime_session_info.txt")
sessionInfo()
sink()