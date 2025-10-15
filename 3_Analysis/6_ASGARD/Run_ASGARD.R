library('Asgard')
library(Seurat)
library(biomaRt)
library(dplyr)
library(limma)
library(cmapR)

rm(list=ls())
set.seed(123456)

################################################################################
# Integrative ovarian aging analysis
# Run ASGARD
################################################################################

# Load drug reference data (shared across all datasets)
my_gene_info <- read.table(file="~/Data/01_Reference/ASGARD/DrugReference/ovary_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info <- read.table(file="~/Data/01_Reference/ASGARD/DrugReference/ovary_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles = GetDrugRef(drug.response.path = '~/Data/01_Reference/ASGARD/DrugReference/ovary_rankMatrix.txt',
                               probe.to.genes = my_gene_info, drug.info = my_drug_info)

GSE92742.gctx.path="~/Data/01_Reference/ASGARD/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"

# Drug score calculation function (shared across all datasets)
CombineP = function (p){
  keep <- (p > 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    res <- list(chisq = NA_real_, df = NA_integer_, p = NA_real_, 
                validp = p[keep])
  }
  else {
    lnp <- log(p[keep])
    chisq <- (-2) * sum(lnp)
    df <- 2 * length(lnp)
    if (length(lnp) != length(p)) {
      warning("Some studies omitted")
    }
    res <- pchisq(chisq,df, lower.tail = FALSE)
  }
  return(res)
}

# Function to calculate drug scores (shared across all datasets)
calculate_drug_scores <- function(Gene.list, Drug.ident.res, cell_metadata, Case, dataset_name) {
  # Calculate cluster proportions
  if (length(Case) > 0) {
    cell_metadata <- subset(cell_metadata, sample %in% Case)
  }
  clustering <- cell_metadata$cluster
  cluster_sizes <- table(clustering)
  cluster_sizes <- cluster_sizes[which(cluster_sizes > 3)]
  cluster_prop <- round(100*cluster_sizes/nrow(cell_metadata), 2) 
  
  # Combine cluster drugs
  drug_list <- data.frame()
  fda_drugs_only <- TRUE
  for (i in names(Drug.ident.res)) {
    ith_cluster_drugs <- Drug.ident.res[[i]]
    drug_names <- ith_cluster_drugs$Drug.name
    ith_cluster_drugs <- ith_cluster_drugs[!duplicated(drug_names), ]
    
    if (fda_drugs_only) {
      drug_names <- intersect(drug_names, FDA.drug)
    }
    
    if (length(drug_names)>0) {
      ith_cluster_drugs <- subset(ith_cluster_drugs, Drug.name %in% drug_names)
      fdrs <- ith_cluster_drugs$FDR
      p_values <- ith_cluster_drugs$P.value
      
      temp <- data.frame(
        drug = drug_names, 
        cluster = i,
        cluster_prop = cluster_prop[i],
        p_value = p_values,
        fdr = fdrs,
        row.names = NULL
      )
      drug_list <- rbind(drug_list, temp)
    }
  }
  drug_list <- unique(drug_list)
  drug_list$weighted_prop <- drug_list$cluster_prop*(-log10(drug_list$fdr))
  drug_list[is.na(drug_list)] <- 0

  drug_coverage <- tapply(drug_list$weighted_prop, drug_list$drug, sum)
  drugs <- rownames(drug_coverage)

  # Combine p-values
  if(length(unique(names(Drug.ident.res)))>1){
    combined_p_values <- tapply(drug_list$p_value, drug_list$drug, CombineP)
  }else{
    combined_p_values <- drug_list$p_value
    names(combined_p_values) <- drug_list$drug
  }

  # Load drug response data
  tissue <- "ovary"
  cell_lines <- subset(cell_data, primary_site == tissue)$cell_id
  drug_metadata_92742 <- col_meta_GSE92742[, c("sig_id", "pert_iname")]
  row.names(drug_metadata_92742) <- drug_metadata_92742$sig_id
  idx <- which(col_meta_GSE92742$cell_id %in% cell_lines & 
                 col_meta_GSE92742$pert_iname %in% drugs)
  sig_ids <- col_meta_GSE92742$sig_id[idx]
  drug_metadata_92742 <- drug_metadata_92742[sig_ids, ]

  exprs <- as.data.frame(parse_gctx(GSE92742.gctx.path, cid=sig_ids)@mat)
  treatments <- colnames(exprs)
  exprs$gene_id <- row.names(exprs)
  tmp <- merge(exprs, gene_meta, by.x="gene_id", by.y="pr_gene_id")
  drug_responses_92742 <- tmp[, c("pr_gene_symbol", treatments)]

  drug_responses <- drug_responses_92742
  row.names(drug_responses) <- drug_responses[, 1]
  drug_responses <- drug_responses[, -1]
  drug_metadata <- drug_metadata_92742

  # Find common DEGs
  common_degs <- list()
  for (i in names(Gene.list)) {
    ith_cluster_degs <- Gene.list[[i]]
    ith_cluster_degs <- subset(ith_cluster_degs, adj.P.Val < 0.1)
    if (length(ith_cluster_degs) > 0) {
      common_degs[[i]] <- rownames(ith_cluster_degs)
    }
  }
  common_degs <- Reduce(intersect, common_degs)

  # Calculate drug scores
  drug_scores <- list()
  for (drug in drugs) {
    treatments <- subset(drug_metadata, pert_iname == drug)$sig_id
    if (length(treatments) > 1) {
      curr_drug_response <- drug_responses[, grepl(paste(treatments, collapse = "|"), colnames(drug_responses))]
      mean_response <- apply(curr_drug_response, 1, mean)
    } else {
      curr_drug_response <- drug_responses[, treatments]
      mean_response <- curr_drug_response
    }
    
    drug_stats <- drug_list[drug_list$drug == drug, ]
    drug_score <- 0
    for (i in names(Gene.list)) {
      cluster_prop <- drug_stats[drug_stats$cluster == i, "cluster_prop"]
      fdr <- drug_stats[drug_stats$cluster == i, "fdr"]
      p_value <- drug_stats[drug_stats$cluster == i, "p_value"]
      
      ith_cluster_degs <- Gene.list[[i]]
      ith_cluster_degs <- subset(ith_cluster_degs, adj.P.Val < 0.05)
      
      treatable_degs <- intersect(row.names(ith_cluster_degs), names(mean_response))
      if (length(treatable_degs > 0)) {
        deg_scores <- ith_cluster_degs[treatable_degs, "score"]
        
        treated_degs <- -deg_scores*mean_response[treatable_degs]
        treated_degs <- treated_degs[which(treated_degs > 0)]
        
        treated_degs_ratio <- length(treated_degs)/length(treatable_degs)
        drug_score <- drug_score +
          (cluster_prop/100)*(-log10(fdr))*treated_degs_ratio
      }
    }
    
    drug_scores[[drug]] <- drug_score
  }
  drug_scores <- t(as.data.frame(drug_scores))

  out <- data.frame(
    Drug.therapeutic.score = drug_scores,
    P.value = combined_p_values[drugs],
    FDR = p.adjust(combined_p_values[drugs], method = "BH")
  )

  Drug.score <- out
  write.table(Drug.score, file = paste0("CZI_drug_score_output_", dataset_name, ".txt"), quote = FALSE, sep = "\t")
  
  return(Drug.score)
}

# Function for limma DEG analysis (shared across all datasets)
run_limma_analysis <- function(seurat_obj, Case, Control, sample_col, group_col) {
  DefaultAssay(seurat_obj) <- "RNA"
  min.cells=3
  Gene.list <- list()
  C_names <- NULL
  
  for(i in unique(seurat_obj@meta.data$celltype.level2)){
    Idents(seurat_obj) <- "celltype.level2"
    c_cells <- subset(seurat_obj, celltype.level2 == i)
    Idents(c_cells) <- group_col
    Samples=c_cells@meta.data
    Controlsample <- row.names(subset(Samples, get(sample_col) %in% Control))
    Casesample <- row.names(subset(Samples, get(sample_col) %in% Case))
    if(length(Controlsample)>min.cells & length(Casesample)>min.cells){
      expr <- as.matrix(c_cells@assays$RNA@data)
      new_expr <- as.matrix(expr[,c(Casesample,Controlsample)])
      new_sample <- data.frame(Samples=c(Casesample,Controlsample),type=c(rep("Case",length(Casesample)),rep("Control",length(Controlsample))))
      row.names(new_sample) <- paste(new_sample$Samples,row.names(new_sample),sep="_")
      expr <- new_expr
      bad <- which(rowSums(expr>0)<3)
      expr <- expr[-bad,]
      mm <- model.matrix(~0 + type, data = new_sample)
      fit <- lmFit(expr, mm)
      contr <- makeContrasts(typeCase - typeControl, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrasts = contr)
      tmp <- eBayes(tmp)
      C_data <- topTable(tmp, sort.by = "P",n = nrow(tmp))
      C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$t,adj.P.Val=C_data$adj.P.Val,P.Value=C_data$P.Value)
      Gene.list[[i]] <- C_data_for_drug
      C_names <- c(C_names,i)
    }
  }
  return(Gene.list)
}

###########################################################
# 1. HUMAN GSE202601 DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/3_Celltype_annotation/Human/Human_GSE202601/2025-06-26/2025-06-26_10x_ovary_Human_GSE202601_celltype_annotated_Seurat_object_combined_final.RData")
ovary.Human.GSE202601.cl@meta.data$sample <- ovary.Human.GSE202601.cl@meta.data$Library

# Sample groups
Old <- c("Human49", "Human51", "Human52", "Human54")
Young <- setdiff(unique(ovary.Human.GSE202601.cl@meta.data$Library), Old)

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.Human.GSE202601.cl, Old, Young, "Library", "Group")

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = Gene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.Human.GSE202601.cl@meta.data
cell_metadata$cluster <- ovary.Human.GSE202601.cl@meta.data$celltype.level2
cell_metadata$sample <- ovary.Human.GSE202601.cl@meta.data$Library

Drug.score.GSE202601 <- calculate_drug_scores(Gene.list, Drug.ident.res, cell_metadata, Old, "Human_GSE202601")

###########################################################
# 2. HUMAN GSE255690 DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-06-26_10x_ovary_Human_GSE255690_celltype_annotated_Seurat_object_combined_final.RData")
ovary.Human.GSE255690.cl@meta.data$sample <- ovary.Human.GSE255690.cl@meta.data$Library
ovary.Human.GSE255690.cl$Group <- "Case"
ovary.Human.GSE255690.cl$Group[ovary.Human.GSE255690.cl$Age != "Old"] <- "Control"

# Sample groups
Old <- c("Old_1", "Old_2", "Old_3")
Young <- c("Young_1", "Young_3", "Young_4", "Middle_2", "Middle_3", "Middle_4")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.Human.GSE255690.cl, Old, Young, "Library", "Group")

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = Gene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.Human.GSE255690.cl@meta.data
cell_metadata$cluster <- ovary.Human.GSE255690.cl@meta.data$celltype.level2
cell_metadata$sample <- ovary.Human.GSE255690.cl@meta.data$Library

Drug.score.GSE255690 <- calculate_drug_scores(Gene.list, Drug.ident.res, cell_metadata, Old, "Human_GSE255690")

###########################################################
# 3. MONKEY DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-07-02_STRTseq_ovary_Monkey_GSE130664_Seurat_object_with_final_annotation.RData")
ovary.Monkey.GSE130664.cl@meta.data$sample <- ovary.Monkey.GSE130664.cl@meta.data$individual
ovary.Monkey.GSE130664.cl$Group <- "Case"
ovary.Monkey.GSE130664.cl$Group[ovary.Monkey.GSE130664.cl$age_group == "Y"] <- "Control"

# Sample groups
Old <- c("OF1", "OF2", "OF3", "OF4")
Young <- c("YF1", "YF2", "YF3", "YF4")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.Monkey.GSE130664.cl, Old, Young, "individual", "Group")

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = Gene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.Monkey.GSE130664.cl@meta.data
cell_metadata$cluster <- ovary.Monkey.GSE130664.cl@meta.data$celltype.level2
cell_metadata$sample <- ovary.Monkey.GSE130664.cl@meta.data$individual

Drug.score.Monkey <- calculate_drug_scores(Gene.list, Drug.ident.res, cell_metadata, Old, "Monkey")

###########################################################
# 4. GOAT DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-07-01_10x_ovary_Goat_PRJNA1010653_Seurat_object_with_final_annotation.RData")
ovary.Goat@meta.data$sample <- ovary.Goat@meta.data$Library

# Sample groups
Old <- c("aging")
Young <- c("young")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.Goat, Old, Young, "Library", "Group")

# Humanize gene names for Goat
goat_symbols <- unique(unlist(lapply(Gene.list, rownames)))
goat_symbols <- goat_symbols[!is.na(goat_symbols) & nzchar(goat_symbols)]

ensembl_goat <- useEnsembl(
  biomart = "genes",
  dataset = "chircus_gene_ensembl",
  mirror  = "useast"
)

attrs <- c("ensembl_gene_id",
           "external_gene_name",
           "hsapiens_homolog_ensembl_gene",
           "hsapiens_homolog_associated_gene_name",
           "hsapiens_homolog_orthology_type")

mapping_raw <- getBM(
  mart       = ensembl_goat,
  filters    = "external_gene_name",
  values     = goat_symbols,
  attributes = attrs
)

map_clean <- mapping_raw %>%
  mutate(
    external_gene_name = trimws(external_gene_name),
    hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name)
  ) %>%
  filter(
    hsapiens_homolog_orthology_type == "ortholog_one2one",
    hsapiens_homolog_associated_gene_name != ""
  ) %>%
  distinct(external_gene_name, .keep_all = TRUE)

goat2human <- map_clean %>%
  select(goat_symbol = external_gene_name,
         human_symbol = hsapiens_homolog_associated_gene_name)

humanize_limma_table <- function(df, goat2human, duplicate_strategy = c("best_padj", "max_abs_score")) {
  duplicate_strategy <- match.arg(duplicate_strategy)
  if (!is.data.frame(df) || !nrow(df)) return(NULL)
  
  for (nm in c("score", "adj.P.Val", "P.Value")) if (!nm %in% colnames(df)) df[[nm]] <- NA_real_
  
  df$goat_symbol <- rownames(df)
  
  joined <- df %>%
    left_join(goat2human, by = c("goat_symbol" = "goat_symbol")) %>%
    filter(!is.na(human_symbol) & human_symbol != "")
  
  if (!nrow(joined)) return(NULL)
  
  resolved <- switch(
    duplicate_strategy,
    best_padj = joined %>%
      group_by(human_symbol) %>%
      arrange(adj.P.Val, .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup(),
    max_abs_score = joined %>%
      mutate(.abs = abs(score)) %>%
      group_by(human_symbol) %>%
      arrange(desc(.abs), .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(-.abs)
  )
  
  out <- resolved %>%
    select(score, adj.P.Val, P.Value, everything()) %>%
    as.data.frame()
  rownames(out) <- out$human_symbol
  out$human_symbol <- NULL
  out$goat_symbol  <- NULL
  out
}

duplicate_strategy <- "best_padj"
hGene.list <- lapply(Gene.list, humanize_limma_table,
                     goat2human = goat2human,
                     duplicate_strategy = duplicate_strategy)

kept <- !vapply(hGene.list, is.null, logical(1))
hGene.list <- hGene.list[kept]

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = hGene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.Goat@meta.data
cell_metadata$cluster <- ovary.Goat@meta.data$celltype.level2
cell_metadata$sample <- ovary.Goat@meta.data$Library

Drug.score.Goat <- calculate_drug_scores(hGene.list, Drug.ident.res, cell_metadata, Old, "Goat")

###########################################################
# 5. MOUSE PRJNA863443 AC DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-06-26_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")
ovary.AC@meta.data$sample <- ovary.AC@meta.data$Library

# Sample groups
Old <- c("OF1", "OF2")
Young <- c("YF1", "YF2")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.AC, Old, Young, "Library", "Group")

# Humanize gene names for Mouse
mouse_symbols <- unique(unlist(lapply(Gene.list, rownames)))
mouse_symbols <- mouse_symbols[!is.na(mouse_symbols) & nzchar(mouse_symbols)]

ensembl_mouse <- useEnsembl(biomart = "genes",
                            dataset = "mmusculus_gene_ensembl",
                            mirror  = "useast")

attrs <- c("ensembl_gene_id",
           "external_gene_name",
           "hsapiens_homolog_ensembl_gene",
           "hsapiens_homolog_associated_gene_name",
           "hsapiens_homolog_orthology_type")

mapping_raw <- getBM(
  mart      = ensembl_mouse,
  filters   = "external_gene_name",
  values    = mouse_symbols,
  attributes = attrs
)

map_clean <- mapping_raw %>%
  mutate(external_gene_name = trimws(external_gene_name),
         hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name)) %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

mouse2human <- map_clean %>%
  select(mouse_symbol = external_gene_name,
         human_symbol = hsapiens_homolog_associated_gene_name)

humanize_limma_table_mouse <- function(df, mouse2human, duplicate_strategy = c("best_padj", "max_abs_score")) {
  duplicate_strategy <- match.arg(duplicate_strategy)
  if (!is.data.frame(df) || !nrow(df)) return(NULL)
  
  for (nm in c("score", "adj.P.Val", "P.Value")) {
    if (!nm %in% colnames(df)) df[[nm]] <- NA_real_
  }
  
  df$mouse_symbol <- rownames(df)
  
  joined <- df %>%
    left_join(mouse2human, by = c("mouse_symbol" = "mouse_symbol")) %>%
    filter(!is.na(human_symbol) & human_symbol != "")
  
  if (!nrow(joined)) return(NULL)
  
  resolved <- switch(
    duplicate_strategy,
    best_padj = joined %>%
      group_by(human_symbol) %>%
      arrange(adj.P.Val, .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup(),
    max_abs_score = joined %>%
      mutate(.abs = abs(score)) %>%
      group_by(human_symbol) %>%
      arrange(desc(.abs), .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(-.abs)
  )
  
  out <- resolved %>%
    select(score, adj.P.Val, P.Value, everything()) %>%
    as.data.frame()
  rownames(out) <- out$human_symbol
  out$human_symbol <- NULL
  out$mouse_symbol <- NULL
  out
}

duplicate_strategy <- "best_padj" 
hGene.list <- lapply(Gene.list, humanize_limma_table_mouse, mouse2human = mouse2human,
                     duplicate_strategy = duplicate_strategy)

kept <- !vapply(hGene.list, is.null, logical(1))
if (any(!kept)) {
  message("No orthologs retained for: ", paste(names(hGene.list)[!kept], collapse = ", "))
  hGene.list <- hGene.list[kept]
}

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = hGene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.AC@meta.data
cell_metadata$cluster <- ovary.AC@meta.data$celltype.level2
cell_metadata$sample <- ovary.AC@meta.data$individual

Drug.score.Mouse_Aging <- calculate_drug_scores(hGene.list, Drug.ident.res, cell_metadata, Old, "Mouse_Aging")

###########################################################
# 6. PRJNA863443 F+ DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-06-27_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")
ovary.Foxl2@meta.data$sample <- ovary.Foxl2@meta.data$Library

# Sample groups
Old <- c("Foxl2_wt_old_0", "Foxl2_wt_old_1", "Foxl2_wt_old_2")
Young <- c("Foxl2_wt_young_1", "Foxl2_wt_young_2")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.Foxl2, Old, Young, "Library", "Group")

# Humanize gene names (reuse function from above)
mouse_symbols <- unique(unlist(lapply(Gene.list, rownames)))
mouse_symbols <- mouse_symbols[!is.na(mouse_symbols) & nzchar(mouse_symbols)]

mapping_raw <- getBM(
  mart      = ensembl_mouse,
  filters   = "external_gene_name",
  values    = mouse_symbols,
  attributes = attrs
)

map_clean <- mapping_raw %>%
  mutate(external_gene_name = trimws(external_gene_name),
         hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name)) %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

mouse2human <- map_clean %>%
  select(mouse_symbol = external_gene_name,
         human_symbol = hsapiens_homolog_associated_gene_name)

hGene.list <- lapply(Gene.list, humanize_limma_table_mouse, mouse2human = mouse2human,
                     duplicate_strategy = duplicate_strategy)

kept <- !vapply(hGene.list, is.null, logical(1))
if (any(!kept)) {
  message("No orthologs retained for: ", paste(names(hGene.list)[!kept], collapse = ", "))
  hGene.list <- hGene.list[kept]
}

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = hGene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.Foxl2@meta.data
cell_metadata$cluster <- ovary.Foxl2@meta.data$celltype.level2
cell_metadata$sample <- ovary.Foxl2@meta.data$individual

Drug.score.Mouse_Foxl2 <- calculate_drug_scores(hGene.list, Drug.ident.res, cell_metadata, Old, "Mouse_Foxl2")

###########################################################
# 7. MOUSE PRJNA863443 VCD DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")
ovary.VCD@meta.data$sample <- ovary.VCD@meta.data$Library

# Sample groups
Old <- c("CTL_10m_90d_1", "CTL_10m_90d_2")
Young <- c("CTL_3m_30d_1", "CTL_3m_90d_1", "CTL_10m_30d_1",  "CTL_3m_30d_2", "CTL_10m_30d_2", "CTL_3m_90d_2")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.VCD, Old, Young, "Library", "Group")

# Humanize gene names
mouse_symbols <- unique(unlist(lapply(Gene.list, rownames)))
mouse_symbols <- mouse_symbols[!is.na(mouse_symbols) & nzchar(mouse_symbols)]

mapping_raw <- getBM(
  mart      = ensembl_mouse,
  filters   = "external_gene_name",
  values    = mouse_symbols,
  attributes = attrs
)

map_clean <- mapping_raw %>%
  mutate(external_gene_name = trimws(external_gene_name),
         hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name)) %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

mouse2human <- map_clean %>%
  select(mouse_symbol = external_gene_name,
         human_symbol = hsapiens_homolog_associated_gene_name)

hGene.list <- lapply(Gene.list, humanize_limma_table_mouse, mouse2human = mouse2human,
                     duplicate_strategy = duplicate_strategy)

kept <- !vapply(hGene.list, is.null, logical(1))
if (any(!kept)) {
  message("No orthologs retained for: ", paste(names(hGene.list)[!kept], collapse = ", "))
  hGene.list <- hGene.list[kept]
}

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = hGene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.VCD@meta.data
cell_metadata$cluster <- ovary.VCD@meta.data$celltype.level2
cell_metadata$sample <- ovary.VCD@meta.data$individual

Drug.score.Mouse_VCD <- calculate_drug_scores(hGene.list, Drug.ident.res, cell_metadata, Old, "Mouse_VCD")

###########################################################
# 8. MOUSE EMTAB12889 DATASET
###########################################################

cat("Processing Mouse EMTAB12889 dataset...\n")

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-07-23_10x_ovary_Mouse_EMTAB12889_Seurat_object_with_final_annotation.RData")
ovary.EMTAB12889@meta.data$sample <- ovary.EMTAB12889@meta.data$Library

# Sample groups
Old <- c("Ovary_15M_1", "Ovary_15M_3", "Ovary_15M_4")
Young <- c("Ovary_9M_1", "Ovary_9M_2", "Ovary_9M_3", "Ovary_9M_4", "Ovary_12M_1", "Ovary_12M_2", "Ovary_12M_3", "Ovary_12M_4", "Ovary_12M_5")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.EMTAB12889, Old, Young, "Library", "Group")

# Humanize gene names
mouse_symbols <- unique(unlist(lapply(Gene.list, rownames)))
mouse_symbols <- mouse_symbols[!is.na(mouse_symbols) & nzchar(mouse_symbols)]

mapping_raw <- getBM(
  mart      = ensembl_mouse,
  filters   = "external_gene_name",
  values    = mouse_symbols,
  attributes = attrs
)

map_clean <- mapping_raw %>%
  mutate(external_gene_name = trimws(external_gene_name),
         hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name)) %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

mouse2human <- map_clean %>%
  select(mouse_symbol = external_gene_name,
         human_symbol = hsapiens_homolog_associated_gene_name)

hGene.list <- lapply(Gene.list, humanize_limma_table_mouse, mouse2human = mouse2human,
                     duplicate_strategy = duplicate_strategy)

kept <- !vapply(hGene.list, is.null, logical(1))
if (any(!kept)) {
  message("No orthologs retained for: ", paste(names(hGene.list)[!kept], collapse = ", "))
  hGene.list <- hGene.list[kept]
}

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = hGene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.EMTAB12889@meta.data
cell_metadata$cluster <- ovary.EMTAB12889@meta.data$celltype.level2
cell_metadata$sample <- ovary.EMTAB12889@meta.data$individual

Drug.score.Mouse_EMTAB12889 <- calculate_drug_scores(hGene.list, Drug.ident.res, cell_metadata, Old, "Mouse_EMTAB12889")

###########################################################
# 9. MOUSE GSE232309 DATASET
###########################################################

# Load data
load("/Volumes/OIProject_II/1_R/0_Annotated_Seurat_objects/2025-07-23_10x_ovary_Mouse_GSE232309_Seurat_object_with_final_annotation.RData")
ovary.GSE232309@meta.data$sample <- ovary.GSE232309@meta.data$Library

# Sample groups
Old <- c("9M_5", "9M_6", "9M_7", "9M_8")
Young <- c("3M_1", "3M_2", "3M_3", "3M_4")

# Limma DEG analysis
Gene.list <- run_limma_analysis(ovary.GSE232309, Old, Young, "Library", "Group")

# Humanize gene names
mouse_symbols <- unique(unlist(lapply(Gene.list, rownames)))
mouse_symbols <- mouse_symbols[!is.na(mouse_symbols) & nzchar(mouse_symbols)]

mapping_raw <- getBM(
  mart      = ensembl_mouse,
  filters   = "external_gene_name",
  values    = mouse_symbols,
  attributes = attrs
)

map_clean <- mapping_raw %>%
  mutate(external_gene_name = trimws(external_gene_name),
         hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name)) %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         hsapiens_homolog_associated_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE)

mouse2human <- map_clean %>%
  select(mouse_symbol = external_gene_name,
         human_symbol = hsapiens_homolog_associated_gene_name)

hGene.list <- lapply(Gene.list, humanize_limma_table_mouse, mouse2human = mouse2human,
                     duplicate_strategy = duplicate_strategy)

kept <- !vapply(hGene.list, is.null, logical(1))
if (any(!kept)) {
  message("No orthologs retained for: ", paste(names(hGene.list)[!kept], collapse = ", "))
  hGene.list <- hGene.list[kept]
}

# Drug repurposing
Drug.ident.res <- GetDrug(
  gene.data = hGene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "all"
)

# Drug score calculation
cell_metadata <- ovary.GSE232309@meta.data
cell_metadata$cluster <- ovary.GSE232309@meta.data$celltype.level2
cell_metadata$sample <- ovary.GSE232309@meta.data$individual

Drug.score.Mouse_GSE232309 <- calculate_drug_scores(hGene.list, Drug.ident.res, cell_metadata, Old, "Mouse_GSE232309")

###########################################################
sink(file = paste0(Sys.Date(), "_Run_ASGARD_session_info.txt"))
sessionInfo()
sink()