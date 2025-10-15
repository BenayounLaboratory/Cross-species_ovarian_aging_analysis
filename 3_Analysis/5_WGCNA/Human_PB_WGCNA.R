options(stringsAsFactors = FALSE)

library(WGCNA)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

rm(list = ls())

################################################################################
# Integrative ovarian aging analysis
# Run WGCNA
# Process human datasets
################################################################################

################################################################################
# 1. Load PB counts table
################################################################################

load("/Volumes/OIProject_II/1_R/4_Analysis/7_PB_datasets/2025-08-06_10x_ovary_Human_GSE202601_PB_VST_counts.RData")
vst.counts.GSE202601 <- vst.counts
rm(vst.counts)

load("/Volumes/OIProject_II/1_R/4_Analysis/7_PB_datasets/2025-08-06_PB_VST_counts.RData")
vst.counts.GSE255690 <- vst.counts
rm(vst.counts)

# Assess detected cell types
all_cell_types < c("Granulosa", "Theca", "Stroma")

cell_types_both  <- intersect(cts1, cts2)    # consensus
cell_types_only1 <- setdiff(cts1, cts2)      # single on GSE202601
cell_types_only2 <- setdiff(cts2, cts1)      # single on GSE255690

# Create folders
for (cell in all_cell_types) {
  dir.create(file.path(getwd(), cell), showWarnings = FALSE, recursive = TRUE)
}

################################################################################
# 2. Build traits
################################################################################

# Grab sample names
smp2 <- colnames(vst.counts.GSE255690[[names(vst.counts.GSE255690)[1]]])
smp1 <- colnames(vst.counts.GSE202601[[names(vst.counts.GSE202601)[1]]])

# GSE255690: "Young_*", "Middle_*", "Old_*"
group2 <- ifelse(grepl("^Young",  smp2, ignore.case = TRUE), "Young",
                 ifelse(grepl("^Middle", smp2, ignore.case = TRUE), "Middle",
                        ifelse(grepl("^Old",    smp2, ignore.case = TRUE), "Old", NA)))
traits.GSE255690 <- data.frame(
  AgeGroup = factor(group2, levels = c("Young","Middle","Old")),
  AgeScore = as.numeric(factor(group2, levels = c("Young","Middle","Old"))) - 1, # 0,1,2
  row.names = smp2,
  check.names = FALSE
)

# GSE202601: names containing "Human5" are Old; others Young
group1 <- ifelse(grepl("Human5", smp1, ignore.case = TRUE), "Old", "Young")
traits.GSE202601 <- data.frame(
  AgeGroup = factor(group1, levels = c("Young","Old")),
  AgeScore = ifelse(group1 == "Young", 0, 2),  # 0,2 to roughly align with GSE255690 scaling
  row.names = smp1,
  check.names = FALSE
)

################################################################################
# 3. Define functions
################################################################################

meta_beta <- function(beta, se){
  w <- 1/(se^2); b <- sum(w*beta)/sum(w)
  se_b <- sqrt(1/sum(w))
  z <- b/se_b; p <- 2*pnorm(-abs(z))
  data.frame(meta_beta=b, meta_se=se_b, meta_z=z, meta_p=p)
}

align_traits <- function(traits_master, sample_names, label){

  tr <- traits_master[intersect(rownames(traits_master), sample_names), , drop = FALSE]
  missing <- setdiff(sample_names, rownames(tr))
  if (length(missing) > 0) {
    add <- as.data.frame(matrix(NA, nrow = length(missing), ncol = ncol(traits_master)),
                         stringsAsFactors = FALSE)
    colnames(add) <- colnames(traits_master)
    rownames(add) <- missing
    tr <- rbind(tr, add)
  }
  tr <- tr[sample_names, , drop = FALSE]

  return(tr)
}

fit_ME_age_safe <- function(ME, AgeScore, covars=NULL, tag=""){

  df <- data.frame(ME = as.numeric(ME), AgeScore = as.numeric(AgeScore))
  if (!is.null(covars) && ncol(covars) > 0) {
    covars <- covars[rownames(df), , drop = FALSE]
    df <- cbind(df, covars)
  }
  m <- lm(ME ~ AgeScore + ., data = df)
  co <- summary(m)$coef

  c(beta = co["AgeScore","Estimate"],
    se   = co["AgeScore","Std. Error"],
    p    = co["AgeScore","Pr(>|t|)"])
}

safe_goodSamplesGenes <- function(datExpr){
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  datExpr[gsg$goodSamples, gsg$goodGenes, drop=FALSE]
}

pick_power_multiple_safe <- function(datExpr.H1, datExpr.H2, powers = 1:30){
  sft1 <- pickSoftThreshold(datExpr.H1, powerVector=powers,
                            networkType="signed", corFnc="bicor", verbose=5)
  sft2 <- pickSoftThreshold(datExpr.H2, powerVector=powers,
                            networkType="signed", corFnc="bicor", verbose=5)
  choose_power <- function(sft){
    pe <- sft$powerEstimate
    if(!is.na(pe)) return(pe)
    fi <- sft$fitIndices
    idx <- which(fi[, "SFT.R.sq"] >= 0.8)
    if(length(idx) > 0) return(fi[min(idx), "Power"])
    8
  }
  p1 <- choose_power(sft1)
  p2 <- choose_power(sft2)
  min(p1, p2)
}

pick_power_single_safe <- function(datExpr, powers = 1:30){
  sft <- pickSoftThreshold(datExpr, powerVector=powers,
                           networkType="signed", corFnc="bicor", verbose=5)
  p <- sft$powerEstimate
  if(is.na(p)){
    fit <- sft$fitIndices
    idx <- which(fit[, "SFT.R.sq"] >= 0.8)
    p <- if(length(idx)>0) fit[min(idx), "Power"] else 8
  }
  p
}

plot_ME_age_heatmap <- function(MEs, AgeScore, outfile, title){
  
  # Correlate each module eigengene to AgeScore
  r <- as.numeric(cor(MEs, AgeScore, use = "p"))
  names(r) <- colnames(MEs)
  
  # P-values per module
  p <- corPvalueStudent(r, nrow(MEs))
  
  # Make a 1 x K matrix (row = "AgeScore", columns = modules)
  corMat  <- t(as.matrix(r))
  rownames(corMat) <- "AgeScore"
  
  textMat <- matrix(paste0(round(r, 2), "\n(", signif(p, 2), ")"),
                    nrow = 1, dimnames = list("AgeScore", colnames(MEs)))
  
  cluster_rows_flag <- nrow(corMat) > 1 
  cluster_cols_flag <- ncol(corMat) > 1
  
  pdf(outfile)
  pheatmap(corMat,
           color = colorRampPalette(c("blue","white","red"))(100),
           display_numbers = textMat,
           number_color = "black",
           cluster_rows = cluster_rows_flag,
           cluster_cols = cluster_cols_flag,
           main = title)
  dev.off()
}

empty_consensus_df <- function(){
  data.frame(Module=character(),
             beta_H1=numeric(), se_H1=numeric(), p_H1=numeric(),
             beta_H2=numeric(), se_H2=numeric(), p_H2=numeric(),
             meta_beta=numeric(), meta_se=numeric(),
             meta_z=numeric(), meta_p=numeric(),
             sign_concordant=logical(), meta_FDR=numeric())
}

empty_single_df <- function(){
  data.frame(Module=character(), beta=numeric(), se=numeric(), p=numeric(),
             FDR=numeric(), dataset=character())
}

prep_datExpr <- function(datExpr, keep_top_mad = 10000, min_genes = 50, label=""){
  v <- apply(datExpr, 2, stats::var, na.rm = TRUE)
  keep <- which(v > 0 & !is.na(v))

  datExpr2 <- datExpr[, keep, drop = FALSE]
  if (!is.null(keep_top_mad) && ncol(datExpr2) > keep_top_mad) {
    madv <- apply(datExpr2, 2, mad, na.rm = TRUE)
    ord <- order(madv, decreasing = TRUE)
    datExpr2 <- datExpr2[, ord[seq_len(keep_top_mad)], drop = FALSE]
  }
  if (ncol(datExpr2) < min_genes) {
    message("[", label, "] Too few genes after filtering (", ncol(datExpr2),
            " < ", min_genes, "); skipping clustering.")
    return(NULL)
  }
  return(datExpr2)
}

choose_minModuleSize <- function(nGenes){
  max(20, round(0.025 * nGenes))  # ~2.5% of genes, min 20
}


################################################################################
# 4) Main loop: consensus for shared, single for unique
################################################################################

enriched_GO_list <- list()
all_module_tables <- list()
all_ME_age_results <- list()

for (cell in all_cell_types) {
  
  outdir <- file.path(getwd(), cell)
  
  in1 <- cell %in% cts1
  in2 <- cell %in% cts2
  
  if (in1 && in2) {
    ############################################################################
    # CONSENSUS WGCNA
    ############################################################################
    mode_str <- "consensus"
    
    # Prepare overlapping genes and expression (samples x genes)
    X1_raw <- as.data.frame(vst.counts.GSE202601[[cell]])
    X2_raw <- as.data.frame(vst.counts.GSE255690[[cell]])
    overlapping.genes <- intersect(rownames(X1_raw), rownames(X2_raw))
    
    X1 <- X1_raw[overlapping.genes, , drop=FALSE]
    X2 <- X2_raw[overlapping.genes, , drop=FALSE]
    
    datExpr.H1 <- safe_goodSamplesGenes(as.data.frame(t(X1)))
    datExpr.H2 <- safe_goodSamplesGenes(as.data.frame(t(X2)))
    
    # Pre-filter per dataset (variance/MAD/min genes)
    datExpr.H1 <- prep_datExpr(datExpr.H1, keep_top_mad = 10000, min_genes = 50, label = paste(cell, "H1"))
    datExpr.H2 <- prep_datExpr(datExpr.H2, keep_top_mad = 10000, min_genes = 50, label = paste(cell, "H2"))
       
    traits1 <- align_traits(traits.GSE202601, rownames(datExpr.H1), "H1/GSE202601")
    traits2 <- align_traits(traits.GSE255690, rownames(datExpr.H2), "H2/GSE255690")
    
    power <- pick_power_multiple_safe(datExpr.H1, datExpr.H2, powers=1:30)
    minModuleSize <- choose_minModuleSize(ncol(datExpr.H1))
    
    cons <- blockwiseConsensusModules(
      multiExpr = list(H1=list(data=datExpr.H1), H2=list(data=datExpr.H2)),
      power = power,
      networkType = "signed",
      corType = "bicor",
      maxBlockSize = 40000,
      TOMType = "signed",
      minModuleSize = minModuleSize,
      deepSplit = 2,
      mergeCutHeight = 0.20,
      pamRespectsDendro = FALSE,
      verbose = 5
    )
    
    moduleColors <- cons$colors
    MElist       <- cons$multiMEs
    
    save(moduleColors, MElist,
         file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_Consensus_Modules.RData")))
    
    # Dendrogram
    pdf(file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_Gene_Dendrogram_with_module_colors.pdf")),
        width = 10, height = 6)
    suppressWarnings({
      adj1 <- adjacency(datExpr.H1, power = power, type="signed")
      TOM1 <- TOMsimilarity(adj1, TOMType="signed")
      geneTree <- hclust(as.dist(1 - TOM1), method = "average")
      plotDendroAndColors(geneTree, moduleColors,
                          groupLabels = "Consensus Modules",
                          main = paste("Gene Clustering -", cell),
                          dendroLabels = FALSE)
    })
    dev.off()
    
    # Eigengenes per dataset
    MEs.H1 <- tryCatch(MElist$H1$data, error=function(e) NULL)
    MEs.H2 <- tryCatch(MElist$H2$data, error=function(e) NULL)
    
    if (!is.null(MEs.H1) && ncol(MEs.H1) > 0 && !all(grepl("^ME", colnames(MEs.H1)))) {
      colnames(MEs.H1) <- paste0("ME", colnames(MEs.H1))
    }
    if (!is.null(MEs.H2) && ncol(MEs.H2) > 0 && !all(grepl("^ME", colnames(MEs.H2)))) {
      colnames(MEs.H2) <- paste0("ME", colnames(MEs.H2))
    }
    
    MEs.H1 <- MEs.H1[rownames(datExpr.H1), , drop=FALSE]
    MEs.H2 <- MEs.H2[rownames(datExpr.H2), , drop=FALSE]
    
    plot_ME_age_heatmap(MEs.H1, traits1$AgeScore,
                        file.path(outdir, paste0(Sys.Date(), "_", cell, "_H1_ME_AgeScore_Heatmap.pdf")),
                        paste(cell, "- H1 (GSE202601): ME ~ AgeScore"))
    plot_ME_age_heatmap(MEs.H2, traits2$AgeScore,
                        file.path(outdir, paste0(Sys.Date(), "_", cell, "_H2_ME_AgeScore_Heatmap.pdf")),
                        paste(cell, "- H2 (GSE255690): ME ~ AgeScore"))
    
    mods <- colnames(MEs.H1)
    if (length(mods) == 0) {
      res_df <- empty_consensus_df()
    } else {
      res_list <- lapply(mods, function(m){
        cov1 <- data.frame(row.names = rownames(MEs.H1))
        cov2 <- data.frame(row.names = rownames(MEs.H2))
        r1 <- fit_ME_age_safe(MEs.H1[,m], traits1$AgeScore, cov1, tag=paste(cell,"H1",m))
        r2 <- fit_ME_age_safe(MEs.H2[,m], traits2$AgeScore, cov2, tag=paste(cell,"H2",m))
        meta <- meta_beta(beta=c(r1["beta"], r2["beta"]),
                          se  =c(r1["se"],   r2["se"]))
        data.frame(Module=m,
                   beta_H1=r1["beta"], se_H1=r1["se"], p_H1=r1["p"],
                   beta_H2=r2["beta"], se_H2=r2["se"], p_H2=r2["p"],
                   meta_beta=meta$meta_beta, meta_se=meta$meta_se,
                   meta_z=meta$meta_z, meta_p=meta$meta_p,
                   sign_concordant = (sign(as.numeric(r1["beta"])) == sign(as.numeric(r2["beta"]))),
                   stringsAsFactors = FALSE)
      })
      res_df <- bind_rows(res_list)
      if (nrow(res_df) == 0) res_df <- empty_consensus_df()
      res_df$meta_FDR <- p.adjust(res_df$meta_p, method="fdr")
    }
    
    write.csv(res_df, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_Consensus_ME_Age_Associations.csv")),
              row.names = FALSE)
    all_ME_age_results[[cell]] <- res_df
    
    # Significant modules
    sig_tab <- tryCatch(res_df %>% dplyr::filter(meta_FDR < 0.1 & sign_concordant), error=function(e) res_df[0,])
    significant_modules <- if (nrow(sig_tab)>0) sig_tab$Module else character(0)
    
    # Gene lists + GO
    if (length(significant_modules) > 0) {
      module_gene_list <- list()
      genes_all <- colnames(datExpr.H1)
      names(moduleColors) <- genes_all
      
      for (m in significant_modules) {
        mcol <- sub("^ME", "", m)
        module_genes <- names(moduleColors)[moduleColors == mcol]
        module_gene_list[[m]] <- module_genes
        
        enriched_GO <- tryCatch({
          enrichGO(gene          = module_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   universe      = overlapping.genes,
                   ont           = "ALL",
                   pvalueCutoff  = 0.1)
        }, error=function(e) NULL)
        
        if(!is.list(enriched_GO_list[[cell]])) enriched_GO_list[[cell]] <- list()
        enriched_GO_list[[cell]][[paste0(m, "_CONSENSUS")]] <- enriched_GO
      }
      
      save(module_gene_list, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_Consensus_Significant_Module_Gene_Lists.RData")))
      write.csv(significant_modules, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_Consensus_Significant_Modules.csv")),
                row.names = FALSE)
    } 
    
    keep_cols <- c("Module","beta_H1","p_H1","beta_H2","p_H2","meta_beta","meta_p","meta_FDR","sign_concordant")
    if (nrow(res_df) > 0) {
      all_module_tables[[cell]] <- res_df[, intersect(keep_cols, colnames(res_df)), drop=FALSE]
    } else {
      all_module_tables[[cell]] <- res_df
    }
    
  } else {
    ############################################################################
    # SINGLE-DATASET WGCNA
    ############################################################################
    src_tag <- if (in1) "H1_GSE202601" else "H2_GSE255690"
    mode_str <- paste0("single_", src_tag)
    
    if (in1) {
      X_raw <- as.data.frame(vst.counts.GSE202601[[cell]])
      traits_src_master <- traits.GSE202601
    } else {
      X_raw <- as.data.frame(vst.counts.GSE255690[[cell]])
      traits_src_master <- traits.GSE255690
    }
    
    datExpr <- safe_goodSamplesGenes(as.data.frame(t(X_raw)))
    datExpr <- prep_datExpr(datExpr, keep_top_mad = 10000, min_genes = 50, label = paste(cell, src_tag))
        
    save(datExpr, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_cleaned_WGCNA_input.RData")))
    
    traits_src <- align_traits(traits_src_master, rownames(datExpr), src_tag)
    
    softPower <- pick_power_single_safe(datExpr, powers=1:30)
    minModuleSize <- choose_minModuleSize(ncol(datExpr))
    
    adjacency_mat <- adjacency(datExpr, power=softPower, type="signed")
    TOM <- TOMsimilarity(adjacency_mat, TOMType="signed")
    dissTOM <- 1 - TOM
    save(TOM, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_TOM_Matrix.RData")))
    
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, minClusterSize = minModuleSize)
    moduleColors <- labels2colors(dynamicMods)
    save(moduleColors, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_Module_Colors.RData")))
    
    pdf(file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_Gene_Dendrogram.pdf")), width = 10, height = 6)
    plot(geneTree, main = paste("Gene Clustering Dendrogram -", cell, "-", src_tag), xlab = "", sub = "", labels = FALSE)
    dev.off()
    
    pdf(file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_Gene_Dendrogram_with_module_colors.pdf")), width = 10, height = 6)
    plotDendroAndColors(geneTree, moduleColors, groupLabels = "Modules",
                        main = paste("Gene Clustering Dendrogram -", cell, "-", src_tag),
                        cex.axis = 0.5, cex.lab = 0.7, dendroLabels = FALSE)
    dev.off()
    
    MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
    if (!is.null(MEs) && ncol(MEs) > 0 && !all(grepl("^ME", colnames(MEs)))) {
      colnames(MEs) <- paste0("ME", colnames(MEs))
    }
    
    plot_ME_age_heatmap(MEs, traits_src$AgeScore,
                        file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_ME_AgeScore_Heatmap.pdf")),
                        paste(cell, "-", src_tag, ": ME ~ AgeScore"))
    
    mods <- colnames(MEs)
    if (length(mods) == 0) {
      res_df <- empty_single_df()
    } else {
      cov_src <- data.frame(row.names = rownames(MEs))
      res_list <- lapply(mods, function(m){
        r <- fit_ME_age_safe(MEs[,m], traits_src$AgeScore, cov_src, tag=paste(cell, src_tag, m))
        data.frame(Module=m, beta=r["beta"], se=r["se"], p=r["p"], stringsAsFactors = FALSE)
      })
      res_df <- bind_rows(res_list)
      if (nrow(res_df) == 0) res_df <- empty_single_df()
      res_df$FDR <- p.adjust(res_df$p, method="fdr")
      res_df$dataset <- src_tag
    }
    
    write.csv(res_df, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_ME_Age_Associations.csv")),
              row.names = FALSE)
    all_ME_age_results[[cell]] <- res_df
    
    significant_modules <- if (nrow(res_df)>0) {
      res_df %>% dplyr::filter(FDR < 0.1) %>% dplyr::pull(Module)
    } else character(0)
    
    if (length(significant_modules) > 0) {
      module_gene_list <- list()
      genes_all <- colnames(datExpr)
      names(moduleColors) <- genes_all
      
      for (m in significant_modules) {
        mcol <- sub("^ME", "", m)
        module_genes <- names(moduleColors)[moduleColors == mcol]
        module_gene_list[[m]] <- module_genes
        
        enriched_GO <- tryCatch({
          enrichGO(gene          = module_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   universe      = genes_all,
                   ont           = "ALL",
                   pvalueCutoff  = 0.1)
        }, error=function(e) NULL)
        
        if(!is.list(enriched_GO_list[[cell]])) enriched_GO_list[[cell]] <- list()
        enriched_GO_list[[cell]][[paste0(m, "_", src_tag)]] <- enriched_GO
      }
      
      save(module_gene_list, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_Significant_Module_Gene_Lists.RData")))
      write.csv(significant_modules, file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_Significant_Modules.csv")),
                row.names = FALSE)
    } else {
      write.csv(data.frame(Module=character(0)),
                file = file.path(outdir, paste0(Sys.Date(), "_", cell, "_", src_tag, "_Significant_Modules.csv")),
                row.names = FALSE)
    }
  }
  
  writeLines(mode_str, con = file.path(outdir, paste0(Sys.Date(), "_", cell, "_mode.txt")))
}

################################################################################
# Aggregate aging-relevant module gene sets into ONE RData
################################################################################

# Define functions

write_gmt <- function(gene_sets, file){
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(gene_sets)) {
    # description field set to 'NA' (second column), adjust as needed
    cat(nm, "NA", paste(gene_sets[[nm]], collapse = "\t"), sep = "\t", file = con)
    cat("\n", file = con)
  }
}

# Discover cell-type folders 
cell_dirs <- list.dirs(getwd(), full.names = TRUE, recursive = FALSE)

aging_gene_sets <- list()  
set_metadata    <- list()  

for (cell_dir in cell_dirs) {
  cell <- basename(cell_dir)
  
  # --------- CONSENSUS (shared cell types) ---------
  cons_csv  <- list.files(cell_dir, pattern = "_Consensus_Significant_Modules\\.csv$", full.names = TRUE)
  cons_rdat <- list.files(cell_dir, pattern = "_Consensus_Significant_Module_Gene_Lists\\.RData$", full.names = TRUE)
  
  if (length(cons_csv) > 0 && length(cons_rdat) > 0) {
    sig_mods <- tryCatch({
      x <- read.csv(cons_csv[1], stringsAsFactors = FALSE)
      if ("ME" %in% colnames(x)) x$Module else character(0)
    }, error = function(e) character(0))
    
    if (length(sig_mods) > 0) {
      # Load module_gene_list object
      env_cons <- new.env()
      ok <- tryCatch({ load(cons_rdat[1], envir = env_cons); TRUE }, error = function(e) FALSE)
      if (ok && exists("module_gene_list", envir = env_cons)) {
        mgl <- get("module_gene_list", envir = env_cons)
        for (mod in intersect(names(mgl), sig_mods)) {
          genes <- unique(na.omit(as.character(mgl[[mod]])))
          if (length(genes) == 0) next
          set_name <- paste("Human", "CONSENSUS", cell, mod, sep = "__")
          aging_gene_sets[[set_name]] <- genes
          set_metadata[[set_name]] <- list(cell = cell, tag = "CONSENSUS", module = mod, n_genes = length(genes))
        }
      }
    }
  }
  
  # --------- SINGLE-DATASET (unique cell types) ---------
  # GSE202601
  s_csv_h1  <- list.files(cell_dir, pattern = "_H1_GSE202601_ME_Age_Associations\\.csv$", full.names = TRUE)
  s_rdat_h1 <- list.files(cell_dir, pattern = "_H1_GSE202601_Significant_Module_Gene_Lists\\.RData$", full.names = TRUE)
  if (length(s_csv_h1) > 0 && length(s_rdat_h1) > 0) {
    sig_mods <- tryCatch({
      x <- read.csv(s_csv_h1[1], stringsAsFactors = FALSE)
      if ("Module" %in% colnames(x) && "FDR" %in% colnames(x)) x$Module[x$FDR < 0.1] else character(0)
    }, error = function(e) character(0))
    
    env_h1 <- new.env()
    ok <- tryCatch({ load(s_rdat_h1[1], envir = env_h1); TRUE }, error = function(e) FALSE)
    if (ok && exists("module_gene_list", envir = env_h1)) {
      mgl <- get("module_gene_list", envir = env_h1)
      for (mod in intersect(names(mgl), sig_mods)) {
        genes <- unique(na.omit(as.character(mgl[[mod]])))
        if (length(genes) == 0) next
        set_name <- paste("Human", "H1_GSE202601", cell, mod, sep = "__")
        aging_gene_sets[[set_name]] <- genes
        set_metadata[[set_name]] <- list(cell = cell, tag = "H1_GSE202601", module = mod, n_genes = length(genes))
      }
    }
  }
  
  # GSE255690
  s_csv_h2  <- list.files(cell_dir, pattern = "_H2_GSE255690_ME_Age_Associations\\.csv$", full.names = TRUE)
  s_rdat_h2 <- list.files(cell_dir, pattern = "_H2_GSE255690_Significant_Module_Gene_Lists\\.RData$", full.names = TRUE)
  if (length(s_csv_h2) > 0 && length(s_rdat_h2) > 0) {
    sig_mods <- tryCatch({
      x <- read.csv(s_csv_h2[1], stringsAsFactors = FALSE)
      if ("Module" %in% colnames(x) && "FDR" %in% colnames(x)) x$Module[x$FDR < 0.1] else character(0)
    }, error = function(e) character(0))
    
    env_h2 <- new.env()
    ok <- tryCatch({ load(s_rdat_h2[1], envir = env_h2); TRUE }, error = function(e) FALSE)
    if (ok && exists("module_gene_list", envir = env_h2)) {
      mgl <- get("module_gene_list", envir = env_h2)
      for (mod in intersect(names(mgl), sig_mods)) {
        genes <- unique(na.omit(as.character(mgl[[mod]])))
        if (length(genes) == 0) next
        set_name <- paste("Human", "H2_GSE255690", cell, mod, sep = "__")
        aging_gene_sets[[set_name]] <- genes
        set_metadata[[set_name]] <- list(cell = cell, tag = "H2_GSE255690", module = mod, n_genes = length(genes))
      }
    }
  }
}

# Save all data
out_rdata <- file.path(getwd(), ))
save(aging_gene_sets, set_metadata, file = paste0(Sys.Date(), "_Human_AgingRelevantModules_geneSets.RData")

# Save as GMT for GSEA
out_gmt <- file.path(getwd(), paste0(Sys.Date(), "_Human_AgingRelevantModules_geneSets.gmt"))
if (length(aging_gene_sets) > 0) {
  write_gmt(aging_gene_sets, out_gmt)
}

################################################################################
# 5) Generate ORA plot
################################################################################

all_go_terms_list <- list()
for (cell in names(enriched_GO_list)) {
  for (module in names(enriched_GO_list[[cell]])) {
    enriched_result <- enriched_GO_list[[cell]][[module]]
    if (!is.null(enriched_result) && nrow(enriched_result@result) > 0) {
      go_data <- enriched_result@result %>%
        dplyr::select(ID, Description, p.adjust, Count) %>%
        mutate(CellType = cell, Module = module)
      all_go_terms_list[[paste0(cell, "_", module)]] <- go_data
    }
  }
}

if (length(all_go_terms_list) > 0) {
  all_go_terms_df <- dplyr::bind_rows(all_go_terms_list)
  all_go_terms_df$logP <- -log10(all_go_terms_df$p.adjust)
  
  pdf(paste0(Sys.Date(), "_Human_WGCNA_ORA_bubbleplots_consensus_and_single.pdf"), width = 20, height = 30)
  print(
    ggplot(all_go_terms_df, aes(x = Module, y = Description, size = Count, color = p.adjust)) +
      geom_point(alpha = 0.8) +
      scale_size_continuous(range = c(3, 10)) +
      scale_color_gradientn(colors = c("darkorange2", "purple"), limits = c(0, 0.1)) +
      facet_wrap(~ CellType, scales = "free_x") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(size = "Gene Count", color = "p.adjust", x = "Module", y = "GO Term",
           title = "Human WGCNA: Consensus (both) + Single-dataset (unique) Modules")
  )
  dev.off()
}

################################################################################
sink(file = paste0(Sys.Date(), "_Human_PB_WGCNA_session_info.txt"))
sessionInfo()
sink()
