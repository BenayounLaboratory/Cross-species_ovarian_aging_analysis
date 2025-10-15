options(stringsAsFactors = FALSE)

library(DESeq2)
library(clusterProfiler)
library(fgsea)    
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)

theme_set(theme_bw())

################################################################################
# Integrative ovarian aging analysis
# Run GSEA - age significant WGCNA module gene sets
################################################################################

################################################################################
# 1. Load data
################################################################################

# Gene sets
gene_sets <- fgsea::gmtPathways("/Volumes/jinho01/Benayoun_lab/Projects/CZI/5_WGCNA/Human_gene_list_objects/2025-09-23_Human_AgingRelevantModules_geneSets.gmt")  # WGCNA modules as GMT gene sets

# DESeq2 results
load(/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/0_Humanized_DESeq2_objects/2025-08-18_all_combined_cell_type_DESeq2_results.RData) 

################################################################################
# 2. Define parameters
################################################################################

fdr_cutoff     <- 0.10     # BH FDR threshold
minGSSize      <- 5        # pathway size lower bound
maxGSSize      <- 10000    # pathway size upper bound
top_per_ds     <- 10       # terms per dataset in combined plot
balance_posneg <- TRUE     # pick half pos / half neg NES in per-dataset top

# Desired dataset order for plotting
DATASET_ORDER <- c(
  "Human_GSE202601", "Human_GSE255690",
  "Monkey_STRTseq_GSE130664",
  "Mouse_Aging", "Mouse_Foxl2", "Mouse_VCD", "Mouse_GSE232309", "Mouse_EMTAB12889",
  "Goat"
)

###############################################################################
# 3. Define functions
###############################################################################

# Extract cell names
get_cell_from_setname <- function(s) {
  p <- strsplit(s, "__", fixed = TRUE)[[1]]
  if (length(p) >= 3) p[3] else NA_character_
}

# Build TERM2GENE for a cell type
term2gene_for_cell <- function(gene_sets, cell) {
  keep_sets <- names(gene_sets)[vapply(names(gene_sets), function(nm) get_cell_from_setname(nm) == cell, logical(1))]
  if (!length(keep_sets)) return(data.frame(gs_name=character(0), gene_symbol=character(0)))
  do.call(rbind, lapply(keep_sets, function(nm) {
    gs <- gene_sets[[nm]]
    data.frame(gs_name = nm, gene_symbol = as.character(gs), stringsAsFactors = FALSE)
  }))
}

# Build ranked gene list from a DESeq2-like table using 'stat'
build_geneList_from_stat <- function(res_df) {
  df <- as.data.frame(res_df)
  if (is.null(rownames(df))) {
    gc <- intersect(c("gene","symbol","Gene","SYMBOL","GeneSymbol"), colnames(df))
    if (length(gc) >= 1) rownames(df) <- as.character(df[[gc[1]]])
  }
  gl <- df$stat
  names(gl) <- rownames(df)
  gl <- gl[!is.na(gl) & is.finite(gl)]
  sort(gl, decreasing = TRUE)
}

# GSEA function
perform_gsea <- function(geneList, term2gene,
                         fdr_filter = 0.1,
                         minGSSize = 5, maxGSSize = 10000,
                         verbose = TRUE) {
  
  gsea_obj <- clusterProfiler::GSEA(
        geneList     = geneList,
        TERM2GENE    = term2gene,
        minGSSize    = minGSSize,
        maxGSSize    = maxGSSize,
        pvalueCutoff = 1,
        verbose      = FALSE
  )
  
  return (gsea_obj)
}

is_human_ds <- function(dsname) grepl("^Human_", dsname)

run_wgcna_gsea_all <- function(combined_cell_type_objects, gene_sets,
                               fdr_cutoff = 0.1,
                               include_human = TRUE,
                               minGSSize = 5, maxGSSize = 10000) {
  results <- list()
  cells_common <- intersect(
    names(combined_cell_type_objects),
    unique(vapply(names(gene_sets), get_cell_from_setname, "", USE.NAMES = FALSE))
  )
  
  for (cell in cells_common) {
    results[[cell]] <- list()
    
    t2g <- term2gene_for_cell(gene_sets, cell)
    
    ds_names <- names(combined_cell_type_objects[[cell]])
    if (!include_human) ds_names <- ds_names[!is_human_ds(ds_names)]
    
    for (ds in ds_names) {
      
      res_df <- combined_cell_type_objects[[cell]][[ds]]
      
      geneList <- build_geneList_from_stat(res_df)
      
      gsea_res <- perform_gsea(geneList, t2g,
                     fdr_filter = fdr_cutoff,
                     minGSSize = minGSSize, maxGSSize = maxGSSize,
                     verbose = TRUE)
      
      tab <- gsea_res@result
      
      out_csv <- paste0(Sys.Date(), "_WGCNA_GSEA_", cell, "_", ds, ".csv")
      write.csv(tab, out_csv, row.names = FALSE)
      
      generate_bubble_plot(cell, tab, ds)
      results[[cell]][[ds]] <- tab
    }
  }
  results
}

# For plotting

NES_LIMITS <- c(-3, 3)
NES_VALUES <- scales::rescale(c(-3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3))
NES_COLORS <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1",
                "white","lightcoral","brown1","firebrick2","firebrick4")

generate_bubble_plot <- function(my.cell.type, my.gsea.data, ds_name, max.path.plot = 20) {

  my.gsea.data <- my.gsea.data[is.finite(my.gsea.data$NES), , drop = FALSE]
  
  # split by direction
  pos <- my.gsea.data[my.gsea.data$NES > 0, ]
  neg <- my.gsea.data[my.gsea.data$NES < 0, ]
  
  pos_ix <- sort(pos$NES, index.return = TRUE, decreasing = TRUE)$ix
  neg_ix <- sort(neg$NES, index.return = TRUE, decreasing = FALSE)$ix
  
  if ((nrow(pos) > round(max.path.plot/2)) && (nrow(neg) > round(max.path.plot/2))) {
    sel <- rbind(pos[pos_ix[1:round(max.path.plot/2)], ],
                 neg[neg_ix[1:round(max.path.plot/2)], ])
  } else {
    sel <- rbind(pos[pos_ix[1:min(round(max.path.plot/2), nrow(pos))], ],
                 neg[neg_ix[1:min(round(max.path.plot/2), nrow(neg))], ])
    sel <- sel[!is.na(sel$ID), ]
  }
  
  sel$minlog10fdr <- -log10(sel$p.adjust)
  sel$PathName    <- paste0(sel$ID, " ", sel$Description)
  
  ord <- order(sel$minlog10fdr, decreasing = TRUE)
  plotdf <- sel[ord, ]
  plotdf$age <- ifelse(plotdf$NES < 0, "YF", "OF")
  plotdf     <- plotdf[order(plotdf$age), ]
  plotdf$PathName  <- factor(plotdf$PathName, levels = rev(unique(plotdf$PathName)))
  plotdf$CellType  <- factor(rep("PB", nrow(plotdf)))
  
  p <- ggplot(plotdf, aes(x = CellType, y = PathName, colour = NES, size = minlog10fdr)) +
    geom_point(shape = 16, alpha = 0.95) +  
    ggtitle(paste("GSEA (WGCNA modules) for", my.cell.type, "-", ds_name)) +
    labs(x = "-log10(FDR, BH)", y = "") +
    scale_colour_gradientn(
      colours = NES_COLORS,
      values  = NES_VALUES,
      limits  = NES_LIMITS,
      oob     = scales::squish,
      na.value = "grey50"
    ) +
    theme(axis.text.y = element_text(size = 6),
          plot.margin = margin(10, 10, 10, 50))
  
  pdf(paste0(Sys.Date(), "_WGCNA_GSEA_BubblePlot_", my.cell.type, "_", ds_name,
             "_top_", nrow(plotdf), "_terms.pdf"),
      onefile = TRUE, height = 6, width = 30)
  print(p)
  dev.off()
}

###############################################################################
# 4. Run GSEA
###############################################################################

wgcna_gsea_results <- run_wgcna_gsea_all(
  combined_cell_type_objects, gene_sets,
  fdr_cutoff = fdr_cutoff,
  include_human = TRUE,
  minGSSize = minGSSize, maxGSSize = maxGSSize
)

save(wgcna_gsea_results, file = paste0(Sys.Date(), "_WGCNA_GSEA_results_list.RData"))

###############################################################################
sink(file = paste0(Sys.Date(), "_WGCNA_module_GSEA_session_info.txt"))
sessionInfo()
sink()
