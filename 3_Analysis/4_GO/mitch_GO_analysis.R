setwd("/Volumes/jinho01/Benayoun_lab/Projects/CZI/2_PB_GO/mitch/2025-09-24")
options(stringsAsFactors = FALSE)

# Load necessary libraries
library(mitch)
library(DESeq2)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(readr)
library(stringr)

theme_set(theme_bw())

rm(list = ls())

################################################################################
# Load data
################################################################################

# Load DESeq2 results
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/0_Humanized_DESeq2_objects/2025-08-18_all_combined_cell_type_DESeq2_results.RData")

# Create species mapping based on the actual dataset names
species_mapping <- data.frame(
  Dataset = c("Human_GSE202601", "Human_GSE255690", "Monkey", "Mouse_Aging", "Mouse_VCD", "Mouse_Foxl2", "Mouse_EMTAB12889", "Mouse_GSE232309", "Goat"),
  Species = c("Human", "Human", "Monkey", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Goat"),
  stringsAsFactors = FALSE
)

# Get all cell types
all_cell_types <- names(combined_cell_type_objects)
all_cell_types

# Get all datasets
first_cell_type <- all_cell_types[1]
names(combined_cell_type_objects[[first_cell_type]])

################################################################################
# Load gene sets for mitch analysis
################################################################################

load_custom_genesets <- function(path,
                                 min_set_size = 5,
                                 max_set_size = 2000,
                                 join_id_and_name = TRUE,
                                 verbose = TRUE) {
  
  # read as character, keep everything
  df <- readr::read_tsv(
    path,
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  
  getcol <- function(d, candidates) {
    for (nm in candidates) if (nm %in% names(d)) return(d[[nm]])
  }
  
  gene_col <- getcol(df, c("Gene name","Gene symbol","gene_symbol","Gene", "HGNC symbol"))
  go_id    <- getcol(df, c("GO term accession","GO.ID","GO ID","go_id"))
  go_name  <- getcol(df, c("GO term name","Term","GO term","go_name"))
  
  # trim whitespace & normalize
  gene_col <- str_trim(gene_col)
  go_id    <- str_trim(go_id)
  go_name  <- str_squish(go_name)
  
  # drop empty rows
  keep <- !(is.na(gene_col) | gene_col == "" | is.na(go_id) | go_id == "")
  df2 <- tibble(
    gene = gene_col[keep],
    go_id = go_id[keep],
    go_name = go_name[keep]
  ) %>%
    distinct()  # remove duplicate rows
  
  # build set name
  set_name <- if (join_id_and_name) {
    paste0(df2$go_id, " | ", df2$go_name)
  } else {
    df2$go_id
  }
  
  df2$set <- set_name
  
  # split to list: one element per GO term set
  # ensure unique, sorted genes per set
  sets <- split(df2$gene, df2$set)
  sets <- lapply(sets, function(v) sort(unique(v)))
  
  # filter by size
  sizes <- vapply(sets, length, integer(1))
  keep_sets <- which(sizes >= min_set_size & sizes <= max_set_size)
  sets_filt <- sets[keep_sets]
  
  return(sets_filt)
}

write_gmt <- function(genesets, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(genesets)) {
    line <- c(nm, "na", genesets[[nm]])
    writeLines(paste(line, collapse = "\t"), con = con)
  }
}

# Load custom GO gene sets
custom_genesets <- load_custom_genesets(
  "/Volumes/jinho01/Benayoun_lab/Projects/CZI/2_PB_GO/mitch/mart_export.txt",
  min_set_size = 5,   
  max_set_size = 5000
)

write_gmt(custom_genesets, file = paste0(Sys.Date(), "_ensembl_custom_GO_sets.gmt"))

################################################################################
# Define functions for mitch
################################################################################

# Function to prepare DESeq2 results for mitch
prepare_deseq_for_mitch <- function(deseq_obj) {
  if (inherits(deseq_obj, "DESeqResults")) {
    # Extract required columns for mitch
    result <- data.frame(
      baseMean = deseq_obj$baseMean,
      log2FoldChange = deseq_obj$log2FoldChange,
      lfcSE = deseq_obj$lfcSE,
      stat = deseq_obj$stat,
      pvalue = deseq_obj$pvalue,
      padj = deseq_obj$padj,
      row.names = rownames(deseq_obj)
    )
  } else if (is.data.frame(deseq_obj) && "stat" %in% colnames(deseq_obj)) {
    result <- deseq_obj
  } 
  
  # Remove rows with NA values
  result <- result[complete.cases(result), ]
  
  return(result)
}

# Function to add small random noise to break identical contrasts
add_minimal_noise <- function(mitch_input_list) {
  lapply(mitch_input_list, function(dataset) {
    # Add tiny random noise to log2FoldChange to break perfect correlations
    noise_factor <- 1e-10
    dataset$log2FoldChange <- dataset$log2FoldChange + rnorm(nrow(dataset), 0, noise_factor)
    return(dataset)
  })
}

# Function to prioritize human datasets in duplicate removal
prioritize_human_datasets <- function(mitch_input_list, species_mapping) {
  # If force mode is enabled, skip human dataset prioritization
  if (FORCE_KEEP_ALL_DATASETS) {
    message(" - FORCE MODE: Skipping human dataset prioritization")
    return(mitch_input_list)
  }
  
  # Identify human datasets
  human_datasets <- names(mitch_input_list)[grepl("Human", names(mitch_input_list))]
  
  if (length(human_datasets) > 1) {
    message(" - Found ", length(human_datasets), " human datasets, prioritizing all of them")
    
    # For human datasets, use even stricter correlation threshold (0.99)
    # This ensures we only remove truly identical datasets
    common_genes <- Reduce(intersect, lapply(mitch_input_list[human_datasets], rownames))
    
    if (length(common_genes) >= 100) {
      lfc_matrix <- do.call(cbind, lapply(mitch_input_list[human_datasets], function(x) {
        x[common_genes, "log2FoldChange"]
      }))
      colnames(lfc_matrix) <- human_datasets
      
      cor_matrix_human <- cor(lfc_matrix, use = "pairwise.complete.obs")
      
      # Only remove if correlation is > 0.99 (essentially identical)
      high_cor_human <- which(cor_matrix_human > 0.99 & upper.tri(cor_matrix_human), arr.ind = TRUE)
      
      if (nrow(high_cor_human) > 0) {
        message(" - Found ", nrow(high_cor_human), " essentially identical human dataset pairs (cor > 0.99)")
        datasets_to_remove <- unique(high_cor_human[, 2])
        datasets_to_keep <- setdiff(human_datasets, human_datasets[datasets_to_keep])
        
        message(" - Keeping human datasets: ", paste(datasets_to_keep, collapse = ", "))
        return(mitch_input_list[datasets_to_keep])
      }
    }
  }
  
  return(mitch_input_list)
}

################################################################################
# Run mitch
################################################################################

# Create mitch input data structure
mitch_data <- list()
mitch_results <- list()

# Process each cell type
for (cell_type in all_cell_types) {
  
  # Get datasets for this cell type
  cell_datasets <- combined_cell_type_objects[[cell_type]]
  
  if (length(cell_datasets) > 0) {
    # Prepare data for mitch
    mitch_input_list <- list()
    
    for (dataset_name in names(cell_datasets)) {
      dataset <- cell_datasets[[dataset_name]]
      
      prepared_data <- prepare_deseq_for_mitch(dataset)
      mitch_input_list[[dataset_name]] <- prepared_data
      
    }
    
    # Only proceed if we have at least 2 datasets
    if (length(mitch_input_list) >= 2) {
      
      mitch_input_list <- add_minimal_noise(mitch_input_list)
      mitch_data[[cell_type]] <- mitch_import(mitch_input_list, DEtype = "DESeq2")
      
      mitch_results[[cell_type]] <- mitch_calc(
            mitch_data[[cell_type]], 
            genesets = custom_genesets,
            priority = "significance", 
            minsetsize = 10, 
            resrows = 50, 
            cores = 1
      )
    }     
  }
}

# Save all mitch data and results
save(mitch_data, mitch_results, 
     file = paste0(Sys.Date(), "_mitch_Comparative_Analysis_Results.RData"))

################################################################################
# Create heatmaps of top pathways
################################################################################

# Function to create heatmap for top pathways
create_pathway_heatmap <- function(mitch_results, cell_type, analysis_type, top_n = 15) {
  if (!is.null(mitch_results) && "enrichment_result" %in% names(mitch_results)) {
    # Get top significant results
    top_results <- head(mitch_results$enrichment_result, top_n)

    if (nrow(top_results) > 0) {
      # Extract effect sizes (s values)
      effect_cols <- grep("^s\\.", colnames(top_results), value = TRUE)
      effect_cols <- effect_cols[effect_cols != "s.dist"]
      
      if (length(effect_cols) > 0) {
        # Create matrix for heatmap
        heatmap_data <- as.matrix(top_results[, effect_cols, drop = FALSE])
        rownames(heatmap_data) <- top_results$set
        
        # Create heatmap
        hm_clustered <- Heatmap(
          heatmap_data,
          name = "Effect Size",
          col = colorRamp2(c(-0.75, 0, 0.75), c("darkblue", "white", "#CC3333")),
          border = TRUE,
          rect_gp = gpar(col = "grey", lwd = 0.5),
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          column_title = paste(cell_type, "-", analysis_type, "(Clustered)"),
          width = length(effect_cols) * unit(7, "mm"),
          height = nrow(heatmap_data) * unit(7, "mm")
        )
        
        # Save heatmap
        pdf(paste0(Sys.Date(), "_", cell_type, "_", analysis_type, "_Top", top_n, "_Pathways_Heatmap_Clustered.pdf"), 
            width = 10, height = 15)
        print(hm_clustered)
        dev.off()
        
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

# Create heatmaps for each cell type
for (cell_type in names(mitch_results)) {
  if (!is.null(mitch_results[[cell_type]])) {
    success <- create_pathway_heatmap(mitch_results[[cell_type]], cell_type, "GO", top_n = 15)
    if (success) message(" - Generated GO heatmap for ", cell_type)
  }
}

################################################################################
# Session info and cleanup
################################################################################
sink(file = paste0(Sys.Date(), "_mitch_GO_analysis_Session_Info.txt"))
sessionInfo()
sink()