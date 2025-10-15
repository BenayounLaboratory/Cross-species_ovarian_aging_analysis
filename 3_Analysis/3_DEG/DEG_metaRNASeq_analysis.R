library(ggplot2)
library(metaRNASeq)
library(ggVennDiagram)
  
rm(list = ls())
  
###############################################################################
# Integrative ovarian aging analysis
# Assess overlapping DEGs (FDR < 0.1) - pairwise
# Process species with more than 2 datasets via metaRNASeq to combine pvalues
# Generate bubble plots
###############################################################################

###############################################################################
# 1. Load data
###############################################################################

# DESeq2 results
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/0_Humanized_DESeq2_objects/2025-08-18_all_combined_cell_type_DESeq2_results.RData")

###############################################################################
# 2. Define parameters and functions
###############################################################################

# Create directories
dir.create("pairwise_bubbles", showWarnings = FALSE)
dir.create("venn_plots", showWarnings = FALSE)

# Define parameters
alpha_fdr  <- 0.10     # BOTH species must pass this FDR
log10_cap  <- 10       # cap for -log10(FDR) for plotting
n_per_dir  <- 10       # up to 10 Up and 10 Down per plot

# Define dataset names - for human & mouse
human_ds_names <- c("Human_GSE202601", "Human_GSE255690")
mouse_ds_names <- c("Mouse_Aging","Mouse_VCD","Mouse_Foxl2",
                    "Mouse_EMTAB12889","Mouse_GSE232309")

# Define colors for plotting
species_colors <- c(
  Human  = "#aed8e6",
  Monkey = "#8dc63f",
  Mouse  = "#ec008c",
  Goat   = "#ffde17"
)

# Functions
# -log10(pval) calculation & cap max value
mlog <- function(p) {
  x <- -log10(p)
  x[!is.finite(x)] <- log10_cap
  pmin(pmax(x, 0), log10_cap)
}

# Perform metaRNASeq for Human datasets
human_meta_one_cell <- function(ct_list) {
  present <- intersect(human_ds_names, names(ct_list))
  
  # If cell type is detected in one dataset
  if (length(present) == 1) {
    df <- as.data.frame(ct_list[[present]])
    out <- df[, c("padj","log2FoldChange"), drop = FALSE]
    colnames(out) <- c("Human_padj","Human_FC")
    return(out)
  }
  
  # If cell type is detected in two datasets
  ds1 <- as.data.frame(ct_list[["Human_GSE202601"]])
  ds2 <- as.data.frame(ct_list[["Human_GSE255690"]])
  common <- intersect(rownames(ds1), rownames(ds2))
  ds1 <- ds1[common, , drop = FALSE]
  ds2 <- ds2[common, , drop = FALSE]
  
  rawp_list <- list(ds1$pvalue, ds2$pvalue)
  inv <- invnorm(rawp_list, nrep = c(8, 9), BHth = 1)
  
  data.frame(
    Human_padj = p.adjust(inv$rawpval, method = "BH"),
    Human_FC   = ds1$log2FoldChange,
    row.names  = common,
    check.names = FALSE
  )
}

# Perform metaRNASeq for Mouse datasets
mouse_meta_one_cell <- function(ct_list) {
  present <- intersect(mouse_ds_names, names(ct_list))
  
  present <- present[vapply(present, function(ds) {
    df <- as.data.frame(ct_list[[ds]])
    !is.null(rownames(df)) && all(c("pvalue","log2FoldChange") %in% colnames(df))
  }, logical(1))]
  
  # If cell type is detected in one dataset
  if (length(present) == 1) {
    df <- as.data.frame(ct_list[[present]])
    out <- df[, c("padj","log2FoldChange"), drop = FALSE]
    colnames(out) <- c("Mouse_padj","Mouse_FC")
    rownames(out) <- rownames(df)
    return(out)
  }
  
  # If cell type is detected in more than one dataset
  common <- Reduce(intersect, lapply(ct_list[present], rownames))
  
  rawp_list <- lapply(present, function(ds) as.data.frame(ct_list[[ds]])[common, "pvalue", drop = TRUE])
  inv <- invnorm(rawp_list, nrep = rep(4, length(present)), BHth = 1)
  
  first_fc <- as.data.frame(ct_list[[present[1]]])[common, "log2FoldChange", drop = TRUE]
  data.frame(
    Mouse_padj = p.adjust(inv$rawpval, method = "BH"),
    Mouse_FC   = first_fc,
    row.names  = common,
    check.names = FALSE
  )
}

# Process monkey and goat datasets
pull_one_ds <- function(df, padj_col = "padj", fc_col = "log2FoldChange", prefix = "X") {
  df <- as.data.frame(df)
  out <- df[, c(padj_col, fc_col), drop = FALSE]
  colnames(out) <- c(paste0(prefix, "_padj"), paste0(prefix, "_FC"))
  rownames(out) <- rownames(df)
  out
}

# Generate bubble plots
make_bubble_long <- function(human_tbl, other_tbl, other_label,
                             fdr_thr = 0.1, n_per_direction = 10) {
  common <- intersect(rownames(human_tbl), rownames(other_tbl))
  
  H <- human_tbl[common, , drop = FALSE]
  O <- other_tbl[common, , drop = FALSE]
  
  keep <- !is.na(H$Human_padj) & !is.na(O[[paste0(other_label, "_padj")]]) &
    (H$Human_padj < fdr_thr) & (O[[paste0(other_label, "_padj")]] < fdr_thr)
  
  H <- H[keep, , drop = FALSE]
  O <- O[rownames(H), , drop = FALSE]
  
  sign_lab <- ifelse(H$Human_FC > 0, "Up",
                     ifelse(H$Human_FC < 0, "Down", "Zero"))
  nz <- sign_lab != "Zero"
  H <- H[nz, , drop = FALSE]
  O <- O[rownames(H), , drop = FALSE]
  H$Sign <- sign_lab[nz]
  
  up_idx   <- which(H$Sign == "Up")
  down_idx <- which(H$Sign == "Down")
  
  ord_up   <- if (length(up_idx))   up_idx[order(H$Human_padj[up_idx],   decreasing = FALSE)] else integer(0)
  ord_down <- if (length(down_idx)) down_idx[order(H$Human_padj[down_idx], decreasing = FALSE)] else integer(0)
  
  sel_up   <- head(ord_up,   n_per_direction)
  sel_down <- head(ord_down, n_per_direction)
  
  keep_rows <- c(sel_up, sel_down)
  
  H <- H[keep_rows, , drop = FALSE]
  O <- O[rownames(H), , drop = FALSE]
  
  gene_order <- rownames(H)
  
  df <- data.frame(
    Gene        = gene_order,
    Human_padj  = H$Human_padj,
    Human_FC    = H$Human_FC,
    Other_padj  = O[[paste0(other_label, "_padj")]],
    Other_FC    = O[[paste0(other_label, "_FC")]],
    Sign        = H$Sign,
    stringsAsFactors = FALSE
  )
  
  long <- rbind(
    data.frame(Gene = df$Gene, Species = "Human",
               mlog10 = mlog(df$Human_padj), Sign = df$Sign,
               stringsAsFactors = FALSE),
    data.frame(Gene = df$Gene, Species = other_label,
               mlog10 = mlog(df$Other_padj), Sign = df$Sign,
               stringsAsFactors = FALSE)
  )
  long$signed_mlog10 <- ifelse(long$Sign == "Up", long$mlog10, -long$mlog10)
  long$Gene <- factor(long$Gene, levels = rev(gene_order))
  long
}

plot_bubbles <- function(long_df, file, title) {
  p <- ggplot(long_df, aes(x = Species, y = Gene)) +
    geom_point(aes(size = mlog10, fill = signed_mlog10), shape = 21, color = "grey30") +
    scale_size(range = c(2.2, 8), name = expression(-log[10]("FDR"))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-log10_cap, log10_cap),
                         name = "Direction\n(by Human FC)") +
    labs(title = title, x = NULL, y = "Genes (Up top, Down bottom; sorted by Human FDR)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          legend.box = "vertical",
          legend.margin = margin(2, 2, 2, 2))
  ggsave(filename = file, plot = p, width = 5.6,
         height = max(4, 0.32 * length(levels(long_df$Gene))), units = "in")
}

###############################################################################
# 3. Run analysis
###############################################################################

available_cell_types <- c("Granulosa", "Theca", "Stroma")

human_meta  <- setNames(vector("list", length(available_cell_types)), available_cell_types)
mouse_meta  <- setNames(vector("list", length(available_cell_types)), available_cell_types)
monkey_tbl  <- setNames(vector("list", length(available_cell_types)), available_cell_types)
goat_tbl    <- setNames(vector("list", length(available_cell_types)), available_cell_types)

for (ct in available_cell_types) {
  ct_list <- combined_cell_type_objects[[ct]]
  
  human_meta[[ct]] <- human_meta_one_cell(ct_list)
  mouse_meta[[ct]] <- mouse_meta_one_cell(ct_list)
  
  mname <- intersect(c("Monkey","Monkey_STRTseq_GSE130664"), names(ct_list))
  if (length(mname)) monkey_tbl[[ct]] <- pull_one_ds(ct_list[[mname[1]]], prefix = "Monkey")
  
  if ("Goat" %in% names(ct_list)) goat_tbl[[ct]] <- pull_one_ds(ct_list[["Goat"]], prefix = "Goat")
}

###############################################################################
# 4. Generate bubble plots
###############################################################################

for (ct in available_cell_types) {
  
  H <- human_meta[[ct]]
  
  # Human vs Monkey
  if (!is.null(H) && !is.null(monkey_tbl[[ct]])) {
    L <- make_bubble_long(H, monkey_tbl[[ct]], "Monkey",
                          fdr_thr = alpha_fdr, n_per_direction = n_per_dir)
    if (!is.null(L) && nrow(L)) {
      plot_bubbles(L,
                   file  = file.path("pairwise_bubbles", paste0(ct, "_Human_vs_Monkey_FDRlt0.1_Up10_Down10.pdf")),
                   title = paste0("Human vs Monkey (", ct, ")  |  FDR<", alpha_fdr, " in BOTH"))
    } 
  }
  
  # Human vs Mouse
  if (!is.null(H) && !is.null(mouse_meta[[ct]])) {
    L <- make_bubble_long(H, mouse_meta[[ct]], "Mouse",
                          fdr_thr = alpha_fdr, n_per_direction = n_per_dir)
    if (!is.null(L) && nrow(L)) {
      plot_bubbles(L,
                   file  = file.path("pairwise_bubbles", paste0(ct, "_Human_vs_Mouse_FDRlt0.1_Up10_Down10.pdf")),
                   title = paste0("Human vs Mouse (", ct, ")  |  FDR<", alpha_fdr, " in BOTH"))
    } 
  }
  
  # Human vs Goat
  if (!is.null(H) && !is.null(goat_tbl[[ct]])) {
    L <- make_bubble_long(H, goat_tbl[[ct]], "Goat",
                          fdr_thr = alpha_fdr, n_per_direction = n_per_dir)
    if (!is.null(L) && nrow(L)) {
      plot_bubbles(L,
                   file  = file.path("pairwise_bubbles", paste0(ct, "_Human_vs_Goat_FDRlt0.1_Up10_Down10.pdf")),
                   title = paste0("Human vs Goat (", ct, ")  |  FDR<", alpha_fdr, " in BOTH"))
    } 
  }
}

###############################################################################
sink(file = paste0(Sys.Date(), "_DEG_metaRNASeq_analysis_session_info.txt"))
sessionInfo()
sink()
