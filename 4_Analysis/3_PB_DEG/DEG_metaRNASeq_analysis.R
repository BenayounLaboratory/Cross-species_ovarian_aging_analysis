# Load libraries

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DESeq2)
  
# rm(list = ls())
  
###############################################################################
# Integrative ovarian aging analysis
# Within-species meta-analysis (one-sided Fisher's method) of age-associated
# DEGs across ovarian scRNA-seq datasets, followed by cross-species integration.
###############################################################################

###############################################################################
# 1. Load data
###############################################################################

# DESeq2 results
load("/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/0_Humanized_DESeq2_objects/2026-08-06_all_combined_cell_type_DESeq2_results.RData")

# Set parameters
HUMANIZED    <- "/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/Humanized_DESeq2"
GOAT_TSV_DIR <- "/Volumes/jinho01/Benayoun_lab/Projects/CZI/1_PB_DESeq2/SVA_DESeq2_Final/Goat_PRJNA1010653/DESeq2_results"
GOAT_ORTH    <- file.path(HUMANIZED, "2026-05-27_ortholog_goat_to_human.csv")

FDR_THRESHOLD    <- 0.1
MIN_STUDIES_HUM  <- 2L   # out of 3
MIN_STUDIES_MOUS <- 3L   # out of 7

HUMAN_DATASETS  <- c("Human_GSE202601", "Human_GSE255690", "Human_TS")
MOUSE_DATASETS  <- c("Mouse_AC", "Mouse_VCD", "Mouse_EMTAB128899",
                     "Mouse_GSE232309", "Mouse_EMTAB11491",
                     "Mouse_Foxl2_wt", "Mouse_GSE290742")
MONKEY_DATASET  <- "Monkey_GSE130664"

cell_type_order <- c(
  "Oocyte", "Granulosa", "Theca", "Stroma", "SMC", "BEC", "LEC",
  "Epithelial", "Myeloid", "DC",
  "ILC", "NK", "NKT", "CD8NKT", "CD8T", "CD4T", "DNT", "DPT", "B"
)

################################################################################
# 1. Load helpers
################################################################################

# Extract the results
extract_age_results <- function(rds_path) {
  obj <- readRDS(rds_path)
  lapply(obj, function(ct_data) {
    res <- ct_data[["age"]]
    # Ensure required columns
    req <- c("pvalue", "stat", "log2FoldChange", "padj")
    if (!all(req %in% colnames(res))) return(NULL)
    res[!is.na(res$pvalue) & !is.na(res$stat), req, drop = FALSE]
  })
}

# Load all datasets for one species group
load_species_group <- function(datasets) {
  out <- list()
  for (ds in datasets) {
    f <- list.files(HUMANIZED, pattern = paste0(ds, "_DESeq2_results\\.rds$"),
                    full.names = TRUE)
    if (length(f) == 0) { message("  No RDS for ", ds); next }
    message("  Loading ", basename(f[1]))
    out[[ds]] <- extract_age_results(f[1])
  }
  out
}

################################################################################
# 2. Fisher's one-sided meta-analysis per cell type
################################################################################

# One gene, one direction: Fisher's method on a numeric vector of p-values
# (only tested studies passed in; df = 2 * length(pvalues))
fisher_p <- function(pvalues) {
  pvalues <- pmax(pvalues, 1e-300)
  pchisq(-2 * sum(log(pvalues)), df = 2L * length(pvalues), lower.tail = FALSE)
}

# Returns list[cell_type] -> data.frame (one row per gene)
run_species_meta <- function(species_data, min_studies) {

  datasets   <- names(species_data)
  cell_types <- unique(unlist(lapply(species_data, names)))

  lapply(setNames(cell_types, cell_types), function(ct) {

    # Gather results for this cell type from each dataset
    ct_list <- Filter(Negate(is.null),
                      lapply(species_data, function(ds_res) ds_res[[ct]]))

    if (length(ct_list) < min_studies) return(NULL)

    all_genes <- unique(unlist(lapply(ct_list, rownames)))
    n_ds      <- length(ct_list)

    # Matrices: rows = genes, cols = datasets
    p_mat  <- matrix(NA_real_, nrow = length(all_genes), ncol = n_ds,
                     dimnames = list(all_genes, names(ct_list)))
    st_mat <- matrix(NA_real_, nrow = length(all_genes), ncol = n_ds,
                     dimnames = list(all_genes, names(ct_list)))
    lf_mat <- matrix(NA_real_, nrow = length(all_genes), ncol = n_ds,
                     dimnames = list(all_genes, names(ct_list)))

    for (ds in names(ct_list)) {
      df   <- ct_list[[ds]]
      gg   <- intersect(rownames(df), all_genes)
      p_mat[gg, ds]  <- df[gg, "pvalue"]
      st_mat[gg, ds] <- df[gg, "stat"]
      lf_mat[gg, ds] <- df[gg, "log2FoldChange"]
    }

    n_tested <- rowSums(!is.na(p_mat))
    keep     <- n_tested >= min_studies
    if (!any(keep)) return(NULL)

    p_mat  <- p_mat[keep, , drop = FALSE]
    st_mat <- st_mat[keep, , drop = FALSE]
    lf_mat <- lf_mat[keep, , drop = FALSE]
    n_tested <- n_tested[keep]

    # One-sided p-values from Wald stat
    # p_up   = P(Z > stat_i) = pnorm(-stat_i)   [small when gene is upregulated]
    # p_down = P(Z < stat_i) = pnorm(stat_i)    [small when gene is downregulated]
    p_up_mat   <- p_down_mat <- matrix(NA_real_, nrow = nrow(p_mat), ncol = ncol(p_mat),
                                        dimnames = dimnames(p_mat))
    mask <- !is.na(st_mat)
    p_up_mat[mask]   <- pnorm(-st_mat[mask])
    p_down_mat[mask] <- pnorm( st_mat[mask])

    # Fisher's combined p (one-sided, per gene, over tested studies only)
    fp_up <- vapply(seq_len(nrow(p_up_mat)), function(i) {
      pv <- p_up_mat[i, ]; pv <- pv[!is.na(pv)]; fisher_p(pv)
    }, numeric(1))
    fp_down <- vapply(seq_len(nrow(p_down_mat)), function(i) {
      pv <- p_down_mat[i, ]; pv <- pv[!is.na(pv)]; fisher_p(pv)
    }, numeric(1))

    # BH correction separately for up and down
    meta_padj_up   <- p.adjust(fp_up,   method = "BH")
    meta_padj_down <- p.adjust(fp_down, method = "BH")

    # Summary log2FC and direction; unname() prevents named vectors becoming list columns
    median_l2fc    <- unname(rowMedians(lf_mat, na.rm = TRUE))
    n_up           <- unname(rowSums(lf_mat > 0, na.rm = TRUE))
    n_down         <- unname(rowSums(lf_mat < 0, na.rm = TRUE))
    n_tested_v     <- unname(n_tested)

    # Direction: which one-sided test is stronger
    meta_direction <- ifelse(fp_up <= fp_down, "up", "down")

    meta_padj_up_sc   <- ifelse(median_l2fc > 0, meta_padj_up,   1)
    meta_padj_down_sc <- ifelse(median_l2fc < 0, meta_padj_down, 1)

    tibble(
      gene           = rownames(p_mat),
      n_tested       = n_tested_v,
      n_up           = n_up,
      n_down         = n_down,
      median_log2FC  = median_l2fc,
      fisher_p_up    = fp_up,
      fisher_p_down  = fp_down,
      meta_padj_up   = meta_padj_up_sc,
      meta_padj_down = meta_padj_down_sc,
      meta_direction = meta_direction,
      meta_sig       = (meta_padj_up_sc < FDR_THRESHOLD) | (meta_padj_down_sc < FDR_THRESHOLD)
    )
  })
}

################################################################################
# 3. Load Goat data (TSV + manual humanization)
################################################################################

load_goat_humanized <- function() {
  # Read all-celltypes TSV
  tsv_f <- list.files(GOAT_TSV_DIR, pattern = "_age_all_celltypes\\.tsv$",
                      full.names = TRUE)
  if (length(tsv_f) == 0) stop("Goat all-celltypes TSV not found.")
  goat_df <- read.table(tsv_f[1], header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
  # Normalize cell type names
  goat_df$celltype <- gsub("CD8 NKT", "CD8NKT", goat_df$celltype)

  # Load ortholog table
  orth <- read.csv(GOAT_ORTH, stringsAsFactors = FALSE, row.names = 1)

  # Humanize: map goat gene names to human symbols
  goat_df <- goat_df %>%
    filter(gene %in% rownames(orth)) %>%
    mutate(human_gene = orth[gene, "hsapiens_homolog_associated_gene_name"]) %>%
    filter(human_gene != "") %>%
    # keep first occurrence per (cell type, human gene)
    group_by(celltype, human_gene) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()

  ct_list <- split(goat_df, goat_df$celltype)
  lapply(ct_list, function(d) {
    df <- d %>% select(human_gene, log2FoldChange, pvalue) %>%
      as.data.frame()
    rownames(df) <- df$human_gene
    df
  })
}

################################################################################
# 4. Main: load all data and run meta-analyses
################################################################################

human_data  <- load_species_group(HUMAN_DATASETS)

mouse_data  <- load_species_group(MOUSE_DATASETS)

monkey_data <- load_species_group(MONKEY_DATASET)
monkey_res  <- monkey_data[[MONKEY_DATASET]] 

message("Loading Goat data (from TSV + ortholog table)...")
goat_res <- tryCatch(load_goat_humanized(),
                     error = function(e) {
                       message("  Goat load failed: ", conditionMessage(e))
                       list()
                     })

human_meta <- run_species_meta(human_data, min_studies = MIN_STUDIES_HUM)
mouse_meta <- run_species_meta(mouse_data, min_studies = MIN_STUDIES_MOUS)

################################################################################
# 5. Cross-species integration
################################################################################

# All cell types present in human or mouse meta
all_ct <- union(names(human_meta), names(mouse_meta))
all_ct <- all_ct[!sapply(all_ct, function(ct)
  is.null(human_meta[[ct]]) && is.null(mouse_meta[[ct]]))]

cross_species_list <- lapply(setNames(all_ct, all_ct), function(ct) {

  hm <- human_meta[[ct]]   
  mm <- mouse_meta[[ct]]  
  mk <- monkey_res[[ct]] 
  gk <- goat_res[[ct]]   

  # Collect all genes appearing in any species
  genes <- unique(c(
    if (!is.null(hm)) hm$gene,
    if (!is.null(mm)) mm$gene
  ))
  if (length(genes) == 0) return(NULL)

  # Helpers to pull columns safely
  pull_col <- function(df, g, col, default = NA_real_) {
    if (is.null(df)) return(rep(default, length(g)))
    vals <- df[[col]][match(g, df$gene)]
    ifelse(is.na(vals), default, vals)
  }
  pull_chr <- function(df, g, col, default = NA_character_) {
    if (is.null(df)) return(rep(default, length(g)))
    vals <- df[[col]][match(g, df$gene)]
    ifelse(is.na(vals), default, vals)
  }

  # Monkey: single study, use raw padj and log2FC direction
  mk_padj <- if (!is.null(mk)) mk[match(genes, rownames(mk)), "padj"]   else rep(NA_real_, length(genes))
  mk_l2fc <- if (!is.null(mk)) mk[match(genes, rownames(mk)), "log2FoldChange"] else rep(NA_real_, length(genes))
  mk_dir  <- ifelse(is.na(mk_l2fc), NA_character_,
                    ifelse(mk_l2fc > 0, "up", "down"))

  # Goat: direction and nominal pvalue only
  gk_pval <- if (!is.null(gk)) gk[match(genes, rownames(gk)), "pvalue"] else rep(NA_real_, length(genes))
  gk_l2fc <- if (!is.null(gk)) gk[match(genes, rownames(gk)), "log2FoldChange"] else rep(NA_real_, length(genes))
  gk_dir  <- ifelse(is.na(gk_l2fc), NA_character_,
                    ifelse(gk_l2fc > 0, "up", "down"))

  df <- tibble(
    gene               = genes,
    cell_type          = ct,

    # Human meta
    human_n_studies    = pull_col(hm, genes, "n_tested",       0),
    human_median_l2fc  = pull_col(hm, genes, "median_log2FC"),
    human_padj_up      = pull_col(hm, genes, "meta_padj_up",   1),
    human_padj_down    = pull_col(hm, genes, "meta_padj_down", 1),
    human_direction    = pull_chr(hm, genes, "meta_direction"),
    human_meta_sig     = (human_padj_up < FDR_THRESHOLD) | (human_padj_down < FDR_THRESHOLD),

    # Mouse meta
    mouse_n_studies    = pull_col(mm, genes, "n_tested",       0),
    mouse_median_l2fc  = pull_col(mm, genes, "median_log2FC"),
    mouse_padj_up      = pull_col(mm, genes, "meta_padj_up",   1),
    mouse_padj_down    = pull_col(mm, genes, "meta_padj_down", 1),
    mouse_direction    = pull_chr(mm, genes, "meta_direction"),
    mouse_meta_sig     = (mouse_padj_up < FDR_THRESHOLD) | (mouse_padj_down < FDR_THRESHOLD),

    # Monkey (single study)
    monkey_padj        = mk_padj,
    monkey_l2fc        = mk_l2fc,
    monkey_direction   = mk_dir,
    monkey_sig         = !is.na(mk_padj) & mk_padj < FDR_THRESHOLD,

    # Goat (directional annotation)
    goat_pvalue        = gk_pval,
    goat_l2fc          = gk_l2fc,
    goat_direction     = gk_dir
  ) %>%
    mutate(
      # Direction consistent between human and mouse (both must be sig and agree)
      hum_mou_same_dir = !is.na(human_direction) & !is.na(mouse_direction) &
                         human_direction == mouse_direction,

      # Directional human padj (the significant side)
      human_meta_padj = case_when(
        human_direction == "up"   ~ human_padj_up,
        human_direction == "down" ~ human_padj_down,
        TRUE ~ pmin(human_padj_up, human_padj_down)
      ),
      mouse_meta_padj = case_when(
        mouse_direction == "up"   ~ mouse_padj_up,
        mouse_direction == "down" ~ mouse_padj_down,
        TRUE ~ pmin(mouse_padj_up, mouse_padj_down)
      ),

      # Tier assignment
      tier = case_when(
        # Tier 1: both species significant + same direction
        human_meta_sig & mouse_meta_sig & hum_mou_same_dir
          ~ "Tier1_conserved",

        # Tier 2: one species significant + monkey agrees in direction
        (human_meta_sig | mouse_meta_sig) & monkey_sig &
          ((!is.na(human_direction) & !is.na(monkey_direction) & human_direction == monkey_direction) |
           (!is.na(mouse_direction) & !is.na(monkey_direction) & mouse_direction == monkey_direction))
          ~ "Tier2_partial",

        # Significant in at least one species
        human_meta_sig | mouse_meta_sig
          ~ "Tier3_single_species",

        TRUE ~ "NS"
      )
    )

  df
})

cross_species_df <- bind_rows(cross_species_list)

################################################################################
# 6. Save outputs
################################################################################

# Human meta results (all cell types, one file)
human_meta_df <- bind_rows(
  lapply(names(human_meta), function(ct) {
    if (is.null(human_meta[[ct]])) return(NULL)
    mutate(human_meta[[ct]], cell_type = ct)
  })
) %>% select(cell_type, everything())

mouse_meta_df <- bind_rows(
  lapply(names(mouse_meta), function(ct) {
    if (is.null(mouse_meta[[ct]])) return(NULL)
    mutate(mouse_meta[[ct]], cell_type = ct)
  })
) %>% select(cell_type, everything())

readr::write_tsv(human_meta_df,
  paste0(Sys.Date(), "_human_meta_DEGs.tsv"))
readr::write_tsv(mouse_meta_df,
  paste0(Sys.Date(), "_mouse_meta_DEGs.tsv"))
readr::write_tsv(cross_species_df,
  paste0(Sys.Date(), "_cross_species_DEGs.tsv"))
message("Main tables written.")

human_meta_df    <- readr::read_tsv(paste0(Sys.Date(), "_human_meta_DEGs.tsv"),
                                     show_col_types = FALSE)
mouse_meta_df    <- readr::read_tsv(paste0(Sys.Date(), "_mouse_meta_DEGs.tsv"),
                                     show_col_types = FALSE)
cross_species_df <- readr::read_tsv(paste0(Sys.Date(), "_cross_species_DEGs.tsv"),
                                     show_col_types = FALSE)

###############################################################################
sink(file = paste0(Sys.Date(), "_DEG_metaRNASeq_analysis_session_info.txt"))
sessionInfo()
sink()
