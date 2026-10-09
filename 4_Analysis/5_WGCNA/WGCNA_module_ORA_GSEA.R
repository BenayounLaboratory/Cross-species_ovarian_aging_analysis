options(stringsAsFactors = FALSE)

library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(scales)
library(dplyr)
library(clusterProfiler)
library(patchwork)

theme_set(theme_bw())

################################################################################
# Integrative ovarian aging analysis
# Run ORA and GSEA of WGCNA modules
# Cross-species analysis
################################################################################

################################################################################
# 0. Configuration
################################################################################

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE         <- "/Volumes/jinho01/Benayoun_lab/Projects/CZI"
RESULTS_DIR  <- file.path(BASE, "5_WGCNA/CrossSpecies_WGCNA_Results")
GMT_FILE     <- file.path(RESULTS_DIR, "GMT_Export/Human_AgingModules.gmt")
WGCNA_RDS    <- file.path(RESULTS_DIR, "wgcna_results_all_species.rds")
DESEQ_DIR    <- file.path(BASE, "1_PB_DESeq2/Humanized_DESeq2")

ORA_DIR      <- file.path(RESULTS_DIR, "ORA_Human")
GSEA_DIR     <- file.path(RESULTS_DIR, "GSEA_Human_vs_OtherSpecies")
dir.create(ORA_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(GSEA_DIR, showWarnings = FALSE, recursive = TRUE)

# All datasets to test (human + non-human; all have humanized t-statistics)
NON_HUMAN_DATASETS <- c(
  "Mouse_AC", "Mouse_EMTAB11491", "Mouse_EMTAB12889",
  "Mouse_Foxl2_wt", "Mouse_GSE232309", "Mouse_GSE290742", "Mouse_VCD",
  "Monkey_GSE130664", "Goat_PRJNA1010653"
)
HUMAN_DATASETS <- c("Human_GSE202601", "Human_GSE255690", "Human_TS")

# Cell types with human aging modules (Stroma has 0; skip it)
TARGET_CELL_TYPES <- c("Granulosa", "Theca")

# Fixed NES color scale for GSEA bubble plots
NES_LIMITS <- c(-3, 3)
NES_VALUES <- rescale(c(-3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3))
NES_COLORS <- c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1",
                "white", "lightcoral", "brown1", "firebrick2", "firebrick4")

# Human datasets are listed first; non-human follow
DATASET_ORDER <- c(
  "Human_GSE202601", "Human_GSE255690", "Human_TS",
  "Mouse_AC", "Mouse_EMTAB11491", "Mouse_EMTAB12889",
  "Mouse_Foxl2_wt", "Mouse_GSE232309", "Mouse_GSE290742", "Mouse_VCD",
  "Monkey_GSE130664", "Goat_PRJNA1010653"
)

###############################################################################
# Helper: parse GMT → named list of character vectors
###############################################################################
read_gmt_list <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  sets  <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    list(name = parts[1], genes = parts[-c(1, 2)])
  })
  out        <- lapply(sets, `[[`, "genes")
  names(out) <- sapply(sets, `[[`, "name")
  out
}

# Extract cell type token from set name "Human__Granulosa__MEturquoise"
set_cell_type <- function(nm) strsplit(nm, "__", fixed = TRUE)[[1]][2]

###############################################################################
# Helper: build ranked gene list from DESeq2 stat column
###############################################################################
build_gene_list <- function(res) {
  if (is.null(res)) return(NULL)
  df <- as.data.frame(res)
  if (!"stat" %in% colnames(df)) return(NULL)
  gl <- setNames(df$stat, rownames(df))
  gl <- gl[!is.na(gl) & is.finite(gl)]
  if (any(duplicated(names(gl)))) {
    gl <- gl[order(abs(gl), decreasing = TRUE)]
    gl <- gl[!duplicated(names(gl))]
  }
  sort(gl, decreasing = TRUE)
}

###############################################################################
# 1. Load inputs
###############################################################################

wgcna_all  <- readRDS(WGCNA_RDS)
human_sets <- read_gmt_list(GMT_FILE)

# Per cell type: which human datasets actually contributed to the WGCNA network.
human_valid_ds <- lapply(wgcna_all[["Human"]], function(x) x$contributing_datasets)

###############################################################################
# Part 1: ORA — enrichGO (GO:BP) for each human aging module
###############################################################################

ora_results <- list()

for (ct in TARGET_CELL_TYPES) {

  genes_universe <- wgcna_all[["Human"]][[ct]]$genes
  if (is.null(genes_universe)) {
    next
  }

  # Modules for this cell type
  ct_sets <- human_sets[sapply(names(human_sets), set_cell_type) == ct]
  if (!length(ct_sets)) {
    next
  }

  ora_results[[ct]] <- list()

  for (mod_name in names(ct_sets)) {
    mod_genes <- unique(ct_sets[[mod_name]])
    message("  Module: ", mod_name, " (", length(mod_genes), " genes)")

    ego <- tryCatch(
      enrichGO(
        gene          = mod_genes,
        universe      = genes_universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,
        qvalueCutoff  = 1,
        minGSSize     = 5,
        maxGSSize     = 2000,
        readable      = FALSE
      ),
      error = function(e) { message("    enrichGO error: ", e$message); NULL }
    )

    if (is.null(ego) || nrow(ego@result) == 0) {
      next
    }

    res_df <- ego@result %>% arrange(p.adjust, pvalue)
    n_sig  <- sum(res_df$p.adjust < 0.05, na.rm = TRUE)

    # Save CSV
    safe_name <- gsub("[^A-Za-z0-9_]", "_", mod_name)
    write.csv(res_df,
              file.path(ORA_DIR, paste0(Sys.Date(), "_", safe_name, "_ORA_GO_BP.csv")),
              row.names = FALSE)

    ora_results[[ct]][[mod_name]] <- res_df
  }

  # Per-cell-type: per-module bubble plots
  if (!length(ora_results[[ct]])) next

  for (mod_name in names(ora_results[[ct]])) {
    df <- ora_results[[ct]][[mod_name]] %>%
      filter(p.adjust < 0.2) %>%
      arrange(p.adjust) %>%
      slice_head(n = 20) %>%
      mutate(
        Term   = stringr::str_trunc(paste0(ID, " | ", Description), 55),
        log10p = -log10(p.adjust)
      )
    if (!nrow(df)) next

    p <- ggplot(df, aes(x = log10p, y = reorder(Term, log10p))) +
      geom_point(aes(size = Count, color = log10p)) +
      scale_size_area(max_size = 8) +
      scale_color_gradient(low = "grey80", high = "purple4") +
      labs(
        title = mod_name,
        x     = expression(-log[10]("FDR")),
        y     = NULL, size = "Overlap", color = expression(-log[10]("FDR"))
      ) +
      theme(axis.text.y = element_text(size = 8),
            plot.title  = element_text(size = 9, face = "bold"))

    safe_name <- gsub("[^A-Za-z0-9_]", "_", mod_name)
    ggsave(file.path(ORA_DIR, paste0(Sys.Date(), "_", safe_name, "_ORA_bubble.pdf")),
           p, width = 8, height = 5, useDingbats = FALSE)
  }

  # Combined faceted plot across all modules in this cell type
  all_top <- bind_rows(
    lapply(names(ora_results[[ct]]), function(nm) {
      ora_results[[ct]][[nm]] %>%
        filter(p.adjust < 0.2) %>%
        arrange(p.adjust) %>%
        slice_head(n = 10) %>%
        mutate(
          Module = nm,
          Term   = stringr::str_trunc(paste0(ID, " | ", Description), 50),
          log10p = -log10(p.adjust)
        )
    })
  )

  if (nrow(all_top)) {
    p_comb <- ggplot(all_top, aes(x = log10p, y = reorder(Term, log10p))) +
      geom_point(aes(size = Count, color = log10p)) +
      scale_size_area(max_size = 6) +
      scale_color_gradient(low = "grey80", high = "purple4") +
      facet_wrap(~ Module, scales = "free_y", ncol = 2) +
      labs(
        title = paste0("ORA GO:BP — Human ", ct, " aging modules"),
        x     = expression(-log[10]("FDR")),
        y     = NULL, size = "Overlap", color = expression(-log[10]("FDR"))
      ) +
      theme(axis.text.y = element_text(size = 7), strip.text = element_text(size = 8))

    ggsave(file.path(ORA_DIR, paste0(Sys.Date(), "_Human_", ct, "_ORA_combined.pdf")),
           p_comb,
           width  = 14,
           height = max(6, 0.4 * length(unique(all_top$Term))),
           useDingbats = FALSE)
  }
}

# Save all ORA results as a single RDS
saveRDS(ora_results, file.path(ORA_DIR, paste0(Sys.Date(), "_ORA_results.rds")))

###############################################################################
# Part 2: GSEA — human aging modules vs non-human DESeq2 t-statistics
###############################################################################

gsea_results <- list()

for (ct in TARGET_CELL_TYPES) {
  message("\n  Cell type: ", ct)

  # Build TERM2GENE for this cell type
  ct_sets <- human_sets[sapply(names(human_sets), set_cell_type) == ct]
  if (!length(ct_sets)) {
    next
  }
  t2g <- bind_rows(
    lapply(names(ct_sets), function(nm)
      data.frame(gs_name = nm, gene_symbol = ct_sets[[nm]], stringsAsFactors = FALSE))
  )

  gsea_results[[ct]] <- list()

  # Combine non-human datasets with human datasets valid for this cell type
  valid_human <- intersect(HUMAN_DATASETS, human_valid_ds[[ct]])
  all_datasets <- c(NON_HUMAN_DATASETS, valid_human)

  for (ds in all_datasets) {
    rds_file <- file.path(DESEQ_DIR, paste0("2026-05-26_", ds, "_DESeq2_results.rds"))
    if (!file.exists(rds_file)) {
      next
    }

    deseq_obj <- tryCatch(readRDS(rds_file),
                          error = function(e) { message("  Load error: ", e$message); NULL })
    if (is.null(deseq_obj)) next

    # Extract this cell type's age contrast
    ct_res <- deseq_obj[[ct]]
    if (is.null(ct_res)) {
      next
    }

    # Get age contrast (first one with a stat column)
    age_res <- ct_res[["age"]]
    if (is.null(age_res)) {
      age_res <- NULL
      for (contrast_name in names(ct_res)) {
        cand <- as.data.frame(ct_res[[contrast_name]])
        if ("stat" %in% colnames(cand)) { age_res <- cand; break }
      }
    }

    gene_list <- build_gene_list(age_res)
    if (is.null(gene_list) || length(gene_list) < 10) {
      next
    }

    gsea_obj <- tryCatch(
      suppressMessages(
        GSEA(
          geneList     = gene_list,
          TERM2GENE    = t2g,
          minGSSize    = 5,
          maxGSSize    = 10000,
          pvalueCutoff = 1,
          verbose      = FALSE
        )
      ),
      error = function(e) { message("    GSEA error: ", e$message); NULL }
    )

    if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
      gsea_results[[ct]][[ds]] <- NULL
      next
    }

    res_tab <- gsea_obj@result %>%
      arrange(p.adjust, desc(abs(NES)))
    n_sig <- sum(res_tab$p.adjust < 0.1, na.rm = TRUE)

    # Save CSV
    write.csv(res_tab,
              file.path(GSEA_DIR,
                        paste0(Sys.Date(), "_GSEA_", ct, "_", ds, ".csv")),
              row.names = FALSE)

    gsea_results[[ct]][[ds]] <- res_tab

    # Per-dataset bubble plot (top 5 pos + top 5 neg by NES, among FDR<0.1)
    plot_df <- res_tab %>%
      filter(is.finite(NES), p.adjust < 0.1) %>%
      { bind_rows(
          arrange(filter(., NES > 0), desc(NES)) %>% slice_head(n = 5),
          arrange(filter(., NES < 0), NES)       %>% slice_head(n = 5)
        ) } %>%
      mutate(
        PathName   = paste0(ID, " ", Description),
        log10fdr   = -log10(p.adjust)
      )

    if (nrow(plot_df)) {
      p <- ggplot(plot_df,
                  aes(x = 1, y = reorder(PathName, log10fdr),
                      color = NES, size = log10fdr)) +
        geom_point(shape = 16, alpha = 0.95) +
        scale_color_gradientn(colors = NES_COLORS, values = NES_VALUES,
                              limits = NES_LIMITS, oob = squish) +
        labs(title = paste("GSEA |", ct, "|", ds),
             x = NULL, y = NULL, color = "NES",
             size = expression(-log[10]("FDR"))) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 8),
              plot.title  = element_text(size = 9, face = "bold"))
      ggsave(file.path(GSEA_DIR,
                       paste0(Sys.Date(), "_GSEA_bubble_", ct, "_", ds, ".pdf")),
             p, width = 10, height = 5, useDingbats = FALSE)
    }
  }
}

# Save all GSEA results
saveRDS(gsea_results, file.path(GSEA_DIR, paste0(Sys.Date(), "_GSEA_results.rds")))

###############################################################################
# Session info
###############################################################################
sink("./WGCNA_module_ORA_GSEA_session_info")
sessionInfo()
sink()