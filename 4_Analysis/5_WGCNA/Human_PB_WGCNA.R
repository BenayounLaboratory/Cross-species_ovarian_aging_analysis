library(WGCNA)

################################################################################
# Integrative ovarian aging analysis
# Run WGCNA
# Process human datasets
################################################################################

################################################################################
# 0. Configuration
################################################################################

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ---- Paths ----
INPUT_DIR  <- "/Volumes/OIProject_II/1_R/4_Analysis/4_WGCNA/WGCNA_input"
BASE_DIR   <- "/Volumes/OIProject_II/1_R/4_Analysis/4_WGCNA/"
OUTPUT_DIR <- file.path(BASE_DIR, paste0(Sys.Date(), "_CrossSpecies_WGCNA_Results"))
HUMAN_DIR  <- file.path(OUTPUT_DIR, "WGCNA_Human")

# ---- WGCNA parameters ----
NETWORK_TYPE     <- "signed"
COR_TYPE         <- "bicor"
MIN_MODULE_SIZE  <- 30
MERGE_CUT_HEIGHT <- 0.25
DEEP_SPLIT       <- 2
MAX_BLOCK_SIZE   <- 30000
MIN_GENES_QC     <- 500
MIN_SAMPLES_QC   <- 3
CELL_TYPES_USE   <- c("Granulosa", "Theca", "Stroma")
SOFT_POW_DEFAULT <- 12

# ---- Datasets ----
DATASET_FILES <- c(
  Human_GSE202601 = "Human_GSE202601_PB_WGCNA_log2_counts.rds",
  Human_GSE255690 = "Human_GSE255690_PB_WGCNA_log2_counts.rds",
  Human_TS        = "Human_TS_PB_WGCNA_log2_counts.rds"
)

# ---- Age metadata (years) ----
AGE_METADATA <- list(
  Human_GSE202601 = c(
    Human23 = 23, Human27 = 27, Human28 = 28, Human29 = 29,
    Human49 = 49, Human51 = 51, Human52 = 52, Human54 = 54
  ),
  Human_GSE255690 = c(
    Young_1 = 23, Young_3 = 23, Young_4 = 23,
    Middle_2 = 38, Middle_3 = 38, Middle_4 = 38,
    Old_1 = 48, Old_2 = 48, Old_3 = 48
  ),
  Human_TS = c(
    TSP27 = 26, TSP28 = 55, TSP30 = 56
  )
)

################################################################################
# 1 : Load Data
################################################################################

dir.create(HUMAN_DIR, showWarnings = FALSE, recursive = TRUE)

.orient_matrix <- function(mat, dataset_name) {
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (nrow(mat) > ncol(mat) * 2) mat <- t(mat)
  mat
}

raw_data <- list()
for (ds_name in names(DATASET_FILES)) {
  obj <- readRDS(file.path(INPUT_DIR, DATASET_FILES[[ds_name]]))
  if (is.list(obj) && !is.data.frame(obj)) {
    obj <- lapply(obj, .orient_matrix, dataset_name = ds_name)
  } else {
    obj <- list(all_cells = .orient_matrix(obj, ds_name))
  }
  raw_data[[ds_name]] <- obj
}

################################################################################
# 2 : QC
################################################################################

qc_data <- list()   # qc_data[[dataset]][[cell_type]] = cleaned matrix

for (ds_name in names(raw_data)) {
  qc_data[[ds_name]] <- list()
  for (ct in intersect(names(raw_data[[ds_name]]), CELL_TYPES_USE)) {
    mat <- raw_data[[ds_name]][[ct]]

    # goodSamplesGenes needs >= 4 samples
    gsg <- tryCatch(
      goodSamplesGenes(mat, verbose = 0),
      error = function(e) {
        good_g <- apply(mat, 2, function(x) all(is.finite(x)) && var(x) > 0)
        good_s <- apply(mat, 1, function(x) all(is.finite(x)))
        list(allOK = all(good_g) && all(good_s),
             goodGenes = good_g, goodSamples = good_s)
      }
    )
    if (!gsg$allOK) mat <- mat[gsg$goodSamples, gsg$goodGenes]

    if (ncol(mat) < MIN_GENES_QC || nrow(mat) < MIN_SAMPLES_QC) {
      next
    }
    qc_data[[ds_name]][[ct]] <- mat
  }
}

################################################################################
# 3 : Soft-Threshold Selection
################################################################################

soft_powers <- list()

for (ds_name in names(qc_data)) {
  soft_powers[[ds_name]] <- list()
  for (ct in names(qc_data[[ds_name]])) {
    sft <- pickSoftThreshold(qc_data[[ds_name]][[ct]],
                             powerVector = seq(4, 20, by = 2),
                             networkType = NETWORK_TYPE,
                             corFnc      = COR_TYPE,
                             verbose     = 0)
    sp <- sft$powerEstimate
    if (is.na(sp)) sp <- SOFT_POW_DEFAULT
    soft_powers[[ds_name]][[ct]] <- sp
  }
}

################################################################################
# 4 : Consensus WGCNA per Cell Type
################################################################################

wgcna_results <- list()

for (ct in CELL_TYPES_USE) {

  available_ds <- names(qc_data)[sapply(names(qc_data),
                                        function(ds) !is.null(qc_data[[ds]][[ct]]))]
  if (length(available_ds) < 2) {
    next
  }

  # Intersect genes across datasets
  mats         <- lapply(available_ds, function(ds) qc_data[[ds]][[ct]])
  common_genes <- Reduce(intersect, lapply(mats, colnames))
  if (length(common_genes) < MIN_GENES_QC) {
    next
  }
  multiExpr <- lapply(mats, function(m) list(data = m[, common_genes, drop = FALSE]))
  names(multiExpr) <- available_ds

  # Most stringent (highest) power across datasets
  sp <- max(sapply(available_ds, function(ds) soft_powers[[ds]][[ct]]))

  consNet <- blockwiseConsensusModules(
    multiExpr         = multiExpr,
    power             = sp,
    networkType       = NETWORK_TYPE,
    corType           = COR_TYPE,
    maxBlockSize      = MAX_BLOCK_SIZE,
    minModuleSize     = MIN_MODULE_SIZE,
    mergeCutHeight    = MERGE_CUT_HEIGHT,
    deepSplit         = DEEP_SPLIT,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE,
    saveConsensusTOMs = FALSE,
    verbose           = 2
  )
  colors <- consNet$colors

  ct_out <- file.path(HUMAN_DIR, ct)
  dir.create(ct_out, showWarnings = FALSE, recursive = TRUE)

  pdf(file.path(ct_out, paste0(Sys.Date(), "_Human_", ct, "_consensus_dendrogram.pdf")),
      width = 12, height = 6)
  plotDendroAndColors(consNet$dendrograms[[1]],
                      colors[consNet$blockGenes[[1]]],
                      "Consensus Modules",
                      dendroLabels = FALSE,
                      hang         = 0.03,
                      addGuide     = TRUE,
                      main = paste("Consensus – Human", ct,
                                   "(", paste(available_ds, collapse = "+"), ")"))
  dev.off()

  result_obj <- list(
    colors                = colors,
    MEs                   = consNet$multiMEs,
    genes                 = common_genes,
    soft_power            = sp,
    contributing_datasets = available_ds,
    n_modules             = length(unique(colors[colors != "grey"]))
  )
  saveRDS(result_obj,
          file.path(ct_out, paste0(Sys.Date(), "_Human_", ct, "_wgcna_result.rds")))
  wgcna_results[[ct]] <- result_obj
}

saveRDS(wgcna_results,
        file.path(HUMAN_DIR, paste0(Sys.Date(), "_Human_all_cellTypes_wgcna.rds")))

################################################################################
# 5 : Module-Eigengene ~ Age Correlation
################################################################################
# Pearson per dataset (datasets with n < 4 give NA), combined across datasets
# by Fisher-z Stouffer's Z weighted by sqrt(n - 3). BH FDR < 0.1 = aging module.

.cor_age <- function(me_vec, age_vec) {
  ok <- complete.cases(me_vec, age_vec)
  if (sum(ok) < 4) return(c(cor = NA_real_, pval = NA_real_))
  ct <- cor.test(me_vec[ok], age_vec[ok], method = "pearson")
  c(cor = unname(ct$estimate), pval = ct$p.value)
}

.stouffer_combine <- function(cor_list, n_list) {
  rs <- sapply(cor_list, `[[`, "cor")
  ns <- unlist(n_list)
  ok <- !is.na(rs) & ns >= 4
  if (sum(ok) == 0) return(c(meta_cor = NA_real_, meta_pval = NA_real_))
  zs <- atanh(rs[ok])
  w  <- sqrt(ns[ok] - 3)
  Z  <- sum(w * zs) / sqrt(sum(w^2))
  c(meta_cor = tanh(sum(w * zs) / sum(w)), meta_pval = 2 * pnorm(-abs(Z)))
}

for (ct in names(wgcna_results)) {
  res          <- wgcna_results[[ct]]
  module_names <- setdiff(unique(res$colors), "grey")
  if (length(module_names) == 0) next

  cor_df <- do.call(rbind, lapply(module_names, function(mod) {
    cor_by_ds <- list()
    n_by_ds   <- list()
    for (ds in res$contributing_datasets) {
      me_mat <- res$MEs[[ds]]$data
      me_col <- paste0("ME", mod)
      if (!me_col %in% colnames(me_mat)) next
      age_vec <- AGE_METADATA[[ds]][rownames(me_mat)]
      cor_by_ds[[ds]] <- .cor_age(me_mat[, me_col], age_vec)
      n_by_ds[[ds]]   <- sum(!is.na(age_vec))
    }

    if (length(cor_by_ds) == 0)
      return(data.frame(module = mod, cor = NA, pval = NA))
    if (length(cor_by_ds) == 1)
      return(data.frame(module = mod,
                        cor    = unname(cor_by_ds[[1]]["cor"]),
                        pval   = unname(cor_by_ds[[1]]["pval"])))
    meta <- .stouffer_combine(cor_by_ds, n_by_ds)
    data.frame(module = mod,
               cor    = unname(meta["meta_cor"]),
               pval   = unname(meta["meta_pval"]))
  }))
  rownames(cor_df) <- NULL
  cor_df$padj           <- p.adjust(cor_df$pval, method = "BH")
  cor_df$aging_relevant <- !is.na(cor_df$padj) & cor_df$padj < 0.1

  write.csv(cor_df,
            file.path(HUMAN_DIR, ct,
                      paste0(Sys.Date(), "_Human_", ct, "_ModuleAge_correlation.csv")),
            row.names = FALSE)
}

################################################################################
sink(file = "./Human_PB_WGCNA_session_info.txt")
sessionInfo()
sink()