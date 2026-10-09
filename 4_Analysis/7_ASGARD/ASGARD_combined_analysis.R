library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
  
rm(list = ls())

################################################################################
# Integrative ovarian aging analysis
# Combined analysis of ASGARD outputs
################################################################################

############################################
# 1. Load data
############################################

load("./CZI_ASGARD_results_combined.RData")

## Objects present:
##   asgard.human.GSE202601, asgard.human.GSE255690, asgard.human.TS
##   asgard.monkey
##   asgard.mouse.aging, asgard.mouse.Foxl2, asgard.mouse.VCD,
##   asgard.mouse.GSE232309, asgard.mouse.EMTAB128899,
##   asgard.mouse.EMTAB11491, asgard.mouse.GSE290742
##   asgard.goat

############################################
# 2. Helper functions
############################################

.zerop <- function(p) pmin(pmax(p, .Machine$double.eps), 1 - 1e-15)

combine_fisher <- function(p) {
  p <- .zerop(p); k <- length(p)
  if (!k) return(NA_real_)
  stat <- -2 * sum(log(p))
  pchisq(stat, df = 2 * k, lower.tail = FALSE)
}

## Signed Stouffer: aligns Z to the beneficial (positive score) direction
combine_stouffer_signed <- function(p, score, w = NULL, target_sign = +1L) {
  p <- .zerop(p)
  s <- sign(score) * target_sign
  z <- qnorm(1 - p) * s
  if (is.null(w)) w <- rep(1, length(z))
  zc <- sum(w * z) / sqrt(sum(w^2))
  2 * pnorm(-abs(zc))
}

binom_overrep <- function(sig_vec, p0) {
  k <- sum(sig_vec, na.rm = TRUE)
  n <- sum(!is.na(sig_vec))
  if (n == 0) return(NA_real_)
  pbinom(q = k - 1, size = n, prob = p0, lower.tail = FALSE)
}

mlog <- function(p, cap = 10) {
  x <- -log10(p)
  x[!is.finite(x)] <- cap
  pmin(pmax(x, 0), cap)
}

############################################
# 3. Assemble tidy data frame
############################################

datasets <- list(
  Human_GSE202601   = asgard.human.GSE202601,
  Human_GSE255690   = asgard.human.GSE255690,
  Human_TS          = asgard.human.TS,
  Monkey            = asgard.monkey,
  Mouse_Aging       = asgard.mouse.aging,
  Mouse_Foxl2       = asgard.mouse.Foxl2,
  Mouse_VCD         = asgard.mouse.VCD,
  Mouse_GSE232309   = asgard.mouse.GSE232309,
  Mouse_EMTAB128899 = asgard.mouse.EMTAB128899,
  Mouse_EMTAB11491  = asgard.mouse.EMTAB11491,
  Mouse_GSE290742   = asgard.mouse.GSE290742,
  Goat              = asgard.goat
)

species_map <- c(
  Human_GSE202601   = "Human",
  Human_GSE255690   = "Human",
  Human_TS          = "Human",
  Monkey            = "Monkey",
  Mouse_Aging       = "Mouse",
  Mouse_Foxl2       = "Mouse",
  Mouse_VCD         = "Mouse",
  Mouse_GSE232309   = "Mouse",
  Mouse_EMTAB128899 = "Mouse",
  Mouse_EMTAB11491  = "Mouse",
  Mouse_GSE290742   = "Mouse",
  Goat              = "Goat"
)

dataset_order <- names(datasets)

extract_asgard <- function(df, dataset_name) {
  if (is.null(df) || !nrow(df)) return(NULL)
  out <- df %>% as.data.frame()

  if (!("drug" %in% names(out))) {
    cand <- intersect(c("Drug", "DRUG", "drug", "Drug.name", "Drug.Name",
                        "compound", "Compound"), names(out))
    if (length(cand)) {
      out$drug <- out[[cand[1]]]
    } else {
      out$drug <- rownames(out)
    }
  }

  pcol <- intersect(c("P.value", "pvalue", "p_value", "P", "p", "P.Value"), names(out))
  scol <- intersect(c("Drug.therapeutic.score", "score", "therapeutic_score",
                       "tau", "effect"), names(out))
  if (!length(pcol) || !length(scol)) return(NULL)

  out %>%
    transmute(
      dataset = dataset_name,
      species = species_map[[dataset_name]],
      drug    = as.character(drug),
      pval    = suppressWarnings(as.numeric(.data[[pcol[1]]])),
      score   = suppressWarnings(as.numeric(.data[[scol[1]]]))
    ) %>%
    filter(!is.na(drug), is.finite(pval), is.finite(score)) %>%
    mutate(pval = .zerop(pval))
}

df <- purrr::imap_dfr(datasets, extract_asgard)
stopifnot(nrow(df) > 0)
df$dataset <- factor(df$dataset, levels = dataset_order, ordered = TRUE)

cat("Datasets loaded:", nlevels(df$dataset), "\n")
cat("Unique drugs   :", n_distinct(df$drug), "\n")

############################################
# 4. Per-dataset BH correction
############################################

df <- df %>%
  group_by(dataset) %>%
  mutate(
    padj = p.adjust(pval, "BH"),
    sig  = as.integer(padj <= 0.10)
  ) %>%
  ungroup()

## Background significance rate (median across datasets, used for over-rep test)
p0 <- df %>%
  group_by(dataset) %>%
  summarise(rate = mean(sig, na.rm = TRUE), .groups = "drop") %>%
  pull(rate) %>%
  median(na.rm = TRUE)
if (!is.finite(p0) || p0 == 0) p0 <- 0.05

############################################
# 5. Per-species meta-analysis (Fisher)
############################################

per_species <- df %>%
  group_by(drug, species) %>%
  summarise(
    fisher_p        = combine_fisher(pval),
    stouffer_signed = combine_stouffer_signed(pval, score),
    n_datasets      = n(),
    mean_score      = mean(score, na.rm = TRUE),
    .groups = "drop"
  )

############################################
# 6. Pan-species meta-analysis
############################################

pan <- per_species %>%
  group_by(drug) %>%
  summarise(
    pan_fisher_p    = combine_fisher(fisher_p),
    pan_signed_p    = combine_fisher(stouffer_signed),
    species_covered = n_distinct(species),
    total_datasets  = sum(n_datasets),
    .groups = "drop"
  )

############################################
# 7. Over-representation & directional consistency
############################################

overrep <- df %>%
  group_by(drug) %>%
  summarise(
    k_sig     = sum(sig, na.rm = TRUE),
    n         = n(),
    overrep_p = binom_overrep(sig, p0),
    .groups   = "drop"
  )

consistency <- df %>%
  mutate(sdir = sign(score)) %>%
  group_by(drug) %>%
  summarise(
    n_nonzero     = sum(sdir != 0, na.rm = TRUE),
    majority_sign = if (n_nonzero > 0)
                      as.integer(names(which.max(table(factor(sdir[sdir != 0],
                                                              levels = c(-1, 1))))))
                    else NA_integer_,
    consistency   = if (n_nonzero > 0)
                      max(table(sdir[sdir != 0])) / n_nonzero
                    else NA_real_,
    .groups = "drop"
  )

############################################
# 8. Assemble master table + BH adjust
############################################

res <- pan %>%
  left_join(overrep,      by = "drug") %>%
  left_join(consistency,  by = "drug") %>%
  mutate(
    q_fisher  = p.adjust(pan_fisher_p, "BH"),
    q_signed  = p.adjust(pan_signed_p, "BH"),
    q_overrep = p.adjust(overrep_p,   "BH")
  ) %>%
  arrange(q_fisher, q_signed)

## Output directory
out_dir <- "2026-09-17/asgard_meta"
dir.create(out_dir, showWarnings = FALSE)

write.csv(res,
          file = file.path(out_dir, "ASGARD_meta_ranked_drugs.csv"),
          row.names = FALSE)

############################################
sink(file = "./ASGARD_combined_analysis_session_info.txt")
sessionInfo()
sink()
