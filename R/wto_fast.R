# Vectorised re-implementation of wTO::wTO.fast() (wTO 2.1).
#
# It follows the package step by step: Spearman/Pearson correlation with a zero
# diagonal, signed wTO = (A %*% A + A) / (min(k_i, k_j) + 1 - |A|) rounded to 3
# decimals, bootstrap resampling of samples, and a p-value equal to the share
# of resamples in which wTO moved by more than `delta` (BH-adjusted).
# The only change is that min(k_i, k_j) is computed with pmin(outer()) instead
# of a double for-loop, which makes 2,000-gene networks practical.
# Returns the same columns as wTO.fast(): Node.1, Node.2, wTO, pval, pval.adj.

wto_matrix <- function(cor_mat) {
  A <- cor_mat
  k <- rowSums(abs(A), na.rm = TRUE)
  K <- pmin(outer(k, rep(1, length(k))), outer(rep(1, length(k)), k))
  round((A %*% A + A) / (K + 1 - abs(A)), 3)
}

signed_cor <- function(expr, method = "spearman") {
  if (anyNA(expr)) {
    C <- suppressWarnings(stats::cor(t(expr), method = method, use = "pairwise.complete.obs"))
  } else {
    # no missing values: Spearman = Pearson on (average-tie) ranks, one BLAS call
    if (method == "spearman") expr <- t(apply(expr, 1, rank, ties.method = "average"))
    C <- suppressWarnings(stats::cor(t(expr)))
  }
  diag(C) <- 0
  C[is.na(C)] <- 0
  C
}

wto_fast <- function(Data, n = 100, method = "spearman", delta = 0.2) {
  Data <- as.matrix(Data)
  w0 <- wto_matrix(signed_cor(Data, method))
  changed <- matrix(0, nrow(w0), ncol(w0))
  for (i in seq_len(n)) {
    boot <- Data[, sample(ncol(Data), replace = TRUE), drop = FALSE]
    wb <- wto_matrix(signed_cor(boot, method))
    changed <- changed + (wb < w0 - delta) + (wb > w0 + delta)
  }
  ut <- which(upper.tri(w0), arr.ind = TRUE)
  pval <- changed[ut] / n
  data.table::data.table(
    Node.1 = rownames(w0)[ut[, 1]], Node.2 = colnames(w0)[ut[, 2]],
    wTO = w0[ut], pval = pval, pval.adj = stats::p.adjust(pval, method = "BH")
  )
}

# Build (or load) the cached wTO networks used by the slides and the hands-on:
# one network per group (controls, bipolar patients off lithium, bipolar
# patients on lithium) on all selected genes. The cache records the gene list and
# the samples and is rebuilt automatically if either changes.
load_wto_networks <- function(network_matrix, sample_info, n = 100,
                              path = "results/wto_cache_2k.rds", seed = 20261109) {
  genes <- rownames(network_matrix)
  if (file.exists(path)) {
    cached <- readRDS(path)
    if (identical(cached$genes, genes) && identical(cached$n, n) &&
        identical(cached$samples, colnames(network_matrix))) return(cached)
  }
  message("Computing wTO networks for ", length(genes), " genes (one-time; ~15 min)...")
  stopifnot(identical(sample_info$sample, colnames(network_matrix)))
  cols <- function(g) sample_info$sample[sample_info$treatment == g]
  set.seed(seed)
  out <- list(
    genes = genes, n = n, samples = colnames(network_matrix),
    control    = wto_fast(network_matrix[, cols("Control"), drop = FALSE], n = n),
    no_lithium = wto_fast(network_matrix[, cols("no Lithium"), drop = FALSE], n = n),
    lithium    = wto_fast(network_matrix[, cols("Lithium"), drop = FALSE], n = n)
  )
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, path)
  out
}
