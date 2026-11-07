# One-time build of results/wto_cache_2k.rds (run from the workshop folder):
#   Rscript R/prepare_gse124326.R   # if data/gse124326_bd_lithium.rds is missing
#   Rscript R/build_wto_cache.R
# Uses the same preprocessing and gene selection as the slides and the hands-on
# notebook: plate/sex regressed out, genes with mean log-CPM >= 2, then the top 2,000 by variance (338 samples).
source("R/wto_fast.R")
course_data <- readRDS("data/gse124326_bd_lithium.rds")
expression_matrix <- course_data$expression
sample_info <- course_data$sample_info
treatment_design <- model.matrix(~ factor(treatment, levels = c("Control", "no Lithium", "Lithium")),
                                 data = sample_info)
adjusted_matrix <- limma::removeBatchEffect(
  expression_matrix, batch = sample_info$plate,
  covariates = model.matrix(~ sex, data = sample_info)[, -1, drop = FALSE],
  design = treatment_design)
expressed <- rowMeans(expression_matrix) >= 2   # drop low, noisy genes first
gene_variance <- apply(adjusted_matrix[expressed, ], 1, var)
selected_genes <- names(sort(gene_variance, decreasing = TRUE))[1:2000]
network_matrix <- adjusted_matrix[selected_genes, , drop = FALSE]
t0 <- Sys.time()
nets <- load_wto_networks(network_matrix, sample_info, n = 100)
cat("Done in", format(Sys.time() - t0), "\n")
print(sapply(nets[c("control", "no_lithium", "lithium")],
             function(d) c(pairs = nrow(d), stable = sum(d$pval.adj < 0.05))))
