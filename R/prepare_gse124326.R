# Build the workshop input from GEO GSE124326 (Krebs et al. 2020, whole blood,
# bipolar disorder; "reveals strong lithium effect").
#   Rscript R/prepare_gse124326.R
# Output: data/gse124326_bd_lithium.rds, a list with
#   expression  genes x samples, log2-CPM (HGNC symbols, filtered)
#   sample_info one row per sample: treatment ("Control" / "no Lithium" /
#               "Lithium"; the last two are bipolar patients), diagnosis, sex,
#               sequencing plate, assessment group, estimated blood-cell fractions
#   provenance  source files, md5 sums, filters, date
# Samples: everyone the authors kept in their final analysis and assessed in
# group A (all patients are group A; controls from group B are left out so that
# diagnosis is not confounded with assessment group). Patients: bipolar I only
# (the 13 bipolar II patients were all off lithium).

suppressPackageStartupMessages(library(AnnotationDbi))

geo <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124326"
files <- c(counts = "GSE124326_count_matrix.txt.gz", series = "GSE124326_series_matrix.txt.gz")
urls  <- c(counts = file.path(geo, "suppl", files[["counts"]]),
           series = file.path(geo, "matrix", files[["series"]]))
raw_dir <- "data/raw"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(raw_dir, files)
names(paths) <- names(files)
for (k in names(paths)) if (!file.exists(paths[[k]])) {
  options(timeout = 600)
  download.file(urls[[k]], paths[[k]], mode = "wb")
}

# ---- sample metadata from the series matrix ------------------------------
lines <- readLines(gzfile(paths[["series"]]))
split_row <- function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1])
titles <- split_row(grep("^!Sample_title", lines, value = TRUE))
gsm    <- split_row(grep("^!Sample_geo_accession", lines, value = TRUE))
cells  <- do.call(rbind, lapply(grep("^!Sample_characteristics_ch1", lines, value = TRUE), split_row))
# characteristics are not always in the same row order, so look each key up per sample
field <- function(key) apply(cells, 2, function(col) {
  hit <- col[startsWith(col, paste0(key, ":"))]
  if (length(hit)) sub("^[^:]*: ", "", hit[1]) else NA_character_
})
meta <- data.frame(
  sample           = titles,
  geo_accession    = gsm,
  diagnosis        = field("bipolar disorder diagnosis"),
  lithium          = field("lithium use (non-user=0, user = 1)"),
  sex              = field("Sex"),
  age              = suppressWarnings(as.numeric(field("age"))),
  plate            = paste0("plate", field("sequencing plate")),
  assessment_group = field("assessment group"),
  neutrophils      = as.numeric(field("neutrophils")),
  monocytes        = as.numeric(field("monocytes")),
  t_cd4_naive      = as.numeric(field("t.cells.cd4.naive")),
  included         = field("included in final analysis"),
  stringsAsFactors = FALSE
)
meta <- subset(meta, included == "TRUE" & assessment_group == "A" & lithium %in% c("0", "1") &
                 diagnosis %in% c("Control", "BP1"))
meta$treatment <- ifelse(meta$diagnosis == "Control", "Control",
                         ifelse(meta$lithium == "1", "Lithium", "no Lithium"))
meta$included <- meta$lithium <- NULL

# ---- counts -> filtered log2-CPM with HGNC symbols -------------------------
counts <- read.delim(gzfile(paths[["counts"]]), row.names = 1, check.names = FALSE)
colnames(counts) <- sub("\\.counts$", "", colnames(counts))
counts <- as.matrix(counts[!grepl("^__", rownames(counts)), meta$sample])
lib_size <- colSums(counts)
cpm <- t(t(counts) / lib_size) * 1e6
min_n <- ceiling(0.5 * nrow(meta))
counts <- counts[rowSums(cpm >= 1) >= min_n, ]          # expressed in at least half the samples

ensembl <- sub("\\..*$", "", rownames(counts))
symbol <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ensembl,
                                                 column = "SYMBOL", keytype = "ENSEMBL"))
logcpm <- log2(t((t(counts) + 0.5) / (lib_size + 1)) * 1e6)  # voom's log-CPM, original library sizes
keep <- !is.na(symbol)
logcpm <- logcpm[keep, ]; symbol <- symbol[keep]
# one row per symbol: keep the Ensembl gene with the highest mean expression
ord <- order(symbol, -rowMeans(logcpm))
logcpm <- logcpm[ord, ]; symbol <- symbol[ord]
first <- !duplicated(symbol)
expression <- round(logcpm[first, ], 3)
rownames(expression) <- symbol[first]
expression <- expression[order(rownames(expression)), ]

out <- list(
  expression = expression,
  sample_info = meta[match(colnames(expression), meta$sample), ],
  provenance = list(
    geo = "GSE124326", reference = "Krebs et al. (2020) Psychol Med 50:2575-2586",
    files = basename(paths), md5 = unname(tools::md5sum(paths)),
    samples = "bipolar I patients and assessment-group-A controls included in the authors' final analysis",
    filter = sprintf("CPM >= 1 in >= %d samples; symbols via org.Hs.eg.db; duplicates -> highest mean", min_n),
    transform = "log2((count + 0.5) / (library size + 1) * 1e6)",
    created = format(Sys.time(), "%Y-%m-%d")
  )
)
stopifnot(identical(out$sample_info$sample, colnames(out$expression)))
dir.create("data", showWarnings = FALSE)
saveRDS(out, "data/gse124326_bd_lithium.rds", compress = "xz")
cat(sprintf("%d genes x %d samples\n", nrow(expression), ncol(expression)))
print(table(out$sample_info$treatment))
print(table(out$sample_info$treatment, out$sample_info$diagnosis))
