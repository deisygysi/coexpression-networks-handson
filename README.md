# Co-expression Networks · Hands-on

Companion notebook for **Tutorial VT1 · Co-expression Networks: Construction and Analysis**,
ISCB-Latin America SoIBio SPBBC Conference on Bioinformatics 2026.
Dr. Deisy Morselli Gysi · Federal University of Paraná (UFPR), Brazil.

You will build and compare gene co-expression networks on whole-blood RNA-seq
from people with bipolar I disorder, on and off lithium, and from controls
(GEO [GSE124326](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124326);
Krebs *et al.*, 2020, *Psychol Med* 50:2575–2586).

| Step | Method |
|:--|:--|
| Data and differential expression | limma (lithium and diagnosis contrasts) |
| Networks | WGCNA, ARACNE, wTO |
| Comparing networks | wTO consensus, CoDiNA |
| Function | ORA and GSEA with clusterProfiler |

## Get started

1. Clone or download this repository and open the folder in RStudio or Positron.
2. Restore the pinned R packages once:

   ```r
   install.packages("renv")
   renv::restore()
   ```

3. Open `coexpression_handson_iscb_la_2026.qmd` and run the chunks in order.
   Each **Exercise** has a chunk with `___` blanks to fill and a collapsed
   **Show solution**.

To just read the notebook with all outputs, open
`coexpression_handson_iscb_la_2026.html` in a browser.

## What is in this folder

| Path | Contents |
|:--|:--|
| `coexpression_handson_iscb_la_2026.qmd` | The notebook you work in |
| `coexpression_handson_iscb_la_2026.html` | The same notebook, rendered |
| `data/gse124326_bd_lithium.rds` | Prepared input: log2-CPM matrix and sample information (338 samples) |
| `results/wto_cache_2k.rds` | Pre-computed wTO networks (2,000 genes, 3 groups), so they load in seconds |
| `R/wto_fast.R` | Fast wTO, identical output to `wTO::wTO.fast()` |
| `R/prepare_gse124326.R` | Rebuilds the input from GEO (optional) |
| `R/build_wto_cache.R` | Rebuilds the wTO cache, about 15 min (optional) |
| `renv.lock`, `renv/`, `.Rprofile` | Pinned package versions |

## Rebuild from scratch (optional)

```r
source("R/prepare_gse124326.R")  # downloads the GEO files into data/raw/
source("R/build_wto_cache.R")    # recomputes results/wto_cache_2k.rds
```

Then render the notebook from a terminal:

```sh
quarto render coexpression_handson_iscb_la_2026.qmd
```
