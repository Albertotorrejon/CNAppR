# CNAppR

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg)](https://cran.r-project.org/)
[![Platform](https://img.shields.io/badge/platform-WES%20%7C%20WGS%20%7C%20Array-brightgreen)]()

CNAppR is an R package for quantifying somatic copy number alterations (CNAs) in cancer. It extends the [CNApp](https://github.com/ait5/CNApp) framework with a native BAM processing pipeline, Nanopore support, GSEA-based pathway enrichment, and survival analysis. The package works with pre-segmented data from any upstream tool (ASCAT, CNVkit, SNP arrays) as well as directly from aligned BAM files.

## Installation

```r
devtools::install_github("Albertotorrejon/CNAppR")
```

For BAM processing (optional — install only what you need):

```r
# Read counting and CBS segmentation
BiocManager::install(c("Rsamtools", "GenomicRanges", "DNAcopy"))

# GC correction for WES (~800 MB each, install only the build you need)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

## Quick run

```r
library(CNAppR)

# Load a pre-segmented CNA table
data <- read_cna_file(system.file("models", "LIHC_seg.txt", package = "CNAppR"))
head(data)
#>          ID chr loc.start  loc.end   seg.mean
#> 1 TCGA-BC-A... 1    554268 11063088  0.0412
#> 2 TCGA-BC-A... 1  11063088 12263592 -0.9871

# Re-segment and classify all samples
sample_ids <- unique(data$ID)
reseg_list <- lapply(sample_ids, function(id) resegment_sample(data, sample_id = id))
names(reseg_list) <- sample_ids

# Compute FCS, BCS and GCS scores
scores <- calculate_cna_scores(reseg_list)
head(scores)
#>            ID      FCS      BCS      GCS
#> 1 TCGA-BC-A...  0.1823   4.0000  0.3241
#> 2 TCGA-DD-A...  0.0541   2.0000  0.1042
```

Visualise the cohort-level CNA frequency:

```r
plot_cn_frequency(reseg_list)
```

## BAM pipeline

For samples where you have aligned reads:

```r
# WES — Illumina
result <- run_bam_pipeline(
  bam_path        = "tumor.bam",
  sequencing_type = "illumina",
  targets_bed     = "capture_kit.bed",
  genome_build    = "hg19"
)

# WGS — Nanopore low-pass cfDNA
result <- run_bam_pipeline(
  bam_path        = "cfDNA_nanopore.bam",
  sample_id       = "sample_01",
  sequencing_type = "nanopore",
  genome_build    = "hg38"
)

result$scores    # data.frame with FCS, BCS, GCS
result$segments  # segmented copy number profile
```

## Pathway enrichment (optional)

Requires `fgsea` and `msigdbr`:

```r
BiocManager::install("fgsea")
install.packages("msigdbr")

annotated <- annotate_cna_genes(reseg_list[["sample_1"]], genome_build = "hg19")
ranks     <- build_gene_ranks(annotated, score_col = "seg.mean")

gsea_res  <- run_gsea(ranks, collections = c("H", "C2"))
gsea_res$plot         # NES barplot
gsea_res$significant  # pathways with FDR < 0.05
```

## Survival analysis (optional)

Requires `survival` and `survminer`:

```r
install.packages(c("survival", "survminer"))

res <- run_survival_analysis(scores, survival_data = surv_df, score_col = "GCS")
res$plot    # Kaplan-Meier curve
res$pvalue  # log-rank p-value
```

## Input format

Minimum required columns for pre-segmented data:

| Column | Type | Description |
|---|---|---|
| `ID` | character | Sample identifier |
| `chr` | integer | Chromosome (1–22) |
| `loc.start` | integer | Segment start (bp) |
| `loc.end` | integer | Segment end (bp) |
| `seg.mean` | numeric | Log2 copy number ratio |

## CNA scores

| Score | Description |
|---|---|
| **FCS** | Focal CNA Score — weighted sub-arm alteration burden |
| **BCS** | Broad CNA Score — count of arm and chromosomal alterations |
| **GCS** | Global CNA Score — normalised combination of FCS and BCS |

## Citation

If you use CNAppR, please cite the original CNApp paper:

> Franch-Exposito S, Bassaganyas L, Vila-Casadesus M, Hernandez-Illan E, Esteban-Fabro R, Diaz-Gay M, Lozano JJ, Castells A, Llovet JM, Castellvi-Bel S, Camps J.
> **CNApp, a tool for the quantification of copy number alterations and integrative analysis revealing clinical implications.**
> *eLife* 2020;9:e50267. https://doi.org/10.7554/eLife.50267

---

MIT © 2026 Alberto Torrejón Aquino · UOC
