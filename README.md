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

# Load a pre-segmented CNA table (your own file or the bundled example)
data <- read_cna_file("path/to/your/segments.txt")

# Example with the bundled LIHC dataset (354 samples)
data <- read_cna_file(
  system.file("models", "LIHC_354_cnvsegments_input_scores_nopurity.txt",
              package = "CNAppR")
)

head(data)
#>            ID chr loc.start   loc.end  seg.mean
#> 1 TCGA-2Y-A9GS   1 205389460 205391959 -0.9772
#> 2 TCGA-2Y-A9GS   1 205392408 247813706  0.1764
```

## Re-segmentation and scoring

`resegment_sample()` processes one sample at a time. To run all samples in a dataset, loop over the sample IDs and collect results in a named list:

```r
sample_ids <- unique(data$ID)

seg <- lapply(sample_ids, function(id) {
  tryCatch(
    resegment_sample(data, sample_id = id),
    error = function(e) {
      warning(id, ": ", conditionMessage(e))
      NULL
    }
  )
})
names(seg) <- sample_ids
seg <- Filter(Negate(is.null), seg)  # drop any samples that failed

# Compute FCS, BCS and GCS scores
scores <- calculate_cna_scores(seg)
head(scores)
#>              FCS BCS    GCS
#> TCGA-2Y-A9GS   2   3  0.82
#> TCGA-DD-A1QA   0   1 -0.54
```

Visualise the cohort-level CNA frequency:

```r
plot_cn_frequency(seg)
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

Optional columns preserved throughout the pipeline: `BAF`, `purity`.

## BAM pipeline

For samples where you have aligned reads, use `read_bam_qc()` to count reads per bin and `run_bam_pipeline()` as a single-call convenience wrapper:

```r
# WES — Illumina
seg_wes <- run_bam_pipeline(
  bam_path        = "tumor.bam",
  sample_id       = "sample_01",
  sequencing_type = "illumina",
  targets_bed     = "capture_kit.bed",
  genome_build    = "hg19"
)

# WGS — Nanopore low-pass cfDNA
seg_wgs <- run_bam_pipeline(
  bam_path        = "cfDNA_nanopore.bam",
  sample_id       = "sample_01",
  sequencing_type = "nanopore",
  bin_size        = 1000000L,
  genome_build    = "hg19"
)
```

Both functions return a data frame in the same format as `resegment_sample()`, ready to be passed directly to `calculate_cna_scores()`.

### Panel of Normals (PoN)

For cohort-level normalisation, build a PoN from a set of matched normal BAMs before running the tumour pipeline:

```r
pon <- build_pon(
  bam_paths       = c("normal1.bam", "normal2.bam", "normal3.bam"),
  genome_build    = "hg19",
  experiment_type = "wes"
)

# Pass the PoN to read_bam_qc() or run_bam_pipeline()
seg <- run_bam_pipeline("tumor.bam", sample_id = "sample_01",
                         sequencing_type = "illumina",
                         genome_build = "hg19", pon = pon)
```

### Combining WES and WGS samples

Use `harmonize_segments()` to merge segment tables from different platforms or sequencing types into a single cohort table:

```r
seg_all <- harmonize_segments(
  segments_list = c(seg_wes_list, seg_wgs_list),
  source        = c(rep("wes", length(seg_wes_list)),
                    rep("wgs", length(seg_wgs_list)))
)

scores <- calculate_cna_scores(split(seg_all, seg_all$ID))
```

## Pathway enrichment (optional)

Requires `fgsea` and `msigdbr`:

```r
BiocManager::install("fgsea")
install.packages("msigdbr")

sample_id <- sample_ids[1]
annotated <- annotate_cna_genes(seg[[sample_id]], genome_build = "hg19")
ranks     <- build_gene_ranks(annotated, score_col = "seg.mean")

gsea_res  <- run_gsea(ranks, collections = c("H", "C2"))
gsea_res$plot         # NES barplot
gsea_res$significant  # pathways with FDR < 0.05
```

## Survival analysis (optional)

Requires `survival` and `survminer`:

```r
install.packages(c("survival", "survminer"))

res <- run_survival_analysis(scores, survival_data = surv_df, score_var = "GCS")
res$plot    # Kaplan-Meier curve
res$pvalue  # log-rank p-value
```

`surv_df` must have columns `time` (numeric, months/days) and `event` (0/1), with row names matching the sample IDs in `scores`.

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
