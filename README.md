# CNAppR

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg)
![Platform](https://img.shields.io/badge/platform-WES%20%7C%20WGS%20%7C%20Array-brightgreen)

R package for somatic Copy Number Alteration (CNA) analysis, developed as part of a Master's thesis in Bioinformatics (UOC). Extends the original [CNApp](https://github.com/ait5/CNApp) framework with native BAM processing, Nanopore support, and multi-platform cohort harmonisation.

---

## Installation

```r
devtools::install_github("Albertotorrejon/CNAppR")
```

For BAM processing (WES/WGS from aligned reads):

```r
BiocManager::install(c("Rsamtools", "GenomicRanges", "DNAcopy"))

# WES GC correction — large download (~800 MB each, install only the build you need)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

---

## Workflow 1 — From BAM file

```r
library(CNAppR)

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

result$scores   # data.frame with FCS, BCS, GCS per sample
result$segments # segmented copy number profile
```

---

## Workflow 2 — From pre-segmented data

Works with output from ASCAT, CNVkit, SNP arrays, aCGH, or any tool that
produces a standard segment table.

```r
# Load and re-segment
data       <- read_cna_file("segments.tsv")
sample_ids <- unique(data$ID)
reseg_list <- lapply(sample_ids, function(id) resegment_sample(data, sample_id = id))
names(reseg_list) <- sample_ids

# Score
scores <- calculate_cna_scores(reseg_list)

# Visualise
plot_cn_frequency(reseg_list)
plot_segmentation(data[data$ID == "sample_1", ], reseg_list[["sample_1"]], "sample_1")

# Gene annotation and clinical association
annotated <- annotate_cna_genes(reseg_list[["sample_1"]], genome_build = "hg19")
test_clinical_association(scores, clinical_data = clinical)
```

### Required input columns

| Column | Type | Description |
|---|---|---|
| `ID` | character | Sample identifier |
| `chr` | integer | Chromosome (1–22) |
| `loc.start` | integer | Segment start (bp) |
| `loc.end` | integer | Segment end (bp) |
| `seg.mean` | numeric | Log2 copy number ratio |

---

## CNA scores

| Score | Description |
|---|---|
| **FCS** | Focal CNA Score — sub-arm alterations |
| **BCS** | Broad CNA Score — arm and chromosomal alterations |
| **GCS** | Global CNA Score — combined FCS + BCS |

---

## Optional extensions

### Pathway enrichment (GSEA)

Requires `fgsea` and `msigdbr`.

```r
BiocManager::install("fgsea"); install.packages("msigdbr")

ranks    <- build_gene_ranks(annotated, score_col = "seg.mean")
gsea_res <- run_gsea(ranks, collections = c("H", "C2"))
gsea_res$plot            # NES barplot
gsea_res$significant     # significant pathways (FDR < 0.05)
```

### Survival analysis

Requires `survival` and `survminer`.

```r
install.packages(c("survival", "survminer"))

res <- run_survival_analysis(scores, survival_data = surv_df, score_col = "GCS")
res$plot    # Kaplan-Meier curve
res$pvalue  # log-rank p-value
```

---

## Citation

> Franch-Exposito S, et al. **CNApp, a tool for the quantification of copy number alterations and integrative analysis revealing clinical implications.** *eLife* 2020;9:e50267. https://doi.org/10.7554/eLife.50267

---

MIT © 2026 Alberto Torrejón Aquino · UOC
