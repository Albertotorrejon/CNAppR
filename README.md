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

---

## Workflow 1 — From pre-segmented data

Works with output from ASCAT, CNVkit, SNP arrays, aCGH, or any tool that
produces a standard segment table.

```r
library(CNAppR)

data       <- read_cna_file("segments.tsv")
sample_ids <- unique(data$ID)
reseg_list <- lapply(sample_ids, function(id) resegment_sample(data, sample_id = id))
names(reseg_list) <- sample_ids

scores <- calculate_cna_scores(reseg_list)
```

### Required input columns

| Column | Type | Description |
|---|---|---|
| `ID` | character | Sample identifier |
| `chr` | integer | Chromosome (1–22) |
| `loc.start` | integer | Segment start (bp) |
| `loc.end` | integer | Segment end (bp) |
| `seg.mean` | numeric | Log2 copy number ratio |

### Multi-platform cohorts

When samples come from different assay types (WES, WGS, SNP array), use
`harmonize_segments()` before scoring to normalise the schema and add a
`source` column:

```r
all_segs <- harmonize_segments(
  c(wes_segs, wgs_segs),
  source = c(rep("wes", length(wes_segs)), rep("wgs", length(wgs_segs)))
)
scores <- calculate_cna_scores(split(all_segs, all_segs$ID))
```

---

## CNA scores

| Score | Description |
|---|---|
| **FCS** | Focal CNA Score — sub-arm alterations |
| **BCS** | Broad CNA Score — arm and chromosomal alterations |
| **GCS** | Global CNA Score — combined FCS + BCS |

---

## Visualisation & clinical association

```r
# Genome-wide frequency plot
plot_cn_frequency(reseg_list)

# Per-sample segmentation QC
plot_segmentation(data[data$ID == "sample_1", ], reseg_list[["sample_1"]], "sample_1")

# Gene annotation
annotated <- annotate_cna_genes(reseg_list[["sample_1"]], genome_build = "hg19")

# Association with clinical variables
test_clinical_association(scores, clinical_data = clinical)
```

For a pathway-level view of CNA burden, `run_gsea()` translates CNA segments
into a gene rank vector and runs pre-ranked GSEA against MSigDB (requires
`fgsea` and `msigdbr`):

```r
BiocManager::install("fgsea"); install.packages("msigdbr")

ranks    <- build_gene_ranks(annotated, score_col = "seg.mean")
gsea_res <- run_gsea(ranks, collections = c("H", "C2"))
gsea_res$plot        # NES barplot
gsea_res$significant # FDR < 0.05
```

---

## Optional — Survival analysis

Requires `survival` and `survminer`.

```r
install.packages(c("survival", "survminer"))

res <- run_survival_analysis(scores, survival_data = surv_df, score_col = "GCS")
res$plot    # Kaplan-Meier curve
res$pvalue  # log-rank p-value
```

---

## Workflow 2 — From BAM file

For samples where you have aligned reads and want to go directly to CNA scores.
Requires Bioconductor packages for read counting and GC correction.

```r
BiocManager::install(c("Rsamtools", "GenomicRanges", "DNAcopy"))

# WES GC correction — large download (~800 MB each)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
```

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

result$scores   # FCS, BCS, GCS per sample
result$segments # segmented copy number profile
```

---

## Citation

> Franch-Exposito S, et al. **CNApp, a tool for the quantification of copy number alterations and integrative analysis revealing clinical implications.** *eLife* 2020;9:e50267. https://doi.org/10.7554/eLife.50267

---

MIT © 2026 Alberto Torrejón Aquino · UOC
