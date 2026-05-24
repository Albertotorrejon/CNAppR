# CNAppR

<!-- badges -->
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg)
![Platform](https://img.shields.io/badge/platform-WES%20%7C%20WGS%20%7C%20Array-brightgreen)

**CNAppR** is an R package for the analysis of somatic Copy Number Alterations (CNAs) in cancer genomics. It extends the original [CNApp](https://github.com/ait5/CNApp) framework with a **native BAM processing pipeline**, **pathway enrichment analysis (GSEA)**, and **survival analysis** — enabling end-to-end analysis directly from sequencing files with no external tools required.

Designed for translational research, CNAppR supports whole-exome sequencing (WES), whole-genome sequencing (WGS, including low-pass Nanopore), and pre-segmented data from SNP arrays or aCGH.

---

## Installation

```r
# Requires devtools
devtools::install_github("Albertotorrejon/CNAppR")
```

### Optional dependencies

Install only the modules you need:

```r
# BAM pipeline (WES/WGS from aligned reads)
BiocManager::install(c("Rsamtools", "GenomicRanges", "DNAcopy"))

# WES GC correction (BSgenome-based)
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19",
                       "BSgenome.Hsapiens.UCSC.hg38", "Biostrings"))

# Faster WES binning (replaces internal loops with BED-standard operations)
install.packages("valr")

# Pathway enrichment (GSEA)
BiocManager::install("fgsea")
install.packages("msigdbr")

# Survival analysis
install.packages(c("survival", "survminer"))
```

---

## Workflows

CNAppR supports two entry points depending on your data:

| Workflow | Input | Key function |
|---|---|---|
| **BAM pipeline** | `.bam` file (WES or WGS) | `run_bam_pipeline()` |
| **Pre-segmented** | Segment table (`.txt`/`.tsv`) | `read_cna_file()` |

---

## Workflow 1 — BAM pipeline

For samples where you have a BAM file and want to go directly to CNA scores.

### WES (whole-exome sequencing)

```r
library(CNAppR)

result <- run_bam_pipeline(
  bam_path        = "tumor.bam",
  sequencing_type = "illumina",
  experiment_type = "wes",
  targets_bed     = "capture_kit.bed",   # required for WES
  genome_build    = "hg19",
  gc_correction   = TRUE
)
```

### WGS (whole-genome sequencing, including low-pass Nanopore)

```r
result <- run_bam_pipeline(
  bam_path        = "cfDNA_sample.bam",
  sequencing_type = "nanopore",          # or "illumina"
  experiment_type = "wgs",
  bin_size        = 500e3,               # 500 kb bins recommended
  genome_build    = "hg38",
  gc_correction   = TRUE
)
```

### Step-by-step pipeline

`run_bam_pipeline()` is a wrapper around three composable functions:

```r
# 1. Read BAM, count reads per bin, apply GC correction
qc <- read_bam_qc(
  "tumor.bam",
  sequencing_type = "illumina",
  experiment_type = "wes",
  targets_bed     = "capture_kit.bed"
)

# 2. CBS segmentation
seg <- run_cbs_segmentation(qc, sample_id = "sample_01")

# 3. Re-segment and classify CNAs
reseg <- resegment_sample(seg, sample_id = "sample_01")
```

### WES binning with valr

When the `valr` package is installed, `read_bam_qc()` automatically uses
`valr::bed_makewindows()` and `valr::bed_complement()` for bin generation
(CNVkit-style target splitting and antitarget computation). If valr is
not available, the original loop-based implementation is used as fallback.

```r
install.packages("valr")   # recommended for WES
```

### Panel of Normals (optional)

```r
pon <- build_pon(
  bam_paths       = c("normal1.bam", "normal2.bam", "normal3.bam"),
  sequencing_type = "illumina",
  experiment_type = "wes",
  targets_bed     = "capture_kit.bed"
)

result <- run_bam_pipeline("tumor.bam", pon = pon,
                            experiment_type = "wes",
                            targets_bed     = "capture_kit.bed")
```

---

## Workflow 2 — Pre-segmented input

For samples from SNP arrays, aCGH, or any upstream segmentation tool.

```r
library(CNAppR)

# Load data (tab-separated segment table)
data <- read_cna_file(system.file("models", "demo.txt", package = "CNAppR"))

# Re-segment and score all samples
sample_ids  <- unique(data$ID)
reseg_list  <- lapply(sample_ids, function(id) resegment_sample(data, sample_id = id))
names(reseg_list) <- sample_ids

scores <- calculate_cna_scores(reseg_list)

# Annotate CNA segments with overlapping genes
annotated <- annotate_cna_genes(reseg_list[["sample_1"]], genome_build = "hg19")

# Visualise
plot_segmentation(
  original_data    = data[data$ID == "sample_1", ],
  resegmented_data = reseg_list[["sample_1"]],
  sample_id        = "sample_1"
)

plot_cn_frequency(reseg_list)

# Test association with clinical variables
results <- test_clinical_association(scores, clinical_data = clinical)
```

---

## Pathway Enrichment Analysis (GSEA)

Requires `fgsea` (Bioconductor) and `msigdbr` (CRAN).

```r
# 1. Annotate segments with gene symbols
reseg <- resegment_sample(data, sample_id = "sample_1")
ann   <- annotate_cna_genes(reseg, genome_build = "hg19")

# 2. Build a gene rank vector (positive = amplified, negative = deleted)
ranks <- build_gene_ranks(ann, score_col = "seg.mean")

# 3. Run GSEA (Hallmark + curated gene sets)
gsea_res <- run_gsea(
  ranks,
  collections = c("H", "C2"),
  alpha       = 0.25,
  top_n       = 20
)

gsea_res$plot           # NES barplot
head(gsea_res$significant)  # significant pathways

# 4. Stratified GSEA by clinical group
ranks_list <- lapply(sample_ids, function(id) {
  build_gene_ranks(annotate_cna_genes(resegment_sample(data, id)))
})
names(ranks_list) <- sample_ids

groups   <- c(sample_1 = "MSI", sample_2 = "MSS")   # example
gsea_str <- run_gsea(ranks_list, clinical_groups = groups, collections = "H")
```

---

## Survival Analysis

Requires `survival` and `survminer`.

```r
# Median split on GCS
surv_data <- data.frame(
  time   = c(24, 36, 12, 60, 48),
  event  = c(1,  0,  1,  0,  1),
  row.names = sample_ids
)

res <- run_survival_analysis(
  scores       = scores,
  survival_data = surv_data,
  score_var    = "GCS",      # or "FCS" / "BCS"
  n_groups     = 2,          # 2 = median split; 4 = quartiles
  time_col     = "time",
  event_col    = "event"
)

res$plot     # Kaplan-Meier curve
res$pvalue   # log-rank p-value

# Stratified by MSI status
res_msi <- run_survival_analysis(scores, surv_data,
                                  score_var    = "FCS",
                                  n_groups     = 4,
                                  clinical_var = "msi_status")
```

---

## Clinical Annotation

```r
# Merge external annotation file with CNA data
data_ann <- prepare_annotation_data(data, "annotations.tsv")

# Extract and classify clinical variables (numeric / categorical)
var_prep <- prepare_clinical_variables(data_ann, exclude_cols = c("BAF", "purity"))

# Test CNA score associations with all clinical variables
assoc <- test_clinical_association(scores, var_prep$mat_variables)
assoc$pval_nonparametric   # Spearman / Kruskal-Wallis
assoc$pval_parametric      # Pearson / ANOVA
```

---

## Input formats

### BAM pipeline

| Parameter | Required | Description |
|---|---|---|
| `.bam` file | Yes | Aligned reads, indexed (`.bai` must exist alongside) |
| `.bed` file | WES only | Capture kit target regions |

### Pre-segmented table

Minimum required columns:

| Column | Type | Description |
|---|---|---|
| `ID` | character | Sample identifier |
| `chr` | integer | Chromosome (1–22) |
| `loc.start` | integer | Segment start (bp) |
| `loc.end` | integer | Segment end (bp) |
| `seg.mean` | numeric | Log2 copy number ratio |

Optional: `BAF` (B-allele frequency), `purity` (tumour purity 0–1).

---

## CNA score definitions

| Score | Description |
|---|---|
| **FCS** | Focal CNA Score — weighted sum of sub-arm alteration scores |
| **BCS** | Broad CNA Score — count of arm- and chromosomal-level alterations |
| **GCS** | Global CNA Score — normalised combination of FCS and BCS |

---

## CNA classification

Each segment is classified based on its length relative to the chromosome:

| Level | Criterion | Description |
|---|---|---|
| **chromosomal** | > 90% chromosome length | Whole-chromosome alteration |
| **arm** | > 50% arm length | Arm-level alteration (p or q) |
| **focal** | < arm threshold | Sub-arm alteration |
| **diploid** | No alteration | Segment within normal range |

Alterations are typed as **Gain**, **Loss** or **CN-LOH** (requires BAF column), with intensity grades **Low**, **Medium** or **High**.

---

## Key parameters — `resegment_sample()`

| Parameter | Default | Description |
|---|---|---|
| `min_length` | 100,000 bp | Minimum segment length to retain |
| `dev_btw_segs` | 0.16 | Max log2 difference to merge adjacent segments |
| `dev_tozero` | 0.16 | Segments within this range of 0 are set to neutral |
| `chrom_percent` | 0.90 | Coverage threshold for chromosomal classification |
| `arm_percent` | 0.50 | Coverage threshold for arm-level classification |

---

## Citation

If you use CNAppR, please cite the original CNApp publication:

> Franch-Exposito S, Bassaganyas L, Vila-Casadesus M, Hernandez-Illan E, Esteban-Fabro R, Diaz-Gay M, Lozano JJ, Castells A, Llovet JM, Castellvi-Bel S, Camps J.
> **CNApp, a tool for the quantification of copy number alterations and integrative analysis revealing clinical implications.**
> *eLife* 2020;9:e50267. https://doi.org/10.7554/eLife.50267

---

## License

MIT © 2026 Alberto Torrejón Aquino · UOC
