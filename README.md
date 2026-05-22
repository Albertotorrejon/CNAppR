# CNAppR

<!-- badges -->
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg)
![Platform](https://img.shields.io/badge/platform-WES%20%7C%20WGS%20%7C%20Array-brightgreen)

**CNAppR** is an R package for the analysis of somatic Copy Number Alterations (CNAs) in cancer genomics. It extends the original [CNApp](https://github.com/ait5/CNApp) framework with a **native BAM processing pipeline**, enabling end-to-end analysis directly from sequencing files — no FASTA, no external tools required.

Designed for translational research, CNAppR is built around next-generation sequencing data — primarily whole-exome sequencing (WES) and whole-genome sequencing (WGS, including low-pass Nanopore). Pre-segmented data from SNP arrays or aCGH can also be used through the second workflow.

---

## Installation

```r
# Requires devtools
devtools::install_github("Albertotorrejon/CNAppR")
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

result$scores   # FCS, BCS, GCS
result$segments # segmented profile
```

### WGS (whole-genome sequencing, including low-pass Nanopore)

```r
result <- run_bam_pipeline(
  bam_path        = "cfDNA_sample.bam",
  sequencing_type = "nanopore",          # or "illumina"
  experiment_type = "wgs",
  bin_size        = 500e3,               # 500 kb bins recommended for low-pass
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
seg <- run_cbs_segmentation(qc)

# 3. Compute CNA scores
scores <- calculate_cna_scores(seg$segments)
```

### Panel of Normals (optional)

For tumour-only workflows, build a PoN from normal samples to remove systematic biases:

```r
# Build once from matched normal BAMs
pon <- build_pon(
  bam_paths       = c("normal1.bam", "normal2.bam", "normal3.bam"),
  sequencing_type = "illumina",
  experiment_type = "wes",
  targets_bed     = "capture_kit.bed"
)

# Apply at inference time
result <- run_bam_pipeline("tumor.bam", pon = pon, ...)
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

A demo dataset of 160 colorectal cancer samples (TCGA-COAD) is bundled with the package.

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
| **FCS** | Focal CNA Score — count of sub-arm alterations weighted by size |
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
