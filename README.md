# CNAppR

**CNAppR** is an R package for integrative analysis of somatic Copy Number Alterations (CNAs) in cancer genomics. It implements the analytical core of the original [CNApp](https://github.com/ait5/CNApp) tool as a standalone R library.

## Features

- **Re-segmentation**: iterative fusion of similar adjacent segments with optional purity correction
- **CNA scoring**: computation of Focal (FCS), Broad (BCS) and Global (GCS) CNA scores
- **Classification**: chromosomal-, arm- and focal-level alteration detection
- **Visualisation**: genome-wide segmentation plots and copy number frequency heatmaps
- **Clinical association**: parametric and non-parametric tests linking CNA scores to clinical variables
- **ML classification**: random forest-based sample classification from CNA scores

## Installation

```r
# From GitHub (requires devtools)
devtools::install_github("Albertotorrejon/CNAppR")
```

## Quick start

```r
library(CNAppR)

# 1. Load and validate your CNA data (TSV with columns: ID, chr, loc.start, loc.end, seg.mean)
data <- read_cna_file("your_data.tsv")
validate_cna_data(data)

# 2. Re-segment and classify alterations for each sample
sample_ids <- unique(data$ID)
reseg_list <- lapply(sample_ids, function(id) resegment_sample(data, sample_id = id))
names(reseg_list) <- sample_ids

# 3. Calculate CNA scores
scores <- calculate_cna_scores(reseg_list)
head(scores)
#>            FCS BCS        GCS
#> sample_1   12   4  1.2043...
#> sample_2    3  11 -1.2043...

# 4. Annotate CNA segments with overlapping genes
annotated_list <- lapply(sample_ids, function(id) {
  annotate_cna_genes(reseg_list[[id]], genome_build = "hg19")
})
names(annotated_list) <- sample_ids

# 5. Visualise
plot_segmentation(
  original_data    = data[data$ID == "sample_1", ],
  resegmented_data = reseg_list[["sample_1"]],
  sample_id        = "sample_1"
)

plot_cn_frequency(reseg_list)

# 6. Test association with clinical variables (optional)
# clinical must be a data.frame with rownames = sample IDs
results <- test_clinical_association(scores, clinical_data = clinical)
```

## Input format

The main data file must be a tab-separated (or CSV) file with at least these columns:

| Column | Type | Description |
|---|---|---|
| `ID` | character | Sample identifier |
| `chr` | integer | Chromosome (1–22) |
| `loc.start` | integer | Segment start position (bp) |
| `loc.end` | integer | Segment end position (bp) |
| `seg.mean` | numeric | Log2 copy number ratio |

Optional columns: `BAF` (B-allele frequency), `purity` (tumour purity 0–1).

## CNA classification

Each segment is classified into one of four levels based on its length relative to the chromosome:

| Level | Criterion | Description |
|---|---|---|
| **chromosomal** | > 90% of chromosome length | Whole-chromosome alteration |
| **arm** | > 50% of a chromosome arm | Arm-level alteration (p or q) |
| **focal** | < arm threshold | Sub-arm alteration, weighted by size |
| **diploid** | No alteration detected | Segment within normal copy number range |

Within each level, alterations are typed as **Gain**, **Loss** or **CN-LOH** (copy-neutral loss of heterozygosity, requires BAF column), and assigned an intensity (**Low**, **Medium**, **High**) based on the seg.mean thresholds.

## Gene annotation

`annotate_cna_genes()` maps each CNA segment to the genes it overlaps using RefSeq (UCSC) reference files bundled with the package for hg19 and hg38:

```r
reseg_annotated <- annotate_cna_genes(
  reseg_list[["sample_1"]],
  genome_build = "hg19",   # or "hg38"
  only_cna     = TRUE      # skip diploid segments
)

# Gene names are returned as a semicolon-separated string in the `genes` column
head(reseg_annotated[!is.na(reseg_annotated$genes),
                     c("chr", "loc.start", "loc.end", "type", "genes")])
```

## CNA score definitions

| Score | Description |
|---|---|
| **FCS** | Focal CNA Score — count of focal alterations |
| **BCS** | Broad CNA Score — count of arm- and chromosomal-level alterations |
| **GCS** | Global CNA Score — normalised sum of FCS and BCS |

## Re-segmentation parameters

The key parameters for `resegment_sample()`:

| Parameter | Default | Description |
|---|---|---|
| `min_length` | 100,000 bp | Minimum segment length to retain |
| `dev_btw_segs` | 0.16 | Max seg.mean difference to merge two segments |
| `dev_tozero` | 0.16 | Segments within this range of 0 are set to neutral |
| `chrom_percent` | 0.9 | Min chromosome coverage for chromosomal classification |
| `arm_percent` | 0.5 | Min arm coverage for arm-level classification |

## Dependencies

Core (required): `ggplot2`, `rlang`, `RColorBrewer`

Suggested: `randomForest`, `caret`, `survival`, `survminer`, `plotly`, `heatmaply`

Bioconductor (optional): `GenomicFeatures`, `GenomicAlignments`, `GenVisR`

## License

MIT
