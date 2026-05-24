# BAM → CBS → CNAppR pipeline
#
# Full workflow: BAM file → QC + read counting → CBS segmentation →
#   CNAppR format → resegment_sample()
#
# Required Bioconductor packages: Rsamtools, GenomicRanges, DNAcopy
# Install with: BiocManager::install(c("Rsamtools", "GenomicRanges", "DNAcopy"))
#
# Illumina (WES/exome): paired-end reads, smaller bins (~500 kb)
# Nanopore (low-pass WGS): single long reads, larger bins (~1 Mb)


# ─── internal helpers ────────────────────────────────────────────────────────

.check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required package not installed: ", pkg,
         "\nInstall with: BiocManager::install('", pkg, "')")
}

.check_bam <- function(bam_path) {
  if (!file.exists(bam_path))
    stop("BAM file not found: ", bam_path)
  bai_path <- paste0(bam_path, ".bai")
  if (!file.exists(bai_path))
    stop("BAM index (.bai) not found.\nRun: samtools index ",
         basename(bam_path))
}

# Loads a WES capture kit BED file and converts it to GRanges
# with chromosome names matching the BAM header.
.load_targets_bed <- function(targets_bed, chr_lens) {
  if (!file.exists(targets_bed))
    stop("Targets BED file not found: ", targets_bed)

  bed_raw <- readLines(targets_bed)
  bed_raw <- bed_raw[!grepl("^(#|track|browser)", bed_raw)]
  bed <- read.table(text = paste(bed_raw, collapse = "\n"), header = FALSE,
                    sep = "\t", stringsAsFactors = FALSE, fill = TRUE)[, 1:3]
  colnames(bed) <- c("chr", "start", "end")
  bed$start <- suppressWarnings(as.integer(bed$start))
  bed$end   <- suppressWarnings(as.integer(bed$end))
  bed <- bed[!is.na(bed$start) & !is.na(bed$end) & bed$end > bed$start, ]
  if (nrow(bed) == 0)
    stop("BED file contains no valid regions: ", targets_bed)

  bam_chrs   <- names(chr_lens)
  has_prefix <- any(grepl("^chr", bam_chrs))
  bed$chr <- if (has_prefix) {
    ifelse(grepl("^chr", bed$chr), bed$chr, paste0("chr", bed$chr))
  } else {
    gsub("^chr", "", bed$chr)
  }

  bed <- bed[bed$chr %in% bam_chrs, ]
  if (nrow(bed) == 0)
    stop("No BED targets match the BAM chromosomes.\n",
         "BAM chromosomes: ", paste(head(bam_chrs, 5), collapse = ", "), "...")

  # BED is 0-based half-open → GRanges is 1-based closed
  GenomicRanges::GRanges(
    seqnames = bed$chr,
    ranges   = IRanges::IRanges(start = bed$start + 1L, end = bed$end)
  )
}

# Loads the pre-computed bin annotation table (GC + mappability)
# and joins it with the BAM bins by chromosomal position.
.load_bin_annotations <- function(bins, genome_build, bin_size) {
  bin_kb  <- round(bin_size / 1000)
  fname   <- sprintf("bins_%s_%dkb.rds", genome_build, bin_kb)
  fpath   <- system.file("extdata", fname, package = "CNAppR")

  if (!nzchar(fpath) || !file.exists(fpath))
    stop("Annotation table not found: ", fname,
         "\nAvailable bin sizes: 500kb, 1000kb, 2000kb, 10000kb",
         "\nAvailable genomes: hg19, hg38",
         "\nOr generate the table with inst/scripts/generate_bin_annotations.R")

  ann <- readRDS(fpath)

  # Join key: chr + loc.start
  bins_df <- data.frame(
    chr       = as.character(GenomicRanges::seqnames(bins)),
    loc.start = GenomicRanges::start(bins),
    stringsAsFactors = FALSE
  )
  merged <- merge(bins_df, ann[, c("chr", "loc.start", "gc", "mappability")],
                  by = c("chr", "loc.start"), all.x = TRUE, sort = FALSE)

  # Restore original order (merge may reorder)
  idx <- match(paste(bins_df$chr, bins_df$loc.start),
               paste(merged$chr, merged$loc.start))
  merged <- merged[idx, ]

  list(gc = merged$gc, mappability = merged$mappability)
}

# HMMcopy-style LOESS GC correction: outlier filtering + reads ~ GC regression
.correct_gc_loess <- function(reads, gc) {
  # Step 1: base filter matching HMMcopy (valid GC and bin with reads)
  base_valid <- !is.na(gc) & gc > 0.2 & gc < 0.8 & reads > 0
  if (sum(base_valid) < 10) {
    warning("Too few valid bins for GC correction, skipping.")
    return(reads)
  }

  # Step 2: remove read outliers by IQR (matching HMMcopy, threshold 2.5)
  r      <- reads[base_valid]
  iqr    <- stats::IQR(r)
  q      <- stats::quantile(r, c(0.25, 0.75))
  in_iqr <- r >= (q[1] - 2.5 * iqr) & r <= (q[2] + 2.5 * iqr)

  fit_valid <- base_valid
  fit_valid[base_valid] <- in_iqr

  if (sum(fit_valid) < 10) {
    warning("Too few bins after outlier filtering, skipping GC correction.")
    return(reads)
  }

  # Step 3: LOESS with span = 0.03 as in HMMcopy, fitted on clean data
  gc_fit    <- gc[fit_valid]
  reads_fit <- reads[fit_valid]
  fit       <- stats::loess(reads_fit ~ gc_fit, span = 0.03)
  predicted <- stats::predict(fit, newdata = data.frame(gc_fit = gc), se = FALSE)

  median_pred <- stats::median(predicted[fit_valid], na.rm = TRUE)

  corrected   <- reads / predicted * median_pred

  # Bins without a valid prediction retain their original value
  bad <- is.na(corrected) | !is.finite(corrected) | corrected < 0
  corrected[bad] <- reads[bad]
  return(corrected)
}

# WES-specific: GC normalization in log2 space (CNVkit flat-reference style).
# Fits LOESS on log2(reads) ~ GC using ONLY the bins indicated by use_for_fit
# (target bins), then predicts for ALL bins.  This avoids the two-population
# problem when target (~100x) and antitarget (~0.1x) bins have very different
# depths: the LOESS learns the GC-depth curve from the reliable target bins and
# applies it to antitarget bins as well.
# Returns NULL if there are not enough fit bins (caller falls back to median).
.wes_loess_log2 <- function(reads, gc, span = 0.3, use_for_fit = NULL) {
  log2_reads <- log2(reads)
  log2_reads[!is.finite(log2_reads)] <- NA_real_

  # Bins eligible for fitting: valid GC, positive reads, within GC range
  eligible <- !is.na(reads) & reads > 0 & !is.na(gc) & gc > 0.2 & gc < 0.8

  # Restrict fitting to the subset indicated by use_for_fit (e.g. target bins)
  fit_eligible <- if (!is.null(use_for_fit)) eligible & use_for_fit else eligible

  if (sum(fit_eligible) < 30L) {
    warning("Too few valid WES bins for GC-LOESS (need >= 30); ",
            "falling back to median normalization.")
    return(NULL)
  }

  # IQR outlier filter on the fit-eligible subset
  r        <- log2_reads[fit_eligible]
  iqr      <- stats::IQR(r)
  q        <- stats::quantile(r, c(0.25, 0.75))
  fit_mask <- fit_eligible
  fit_mask[fit_eligible] <- r >= (q[1] - 2.5 * iqr) & r <= (q[2] + 2.5 * iqr)

  if (sum(fit_mask) < 30L) {
    warning("Too few WES bins after IQR filter; falling back to median normalization.")
    return(NULL)
  }

  # Fit LOESS on target bins, predict for ALL bins
  fit <- stats::loess(log2_reads[fit_mask] ~ gc[fit_mask], span = span)
  predicted <- tryCatch(
    stats::predict(fit, newdata = data.frame(`gc[fit_mask]` = gc,
                                              check.names  = FALSE), se = FALSE),
    error = function(e) rep(NA_real_, length(gc))
  )

  # Residual = log2(observed / GC-expected) — diploid bins cluster around 0
  log2_ratio <- log2_reads - predicted
  log2_ratio[is.na(reads) | reads == 0L] <- NA_real_
  log2_ratio
}

# Splits large target bins into sub-bins of approximately bin_size bp,
# mirroring CNVkit's `target --split` command.
# Strategy: try valr::bed_makewindows() first (BED-standard, no edge-case
# arithmetic); fall back to the original loop implementation if valr is
# unavailable or raises an error.
.split_target_bins <- function(targets_gr, bin_size = 267L) {
  bin_size <- as.integer(bin_size)

  # ── valr path ──────────────────────────────────────────────────────────────
  if (requireNamespace("valr", quietly = TRUE)) {
    result <- tryCatch({
      bed_df <- data.frame(
        chrom = as.character(GenomicRanges::seqnames(targets_gr)),
        start = GenomicRanges::start(targets_gr) - 1L,   # GRanges→0-based
        end   = GenomicRanges::end(targets_gr),
        stringsAsFactors = FALSE
      )
      split_df <- valr::bed_makewindows(bed_df, win_size = bin_size)
      .valr_to_granges(split_df)
    }, error = function(e) {
      warning("valr::bed_makewindows failed in .split_target_bins (",
              conditionMessage(e), "); using loop fallback.")
      NULL
    })
    if (!is.null(result)) return(result)
  }

  # ── original loop fallback ─────────────────────────────────────────────────
  starts   <- GenomicRanges::start(targets_gr)
  ends     <- GenomicRanges::end(targets_gr)
  chrs     <- as.character(GenomicRanges::seqnames(targets_gr))
  widths   <- ends - starts + 1L
  n_bins   <- pmax(1L, as.integer(round(widths / bin_size)))

  total   <- sum(n_bins)
  out_chr <- character(total)
  out_s   <- integer(total)
  out_e   <- integer(total)

  pos <- 1L
  for (i in seq_along(targets_gr)) {
    n <- n_bins[i]
    if (n == 1L) {
      out_chr[pos] <- chrs[i]
      out_s[pos]   <- starts[i]
      out_e[pos]   <- ends[i]
      pos <- pos + 1L
    } else {
      bp <- as.integer(round(seq(starts[i], ends[i] + 1L, length.out = n + 1L)))
      for (j in seq_len(n)) {
        out_chr[pos] <- chrs[i]
        out_s[pos]   <- bp[j]
        out_e[pos]   <- max(bp[j + 1L] - 1L, bp[j])   # guard against out_e < out_s
        pos <- pos + 1L
      }
    }
  }

  GenomicRanges::GRanges(
    seqnames = out_chr,
    ranges   = IRanges::IRanges(start = out_s, end = out_e)
  )
}

# Computes antitarget bins for WES: genomic gaps between captured exons,
# split into fixed-size windows.  Mirrors the CNVkit antitarget strategy.
# min_size  — minimum gap width to keep as antitarget (default 10 kb)
# bin_size  — split large gaps into bins of this size (default 100 kb)
#
# Strategy: try valr (bed_complement + bed_makewindows) first; fall back to
# the original nested-loop implementation if valr is unavailable or errors.
.compute_antitarget_bins <- function(targets_gr, chr_lens,
                                      min_size = 10000L,
                                      bin_size  = 100000L) {
  min_size <- as.integer(min_size)
  bin_size <- as.integer(bin_size)

  # ── valr path ──────────────────────────────────────────────────────────────
  if (requireNamespace("valr", quietly = TRUE)) {
    result <- tryCatch({
      genome_df <- data.frame(
        chrom = names(chr_lens),
        size  = as.integer(chr_lens),
        stringsAsFactors = FALSE
      )
      bed_df <- data.frame(
        chrom = as.character(GenomicRanges::seqnames(targets_gr)),
        start = GenomicRanges::start(targets_gr) - 1L,   # GRanges→0-based
        end   = GenomicRanges::end(targets_gr),
        stringsAsFactors = FALSE
      )
      merged  <- valr::bed_merge(bed_df)
      gaps    <- valr::bed_complement(merged, genome_df)
      gaps    <- gaps[gaps$end - gaps$start >= min_size, , drop = FALSE]
      if (nrow(gaps) == 0L) return(NULL)
      anti_df <- valr::bed_makewindows(gaps, win_size = bin_size)
      anti_df <- anti_df[anti_df$end - anti_df$start >= min_size, , drop = FALSE]
      if (nrow(anti_df) == 0L) return(NULL)
      .valr_to_granges(anti_df)
    }, error = function(e) {
      warning("valr failed in .compute_antitarget_bins (",
              conditionMessage(e), "); using loop fallback.")
      NULL
    })
    if (!is.null(result) || inherits(result, "GRanges")) return(result)
    # result == NULL returned from tryCatch means 0 antitarget bins (not an error)
    if (is.null(result)) return(NULL)
  }

  # ── original nested-loop fallback ─────────────────────────────────────────
  targets_merged  <- GenomicRanges::reduce(targets_gr)
  antitarget_list <- list()

  for (ch in names(chr_lens)) {
    chr_len  <- as.integer(chr_lens[[ch]])
    chr_tgts <- targets_merged[
      as.character(GenomicRanges::seqnames(targets_merged)) == ch
    ]

    t_starts <- sort(GenomicRanges::start(chr_tgts))
    t_ends   <- sort(GenomicRanges::end(chr_tgts))

    if (length(t_starts) == 0L) {
      gap_s <- 1L
      gap_e <- chr_len
    } else {
      gap_s <- c(1L,            t_ends   + 1L)
      gap_e <- c(t_starts - 1L, chr_len)
    }

    gap_len <- gap_e - gap_s + 1L
    keep    <- gap_len >= min_size & gap_s <= gap_e & gap_s >= 1L
    gap_s   <- gap_s[keep]
    gap_e   <- gap_e[keep]

    if (length(gap_s) == 0L) next

    bins_s <- bins_e <- integer(0)
    for (j in seq_along(gap_s)) {
      s      <- gap_s[j]
      e      <- gap_e[j]
      starts <- seq(s, e, by = bin_size)
      ends   <- pmin(starts + bin_size - 1L, e)
      ok     <- (ends - starts + 1L) >= min_size
      bins_s <- c(bins_s, starts[ok])
      bins_e <- c(bins_e, ends[ok])
    }

    if (length(bins_s) == 0L) next

    antitarget_list[[ch]] <- GenomicRanges::GRanges(
      seqnames = ch,
      ranges   = IRanges::IRanges(start = bins_s, end = bins_e)
    )
  }

  if (length(antitarget_list) == 0L) return(NULL)
  do.call(c, unname(antitarget_list))
}

# Computes GC content per bin directly from BSgenome (WES mode).
# BSgenome packages always use UCSC chr-prefixed names; the function
# temporarily normalises seqnames for the lookup regardless of BAM naming.
.compute_wes_gc <- function(targets_gr, genome_build) {
  bsgenome_pkg <- switch(genome_build,
    hg19 = "BSgenome.Hsapiens.UCSC.hg19",
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    stop("genome_build must be 'hg19' or 'hg38' for WES GC correction.")
  )
  for (pkg in c(bsgenome_pkg, "Biostrings")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("WES GC correction requires ", pkg, ".\n",
           "Install with: BiocManager::install('", pkg, "')")
  }

  bsg <- getExportedValue(bsgenome_pkg, bsgenome_pkg)

  # BSgenome uses chr prefix; rename seqlevels if BAM-normalised GRanges lacks it
  lookup_gr  <- targets_gr
  cur_levels <- GenomeInfoDb::seqlevels(lookup_gr)
  if (!all(grepl("^chr", cur_levels))) {
    new_levels <- ifelse(grepl("^chr", cur_levels),
                         cur_levels, paste0("chr", cur_levels))
    lookup_gr <- GenomeInfoDb::renameSeqlevels(
      lookup_gr, stats::setNames(new_levels, cur_levels)
    )
  }

  # Restrict to chromosomes present in the BSgenome object
  bsg_chrs <- as.character(GenomeInfoDb::seqnames(bsg))
  valid    <- as.character(GenomicRanges::seqnames(lookup_gr)) %in% bsg_chrs
  gc_out   <- rep(NA_real_, length(lookup_gr))

  if (sum(valid) > 0) {
    seqs <- BSgenome::getSeq(bsg, lookup_gr[valid])
    gc_out[valid] <- as.numeric(
      Biostrings::letterFrequency(seqs, "GC", as.prob = TRUE)[, 1]
    )
  }
  gc_out
}

# Normalizes chromosome names to UCSC style "chrN" (with chr prefix)
.norm_chr_to_ucsc <- function(x) {
  x <- as.character(x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

# Normalizes chromosome names to integer (no prefix; X→23, Y→24)
.norm_chr_to_int <- function(x) {
  x <- gsub("^chr", "", as.character(x))
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  suppressWarnings(as.integer(x))
}

# Converts a valr-style BED data.frame (chrom, start, end; 0-based half-open)
# to a GRanges object (1-based closed) required by Rsamtools::countBam().
# Extra columns (e.g. .win_id from bed_makewindows) are silently ignored.
.valr_to_granges <- function(bed_df) {
  GenomicRanges::GRanges(
    seqnames = bed_df[["chrom"]],
    ranges   = IRanges::IRanges(start = bed_df[["start"]] + 1L,
                                end   = bed_df[["end"]])
  )
}


# ─── read_bam_qc ─────────────────────────────────────────────────────────────

#' Per-bin read counting from a BAM file
#'
#' Divides the genome into fixed-size bins, counts reads per bin using
#' Rsamtools, and computes the median-normalised log2 ratio.
#' This is the required input for run_cbs_segmentation().
#'
#' @param bam_path Character path to the BAM file (must be indexed, .bai).
#' @param bin_size Integer bin size in bp (default 500000 for WES,
#'   use 1000000 for low-pass Nanopore).
#' @param genome_build Character "hg19" or "hg38" (default "hg19"). Used
#'   only as metadata; chromosome sizes are read from the BAM header.
#' @param min_mapq Integer minimum mapping quality to include (default 20).
#' @param gc_correction Logical, apply HMMcopy-style GC content correction
#'   using pre-computed tables included in the package (default FALSE).
#'   Recommended for low-pass WGS and WES.
#' @param sequencing_type Character "illumina" or "nanopore". Controls
#'   mappability filtering (Illumina only) and smooth.region in CBS.
#'   Does NOT determine WGS vs WES — use experiment_type for that.
#' @param experiment_type Character "wes" or "wgs" (default "wgs").
#'   "wes": bins are created from targets_bed (or the standard exome BED);
#'   bins with zero reads become NA instead of a large negative log2.
#'   "wgs": fixed genome-wide bins via tileGenome; zero-read bins keep
#'   their log2 value (needed to detect deep deletions).
#' @param min_mappability Numeric minimum mappability threshold for filtering
#'   bins in Illumina data (default 0.9). Ignored for Nanopore.
#' @param targets_bed Character path to a BED file with capture kit regions.
#'   Only used when experiment_type = "wes". When NULL, the standard exome BED
#'   bundled in the package is used automatically (requires running
#'   generate_standard_exome_bed.R once). Default NULL.
#' @param min_reads Integer minimum number of reads for a bin to be considered
#'   valid. Applied only to target bins in WES mode; antitarget bins are
#'   filtered by the LOESS IQR step instead (default 0, no filtering).
#' @param use_antitarget Logical, only used when experiment_type = "wes".
#'   If TRUE, computes off-target bins (complement of the capture kit regions)
#'   and includes them in read counting and GC-LOESS normalisation, mirroring
#'   CNVkit's antitarget strategy.  Improves detection in gene-poor regions.
#'   Requires gc_correction = TRUE (BSgenome).  Default FALSE.
#' @param antitarget_bin_size Integer bin size for antitarget regions in bp
#'   (default 100000 = 100 kb).  Large gaps are split into bins of this size.
#' @param split_targets Logical or NULL. Only used when experiment_type = "wes".
#'   If TRUE, splits each capture region into sub-bins of approximately
#'   target_bin_size bp, mirroring CNVkit's \code{target --split} command.
#'   This increases resolution within large exons and gives the GC-LOESS fit
#'   more data points, improving GC correction accuracy.
#'   NULL (default) means automatic: TRUE for WES, FALSE for WGS.
#' @param target_bin_size Integer target sub-bin size in bp when
#'   split_targets = TRUE (default 267, the CNVkit default).  Exons smaller
#'   than this value are kept as single bins.
#' @param pon Optional Panel of Normals object produced by build_pon(). When
#'   provided, the per-bin PoN median log2 is subtracted from the tumour log2
#'   ratios after GC correction, removing systematic protocol-specific
#'   artefacts that are consistent across samples.  Must have been built with
#'   the same genome_build, experiment_type, targets_bed, and gc_correction
#'   settings.  Default NULL (no PoN normalization).
#'
#' @return Data frame with columns: chr, loc.start, loc.end, reads, log2_ratio,
#'   bin_type.  If gc_correction = TRUE also includes gc.
#'   If pon is provided, log2_ratio reflects PoN-normalised values.
#'
#' @details
#' Requires: Rsamtools, GenomicRanges (BiocManager::install(...))
#'
#' **WES vs WGS mode**
#' Use experiment_type = "wes" for exome/targeted data and "wgs" for
#' whole-genome. sequencing_type only affects mappability filtering and
#' CBS smoothing — you can have Illumina WGS or Nanopore WES independently.
#'
#' **Automatic standard exome BED**
#' When experiment_type = "wes" and targets_bed = NULL, the function
#' automatically uses the standard exome BED bundled in the package
#' (inst/extdata/standard_exome_hg19.bed or standard_exome_hg38.bed).
#' Generate it once with: source("inst/scripts/generate_standard_exome_bed.R")
#'
#' **GC correction — WGS vs WES**
#' For WGS, GC correction follows the HMMcopy approach: IQR-based outlier
#' filtering + LOESS regression (span = 0.03) on raw read counts, then
#' median normalisation. Pre-computed GC/mappability tables are used
#' (500kb, 1Mb, 2Mb, 10Mb; hg19/hg38).
#'
#' For WES with gc_correction = TRUE, a different (CNVkit flat-reference style)
#' approach is used: LOESS is fitted on log2(reads) ~ GC in log2 space
#' (span = 0.3). The residuals are the GC-corrected log2 ratios with a diploid
#' baseline of 0, making them directly comparable to CNVkit output and
#' robust to tumour-wide ploidy shifts that would bias median normalisation.
#'
#' @examples
#' \dontrun{
#'   # Nanopore low-pass WGS
#'   qc <- read_bam_qc("sample.bam", bin_size = 1e6,
#'                      sequencing_type = "nanopore", experiment_type = "wgs")
#'
#'   # Illumina WGS (NOT exome)
#'   qc <- read_bam_qc("illumina_wgs.bam", bin_size = 5e5,
#'                      sequencing_type = "illumina", experiment_type = "wgs")
#'
#'   # Illumina WES — standard exome BED used automatically
#'   qc <- read_bam_qc("sample_wes.bam",
#'                      sequencing_type = "illumina", experiment_type = "wes",
#'                      min_reads = 10L)
#'
#'   # Illumina WES — override with your own capture kit BED
#'   qc <- read_bam_qc("sample_wes.bam", targets_bed = "capture_kit.bed",
#'                      sequencing_type = "illumina", experiment_type = "wes",
#'                      min_reads = 10L)
#' }
#'
#' @export
read_bam_qc <- function(bam_path, bin_size = 500000L, genome_build = "hg19",
                         min_mapq = 20L, gc_correction = FALSE,
                         sequencing_type = c("illumina", "nanopore"),
                         experiment_type  = c("wgs", "wes"),
                         min_mappability = 0.9,
                         targets_bed = NULL,
                         min_reads = 0L,
                         use_antitarget = FALSE,
                         antitarget_bin_size = 100000L,
                         split_targets = NULL,
                         target_bin_size = 267L,
                         pon = NULL) {
  sequencing_type <- match.arg(sequencing_type)
  experiment_type  <- match.arg(experiment_type)

  # split_targets = NULL → auto: TRUE for WES (best default), FALSE for WGS
  if (is.null(split_targets))
    split_targets <- (experiment_type == "wes")

  # Auto-use standard exome BED in WES mode when no BED is provided.
  # Works regardless of sequencing_type (handles Illumina WES and any other).
  if (experiment_type == "wes" && is.null(targets_bed)) {
    std_bed <- system.file("extdata",
                            sprintf("standard_exome_%s.bed", genome_build),
                            package = "CNAppR")
    if (nzchar(std_bed) && file.exists(std_bed)) {
      targets_bed <- std_bed
      message("  Using standard exome BED (", genome_build,
              "). Override with the targets_bed parameter.")
    } else {
      warning("Standard exome BED not found for ", genome_build, ".\n",
              "Run inst/scripts/generate_standard_exome_bed.R once ",
              "to generate it.")
    }
  }

  .check_bam(bam_path)
  .check_pkg("Rsamtools")
  .check_pkg("GenomicRanges")

  # Chromosome sizes from the BAM header (no BSgenome dependency)
  header   <- Rsamtools::scanBamHeader(bam_path)[[1]]
  chr_lens <- header$targets   # named integer: chr_name → length

  # Keep only standard chromosomes (1-22, X, Y), with or without "chr" prefix
  keep     <- grepl("^(chr)?([0-9]{1,2}|X|Y)$", names(chr_lens))
  chr_lens <- chr_lens[keep]
  if (length(chr_lens) == 0)
    stop("No standard chromosomes found in the BAM header.")

  # Create bins: from targets BED (WES) or fixed genome-wide windows (WGS)
  if (!is.null(targets_bed)) {
    message("  WES mode: loading targets from BED...")
    original_targets <- .load_targets_bed(targets_bed, chr_lens)

    if (split_targets && experiment_type == "wes") {
      message("  Splitting target bins (~", target_bin_size, " bp each, ",
              "CNVkit style)...")
      bins <- .split_target_bins(original_targets,
                                  bin_size = as.integer(target_bin_size))
      message("  ", length(original_targets), " exons -> ",
              length(bins), " bins after split")
    } else {
      bins <- original_targets
    }
  } else {
    original_targets <- NULL
    # WGS tiling: try valr::bed_makewindows() first; fall back to tileGenome.
    if (requireNamespace("valr", quietly = TRUE)) {
      bins <- tryCatch({
        genome_df <- data.frame(
          chrom = names(chr_lens),
          size  = as.integer(chr_lens),
          stringsAsFactors = FALSE
        )
        .valr_to_granges(
          valr::bed_makewindows(genome_df, win_size = as.integer(bin_size))
        )
      }, error = function(e) {
        warning("valr::bed_makewindows failed for WGS tiling (",
                conditionMessage(e), "); using tileGenome fallback.")
        GenomicRanges::tileGenome(
          seqlengths             = chr_lens,
          tilewidth              = as.integer(bin_size),
          cut.last.tile.in.chrom = TRUE
        )
      })
    } else {
      bins <- GenomicRanges::tileGenome(
        seqlengths             = chr_lens,
        tilewidth              = as.integer(bin_size),
        cut.last.tile.in.chrom = TRUE
      )
    }
  }

  # WES antitarget bins: genomic gaps between captured exons.
  # Computed from the ORIGINAL (pre-split) exon regions so that the
  # complement is correct regardless of whether split_targets is TRUE.
  n_target <- length(bins)
  bin_type <- rep("target", n_target)

  if (experiment_type == "wes" && use_antitarget) {
    if (!gc_correction)
      stop("use_antitarget = TRUE requires gc_correction = TRUE ",
           "(BSgenome is needed to compute GC content for antitarget bins).")
    anti_source <- if (!is.null(original_targets)) original_targets else bins
    message("  Computing antitarget bins (bin size: ",
            round(antitarget_bin_size / 1000), " kb, min gap: 10 kb)...")
    anti_gr <- .compute_antitarget_bins(anti_source, chr_lens,
                                         min_size = 10000L,
                                         bin_size  = as.integer(antitarget_bin_size))
    if (!is.null(anti_gr)) {
      message("  Antitarget bins: ", length(anti_gr),
              "  (target bins: ", n_target, ")")
      bins     <- c(bins, anti_gr)
      bin_type <- c(rep("target",     n_target),
                    rep("antitarget", length(anti_gr)))
    } else {
      message("  No antitarget regions found.")
    }
  }

  # Count reads per bin (only reads with MAPQ >= min_mapq)
  param  <- Rsamtools::ScanBamParam(
    which       = bins,
    mapqFilter  = as.integer(min_mapq)
  )
  counts <- Rsamtools::countBam(bam_path, param = param)

  reads <- counts$records
  total <- sum(reads)
  if (total == 0)
    stop("No reads counted in '", basename(bam_path),
         "'. Check that the BAM is not empty and the index is valid.")

  gc          <- NULL
  mappability <- NULL
  log2_ratio  <- NULL   # set directly for WES + gc_correction (LOESS-log2 path)

  if (gc_correction) {
    if (experiment_type == "wes") {
      # WES: compute GC per exon from BSgenome, then normalise in log2 space
      # (CNVkit flat-reference style: residuals of log2(reads) ~ GC LOESS).
      # min_reads filter is applied before the fit so low-coverage bins don't
      # bias the regression.
      message("  Computing GC content per exon from BSgenome (", genome_build, ")...")
      gc          <- .compute_wes_gc(bins, genome_build)
      mappability <- NULL

      reads_filtered <- reads
      if (min_reads > 0L) {
        # Apply min_reads only to target bins; antitarget bins have inherently
        # lower coverage and are cleaned by the LOESS IQR filter instead.
        target_idx <- bin_type == "target"
        reads_filtered[target_idx & !is.na(reads_filtered) &
                         reads_filtered < min_reads] <- NA_real_
      }

      message("  Applying GC-LOESS in log2 space (WES flat-reference normalization)...")
      log2_ratio <- .wes_loess_log2(reads_filtered, gc,
                                     use_for_fit = (bin_type == "target"))

      if (!is.null(log2_ratio)) {
        reads <- reads_filtered   # store filtered counts in output
      } else {
        message("  GC-LOESS fallback: using median normalization.")
        reads <- reads_filtered
      }
    } else {
      message("  Loading bin annotations (GC + mappability)...")
      ann         <- .load_bin_annotations(bins, genome_build, bin_size)
      gc          <- ann$gc
      mappability <- ann$mappability

      # Filter bins with low mappability (Illumina only; Nanopore maps long reads)
      if (sequencing_type == "illumina" && !all(is.na(mappability))) {
        low_map        <- !is.na(mappability) & mappability < min_mappability
        reads[low_map] <- NA_real_
        message("  Bins filtered by mappability < ", min_mappability, ": ", sum(low_map))
      }

      message("  Applying LOESS GC correction...")
      reads_for_loess        <- reads
      reads_for_loess[is.na(reads_for_loess)] <- 0
      reads_corrected        <- .correct_gc_loess(reads_for_loess, gc)
      reads_corrected[is.na(reads)] <- NA_real_
      reads <- reads_corrected
    }
  }

  # ── Median normalization (WGS, or WES without gc_correction, or LOESS fallback)
  if (is.null(log2_ratio)) {
    # Minimum reads filter (target bins only in WES mode)
    if (min_reads > 0L) {
      target_idx <- bin_type == "target"
      reads[target_idx & !is.na(reads) & reads < min_reads] <- NA_real_
    }

    valid_reads <- reads[!is.na(reads) & reads > 0]
    if (length(valid_reads) == 0)
      stop("No valid reads after filtering.")
    median_reads <- stats::median(valid_reads)

    if (experiment_type == "wes") {
      # WES: bins with 0 reads are uncaptured regions → NA instead of large negative log2
      log2_ratio <- ifelse(is.na(reads) | reads == 0L, NA_real_,
                           log2(reads / median_reads))
    } else {
      # WGS: pmax avoids log2(0); zero-read bins keep a large negative log2
      # so deep deletions are visible to CBS
      log2_ratio <- log2(pmax(reads, 0.5) / median_reads)
      log2_ratio[is.na(reads)] <- NA_real_
    }
  }

  out <- data.frame(
    chr        = as.character(GenomicRanges::seqnames(bins)),
    loc.start  = GenomicRanges::start(bins),
    loc.end    = GenomicRanges::end(bins),
    reads      = reads,
    log2_ratio = log2_ratio,
    bin_type   = bin_type,
    stringsAsFactors = FALSE
  )
  if (!is.null(gc))          out$gc          <- gc
  if (!is.null(mappability)) out$mappability <- mappability

  # ── PoN normalization ──────────────────────────────────────────────────────
  if (!is.null(pon)) {
    if (!is.list(pon) || is.null(pon$median_log2) || is.null(pon$bins))
      stop("pon must be a PoN object produced by build_pon().")

    if (!isTRUE(pon$gc_correction == gc_correction))
      warning("PoN was built with gc_correction = ", pon$gc_correction,
              " but this call uses gc_correction = ", gc_correction,
              ". Log2 values may not be on the same scale.")

    tumor_key <- paste(out$chr, out$loc.start)
    pon_key   <- paste(pon$bins$chr, pon$bins$loc.start)
    pon_idx   <- match(tumor_key, pon_key)

    n_matched <- sum(!is.na(pon_idx))
    if (n_matched == 0L)
      stop("No bins matched between tumour and PoN. ",
           "Ensure the same genome_build, experiment_type, and targets_bed.")
    if (n_matched < nrow(out))
      warning(nrow(out) - n_matched,
              " tumour bins have no PoN entry; their log2_ratio set to NA.")

    pon_log2          <- rep(NA_real_, nrow(out))
    matched           <- !is.na(pon_idx)
    pon_log2[matched] <- pon$median_log2[pon_idx[matched]]

    out$log2_ratio <- out$log2_ratio - pon_log2
    message("  PoN normalization: ", n_matched, "/", nrow(out),
            " bins matched (", pon$n_samples, " normals in PoN).")
  }

  out
}


# ─── run_cbs_segmentation ────────────────────────────────────────────────────

#' CBS segmentation (Circular Binary Segmentation)
#'
#' Runs CBS on the log2 ratios produced by read_bam_qc() and returns a
#' table in the input format expected by resegment_sample().
#'
#' @param qc_data Data frame output from read_bam_qc().
#' @param sample_id Character sample identifier.
#' @param alpha Numeric CBS significance level (default 0.01).
#' @param nperm Integer number of CBS permutations (default 10000).
#' @param sequencing_type Character "illumina" or "nanopore". Controls
#'   outlier smoothing: more aggressive for Nanopore long reads.
#'
#' @return Data frame with columns: ID, chr, loc.start, loc.end, seg.mean.
#'   Ready to pass directly to resegment_sample().
#'
#' @details
#' Requires: DNAcopy (BiocManager::install("DNAcopy"))
#'
#' CBS is the most widely used segmentation algorithm for detecting
#' log ratios in copy number data. Implemented in the R package DNAcopy.
#'
#' @examples
#' \dontrun{
#'   qc  <- read_bam_qc("sample.bam", bin_size = 1e6)
#'   seg <- run_cbs_segmentation(qc, sample_id = "sample_01")
#'   head(seg)
#' }
#'
#' @export
run_cbs_segmentation <- function(qc_data, sample_id,
                                  alpha = 0.01, nperm = 10000L,
                                  sequencing_type = c("illumina", "nanopore")) {
  sequencing_type <- match.arg(sequencing_type)
  .check_pkg("DNAcopy")

  if (!all(c("chr", "loc.start", "loc.end", "log2_ratio") %in% colnames(qc_data)))
    stop("qc_data must have columns: chr, loc.start, loc.end, log2_ratio.",
         "\nUse read_bam_qc() to generate it.")

  # Normalize chromosomes to integer for DNAcopy
  chr_int <- .norm_chr_to_int(qc_data$chr)
  valid   <- !is.na(chr_int)
  if (!any(valid)) stop("No valid chromosomes in qc_data.")

  qc_data <- qc_data[valid, ]
  chr_int <- chr_int[valid]

  # Bin midpoint as coordinate for CBS
  pos <- as.integer((qc_data$loc.start + qc_data$loc.end) / 2L)

  # CNA object for DNAcopy
  cna_obj <- DNAcopy::CNA(
    genomdat  = qc_data$log2_ratio,
    chrom     = chr_int,
    maploc    = pos,
    data.type = "logratio",
    sampleid  = sample_id
  )

  # Outlier smoothing (more aggressive for Nanopore due to long-read noise)
  smooth_region <- if (sequencing_type == "nanopore") 5L else 2L
  cna_smooth    <- DNAcopy::smooth.CNA(cna_obj, smooth.region = smooth_region)

  # Segmentation
  seg    <- DNAcopy::segment(cna_smooth,
                              alpha   = alpha,
                              nperm   = as.integer(nperm),
                              verbose = 0)
  seg_df <- seg$output

  if (nrow(seg_df) == 0)
    stop("CBS produced no segments for sample '", sample_id, "'.")

  # CNAppR format: ID, chr, loc.start, loc.end, seg.mean
  data.frame(
    ID        = sample_id,
    chr       = as.integer(seg_df$chrom),
    loc.start = as.integer(seg_df$loc.start),
    loc.end   = as.integer(seg_df$loc.end),
    seg.mean  = as.numeric(seg_df$seg.mean),
    stringsAsFactors = FALSE
  )
}


# ─── build_pon ───────────────────────────────────────────────────────────────

#' Build a Panel of Normals (PoN) for copy-number normalisation
#'
#' Processes a set of normal BAM files through the same read-counting and
#' GC-correction pipeline used by read_bam_qc(), then computes a per-bin
#' median log2 that represents the expected diploid baseline for the protocol.
#' The resulting PoN object can be passed to read_bam_qc(pon = ...) or
#' run_bam_pipeline(pon = ...) to remove systematic protocol-specific
#' artefacts — wave patterns, GC bias, consistently poorly-captured exons —
#' that are not removed by single-sample GC correction alone.
#'
#' @param bam_paths Character vector of paths to normal BAM files.  Each must
#'   have a corresponding .bai index.  Use at least 5 normals for stable
#'   medians; 10–20 is typical.
#' @param genome_build Character "hg19" or "hg38" (default "hg19").
#' @param sequencing_type Character "illumina" or "nanopore".
#' @param experiment_type Character "wgs" or "wes".
#' @param bin_size Integer bin size in bp (WGS only; default 500000).
#' @param targets_bed Character path to the capture kit BED (WES only).
#'   When NULL the standard exome BED bundled in the package is used.
#' @param min_mapq Integer minimum mapping quality (default 20).
#' @param gc_correction Logical, apply GC correction when counting each normal
#'   (default TRUE).  Must match the gc_correction argument used in the tumour
#'   call so that log2 values are on the same scale.
#' @param use_antitarget Logical, include antitarget bins (WES only, requires
#'   gc_correction = TRUE; default FALSE).
#' @param antitarget_bin_size Integer bin size for antitarget regions in bp
#'   (default 100000).
#' @param output_file Character path to save the PoN as an RDS file (optional).
#'   Reload later with pon <- readRDS(output_file).
#'
#' @return A named list (invisibly) containing:
#'   \describe{
#'     \item{bins}{Data frame: chr, loc.start, loc.end, bin_type}
#'     \item{median_log2}{Numeric vector, per-bin median log2 across normals}
#'     \item{n_samples}{Number of normal BAMs processed}
#'     \item{genome_build, experiment_type, sequencing_type, gc_correction}{
#'       Parameters used — must match the tumour call}
#'     \item{created}{POSIXct timestamp}
#'   }
#'
#' @examples
#' \dontrun{
#'   normals <- c("normal1.bam", "normal2.bam", "normal3.bam")
#'   pon <- build_pon(normals,
#'                    genome_build    = "hg19",
#'                    experiment_type = "wes",
#'                    gc_correction   = TRUE,
#'                    output_file     = "my_pon.rds")
#'
#'   # Later:
#'   pon <- readRDS("my_pon.rds")
#'   qc  <- read_bam_qc("tumour.bam",
#'                       experiment_type = "wes",
#'                       gc_correction   = TRUE,
#'                       pon             = pon)
#' }
#'
#' @export
build_pon <- function(bam_paths,
                      genome_build        = "hg19",
                      sequencing_type     = c("illumina", "nanopore"),
                      experiment_type     = c("wgs", "wes"),
                      bin_size            = 500000L,
                      targets_bed         = NULL,
                      min_mapq            = 20L,
                      gc_correction       = TRUE,
                      use_antitarget      = FALSE,
                      antitarget_bin_size = 100000L,
                      split_targets       = FALSE,
                      target_bin_size     = 267L,
                      output_file         = NULL) {
  sequencing_type <- match.arg(sequencing_type)
  experiment_type <- match.arg(experiment_type)

  if (length(bam_paths) == 0L)
    stop("bam_paths must contain at least one BAM file.")
  if (length(bam_paths) < 2L)
    warning("A PoN with a single sample has limited statistical power; ",
            "use >= 5 normals for reliable medians.")

  message("Building Panel of Normals from ", length(bam_paths), " sample(s)...")

  log2_list <- vector("list", length(bam_paths))
  bins_ref  <- NULL

  for (i in seq_along(bam_paths)) {
    message("  [", i, "/", length(bam_paths), "] ", basename(bam_paths[i]))
    qc_i <- read_bam_qc(
      bam_path            = bam_paths[i],
      bin_size            = bin_size,
      genome_build        = genome_build,
      min_mapq            = min_mapq,
      gc_correction       = gc_correction,
      sequencing_type     = sequencing_type,
      experiment_type     = experiment_type,
      targets_bed         = targets_bed,
      min_reads           = 0L,
      use_antitarget      = use_antitarget,
      antitarget_bin_size = antitarget_bin_size,
      split_targets       = split_targets,
      target_bin_size     = target_bin_size
    )

    if (is.null(bins_ref)) {
      bins_ref <- qc_i[, c("chr", "loc.start", "loc.end", "bin_type")]
    } else if (nrow(qc_i) != nrow(bins_ref)) {
      stop("BAM [", i, "] produced ", nrow(qc_i), " bins but expected ",
           nrow(bins_ref), ". All normals must use the same genome_build ",
           "and targets_bed.")
    }

    log2_list[[i]] <- qc_i$log2_ratio
  }

  # Per-bin median log2 across all normal samples (NA-aware)
  log2_mat    <- do.call(cbind, log2_list)
  median_log2 <- apply(log2_mat, 1L,
                       function(x) stats::median(x, na.rm = TRUE))

  n_valid <- sum(!is.na(median_log2))
  message("PoN complete: ", nrow(bins_ref), " bins, ",
          n_valid, " with valid median log2.")

  pon <- list(
    bins            = bins_ref,
    median_log2     = median_log2,
    n_samples       = length(bam_paths),
    genome_build    = genome_build,
    experiment_type = experiment_type,
    sequencing_type = sequencing_type,
    gc_correction   = gc_correction,
    bin_size        = bin_size,
    created         = Sys.time()
  )

  if (!is.null(output_file)) {
    saveRDS(pon, output_file)
    message("PoN saved to: ", output_file)
  }

  invisible(pon)
}


# ─── run_bam_pipeline ────────────────────────────────────────────────────────

#' Full BAM → CBS → CNAppR pipeline
#'
#' Convenience wrapper that runs read_bam_qc(), run_cbs_segmentation()
#' and resegment_sample() in a single call.
#'
#' @param bam_path Character path to the indexed BAM file.
#' @param sample_id Character sample identifier.
#' @param sequencing_type Character "illumina" or "nanopore".
#' @param bin_size Integer bin size in bp for read counting.
#'   Default 500000 (Illumina) or 1000000 (Nanopore).
#' @param genome_build Character "hg19" or "hg38" (default "hg19").
#' @param min_mapq Integer minimum mapping quality (default 20).
#' @param alpha Numeric CBS significance level (default 0.01).
#' @param nperm Integer CBS permutations (default 10000).
#' @param experiment_type Character "wgs" or "wes". Controls whether bins come
#'   from fixed genome-wide windows ("wgs") or a capture kit BED ("wes").
#'   Independent of sequencing_type: you can have Illumina WGS or Nanopore WES.
#'   Default "wgs".
#' @param targets_bed Character path to a BED file with capture kit regions.
#'   Only used when experiment_type = "wes". When NULL, the standard exome BED
#'   bundled in the package is used automatically. Default NULL.
#' @param min_reads Integer minimum number of reads for a bin to be considered
#'   valid. Bins below this threshold are set to NA. Recommended min_reads = 10
#'   for WES (default 0, no filtering).
#' @param pon Optional Panel of Normals object produced by build_pon().
#'   Passed through to read_bam_qc(). Default NULL.
#' @param ... Additional arguments passed to resegment_sample()
#'   (e.g. purity, acrocentric, skip_resegmentation).
#'
#' @return Output of resegment_sample(): data frame with classified CNA
#'   segments (chromosomal, arm, focal, diploid).
#'
#' @examples
#' \dontrun{
#'   # Nanopore low-pass WGS
#'   result <- run_bam_pipeline(
#'     "sample_nanopore.bam",
#'     sample_id       = "sample_01",
#'     sequencing_type = "nanopore",
#'     bin_size        = 1000000
#'   )
#'
#'   # Illumina WES — standard exome BED used automatically
#'   result <- run_bam_pipeline(
#'     "sample_wes.bam",
#'     sample_id       = "sample_02",
#'     sequencing_type = "illumina",
#'     min_reads       = 10L
#'   )
#'
#'   # Illumina WES — override with your own capture kit BED
#'   result <- run_bam_pipeline(
#'     "sample_wes.bam",
#'     sample_id       = "sample_02",
#'     sequencing_type = "illumina",
#'     targets_bed     = "capture_kit.bed",
#'     min_reads       = 10L
#'   )
#' }
#'
#' @export
run_bam_pipeline <- function(bam_path,
                              sample_id,
                              sequencing_type     = c("illumina", "nanopore"),
                              experiment_type      = c("wgs", "wes"),
                              bin_size            = NULL,
                              genome_build        = "hg19",
                              min_mapq            = 20L,
                              gc_correction       = FALSE,
                              alpha               = 0.01,
                              nperm               = 10000L,
                              targets_bed         = NULL,
                              min_reads           = 0L,
                              use_antitarget      = FALSE,
                              antitarget_bin_size = 100000L,
                              split_targets       = NULL,
                              target_bin_size     = 267L,
                              pon                 = NULL,
                              ...) {
  sequencing_type <- match.arg(sequencing_type)
  experiment_type  <- match.arg(experiment_type)

  # Default bin size based on sequencing type
  if (is.null(bin_size)) {
    bin_size <- if (sequencing_type == "nanopore") 1000000L else 500000L
  }

  message("Step 1/3: Counting reads per bin from BAM...")
  qc_data <- read_bam_qc(
    bam_path            = bam_path,
    bin_size            = bin_size,
    genome_build        = genome_build,
    min_mapq            = min_mapq,
    gc_correction       = gc_correction,
    sequencing_type     = sequencing_type,
    experiment_type     = experiment_type,
    targets_bed         = targets_bed,
    min_reads           = min_reads,
    use_antitarget      = use_antitarget,
    antitarget_bin_size = antitarget_bin_size,
    split_targets       = split_targets,
    target_bin_size     = target_bin_size,
    pon                 = pon
  )

  message("Step 2/3: CBS segmentation (", sequencing_type, ")...")
  seg_data <- run_cbs_segmentation(
    qc_data         = qc_data,
    sample_id       = sample_id,
    alpha           = alpha,
    nperm           = nperm,
    sequencing_type = sequencing_type
  )

  message("Step 3/3: CNAppR re-segmentation and classification...")
  resegment_sample(seg_data, sample_id = sample_id, ...)
}
