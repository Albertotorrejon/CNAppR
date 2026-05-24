#' Compute CNA Frequency Profiles by Reference Regions
#'
#' Computes per-region copy-number gain and loss frequencies across a cohort
#' using the weighted-mean methodology of the original CNApp tool
#' (Franch-Expósito et al. 2020, eLife 50267).
#'
#' @param resegmented_list List of re-segmented data frames, one per sample
#'   (output from \code{resegment_sample()}). Each data frame must have at
#'   minimum the columns \code{chr} (integer or character), \code{loc.start},
#'   \code{loc.end}, and \code{seg.mean}.
#' @param reference_file Either a data frame with at least three columns
#'   (chr, start, end) or a path to a tab-delimited file with the same
#'   structure (no header). A fourth column \code{name} is preserved in the
#'   output if present.
#' @param gain_thr Numeric gain threshold in log2 copy-number ratio
#'   (default 0.2, corresponding to the CNApp "default" thresholds).
#' @param loss_thr Numeric loss threshold in log2 copy-number ratio
#'   (default -0.2).
#' @param sample_names Character vector of sample identifiers. If \code{NULL}
#'   (default) the names of \code{resegmented_list} are used; if the list is
#'   unnamed, samples are labelled \code{"sample_1"}, \code{"sample_2"}, etc.
#'
#' @return A data frame with the same rows as \code{reference_file} plus two
#'   additional columns:
#'   \describe{
#'     \item{gain_freq}{Proportion of samples (0–1) with a weighted-mean
#'       copy-number ratio \eqn{\ge} \code{gain_thr} in that region.}
#'     \item{loss_freq}{Proportion of samples (0–1) with a weighted-mean
#'       copy-number ratio \eqn{\le} \code{loss_thr} in that region.}
#'   }
#'
#' @details
#' For each sample × reference region, the function computes a length-weighted
#' mean copy-number ratio from all overlapping copy-number segments:
#' \deqn{
#'   \bar{m} = \frac{\sum_j m_j \cdot w_j}{\sum_j w_j}, \quad
#'   w_j = \min\!\left(\frac{\ell_j^{\cap}}{L_{\text{region}}},\, 1\right)
#' }
#' where \eqn{\ell_j^{\cap}} is the overlap length between segment \eqn{j}
#' and the reference region and \eqn{L_{\text{region}}} is the region length.
#'
#' The weighted mean is then classified against \code{gain_thr} and
#' \code{loss_thr}. A region with no overlapping segments receives a mean of
#' zero (neutral).
#'
#' This replicates the pipeline used in the original CNApp Shiny application
#' for cohorts with \eqn{\ge 150} samples (functions
#' \code{scoring.segments}, \code{compare.locations.extract.means}, and
#' \code{classify.values}).
#'
#' @examples
#' \dontrun{
#'   # Arm-level frequencies with default thresholds
#'   arm_ref <- system.file(
#'     "aux_files/segmented_files_hg19/autosomes_hg19_by_arms.txt",
#'     package = "CNAppR"
#'   )
#'   freq <- compute_region_frequencies(reseg_list, arm_ref)
#'   head(freq)
#' }
#'
#' @export
compute_region_frequencies <- function(resegmented_list,
                                       reference_file,
                                       gain_thr    = 0.2,
                                       loss_thr    = -0.2,
                                       sample_names = NULL) {

  # ── validate inputs ──────────────────────────────────────────────────────────
  if (!is.list(resegmented_list))
    stop("'resegmented_list' must be a list of data frames.")

  if (is.null(sample_names)) {
    sample_names <- names(resegmented_list)
    if (is.null(sample_names))
      sample_names <- paste0("sample_", seq_along(resegmented_list))
  }
  n <- length(resegmented_list)

  # ── load reference file ──────────────────────────────────────────────────────
  if (is.character(reference_file) && length(reference_file) == 1L) {
    reference_file <- utils::read.table(reference_file, header = FALSE,
                                        sep = "\t", stringsAsFactors = FALSE)
  }
  if (!is.data.frame(reference_file))
    stop("'reference_file' must be a data frame or a path to a tab-delimited file.")
  if (ncol(reference_file) < 3L)
    stop("'reference_file' must have at least 3 columns: chr, start, end.")

  colnames(reference_file)[1:3] <- c("chr", "start", "end")
  reference_file$start <- as.numeric(reference_file$start)
  reference_file$end   <- as.numeric(reference_file$end)

  # Ensure "chr" prefix in the reference chr column
  reference_file$chr <- ifelse(
    grepl("^chr", reference_file$chr),
    reference_file$chr,
    paste0("chr", reference_file$chr)
  )

  nr          <- nrow(reference_file)
  gain_counts <- integer(nr)
  loss_counts <- integer(nr)

  # ── per-sample loop ──────────────────────────────────────────────────────────
  for (si in seq_len(n)) {
    seg <- resegmented_list[[si]]
    if (is.null(seg) || nrow(seg) == 0L) next

    # Add chr column with "chr" prefix (input is integer from resegment_sample)
    seg$chr_c <- paste0("chr", as.character(seg$chr))

    # Remove degenerate segments
    seg <- seg[!is.na(seg$seg.mean) & seg$loc.end > seg$loc.start, , drop = FALSE]
    if (nrow(seg) == 0L) next

    # ── per-region weighted mean ─────────────────────────────────────────────
    for (ri in seq_len(nr)) {
      rch    <- reference_file$chr[ri]
      rstart <- reference_file$start[ri]
      rend   <- reference_file$end[ri]
      rlen   <- rend - rstart
      if (rlen <= 0) next

      s <- seg[seg$chr_c == rch, , drop = FALSE]
      if (nrow(s) == 0L) next

      ov_idx <- which(s$loc.start <= rend & s$loc.end >= rstart)
      if (length(ov_idx) == 0L) next

      ss <- s$loc.start[ov_idx]
      se <- s$loc.end[ov_idx]
      sv <- s$seg.mean[ov_idx]

      # Overlap length within the reference region
      ov_len <- pmax(pmin(se, rend) - pmax(ss, rstart), 0)

      # Weight = overlap_length / region_length, capped at 1
      wt     <- pmin(ov_len / rlen, 1)
      wt_sum <- sum(wt)

      if (wt_sum == 0) next

      wmean <- sum(sv * wt) / wt_sum

      if (!is.na(wmean)) {
        if (wmean >= gain_thr)
          gain_counts[ri] <- gain_counts[ri] + 1L
        else if (wmean <= loss_thr)
          loss_counts[ri] <- loss_counts[ri] + 1L
      }
    }
  }

  # ── assemble output ──────────────────────────────────────────────────────────
  out           <- reference_file
  out$gain_freq <- gain_counts / n
  out$loss_freq <- loss_counts / n

  out
}
