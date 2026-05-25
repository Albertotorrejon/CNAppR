#' Re-segment Copy Number Data
#'
#' @param data Data frame with columns: ID, chr, loc.start, loc.end, seg.mean, [BAF, purity].
#' @param sample_id Character sample identifier to process.
#' @param min_length Numeric minimum segment length in bp (default 100000).
#' @param max_dist_segm Numeric maximum distance between segments to consider merging (bp, default 1000000).
#' @param percent_dist Numeric percentage-based distance threshold (default 2).
#' @param dev_btw_segs Numeric max deviation between segment means for merging (default 0.16).
#' @param dev_tozero Numeric deviation from zero threshold (default 0.16).
#' @param dev_baf Numeric max BAF deviation between segments (default 0.1).
#' @param low_gain Numeric low gain threshold in log2 (default 0.2).
#' @param normal_gain Numeric normal gain threshold (default 0.58, i.e. round(log2(3/2), 2)).
#' @param high_gain Numeric high gain threshold (default 1.0, i.e. round(log2(4/2), 2)).
#' @param low_loss Numeric low loss threshold (default -0.2).
#' @param normal_loss Numeric normal loss threshold (default -1.0, i.e. round(log2(1/2), 2)).
#' @param high_loss Numeric high loss threshold (default -1.74, i.e. round(log2(0.6/2), 2)).
#' @param chrom_percent Numeric minimum coverage for chromosomal alterations (default 0.9).
#' @param arm_percent Numeric minimum coverage for arm-level alterations (default 0.5).
#' @param focal_percent_low Numeric focal alteration low threshold (default 0.05).
#' @param focal_percent_medium Numeric focal alteration medium threshold (default 0.15).
#' @param focal_percent_high Numeric focal alteration high threshold (default 0.3).
#' @param min_baf Numeric minimum BAF for CN-LOH detection (default 0.2).
#' @param acrocentric Logical whether to ignore p arms of chromosomes 13,14,15,21,22 (default FALSE).
#' @param skip_resegmentation Logical whether to skip re-segmentation (default FALSE).
#'
#' @return Data frame with additional columns: score, classified, type, intensity, weight, comments
#'
#' @export
resegment_sample <- function(data,
                             sample_id,
                             min_length = 100000,
                             max_dist_segm = 1000000,
                             percent_dist = 2,
                             dev_btw_segs = 0.16,
                             dev_tozero = 0.16,
                             dev_baf = 0.1,
                             low_gain = 0.2,
                             normal_gain = round(log2(3 / 2), 2),
                             high_gain = round(log2(4 / 2), 2),
                             low_loss = -0.2,
                             normal_loss = round(log2(1 / 2), 2),
                             high_loss = round(log2(0.6 / 2), 2),
                             chrom_percent = 0.9,
                             arm_percent = 0.5,
                             focal_percent_low = 0.05,
                             focal_percent_medium = 0.15,
                             focal_percent_high = 0.3,
                             min_baf = 0.2,
                             acrocentric = FALSE,
                             skip_resegmentation = FALSE) {

  tryCatch({

  file <- data[data$ID == sample_id, ]
  if (nrow(file) == 0) stop("Sample '", sample_id, "' not found in data.")

  file$length <- file$loc.end - file$loc.start

  # Normalize chr column: remove "chr" prefix, map X->23, Y->24, drop unknowns
  file$chr <- gsub("^chr", "", as.character(file$chr))
  file$chr[file$chr == "X"] <- "23"
  file$chr[file$chr == "Y"] <- "24"
  file$chr <- suppressWarnings(as.numeric(file$chr))
  file <- file[!is.na(file$chr), , drop = FALSE]
  if (nrow(file) == 0) stop("No valid chromosomes found in sample '", sample_id, "'.")

  # If BAF column is absent, set to 0.5 (neutral value, does not affect analysis)
  if (is.null(file$BAF)) file$BAF <- 0.5

  l3 <- get_cytobands_data("level3")
  l4 <- get_cytobands_data("level4")

  if (acrocentric) {
    acro_chrs <- c(13, 14, 15, 21, 22)
    for (chr in acro_chrs) {
      idx_l3 <- which(l3$chr == chr & l3$label == "p")
      if (length(idx_l3) > 0) {
        idx_l4 <- which(l4$chr == chr)
        if (length(idx_l4) > 0) {
          l4$length[idx_l4] <- l4$length[idx_l4] - l3$length[idx_l3]
        }
      }
    }
  }

  # ===== Purity Correction (disabled) =====
  # The original CNApp applied purity correction using a minimum floor of 0.4.
  # This floor is inappropriate for liquid biopsy (tumor fraction 5-15%) and
  # would overcorrect seg.mean values in that context. Users requiring purity
  # correction should apply it externally (e.g. using ichorCNA tumor fraction
  # estimates) before passing segments to resegment_sample().
  #
  # r_lim <- 0.4
  # if (!is.null(file$purity) && !is.na(file$purity[1])) {
  #   r <- file$purity[1]
  #   if (r < r_lim) r <- r_lim
  #   v <- file$seg.mean
  #   file$seg.mean <- sapply(v, function(n, r) {
  #     inside_log <- (2^n + (r - 1)) / r
  #     if (inside_log <= 2^log2((1/2) * r_lim)) inside_log <- 2^log2((1/2) * r_lim)
  #     log2(inside_log)
  #   }, r = r)
  #   file$seg.mean <- ifelse(is.na(file$seg.mean), 0, file$seg.mean)
  # }

  loss_lim <- log2(0.2)  # cap floor: log2((1/2) * 0.4)

  # ===== Cap seg.mean always =====
  file$seg.mean <- sapply(file$seg.mean, function(n) {
    if (is.na(n)) n <- 0
    if (n < loss_lim) loss_lim else n
  })

  # ===== Re-segmentation =====
  if (!skip_resegmentation) {
    filt <- .resegment_iterative(
      file,
      min_length = min_length,
      max_dist_segm = max_dist_segm,
      percent_dist = percent_dist,
      dev_btw_segs = dev_btw_segs,
      dev_tozero = dev_tozero,
      dev_baf = dev_baf
    )
  } else {
    filt <- file
    filt <- filt[filt$length > min_length, ]
    if (nrow(filt) == 0) stop("Sample '", sample_id, "' has no segments after minimum length filtering.")
  }

  # ===== Classify CNAs =====
  filt <- .classify_cna_segments(
    filt,
    l3 = l3, l4 = l4,
    low_gain = low_gain, normal_gain = normal_gain, high_gain = high_gain,
    low_loss = low_loss, normal_loss = normal_loss, high_loss = high_loss,
    chrom_percent = chrom_percent, arm_percent = arm_percent,
    focal_percent_low = focal_percent_low,
    focal_percent_medium = focal_percent_medium,
    focal_percent_high = focal_percent_high,
    min_baf = min_baf
  )

  return(filt)

  }, error = function(e) {
    stop("Error processing sample '", sample_id, "': ", conditionMessage(e),
         "\nCheck that the input data has the correct format (columns: ID, chr, loc.start, loc.end, seg.mean).")
  })
}


#' @keywords internal
.resegment_iterative <- function(file, min_length, max_dist_segm, percent_dist,
                                 dev_btw_segs, dev_tozero, dev_baf) {
  n_loops <- 100

  filt <- file[file$length > min_length, , drop = FALSE]
  if (nrow(filt) == 0) stop("No segments remain after minimum length filtering.")

  near_zero <- filt$seg.mean > -dev_tozero & filt$seg.mean < dev_tozero
  filt$seg.mean[near_zero] <- 0

  # Remove chromosome Y (chr==24) but keep X to avoid losing relevant CNAs
  filt <- filt[filt$chr != 24, , drop = FALSE]
  if (nrow(filt) == 0) stop("No segments remain after Y chromosome removal.")

  filt <- filt[order(filt$chr, filt$loc.start), , drop = FALSE]
  rownames(filt) <- NULL

  prev_nrow <- nrow(filt)

  for (zz in 1:n_loops) {
    merged_any <- FALSE
    result_list <- list()
    i <- 1

    while (i <= nrow(filt)) {
      if (i == nrow(filt)) {
        result_list[[length(result_list) + 1]] <- filt[i, , drop = FALSE]
        i <- i + 1
      } else {
        s1 <- filt[i, , drop = FALSE]
        s2 <- filt[i + 1, , drop = FALSE]
        should_merge <- FALSE

        if (s1$chr == s2$chr) {
          dist     <- s2$loc.start - s1$loc.end
          mean_diff <- abs(s1$seg.mean - s2$seg.mean)
          baf_diff  <- abs(s1$BAF - s2$BAF)
          merge_dist <- dist < max_dist_segm
          merge_dev  <- mean_diff < dev_btw_segs
          merge_baf  <- baf_diff < dev_baf
          # Gap between segments must be less than percent_dist * combined length
          merge_pct  <- dist < (s1$length + s2$length) * percent_dist
          should_merge <- merge_dist && merge_dev && merge_baf && merge_pct
        }

        if (should_merge) {
          merged_segment <- s1
          merged_segment$loc.end <- s2$loc.end
          w1 <- s1$length
          w2 <- s2$length
          merged_segment$seg.mean <- (s1$seg.mean * w1 + s2$seg.mean * w2) /
            (w1 + w2)
          merged_segment$BAF <- (s1$BAF * w1 + s2$BAF * w2) / (w1 + w2)
          # Length is recalculated after each full pass, not during merging,
          # to avoid inflating the merge_pct denominator in subsequent steps
          result_list[[length(result_list) + 1]] <- merged_segment
          merged_any <- TRUE
          i <- i + 2
        } else {
          result_list[[length(result_list) + 1]] <- s1
          i <- i + 1
        }
      }
    }

    if (length(result_list) > 0) {
      filt <- as.data.frame(do.call(rbind, result_list), stringsAsFactors = FALSE)
      filt$chr       <- as.integer(filt$chr)
      filt$loc.start <- as.integer(filt$loc.start)
      filt$loc.end   <- as.integer(filt$loc.end)
      filt$seg.mean  <- as.numeric(filt$seg.mean)
      filt$BAF       <- as.numeric(filt$BAF)
      if (!is.null(filt$purity)) filt$purity <- as.numeric(filt$purity)
      # Recalculate length after each full pass (matching original behaviour)
      filt$length    <- filt$loc.end - filt$loc.start
      rownames(filt) <- NULL
    }

    curr_nrow <- nrow(filt)
    if (!merged_any || curr_nrow >= prev_nrow) break
    prev_nrow <- curr_nrow
  }

  filt$length <- filt$loc.end - filt$loc.start
  return(filt)
}


#' @keywords internal
.classify_cna_segments <- function(filt, l3, l4,
                                   low_gain, normal_gain, high_gain,
                                   low_loss, normal_loss, high_loss,
                                   chrom_percent, arm_percent,
                                   focal_percent_low, focal_percent_medium, focal_percent_high,
                                   min_baf) {

  filt$score      <- NA_integer_
  filt$classified <- NA_character_
  filt$type       <- NA_character_
  filt$intensity  <- NA_character_
  filt$weight     <- NA_integer_
  filt$comments   <- NA_character_

  n_cna   <- 0
  n_arm   <- 0
  n_focal <- 0

  for (i in seq_len(nrow(filt))) {
    chr_val <- gsub("chr", "", filt$chr[i])
    if (chr_val %in% c("Y", "24")) next
    chr <- suppressWarnings(as.numeric(chr_val))
    if (is.na(chr)) next

    chr_len <- l4[l4$chr == chr, "length"]
    if (length(chr_len) == 0) { chr_len <- 1e9 } else { chr_len <- chr_len[1] }

    # ===== CHROMOSOMAL-LEVEL =====
    if (filt$length[i] > chrom_percent * chr_len) {

      w <- 0
      # Gain/loss first, then CN-LOH overwrites if applicable
      w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain, high_gain, low_loss, normal_loss, high_loss)
      if (!is.na(filt$BAF[i]) &&
          filt$seg.mean[i] >= low_loss && filt$seg.mean[i] <= low_gain &&
          (filt$BAF[i] > (0.5 + min_baf) | filt$BAF[i] < (0.5 - min_baf))) {
        w <- 2
        filt$type[i] <- "CN-LOH"
      }

      if (w != 0) {
        filt$classified[i] <- "chromosomal"
        if (is.na(filt$type[i])) {
          filt$type[i] <- ifelse(filt$seg.mean[i] > 0, "Gain", "Loss")
        }
        filt$intensity[i] <- .classify_intensity(w)
        filt$score[i] <- w
        n_cna <- n_cna + w
      }
      next
    }

    # ===== ARM-LEVEL =====
    arms <- l3[l3$chr == chr, ]
    if (nrow(arms) >= 2) {
      centromere <- arms[1, "end"]
      l_p <- (min(c(centromere, filt$loc.end[i])) - filt$loc.start[i]) / arms[1, "length"]
      l_q <- (filt$loc.end[i] - max(c(centromere, filt$loc.start[i]))) / arms[2, "length"]

      if (l_p > arm_percent && l_q > arm_percent) {
        filt$comments[i] <- "Both arms significant"
      }

      if (l_p > arm_percent) {
        w <- 0
        # Gain/loss first, then CN-LOH overwrites if applicable
        w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain,
                          high_gain, low_loss, normal_loss, high_loss)
        if (!is.na(filt$BAF[i]) &&
            filt$seg.mean[i] >= low_loss && filt$seg.mean[i] <= low_gain &&
            (filt$BAF[i] > (0.5 + min_baf) | filt$BAF[i] < (0.5 - min_baf))) {
          w <- 2
          filt$type[i] <- "CN-LOH"
        }

        if (w != 0) {
          filt$classified[i] <- "arm"
          if (is.na(filt$type[i])) {
            filt$type[i] <- ifelse(filt$seg.mean[i] > 0, "Gain", "Loss")
          }
          filt$intensity[i] <- .classify_intensity(w)
          filt$score[i]     <- w
          # Arm-level alterations do not have a weight (only focal ones do)
          n_arm <- n_arm + w
        }
        next
      }

      if (l_q > arm_percent) {
        w <- 0
        # Gain/loss first, then CN-LOH overwrites if applicable
        w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain,
                          high_gain, low_loss, normal_loss, high_loss)
        if (!is.na(filt$BAF[i]) &&
            filt$seg.mean[i] >= low_loss && filt$seg.mean[i] <= low_gain &&
            (filt$BAF[i] > (0.5 + min_baf) | filt$BAF[i] < (0.5 - min_baf))) {
          w <- 2
          filt$type[i] <- "CN-LOH"
        }

        if (w != 0) {
          filt$classified[i] <- "arm"
          if (is.na(filt$type[i])) {
            filt$type[i] <- ifelse(filt$seg.mean[i] > 0, "Gain", "Loss")
          }
          filt$intensity[i] <- .classify_intensity(w)
          filt$score[i]     <- w
          # Arm-level alterations do not have a weight (only focal ones do)
          n_arm <- n_arm + w
        }
        next
      }
    }

    # ===== FOCAL-LEVEL =====
    arms <- l3[l3$chr == chr, ]
    if (nrow(arms) >= 2) {
      centromere <- arms[1, "end"]
      l_p <- (min(c(centromere, filt$loc.end[i])) - filt$loc.start[i]) / arms[1, "length"]
      l_q <- (filt$loc.end[i] - max(c(centromere, filt$loc.start[i]))) / arms[2, "length"]
      l_tot <- max(l_p, l_q)

      if (l_p > 0 && l_q > 0) {
        filt$comments[i] <- paste0("Spans p and q arms: p=", round(l_p, 3), " q=", round(l_q, 3))
      }

      ww <- 0
      if (l_tot <= focal_percent_low) ww <- 1
      if (l_tot <= focal_percent_medium && l_tot > focal_percent_low) ww <- 2
      if (l_tot <= focal_percent_high && l_tot > focal_percent_medium) ww <- 3
      if (l_tot > focal_percent_high) ww <- 4

      w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain, high_gain, low_loss, normal_loss, high_loss)

      if (w != 0) {
        filt$classified[i] <- "focal"
        filt$type[i] <- ifelse(filt$seg.mean[i] > 0, "Gain", "Loss")
        filt$intensity[i] <- .classify_intensity(w)
        filt$score[i] <- w * ww
        filt$weight[i] <- ww
        n_focal <- n_focal + w * ww
        next
      }

      # ===== CN-LOH focal =====
      if (!is.na(filt$BAF[i])) {
        w <- 0
        if (filt$seg.mean[i] >= low_loss && filt$seg.mean[i] <= low_gain &&
            (filt$BAF[i] > (0.5 + min_baf) | filt$BAF[i] < (0.5 - min_baf))) {
          w <- 2
        }
        if (w != 0) {
          filt$classified[i] <- "focal"
          filt$type[i] <- "CN-LOH"
          filt$score[i] <- w * ww
          filt$intensity[i] <- .classify_intensity(w)
          filt$weight[i] <- ww
          n_focal <- n_focal + w * ww
        }
      }
    }

    # If no condition classified this segment, mark it as diploid
    if (is.na(filt$classified[i])) {
      filt$classified[i] <- "diploid"
      filt$type[i]       <- "diploid"
      filt$score[i]      <- 0
      filt$intensity[i]  <- "None"
    }

  } # cierra for

  # Segments that were skipped early (e.g. Y chr, invalid chr) remain NA —
  # they are regions smoothed back to diploid by the resegmentation algorithm.
  na_cls <- is.na(filt$classified)
  if (any(na_cls)) {
    filt$classified[na_cls] <- "diploid"
    filt$type[na_cls]       <- "diploid"
    filt$score[na_cls]      <- 0L
    filt$intensity[na_cls]  <- "None"
    filt$weight[na_cls]     <- 0L
  }

  attr(filt, "n_chromosomal") <- n_cna
  attr(filt, "n_arm") <- n_arm
  attr(filt, "n_focal") <- n_focal
  return(filt)
} # cierra .classify_cna_segments


#' @keywords internal
.classify_intensity <- function(w) {
  if (w == 1) return("Low")
  if (w == 2) return("Medium")
  if (w == 3) return("High")
  return(NA_character_)
}


#' Classify Single CNA Based on seg.mean
#'
#' @param seg_mean Numeric log2 ratio.
#' @param low_gain Numeric low gain cutoff.
#' @param medium_gain Numeric medium gain cutoff.
#' @param high_gain Numeric high gain cutoff.
#' @param low_loss Numeric low loss cutoff.
#' @param medium_loss Numeric medium loss cutoff.
#' @param high_loss Numeric high loss cutoff.
#'
#' @return Integer (0=no alteration, 1=low, 2=medium, 3=high).
#'
#' @export
classify_cna <- function(seg_mean, low_gain, medium_gain, high_gain,
                         low_loss, medium_loss, high_loss) {
  w <- 0
  if (seg_mean >= low_gain && seg_mean < medium_gain) w <- 1
  if (seg_mean >= medium_gain && seg_mean < high_gain) w <- 2
  if (seg_mean >= high_gain) w <- 3
  if (seg_mean <= low_loss && seg_mean > medium_loss) w <- 1
  if (seg_mean <= medium_loss && seg_mean > high_loss) w <- 2
  if (seg_mean <= high_loss) w <- 3
  return(w)
}
