#' Re-segment Copy Number Data
#'
#' Performs iterative re-segmentation, fusion of similar segments, and classification
#' of copy number alterations (CNA).
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
#' @param normal_gain Numeric normal gain threshold (default log2(3/2) ≈ 0.585).
#' @param high_gain Numeric high gain threshold (default log2(4/2) = 1.0).
#' @param low_loss Numeric low loss threshold (default -0.2).
#' @param normal_loss Numeric normal loss threshold (default log2(1/2) ≈ -1.0).
#' @param high_loss Numeric high loss threshold (default log2(0.6/2) ≈ -1.737).
#' @param chrom_percent Numeric minimum coverage for chromosomal alterations (0-1, default 0.9).
#' @param arm_percent Numeric minimum coverage for arm-level alterations (0-1, default 0.5).
#' @param focal_percent_low Numeric focal alteration low threshold (default 0.05).
#' @param focal_percent_medium Numeric focal alteration medium threshold (default 0.15).
#' @param focal_percent_high Numeric focal alteration high threshold (default 0.3).
#' @param min_baf Numeric minimum BAF for CN-LOH detection (default 0.2).
#' @param acrocentric Logical whether to ignore p arms of chromosomes 13,14,15,21,22 (default FALSE).
#' @param skip_resegmentation Logical whether to skip re-segmentation (default FALSE).
#'
#' @return Data frame (re-segmented) with additional columns:
#'   - score, classified, type, intensity, weight, comments
#'
#' @details
#' This is the main CNA processing pipeline:
#' 1. Purity correction (re-centralization if purity provided)
#' 2. Minimum length filtering
#' 3. Iterative segment fusion (merging similar, nearby segments)
#' 4. Classification into chromosomal, arm-level, focal, and CN-LOH events
#'
#' @examples
#' \dontrun{
#'   data <- read_cna_file("sample.tsv")
#'   result <- resegment_sample(data, sample_id = "sample_001")
#' }
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
                             normal_gain = log2(3 / 2),
                             high_gain = log2(4 / 2),
                             low_loss = -0.2,
                             normal_loss = log2(1 / 2),
                             high_loss = log2(0.6 / 2),
                             chrom_percent = 0.9,
                             arm_percent = 0.5,
                             focal_percent_low = 0.05,
                             focal_percent_medium = 0.15,
                             focal_percent_high = 0.3,
                             min_baf = 0.2,
                             acrocentric = FALSE,
                             skip_resegmentation = FALSE) {

  # Filter data for sample
  file <- data[data$ID == sample_id, ]
  if (nrow(file) == 0) {
    stop("Sample '", sample_id, "' not found in data.")
  }

  # Initialize columns
  file$length <- file$loc.end - file$loc.start
  file$BAF <- ifelse(is.null(file$BAF), 0.5, file$BAF)

  # Load cytoband reference
  l3 <- get_cytobands_data("level3")
  l4 <- get_cytobands_data("level4")

  # Adjust acrocentrics if needed
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

  # ===== Purity Correction =====
  r_lim <- 0.4
  loss_lim <- log2((1 / 2) * r_lim)

  if (!is.null(file$purity) && !is.na(file$purity[1])) {
    r <- file$purity[1]
    if (r < r_lim) r <- r_lim

    v <- file$seg.mean
    file$seg.mean <- sapply(v, function(n, r) {
      inside_log <- (2^n + (r - 1)) / r
      if (inside_log <= 2^loss_lim) inside_log <- 2^loss_lim
      log2(inside_log)
    }, r = r)
  }

  # Cap seg.mean values when NO purity
  file$seg.mean <- pmax(file$seg.mean, loss_lim)

  # ===== Skip or Perform Re-segmentation =====
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
    # Still remove segments below min length if skipping
    filt <- filt[filt$length > min_length, ]
    if (nrow(filt) == 0) {
      stop("Sample '", sample_id, "' has no segments after minimum length filtering.")
    }
  }

  # ===== Classify CNAs =====
  filt <- .classify_cna_segments(
    filt,
    l3 = l3,
    l4 = l4,
    low_gain = low_gain,
    normal_gain = normal_gain,
    high_gain = high_gain,
    low_loss = low_loss,
    normal_loss = normal_loss,
    high_loss = high_loss,
    chrom_percent = chrom_percent,
    arm_percent = arm_percent,
    focal_percent_low = focal_percent_low,
    focal_percent_medium = focal_percent_medium,
    focal_percent_high = focal_percent_high,
    min_baf = min_baf
  )

  return(filt)
}


#' Internal: Iterative Re-segmentation
#'
#' @keywords internal
.resegment_iterative <- function(file,
                                 min_length,
                                 max_dist_segm,
                                 percent_dist,
                                 dev_btw_segs,
                                 dev_tozero,
                                 dev_baf) {
  n_loops <- 100

  # Step 1: Filter by minimum length
  new_file <- list()
  new_file[[1]] <- file[file$length > min_length, ]

  if (nrow(new_file[[1]]) == 0) {
    stop("No segments remain after minimum length filtering.")
  }

  # Step 2: Set near-zero segments to 0
  near_zero <- new_file[[1]]$seg.mean > -dev_tozero & new_file[[1]]$seg.mean < dev_tozero
  new_file[[1]]$seg.mean[near_zero] <- 0

  # Step 3: Remove sex chromosomes (keep only 1-22)
  new_file[[1]] <- new_file[[1]][new_file[[1]]$chr != 24 & new_file[[1]]$chr != 23, ]

  # Step 4: Iterative fusion of segments
  for (zz in 2:n_loops) {
    new_file[[zz]] <- data.frame(
      ID = NA_character_, chr = NA_integer_, loc.start = NA_integer_,
      loc.end = NA_integer_, seg.mean = NA_real_, length = NA_integer_,
      BAF = NA_real_, stringsAsFactors = FALSE
    )

    k <- 1
    flag <- 0

    for (i in 1:(nrow(new_file[[zz - 1]]) - 1)) {
      if (flag == 1) {
        flag <- 0
        next
      }

      s1 <- new_file[[zz - 1]][i, ]
      new_file[[zz]][k, ] <- s1
      k <- k + 1

      s2 <- new_file[[zz - 1]][i + 1, ]

      # Check if segments should be merged
      if (s1$chr == s2$chr) {
        dist <- s2$loc.start - s1$loc.end
        mean_diff <- abs(s1$seg.mean - s2$seg.mean)
        baf_diff <- abs(s1$BAF - s2$BAF)

        # Condition for merging
        merge_dist <- dist < max_dist_segm
        merge_dev <- mean_diff < dev_btw_segs
        merge_baf <- baf_diff < dev_baf
        merge_pct <- (dist / (s1$length + s2$length)) * 100 < percent_dist

        if (merge_dist && merge_dev && merge_baf && merge_pct) {
          # Merge segments
          new_file[[zz]][k, "loc.end"] <- s2$loc.end
          new_file[[zz]][k, "seg.mean"] <- mean(c(s1$seg.mean, s2$seg.mean))
          new_file[[zz]][k, "BAF"] <- mean(c(s1$BAF, s2$BAF))
          new_file[[zz]][k, "length"] <- new_file[[zz]][k, "loc.end"] - new_file[[zz]][k, "loc.start"]
          flag <- 1
        }
      }
    }

    # Add last segment if not merged
    if (flag == 0) {
      new_file[[zz]][k, ] <- new_file[[zz - 1]][nrow(new_file[[zz - 1]]), ]
    }

    # Check convergence (no changes from previous iteration)
    if (zz > 2) {
      if (nrow(new_file[[zz]]) == nrow(new_file[[zz - 1]])) {
        # Check if data is identical (convergence)
        if (all.equal(new_file[[zz]], new_file[[zz - 1]]) == TRUE) {
          break
        }
      }
    }
  }

  filt <- new_file[[length(new_file)]]
  filt <- filt[!is.na(filt$ID), ]
  filt$length <- filt$loc.end - filt$loc.start

  return(filt)
}


#' Internal: Classify CNA Segments
#'
#' @keywords internal
.classify_cna_segments <- function(filt,
                                   l3,
                                   l4,
                                   low_gain,
                                   normal_gain,
                                   high_gain,
                                   low_loss,
                                   normal_loss,
                                   high_loss,
                                   chrom_percent,
                                   arm_percent,
                                   focal_percent_low,
                                   focal_percent_medium,
                                   focal_percent_high,
                                   min_baf) {

  # Initialize classification columns
  filt$score <- NA_integer_
  filt$classified <- NA_character_
  filt$type <- NA_character_
  filt$intensity <- NA_character_
  filt$weight <- NA_integer_
  filt$comments <- NA_character_

  n_cna <- 0
  n_arm <- 0
  n_focal <- 0

  for (i in seq_len(nrow(filt))) {
    chr <- as.numeric(gsub("chr", "", filt$chr[i]))

    # Get chromosome length from l4
    chr_len <- l4[l4$chr == chr, "length"]
    if (length(chr_len) == 0) chr_len <- 1e9 # fallback

    # ===== CHROMOSOMAL-LEVEL =====
    if (filt$length[i] > chrom_percent * chr_len) {
      w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain, high_gain, low_loss, normal_loss, high_loss)

      # CN-LOH check
      if (!is.na(filt$BAF[i]) && filt$BAF[i] > (1 - min_baf) && filt$BAF[i] < min_baf) {
        filt$type[i] <- "CN-LOH"
      }

      if (w != 0) {
        filt$classified[i] <- "Chromosomal"
        filt$type[i] <- ifelse(filt$seg.mean[i] > 0, "Gain", "Loss")
        filt$intensity[i] <- .classify_intensity(w)
        filt$score[i] <- w
        n_cna <- n_cna + 1
      }
      next
    }

    # ===== ARM-LEVEL =====
    arms <- l3[l3$chr == chr, ]
    if (nrow(arms) >= 2) {
      centromere <- arms[1, "end"]

      # p-arm coverage
      l_p <- (min(c(centromere, filt$loc.end[i])) - filt$loc.start[i]) / arms[1, "length"]

      # q-arm coverage
      l_q <- (filt$loc.end[i] - max(c(centromere, filt$loc.start[i]))) / arms[2, "length"]

      if (l_p > arm_percent && l_q > arm_percent) {
        filt$comments[i] <- "Both arms significant"
      }

      if (l_p > arm_percent) {
        w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain, high_gain, low_loss, normal_loss, high_loss)

        if (!is.na(filt$BAF[i]) && filt$BAF[i] > (1 - min_baf) && filt$BAF[i] < min_baf) {
          filt$type[i] <- "CN-LOH"
        }

        if (w != 0) {
          filt$classified[i] <- "Arm"
          filt$type[i] <- ifelse(is.na(filt$type[i]), ifelse(filt$seg.mean[i] > 0, "Gain", "Loss"), filt$type[i])
          filt$intensity[i] <- .classify_intensity(w)
          filt$score[i] <- w
          filt$weight[i] <- 1 # p-arm
          n_arm <- n_arm + 1
        }
        next
      }

      if (l_q > arm_percent) {
        w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain, high_gain, low_loss, normal_loss, high_loss)

        if (!is.na(filt$BAF[i]) && filt$BAF[i] > (1 - min_baf) && filt$BAF[i] < min_baf) {
          filt$type[i] <- "CN-LOH"
        }

        if (w != 0) {
          filt$classified[i] <- "Arm"
          filt$type[i] <- ifelse(is.na(filt$type[i]), ifelse(filt$seg.mean[i] > 0, "Gain", "Loss"), filt$type[i])
          filt$intensity[i] <- .classify_intensity(w)
          filt$score[i] <- w
          filt$weight[i] <- 2 # q-arm
          n_arm <- n_arm + 1
        }
        next
      }
    }

    # ===== FOCAL-LEVEL =====
    if (filt$seg.mean[i] > low_gain || filt$seg.mean[i] < low_loss) {
      arms <- l3[l3$chr == chr, ]
      if (nrow(arms) >= 2) {
        centromere <- arms[1, "end"]
        l_p <- (min(c(centromere, filt$loc.end[i])) - filt$loc.start[i]) / arms[1, "length"]
        l_q <- (filt$loc.end[i] - max(c(centromere, filt$loc.start[i]))) / arms[2, "length"]
        l_tot <- max(l_p, l_q)

        if (l_p > 0 && l_q > 0) {
          filt$comments[i] <- paste0("Spans p and q arms: p=", round(l_p, 3), " q=", round(l_q, 3))
        }

        # Weight focal alterations
        ww <- 0
        if (l_tot <= focal_percent_low) ww <- 1
        if (l_tot <= focal_percent_medium && l_tot > focal_percent_low) ww <- 2
        if (l_tot <= focal_percent_high && l_tot > focal_percent_medium) ww <- 3
        if (l_tot > focal_percent_high) ww <- 4

        w <- classify_cna(filt$seg.mean[i], low_gain, normal_gain, high_gain, low_loss, normal_loss, high_loss)

        if (w != 0) {
          filt$classified[i] <- "Focal"
          filt$type[i] <- ifelse(filt$seg.mean[i] > 0, "Gain", "Loss")
          filt$intensity[i] <- .classify_intensity(w)
          filt$score[i] <- w
          filt$weight[i] <- ww
          n_focal <- n_focal + 1
        }
      }
    }
  }

  # Store summary counts as attributes
  attr(filt, "n_chromosomal") <- n_cna
  attr(filt, "n_arm") <- n_arm
  attr(filt, "n_focal") <- n_focal

  return(filt)
}


#' Classify CNA by Intensity
#'
#' Converts numeric intensity code (1-3) to string label.
#'
#' @param w Numeric intensity (1=low, 2=medium, 3=high).
#'
#' @return Character intensity label.
#'
#' @keywords internal
.classify_intensity <- function(w) {
  if (w == 1) return("Low")
  if (w == 2) return("Medium")
  if (w == 3) return("High")
  return(NA_character_)
}


#' Classify Single CNA Based on seg.mean
#'
#' Internal function that classifies a segment as gain/loss and returns intensity.
#'
#' @param seg_mean Numeric log2 ratio.
#' @param low_gain Numeric low gain cutoff.
#' @param medium_gain Numeric medium gain cutoff.
#' @param high_gain Numeric high gain cutoff.
#' @param low_loss Numeric low loss cutoff.
#' @param medium_loss Numeric medium loss cutoff.
#' @param high_loss Numeric high loss cutoff.
#'
#' @return Integer (0=no alteration, 1=low intensity, 2=medium, 3=high).
#'
#' @export
classify_cna <- function(seg_mean,
                         low_gain,
                         medium_gain,
                         high_gain,
                         low_loss,
                         medium_loss,
                         high_loss) {
  w <- 0

  # Classify gains
  if (seg_mean >= low_gain && seg_mean < medium_gain) w <- 1
  if (seg_mean >= medium_gain && seg_mean < high_gain) w <- 2
  if (seg_mean >= high_gain) w <- 3

  # Classify losses
  if (seg_mean <= low_loss && seg_mean > medium_loss) w <- 1
  if (seg_mean <= medium_loss && seg_mean > high_loss) w <- 2
  if (seg_mean <= high_loss) w <- 3

  return(w)
}
