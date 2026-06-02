#' @importFrom rlang .data
NULL

#' Plot Segmentation Before and After Re-segmentation
#'
#' Visualizes original segments and re-segmented results on a genome-wide plot.
#'
#' @param original_data Data frame with original segments (ID, chr, loc.start, loc.end, seg.mean).
#' @param resegmented_data Data frame with re-segmented segments (same columns).
#' @param sample_id Character sample identifier (for title).
#' @param cytobands Level3 cytoband reference (from get_cytobands_data()).
#' @param centromeres Data frame with centromere positions (chr, position).
#' @param thresholds List or numeric vector with gain/loss thresholds:
#'   - low_gain, medium_gain, high_gain, low_loss, medium_loss, high_loss
#'
#' @return A ggplot2 object.
#'
#'
#' @details
#' Creates a genome-wide plot with:
#' - Chromosomes separated by vertical lines
#' - Centromeres marked with dashed lines
#' - Horizontal lines for gain/loss thresholds
#' - Segments colored by type (original=black, resegmented=green)
#'
#' @examples
#' \dontrun{
#'   plot <- plot_segmentation(
#'     original = orig_list[[1]],
#'     resegmented = filt_list[[1]],
#'     sample_id = "sample_001",
#'     cytobands = l3
#'   )
#'   print(plot)
#' }
#'
#' @export
plot_segmentation <- function(original_data,
                              resegmented_data,
                              sample_id = "Sample",
                              cytobands = NULL,
                              centromeres = NULL,
                              thresholds = NULL) {

  if (is.null(cytobands)) {
    cytobands <- get_cytobands_data("level3")
  }

  # Prepare data for plotting
  l4 <- get_cytobands_data("level4")

  # Add cumulative positions for genome-wide plot
  l4_sort <- l4[order(l4$chr), ]
  l4_sort$cum <- cumsum(l4_sort$length)
  l4_sort$to_sum <- l4_sort$cum - l4_sort$length

  # Function to add adjusted coordinates
  adjust_coords <- function(data, ref) {
    if (is.null(data) || nrow(data) == 0) return(data)
    data$adj_start <- rep(NA_real_, nrow(data))

    for (i in seq_len(nrow(data))) {
      chr_code <- as.numeric(data$chr[i])
      chr_match <- which(ref$chr == chr_code)

      if (length(chr_match) > 0) {
        data$adj_start[i] <- data$loc.start[i] + ref$to_sum[chr_match]
      }
    }

    data$adj_end <- data$adj_start + (data$loc.end - data$loc.start)
    return(data)
  }

  # Adjust coordinates
  original_data <- adjust_coords(original_data, l4_sort)
  resegmented_data <- adjust_coords(resegmented_data, l4_sort)

  # Remove NA rows (unmapped segments)
  original_data <- original_data[!is.na(original_data$adj_start), ]
  resegmented_data <- resegmented_data[!is.na(resegmented_data$adj_start), ]

  # Create ggplot
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = original_data,
      ggplot2::aes(x = .data$adj_start, xend = .data$adj_end, y = .data$seg.mean, yend = .data$seg.mean),
      color = "black", linewidth = 0.7, alpha = 0.6
    ) +
    ggplot2::geom_segment(
      data = resegmented_data,
      ggplot2::aes(x = .data$adj_start, xend = .data$adj_end, y = .data$seg.mean, yend = .data$seg.mean),
      color = "forestgreen", linewidth = 1.0, alpha = 0.8
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray30", linewidth = 0.5) +
    ggplot2::facet_wrap(~factor(1), ncol = 1) +
    ggplot2::labs(
      title = paste("Segmentation Comparison:", sample_id),
      x = "Genomic Position", y = "log2(ratio)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "lightgray", linewidth = 0.3),
      axis.text.x = ggplot2::element_blank()
    )

  # Add threshold lines if provided
  if (!is.null(thresholds)) {
    threshold_values <- c(
      thresholds[["low_gain"]], thresholds[["normal_gain"]], thresholds[["high_gain"]],
      thresholds[["low_loss"]], thresholds[["normal_loss"]], thresholds[["high_loss"]]
    )
    threshold_values <- threshold_values[!is.na(threshold_values)]

    p <- p +
      ggplot2::geom_hline(
        yintercept = threshold_values,
        linetype = "dashed", color = "orange", linewidth = 0.4, alpha = 0.7
      )
  }

  return(p)
}


#' Plot Copy Number Frequency Heatmap
#'
#' Creates a frequency heatmap showing gains/losses across samples by chromosome.
#'
#' @param resegmented_list List of re-segmented data frames (one per sample).
#' @param sample_names Character vector of sample IDs.
#' @param gain_threshold Numeric seg.mean cutoff for gains (default 0.23).
#' @param loss_threshold Numeric seg.mean cutoff for losses (default -0.23).
#' @param group_var (optional) Categorical variable for stratification.
#' @param group_values (optional) Subset of group_var values to plot.
#'
#' @return A ggplot2 object (heatmap).
#'
#' @details
#' Counts gains/losses per chromosome across samples and displays as heatmap.
#' Can be stratified by clinical group.
#'
#' @examples
#' \dontrun{
#'   plot <- plot_cn_frequency(filt_list, sample_names = names(filt_list))
#'   print(plot)
#' }
#'
#' @export
plot_cn_frequency <- function(resegmented_list,
                              sample_names = NULL,
                              gain_threshold = 0.23,
                              loss_threshold = -0.23,
                              group_var = NULL,
                              group_values = NULL) {

  if (is.null(sample_names)) {
    sample_names <- names(resegmented_list)
    if (is.null(sample_names)) {
      sample_names <- paste0("sample_", seq_along(resegmented_list))
    }
  }

  # Count gains/losses per chromosome per sample
  freq_matrix <- matrix(0, nrow = 24, ncol = length(sample_names))
  rownames(freq_matrix) <- 1:24
  colnames(freq_matrix) <- sample_names

  for (i in seq_along(resegmented_list)) {
    reseg <- resegmented_list[[i]]
    if (is.null(reseg) || nrow(reseg) == 0L) next

    chr_num <- suppressWarnings(
      as.integer(gsub("^chr", "", as.character(reseg$chr)))
    )
    valid <- !is.na(chr_num) & chr_num >= 1L & chr_num <= 24L
    if (!any(valid)) next

    chr_v <- chr_num[valid]
    seg_v <- reseg$seg.mean[valid]

    gain_idx <- which(seg_v > gain_threshold)
    loss_idx <- which(seg_v < loss_threshold)

    if (length(gain_idx) > 0L)
      freq_matrix[, i] <- freq_matrix[, i] +
        tabulate(chr_v[gain_idx], nbins = 24L)
    if (length(loss_idx) > 0L)
      freq_matrix[, i] <- freq_matrix[, i] -
        tabulate(chr_v[loss_idx], nbins = 24L)
  }

  # Convert to long format for ggplot
  freq_df <- as.data.frame(as.table(freq_matrix))
  colnames(freq_df) <- c("chr", "sample", "frequency")
  freq_df$chr <- as.numeric(as.character(freq_df$chr))
  freq_df$sample <- as.character(freq_df$sample)

  # Create heatmap
  p <- ggplot2::ggplot(freq_df, ggplot2::aes(x = .data$sample, y = .data$chr, fill = .data$frequency)) +
    ggplot2::geom_tile(color = "white", size = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", name = "Gain/Loss Freq"
    ) +
    ggplot2::labs(x = "Sample", y = "Chromosome", title = "Copy Number Frequency") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  return(p)
}


#' Plot Genome-Wide CNA Profile (ASCAT-style)
#'
#' Generates a genome-wide CNA visualization showing per-bin log2 ratio dots
#' colored by CNA call (NEUT=blue, HETD=green, GAIN/AMP/HLAMP=red) with
#' segment mean horizontal lines superimposed. Replicates the style produced
#' by ASCAT/ichorCNA genomeWide plots.
#'
#' @param bins_data Data frame from \code{read_bam_qc()}. Required columns:
#'   \code{chr}, \code{loc.start}, \code{loc.end}, \code{log2_ratio}.
#' @param seg_data Data frame of segments. Required columns: \code{chr}
#'   (integer or character), \code{loc.start}, \code{loc.end}, \code{call}.
#'   The Y value is taken from \code{seg.mean} if present, otherwise
#'   \code{seg.median.logR}.
#' @param sample_id Character sample label shown in the title.
#' @param tumor_fraction Numeric tumor fraction (optional, shown in subtitle).
#' @param ploidy Numeric ploidy (optional, shown in subtitle).
#' @param ylim Numeric vector of length 2 for Y-axis limits (default
#'   \code{c(-2, 2)}).
#' @param call_colors Named character vector to override default call colors.
#'   Recognized names: NEUT, HETD, HOMD, GAIN, AMP, HLAMP, HLAMP2, HLAMP3.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#'   bins <- read_bam_qc("sample.bam", bin_size = 1000000,
#'                        sequencing_type = "nanopore", experiment_type = "wgs")
#'   segs <- read.table("sample.seg.txt", header = TRUE, sep = "\t")
#'   # segs must have: chr, loc.start, loc.end, seg.mean (or seg.median.logR), call
#'   p <- plot_genome_wide_cna(bins, segs, sample_id = "barcode66",
#'                              tumor_fraction = 0.48, ploidy = 1.96)
#'   print(p)
#' }
#'
#' @export
plot_genome_wide_cna <- function(bins_data,
                                  seg_data,
                                  sample_id      = "Sample",
                                  tumor_fraction = NULL,
                                  ploidy         = NULL,
                                  ylim           = c(-2, 2),
                                  call_colors    = NULL) {

  # ── color scheme ────────────────────────────────────────────────
  cna_colors <- c(
    NEUT   = "steelblue",
    HETD   = "forestgreen",
    HOMD   = "darkgreen",
    GAIN   = "red3",
    AMP    = "red3",
    HLAMP  = "red3",
    HLAMP2 = "red3",
    HLAMP3 = "red3"
  )
  if (!is.null(call_colors)) cna_colors[names(call_colors)] <- call_colors

  # ── cytoband cumulative positions (autosomes only) ───────────────
  l4 <- get_cytobands_data("level4")
  l4 <- l4[order(l4$chr), ]
  l4 <- l4[l4$chr <= 22, ]
  l4$cum_start <- c(0, cumsum(l4$length[-nrow(l4)]))
  offset_map   <- stats::setNames(l4$cum_start, as.character(l4$chr))

  add_offset <- function(chr_int, pos) {
    off <- offset_map[as.character(chr_int)]
    pos + ifelse(is.na(off), 0, off)
  }

  norm_chr <- function(x) suppressWarnings(as.integer(gsub("^chr", "", as.character(x))))

  # ── prepare bins ─────────────────────────────────────────────────
  bins_data$chr_int <- norm_chr(bins_data$chr)
  bins_data <- bins_data[!is.na(bins_data$chr_int) & bins_data$chr_int <= 22, ]
  bins_data <- bins_data[!is.na(bins_data$log2_ratio), ]
  bins_data$x_mid  <- add_offset(bins_data$chr_int,
                                  (bins_data$loc.start + bins_data$loc.end) / 2)
  # Clamp slightly outside ylim so clipped dots still appear at the boundary
  bins_data$y_plot <- pmax(pmin(bins_data$log2_ratio, ylim[2] + 0.25), ylim[1] - 0.25)

  # ── prepare segments ─────────────────────────────────────────────
  seg_data$chr_int <- norm_chr(seg_data$chr)
  seg_data <- seg_data[!is.na(seg_data$chr_int) & seg_data$chr_int <= 22, ]

  seg_mean_col  <- if ("seg.mean" %in% colnames(seg_data)) "seg.mean" else "seg.median.logR"
  seg_data$y        <- seg_data[[seg_mean_col]]
  seg_data$x_start  <- add_offset(seg_data$chr_int, seg_data$loc.start)
  seg_data$x_end    <- add_offset(seg_data$chr_int, seg_data$loc.end)

  # Derive 'call' from CNAppR 'type' column when ASCAT 'call' is absent
  if (!"call" %in% colnames(seg_data)) {
    type_col <- if ("type" %in% colnames(seg_data)) seg_data$type else rep(NA_character_, nrow(seg_data))
    seg_data$call <- "NEUT"
    seg_data$call[!is.na(type_col) & type_col == "Gain" & seg_data$y < 1.0]  <- "GAIN"
    seg_data$call[!is.na(type_col) & type_col == "Gain" & seg_data$y >= 1.0] <- "AMP"
    seg_data$call[!is.na(type_col) & type_col == "Loss" & seg_data$y > -1.5] <- "HETD"
    seg_data$call[!is.na(type_col) & type_col == "Loss" & seg_data$y <= -1.5] <- "HOMD"
  }

  # ── assign CNA call to each bin (by overlapping segment) ─────────
  bins_data$call <- "NEUT"
  for (c in unique(bins_data$chr_int)) {
    idx_b <- which(bins_data$chr_int == c)
    seg_c <- seg_data[seg_data$chr_int == c, , drop = FALSE]
    if (nrow(seg_c) == 0L) next
    for (i in seq_len(nrow(seg_c))) {
      hit <- idx_b[bins_data$loc.start[idx_b] >= seg_c$loc.start[i] &
                   bins_data$loc.end[idx_b]   <= seg_c$loc.end[i]]
      if (length(hit) > 0L) bins_data$call[hit] <- seg_c$call[i]
    }
  }

  # ── chromosome label colors (dominant call per chromosome) ────────
  chr_mids      <- l4$cum_start + l4$length / 2
  chr_label_col <- vapply(l4$chr, function(c) {
    sc <- seg_data[seg_data$chr_int == c, , drop = FALSE]
    if (nrow(sc) == 0L) return("steelblue")
    sc$seg_len <- sc$x_end - sc$x_start
    dom <- sc$call[which.max(sc$seg_len)]
    if (is.na(dom) || dom == "NEUT")                              return("steelblue")
    if (dom %in% c("GAIN", "AMP", "HLAMP", "HLAMP2", "HLAMP3")) return("red3")
    if (dom %in% c("HETD", "HOMD"))                              return("forestgreen")
    "steelblue"
  }, character(1L))

  # ── title ─────────────────────────────────────────────────────────
  sub_parts <- character(0)
  if (!is.null(tumor_fraction))
    sub_parts <- c(sub_parts, paste0("Tumor Fraction: ", round(tumor_fraction, 4)))
  if (!is.null(ploidy))
    sub_parts <- c(sub_parts, paste0("Ploidy: ", round(ploidy, 2)))
  title_str <- if (length(sub_parts) > 0L) {
    paste0(sample_id, "\n", paste(sub_parts, collapse = ", "))
  } else {
    sample_id
  }

  x_max <- max(l4$cum_start + l4$length)

  # ── build plot ───────────────────────────────────────────────────
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = bins_data,
      ggplot2::aes(x = .data$x_mid, y = .data$y_plot, color = .data$call),
      size = 0.4, alpha = 0.55, shape = 16
    ) +
    ggplot2::geom_segment(
      data = seg_data,
      ggplot2::aes(
        x     = .data$x_start, xend = .data$x_end,
        y     = .data$y,       yend = .data$y,
        color = .data$call
      ),
      linewidth = 1.1, alpha = 0.9
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
    ggplot2::geom_vline(
      xintercept = l4$cum_start,
      linetype   = "dashed", color = "gray55", linewidth = 0.25
    ) +
    ggplot2::scale_color_manual(values = cna_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max),
      labels = NULL, breaks = NULL,
      expand = c(0.005, 0)
    ) +
    ggplot2::coord_cartesian(ylim = ylim, clip = "off") +
    ggplot2::labs(title = title_str, x = NULL, y = "Copy Number (log2 ratio)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(hjust = 0.5, size = 9),
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      axis.line.x        = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 8),
      panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.25),
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(t = 5, r = 5, b = 22, l = 5)
    ) +
    ggplot2::annotate(
      "text",
      x     = chr_mids,
      y     = ylim[1] - (ylim[2] - ylim[1]) * 0.12,
      label = as.character(l4$chr),
      color = chr_label_col,
      size  = 2.8, hjust = 0.5
    )

  p
}


#' Plot Copy Number Frequency Profile (Genomic Position on X)
#'
#' Creates the classic copy-number frequency bar plot: gain frequency shown as
#' positive bars (red) and loss frequency as negative bars (blue), with
#' genomic position on the X axis. Accepts the output of
#' \code{compute_region_frequencies()}.
#'
#' @param freq_df Data frame output of \code{compute_region_frequencies()}.
#'   Required columns: \code{chr}, \code{start}, \code{end},
#'   \code{gain_freq}, \code{loss_freq} (all as proportions 0–1).
#' @param title Character plot title.
#' @param gain_color Character hex/name color for gain bars (default red).
#' @param loss_color Character hex/name color for loss bars (default blue).
#' @param alpha Numeric transparency for bars (default 0.85).
#'
#' @return A ggplot2 object.
#'
#' @details
#' The Y axis shows percentage of samples (0–100\%). Chromosome boundaries
#' are drawn as dashed vertical lines; labels are placed below the axis.
#' Works with any reference granularity: arm-level (44 bars), cytoband
#' (~800 bars), or WES targets (246k tiny bars, rendering as a smooth curve).
#'
#' @examples
#' \dontrun{
#'   arm_ref <- "inst/aux_files/segmented_files_hg19/autosomes_hg19_by_arms.txt"
#'   freq    <- compute_region_frequencies(lihc_list, arm_ref)
#'   plot_frequency_profile(freq, title = "TCGA-LIHC  |  n = 354")
#' }
#'
#' @export
plot_frequency_profile <- function(freq_df,
                                   title      = "Copy Number Frequency Profile",
                                   gain_color = "#d73027",
                                   loss_color = "#4575b4",
                                   alpha      = 0.85) {

  # ── normalise chr column ─────────────────────────────────────────────
  freq_df$chr_int <- suppressWarnings(
    as.integer(gsub("^chr", "", as.character(freq_df$chr)))
  )
  freq_df <- freq_df[!is.na(freq_df$chr_int) & freq_df$chr_int <= 22L, ]
  freq_df <- freq_df[order(freq_df$chr_int, freq_df$start), ]
  if (nrow(freq_df) == 0L) stop("No autosomal regions found in freq_df.")

  # ── cumulative genomic offsets ───────────────────────────────────────
  chrs      <- sort(unique(freq_df$chr_int))
  chr_max   <- vapply(chrs, function(c) max(freq_df$end[freq_df$chr_int == c]),
                      numeric(1L))
  cum_off   <- c(0, cumsum(chr_max[-length(chr_max)]))
  names(cum_off) <- as.character(chrs)

  freq_df$x_start <- freq_df$start + cum_off[as.character(freq_df$chr_int)]
  freq_df$x_end   <- freq_df$end   + cum_off[as.character(freq_df$chr_int)]

  # Boundaries and midpoints for chromosome labels
  chr_bounds <- vapply(chrs, function(c) max(freq_df$x_end[freq_df$chr_int == c]),
                       numeric(1L))
  chr_mids   <- vapply(chrs, function(c)
    mean(range(freq_df$x_start[freq_df$chr_int == c],
               freq_df$x_end[freq_df$chr_int == c])),
    numeric(1L))

  x_max <- max(chr_bounds)

  # ── build plot ───────────────────────────────────────────────────────
  p <- ggplot2::ggplot(freq_df) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$x_start, xmax = .data$x_end,
                   ymin = 0, ymax = .data$gain_freq * 100),
      fill = gain_color, alpha = alpha
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$x_start, xmax = .data$x_end,
                   ymin = -.data$loss_freq * 100, ymax = 0),
      fill = loss_color, alpha = alpha
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_vline(
      xintercept = c(0, chr_bounds),
      color = "gray65", linewidth = 0.25, linetype = "dashed"
    ) +
    ggplot2::annotate(
      "text",
      x     = chr_mids,
      y     = -108,
      label = as.character(chrs),
      size  = 2.6, hjust = 0.5, color = "gray25"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max),
      labels = NULL, breaks = NULL, expand = c(0.005, 0)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(abs(x), "%"),
      limits = c(-115, 105),
      breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100)
    ) +
    ggplot2::labs(title = title, x = NULL, y = "Samples (%)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(hjust = 0.5, size = 11),
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      axis.line.x        = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.25),
      plot.margin        = ggplot2::margin(t = 5, r = 5, b = 22, l = 5)
    )

  p
}
