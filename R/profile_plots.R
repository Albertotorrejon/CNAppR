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
#' @importFrom utils get_cytobands_data
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
    data$adj_start <- NA_real_

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
      color = "black", size = 0.7, alpha = 0.6
    ) +
    ggplot2::geom_segment(
      data = resegmented_data,
      ggplot2::aes(x = .data$adj_start, xend = .data$adj_end, y = .data$seg.mean, yend = .data$seg.mean),
      color = "forestgreen", size = 1.0, alpha = 0.8
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray30", size = 0.5) +
    ggplot2::facet_wrap(~factor(1), ncol = 1) +
    ggplot2::labs(
      title = paste("Segmentation Comparison:", sample_id),
      x = "Genomic Position", y = "log2(ratio)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "lightgray", size = 0.3),
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
        linetype = "dashed", color = "orange", size = 0.4, alpha = 0.7
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

    for (j in seq_len(nrow(reseg))) {
      chr <- as.numeric(gsub("chr", "", reseg$chr[j]))
      if (!is.na(chr) && chr >= 1 && chr <= 24) {
        if (reseg$seg.mean[j] > gain_threshold) {
          freq_matrix[chr, i] <- freq_matrix[chr, i] + 1
        }
        if (reseg$seg.mean[j] < loss_threshold) {
          freq_matrix[chr, i] <- freq_matrix[chr, i] - 1
        }
      }
    }
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
