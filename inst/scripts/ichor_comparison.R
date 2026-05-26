# Comparison: ichorCNA segmentation vs CNAppR resegmentation
# Input:  inst/models/Datos laia/barcode{66,68}.aln.sorted.seg.txt
# Output: figures_ichor/ (PDF plots + comparison statistics)

library(CNAppR)
library(ggplot2)

dir.create("figures_ichor", showWarnings = FALSE)

samples <- c("barcode66.aln.sorted", "barcode68.aln.sorted")

# ── Chromosome order and cumulative offsets for genome-wide x-axis ────────────
chr_lengths <- c(
  `1`=248956422, `2`=242193529, `3`=198295559, `4`=190214555,
  `5`=181538259, `6`=170805979, `7`=159345973, `8`=145138636,
  `9`=138394717, `10`=133797422, `11`=135086622, `12`=133275309,
  `13`=114364328, `14`=107043718, `15`=101991189, `16`=90338345,
  `17`=83257441, `18`=80373285, `19`=58617616, `20`=64444167,
  `21`=46709983, `22`=50818468
)
chr_offsets <- c(0, cumsum(chr_lengths[-length(chr_lengths)]))
names(chr_offsets) <- names(chr_lengths)
chr_mids <- chr_offsets + chr_lengths / 2

# ── Helper: call → colour ─────────────────────────────────────────────────────
call_colour <- function(call) {
  dplyr::case_when(
    call %in% c("AMP", "GAIN", "HLAMP")        ~ "#e74c3c",
    call %in% c("HETD", "HOMD", "DEL", "LOSS") ~ "#27ae60",
    TRUE                                         ~ "#2980b9"
  )
}

# ── Helper: CNAppR classified → colour ───────────────────────────────────────
cnapp_colour <- function(classified) {
  dplyr::case_when(
    grepl("Gain|gain|AMP|amp", classified) ~ "#e74c3c",
    grepl("Loss|loss|DEL|del", classified) ~ "#27ae60",
    TRUE                                    ~ "#2980b9"
  )
}

# ── Helper: genome-wide segment data.frame ────────────────────────────────────
to_genome_coords <- function(seg_df, chr_col, start_col, end_col) {
  seg_df$chr_num <- suppressWarnings(as.integer(gsub("chr", "", seg_df[[chr_col]])))
  seg_df <- seg_df[!is.na(seg_df$chr_num) & seg_df$chr_num %in% 1:22, ]
  offset <- chr_offsets[as.character(seg_df$chr_num)]
  seg_df$gstart <- seg_df[[start_col]] + offset
  seg_df$gend   <- seg_df[[end_col]]   + offset
  seg_df$gmid   <- (seg_df$gstart + seg_df$gend) / 2
  seg_df
}

# ── Stats: segment-level concordance ─────────────────────────────────────────
compare_segs <- function(ichor, cnapp) {
  # For each ichorCNA segment, find CNAppR segments that overlap it
  # and compute weighted correlation of log2 ratio
  ichor$length <- ichor$end - ichor$start
  cnapp$length <- cnapp$loc.end - cnapp$loc.start

  # Build a common grid: split genome into 1 Mb bins, assign log2 from each tool
  bins <- seq(1, 3e9, by = 1e6)
  get_seg_value <- function(pos, chr_val, seg_df, chr_col, start_col, end_col, val_col) {
    hit <- seg_df[[chr_col]] == chr_val &
           seg_df[[start_col]] <= pos &
           seg_df[[end_col]]   >= pos
    if (any(hit)) seg_df[[val_col]][which(hit)[1]] else NA_real_
  }

  results <- lapply(1:22, function(ch) {
    bin_starts <- seq(1, chr_lengths[as.character(ch)], by = 1e6)
    ichor_vals <- sapply(bin_starts, get_seg_value,
                         chr_val   = ch,
                         seg_df    = ichor,
                         chr_col   = "chrom",
                         start_col = "start",
                         end_col   = "end",
                         val_col   = "seg.median.logR")
    cnapp_vals <- sapply(bin_starts, get_seg_value,
                         chr_val   = ch,
                         seg_df    = cnapp,
                         chr_col   = "chr",
                         start_col = "loc.start",
                         end_col   = "loc.end",
                         val_col   = "seg.mean")
    data.frame(ichor = ichor_vals, cnapp = cnapp_vals)
  })
  grid <- do.call(rbind, results)
  grid <- grid[!is.na(grid$ichor) & !is.na(grid$cnapp), ]

  list(
    pearson  = cor(grid$ichor, grid$cnapp, method = "pearson"),
    spearman = cor(grid$ichor, grid$cnapp, method = "spearman"),
    n_bins   = nrow(grid),
    n_ichor  = nrow(ichor),
    n_cnapp  = nrow(cnapp),
    grid     = grid
  )
}

# ═════════════════════════════════════════════════════════════════════════════
all_stats <- list()

for (sid in samples) {

  cat("\n══ Processing:", sid, "══\n")

  seg_path <- file.path(
    system.file("models", package = "CNAppR"),
    "Datos laia",
    paste0(sid, ".seg.txt")
  )

  ichor_raw <- read.table(seg_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)

  # ── 1. ichorCNA genome-wide plot ──────────────────────────────────────────
  ichor_gw <- to_genome_coords(ichor_raw, "chrom", "start", "end")
  ichor_gw$colour <- call_colour(ichor_gw$Corrected_Call)

  ichor_label <- sprintf(
    "ichorCNA  |  Tumor Fraction: %.4f  |  Ploidy: %.2f",
    ichor_raw$seg.median.logR[1] * 0,   # placeholder — read from header if available
    2
  )
  # Extract tumor fraction from first row ID string is not available in seg.txt
  # Use fixed values from the plots shared by Laia
  tf_map   <- c("barcode66.aln.sorted" = 0.4801, "barcode68.aln.sorted" = 0.3355)
  ply_map  <- c("barcode66.aln.sorted" = 1.96,   "barcode68.aln.sorted" = 2.01)
  tf       <- tf_map[sid]
  ploidy   <- ply_map[sid]

  p_ichor <- ggplot(ichor_gw) +
    geom_segment(aes(x = gstart, xend = gend,
                     y = seg.median.logR, yend = seg.median.logR,
                     colour = colour), linewidth = 1.2) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = chr_mids,
      labels = names(chr_mids),
      expand = c(0.01, 0)
    ) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(
      title    = sid,
      subtitle = sprintf("ichorCNA  |  Tumor Fraction: %.4f  |  Ploidy: %.2f",
                         tf, ploidy),
      x        = "Chromosome",
      y        = "Copy Number (log2 ratio)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(size = 8),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
    )

  # ── 2. Convert ichorCNA → CNAppR format and resegment ────────────────────
  cnapp_input <- data.frame(
    ID        = sid,
    chr       = ichor_raw$chrom,
    loc.start = ichor_raw$start,
    loc.end   = ichor_raw$end,
    seg.mean  = ichor_raw$seg.median.logR,
    stringsAsFactors = FALSE
  )
  cnapp_input$chr <- suppressWarnings(as.integer(gsub("chr", "", cnapp_input$chr)))
  cnapp_input <- cnapp_input[!is.na(cnapp_input$chr) & cnapp_input$chr %in% 1:22, ]

  reseg <- tryCatch(
    resegment_sample(cnapp_input, sample_id = sid),
    error = function(e) { cat("resegment error:", e$message, "\n"); NULL }
  )

  if (is.null(reseg)) next

  # ── 3. CNAppR genome-wide plot ────────────────────────────────────────────
  cnapp_gw <- to_genome_coords(reseg, "chr", "loc.start", "loc.end")
  cnapp_gw$colour <- cnapp_colour(cnapp_gw$classified)

  p_cnapp <- ggplot(cnapp_gw) +
    geom_segment(aes(x = gstart, xend = gend,
                     y = seg.mean, yend = seg.mean,
                     colour = colour), linewidth = 1.2) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = chr_mids,
      labels = names(chr_mids),
      expand = c(0.01, 0)
    ) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(
      title    = sid,
      subtitle = "CNAppR resegmentation",
      x        = "Chromosome",
      y        = "Copy Number (log2 ratio)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(size = 8),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
    )

  # ── 4. Save combined plot ─────────────────────────────────────────────────
  out_pdf <- file.path("figures_ichor", paste0(sid, "_comparison.pdf"))
  pdf(out_pdf, width = 14, height = 8)
  gridExtra::grid.arrange(p_ichor, p_cnapp, ncol = 1)
  dev.off()
  cat("Plot saved:", out_pdf, "\n")

  # ── 5. Concordance statistics ─────────────────────────────────────────────
  stats <- compare_segs(
    ichor = cnapp_input,   # use raw ichor values
    cnapp = reseg
  )

  cat(sprintf(
    "  Pearson r  = %.4f\n  Spearman r = %.4f\n  Bins       = %d\n  ichor segs = %d  |  CNAppR segs = %d\n",
    stats$pearson, stats$spearman, stats$n_bins, stats$n_ichor, stats$n_cnapp
  ))

  all_stats[[sid]] <- data.frame(
    sample       = sid,
    pearson_r    = round(stats$pearson,  4),
    spearman_r   = round(stats$spearman, 4),
    n_bins       = stats$n_bins,
    n_segs_ichor = stats$n_ichor,
    n_segs_cnapp = stats$n_cnapp
  )

  # ── 6. Scatter: ichor vs CNAppR log2 ratio per 1 Mb bin ──────────────────
  p_scatter <- ggplot(stats$grid, aes(x = ichor, y = cnapp)) +
    geom_point(alpha = 0.3, size = 0.8, colour = "#2c3e50") +
    geom_abline(slope = 1, intercept = 0, colour = "#e74c3c",
                linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, colour = "#2980b9",
                linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = sprintf("Pearson r = %.3f\nSpearman r = %.3f",
                             stats$pearson, stats$spearman),
             size = 3.5) +
    labs(
      title = sid,
      subtitle = "ichorCNA vs CNAppR — log2 ratio per 1 Mb bin",
      x = "ichorCNA seg.median.logR",
      y = "CNAppR seg.mean"
    ) +
    theme_classic(base_size = 11)

  out_scatter <- file.path("figures_ichor", paste0(sid, "_scatter.pdf"))
  ggsave(out_scatter, p_scatter, width = 6, height = 5)
  cat("Scatter saved:", out_scatter, "\n")
}

# ── Summary table ─────────────────────────────────────────────────────────────
summary_df <- do.call(rbind, all_stats)
cat("\n── Summary ──────────────────────────────────────────────\n")
print(summary_df, row.names = FALSE)
write.table(summary_df,
            file.path("figures_ichor", "comparison_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nStats saved: figures_ichor/comparison_stats.tsv\n")
