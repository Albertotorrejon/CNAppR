# TCGA-LIHC analysis — figures analogous to CNApp paper (Franch-Exposito et al., eLife 2020)
#
# Input:  inst/models/LIHC_354_cnvsegments_input_scores_nopurity.txt
# Output: figures_LIHC/ (PDF + PNG)
#
# Run from the package root:
#   devtools::load_all(".")
#   source("inst/scripts/lihc_paper_figures.R")

suppressPackageStartupMessages({
  library(CNAppR)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# ── Configuration ──────────────────────────────────────────────────────────────

DATA_FILE  <- system.file("models",
                           "LIHC_354_cnvsegments_input_scores_nopurity.txt",
                           package = "CNAppR")
if (!nzchar(DATA_FILE) || !file.exists(DATA_FILE))
  DATA_FILE <- "inst/models/LIHC_354_cnvsegments_input_scores_nopurity.txt"

OUT_DIR    <- "figures_LIHC"
GAIN_THR   <-  0.2    # log2 threshold to call a gain
LOSS_THR   <- -0.2    # log2 threshold to call a loss
BIN_SIZE   <-  1e6    # 1 Mb bins for CNA frequency profile
GENOME     <- "hg19"

dir.create(OUT_DIR, showWarnings = FALSE)

# ── 1. Load data ───────────────────────────────────────────────────────────────

raw <- read.table(DATA_FILE, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
raw$chr <- as.character(raw$chr)
raw$chr <- sub("^chr", "", raw$chr)    # strip chr prefix if present

sample_ids <- unique(raw$ID)
cat(sprintf("Loaded %d samples, %d segments\n", length(sample_ids), nrow(raw)))

# ── 2. Re-segmentation ─────────────────────────────────────────────────────────

cat("Re-segmenting samples...\n")
reseg_raw <- lapply(sample_ids, function(id) {
  tryCatch(
    resegment_sample(raw, sample_id = id,
                     dev_btw_segs  = 0.16,
                     dev_tozero    = 0.16,
                     min_length    = 10000,
                     chrom_percent = 0.9,
                     arm_percent   = 0.5),
    error = function(e) NULL
  )
})
names(reseg_raw) <- sample_ids

failed     <- names(which(sapply(reseg_raw, is.null)))
reseg_list <- Filter(Negate(is.null), reseg_raw)
if (length(failed) > 0)
  cat(sprintf("Skipped %d samples with no usable segments: %s\n",
              length(failed), paste(failed, collapse = ", ")))

# ── 3. CNA scores ──────────────────────────────────────────────────────────────

scores <- calculate_cna_scores(reseg_list)
scores$ID <- rownames(scores)
cat(sprintf("Scores computed for %d samples\n", nrow(scores)))
print(summary(scores[, c("FCS", "BCS", "GCS")]))

# ── Figure 1: Score distributions ─────────────────────────────────────────────
# Equivalent to Figure 2B / 3B in the paper:
# violin + boxplot for FCS, BCS, GCS

scores_long <- tidyr::pivot_longer(scores,
                                    cols      = c("FCS", "BCS", "GCS"),
                                    names_to  = "Score",
                                    values_to = "Value")
scores_long$Score <- factor(scores_long$Score, levels = c("FCS", "BCS", "GCS"))

pal_scores <- c(FCS = "#c0392b", BCS = "#2980b9", GCS = "#27ae60")

fig1 <- ggplot(scores_long, aes(x = Score, y = Value, fill = Score)) +
  geom_violin(alpha = 0.65, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "#222",
               outlier.size = 0.6, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = pal_scores) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title    = "TCGA-LIHC — CNA Score Distributions",
    subtitle = sprintf("n = %d hepatocellular carcinoma samples (no purity correction)",
                       length(sample_ids)),
    x        = NULL,
    y        = "Score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    axis.line        = element_line(color = "#444"),
    axis.ticks       = element_line(color = "#444")
  )

# ── Figure 2: Score correlations ───────────────────────────────────────────────
# Equivalent to Figure 2C in the paper (Spearman r between all score pairs)

cor_label <- function(x, y)
  sprintf("Spearman r = %.2f", cor(x, y, method = "spearman", use = "complete.obs"))

make_cor_panel <- function(df, xvar, yvar, color) {
  lbl <- cor_label(df[[xvar]], df[[yvar]])
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.45, size = 1.0, color = "#444") +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.08, vjust = -0.6,
             label = lbl, size = 3.2, color = "#444") +
    labs(x = xvar, y = yvar) +
    theme_classic(base_size = 11) +
    theme(axis.line = element_line(color = "#444"))
}

fig2 <- (make_cor_panel(scores, "BCS", "FCS", pal_scores["FCS"]) |
         make_cor_panel(scores, "BCS", "GCS", pal_scores["GCS"]) |
         make_cor_panel(scores, "FCS", "GCS", pal_scores["BCS"])) +
  plot_annotation(
    title = "TCGA-LIHC — CNA Score Correlations",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

# ── Figure 3: Genome-wide CNA frequency profile ────────────────────────────────
# Equivalent to Figure 3A in the paper: % of samples with gain or loss per 1 Mb bin

CHR_LENGTHS_HG19 <- c(
  "1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276,
  "5"=180915260, "6"=171115067, "7"=159138663, "8"=146364022,
  "9"=141213431, "10"=135534747,"11"=135006516,"12"=133851895,
  "13"=115169878,"14"=107349540,"15"=102531392,"16"=90354753,
  "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520,
  "21"=48129895, "22"=51304566
)
AUTOSOMES <- as.character(1:22)

# Build 1 Mb bin table
bins <- do.call(rbind, lapply(AUTOSOMES, function(ch) {
  len    <- CHR_LENGTHS_HG19[ch]
  starts <- seq(1, len, by = BIN_SIZE)
  data.frame(chr   = ch,
             b_idx = seq_along(starts),
             start = starts,
             end   = pmin(starts + BIN_SIZE - 1, len),
             stringsAsFactors = FALSE)
}))
bins$global_idx <- seq_len(nrow(bins))

# Cumulative chromosome offsets (for x-axis)
chr_offsets      <- c(0, cumsum(CHR_LENGTHS_HG19[AUTOSOMES[-22]]))
names(chr_offsets) <- AUTOSOMES
bins$genome_pos  <- chr_offsets[bins$chr] + (bins$start + bins$end) / 2

# Stack all re-segmented data
all_segs <- do.call(rbind, lapply(names(reseg_list), function(id) {
  s      <- reseg_list[[id]]
  s$ID   <- id
  s$chr  <- as.character(s$chr)
  s[s$chr %in% AUTOSOMES, ]
}))

# Count gains and losses per bin
n_samp    <- length(sample_ids)
gain_cnt  <- integer(nrow(bins))
loss_cnt  <- integer(nrow(bins))

bins_by_chr <- split(seq_len(nrow(bins)), bins$chr)

cat("Computing CNA frequency profile...\n")
for (ch in AUTOSOMES) {
  ch_segs    <- all_segs[all_segs$chr == ch, ]
  if (nrow(ch_segs) == 0) next
  ch_bin_idx <- bins_by_chr[[ch]]
  ch_bins    <- bins[ch_bin_idx, ]

  for (k in seq_len(nrow(ch_segs))) {
    s_start <- ch_segs$loc.start[k]
    s_end   <- ch_segs$loc.end[k]
    s_val   <- ch_segs$seg.mean[k]
    if (is.na(s_val)) next

    ov <- which(ch_bins$start <= s_end & ch_bins$end >= s_start)
    if (length(ov) == 0) next
    gi <- ch_bin_idx[ov]

    if      (s_val >  GAIN_THR) gain_cnt[gi] <- gain_cnt[gi] + 1L
    else if (s_val <  LOSS_THR) loss_cnt[gi] <- loss_cnt[gi] + 1L
  }
}

freq_df <- data.frame(
  genome_pos = bins$genome_pos,
  chr        = bins$chr,
  gain_pct   =  gain_cnt / n_samp * 100,
  loss_pct   = -loss_cnt / n_samp * 100
)

chr_mids   <- chr_offsets + CHR_LENGTHS_HG19[AUTOSOMES] / 2
chr_breaks <- data.frame(
  pos   = as.numeric(chr_mids),
  label = AUTOSOMES
)

# Alternating background bands per chromosome
chr_bands <- do.call(rbind, lapply(seq_along(AUTOSOMES), function(i) {
  ch <- AUTOSOMES[i]
  data.frame(
    xmin  = chr_offsets[ch],
    xmax  = chr_offsets[ch] + CHR_LENGTHS_HG19[ch],
    fill  = ifelse(i %% 2 == 0, "even", "odd")
  )
}))

fig3 <- ggplot() +
  geom_rect(data = chr_bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.25, inherit.aes = FALSE) +
  scale_fill_manual(values = c(odd = "#f5f5f5", even = "#ffffff"), guide = "none") +
  geom_area(data = freq_df, aes(x = genome_pos, y = gain_pct),
            fill = "#c0392b", alpha = 0.8) +
  geom_area(data = freq_df, aes(x = genome_pos, y = loss_pct),
            fill = "#2980b9", alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "#333") +
  scale_x_continuous(breaks = chr_breaks$pos, labels = chr_breaks$label,
                     expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) paste0(abs(x), "%"),
                     limits = c(-100, 100),
                     breaks = c(-75, -50, -25, 0, 25, 50, 75)) +
  labs(
    title    = "TCGA-LIHC — Genome-wide CNA Frequency Profile",
    subtitle = sprintf("n = %d samples  |  gains (red)  losses (blue)  |  threshold = ±%.1f log2R",
                       n_samp, GAIN_THR),
    x        = "Chromosome",
    y        = "% samples"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 9),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(size = 7.5),
    axis.line          = element_line(color = "#444")
  )

# ── Figure 4: CNA burden — high vs low GCS ─────────────────────────────────────
# Equivalent to the clinical stratification in Figures 3–5 of the paper

med_gcs         <- median(scores$GCS, na.rm = TRUE)
scores$GCS_group <- factor(
  ifelse(scores$GCS >= med_gcs, "High CNA burden", "Low CNA burden"),
  levels = c("Low CNA burden", "High CNA burden")
)

pval <- wilcox.test(GCS ~ GCS_group, data = scores)$p.value
pval_label <- ifelse(pval < 0.001, "p < 0.001",
              ifelse(pval < 0.01,  "p < 0.01",
              sprintf("p = %.3f", pval)))

fig4_gcs <- ggplot(scores, aes(x = GCS_group, y = GCS, fill = GCS_group)) +
  geom_violin(alpha = 0.65, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "#222",
               outlier.size = 0.6, outlier.alpha = 0.4, linewidth = 0.4) +
  annotate("text", x = 1.5, y = max(scores$GCS, na.rm = TRUE) * 1.02,
           label = pval_label, size = 3.5, fontface = "italic") +
  scale_fill_manual(values = c("Low CNA burden"  = "#2980b9",
                               "High CNA burden" = "#c0392b")) +
  labs(
    title    = "TCGA-LIHC — GCS stratification (median split)",
    subtitle = sprintf("Median GCS = %.1f  |  Wilcoxon %s", med_gcs, pval_label),
    x        = NULL,
    y        = "Global CNA Score (GCS)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(color = "grey40", size = 10),
    axis.line       = element_line(color = "#444")
  )

# Same for FCS and BCS side by side
scores$FCS_group <- factor(
  ifelse(scores$FCS >= median(scores$FCS, na.rm = TRUE),
         "High FCS", "Low FCS"),
  levels = c("Low FCS", "High FCS")
)
scores$BCS_group <- factor(
  ifelse(scores$BCS >= median(scores$BCS, na.rm = TRUE),
         "High BCS", "Low BCS"),
  levels = c("Low BCS", "High BCS")
)

fig4_fcs <- ggplot(scores, aes(x = FCS_group, y = FCS, fill = FCS_group)) +
  geom_violin(alpha = 0.65, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "#222",
               outlier.size = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = c("Low FCS" = "#2980b9", "High FCS" = "#c0392b")) +
  labs(x = NULL, y = "Focal CNA Score (FCS)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.line = element_line(color = "#444"))

fig4_bcs <- ggplot(scores, aes(x = BCS_group, y = BCS, fill = BCS_group)) +
  geom_violin(alpha = 0.65, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "#222",
               outlier.size = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = c("Low BCS" = "#2980b9", "High BCS" = "#c0392b")) +
  labs(x = NULL, y = "Broad CNA Score (BCS)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.line = element_line(color = "#444"))

fig4 <- (fig4_fcs | fig4_bcs | fig4_gcs) +
  plot_annotation(
    title = "TCGA-LIHC — Score Stratification (Median Split)",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

# ── Figure 5: Cumulative distribution of GCS ───────────────────────────────────
# ECDF — useful to show heterogeneity across tumours

fig5 <- ggplot(scores, aes(x = GCS)) +
  stat_ecdf(geom = "step", color = "#27ae60", linewidth = 0.8) +
  geom_vline(xintercept = med_gcs, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = med_gcs, y = 0.05, hjust = -0.1,
           label = sprintf("median = %.1f", med_gcs),
           size = 3.2, color = "grey40") +
  labs(
    title    = "TCGA-LIHC — Cumulative Distribution of GCS",
    subtitle = sprintf("n = %d samples", n_samp),
    x        = "Global CNA Score (GCS)",
    y        = "Cumulative proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    axis.line     = element_line(color = "#444")
  )

# ── Save all figures ───────────────────────────────────────────────────────────

save_fig <- function(fig, name, width, height) {
  path_pdf <- file.path(OUT_DIR, paste0(name, ".pdf"))
  path_png <- file.path(OUT_DIR, paste0(name, ".png"))
  ggsave(path_pdf, fig, width = width, height = height)
  ggsave(path_png, fig, width = width, height = height, dpi = 300)
  cat("Saved:", name, "\n")
}

save_fig(fig1, "fig1_score_distributions",    width = 6,  height = 5)
save_fig(fig2, "fig2_score_correlations",     width = 11, height = 4)
save_fig(fig3, "fig3_cna_frequency_profile",  width = 13, height = 4.5)
save_fig(fig4, "fig4_score_stratification",   width = 11, height = 5)
save_fig(fig5, "fig5_gcs_ecdf",               width = 6,  height = 5)

cat(sprintf("\nAll figures written to %s/\n", OUT_DIR))




