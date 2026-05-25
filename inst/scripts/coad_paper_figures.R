# TCGA-COAD analysis — improved figures from CNApp paper (Franch-Exposito et al., eLife 2020)
#
# Input:  inst/models/demo_data_random160COAD_cms_SURVIVAL_14.11.2018.txt
# Output: figures_COAD/ (PDF + PNG)
#
# Run from the package root:
#   devtools::load_all(".")
#   source("inst/scripts/coad_paper_figures.R")

suppressPackageStartupMessages({
  library(CNAppR)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# ── Configuration ──────────────────────────────────────────────────────────────

DATA_FILE <- system.file("models",
                          "demo_data_random160COAD_cms_SURVIVAL_14.11.2018.txt",
                          package = "CNAppR")
if (!nzchar(DATA_FILE) || !file.exists(DATA_FILE))
  DATA_FILE <- "inst/models/demo_data_random160COAD_cms_SURVIVAL_14.11.2018.txt"

OUT_DIR  <- "figures_COAD"
GAIN_THR <-  0.2
LOSS_THR <- -0.2
BIN_SIZE <-  1e6
GENOME   <- "hg19"

dir.create(OUT_DIR, showWarnings = FALSE)

PAL <- c(FCS = "#c0392b", BCS = "#2980b9", GCS = "#27ae60")
PAL_HL <- c(Low = "#2980b9", High = "#c0392b")

# ── 1. Load data ───────────────────────────────────────────────────────────────

raw <- read.table(DATA_FILE, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
raw$chr <- as.character(raw$chr)
raw$chr <- sub("^chr", "", raw$chr)

sample_ids <- unique(raw$ID)
cat(sprintf("Loaded %d samples, %d segments\n", length(sample_ids), nrow(raw)))

# Extract clinical metadata (one row per sample)
clinical_cols <- c("gender", "tumorLocation", "msi", "cimp",
                   "kras_mut", "braf_mut", "surv_status", "surv_time")
clin <- raw[!duplicated(raw$ID), c("ID", clinical_cols)]
rownames(clin) <- clin$ID

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
  cat(sprintf("Skipped %d samples: %s\n", length(failed),
              paste(failed, collapse = ", ")))

# ── 3. CNA scores ──────────────────────────────────────────────────────────────

scores <- calculate_cna_scores(reseg_list)
scores$ID <- rownames(scores)
scores <- merge(scores, clin, by = "ID", all.x = TRUE)
cat(sprintf("Scores computed for %d samples\n", nrow(scores)))
print(summary(scores[, c("FCS", "BCS", "GCS")]))

# Median-split groups
med_gcs <- median(scores$GCS, na.rm = TRUE)
med_fcs <- median(scores$FCS, na.rm = TRUE)
med_bcs <- median(scores$BCS, na.rm = TRUE)
scores$GCS_group <- factor(ifelse(scores$GCS >= med_gcs, "High", "Low"),
                            levels = c("Low", "High"))
scores$FCS_group <- factor(ifelse(scores$FCS >= med_fcs, "High", "Low"),
                            levels = c("Low", "High"))
scores$BCS_group <- factor(ifelse(scores$BCS >= med_bcs, "High", "Low"),
                            levels = c("Low", "High"))

# ── Helper: violin + boxplot panel ─────────────────────────────────────────────

violin_panel <- function(df, xvar, yvar, fill_var, pal, xlab = NULL, ylab = yvar,
                         title = NULL, subtitle = NULL) {
  pval <- tryCatch(
    wilcox.test(df[[yvar]] ~ df[[xvar]])$p.value,
    error = function(e) NA
  )
  pval_lbl <- if (is.na(pval)) ""
              else if (pval < 0.001) "p < 0.001"
              else if (pval < 0.01)  "p < 0.01"
              else if (pval < 0.05)  sprintf("p = %.3f *", pval)
              else                   sprintf("p = %.3f", pval)

  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[fill_var]])) +
    geom_violin(alpha = 0.65, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.12, fill = "white", color = "#222",
                 outlier.size = 0.7, outlier.alpha = 0.4, linewidth = 0.4) +
    annotate("text",
             x     = length(unique(df[[xvar]])) / 2 + 0.5,
             y     = max(df[[yvar]], na.rm = TRUE) * 1.05,
             label = pval_lbl, size = 3.2, color = "#444", fontface = "italic") +
    scale_fill_manual(values = pal) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(color = "grey40", size = 9),
          axis.line     = element_line(color = "#444"))
}

# ── Figure 1: Score distributions ─────────────────────────────────────────────

scores_long <- pivot_longer(scores, cols = c("FCS", "BCS", "GCS"),
                             names_to = "Score", values_to = "Value")
scores_long$Score <- factor(scores_long$Score, levels = c("FCS", "BCS", "GCS"))

fig1 <- ggplot(scores_long, aes(x = Score, y = Value, fill = Score)) +
  geom_violin(alpha = 0.65, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "#222",
               outlier.size = 0.6, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = PAL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.07))) +
  labs(
    title    = "TCGA-COAD — CNA Score Distributions",
    subtitle = sprintf("n = %d colorectal cancer samples", nrow(scores)),
    x        = NULL, y = "Score"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(color = "grey40", size = 10),
        axis.line        = element_line(color = "#444"))

# ── Figure 2: Score correlations ───────────────────────────────────────────────

cor_label <- function(x, y)
  sprintf("Spearman r = %.2f", cor(x, y, method = "spearman", use = "complete.obs"))

cor_panel <- function(df, xv, yv, col) {
  lbl <- cor_label(df[[xv]], df[[yv]])
  ggplot(df, aes(x = .data[[xv]], y = .data[[yv]])) +
    geom_point(alpha = 0.4, size = 1.0, color = "#444") +
    geom_smooth(method = "lm", se = FALSE, color = col, linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6,
             label = lbl, size = 3.2, color = "#444") +
    labs(x = xv, y = yv) +
    theme_classic(base_size = 11) +
    theme(axis.line = element_line(color = "#444"))
}

fig2 <- (cor_panel(scores, "BCS", "FCS", PAL["FCS"]) |
         cor_panel(scores, "BCS", "GCS", PAL["GCS"]) |
         cor_panel(scores, "FCS", "GCS", PAL["BCS"])) +
  plot_annotation(
    title = "TCGA-COAD — CNA Score Correlations",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

# ── Figure 3: Genome-wide CNA frequency profile ────────────────────────────────

CHR_LENGTHS_HG19 <- c(
  "1"=249250621,"2"=243199373,"3"=198022430,"4"=191154276,
  "5"=180915260,"6"=171115067,"7"=159138663,"8"=146364022,
  "9"=141213431,"10"=135534747,"11"=135006516,"12"=133851895,
  "13"=115169878,"14"=107349540,"15"=102531392,"16"=90354753,
  "17"=81195210,"18"=78077248,"19"=59128983,"20"=63025520,
  "21"=48129895,"22"=51304566
)
AUTOSOMES <- as.character(1:22)

bins <- do.call(rbind, lapply(AUTOSOMES, function(ch) {
  len    <- CHR_LENGTHS_HG19[ch]
  starts <- seq(1, len, by = BIN_SIZE)
  data.frame(chr = ch, start = starts, end = pmin(starts + BIN_SIZE - 1, len),
             stringsAsFactors = FALSE)
}))

chr_offsets     <- c(0, cumsum(CHR_LENGTHS_HG19[AUTOSOMES[-22]]))
names(chr_offsets) <- AUTOSOMES
bins$genome_pos <- chr_offsets[bins$chr] + (bins$start + bins$end) / 2

all_segs <- do.call(rbind, lapply(names(reseg_list), function(id) {
  s <- reseg_list[[id]]; s$ID <- id; s$chr <- as.character(s$chr)
  s[s$chr %in% AUTOSOMES, ]
}))

n_samp   <- length(reseg_list)
gain_cnt <- integer(nrow(bins))
loss_cnt <- integer(nrow(bins))
bins_by_chr <- split(seq_len(nrow(bins)), bins$chr)

cat("Computing CNA frequency profile...\n")
for (ch in AUTOSOMES) {
  ch_segs    <- all_segs[all_segs$chr == ch, ]
  if (nrow(ch_segs) == 0) next
  ch_bin_idx <- bins_by_chr[[ch]]
  ch_bins    <- bins[ch_bin_idx, ]
  for (k in seq_len(nrow(ch_segs))) {
    sv <- ch_segs$seg.mean[k]; if (is.na(sv)) next
    ov <- which(ch_bins$start <= ch_segs$loc.end[k] &
                ch_bins$end   >= ch_segs$loc.start[k])
    if (length(ov) == 0) next
    gi <- ch_bin_idx[ov]
    if      (sv >  GAIN_THR) gain_cnt[gi] <- gain_cnt[gi] + 1L
    else if (sv <  LOSS_THR) loss_cnt[gi] <- loss_cnt[gi] + 1L
  }
}

freq_df <- data.frame(
  genome_pos = bins$genome_pos,
  chr        = bins$chr,
  gain_pct   =  gain_cnt / n_samp * 100,
  loss_pct   = -loss_cnt / n_samp * 100
)

chr_mids  <- chr_offsets + CHR_LENGTHS_HG19[AUTOSOMES] / 2
chr_bands <- do.call(rbind, lapply(seq_along(AUTOSOMES), function(i) {
  ch <- AUTOSOMES[i]
  data.frame(xmin = chr_offsets[ch],
             xmax = chr_offsets[ch] + CHR_LENGTHS_HG19[ch],
             fill = ifelse(i %% 2 == 0, "even", "odd"))
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
  scale_x_continuous(breaks = as.numeric(chr_mids),
                     labels = AUTOSOMES, expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) paste0(abs(x), "%"),
                     limits = c(-100, 100),
                     breaks = c(-75, -50, -25, 0, 25, 50, 75)) +
  labs(
    title    = "TCGA-COAD — Genome-wide CNA Frequency Profile",
    subtitle = sprintf("n = %d samples  |  gains (red)  losses (blue)  |  threshold = ±%.1f log2R",
                       n_samp, GAIN_THR),
    x = "Chromosome", y = "% samples"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "grey40", size = 9),
        panel.grid.major.x = element_blank(),
        axis.text.x        = element_text(size = 7.5),
        axis.line          = element_line(color = "#444"))

# ── Figure 4: Clinical associations — MSI, CIMP, KRAS, BRAF ───────────────────
# Equivalent to Figures 4-5 in the eLife paper

# MSI
msi_df <- scores[!is.na(scores$msi), ]
msi_df$msi <- factor(msi_df$msi, levels = c("mss", "msi"))
pal_msi <- c(mss = "#aaaaaa", msi = "#e67e22")

p_msi_fcs <- violin_panel(msi_df, "msi", "FCS", "msi", pal_msi, ylab = "FCS")
p_msi_bcs <- violin_panel(msi_df, "msi", "BCS", "msi", pal_msi, ylab = "BCS")
p_msi_gcs <- violin_panel(msi_df, "msi", "GCS", "msi", pal_msi, ylab = "GCS")

fig4 <- (p_msi_fcs | p_msi_bcs | p_msi_gcs) +
  plot_annotation(
    title    = "TCGA-COAD — CNA Scores by MSI Status",
    subtitle = sprintf("MSS n=%d  |  MSI n=%d",
                       sum(msi_df$msi == "mss", na.rm = TRUE),
                       sum(msi_df$msi == "msi",  na.rm = TRUE)),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(color = "grey40", size = 10))
  )

# CIMP
cimp_df <- scores[!is.na(scores$cimp), ]
cimp_df$cimp <- factor(cimp_df$cimp, levels = c("CIMP.Low", "CIMP.High"))
pal_cimp <- c(CIMP.Low = "#aaaaaa", CIMP.High = "#8e44ad")

p_cimp_fcs <- violin_panel(cimp_df, "cimp", "FCS", "cimp", pal_cimp, ylab = "FCS")
p_cimp_bcs <- violin_panel(cimp_df, "cimp", "BCS", "cimp", pal_cimp, ylab = "BCS")
p_cimp_gcs <- violin_panel(cimp_df, "cimp", "GCS", "cimp", pal_cimp, ylab = "GCS")

fig5 <- (p_cimp_fcs | p_cimp_bcs | p_cimp_gcs) +
  plot_annotation(
    title    = "TCGA-COAD — CNA Scores by CIMP Status",
    subtitle = sprintf("CIMP-Low n=%d  |  CIMP-High n=%d",
                       sum(cimp_df$cimp == "CIMP.Low",  na.rm = TRUE),
                       sum(cimp_df$cimp == "CIMP.High", na.rm = TRUE)),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(color = "grey40", size = 10))
  )

# KRAS / BRAF mutations
mut_df <- scores[!is.na(scores$kras_mut) & !is.na(scores$braf_mut), ]
mut_df$kras_mut <- factor(mut_df$kras_mut, levels = c("wt", "mut"))
mut_df$braf_mut <- factor(mut_df$braf_mut, levels = c("wt", "mut"))
pal_mut <- c(wt = "#aaaaaa", mut = "#c0392b")

p_kras <- violin_panel(mut_df, "kras_mut", "GCS", "kras_mut", pal_mut,
                        xlab = "KRAS", ylab = "GCS")
p_braf <- violin_panel(mut_df, "braf_mut", "GCS", "braf_mut", pal_mut,
                        xlab = "BRAF", ylab = "GCS")

fig6 <- (p_kras | p_braf) +
  plot_annotation(
    title = "TCGA-COAD — GCS by KRAS / BRAF Mutation Status",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

# ── Figure 7: Tumor location ───────────────────────────────────────────────────

loc_df <- scores[!is.na(scores$tumorLocation), ]
pal_loc <- c(right = "#2ecc71", left = "#e74c3c", rectal = "#f39c12")

p_loc_fcs <- violin_panel(loc_df, "tumorLocation", "FCS", "tumorLocation",
                            pal_loc, ylab = "FCS")
p_loc_bcs <- violin_panel(loc_df, "tumorLocation", "BCS", "tumorLocation",
                            pal_loc, ylab = "BCS")
p_loc_gcs <- violin_panel(loc_df, "tumorLocation", "GCS", "tumorLocation",
                            pal_loc, ylab = "GCS")

fig7 <- (p_loc_fcs | p_loc_bcs | p_loc_gcs) +
  plot_annotation(
    title = "TCGA-COAD — CNA Scores by Tumor Location",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

# ── Figure 8: Survival analysis (Kaplan-Meier by GCS) ─────────────────────────

surv_df <- scores[!is.na(scores$surv_time) & !is.na(scores$surv_status), ]

if (nrow(surv_df) >= 20 &&
    requireNamespace("survival",    quietly = TRUE) &&
    requireNamespace("survminer",   quietly = TRUE)) {

  library(survival)
  library(survminer)

  fit <- survfit(Surv(surv_time, surv_status) ~ GCS_group, data = surv_df)

  fig8 <- ggsurvplot(
    fit, data = surv_df,
    pval          = TRUE,
    risk.table    = TRUE,
    risk.table.col= "strata",
    palette       = c("#2980b9", "#c0392b"),
    title         = "TCGA-COAD — Overall Survival by GCS (median split)",
    xlab          = "Days",
    ylab          = "Survival probability",
    legend.labs   = c("Low GCS", "High GCS"),
    ggtheme       = theme_classic(base_size = 12)
  )

  pdf(file.path(OUT_DIR, "fig8_survival_gcs.pdf"), width = 8, height = 7)
  print(fig8)
  dev.off()
  png(file.path(OUT_DIR, "fig8_survival_gcs.png"), width = 8, height = 7,
      units = "in", res = 300)
  print(fig8)
  dev.off()
  cat("Saved: fig8_survival_gcs\n")

} else {
  cat("Skipping survival figure (survminer not installed or too few samples with survival data)\n")
  cat("Install with: install.packages(c('survival','survminer'))\n")
}

# ── Save figures 1-7 ──────────────────────────────────────────────────────────

save_fig <- function(fig, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), fig, width = width, height = height)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), fig, width = width, height = height, dpi = 300)
  cat("Saved:", name, "\n")
}

save_fig(fig1, "fig1_score_distributions",   width = 6,  height = 5)
save_fig(fig2, "fig2_score_correlations",    width = 11, height = 4)
save_fig(fig3, "fig3_cna_frequency_profile", width = 13, height = 4.5)
save_fig(fig4, "fig4_cna_vs_msi",            width = 10, height = 5)
save_fig(fig5, "fig5_cna_vs_cimp",           width = 10, height = 5)
save_fig(fig6, "fig6_cna_vs_mutations",      width = 7,  height = 5)
save_fig(fig7, "fig7_cna_vs_location",       width = 10, height = 5)

# ── Export score table ─────────────────────────────────────────────────────────

score_export <- scores[, c("ID", "FCS", "BCS", "GCS",
                            "msi", "cimp", "kras_mut", "braf_mut",
                            "tumorLocation", "surv_status", "surv_time")]
write.table(score_export,
            file      = file.path(OUT_DIR, "COAD_scores_with_clinical.tsv"),
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)
cat("Saved: COAD_scores_with_clinical.tsv\n")

cat(sprintf("\nAll outputs written to %s/\n", OUT_DIR))
