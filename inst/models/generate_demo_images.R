# ============================================================
# Genera imágenes demo para docs/img/ usando datos reales TFM
# Funciones: plot_frequency_profile, plot_genome_wide_cna, plot_segmentation
# ============================================================
library(CNAppR)
library(ggplot2)

# ── Rutas ────────────────────────────────────────────────────────────
pkg_root  <- system.file(package = "CNAppR")
data_dir  <- file.path(pkg_root, "models/datos TFM")

aux_hg19 <- file.path(
  pkg_root, "..", "..", "..", "Github Sebas",
  "Githubsebas", "aux_files", "segmented_files_hg19"
)

arm_ref  <- file.path(aux_hg19, "autosomes_hg19_by_arms.txt")
wes_bed  <- file.path(pkg_root, "extdata", "standard_exome_hg19.bed")

lihc_path <- file.path(
  data_dir, "LIHC_354_cnvsegments_input_scores_nopurity.txt"
)
bc66_path <- file.path(data_dir, "barcode66.aln.sorted.seg.txt")
bc68_path <- file.path(data_dir, "barcode68.aln.sorted.seg.txt")

out_dir <- normalizePath(
  file.path(pkg_root, "..", "..", "docs", "img"),
  mustWork = FALSE
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Carga y prepara datos LIHC ────────────────────────────────────────
lihc_raw  <- read.table(lihc_path, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
lihc_list <- lapply(
  split(lihc_raw, lihc_raw$ID),
  function(df) df[, c("chr", "loc.start", "loc.end", "seg.mean")]
)
n_lihc <- length(lihc_list)
cat(sprintf("LIHC cargado: %d muestras\n", n_lihc))

# ── 1a. plot_frequency_profile — resolución por brazos ────────────────
cat("Calculando frecuencias (brazos)...\n")
freq_arms <- compute_region_frequencies(lihc_list, arm_ref,
                                        gain_thr = 0.23, loss_thr = -0.23)

p_arms <- plot_frequency_profile(
  freq_df    = freq_arms,
  title      = sprintf("Copy Number Frequency  |  TCGA-LIHC  (n = %d)",
                       n_lihc),
  gain_color = "#d73027",
  loss_color = "#4575b4"
)

ggsave(file.path(out_dir, "lihc_frequency_arms.png"),
       p_arms, width = 12, height = 4, dpi = 150, bg = "white")
cat("  -> lihc_frequency_arms.png\n")

# ── 1b. plot_frequency_profile — resolución WES targets ───────────────
cat("Calculando frecuencias (WES targets, puede tardar unos minutos)...\n")
freq_wes <- compute_region_frequencies(lihc_list, wes_bed,
                                       gain_thr = 0.23, loss_thr = -0.23)

p_wes <- plot_frequency_profile(
  freq_df    = freq_wes,
  title      = sprintf("Copy Number Frequency  |  TCGA-LIHC  WES targets  (n = %d)", n_lihc),
  gain_color = "#d73027",
  loss_color = "#4575b4"
)

ggsave(file.path(out_dir, "lihc_frequency_wes.png"),
       p_wes, width = 12, height = 4, dpi = 150, bg = "white")
cat("  -> lihc_frequency_wes.png\n")

# ── 2. plot_genome_wide_cna — Nanopore barcode66 + barcode68 ──────────
cat("Generando plot_genome_wide_cna...\n")

read_ichor_seg <- function(path) {
  seg <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  seg$chr       <- seg$chrom
  seg$loc.start <- seg$start
  seg$loc.end   <- seg$end
  seg
}

# Simula bins 1 Mb a partir de segmentos + ruido (substituye read_bam_qc)
sim_bins_from_segs <- function(seg_df, bin_size = 1e6, sd_noise = 0.15) {
  set.seed(42)
  rows <- lapply(seq_len(nrow(seg_df)), function(i) {
    s   <- seg_df[i, ]
    pos <- seq(s$loc.start,
               max(s$loc.start, s$loc.end - bin_size),
               by = bin_size)
    data.frame(
      chr        = s$chr,
      loc.start  = pos,
      loc.end    = pos + bin_size - 1L,
      log2_ratio = s$seg.median.logR + rnorm(length(pos), 0, sd_noise)
    )
  })
  do.call(rbind, rows)
}

bc66  <- read_ichor_seg(bc66_path)
bc68  <- read_ichor_seg(bc68_path)
bins66 <- sim_bins_from_segs(bc66)
bins68 <- sim_bins_from_segs(bc68)

p66 <- plot_genome_wide_cna(bins66, bc66,
                              sample_id      = "barcode66  |  Nanopore WGS",
                              tumor_fraction = 0.48, ploidy = 1.96)
p68 <- plot_genome_wide_cna(bins68, bc68,
                              sample_id = "barcode68  |  Nanopore WGS")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  panel_nano <- p66 / p68
} else {
  panel_nano <- gridExtra::arrangeGrob(p66, p68, ncol = 1)
}

ggsave(file.path(out_dir, "nanopore_profiles.png"),
       panel_nano, width = 12, height = 7, dpi = 150, bg = "white")
cat("  -> nanopore_profiles.png\n")

# ── 3. plot_segmentation — comparación antes/después ─────────────────
cat("Generando plot_segmentation...\n")

# Muestra con CNAs visibles
sample_id_demo <- names(lihc_list)[[42]]
orig_df <- lihc_raw[lihc_raw$ID == sample_id_demo,
                    c("ID", "chr", "loc.start", "loc.end", "seg.mean")]

# Resegmentado simulado: fusiona adyacentes con diferencia < 0.15
merge_adjacent_segs <- function(df, tol = 0.15) {
  if (nrow(df) <= 1L) return(df)
  out <- df[1L, , drop = FALSE]
  for (i in seq(2L, nrow(df))) {
    last <- nrow(out)
    if (df$chr[i] == out$chr[last] &&
        abs(df$seg.mean[i] - out$seg.mean[last]) < tol) {
      out$loc.end[last]  <- df$loc.end[i]
      out$seg.mean[last] <- (out$seg.mean[last] + df$seg.mean[i]) / 2
    } else {
      out <- rbind(out, df[i, ])
    }
  }
  rownames(out) <- NULL
  out
}

reseg_df <- merge_adjacent_segs(orig_df)
cat(sprintf("  %s: %d segs -> %d re-segmentados\n",
            sample_id_demo, nrow(orig_df), nrow(reseg_df)))

p_seg <- plot_segmentation(orig_df, reseg_df, sample_id = sample_id_demo)
ggsave(file.path(out_dir, "segmentation_comparison.png"),
       p_seg, width = 12, height = 4, dpi = 150, bg = "white")
cat("  -> segmentation_comparison.png\n")

cat("\nTodas las imágenes guardadas en:", out_dir, "\n")
