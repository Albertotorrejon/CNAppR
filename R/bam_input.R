# BAM → CBS → CNAppR pipeline
#
# Flujo completo: BAM file → QC + conteo de reads → segmentación CBS →
#   formato CNAppR → resegment_sample()
#
# Paquetes Bioconductor requeridos: Rsamtools, GenomicRanges, DNAcopy
# Instalación: BiocManager::install(c("Rsamtools", "GenomicRanges", "DNAcopy"))
#
# Illumina (exoma): paired-end reads, bins más pequeños (~500 kb)
# Nanopore (low-pass WGS): single long reads, bins más grandes (~1 Mb)


# ─── helpers internos ────────────────────────────────────────────────────────

.check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Paquete requerido no instalado: ", pkg,
         "\nInstalar con: BiocManager::install('", pkg, "')")
}

.check_bam <- function(bam_path) {
  if (!file.exists(bam_path))
    stop("Archivo BAM no encontrado: ", bam_path)
  bai_path <- paste0(bam_path, ".bai")
  if (!file.exists(bai_path))
    stop("Índice BAM (.bai) no encontrado.\nEjecuta: samtools index ",
         basename(bam_path))
}
  
# Normaliza nombres de cromosoma al estilo "chrN" (con prefijo chr)
.norm_chr_to_ucsc <- function(x) {
  x <- as.character(x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

# Normaliza nombres de cromosoma a entero (sin prefijo; X→23, Y→24)
.norm_chr_to_int <- function(x) {
  x <- gsub("^chr", "", as.character(x))
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  suppressWarnings(as.integer(x))
}


# ─── read_bam_qc ─────────────────────────────────────────────────────────────

#' Conteo de reads por bin desde un archivo BAM
#'
#' Divide el genoma en bins de tamaño fijo, cuenta los reads de cada bin
#' usando Rsamtools y calcula el log2 ratio normalizado por mediana.
#' Este es el input necesario para run_cbs_segmentation().
#'
#' @param bam_path Character ruta al archivo BAM (debe estar indexado, .bai).
#' @param bin_size Integer tamaño de bin en pb (default 500000 para exoma,
#'   usar 1000000 para low-pass Nanopore).
#' @param genome_build Character "hg19" o "hg38" (default "hg19"). Usado
#'   solo como metadato; los tamaños de cromosoma se leen del header del BAM.
#' @param min_mapq Integer calidad mínima de mapeo a considerar (default 20).
#'
#' @return Data frame con columnas: chr, loc.start, loc.end, reads, log2_ratio.
#'   El log2_ratio está normalizado por la mediana de bins con reads > 0.
#'
#' @details
#' Requiere: Rsamtools, GenomicRanges (BiocManager::install(...))
#'
#' Para Illumina (exoma): bin_size = 500000 (500 kb)
#' Para Nanopore low-pass: bin_size = 1000000 (1 Mb)
#'
#' Nota: La corrección de GC content no está implementada en esta versión.
#' Es recomendable aplicarla en datos con sesgos fuertes de GC.
#'
#' @examples
#' \dontrun{
#'   qc <- read_bam_qc("sample.bam", bin_size = 1e6, genome_build = "hg19")
#'   head(qc)
#' }
#'
#' @export
read_bam_qc <- function(bam_path, bin_size = 500000L, genome_build = "hg19",
                         min_mapq = 20L) {
  .check_bam(bam_path)
  .check_pkg("Rsamtools")
  .check_pkg("GenomicRanges")

  # Tamaños de cromosoma desde el header del BAM (no depende de BSgenome)
  header   <- Rsamtools::scanBamHeader(bam_path)[[1]]
  chr_lens <- header$targets   # named integer: chr_name → length

  # Quedarse solo con cromosomas estándar (1-22, X, Y), con o sin prefijo "chr"
  keep     <- grepl("^(chr)?([0-9]{1,2}|X|Y)$", names(chr_lens))
  chr_lens <- chr_lens[keep]
  if (length(chr_lens) == 0)
    stop("No se encontraron cromosomas estándar en el header del BAM.")

  # Crear bins de tamaño fijo en todo el genoma
  bins <- GenomicRanges::tileGenome(
    seqlengths             = chr_lens,
    tilewidth              = as.integer(bin_size),
    cut.last.tile.in.chrom = TRUE
  )

  # Contar reads por bin (solo reads con MAPQ >= min_mapq)
  param  <- Rsamtools::ScanBamParam(
    which       = bins,
    mapqFilter  = as.integer(min_mapq)
  )
  counts <- Rsamtools::countBam(bam_path, param = param)

  reads <- counts$records
  total <- sum(reads)
  if (total == 0)
    stop("No se contaron reads en '", basename(bam_path),
         "'. Verifica que el BAM no esté vacío y el índice sea válido.")

  # Log2 ratio vs mediana de bins no vacíos (más robusto que vs media total)
  median_reads <- stats::median(reads[reads > 0])
  log2_ratio   <- log2(pmax(reads, 0.5) / median_reads)

  data.frame(
    chr        = as.character(GenomicRanges::seqnames(bins)),
    loc.start  = GenomicRanges::start(bins),
    loc.end    = GenomicRanges::end(bins),
    reads      = reads,
    log2_ratio = log2_ratio,
    stringsAsFactors = FALSE
  )
}


# ─── run_cbs_segmentation ────────────────────────────────────────────────────

#' Segmentación CBS (Circular Binary Segmentation)
#'
#' Ejecuta CBS sobre los log2 ratios producidos por read_bam_qc() y devuelve
#' una tabla en el formato de input de resegment_sample().
#'
#' @param qc_data Data frame output de read_bam_qc().
#' @param sample_id Character identificador de la muestra.
#' @param alpha Numeric nivel de significación CBS (default 0.01).
#' @param nperm Integer número de permutaciones CBS (default 10000).
#' @param sequencing_type Character "illumina" o "nanopore". Ajusta el
#'   suavizado de outliers: más agresivo para long reads de Nanopore.
#'
#' @return Data frame con columnas: ID, chr, loc.start, loc.end, seg.mean.
#'   Listo para pasar directamente a resegment_sample().
#'
#' @details
#' Requiere: DNAcopy (BiocManager::install("DNAcopy"))
#'
#' CBS es el algoritmo de segmentación más usado mundialmente para detectar
#' log ratios en datos de copy number. Implementado en el paquete R DNAcopy.
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
    stop("qc_data debe tener columnas: chr, loc.start, loc.end, log2_ratio.",
         "\nUsa read_bam_qc() para generarlo.")

  # Normalizar cromosomas a entero para DNAcopy
  chr_int <- .norm_chr_to_int(qc_data$chr)
  valid   <- !is.na(chr_int)
  if (!any(valid)) stop("No hay cromosomas válidos en qc_data.")

  qc_data <- qc_data[valid, ]
  chr_int <- chr_int[valid]

  # Posición central de cada bin como coordenada para CBS
  pos <- as.integer((qc_data$loc.start + qc_data$loc.end) / 2L)

  # Objeto CNA para DNAcopy
  cna_obj <- DNAcopy::CNA(
    genomdat  = qc_data$log2_ratio,
    chrom     = chr_int,
    maploc    = pos,
    data.type = "logratio",
    sampleid  = sample_id
  )

  # Suavizado de outliers (más agresivo en Nanopore por ruido de long reads)
  smooth_region <- if (sequencing_type == "nanopore") 5L else 2L
  cna_smooth    <- DNAcopy::smooth.CNA(cna_obj, smooth.region = smooth_region)

  # Segmentación
  seg    <- DNAcopy::segment(cna_smooth,
                              alpha   = alpha,
                              nperm   = as.integer(nperm),
                              verbose = 0)
  seg_df <- seg$output

  if (nrow(seg_df) == 0)
    stop("CBS no generó segmentos para la muestra '", sample_id, "'.")

  # Formato CNAppR: ID, chr, loc.start, loc.end, seg.mean
  data.frame(
    ID        = sample_id,
    chr       = as.integer(seg_df$chrom),
    loc.start = as.integer(seg_df$loc.start),
    loc.end   = as.integer(seg_df$loc.end),
    seg.mean  = as.numeric(seg_df$seg.mean),
    stringsAsFactors = FALSE
  )
}


# ─── run_bam_pipeline ────────────────────────────────────────────────────────

#' Pipeline completo BAM → CBS → CNAppR
#'
#' Wrapper de convenncia que ejecuta read_bam_qc(), run_cbs_segmentation()
#' y resegment_sample() en una sola llamada.
#'
#' @param bam_path Character ruta al archivo BAM indexado.
#' @param sample_id Character identificador de la muestra.
#' @param sequencing_type Character "illumina" o "nanopore".
#' @param bin_size Integer tamaño de bin en pb para el conteo de reads.
#'   Default 500000 (Illumina) o 1000000 (Nanopore).
#' @param genome_build Character "hg19" o "hg38" (default "hg19").
#' @param min_mapq Integer calidad mínima de mapeo (default 20).
#' @param alpha Numeric nivel de significación CBS (default 0.01).
#' @param nperm Integer permutaciones CBS (default 10000).
#' @param ... Argumentos adicionales pasados a resegment_sample()
#'   (p.ej. purity, acrocentric, skip_resegmentation).
#'
#' @return Output de resegment_sample(): data frame con segmentos CNA
#'   clasificados (chromosomal, arm, focal, diploid).
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
#'   # Illumina exoma
#'   result <- run_bam_pipeline(
#'     "sample_illumina.bam",
#'     sample_id       = "sample_02",
#'     sequencing_type = "illumina",
#'     bin_size        = 500000
#'   )
#' }
#'
#' @export
run_bam_pipeline <- function(bam_path,
                              sample_id,
                              sequencing_type = c("illumina", "nanopore"),
                              bin_size        = NULL,
                              genome_build    = "hg19",
                              min_mapq        = 20L,
                              alpha           = 0.01,
                              nperm           = 10000L,
                              ...) {
  sequencing_type <- match.arg(sequencing_type)

  # Bin size por defecto según tipo de secuenciación
  if (is.null(bin_size)) {
    bin_size <- if (sequencing_type == "nanopore") 1000000L else 500000L
  }

  message("Paso 1/3: Contando reads por bin desde BAM...")
  qc_data <- read_bam_qc(
    bam_path     = bam_path,
    bin_size     = bin_size,
    genome_build = genome_build,
    min_mapq     = min_mapq
  )

  message("Paso 2/3: Segmentación CBS (", sequencing_type, ")...")
  seg_data <- run_cbs_segmentation(
    qc_data         = qc_data,
    sample_id       = sample_id,
    alpha           = alpha,
    nperm           = nperm,
    sequencing_type = sequencing_type
  )

  message("Paso 3/3: Resegmentación y clasificación CNAppR...")
  resegment_sample(seg_data, sample_id = sample_id, ...)
}
