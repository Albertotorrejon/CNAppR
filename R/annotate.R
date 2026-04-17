#' Annotate Genes in CNA Segments
#'
#' For each segment classified as a CNA, identifies which reference genes
#' fall within its coordinates.
#'
#' @param segments Data frame with classified segments
#'   (output from resegment_sample()). Must contain: chr, loc.start,
#'   loc.end, classified, type.
#' @param genome_build Character "hg19" or "hg38" (default "hg19").
#' @param only_cna Logical. If TRUE (default), only annotates CNA segments
#'   (excludes segments classified as "diploid").
#'
#' @return The input data frame with one additional column:
#'   - genes: Character with gene names separated by "; ",
#'     or NA if no genes overlap the segment.
#'
#' @details
#' Reference files are RefSeq (UCSC) for hg19 and hg38, bundled with the
#' package. A gene is considered within a segment if its coordinates
#' overlap with the segment start and end positions.
#'
#' @examples
#' \dontrun{
#'   reseg <- resegment_sample(data, sample_id = "sample_01")
#'   reseg_ann <- annotate_cna_genes(reseg, genome_build = "hg19")
#'   head(reseg_ann[!is.na(reseg_ann$genes),
#'                  c("chr", "loc.start", "loc.end", "type", "genes")])
#' }
#'
#' @export
annotate_cna_genes <- function(segments, genome_build = "hg19",
                               only_cna = TRUE) {

  ref_file <- system.file("extdata",
                          paste0("refGene_", genome_build, ".txt"),
                          package = "CNAppR")

  if (!file.exists(ref_file)) {
    stop("Reference file not found for genome_build = '", genome_build,
         "'. Valid options: 'hg19', 'hg38'.")
  }

  ref_genes <- utils::read.table(ref_file, sep = "\t", header = TRUE,
                                 stringsAsFactors = FALSE)

  # Normalize chr column (remove "chr" prefix, convert X/Y to numeric)
  ref_genes$chrom <- gsub("chr", "", ref_genes$chrom)
  ref_genes$chrom[ref_genes$chrom == "X"] <- "23"
  ref_genes$chrom[ref_genes$chrom == "Y"] <- "24"
  ref_genes$chrom <- suppressWarnings(as.numeric(ref_genes$chrom))
  ref_genes <- ref_genes[!is.na(ref_genes$chrom), ]

  segments$genes <- NA_character_

  rows_to_annotate <- seq_len(nrow(segments))
  if (only_cna) {
    rows_to_annotate <- which(!is.na(segments$classified) &
                                segments$classified != "diploid")
  }

  for (i in rows_to_annotate) {
    chr_val <- suppressWarnings(
      as.numeric(gsub("chr", "", segments$chr[i]))
    )
    if (is.na(chr_val)) next

    seg_start <- segments$loc.start[i]
    seg_end   <- segments$loc.end[i]

    genes_in_seg <- ref_genes[
      ref_genes$chrom == chr_val &
        ref_genes$end   >= seg_start &
        ref_genes$start <= seg_end,
      "symbol_name"
    ]

    genes_in_seg <- unique(
      genes_in_seg[!is.na(genes_in_seg) & genes_in_seg != ""]
    )

    if (length(genes_in_seg) > 0) {
      segments$genes[i] <- paste(genes_in_seg, collapse = "; ")
    }
  }

  segments
}
