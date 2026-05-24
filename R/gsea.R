# Gene Set Enrichment Analysis (GSEA) for CNA data
#
# Dependencies (Suggests): fgsea (Bioconductor), msigdbr (CRAN)
# Install: BiocManager::install("fgsea"); install.packages("msigdbr")


# ─── internal helpers ────────────────────────────────────────────────────────

# Retrieves gene sets from MSigDB via msigdbr, handles API differences across
# msigdbr versions (≥7.5.1 uses category=; 2024+ may use collection=).
.build_msigdb_genesets <- function(collections, species) {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required.\nInstall with: install.packages('msigdbr')")

  gs_list <- list()

  for (coll in collections) {
    db <- tryCatch(
      msigdbr::msigdbr(species = species, category = coll),
      error = function(e) {
        warning("Could not fetch MSigDB collection '", coll,
                "': ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(db) || nrow(db) == 0L) {
      warning("No gene sets returned for MSigDB collection '", coll, "'.")
      next
    }

    # Column names differ across msigdbr versions — pick defensively
    name_col <- intersect(c("gs_name", "gs_id"), colnames(db))[1L]
    sym_col  <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(db))[1L]

    if (is.na(name_col) || is.na(sym_col)) {
      warning("Unexpected msigdbr column layout for collection '", coll,
              "'. Skipping.")
      next
    }

    coll_gs <- split(db[[sym_col]], db[[name_col]])
    gs_list <- c(gs_list, coll_gs)
  }

  gs_list
}

# NES barplot: top_n pathways ranked by |NES|, coloured by direction.
.plot_gsea_results <- function(res_df, top_n = 20L, alpha = 0.25,
                               collections = character(0)) {
  sig <- res_df[!is.na(res_df$padj) & res_df$padj <= alpha, , drop = FALSE]

  if (nrow(sig) == 0L) {
    message("No significant pathways at padj <= ", alpha,
            "; showing top ", top_n, " by |NES|.")
    sig <- res_df
  }

  top <- sig[order(abs(sig$NES), decreasing = TRUE), , drop = FALSE]
  top <- top[seq_len(min(as.integer(top_n), nrow(top))), , drop = FALSE]
  top <- top[order(top$NES), , drop = FALSE]

  top$label <- ifelse(
    nchar(top$pathway) > 55L,
    paste0(substr(top$pathway, 1L, 52L), "..."),
    top$pathway
  )
  top$label     <- factor(top$label, levels = top$label)
  top$direction <- ifelse(top$NES > 0, "Amplification", "Deletion")

  coll_str <- if (length(collections) > 0L)
    paste("Collections:", paste(collections, collapse = ", ")) else NULL

  ggplot2::ggplot(top, ggplot2::aes(
    x    = .data$NES,
    y    = .data$label,
    fill = .data$direction
  )) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_fill_manual(
      values = c(Amplification = "#d73027", Deletion = "#4575b4")
    ) +
    ggplot2::labs(
      title    = sprintf("GSEA — top %d pathways by |NES|", nrow(top)),
      subtitle = coll_str,
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL,
      fill     = NULL
    ) +
    ggplot2::theme_bw(base_size = 11L) +
    ggplot2::theme(
      legend.position  = "top",
      axis.text.y      = ggplot2::element_text(size = 9L),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Runs fgseaMultilevel on a single gene rank vector and returns a tidy df.
.run_fgsea_single <- function(gene_ranks, gs_list, min_size, max_size, eps) {
  res <- fgsea::fgseaMultilevel(
    pathways = gs_list,
    stats    = gene_ranks,
    minSize  = as.integer(min_size),
    maxSize  = as.integer(max_size),
    eps      = eps
  )
  res_df <- as.data.frame(res[order(res$padj), , drop = FALSE])
  res_df$leadingEdge <- vapply(
    res$leadingEdge,
    function(x) paste(x, collapse = "; "),
    character(1L)
  )
  res_df
}


# ─── build_gene_ranks ────────────────────────────────────────────────────────

#' Build a gene rank vector from annotated CNA segments
#'
#' Converts the output of annotate_cna_genes() into a named numeric vector
#' suitable as input to run_gsea(). Each gene receives the score of the CNA
#' segment it overlaps; when a gene appears in multiple segments the value
#' with the highest absolute magnitude is kept by default.
#'
#' @param annotated_segs Data frame returned by annotate_cna_genes().
#'   Must contain a "genes" column (semicolon-separated HGNC symbols) and
#'   the score column specified by score_col.
#' @param score_col Character name of the numeric score column to use as
#'   the ranking metric. Default "seg.mean" (log2 copy-number ratio).
#'   Other useful options: "score" (CNAppR weight score), or any numeric
#'   column present in the data.
#' @param agg_fun Function applied when a gene appears in more than one
#'   segment. Receives a numeric vector and returns a scalar. Default:
#'   keep the value with the highest absolute magnitude.
#'
#' @return Named numeric vector sorted in decreasing order. Positive values
#'   indicate amplification; negative values indicate deletion.
#'   Names are HGNC gene symbols.
#'
#' @examples
#' \dontrun{
#'   reseg  <- resegment_sample(data, sample_id = "s1")
#'   ann    <- annotate_cna_genes(reseg)
#'   ranks  <- build_gene_ranks(ann)
#'   gsea_r <- run_gsea(ranks, collections = "H")
#' }
#'
#' @export
build_gene_ranks <- function(annotated_segs, score_col = "seg.mean",
                              agg_fun = function(x) x[which.max(abs(x))]) {
  if (!score_col %in% colnames(annotated_segs))
    stop("Column '", score_col, "' not found in annotated_segs.")
  if (!"genes" %in% colnames(annotated_segs))
    stop("Column 'genes' not found. Run annotate_cna_genes() first.")

  has_genes <- !is.na(annotated_segs$genes) & nzchar(annotated_segs$genes)
  segs_ann  <- annotated_segs[has_genes, , drop = FALSE]

  if (nrow(segs_ann) == 0L)
    stop("No annotated segments found. Check annotate_cna_genes() output.")

  # Expand: one row per gene-segment pair
  gene_rows <- vector("list", nrow(segs_ann))
  for (i in seq_len(nrow(segs_ann))) {
    genes <- trimws(strsplit(segs_ann$genes[i], ";")[[1L]])
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0L) next
    gene_rows[[i]] <- data.frame(
      gene  = genes,
      score = segs_ann[[score_col]][i],
      stringsAsFactors = FALSE
    )
  }

  gene_df <- do.call(rbind, gene_rows[!vapply(gene_rows, is.null, logical(1L))])
  if (is.null(gene_df) || nrow(gene_df) == 0L)
    stop("No gene-score pairs could be built.")

  scores     <- tapply(gene_df$score, gene_df$gene, agg_fun)
  gene_ranks <- sort(unlist(scores, use.names = TRUE), decreasing = TRUE)
  gene_ranks
}


# ─── run_gsea ────────────────────────────────────────────────────────────────

#' Gene Set Enrichment Analysis on CNA-ranked genes
#'
#' Runs fgsea::fgseaMultilevel() using gene sets from MSigDB (via msigdbr)
#' on a vector of genes ranked by their CNA score. Optionally stratifies
#' analysis by clinical group when a list of per-sample rank vectors is
#' provided.
#'
#' @param gene_ranks Named numeric vector (names = HGNC gene symbols, values =
#'   CNA score; use build_gene_ranks() to produce it from annotate_cna_genes()
#'   output). For stratified analysis, a named list of such vectors (one per
#'   sample); see clinical_groups.
#' @param collections Character vector of MSigDB collection codes to test.
#'   Default c("H", "C2"). Common options:
#'   "H" (Hallmarks), "C2" (curated), "C5" (GO), "C6" (oncogenic).
#'   Full list: msigdbr::msigdbr_collections().
#' @param species Character species name for msigdbr (default "Homo sapiens").
#' @param min_size Integer minimum gene set size after filtering (default 15).
#' @param max_size Integer maximum gene set size after filtering (default 500).
#' @param eps Numeric boundary for fgseaMultilevel p-value precision
#'   (default 1e-10). Smaller = more accurate but slower.
#' @param alpha Numeric FDR threshold for significant pathways (default 0.25,
#'   per GSEA convention for exploratory analyses).
#' @param top_n Integer number of top pathways to include in the plot
#'   (by |NES|, default 20).
#' @param plot Logical, generate a NES barplot of the top pathways
#'   (default TRUE).
#' @param clinical_groups Named character vector mapping sample IDs to group
#'   labels (e.g. c(sample_1 = "GroupA", sample_2 = "GroupB")). Only used
#'   when gene_ranks is a named list of per-sample rank vectors. Within each
#'   group, per-gene scores are averaged before GSEA is run.
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{Data frame with columns pathway, NES, pval, padj, size,
#'       leadingEdge (semicolon-separated). In stratified mode, also includes
#'       a "group" column.}
#'     \item{significant}{Subset of results where padj <= alpha.}
#'     \item{plot}{ggplot2 object (NULL if plot = FALSE).}
#'   }
#'
#' @details
#' Requires fgsea (Bioconductor) and msigdbr (CRAN):
#' \code{BiocManager::install("fgsea"); install.packages("msigdbr")}
#'
#' **Recommended pipeline:**
#' \enumerate{
#'   \item resegment_sample() → per-sample CNA segments
#'   \item annotate_cna_genes() → add gene symbols per segment
#'   \item build_gene_ranks() → named numeric rank vector
#'   \item run_gsea() → enrichment results + plot
#' }
#'
#' **Stratified analysis:** pass gene_ranks as a named list (one vector per
#' sample) and provide clinical_groups to run GSEA per clinical subgroup.
#' Within each group, per-gene scores are averaged across samples before
#' fgseaMultilevel is called.
#'
#' @examples
#' \dontrun{
#'   # Single-sample GSEA
#'   reseg <- resegment_sample(data, sample_id = "s1")
#'   ann   <- annotate_cna_genes(reseg)
#'   ranks <- build_gene_ranks(ann)
#'   res   <- run_gsea(ranks, collections = c("H", "C2"), alpha = 0.25)
#'   res$plot
#'   head(res$significant)
#'
#'   # Stratified GSEA by clinical group
#'   sample_ids <- unique(data$ID)
#'   ranks_list <- lapply(sample_ids, function(id) {
#'     build_gene_ranks(annotate_cna_genes(resegment_sample(data, id)))
#'   })
#'   names(ranks_list) <- sample_ids
#'
#'   groups <- c(s1 = "MSI", s2 = "MSS", s3 = "MSI")  # example
#'   res_strat <- run_gsea(ranks_list, clinical_groups = groups, collections = "H")
#' }
#'
#' @export
run_gsea <- function(gene_ranks,
                     collections     = c("H", "C2"),
                     species         = "Homo sapiens",
                     min_size        = 15L,
                     max_size        = 500L,
                     eps             = 1e-10,
                     alpha           = 0.25,
                     top_n           = 20L,
                     plot            = TRUE,
                     clinical_groups = NULL) {

  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("Package 'fgsea' is required.\n",
         "Install with: BiocManager::install('fgsea')")
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required.\n",
         "Install with: install.packages('msigdbr')")

  # ── validate gene_ranks ──────────────────────────────────────────────────
  is_list_input <- is.list(gene_ranks) && !is.data.frame(gene_ranks)

  if (!is_list_input) {
    if (!is.numeric(gene_ranks) || is.null(names(gene_ranks)))
      stop("gene_ranks must be a named numeric vector or a named list of ",
           "named numeric vectors.")
  } else {
    if (is.null(names(gene_ranks)))
      stop("When gene_ranks is a list, it must be named (one entry per sample).")
  }

  # ── fetch gene sets ──────────────────────────────────────────────────────
  message("Fetching gene sets from MSigDB (",
          paste(collections, collapse = ", "), ")...")
  gs_list <- .build_msigdb_genesets(collections, species)

  if (length(gs_list) == 0L)
    stop("No gene sets retrieved from MSigDB for collections: ",
         paste(collections, collapse = ", "))

  message("  Gene sets loaded: ", length(gs_list))

  # ── stratified mode ──────────────────────────────────────────────────────
  if (is_list_input && !is.null(clinical_groups)) {
    groups <- unique(as.character(clinical_groups))
    all_res <- vector("list", length(groups))
    names(all_res) <- groups

    for (grp in groups) {
      grp_samples <- names(clinical_groups)[clinical_groups == grp]
      grp_samples <- intersect(grp_samples, names(gene_ranks))

      if (length(grp_samples) == 0L) {
        warning("Group '", grp, "': no matching samples in gene_ranks. Skipping.")
        next
      }

      # Average per-gene scores across group samples
      grp_ranks_list <- gene_ranks[grp_samples]
      all_genes <- unique(unlist(lapply(grp_ranks_list, names)))
      grp_vec <- vapply(all_genes, function(g) {
        vals <- vapply(grp_ranks_list,
                       function(r) if (g %in% names(r)) r[[g]] else NA_real_,
                       numeric(1L))
        mean(vals, na.rm = TRUE)
      }, numeric(1L))
      grp_vec <- sort(grp_vec[is.finite(grp_vec)], decreasing = TRUE)

      message("Running GSEA for group '", grp, "' (",
              length(grp_samples), " samples, ", length(grp_vec), " genes)...")
      res_i       <- .run_fgsea_single(grp_vec, gs_list, min_size, max_size, eps)
      res_i$group <- grp
      all_res[[grp]] <- res_i
    }

    combined <- do.call(rbind, all_res[!vapply(all_res, is.null, logical(1L))])
    sig      <- combined[!is.na(combined$padj) & combined$padj <= alpha, ,
                         drop = FALSE]
    p <- NULL
    if (plot && nrow(combined) > 0L) {
      first_grp <- combined[combined$group == groups[1L], , drop = FALSE]
      p <- .plot_gsea_results(first_grp, top_n, alpha, collections)
    }
    return(list(results = combined, significant = sig, plot = p))
  }

  # ── single-vector mode ───────────────────────────────────────────────────
  if (is_list_input) {
    # List provided but no clinical_groups: pool all samples (average per gene)
    warning("gene_ranks is a list but clinical_groups is NULL. ",
            "Averaging ranks across all samples for a single GSEA run.")
    all_genes <- unique(unlist(lapply(gene_ranks, names)))
    pooled    <- vapply(all_genes, function(g) {
      vals <- vapply(gene_ranks,
                     function(r) if (g %in% names(r)) r[[g]] else NA_real_,
                     numeric(1L))
      mean(vals, na.rm = TRUE)
    }, numeric(1L))
    gene_ranks <- sort(pooled[is.finite(pooled)], decreasing = TRUE)
  }

  # Remove NAs; resolve duplicates (keep max |score|)
  gene_ranks <- gene_ranks[!is.na(gene_ranks)]
  if (anyDuplicated(names(gene_ranks))) {
    gene_ranks <- tapply(gene_ranks, names(gene_ranks),
                         function(x) x[which.max(abs(x))])
    gene_ranks <- sort(unlist(gene_ranks, use.names = TRUE), decreasing = TRUE)
  } else {
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  }

  message("Running fgseaMultilevel (", length(gs_list), " gene sets, ",
          length(gene_ranks), " genes)...")
  res_df <- .run_fgsea_single(gene_ranks, gs_list, min_size, max_size, eps)

  sig <- res_df[!is.na(res_df$padj) & res_df$padj <= alpha, , drop = FALSE]
  message("Significant pathways (padj <= ", alpha, "): ", nrow(sig))

  p <- NULL
  if (plot && nrow(res_df) > 0L) {
    p <- .plot_gsea_results(res_df, top_n, alpha, collections)
  }

  list(results = res_df, significant = sig, plot = p)
}
