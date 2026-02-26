#' Calculate CNA Scores
#'
#' Computes FCS (Focal Copy Number Score), BCS (Broad Copy Number Score),
#' and GCS (Global Copy Number Score) for samples.
#'
#' @param resegmented_list List of re-segmented data frames, one per sample (output from resegment_sample()).
#' @param sample_names Character vector of sample IDs (optional; inferred from names if list is named).
#'
#' @return Data frame with columns: FCS, BCS, GCS and row names as sample_names.
#'
#' @details
#' Scores are calculated as:
#' - **FCS** = count of focal alterations
#' - **BCS** = count of chromosomal alterations + count of arm-level alterations
#' - **GCS** = (FCS - mean(FCS)) / sd(FCS) + (BCS - mean(BCS)) / sd(BCS)
#'
#' If all scores are 0, normalization is skipped.
#'
#' @examples
#' \dontrun{
#'   reseg_list <- lapply(c("sample_1", "sample_2"), function(id) {
#'     resegment_sample(data, sample_id = id)
#'   })
#'   names(reseg_list) <- c("sample_1", "sample_2")
#'
#'   scores <- calculate_cna_scores(reseg_list)
#'   head(scores)
#' }
#'
#' @export
calculate_cna_scores <- function(resegmented_list, sample_names = NULL) {
  if (is.null(sample_names)) {
    sample_names <- names(resegmented_list)
    if (is.null(sample_names)) {
      sample_names <- paste0("sample_", seq_along(resegmented_list))
    }
  }

  if (length(sample_names) != length(resegmented_list)) {
    stop("Length of sample_names must match length of resegmented_list.")
  }

  # Initialize score vectors
  fcs <- numeric(length(resegmented_list))
  bcs <- numeric(length(resegmented_list))

  # Extract counts from each resegmented sample
  for (i in seq_along(resegmented_list)) {
    reseg <- resegmented_list[[i]]

    # Count alterations by type
    n_focal <- sum(reseg$classified == "Focal", na.rm = TRUE)
    n_chromosomal <- sum(reseg$classified == "Chromosomal", na.rm = TRUE)
    n_arm <- sum(reseg$classified == "Arm", na.rm = TRUE)

    fcs[i] <- n_focal
    bcs[i] <- n_chromosomal + n_arm
  }

  # Normalize scores
  norm_fcs <- .normalize_score(fcs)
  norm_bcs <- .normalize_score(bcs)

  # Global score
  gcs <- norm_fcs + norm_bcs

  scores <- data.frame(
    FCS = fcs,
    BCS = bcs,
    GCS = gcs,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  return(scores)
}


#' Internal: Normalize Score Vector
#'
#' Standardizes a score vector (mean 0, sd 1) if not all zeros.
#'
#' @param score Numeric vector of scores.
#'
#' @return Numeric vector (normalized or original if all zeros).
#'
#' @keywords internal
.normalize_score <- function(score) {
  if (length(which(score == 0)) == length(score)) {
    # All zeros: no normalization
    return(score)
  } else {
    # Standard normalization
    return((score - mean(score)) / stats::sd(score))
  }
}


#' Test Clinical Association
#'
#' Performs parametric (correlation) and non-parametric (Kruskal-Wallis, Wilcoxon)
#' tests between CNA scores and clinical variables.
#'
#' @param scores Data frame with CNA scores (rows = samples, cols = FCS, BCS, GCS, etc.).
#' @param clinical_data Data frame with clinical variables (rows = samples, cols = variables).
#' @param method Character "parametric", "nonparametric", or "both" (default).
#'
#' @return List containing:
#'   - pval_parametric: Data frame of p-values (parametric tests)
#'   - pval_nonparametric: Data frame of p-values (non-parametric tests)
#'
#' @details
#' - **Parametric**: Pearson correlation (numeric vs numeric), ANOVA (categorical vs numeric)
#' - **Non-parametric**: Spearman correlation (numeric vs numeric), Kruskal-Wallis (categorical vs numeric)
#'
#' @examples
#' \dontrun{
#'   clinical_assoc <- test_clinical_association(scores, clinical_data)
#'   print(clinical_assoc$pval_parametric)
#' }
#'
#' @export
test_clinical_association <- function(scores, clinical_data, method = "both") {
  if (is.null(clinical_data) || nrow(clinical_data) == 0) {
    warning("No clinical data provided.")
    return(list(pval_parametric = NULL, pval_nonparametric = NULL))
  }

  # Match rows
  common_samples <- intersect(rownames(scores), rownames(clinical_data))
  scores <- scores[common_samples, , drop = FALSE]
  clinical_data <- clinical_data[common_samples, , drop = FALSE]

  # Initialize p-value matrices
  pval_parametric <- data.frame(matrix(NA, nrow = ncol(clinical_data), ncol = ncol(scores)))
  colnames(pval_parametric) <- colnames(scores)
  rownames(pval_parametric) <- colnames(clinical_data)

  pval_nonparametric <- pval_parametric

  # Test each clinical variable against each score
  for (j in seq_len(ncol(clinical_data))) {
    var_name <- colnames(clinical_data)[j]
    groups <- clinical_data[, j]

    # Determine variable class
    var_class <- class(groups)

    for (k in seq_len(ncol(scores))) {
      score_col <- scores[, k]

      # Remove NA pairs
      valid_idx <- !is.na(groups) & !is.na(score_col)
      if (sum(valid_idx) < 4) next # Require at least 4 observations

      groups_valid <- groups[valid_idx]
      score_valid <- score_col[valid_idx]

      # Numeric vs Numeric: Correlation
      if (var_class %in% c("numeric", "integer")) {
        groups_num <- as.numeric(groups_valid)

        if (method %in% c("parametric", "both")) {
          res <- stats::cor.test(groups_num, score_valid, method = "pearson")
          pval_parametric[j, k] <- res$p.value
        }

        if (method %in% c("nonparametric", "both")) {
          res <- stats::cor.test(groups_num, score_valid, method = "spearman")
          pval_nonparametric[j, k] <- res$p.value
        }
      }

      # Categorical vs Numeric: Kruskal-Wallis / ANOVA
      if (var_class %in% c("character", "factor")) {
        groups_factor <- as.factor(groups_valid)

        if (method %in% c("nonparametric", "both")) {
          res <- stats::kruskal.test(score_valid ~ groups_factor)
          pval_nonparametric[j, k] <- res$p.value
        }

        if (method %in% c("parametric", "both")) {
          # ANOVA
          aov_res <- stats::aov(score_valid ~ groups_factor)
          pval_parametric[j, k] <- stats::anova(aov_res)$`Pr(>F)`[1]
        }
      }
    }
  }

  return(list(
    pval_parametric = pval_parametric,
    pval_nonparametric = pval_nonparametric
  ))
}
