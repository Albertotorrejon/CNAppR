# Survival analysis stratified by CNA scores
#
# Dependencies (Suggests): survival, survminer
# Install: install.packages(c("survival", "survminer"))


#' Survival Analysis Stratified by CNA Scores
#'
#' Performs Kaplan-Meier survival analysis and log-rank testing stratifying
#' samples by quantile groups of a CNA score (FCS, BCS, or GCS). Optionally
#' conditions on a clinical grouping variable.
#'
#' @param scores Data frame of CNA scores returned by calculate_cna_scores()
#'   (rows = samples, columns include FCS, BCS, GCS). Row names must be sample
#'   identifiers matching those in survival_data.
#' @param survival_data Data frame with at least two columns: survival time and
#'   event status. Row names OR a column named by sample_col must correspond to
#'   the row names of scores.
#' @param score_var Character name of the score column to stratify on. One of
#'   "GCS", "FCS", or "BCS" (default "GCS").
#' @param time_col Character name of the survival time column in survival_data
#'   (default "time"). Values must be non-negative numeric.
#' @param event_col Character name of the event/status column in survival_data
#'   (default "event"). Values must be 0 (censored) or 1 (event).
#' @param sample_col Character name of the sample ID column in survival_data.
#'   If NULL (default), row names of survival_data are used for matching.
#' @param n_groups Integer number of quantile groups (default 2 = median split:
#'   High vs Low). Use 4 for quartile groups (Q1–Q4).
#' @param group_labels Character vector of labels for the quantile groups, in
#'   order from lowest to highest score. Length must equal n_groups. Default
#'   NULL uses "Low"/"High" for n_groups = 2, or "Q1"…"Q4" for n_groups = 4.
#' @param alpha Numeric significance level for the log-rank test (default 0.05).
#' @param plot Logical, generate a Kaplan-Meier plot with survminer
#'   (default TRUE).
#' @param palette Character colour palette for survminer (default "jco").
#' @param clinical_var Optional character name of a clinical grouping column in
#'   survival_data. When provided, separate KM curves are generated for each
#'   level of this variable (faceted analysis).
#'
#' @return A list containing:
#'   \describe{
#'     \item{fit}{survfit object (or named list of survfit objects in faceted
#'       mode).}
#'     \item{logrank}{Survival difference test result (survdiff or named list
#'       in faceted mode).}
#'     \item{pvalue}{Numeric log-rank p-value (or named numeric vector in
#'       faceted mode).}
#'     \item{plot}{ggplot2/survminer object (NULL if plot = FALSE).}
#'     \item{data}{Data frame used for analysis with columns: sample, time,
#'       event, score, score_group (and clinical_var if provided).}
#'   }
#'
#' @details
#' Requires the survival and survminer packages:
#' \code{install.packages(c("survival", "survminer"))}
#'
#' Samples present in scores but absent from survival_data (or vice versa)
#' are dropped with a warning. Samples with NA in time, event, or score are
#' also excluded.
#'
#' @examples
#' \dontrun{
#'   scores <- calculate_cna_scores(reseg_list)
#'   surv   <- data.frame(
#'     time   = c(24, 36, 12, 60, 48),
#'     event  = c(1, 0, 1, 0, 1),
#'     row.names = c("s1", "s2", "s3", "s4", "s5")
#'   )
#'
#'   # Median split on GCS
#'   res <- run_survival_analysis(scores, surv, score_var = "GCS")
#'   res$plot
#'   res$pvalue
#'
#'   # Quartile split on FCS, conditioned on MSI status
#'   res_q <- run_survival_analysis(scores, surv,
#'                                   score_var    = "FCS",
#'                                   n_groups     = 4,
#'                                   clinical_var = "msi_status")
#' }
#'
#' @export
run_survival_analysis <- function(scores,
                                   survival_data,
                                   score_var     = "GCS",
                                   time_col      = "time",
                                   event_col     = "event",
                                   sample_col    = NULL,
                                   n_groups      = 2L,
                                   group_labels  = NULL,
                                   alpha         = 0.05,
                                   plot          = TRUE,
                                   palette       = "jco",
                                   clinical_var  = NULL) {

  for (pkg in c("survival", "survminer")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required.\n",
           "Install with: install.packages(c('survival', 'survminer'))")
  }

  n_groups <- as.integer(n_groups)
  if (n_groups < 2L)
    stop("n_groups must be >= 2.")

  # ── validate score column ────────────────────────────────────────────────
  if (!score_var %in% colnames(scores))
    stop("score_var '", score_var, "' not found in scores data frame. ",
         "Available: ", paste(colnames(scores), collapse = ", "))

  # ── resolve sample IDs in survival_data ─────────────────────────────────
  if (!is.null(sample_col)) {
    if (!sample_col %in% colnames(survival_data))
      stop("sample_col '", sample_col, "' not found in survival_data.")
    surv_ids <- as.character(survival_data[[sample_col]])
  } else {
    surv_ids <- rownames(survival_data)
    if (is.null(surv_ids))
      stop("survival_data has no row names. ",
           "Provide sample IDs via the sample_col argument.")
  }

  # ── validate time and event columns ─────────────────────────────────────
  for (col in c(time_col, event_col)) {
    if (!col %in% colnames(survival_data))
      stop("Column '", col, "' not found in survival_data.")
  }

  # ── validate clinical_var ────────────────────────────────────────────────
  if (!is.null(clinical_var) && !clinical_var %in% colnames(survival_data))
    stop("clinical_var '", clinical_var, "' not found in survival_data.")

  # ── build group labels ───────────────────────────────────────────────────
  if (is.null(group_labels)) {
    group_labels <- if (n_groups == 2L) {
      c("Low", "High")
    } else if (n_groups == 4L) {
      paste0("Q", seq_len(4L))
    } else {
      paste0("G", seq_len(n_groups))
    }
  }
  if (length(group_labels) != n_groups)
    stop("group_labels must have length equal to n_groups (", n_groups, ").")

  # ── match samples ────────────────────────────────────────────────────────
  score_ids  <- rownames(scores)
  common_ids <- intersect(score_ids, surv_ids)

  extra_score <- setdiff(score_ids, surv_ids)
  extra_surv  <- setdiff(surv_ids,  score_ids)
  if (length(extra_score) > 0L)
    warning(length(extra_score),
            " sample(s) in scores have no survival data and are dropped.")
  if (length(extra_surv) > 0L)
    warning(length(extra_surv),
            " sample(s) in survival_data have no CNA scores and are dropped.")
  if (length(common_ids) == 0L)
    stop("No samples in common between scores and survival_data.")

  # Build analysis data frame
  surv_idx <- match(common_ids, surv_ids)
  adf <- data.frame(
    sample = common_ids,
    time   = as.numeric(survival_data[[time_col]][surv_idx]),
    event  = as.numeric(survival_data[[event_col]][surv_idx]),
    score  = as.numeric(scores[common_ids, score_var]),
    stringsAsFactors = FALSE
  )
  if (!is.null(clinical_var))
    adf$clinical_var <- as.character(survival_data[[clinical_var]][surv_idx])

  # Drop incomplete rows
  na_rows <- is.na(adf$time) | is.na(adf$event) | is.na(adf$score)
  if (any(na_rows)) {
    warning(sum(na_rows), " sample(s) dropped due to NA in time/event/score.")
    adf <- adf[!na_rows, , drop = FALSE]
  }
  if (nrow(adf) < 4L)
    stop("Too few samples after filtering (need >= 4).")

  # ── quantile-based grouping ──────────────────────────────────────────────
  breaks <- stats::quantile(adf$score,
                            probs = seq(0, 1, length.out = n_groups + 1L),
                            na.rm = TRUE)
  adf$score_group <- cut(adf$score, breaks = breaks,
                          labels = group_labels, include.lowest = TRUE)

  colnames(adf)[colnames(adf) == "score_group"] <- paste0(score_var, "_group")
  group_col <- paste0(score_var, "_group")

  # ── analysis function (single stratum or per-level) ──────────────────────
  .run_km <- function(df) {
    surv_obj <- survival::Surv(time = df$time, event = df$event)
    grp_fac  <- df[[group_col]]
    fit      <- survival::survfit(surv_obj ~ grp_fac, data = df)
    lr_test  <- survival::survdiff(surv_obj ~ grp_fac, data = df)
    pval     <- 1 - stats::pchisq(lr_test$chisq,
                                   df = length(lr_test$n) - 1L)
    list(fit = fit, logrank = lr_test, pvalue = pval)
  }

  # ── run analysis ─────────────────────────────────────────────────────────
  if (is.null(clinical_var)) {
    km <- .run_km(adf)
    result <- list(
      fit      = km$fit,
      logrank  = km$logrank,
      pvalue   = km$pvalue,
      plot     = NULL,
      data     = adf
    )

    if (plot) {
      pval_str <- if (km$pvalue < 0.001) "p < 0.001" else
        sprintf("p = %.3f", km$pvalue)
      result$plot <- survminer::ggsurvplot(
        km$fit,
        data          = adf,
        pval          = FALSE,
        conf.int      = FALSE,
        palette       = palette,
        legend.title  = score_var,
        title         = sprintf("KM — %s (%d groups, log-rank %s)",
                                score_var, n_groups, pval_str),
        xlab          = "Time",
        ylab          = "Survival probability",
        risk.table    = TRUE,
        ggtheme       = ggplot2::theme_bw()
      )
    }

  } else {
    # Faceted: run separately per clinical group
    levels_cv <- unique(na.omit(adf$clinical_var))
    km_list   <- list()
    pval_vec  <- numeric(length(levels_cv))
    names(pval_vec) <- levels_cv

    for (lv in levels_cv) {
      sub_df <- adf[!is.na(adf$clinical_var) & adf$clinical_var == lv, ]
      if (nrow(sub_df) < 4L) {
        warning("Group '", lv, "': too few samples (< 4), skipping.")
        next
      }
      km_list[[lv]] <- .run_km(sub_df)
      pval_vec[[lv]] <- km_list[[lv]]$pvalue
    }

    result <- list(
      fit     = lapply(km_list, `[[`, "fit"),
      logrank = lapply(km_list, `[[`, "logrank"),
      pvalue  = pval_vec,
      plot    = NULL,
      data    = adf
    )

    if (plot && length(km_list) > 0L) {
      # Use the first group for the representative plot
      first_lv  <- names(km_list)[1L]
      sub_df_p  <- adf[!is.na(adf$clinical_var) & adf$clinical_var == first_lv, ]
      pval_str  <- if (pval_vec[[first_lv]] < 0.001) "p < 0.001" else
        sprintf("p = %.3f", pval_vec[[first_lv]])
      result$plot <- survminer::ggsurvplot(
        km_list[[first_lv]]$fit,
        data         = sub_df_p,
        pval         = FALSE,
        conf.int     = FALSE,
        palette      = palette,
        legend.title = score_var,
        title        = sprintf("KM — %s, %s = %s (log-rank %s)",
                               score_var, clinical_var, first_lv, pval_str),
        xlab         = "Time",
        ylab         = "Survival probability",
        risk.table   = TRUE,
        ggtheme      = ggplot2::theme_bw()
      )
    }
  }

  result
}
