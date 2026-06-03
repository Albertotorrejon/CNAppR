#' @importFrom stats predict
NULL

#' Train a CNA-Based Random Forest Classifier
#'
#' Trains a Random Forest model to predict a clinical or molecular label from
#' CNA scores (FCS, BCS, GCS) or any numeric feature matrix derived from
#' \code{calculate_cna_scores()}.
#'
#' @param scores Data frame with at least columns \code{FCS}, \code{BCS},
#'   \code{GCS} and row names matching sample IDs. Additional numeric columns
#'   are included as features.
#' @param labels Named character or factor vector of class labels. Names must
#'   match row names of \code{scores}.
#' @param ntree Integer number of trees (default 500).
#' @param mtry Integer number of features tried at each split. Defaults to
#'   \code{floor(sqrt(ncol(features)))}.
#' @param seed Integer random seed for reproducibility (default 42).
#'
#' @return A list with:
#'   \describe{
#'     \item{model}{The trained \code{randomForest} object.}
#'     \item{oob_accuracy}{Out-of-bag accuracy (proportion correctly classified).}
#'     \item{importance}{Data frame of variable importance scores.}
#'   }
#'
#' @details
#' Requires the \code{randomForest} package (\code{install.packages("randomForest")}).
#' For a proper evaluation split samples into train/test before calling this
#' function:
#' \preformatted{
#'   idx   <- sample(nrow(scores), 0.7 * nrow(scores))
#'   train <- scores[idx, ]; test <- scores[-idx, ]
#'   lab_train <- labels[idx]; lab_test <- labels[-idx]
#'   mod <- train_cna_classifier(train, lab_train)
#'   preds <- predict_cna_class(test, mod$model)
#' }
#'
#' @examples
#' \dontrun{
#'   scores <- calculate_cna_scores(reseg_list)
#'   labels <- setNames(clinical$subtype, rownames(scores))
#'   mod    <- train_cna_classifier(scores, labels)
#'   mod$oob_accuracy
#' }
#'
#' @export
train_cna_classifier <- function(scores,
                                  labels,
                                  ntree = 500L,
                                  mtry  = NULL,
                                  seed  = 42L) {

  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' is required. Install with: install.packages('randomForest')")

  # Align labels to scores
  common <- intersect(rownames(scores), names(labels))
  if (length(common) == 0L)
    stop("No matching sample IDs between 'scores' row names and 'labels' names.")
  if (length(common) < nrow(scores))
    warning(nrow(scores) - length(common), " samples dropped (no matching label).")

  X <- scores[common, vapply(scores, is.numeric, logical(1L)), drop = FALSE]
  y <- factor(labels[common])

  if (is.null(mtry))
    mtry <- max(1L, floor(sqrt(ncol(X))))

  set.seed(seed)
  mod <- randomForest::randomForest(x = X, y = y, ntree = ntree, mtry = mtry,
                                    importance = TRUE)

  oob_acc <- 1 - mod$err.rate[ntree, "OOB"]
  imp_df  <- as.data.frame(randomForest::importance(mod))
  imp_df  <- imp_df[order(-imp_df$MeanDecreaseAccuracy), , drop = FALSE]

  list(model = mod, oob_accuracy = oob_acc, importance = imp_df)
}


#' Predict CNA-Based Sample Classification
#'
#' Applies a trained Random Forest model to new samples and returns predicted
#' class labels and class probabilities.
#'
#' @param scores Data frame of CNA scores for the samples to predict. Must
#'   contain the same feature columns used during training.
#' @param model A \code{randomForest} object (output of
#'   \code{train_cna_classifier()$model} or loaded with
#'   \code{load_prediction_model()}).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{ID}{Sample identifier (from row names of \code{scores}).}
#'     \item{predicted_class}{Predicted class label.}
#'     \item{...}{One column per class with the predicted probability.}
#'   }
#'
#' @examples
#' \dontrun{
#'   preds <- predict_cna_class(test_scores, mod$model)
#'   head(preds)
#' }
#'
#' @export
predict_cna_class <- function(scores, model) {

  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' is required. Install with: install.packages('randomForest')")

  if (!inherits(model, "randomForest"))
    stop("'model' must be a randomForest object.")

  X     <- scores[, vapply(scores, is.numeric, logical(1L)), drop = FALSE]
  preds <- stats::predict(model, X, type = "class")
  probs <- stats::predict(model, X, type = "prob")

  out <- data.frame(
    ID              = rownames(scores),
    predicted_class = as.character(preds),
    stringsAsFactors = FALSE
  )
  cbind(out, as.data.frame(probs))
}


#' Load a Saved Random Forest Model
#'
#' Loads a \code{randomForest} model previously saved with \code{saveRDS()} or
#' \code{save()}.
#'
#' @param path Character path to the model file (\code{.rds} or \code{.RData}).
#'
#' @return A \code{randomForest} object.
#'
#' @examples
#' \dontrun{
#'   # Save after training
#'   saveRDS(mod$model, "my_model.rds")
#'   # Load in a new session
#'   model <- load_prediction_model("my_model.rds")
#'   preds <- predict_cna_class(scores, model)
#' }
#'
#' @export
load_prediction_model <- function(path) {

  if (!file.exists(path))
    stop("Model file not found: ", path)

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    model <- readRDS(path)
  } else {
    env <- new.env()
    load(path, envir = env)
    objs <- ls(env)
    model <- get(objs[1L], envir = env)
  }

  if (!inherits(model, "randomForest"))
    warning("Loaded object does not appear to be a randomForest model.")

  model
}
