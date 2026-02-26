#' Load Prediction Model
#'
#' Loads a pre-trained random forest model for CNA-based classification.
#'
#' @param model_name Character model file name (without .RData extension).
#'   Should be located in inst/models/ directory.
#'
#' @return A random forest model object.
#'
#' @details
#' Models should be trained and saved as .RData files containing a
#' variable named "model_cms" or equivalent.
#'
#' @keywords internal
load_prediction_model <- function(model_name) {
  model_path <- system.file("models", paste0(model_name, ".RData"), package = "CNAppR")

  if (!file.exists(model_path)) {
    warning("Model file not found: ", model_path)
    return(NULL)
  }

  env <- new.env()
  load(model_path, envir = env)

  # Try common model names
  if ("model_cms" %in% ls(env)) {
    return(env$model_cms)
  } else if ("model" %in% ls(env)) {
    return(env$model)
  } else {
    # Return first object that looks like a model
    objs <- ls(env)
    if (length(objs) > 0) {
      return(get(objs[1], envir = env))
    }
  }

  warning("Could not find model in file: ", model_path)
  return(NULL)
}


#' Predict CNA-Based Classification
#'
#' Uses a random forest model to predict sample classification from CNA scores.
#'
#' @param scores Data frame with FCS, BCS, GCS scores (from calculate_cna_scores()).
#' @param model_name Character model file name.
#'
#' @return Data frame with:
#'   - ID: sample identifier
#'   - predicted_class: predicted classification
#'   - probability: prediction probability (if available from model)
#'
#' @details
#' Returns NULL if model cannot be loaded.
#'
#' @keywords internal
predict_cna_class <- function(scores, model_name) {
  model <- load_prediction_model(model_name)

  if (is.null(model)) {
    warning("Could not load model: ", model_name)
    return(NULL)
  }

  tryCatch({
    predictions <- stats::predict(model, scores, type = "class")
    probabilities <- stats::predict(model, scores, type = "prob")

    result <- data.frame(
      ID = rownames(scores),
      predicted_class = as.character(predictions),
      stringsAsFactors = FALSE
    )

    if (!is.null(probabilities)) {
      result <- cbind(result, probabilities)
    }

    return(result)
  }, error = function(e) {
    warning("Error during prediction: ", e$message)
    return(NULL)
  })
}
