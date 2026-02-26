#' Run CNApp Shiny Application
#'
#' Launches the interactive CNApp Shiny interface.
#'
#' @param port Integer port number (default 3838).
#' @param ... Additional arguments passed to shiny::runApp().
#'
#' @return No return value; launches Shiny app in browser.
#'
#' @details
#' This function requires the \code{shiny} package to be installed.
#' It will check for Shiny availability and provide installation instructions
#' if necessary.
#'
#' The app is located in inst/shiny/ and includes:
#' - Data upload and validation
#' - Parameter configuration
#' - Re-segmentation and scoring
#' - Visualization and download
#'
#' @examples
#' \dontrun{
#'   run_cnapp_app()
#' }
#'
#' @export
run_cnapp_app <- function(port = 3838, ...) {
  # Check if Shiny is available
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "The 'shiny' package is required to run CNApp interface.\n",
      "Install it with: install.packages('shiny')\n",
      "Or install all optional dependencies with: install.packages(c('shiny', 'shinyBS', 'shinysky', 'shinyjs', 'shinythemes', 'shinyWidgets', 'shinydashboard', 'plotly'))"
    )
  }

  # Locate Shiny app
  app_path <- system.file("shiny", package = "CNApp")

  if (!dir.exists(app_path)) {
    stop("Shiny app not found in package installation: ", app_path)
  }

  # Launch app
  shiny::runApp(app_path, port = port, ...)
}
