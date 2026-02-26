# CNApp Shiny Server
# Thin wrapper: delegates heavy lifting to R/package core functions

library(shiny)
library(CNApp)
library(DT)

function(input, output, session) {
  useShinyjs()

  # ===== REACTIVE STATE =====
  rv <- reactiveValues(
    data = NULL,
    data_validation = NULL,
    variables_info = NULL,
    sample_names = NULL,
    resegmented_list = NULL,
    scores = NULL,
    status_msg = NULL
  )

  # ===== TAB 1: LOAD DATA =====

  observeEvent(input$btn_load_data, {
    tryCatch({
      # Validate file upload
      req(input$data_file)

      # Read main data
      rv$data <- read_cna_file(
        input$data_file$datapath,
        sep = input$sep,
        dec = input$dec
      )

      # Validate structure
      rv$data_validation <- validate_cna_data(rv$data)

      if (!rv$data_validation$valid) {
        rv$status_msg <- paste("Validation error:", paste(rv$data_validation$errors, collapse = "; "))
        output$load_status <- renderUI({
          div(style = "color: red; padding: 10px; border: 1px solid red; margin: 10px 0;",
              rv$status_msg)
        })
        return()
      }

      # Load annotation if provided
      if (!is.null(input$annot_file)) {
        rv$data <- prepare_annotation_data(
          rv$data,
          input$annot_file$datapath,
          sep = input$sep,
          dec = input$dec
        )
      }

      # Prepare variables
      var_prep <- prepare_clinical_variables(rv$data)
      rv$variables_info <- var_prep$variables_info
      rv$sample_names <- var_prep$sample_names

      rv$status_msg <- paste("Loaded successfully:", nrow(rv$data), "rows,",
                             length(rv$sample_names), "samples")

      output$load_status <- renderUI({
        div(style = "color: green; padding: 10px; border: 1px solid green; margin: 10px 0;",
            rv$status_msg)
      })

      output$data_summary <- renderUI({
        div(
          h4("Data Summary"),
          p("Samples:", length(rv$sample_names)),
          p("Total rows:", nrow(rv$data)),
          p("Columns:", paste(colnames(rv$data)[1:5], collapse = ", "), "...")
        )
      })
    }, error = function(e) {
      rv$status_msg <- paste("Error:", e$message)
      output$load_status <- renderUI({
        div(style = "color: red; padding: 10px;", rv$status_msg)
      })
    })
  })

  observeEvent(input$btn_clear, {
    rv$data <- NULL
    rv$resegmented_list <- NULL
    rv$scores <- NULL
    output$load_status <- renderUI({})
    output$data_summary <- renderUI({})
    output$reseg_results <- renderUI({})
    output$results_panel <- renderUI({})
  })

  # ===== TAB 2: RE-SEGMENTATION =====

  observeEvent(input$btn_run_reseg, {
    tryCatch({
      req(rv$data)

      # Show busy indicator
      shinyjs::show("busy_reseg")

      # Run re-segmentation for all samples
      rv$resegmented_list <- lapply(rv$sample_names, function(sample_id) {
        resegment_sample(
          rv$data,
          sample_id = sample_id,
          min_length = input$min_length,
          max_dist_segm = input$max_dist,
          dev_btw_segs = input$dev_btw,
          dev_tozero = input$dev_zero,
          dev_baf = input$dev_baf,
          low_gain = input$low_gain,
          normal_gain = input$norm_gain,
          high_gain = input$high_gain,
          low_loss = input$low_loss,
          normal_loss = input$norm_loss,
          high_loss = input$high_loss,
          chrom_percent = input$chrom_pct,
          arm_percent = input$arm_pct
        )
      })

      names(rv$resegmented_list) <- rv$sample_names

      # Calculate scores
      rv$scores <- calculate_cna_scores(rv$resegmented_list, rv$sample_names)

      # Hide busy indicator
      shinyjs::hide("busy_reseg")

      # Show results
      output$reseg_results <- renderUI({
        div(
          h4("Re-segmentation Complete"),
          p("Processed", length(rv$sample_names), "samples"),
          dataTableOutput("scores_table"),
          br(),
          downloadButton("dl_reseg_data", "Download Re-segmented Data"),
          HTML("&nbsp;"),
          downloadButton("dl_scores", "Download Scores")
        )
      })

      output$scores_table <- DT::renderDataTable({
        DT::datatable(rv$scores, options = list(pageLength = 10))
      })

      output$dl_reseg_data <- downloadHandler(
        filename = function() paste0("resegmented_", Sys.Date(), ".tsv"),
        content = function(file) {
          # Bind all resegmented data
          all_reseg <- do.call(rbind, rv$resegmented_list)
          utils::write.table(all_reseg, file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      )

      output$dl_scores <- downloadHandler(
        filename = function() paste0("scores_", Sys.Date(), ".tsv"),
        content = function(file) {
          utils::write.table(rv$scores, file, sep = "\t", quote = FALSE)
        }
      )

      output$results_panel <- renderUI({
        div(
          h4("Analysis Results"),
          tabsetPanel(
            tabPanel("Scores Summary", DT::dataTableOutput("scores_table")),
            tabPanel("Download",
                     downloadButton("dl_reseg_data", "Re-segmented Data"),
                     HTML("<br><br>"),
                     downloadButton("dl_scores", "CNA Scores"))
          )
        )
      })
    }, error = function(e) {
      shinyjs::hide("busy_reseg")
      output$reseg_results <- renderUI({
        div(style = "color: red;", p("Error:", e$message))
      })
    })
  })
}
