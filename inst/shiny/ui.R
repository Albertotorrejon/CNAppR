# CNApp Shiny UI
# Minimal interface delegating logic to server.R

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)

dashboardPage(
  dashboardHeader(
    title = "CNApp - Copy Number Analysis",
    titleWidth = 400
  ),

  dashboardSidebar(
    sidebarMenu(
      id = "main_menu",
      menuItem("Load Data", tabName = "tab_load", icon = icon("folder-open"), selected = TRUE),
      menuItem("RE-SEG & SCORE", tabName = "tab_reseg", icon = icon("star")),
      menuItem("Region Profile", tabName = "tab_profile", icon = icon("map")),
      menuItem("Classifier", tabName = "tab_classifier", icon = icon("microchip")),
      menuItem("Results", tabName = "tab_results", icon = icon("bar-chart"))
    )
  ),

  dashboardBody(
    useShinyjs(),

    tabItems(
      # ===== TAB 1: LOAD DATA =====
      tabItem(
        tabName = "tab_load",
        h2("Load Copy Number Data"),
        fluidRow(
          column(6,
                 box(title = "Main Data File", status = "primary", solidHeader = TRUE, width = 12,
                     fileInput("data_file", "Choose CNV data file (TSV/CSV):", accept = c(".tsv", ".csv", ".txt")),
                     selectInput("sep", "Column separator:", choices = list("Tab" = "\t", "Comma" = ","), selected = "\t"),
                     selectInput("dec", "Decimal separator:", choices = list("." = ".", "," = ","), selected = "."),
                     selectInput("genome", "Genome build:", choices = list("GRCh38/hg38" = "hg38", "GRCh37/hg19" = "hg19"), selected = "hg38")
                 )
          ),
          column(6,
                 box(title = "Annotation Data (optional)", status = "info", solidHeader = TRUE, width = 12,
                     fileInput("annot_file", "Choose annotation/clinical data (TSV/CSV):", accept = c(".tsv", ".csv", ".txt")),
                     HTML("<small>Must have ID column matching main data.</small>")
                 )
          )
        ),

        fluidRow(
          column(12,
                 div(style = "text-align: center; padding: 20px;",
                     actionBttn("btn_load_data", "Load & Validate Data", style = "fill", size = "lg", color = "primary"),
                     HTML("&nbsp;&nbsp;"),
                     actionBttn("btn_clear", "Clear", style = "simple", size = "sm", color = "danger")
                 )
          )
        ),

        fluidRow(
          column(12,
                 uiOutput("load_status"),
                 uiOutput("data_summary")
          )
        )
      ),

      # ===== TAB 2: RE-SEG & SCORE =====
      tabItem(
        tabName = "tab_reseg",
        h2("Re-segmentation and CNA Scoring"),

        fluidRow(
          column(6,
                 box(title = "Re-segmentation Parameters", status = "primary", solidHeader = TRUE, width = 12,
                     numericInput("min_length", "Min segment length (bp):", value = 100000, min = 0, step = 10000),
                     numericInput("max_dist", "Max distance between segments (bp):", value = 1000000, min = 0, step = 100000),
                     numericInput("dev_btw", "Max seg.mean deviation:", value = 0.16, min = 0, max = 1, step = 0.01),
                     numericInput("dev_zero", "Min seg.mean deviation from zero:", value = 0.16, min = 0, max = 1, step = 0.01),
                     numericInput("dev_baf", "Max BAF deviation:", value = 0.1, min = 0, max = 1, step = 0.01)
                 )
          ),
          column(6,
                 box(title = "CNA Classification Thresholds", status = "primary", solidHeader = TRUE, width = 12,
                     fluidRow(
                       column(6, h5("Gains"),
                              numericInput("high_gain", "High gain:", value = 1.0, step = 0.01),
                              numericInput("norm_gain", "Normal gain:", value = 0.585, step = 0.01),
                              numericInput("low_gain", "Low gain:", value = 0.2, step = 0.01)
                       ),
                       column(6, h5("Losses"),
                              numericInput("low_loss", "Low loss:", value = -0.2, step = 0.01),
                              numericInput("norm_loss", "Normal loss:", value = -1.0, step = 0.01),
                              numericInput("high_loss", "High loss:", value = -1.737, step = 0.01)
                       )
                     ),
                     sliderInput("chrom_pct", "Chromosomal coverage %:", 0, 1, 0.9, step = 0.05),
                     sliderInput("arm_pct", "Arm-level coverage %:", 0, 1, 0.5, step = 0.05)
                 )
          )
        ),

        fluidRow(
          column(12,
                 div(style = "text-align: center; padding: 20px;",
                     actionBttn("btn_run_reseg", "Run Re-segmentation", style = "fill", size = "lg", color = "primary"),
                     HTML("&nbsp;&nbsp;"),
                     div(id = "busy_reseg", style = "display: none;", HTML("<i class='fas fa-spinner fa-spin'></i> Processing..."))
                 )
          )
        ),

        fluidRow(
          column(12,
                 uiOutput("reseg_results")
          )
        )
      ),

      # ===== TAB 3: REGION PROFILE =====
      tabItem(
        tabName = "tab_profile",
        h2("Region Profile Analysis"),
        fluidRow(
          column(12, "Region profile functionality (future)")
        )
      ),

      # ===== TAB 4: CLASSIFIER =====
      tabItem(
        tabName = "tab_classifier",
        h2("Classifier Models"),
        fluidRow(
          column(12, "Classifier functionality (future)")
        )
      ),

      # ===== TAB 5: RESULTS =====
      tabItem(
        tabName = "tab_results",
        h2("Results and Downloads"),
        fluidRow(
          column(12,
                 uiOutput("results_panel")
          )
        )
      )
    )
  )
)
