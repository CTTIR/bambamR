# bambamR Shiny Application
# =========================

# Check for required packages
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required. Install with: install.packages('shiny')")
}

# Use bslib if available, otherwise shinydashboard
use_bslib <- requireNamespace("bslib", quietly = TRUE)
use_dt <- requireNamespace("DT", quietly = TRUE)

# --- UI ---
if (use_bslib) {
  ui <- bslib::page_navbar(
    title = "bambamR",
    theme = bslib::bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary = "#2166AC"
    ),

    # Home tab
    bslib::nav_panel("Home",
      bslib::layout_column_wrap(
        width = 1,
        bslib::card(
          bslib::card_header("Welcome to bambamR"),
          bslib::card_body(
            shiny::h4("End-to-End RNA-Seq Processing"),
            shiny::p("bambamR provides a streamlined pipeline from FASTQ/BAM
                      files through alignment, counting, differential expression,
                      and publication-ready visualizations."),
            shiny::hr(),
            shiny::h5("Quick Start"),
            shiny::tags$ol(
              shiny::tags$li("Upload your count matrix or BAM files in the Data Input tab"),
              shiny::tags$li("Configure pipeline parameters"),
              shiny::tags$li("Run the analysis"),
              shiny::tags$li("Explore interactive results and export")
            )
          )
        )
      )
    ),

    # Data Input tab
    bslib::nav_panel("Data Input",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          shiny::h4("Upload Data"),
          shiny::fileInput("count_file", "Count Matrix (CSV/TSV/RDS)",
                           accept = c(".csv", ".tsv", ".rds")),
          shiny::fileInput("metadata_file", "Sample Metadata (CSV)",
                           accept = ".csv"),
          shiny::hr(),
          shiny::actionButton("load_example", "Load Example Data",
                               class = "btn-primary")
        ),
        bslib::card(
          bslib::card_header("Data Preview"),
          bslib::card_body(
            shiny::verbatimTextOutput("data_summary"),
            if (use_dt) DT::DTOutput("count_preview")
            else shiny::tableOutput("count_preview_basic")
          )
        )
      )
    ),

    # Pipeline Config tab
    bslib::nav_panel("Pipeline",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          shiny::h4("Configuration"),
          shiny::selectInput("de_method", "DE Method",
                              choices = c("DESeq2", "edgeR", "limma")),
          shiny::selectInput("norm_method", "Normalization",
                              choices = c("cpm", "tpm", "tmm", "rle")),
          shiny::numericInput("p_cutoff", "P-value Cutoff", 0.05,
                               min = 0, max = 1, step = 0.01),
          shiny::numericInput("fc_cutoff", "Log2FC Cutoff", 1,
                               min = 0, step = 0.5),
          shiny::hr(),
          shiny::actionButton("run_pipeline", "Run Analysis",
                               class = "btn-success btn-lg")
        ),
        bslib::card(
          bslib::card_header("Analysis Log"),
          bslib::card_body(
            shiny::verbatimTextOutput("pipeline_log")
          )
        )
      )
    ),

    # Results tab
    bslib::nav_panel("Results",
      bslib::navset_card_tab(
        bslib::nav_panel("Volcano",
          shiny::plotOutput("volcano_plot", height = "600px")
        ),
        bslib::nav_panel("PCA",
          shiny::plotOutput("pca_plot", height = "600px")
        ),
        bslib::nav_panel("Heatmap",
          shiny::plotOutput("heatmap_plot", height = "700px")
        ),
        bslib::nav_panel("MA Plot",
          shiny::plotOutput("ma_plot", height = "600px")
        ),
        bslib::nav_panel("DE Table",
          if (use_dt) DT::DTOutput("de_table")
          else shiny::tableOutput("de_table_basic")
        )
      )
    ),

    # Export tab
    bslib::nav_panel("Export",
      bslib::card(
        bslib::card_header("Export Results"),
        bslib::card_body(
          shiny::downloadButton("dl_de_csv", "Download DE Results (CSV)"),
          shiny::downloadButton("dl_result_rds", "Download Full Result (RDS)"),
          shiny::downloadButton("dl_volcano_pdf", "Download Volcano (PDF)"),
          shiny::downloadButton("dl_heatmap_pdf", "Download Heatmap (PDF)")
        )
      )
    )
  )
} else {
  # Fallback: basic shiny UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("bambamR - RNA-Seq Analysis"),
    shiny::tabsetPanel(
      shiny::tabPanel("Home",
        shiny::h3("Welcome to bambamR"),
        shiny::p("Upload your data in the Data Input tab to get started.")
      ),
      shiny::tabPanel("Data Input",
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::fileInput("count_file", "Count Matrix (CSV/TSV/RDS)"),
            shiny::fileInput("metadata_file", "Sample Metadata (CSV)"),
            shiny::actionButton("load_example", "Load Example Data")
          ),
          shiny::mainPanel(
            shiny::verbatimTextOutput("data_summary")
          )
        )
      ),
      shiny::tabPanel("Pipeline",
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::selectInput("de_method", "DE Method",
                                c("DESeq2", "edgeR", "limma")),
            shiny::numericInput("p_cutoff", "P-value Cutoff", 0.05),
            shiny::numericInput("fc_cutoff", "Log2FC Cutoff", 1),
            shiny::actionButton("run_pipeline", "Run Analysis")
          ),
          shiny::mainPanel(
            shiny::verbatimTextOutput("pipeline_log")
          )
        )
      ),
      shiny::tabPanel("Results",
        shiny::plotOutput("volcano_plot", height = "500px"),
        shiny::plotOutput("pca_plot", height = "500px")
      ),
      shiny::tabPanel("Export",
        shiny::downloadButton("dl_de_csv", "Download DE Results (CSV)")
      )
    )
  )
}

# --- Server ---
server <- function(input, output, session) {

  # Reactive values
  rv <- shiny::reactiveValues(
    counts = NULL,
    metadata = NULL,
    de_results = NULL,
    norm_counts = NULL,
    result = NULL,
    log_text = ""
  )

  # Helper to append log
  log_msg <- function(msg) {
    rv$log_text <- paste0(rv$log_text, "\n", format(Sys.time(), "%H:%M:%S"),
                          " | ", msg)
  }

  # Load example data
  shiny::observeEvent(input$load_example, {
    example_path <- system.file("extdata", "example_counts.rds",
                                 package = "bambamR")
    if (file.exists(example_path)) {
      example_data <- readRDS(example_path)
      rv$counts <- example_data$counts
      rv$metadata <- example_data$metadata
      log_msg("Loaded example data")
    } else {
      log_msg("Example data not found. Generating synthetic data...")
      set.seed(42)
      rv$counts <- matrix(
        rpois(1000, lambda = 100),
        nrow = 100, ncol = 10,
        dimnames = list(paste0("gene", 1:100), paste0("S", 1:10))
      )
      rv$metadata <- data.frame(
        condition = factor(rep(c("ctrl", "treat"), each = 5)),
        row.names = paste0("S", 1:10)
      )
      log_msg("Generated synthetic example data")
    }
  })

  # Upload count matrix
  shiny::observeEvent(input$count_file, {
    req <- shiny::req
    req(input$count_file)
    ext <- tools::file_ext(input$count_file$name)
    tryCatch({
      if (ext == "rds") {
        rv$counts <- readRDS(input$count_file$datapath)
      } else if (ext == "csv") {
        df <- data.table::fread(input$count_file$datapath, data.table = FALSE)
        rownames(df) <- df[[1]]
        rv$counts <- as.matrix(df[, -1])
      } else if (ext == "tsv") {
        df <- data.table::fread(input$count_file$datapath, sep = "\t",
                                 data.table = FALSE)
        rownames(df) <- df[[1]]
        rv$counts <- as.matrix(df[, -1])
      }
      log_msg(paste("Loaded count matrix:", nrow(rv$counts), "genes x",
                     ncol(rv$counts), "samples"))
    }, error = function(e) {
      log_msg(paste("Error loading counts:", e$message))
    })
  })

  # Upload metadata
  shiny::observeEvent(input$metadata_file, {
    req <- shiny::req
    req(input$metadata_file)
    tryCatch({
      df <- data.table::fread(input$metadata_file$datapath, data.table = FALSE)
      rownames(df) <- df[[1]]
      rv$metadata <- df[, -1, drop = FALSE]
      log_msg(paste("Loaded metadata:", nrow(rv$metadata), "samples"))
    }, error = function(e) {
      log_msg(paste("Error loading metadata:", e$message))
    })
  })

  # Data summary
  output$data_summary <- shiny::renderPrint({
    if (is.null(rv$counts)) {
      cat("No data loaded. Upload a count matrix or click 'Load Example Data'.")
    } else {
      cat("Count matrix:", nrow(rv$counts), "genes x", ncol(rv$counts),
          "samples\n")
      if (!is.null(rv$metadata)) {
        cat("Metadata columns:", paste(colnames(rv$metadata), collapse = ", "),
            "\n")
      }
    }
  })

  # Pipeline log
  output$pipeline_log <- shiny::renderText(rv$log_text)

  # Run pipeline
  shiny::observeEvent(input$run_pipeline, {
    req <- shiny::req
    req(rv$counts)

    log_msg("Starting analysis...")

    tryCatch({
      # Normalize
      rv$norm_counts <- bambamR::bb_normalize(rv$counts, method = "cpm")
      log_msg("Normalization complete (CPM)")

      # DE analysis
      if (!is.null(rv$metadata)) {
        rv$de_results <- tryCatch(
          .run_de_app(rv$counts, rv$metadata, input$de_method),
          error = function(e) {
            log_msg(paste("DE failed:", e$message))
            NULL
          }
        )
        if (!is.null(rv$de_results)) {
          n_sig <- sum(rv$de_results$padj < input$p_cutoff, na.rm = TRUE)
          log_msg(paste("DE complete:", n_sig, "significant genes"))
        }
      }

      log_msg("Analysis complete!")
    }, error = function(e) {
      log_msg(paste("Pipeline error:", e$message))
    })
  })

  # Plots
  output$volcano_plot <- shiny::renderPlot({
    req <- shiny::req
    req(rv$de_results)
    bambamR::bb_volcano(rv$de_results, fc_cutoff = input$fc_cutoff,
                         p_cutoff = input$p_cutoff)
  })

  output$pca_plot <- shiny::renderPlot({
    req <- shiny::req
    req(rv$norm_counts, rv$metadata)
    if ("condition" %in% colnames(rv$metadata)) {
      bambamR::bb_pca(rv$norm_counts, rv$metadata, color_by = "condition")
    }
  })

  output$heatmap_plot <- shiny::renderPlot({
    req <- shiny::req
    req(rv$norm_counts)
    bambamR::bb_heatmap(rv$norm_counts, de_result = rv$de_results, n_genes = 50)
  })

  output$ma_plot <- shiny::renderPlot({
    req <- shiny::req
    req(rv$de_results)
    if ("basemean" %in% colnames(rv$de_results)) {
      bambamR::bb_ma_plot(rv$de_results, p_cutoff = input$p_cutoff)
    }
  })

  # DE table
  if (use_dt) {
    output$de_table <- DT::renderDT({
      req <- shiny::req
      req(rv$de_results)
      DT::datatable(rv$de_results, options = list(pageLength = 25))
    })
  }

  # Downloads
  output$dl_de_csv <- shiny::downloadHandler(
    filename = function() paste0("bambamR_DE_", Sys.Date(), ".csv"),
    content = function(file) {
      if (!is.null(rv$de_results)) {
        data.table::fwrite(rv$de_results, file, sep = ",")
      }
    }
  )

  output$dl_result_rds <- shiny::downloadHandler(
    filename = function() paste0("bambamR_result_", Sys.Date(), ".rds"),
    content = function(file) {
      saveRDS(list(
        counts = rv$counts,
        metadata = rv$metadata,
        de_results = rv$de_results,
        norm_counts = rv$norm_counts
      ), file)
    }
  )

  output$dl_volcano_pdf <- shiny::downloadHandler(
    filename = function() "volcano_plot.pdf",
    content = function(file) {
      if (!is.null(rv$de_results)) {
        p <- bambamR::bb_volcano(rv$de_results)
        ggplot2::ggsave(file, p, width = 8, height = 6)
      }
    }
  )

  output$dl_heatmap_pdf <- shiny::downloadHandler(
    filename = function() "heatmap_plot.pdf",
    content = function(file) {
      if (!is.null(rv$norm_counts)) {
        p <- bambamR::bb_heatmap(rv$norm_counts, de_result = rv$de_results)
        ggplot2::ggsave(file, p, width = 10, height = 8)
      }
    }
  )
}

# Helper for DE in app context
.run_de_app <- function(counts, metadata, method) {
  shared <- intersect(colnames(counts), rownames(metadata))
  if (length(shared) < 2) return(NULL)
  counts <- counts[, shared, drop = FALSE]
  metadata <- metadata[shared, , drop = FALSE]

  switch(method,
    DESeq2 = {
      if (!requireNamespace("DESeq2", quietly = TRUE)) return(NULL)
      bambamR::bb_deseq2(counts, metadata)
    },
    edgeR = {
      if (!requireNamespace("edgeR", quietly = TRUE)) return(NULL)
      bambamR::bb_edger(counts, metadata$condition)
    },
    limma = {
      if (!requireNamespace("limma", quietly = TRUE)) return(NULL)
      if (!requireNamespace("edgeR", quietly = TRUE)) return(NULL)
      dm <- stats::model.matrix(~ condition, data = metadata)
      bambamR::bb_limma_voom(counts, dm)
    }
  )
}

shiny::shinyApp(ui = ui, server = server)
