# R/module_upload.R
# ============================================
# UI Function
# ============================================
uploadUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Data Selection"),
    p("Choose whether you would like to use the demo data from the publication associated with this dashboard or upload your own data."),
    
    radioButtons(
      ns("data_source"),
      "Choose your data source:",
      choices = c("Use demo data" = "demo",
                  "Upload my own data" = "upload"),
      selected = "upload"
    ),
    
    # ========================================
    # UPLOAD PATH
    # ========================================
    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
      
      h3("Upload Data"),
      p("Upload your CarveMe analysis outputs. Sensitivities file is required; other files are optional but enable additional analyses."),
      
      wellPanel(
        h4("Required Data"),
        fileInput(
          ns("file_sens"),
          "Sensitivities CSV",
          accept = ".csv",
          buttonLabel = "Browse...",
          placeholder = "No file selected"
        ),
        helpText("Format: columns must include genome_id, nutrient_class, sensitivity_score")
      ),
      
      wellPanel(
        h4("Optional Data (for contextualization)"),
        
        fileInput(ns("file_tax"), "Taxonomy CSV (optional)", accept = ".csv"),
        helpText("Format: genome_id, phylum, class, genus, species"),
        
        fileInput(ns("file_qual"), "Quality Metrics CSV (optional)", accept = ".csv"),
        helpText("Format: genome_id, completeness, contamination"),
        
        fileInput(ns("file_env"), "Environmental Data CSV (optional)", accept = ".csv"),
        helpText("Format: genome_id, temperature_category, depth_category, etc.")
      ),
      
      actionButton(
        ns("process"),
        "Load & Validate Data",
        class = "btn-primary btn-lg",
        icon = icon("upload")
      )
    ),
    
    # ========================================
    # DEMO PATH
    # ========================================
    conditionalPanel(
      condition = sprintf("input['%s'] == 'demo'", ns("data_source")),
      
      h3("Demo Data"),
      p("Load example dataset from the publication to explore the dashboard features."),
      
      actionButton(
        ns("use_demo"),
        "Load Demo Data",
        class = "btn-primary btn-lg",
        icon = icon("database")
      )
    ),
    
    # ========================================
    # SHARED PREVIEWS (shown for both paths)
    # ========================================
    hr(),
    h4("Data Preview"),
    uiOutput(ns("load_status")),
    tableOutput(ns("load_preview")),
    uiOutput(ns("transform_status")),
    tableOutput(ns("transform_preview"))
  )
}

# ============================================
# Server Function
# ============================================
uploadServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    # ========================================
    # VALIDATION FUNCTION
    # ========================================
    
    validate_sensitivities <- function(data) {
      required_cols <- c("genome_id", "nutrient_class", "sensitivity_score")
      
      missing_cols <- setdiff(required_cols, names(data))
      
      if (length(missing_cols) > 0) {
        stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
      }
      
      if (any(is.na(data$genome_id))) {
        stop("genome_id column contains missing values")
      }
      
      if (any(is.na(data$sensitivity_score))) {
        stop("sensitivity_score column contains missing values")
      }
      
      if (!is.numeric(data$sensitivity_score)) {
        stop("sensitivity_score must be numeric")
      }
      
      return(TRUE)
    }
    
    # ========================================
    # UPLOAD PATH: File preview & processing
    # ========================================
    
    # Quick preview of uploaded file (first 5 rows)
    upload_preview <- reactive({
      req(input$file_sens)
      
      tryCatch({
        read_csv(input$file_sens$datapath, n_max = 5, show_col_types = FALSE)
      }, error = function(e) {
        showNotification(paste("Error reading file:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # Full load when button clicked
    uploaded_data <- eventReactive(input$process, {
      req(input$file_sens)
      
      withProgress(message = "Loading data files...", value = 0, {
        
        incProgress(0.2, detail = "Reading sensitivities...")
        sens <- read_csv(input$file_sens$datapath, show_col_types = FALSE)
        
        incProgress(0.4, detail = "Validating format...")
        validate_sensitivities(sens)
        
        incProgress(0.5, detail = "Transforming to feature matrix...")
        feat <- generate_feature_data(sens)
        
        # Load optional files
        tax <- NULL
        if (!is.null(input$file_tax)) {
          incProgress(0.6, detail = "Reading taxonomy...")
          tax <- read_csv(input$file_tax$datapath, show_col_types = FALSE)
        }
        
        qual <- NULL
        if (!is.null(input$file_qual)) {
          incProgress(0.7, detail = "Reading quality metrics...")
          qual <- read_csv(input$file_qual$datapath, show_col_types = FALSE)
        }
        
        env <- NULL
        if (!is.null(input$file_env)) {
          incProgress(0.8, detail = "Reading environmental data...")
          env <- read_csv(input$file_env$datapath, show_col_types = FALSE)
        }
        
        incProgress(1.0, detail = "Complete!")
        
        list(
          sensitivities = sens,
          features = feat,
          taxonomy = tax,
          quality = qual,
          environment = env
        )
      })
    })
    
    # Store uploaded data to rv
    observeEvent(uploaded_data(), {
      data <- uploaded_data()
      
      rv$user_data <- data$sensitivities
      rv$taxonomy <- data$taxonomy
      rv$quality <- data$quality
      rv$environment <- data$environment
      rv$feature_data <- data$features
      
      rv$has_taxonomy <- !is.null(data$taxonomy)
      rv$has_quality <- !is.null(data$quality)
      rv$has_environment <- !is.null(data$environment)
      rv$upload_complete <- TRUE
      
      showNotification(
        paste("Data loaded successfully!", 
              n_distinct(rv$user_data$genome_id), "genomes found."),
        type = "message",
        duration = 5
      )
    })
    
    # ========================================
    # DEMO PATH: Load bundled data
    # ========================================
    
    demo_data <- eventReactive(input$use_demo, {
      
      withProgress(message = "Loading demo data...", value = 0, {
        
        tryCatch({
          incProgress(0.3, detail = "Reading sensitivities...")
          sens <- read_csv("data/processed_sensitivities_long.csv", show_col_types = FALSE)
          
          incProgress(0.6, detail = "Reading quality metrics...")
          qual <- read_csv("data/processed_quality_data.csv", show_col_types = FALSE)
          
          incProgress(0.8, detail = "Transforming to feature matrix...")
          feat <- generate_feature_data(sens)
          
          incProgress(1.0, detail = "Complete!")
          
          list(
            sensitivities = sens,
            features = feat,
            #taxonomy = tax,
            quality = qual
            #environment = env
          )
          
        }, error = function(e) {
          showNotification(paste("Error loading demo data:", e$message), 
                           type = "error", duration = NULL)
          return(NULL)
        })
      })
    })
    
    # Store demo data to rv
    observeEvent(demo_data(), {
      req(demo_data())
      data <- demo_data()
      
      rv$user_data <- data$sensitivities
      #rv$taxonomy <- data$taxonomy
      rv$quality <- data$quality
      #rv$environment <- data$environment
      rv$feature_data <- data$features
      
      #rv$has_taxonomy <- !is.null(data$taxonomy)
      rv$has_quality <- !is.null(data$quality)
      #rv$has_environment <- !is.null(data$environment)
      rv$upload_complete <- TRUE
      
      showNotification(
        paste("Demo data loaded successfully!", 
              n_distinct(rv$user_data$genome_id), "genomes found."),
        type = "message",
        duration = 5
      )
    })
    
    # ========================================
    # SHARED OUTPUTS: Work for both paths
    # ========================================
    
    # Raw data preview (shows upload preview OR demo preview)
    output$load_preview <- renderTable({
      # Show upload preview if in upload mode and file selected
      if (input$data_source == "upload" && !is.null(input$file_sens)) {
        df <- upload_preview()
        if (!is.null(df)) return(head(df, 5))
      }
      
      # Show demo preview if loaded
      if (input$data_source == "demo" && !is.null(demo_data())) {
        return(head(demo_data()$sensitivities, 5))
      }
      
      # Show uploaded data if processed
      if (!is.null(uploaded_data())) {
        return(head(uploaded_data()$sensitivities, 5))
      }
      
      return(NULL)
    })
    
    # Transformed feature data preview
    output$transform_preview <- renderTable({
      req(rv$feature_data)
      head(rv$feature_data, 5)
    })
    
    # Load status message
    output$load_status <- renderUI({
      
      # Check if any data loaded
      data <- if (!is.null(uploaded_data())) {
        uploaded_data()
      } else if (!is.null(demo_data())) {
        demo_data()
      } else {
        return(p("Upload files or load demo data to begin.", class = "text-muted"))
      }
      
      # Success message with summary
      tagList(
        div(
          class = "alert alert-success",
          icon("check-circle"),
          strong("Data loaded successfully!"),
          tags$ul(
            tags$li(sprintf("Sensitivities: %d genomes, %d measurements",
                            n_distinct(data$sensitivities$genome_id),
                            nrow(data$sensitivities))),
            if (!is.null(data$taxonomy)) {
              tags$li(sprintf("Taxonomy: %d genomes", nrow(data$taxonomy)))
            },
            if (!is.null(data$quality)) {
              tags$li(sprintf("Quality metrics: %d genomes", nrow(data$quality)))
            },
            if (!is.null(data$environment)) {
              tags$li(sprintf("Environmental data: %d genomes", nrow(data$environment)))
            }
          )
        )
      )
    })
    
    # Transform status message
    output$transform_status <- renderUI({
      if (is.null(rv$feature_data)) {
        return(p("Awaiting data transformation...", class = "text-muted"))
      }
      
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(sprintf("Feature matrix created: %d genomes × %d metabolites",
                       nrow(rv$feature_data),
                       ncol(rv$feature_data)))
      )
    })
    
  })
}