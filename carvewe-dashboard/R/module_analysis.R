analysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Run SOM Analysis"),
    p("Configure parameters below and click Run to train the SOM."),
    
    numericInput(
      ns("gridsize"), 
      "SOM grid size (square toroid with dimension X by X):", 
      min = 5, max = 50, value = 20
    ),
    
    actionButton(
      ns("run"), 
      "Run SOM Analysis", 
      class = "btn-primary btn-lg",
      icon = icon("play")
    ),
    uiOutput(ns("analysis_status")),
    plotOutput(ns("eval_plot")),
    uiOutput(ns("clustering_info")),
    uiOutput(ns("choose_clusters")),
    hr(),
    plotOutput(ns("grid_plot")),
    uiOutput(ns("results_summary")),
    hr(),
    tableOutput(ns("replicates")),
    tableOutput(ns("genome_clustering"))
  )
}

analysisServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    cat("Analysis module loaded\n")
    
    # Track button clicks
    observeEvent(input$run, {
      cat("=== RUN BUTTON CLICKED ===\n")
      cat("  Grid size:", input$gridsize, "\n")
      cat("  Clusters:", input$clusters, "\n")
    })
    
    # Train SOM (expensive)
    som_result <- eventReactive(input$run, {
      tryCatch({
        cat("Step 1: Starting SOM training\n")
        
        req(rv$feature_data)
        cat("Step 2: Feature data validated\n")
        cat("  - Class:", class(rv$feature_data), "\n")
        cat("  - Dimensions:", dim(rv$feature_data), "\n")
        
        cat("Step 3: Calling train_som()\n")
        result <- train_som(rv$feature_data, grid_size = input$gridsize)
        cat("Step 4: train_som() completed\n")
        
        # Debug structure
        cat("Step 5: Result structure:\n")
        cat("  - Class:", class(result), "\n")
        cat("  - Names:", paste(names(result), collapse = ", "), "\n")
        if ("som_codes" %in% names(result)) {
          cat("  - som_codes dimensions:", dim(result$som_codes), "\n")
        }
        
        return(result)
        
      }, error = function(e) {
        cat("ERROR in train_som:\n")
        cat("  Message:", e$message, "\n")
        showNotification(
          paste("SOM training error:", e$message), 
          type = "error", 
          duration = NULL
        )
        return(NULL)
      })
    })
    
    observeEvent(som_result(), {
      req(som_result())
      
      rv$gridsize <- input$gridsize
      rv$som_result <- som_result()
      rv$som_codes <- som_result()$som_codes
      rv$som_node_assignments <- som_result()$som_model$unit.classif
      rv$som_complete <- TRUE
    })
    
    cluster_testing <- eventReactive(som_result(), {
      req(som_result())
      codes <- rv$som_codes
      
      cluster_eval <- evaluate_cluster_range(codes)
    })
    
    # Output a plot
    output$eval_plot <- renderPlot({
      req(cluster_testing())
      cat("Step X: Rendering plot of intra-cluster distances.\n")
      
      result <- plot_cluster_evaluation(cluster_testing())
      plot(result)
    })
    
    # Render additional UI conditional on SOM finishing and displaying graphics
    output$choose_clusters <- renderUI({
      req(rv$som_complete)
      if (rv$som_complete){
        tagList(
          numericInput(
            ns("clusters"),
            "Number of clusters for segmentation:",
            min = 2, max = 20, value = 8
          ),
          actionButton(
            ns("cluster"),
            "Run Clustering", 
            class = "btn-primary btn-lg",
            icon = icon("play"))
        )
      }
    })
    
    # Cluster SOM nodes (fast - re-runs when k changes)
    cluster_result <- eventReactive(input$cluster,  # Trigger on either change
      {
        req(input$clusters)
        tryCatch({
          cat("Step 6: Starting clustering\n")
          
          req(som_result())
          
          # Extract codes from the reactive VALUE
          codes <- som_result()$som_codes
          cat("Step 7: Extracted codes\n")
          cat("  - Class:", class(codes), "\n")
          cat("  - Dimensions:", dim(codes), "\n")
          
          cat("Step 8: Calling cluster_som_nodes() with k =", input$clusters, "\n")
          result <- cluster_som_nodes(codes, k = input$clusters)
          cat("Step 9: Clustering complete\n")
          cat("  - Result length:", length(result), "\n")
          
          return(result)
          
        }, error = function(e) {
          cat("ERROR in clustering:\n")
          cat("  Message:", e$message, "\n")
          cat("  Call:", deparse(e$call), "\n")
          
          # Debug what we're passing
          cat("  Debug info:\n")
          if (exists("codes")) {
            cat("    - codes class:", class(codes), "\n")
            cat("    - codes is.matrix:", is.matrix(codes), "\n")
            if (is.matrix(codes)) {
              cat("    - codes dim:", dim(codes), "\n")
            }
          }
          
          showNotification(
            paste("Clustering error:", e$message), 
            type = "error", 
            duration = NULL
          )
          return(NULL)
        })
      },
      ignoreNULL = TRUE
    )
    
    # Store results when clustering completes
    observeEvent(cluster_result(), {
      req(cluster_result())
      
      rv$node_clusters <- cluster_result()
      rv$cluster_complete <- TRUE
      
      cat("Step 10: Results stored in rv\n")
      
      showNotification(
        "Analysis complete!", 
        type = "message",
        duration = 3
      )
    })
    
    # Render plot
    output$grid_plot <- renderPlot({
      req(rv$node_clusters)
      cat("Step 11: Rendering plot\n")
      
      cat("  - som_result available:", !is.null(som_result()), "\n")
      cat("  - cluster_result available:", !is.null(cluster_result()), "\n")
      
      # Pass to plotting function (adjust based on what map plotting expects)
      generate_som_map(node_clusters = rv$node_clusters, grid_size = input$gridsize)
    })
    
    check_replicates <- eventReactive(cluster_result(), {
      req(rv$cluster_complete == TRUE)
      
      result <- check_duplicate_replicates(rv$feature_data)
    })
    
    observeEvent(check_replicates(), {
      rv$replicates <- check_replicates()
    })
    
    output$replicates <- renderUI({
      req(check_replicates())
      
      tagList(
        h3("Number of genomes with identical replicate models:"),
        renderTable({
        df <- data.frame(`Number of Duplicates` = rv$replicates[[2]], 
                         `Percent Duplicates` = rv$replicates[[3]])
        head(df, 5)
      })
      )
    })
    
    cluster_genomes <- eventReactive(cluster_result(), {
      req(cluster_result())
      
      result <- assign_genomes_by_plurality(
        feature_matrix = rv$feature_data,
        som_node_assignments = rv$som_node_assignments,
        node_clusters = rv$node_clusters
      )
    })
    
    observeEvent(cluster_genomes(), {
      req(cluster_genomes())
      
      rv$genome_assignments <- cluster_genomes()%>%
        select(c(genome_id, majority_cluster))
    })
    
    output$genome_clustering <- renderUI({
      req(cluster_genomes())
      
      tagList(
        h3("Preview of genome cluster assignments. This is determined by a plurality-based votive threshold if genomes have replicate models."),
        renderTable({
        df <- cluster_genomes()
        head(df, 5)
        })
      )
    })
    
    # Status badges after SOM completes
    output$analysis_status <- renderUI({
      req(rv$som_complete)
      
      tagList(
        h5("SOM Training Complete"),
        span(class = "badge bg-success", 
             icon("check-circle"), "SOM Trained"),
        span(class = "badge bg-primary", 
             icon("th"), sprintf("%d × %d grid (%d nodes)", 
                                 input$gridsize, input$gridsize, 
                                 input$gridsize^2)),
        span(class = "badge bg-info",
             icon("database"), sprintf("%d genomes mapped", 
                                       length(rv$som_node_assignments)))
      )
    })
    
    # Clustering quality info
    output$clustering_info <- renderUI({
      req(cluster_testing())
      
      eval_data <- cluster_testing()
      optimal_k <- eval_data$k[which.min(eval_data$kmeans_dist)]
      
      div(class = "alert alert-info",
          icon("info-circle"),
          strong("Clustering Evaluation Complete"),
          tags$p(sprintf("Tested k = %d to %d clusters", 
                         min(eval_data$k), max(eval_data$k))),
          tags$p(sprintf("Suggested optimal k based on elbow: %d", optimal_k),
                 style = "margin-bottom: 0;")
      )
    })
    
    # Results summary after clustering
    output$results_summary <- renderUI({
      req(rv$cluster_complete)
      
      cluster_sizes <- table(rv$node_clusters)
      
      tagList(
        div(class = "alert alert-success",
            icon("check-circle"),
            strong("Clustering Complete!"),
            tags$ul(
              tags$li(sprintf("Method: k-means with k = %d", input$clusters)),
              tags$li(sprintf("Cluster sizes range: %d - %d nodes",
                              min(cluster_sizes), max(cluster_sizes))),
              tags$li(sprintf("Genome assignments: %d total", 
                              nrow(rv$genome_assignments)))
            )
        ),
        
        # Compact cluster size badges
        div(style = "margin-top: 10px;",
            h6("Cluster Sizes:"),
            lapply(seq_along(cluster_sizes), function(i) {
              span(class = "badge bg-secondary", 
                   style = "margin-right: 5px;",
                   sprintf("C%d: %d", i, cluster_sizes[i]))
            })
        )
      )
    })

  })
}