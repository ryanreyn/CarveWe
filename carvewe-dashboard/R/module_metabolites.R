# R/module_metabolites.R
# ============================================
# UI Function
# ============================================

metabolitesUI <- function(id){
  ns <- NS(id)
  
  tagList(
    h3("Run metabolite analysis"),
    p("Take processed and clustered data and generate analytical output to visualize metabolic patterns."),
    
    actionButton(ns("run"),
                 "Run metabolite analysis",
                 class = "btn-primary btn-lg",
                 icon = icon("play")),
    
    uiOutput(ns("preview_display")),
    hr(),
    h3("SOM prototype plot:"),
    plotOutput(ns("prototype_plot"), height = "800px"),
    hr(),
    h3("Variance versus Consensus:"),
    plotOutput(ns("variance_plot"), height = "800px"),
    hr(),
    h3("High Sensitivity Plot:"),
    plotOutput(ns("high_sensitivity"), height = "800px"),
    hr(),
    h3("Sensitivity Profiles Plot:"),
    plotOutput(ns("profile_plot"), height = "800px"),
    hr(),
    h3("Flux Bubble Plot:"),
    plotOutput(ns("bubble_plot"), height = "800px")
  )
}


# ============================================
# Server Function
# ============================================

metabolitesServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run, {
      cat("Debugging reactive values:\n")
      cat("  genome_assignments: ", class(rv$genome_assignments), "\n")
      cat("  Dimensions: ", dim(rv$genome_assignments), "\n")
      cat("  Column names: ", paste(colnames(rv$genome_assignments), collapse = ", "), "\n")
      cat("  feature_data rownames present: ", !is.null(rownames(rv$feature_data)), "\n")
    })
    
    scale_data <- eventReactive(input$run, {
      cat("Step XX: Scaling the feature data")
      return_scaled <- prepare_scaled_data(
        rv$feature_data,
        rv$genome_assignments)
      
      return_prototype <- prepare_prototype_data(
        rv$som_codes,
        rv$node_clusters,
        rv$gridsize
      )
      
      list(return_scaled, return_prototype)
    })
    
    observeEvent(scale_data(), {
      req(scale_data())
      
      rv$scaled_data <- scale_data()[[1]]
      rv$prototype_data <- scale_data()[[2]]
    })
    
    output$preview_display <- renderUI({
      req(rv$scaled_data, rv$prototype_data)
      
      tagList(
        h4("Scaled Data Preview"),
        tableOutput(ns("scaled_preview")),
        hr(),
        h4("Prototype Data Preview"),
        tableOutput(ns("prototype_preview"))
      )
    })
    
    output$scaled_preview <- renderTable({
      req(rv$scaled_data)
      head(rv$scaled_data, 5)
    })
    
    output$prototype_preview <- renderTable({
      req(rv$prototype_data)
      pivot_data <- rv$prototype_data$pivot
      head(pivot_data, 5)
    })
    
    analyze_prototypes <- eventReactive(scale_data(), {
      req(scale_data())
      
      result <- analyze_prototype_distributions(
        rv$feature_data,
        rv$genome_assignments)
    })
    
    output$prototype_plot <- renderPlot({
      req(analyze_prototypes())
      
      proto_output <- analyze_prototypes()
      
      plot(proto_output$plot)
    })
    
    analyze_variance <- eventReactive(scale_data(), {
      req(rv$quality)
      
      result <- analyze_replicate_variance(
        rv$feature_data,
        rv$genome_assignments,
        rv$quality)
    })
    
    output$variance_plot <- renderPlot({
      output <- analyze_variance()
      
      plot(output$plot)
    })
    
    analyze_high_sens <- eventReactive(scale_data(), {
      result <- analyze_high_sensitivity_metabolites(
        rv$feature_data,
        rv$genome_assignments
      )
    })
    
    output$high_sensitivity <- renderPlot({
      output <- analyze_high_sens()
      
      plot(output$plot)
    })
    
    analyze_profiles <- eventReactive(scale_data(), {
      result <- analyze_cluster_metabolite_profiles(
        rv$feature_data,
        rv$genome_assignments
      )
    })
    
    output$profile_plot <- renderPlot({
      output <- analyze_profiles()
      
      plot(output$plot)
    })
    
    analyze_bubbles <- eventReactive(scale_data(), {
      result <- analyze_metabolite_bubble_plot(
        rv$feature_data,
        rv$genome_assignments
      )
    })
    
    output$bubble_plot <- renderPlot({
      output <- analyze_bubbles()
      
      plot(output$plot)
    })
    
  })
}