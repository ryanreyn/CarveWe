#Main app.R script
ui <- fluidPage(
  titlePanel("CarveWe Dashboard"),
  tabsetPanel(
    tabPanel("Upload", uploadUI("upload")),
    tabPanel("Analysis", analysisUI("analysis")),
    tabPanel("Metabolite Analysis", metabolitesUI("metabolites"))
  )
)

server <- function(input, output, session) {
  
  # Create reactive values
  rv <- reactiveValues(
    user_data = NULL,
    taxonomy = NULL,
    quality = NULL,
    environment = NULL,
    feature_data = NULL,
    has_taxonomy = FALSE,
    has_quality = FALSE,
    has_environment = FALSE,
    upload_complete = FALSE,
    som_complete = FALSE,
    cluster_complete = FALSE
  )
  
  # Run module
  uploadServer("upload", rv)
  analysisServer("analysis", rv)
  metabolitesServer("metabolites", rv)
  
  # Debug: Print when data changes
  observe({
    if (!is.null(rv$user_data)) {
      cat("Data loaded:", nrow(rv$user_data), "rows\n")
    }
  })
}

shinyApp(ui, server)