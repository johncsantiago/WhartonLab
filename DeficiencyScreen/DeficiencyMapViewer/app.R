##library(rsconnect)
##deployapp()

library(shiny)
library(plotly)

source("https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/DeficiencyScreen/DeficiencyMapViewerFunctions.R")

ui <- fluidPage(
  
  ##Sidebar control code 
  titlePanel("Deficiency Map Viewer"),
  sidebarLayout(position="left", sidebarPanel(
    
    ##uiOutput("click"),
    
    fluidRow(
      ##Center on a gene of interest or a position
      column(6,
             radioButtons("locfactor", h5("Select Location Type"), 
                         choices = list("Gene of interest" = "gene",
                                        "Location"     = "pos"),
                         selected = "gene"),),),
    
    conditionalPanel(
      condition = "input.locfactor == 'gene'",
      ##gene.Controls; Gene selection tool. Uses a reactive variable
      uiOutput("geneControls"),
      p("Select Gene: Use gene symbol if available. autocompletes"),
      
      ##url; Link to selected gene flybase summary page. Uses a reactive variable
      uiOutput("url"),
      br(),
      
      fluidRow(
        ##Upstream from gene start
        column(6,
               numericInput("start.range",
                            h5("Upstream Distance from TSS"),
                            min = 0,
                            value = 1000000),),
        
        ##Downstream from gene end
        column(6,
               numericInput("end.range",
                            h5("Downstream Distance from TES"),
                            min = 0,
                            value = 1000000),),
      ),
      ),
    
    conditionalPanel(
      condition = "input.locfactor == 'pos'",
      selectInput("CHR", h5("Chromosome"),
                  choices = list("X"    = "X",
                                 "2L"    = "2L",
                                 "2R"    = "2R",
                                 "3L"    = "3L",
                                 "3R"     = "3R"),
                  selected="X"),
      fluidRow(
        ##Upstream from gene start
        column(6,
               numericInput("start",
                            h5("Start Position"),
                            min = 0,
                            value = 14000000),),
        
        ##Downstream from gene end
        column(6,
               numericInput("end",
                            h5("End Position"),
                            min = 0,
                            value = 15000000),),
      ),),
    
    width = 3),
    
    ##Figure display settings    
    mainPanel(
      h5("Click gene for flybase link"),
      uiOutput("click"),
      ##Display the boxplot
      plotlyOutput(outputId = "plot",
                   height=800, width=1250)
    ),),)



# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  output$geneControls = renderUI({
    selectizeInput("gene", label = "Select Gene", 
                   choices = NULL,  # leave empty for server-side
                   options = list(
                     create = TRUE,
                     placeholder = 'Select a gene name'
                   ),
                   selected = "Arc1")
  })  
  
  observe({
    updateSelectizeInput(
      session, 
      "gene", 
      choices = GeneIDKey$Symbol, 
      server = TRUE
    )
  })
  
  # Debounced reactive inputs
  debounced_gene <- debounce(reactive(input$gene), 500)
  debounced_start <- debounce(reactive(input$start.range), 500)
  debounced_end <- debounce(reactive(input$end.range), 500)
  
  
  output$plot <- renderPlotly({
    if (input$locfactor == "gene") {
      req(debounced_gene())  # Wait until gene input stabilizes
      
      p = plot.gene(
        gene.symbol = debounced_gene(),
        start.range = debounced_start(),
        end.range = debounced_end()
      )
    } else {
      p = plot.pos(
        chr = input$CHR,
        start = input$start,
        end = input$end
      )
    }
    
    p %>%
      layout(dragmode = "select") %>%
      event_register("plotly_selecting") %>%
      event_register("plotly_click")
  })
  
  
  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) {
      "Flybase link will appear here after click (double-click to clear)"
    } else {
      page = a(
        paste0("Flybase Page for ", fbID$Symbol[grep(d$x, fbID$Position)]),
        href = paste0("http://flybase.org/reports/", fbID$FBgn[grep(d$x, fbID$Position)]),
        target = "_blank"
      )
      tagList(page)
    }
  })
}






# Run the application 
shinyApp(ui = ui, server = server)
