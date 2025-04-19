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
    
    actionButton("updatePlot", "Update Plot"),
    
    
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
  
  
  session$onFlushed(function() {
    # Only trigger once to avoid infinite loops
    isolate({
      updateActionButton(session, "updatePlot", label = NULL, icon = NULL)
      session$sendInputMessage("updatePlot", list(value = 1))
    })
  }, once = TRUE)
  
  
  updateData <- eventReactive(input$updatePlot, {
    if (input$locfactor == "gene") {
      req(input$gene, input$start.range, input$end.range)
      list(
        type = "gene",
        gene = input$gene,
        start.range = input$start.range,
        end.range = input$end.range
      )
    } else {
          req(input$CHR, input$start, input$end)
      list(
        type = "pos",
        chr = input$CHR,
        start = input$start,
        end = input$end
      )
    }
  })


  output$plot <- renderPlotly({
    req(updateData())
    params <- updateData()
    
    if (params$type == "gene") {
      plot.gene(params$gene, params$start.range, params$end.range)
    } else {
      plot.pos(params$chr, params$start, params$end)
    }
  })
  
  
  output$geneControls <- renderUI({
    selectizeInput("gene",
                   label = h5("Gene Symbol"),
                   choices = NULL,   # No choices here
                   selected = "Arc1",
                   multiple = FALSE)
  })
  
  
  observe({
    updateSelectizeInput(
      session,
      inputId = "gene",
      choices = sort(unique(GeneIDKey$Symbol)),
      server = TRUE,
      selected = "Arc1"
    )
  })
  
  
  output$url <- renderUI({
    req(input$gene)
    tags$a(href = paste0("http://flybase.org/reports/", 
                         GeneIDKey[GeneIDKey$Symbol == input$gene, "FBgn"]),
           "FlyBase Gene Summary Page",
           target = "_blank")
  })
}







# Run the application 
shinyApp(ui = ui, server = server)
