library(shiny)
library(gplots)
library(DT)
library(gplots)

source('https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/WhartonLabOmicsDataApp/WhartonLabOmicsDataAppFunctions.R')

# UI
ui <- fluidPage(
  titlePanel("Wharton Lab Omics Data Viewer"),
  
  # Move the input above the main panel
  fluidRow(
    column(
      width = 12,
      selectInput(
        "datasetType",
        "Select Data Type:",
        choices = c("G85R-TKT Metabolomics", "G85R-TKT Transcriptomics", "A4V Transcriptomics"),
        selected = "G85R-TKT Metabolomics"
      ),
      uiOutput("dynamicIDSelector")
    )
  ),
  
  fluidRow(
    column(
      width = 8,
      div(
        plotOutput("metabPlot", height = "500px"),
        style = "max-width: 1100px;"  # Example narrower width
      )),
    column(
      width = 4,
      h4("Fold Change and FDR Summary"),
      div(
        DTOutput("metabTable", height = "350px"),
        style = "max-width: 400px;"  # Example narrower width
      )
    )
  )
)


# Server
server <- function(input, output, session) {
  
  output$dynamicIDSelector <- renderUI({
    selectizeInput(
      "metabID",
      "Search and Select ID:",
      choices = NULL,
      selected = "",
      multiple = FALSE,
      options = list(
        placeholder = '🔍 Type a gene symbol or FBgn ID...',
        create = TRUE
      )
    )
  })
  
  observeEvent(input$datasetType, {
    if (input$datasetType == "G85R-TKT Metabolomics") {
      choices <- rownames(norm.data)
    } else {
      symbols <- GeneIDKey$Symbol
      fbgns <- GeneIDKey$FBgn
      display <- paste0(symbols, " (", fbgns, ")")
      choices <- setNames(fbgns, display)
    }
    
    updateSelectizeInput(
      session,
      "metabID",
      choices = c("", choices),
      server = TRUE
    )
  })
  
  
  
  
  
  # Plot
  output$metabPlot <- renderPlot({
    req(input$metabID, input$datasetType)
    id <- input$metabID
    
    if (input$datasetType == "G85R-TKT Metabolomics") {
      plot.metab(id)
    } else if (input$datasetType == "G85R-TKT Transcriptomics") {
      plot.G85R(id)
    } else if (input$datasetType == "A4V Transcriptomics") {
      plot.A4V(id)
    }
  })
  
  
  # FC/FDR Table
  output$metabTable <- renderDT({
    req(input$metabID, input$datasetType)
    id <- input$metabID
    
    df <- NULL
    
    if (input$datasetType == "G85R-TKT Metabolomics") {
      if (id %in% rownames(Metab.FC) && id %in% rownames(Metab.FDR)) {
        df <- data.frame(
          Comparison = colnames(Metab.FC),
          FoldChange = signif(as.numeric(Metab.FC[id, ]), 4),
          FDR = signif(as.numeric(Metab.FDR[id, ]), 4)
        )
      }
    } else {
      match_row <- GeneIDKey[GeneIDKey$FBgn == id | GeneIDKey$Symbol == id, ]
      
      if (nrow(match_row) == 1) {
        FBgn <- match_row$FBgn
        
        if (input$datasetType == "G85R-TKT Transcriptomics") {
          if (FBgn %in% rownames(TKT.FC) && FBgn %in% rownames(TKT.FDR)) {
            df <- data.frame(
              Comparison = colnames(TKT.FC),
              FoldChange = signif(as.numeric(TKT.FC[FBgn, ]), 4),
              FDR = signif(as.numeric(TKT.FDR[FBgn, ]), 4)
            )
          }
        } else if (input$datasetType == "A4V Transcriptomics") {
          if (FBgn %in% rownames(A4V.FC) && FBgn %in% rownames(A4V.FDR)) {
            df <- data.frame(
              Comparison = colnames(A4V.FC),
              FoldChange = signif(as.numeric(A4V.FC[FBgn, ]), 4),
              FDR = signif(as.numeric(A4V.FDR[FBgn, ]), 4)
            )
          }
        }
      }
    }
    
    if (!is.null(df)) {
      datatable(df, options = list(
        pageLength = 6,
        dom = 't',
        paging = FALSE,
        scrollX = TRUE,
        scrollCollapse = TRUE,
        scrollY = "350px",  # adjust height as needed
        scroller = TRUE,
        rowCallback = JS(
          'function(row, data) {
           if (parseFloat(data[2]) < 0.05) {
             $("td:eq(2)", row).css("background-color", "#ffcccc");
           }
         }'
        )
      ),
      rownames = FALSE,  # 👈 this hides the row numbers
      class = 'compact stripe hover',
      selection = 'none',
      style = 'bootstrap')
    }
  })
  
  
}

# Run app
shinyApp(ui = ui, server = server)
