library(shiny)
library(gplots)
library(DT)
library(visNetwork)
#library(rsconnect)
#deployApp("/Users/johncsantiago/Documents/GitHub/WhartonLab/WhartonLabOmicsDataApp/WhartonOmicsDataViewer")

source("https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/WhartonLabOmicsDataApp/WhartonLabOmicsDataAppFunctions.R")

NETWORK_DATASET <- "G85R-TKT Network (Metab + Enzyme)"

ui <- fluidPage(
  titlePanel("Wharton Lab Omics Data Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      selectInput(
        "datasetType",
        "Select Data Type:",
        choices = c(
          "G85R-TKT Metabolomics",
          "G85R-TKT Transcriptomics",
          "A4V Transcriptomics",
          NETWORK_DATASET
        ),
        selected = "G85R-TKT Metabolomics"
      ),
      
      # Dynamic controls: ID selector for omics, network controls for network
      uiOutput("dynamicControls"),
      
      # A little spacing
      br()
    ),
    
    mainPanel(
      width = 9,
      
      # Right side: big viz
      conditionalPanel(
        condition = sprintf("input.datasetType != '%s'", NETWORK_DATASET),
        plotOutput("metabPlot", height = "650px", width = "100%")
      ),
      
      conditionalPanel(
        condition = sprintf("input.datasetType == '%s'", NETWORK_DATASET),
        visNetworkOutput("network", height = "900px", width = "100%")
      ),
      
      br(),
      
      # Table (only for non-network modes)
      conditionalPanel(
        condition = sprintf("input.datasetType != '%s'", NETWORK_DATASET),
        h4("Fold Change and FDR Summary"),
        DTOutput("metabTable")
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$dynamicControls <- renderUI({
    if (input$datasetType == NETWORK_DATASET) {
      tagList(
        h4("Network Controls"),
        selectInput("enzyme_comp", "Transcriptomics (enzyme) comparison", choices = ENZYME_CHOICES, selected = "GRxWT.FC"),
        numericInput("fccutoff", "FC cutoff (abs)", value = 1.25, min = 0, step = 0.1),
        numericInput("meansumcutoff", "Mean CPM sum cutoff", value = 5, min = 0, step = 1)
      )
    } else {
      selectizeInput(
        "metabID",
        "Search and Select ID:",
        choices = NULL,
        selected = "",
        multiple = FALSE,
        options = list(
          placeholder = "🔍 Type a gene symbol or FBgn ID...",
          create = TRUE
        )
      )
    }
  })
  
  # Update ID selector choices when dataset changes (only for non-network modes)
  observeEvent(input$datasetType, {
    if (input$datasetType == NETWORK_DATASET) return()
    
    if (input$datasetType == "G85R-TKT Metabolomics") {
      choices <- rownames(norm.data)
      updateSelectizeInput(session, "metabID", choices = c("", choices), server = TRUE)
    } else {
      symbols <- GeneIDKey$Symbol
      fbgns   <- GeneIDKey$FBgn
      display <- paste0(symbols, " (", fbgns, ")")
      choices <- setNames(fbgns, display)  # value=FBgn, label="Symbol (FBgn)"
      updateSelectizeInput(session, "metabID", choices = c("", choices), server = TRUE)
    }
  }, ignoreInit = FALSE)
  
  # Plot (non-network)
  output$metabPlot <- renderPlot({
    req(input$datasetType)
    validate(need(input$datasetType != NETWORK_DATASET, ""))
    
    req(input$metabID)
    id <- input$metabID
    
    if (input$datasetType == "G85R-TKT Metabolomics") {
      plot.metab(id)
    } else if (input$datasetType == "G85R-TKT Transcriptomics") {
      plot.G85R(id)
    } else if (input$datasetType == "A4V Transcriptomics") {
      plot.A4V(id)
    }
  })
  
  # Network (bigger)
  output$network <- renderVisNetwork({
    req(input$datasetType == NETWORK_DATASET)
    req(input$enzyme_comp, input$fccutoff, input$meansumcutoff)
    
    validate(need(exists("fullnetwork"), "fullnetwork() not found."))
    
    fullnetwork(
      enzyme = input$enzyme_comp,
      fccutoff = input$fccutoff,
      meansumcutoff = input$meansumcutoff
    )
  })
  
  # Table (non-network)
  output$metabTable <- renderDT({
    req(input$datasetType)
    validate(need(input$datasetType != NETWORK_DATASET, ""))
    
    req(input$metabID)
    id <- input$metabID
    
    df <- NULL
    
    if (input$datasetType == "G85R-TKT Metabolomics") {
      if (id %in% rownames(Metab.FC) && id %in% rownames(Metab.FDR)) {
        df <- data.frame(
          Comparison = colnames(Metab.FC),
          FoldChange = signif(as.numeric(Metab.FC[id, ]), 4),
          FDR        = signif(as.numeric(Metab.FDR[id, ]), 4)
        )
      }
    } else {
      match_row <- GeneIDKey[GeneIDKey$FBgn == id | GeneIDKey$Symbol == id, , drop = FALSE]
      if (nrow(match_row) >= 1) {
        FBgn <- match_row$FBgn[1]
        
        if (input$datasetType == "G85R-TKT Transcriptomics") {
          if (FBgn %in% rownames(TKT.FC) && FBgn %in% rownames(TKT.FDR)) {
            df <- data.frame(
              Comparison = colnames(TKT.FC),
              FoldChange = signif(as.numeric(TKT.FC[FBgn, ]), 4),
              FDR        = signif(as.numeric(TKT.FDR[FBgn, ]), 4)
            )
          }
        } else if (input$datasetType == "A4V Transcriptomics") {
          if (FBgn %in% rownames(A4V.FC) && FBgn %in% rownames(A4V.FDR)) {
            df <- data.frame(
              Comparison = colnames(A4V.FC),
              FoldChange = signif(as.numeric(A4V.FC[FBgn, ]), 4),
              FDR        = signif(as.numeric(A4V.FDR[FBgn, ]), 4)
            )
          }
        }
      }
    }
    
    if (is.null(df)) return(NULL)
    
    datatable(
      df,
      options = list(
        pageLength = 6,
        dom = "t",
        paging = FALSE,
        scrollX = TRUE,
        scrollCollapse = TRUE,
        scrollY = "300px",
        scroller = TRUE,
        rowCallback = JS(
          "function(row, data) {
             if (parseFloat(data[2]) < 0.05) {
               $('td:eq(2)', row).css('background-color', '#ffcccc');
             }
           }"
        )
      ),
      rownames = FALSE,
      class = "compact stripe hover",
      selection = "none",
      style = "bootstrap"
    )
  })
}

shinyApp(ui = ui, server = server)
