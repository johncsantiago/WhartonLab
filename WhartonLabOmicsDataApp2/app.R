# app.R  ---------------------------------------------------------------
library(shiny)
library(gplots)
library(DT)
library(visNetwork)

source("https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/WhartonLabOmicsDataApp/WhartonLabOmicsDataAppFunctions.R")

# ---------------- UI ----------------
ui <- fluidPage(
  titlePanel("Wharton Lab Omics Data Viewer"),
  
  fluidRow(
    column(
      width = 12,
      selectInput(
        "datasetType",
        "Select Data Type:",
        choices = c("G85R-TKT Metabolomics", "G85R-TKT Transcriptomics", "A4V Transcriptomics"),
        selected = "G85R-TKT Metabolomics"
      ),
      uiOutput("dynamicViewControls"),
      uiOutput("dynamicIDSelector")
    )
  ),
  
  fluidRow(
    column(
      width = 8,
      div(
        # Show plot OR network, depending on viewMode
        conditionalPanel(
          condition = "input.viewMode == 'Plot'",
          plotOutput("metabPlot", height = "500px")
        ),
        conditionalPanel(
          condition = "input.viewMode == 'Network'",
          visNetworkOutput("network", height = "500px")
        ),
        style = "max-width: 1100px;"
      )
    ),
    column(
      width = 4,
      h4("Fold Change and FDR Summary"),
      div(
        DTOutput("metabTable", height = "300px"),
        style = "max-width: 400px;"
      )
    )
  )
)

# ---------------- Server ----------------
server <- function(input, output, session) {
  
  # View controls (Plot vs Network) + network parameters shown only when relevant
  output$dynamicViewControls <- renderUI({
    # Only allow Network for the G85R-TKT datasets
    allow_network <- input$datasetType %in% c("G85R-TKT Metabolomics", "G85R-TKT Transcriptomics")
    
    tagList(
      selectInput(
        "viewMode",
        "View:",
        choices = if (allow_network) c("Plot", "Network") else c("Plot"),
        selected = "Plot"
      ),
      
      conditionalPanel(
        condition = "input.viewMode == 'Network'",
        fluidRow(
          column(
            6,
            selectInput("metab_comp", "Metabolomics comparison", choices = METAB_CHOICES, selected = "GRxWT.C")
          ),
          column(
            6,
            selectInput("enzyme_comp", "Transcriptomics comparison", choices = ENZYME_CHOICES, selected = "GRxWT.FC")
          )
        ),
        fluidRow(
          column(6, numericInput("fccutoff", "FC cutoff (abs)", value = 0, min = 0, step = 0.1)),
          column(6, numericInput("meansumcutoff", "Mean CPM sum cutoff", value = 0, min = 0, step = 1))
        )
      )
    )
  })
  
  # ID selector
  output$dynamicIDSelector <- renderUI({
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
  })
  
  # Update selector options when dataset changes
  observeEvent(input$datasetType, {
    if (input$datasetType == "G85R-TKT Metabolomics") {
      choices <- rownames(norm.data)
      updateSelectizeInput(session, "metabID", choices = c("", choices), server = TRUE)
    } else {
      symbols <- GeneIDKey$Symbol
      fbgns   <- GeneIDKey$FBgn
      display <- paste0(symbols, " (", fbgns, ")")
      choices <- setNames(fbgns, display) # value=FBgn, label="Symbol (FBgn)"
      updateSelectizeInput(session, "metabID", choices = c("", choices), server = TRUE)
    }
    
    # If user was on Network and switched to A4V, force Plot
    if (input$datasetType == "A4V Transcriptomics") {
      updateSelectInput(session, "viewMode", selected = "Plot")
    }
  }, ignoreInit = FALSE)
  
  # ---- Plot view ----
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
  
  # ---- Network view ----
  output$network <- renderVisNetwork({
    req(input$viewMode == "Network")
    req(input$datasetType %in% c("G85R-TKT Metabolomics", "G85R-TKT Transcriptomics"))
    req(input$metab_comp, input$enzyme_comp)
    
    validate(
      need(exists("fullnetwork"), "fullnetwork() not found. Define it in your Functions.R or paste it into app.R.")
    )
    
    fullnetwork(
      metab = input$metab_comp,
      enzyme = input$enzyme_comp,
      fccutoff = input$fccutoff,
      meansumcutoff = input$meansumcutoff
    )
  })
  
  # ---- FC/FDR Table ----
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
      match_row <- GeneIDKey[GeneIDKey$FBgn == id | GeneIDKey$Symbol == id, , drop = FALSE]
      
      if (nrow(match_row) >= 1) {
        # If multiple matches, take first (rare, but prevents crashes)
        FBgn <- match_row$FBgn[1]
        
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
