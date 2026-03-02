library(shiny)
library(gplots)
library(DT)
library(visNetwork)
#library(rsconnect)
#deployApp("/Users/johncsantiago/Documents/GitHub/WhartonLab/WhartonLabOmicsDataApp/WhartonOmicsDataViewer")

#source("https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/WhartonLabOmicsDataApp/WhartonLabOmicsDataAppFunctions.R")

source("/Users/johncsantiago/Documents/GitHub/WhartonLab/WhartonLabOmicsDataApp/WhartonLabOmicsDataAppFunctions.R")

#METAB_CHOICES <- c(
#  "GRxWT.C","GRxWT.Df","GRxWT.OE","GR.CxDf","WT.CxDf",
#  "GR.CxOE","WT.CxOE","GR.DfxWT.C","GR.OExWT.C","none"
#)

#METAB_CHOICES <- c(
#  "GRxWT","GRxWT.Df","GRxWT.OE","GR.DfxC","WT.DfxC",
#  "GR.OExC","WT.OExC","GR.DfxWT","GR.OExWT","none"
#)

G85R.sigTable = setNames(
  c(
    "GRF.DfxC","GRF.OExC","WTF.DfxC","WTF.OExC",
    "GRxWT.F", "GRxWT.FDf","GRxWT.FOE","GRDfxWT.F","GROExWT.F","GRM.DfxC",
    "GRM.OExC","WTM.DfxC","WTM.OExC","GRxWT.M","GRxWT.MDf",
    "GRxWT.MOE","GRDfxWT.M","GROExWT.M","GR.FxM","WT.FxM"
  ),
  c(
  "GRF.CxDf","GRF.CxOE","WTF.CxDf","WTF.CxOE","GRxWT.FC",
  "GRxWT.FDf","GRxWT.FOE","GRDfxWT.F","GROExWT.F","GRM.CxDf",
  "GRM.CxOE","WTM.CxDf","WTM.CxOE","GRxWT.MC","GRxWT.MDf",
  "GRxWT.MOE","GRDfxWT.M","GROExWT.M","GR.FxM","WT.FxM"
  )
)

Metab.sigTable = setNames(c(
    "GRxWT","GRxWT.Df","GRxWT.OE","GR.DfxC","WT.DfxC",
    "GR.OExC","WT.OExC","GR.DfxWT","GR.OExWT","none"
  ),
  c(
    "GRxWT.C","GRxWT.Df","GRxWT.OE","GR.CxDf","WT.CxDf",
    "GR.CxOE","WT.CxOE","GR.DfxWT.C","GR.OExWT.C","none"
    )
)

MX_GROUPS <- list(
  "Group 1: Within-genotype TKT" = c("WT.DfxC", "WT.OExC",
                                     "GR.DfxC", "GR.OExC"),
  "Group 2: G85R vs Silent" = c("GRxWT", 
                                "GR.DfxWT", "GRxWT.Df",
                                "GR.OExWT", "GRxWT.OE")
)

TX_GROUPS <- list(
  "Group 1: Female within-genotype TKT" = c("WTF.DfxC","WTF.OExC","GRF.DfxC","GRF.OExC"),
  "Group 2: Female G85R vs Silent" = c("GRxWT.F","GRDfxWT.F","GRxWT.FDf","GROExWT.F","GRxWT.FOE"),
  "Group 3: Male within-genotype TKT" = c("WTM.DfxC","WTM.OExC","GRM.DfxC","GRM.OExC"),
  "Group 4: Male G85R vs Silent" = c("GRxWT.M","GRDfxWT.M","GRxWT.MDf","GROExWT.M","GRxWT.MOE"),
  "Group 5: Sex effects" = c("GR.FxM","WT.FxM")
)

AX_GROUPS <- list(
  "Group 1: Female silent across timepoints" = c(
    "S3v9FH", "S3v40FH", "S9v40FH",
    "S3v9FT", "S3v40FT", "S9v40FT",
    "S3v9FA", "S3v40FA", "S9v40FA"),
  "Group 2: Female A4V across timepoints" = c(
    "A3v9FH", "A3v40FH", "A9v40FH",
    "A3v9FT", "A3v40FT", "A9v40FT",
    "A3v9FA", "A3v40FA", "A9v40FA"),
  "Group 3: Female across genotype" = c(
    "AvS3FH", "AvS9FH", "AvS40FH",
    "AvS3FT", "AvS9FT", "AvS40FT",
    "AvS3FA", "AvS9FA", "AvS40FA"),
  "Group 4: Male silent across timepoints" = c(
    "S3v9MH",
    "S3v9MT",
    "S3v9MA"),
  "Group 5: Male A4V across timepoints" = c(
    "A3v9MH",
    "A3v9MT",
    "A3v9MA"),
  "Group 6: Male across genotype" = c(
    "AvS3MH", "AvS9MH",
    "AvS3MT", "AvS9MT",
    "AvS3MA", "AvS9MA"),
  "Group 7: Silent across sex" = c(
    "S3FvMH", "S9FvMH",
    "S3FvMT", "S9FvMT",
    "S3FvMA", "S9FvMA"),
  "Group 8: A4V across sex" = c(
    "A3FvMH", "A9FvMH",
    "A3FvMT", "A9FvMT",
    "A3FvMA", "A9FvMA")
)

G85R.sigTable.order = c(
  "WTF.CxDf",  "WTF.CxOE", "GRF.CxDf",  "GRF.CxOE",
  "GRxWT.FC", "GRDfxWT.F", "GRxWT.FDf","GROExWT.F","GRxWT.FOE",
  "WTM.CxDf",  "WTM.CxOE", "GRM.CxDf",  "GRM.CxOE",
  "GRxWT.MC", "GRDfxWT.M", "GRxWT.MDf", "GROExWT.M", "GRxWT.MOE",
  "GR.FxM",    "WT.FxM"
)


MX_ORDER <- unlist(MX_GROUPS, use.names = FALSE)

TX_ORDER <- unlist(TX_GROUPS, use.names = FALSE)

AX_ORDER <- unlist(AX_GROUPS, use.names = FALSE)

MX_SECTION_FOR <- function(comp) {
  # returns section label for each Comparison value
  out <- rep(NA_character_, length(comp))
  for (nm in names(MX_GROUPS)) {
    out[comp %in% MX_GROUPS[[nm]]] <- nm
  }
  out
}

TX_SECTION_FOR <- function(comp) {
  # returns section label for each Comparison value
  out <- rep(NA_character_, length(comp))
  for (nm in names(TX_GROUPS)) {
    out[comp %in% TX_GROUPS[[nm]]] <- nm
  }
  out
}

AX_SECTION_FOR <- function(comp) {
  # returns section label for each Comparison value
  out <- rep(NA_character_, length(comp))
  for (nm in names(AX_GROUPS)) {
    out[comp %in% AX_GROUPS[[nm]]] <- nm
  }
  out
}

NETWORK_DATASET <- "G85R-TKT Network (Metab + Enzyme)"

ui <- fluidPage(
  titlePanel("Wharton Lab Omics Data Viewer"),
    
    # 👇 Put CSS HERE (top-level inside fluidPage)
    tags$style(HTML("
    tr.dt-rowgroup td.dt-group {
      font-weight: 700;
      font-size: 14px;
      background: #f0f2f5 !important;
      border-top: 3px solid #333 !important;
      border-bottom: 1px solid #bbb !important;
      padding: 8px !important;
    }
  ")),

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
        selectInput("enzyme_comp", "Transcriptomics (enzyme) comparison", choices = as.character(G85R.sigTable[ENZYME_CHOICES]), selected = "GRxWT.FC"),
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
          Comparison = as.character(Metab.sigTable[colnames(Metab.FC)]),
          FoldChange = signif(as.numeric(Metab.FC[id, ]), 4),
          FDR        = signif(as.numeric(Metab.FDR[id, ]), 4),
          stringsAsFactors = FALSE
        )
        
        # keep only comparisons that exist (safe if names change)
        keep <- df$Comparison %in% MX_ORDER
        df <- df[keep, , drop = FALSE]
        
        # order rows to match your group order
        df$Comparison <- factor(df$Comparison, levels = MX_ORDER)
        df <- df[order(df$Comparison), , drop = FALSE]
        df$Comparison <- as.character(df$Comparison)
        
        # add Section column (for RowGroup)
        df$Section <- MX_SECTION_FOR(df$Comparison)
        df$Section <- factor(df$Section, levels = names(MX_GROUPS))        
        
      }
    } else {
      match_row <- GeneIDKey[GeneIDKey$FBgn == id | GeneIDKey$Symbol == id, , drop = FALSE]
      if (nrow(match_row) >= 1) {
        FBgn <- match_row$FBgn[1]
        
        if (input$datasetType == "G85R-TKT Transcriptomics") {
          if (FBgn %in% rownames(TKT.FC) && FBgn %in% rownames(TKT.FDR)) {
            df <- data.frame(
              Comparison = as.character(G85R.sigTable[colnames(TKT.FC)]),
              FoldChange = signif(as.numeric(TKT.FC[FBgn, ]), 4),
              FDR        = signif(as.numeric(TKT.FDR[FBgn, ]), 4),
              stringsAsFactors = FALSE
            )
            
            # keep only comparisons that exist (safe if names change)
            keep <- df$Comparison %in% TX_ORDER
            df <- df[keep, , drop = FALSE]
            
            # order rows to match your group order
            df$Comparison <- factor(df$Comparison, levels = TX_ORDER)
            df <- df[order(df$Comparison), , drop = FALSE]
            df$Comparison <- as.character(df$Comparison)
            
            # add Section column (for RowGroup)
            df$Section <- TX_SECTION_FOR(df$Comparison)
            df$Section <- factor(df$Section, levels = names(TX_GROUPS))
          }
        } else if (input$datasetType == "A4V Transcriptomics") {
          if (FBgn %in% rownames(A4V.FC) && FBgn %in% rownames(A4V.FDR)) {
            df <- data.frame(
              Comparison = colnames(A4V.FC),
              FoldChange = signif(as.numeric(A4V.FC[FBgn, ]), 4),
              FDR        = signif(as.numeric(A4V.FDR[FBgn, ]), 4),
              stringsAsFactors = FALSE
            )
            
            # keep only comparisons that exist (safe if names change)
            keep <- df$Comparison %in% AX_ORDER
            df <- df[keep, , drop = FALSE]
            
            # order rows to match your group order
            df$Comparison <- factor(df$Comparison, levels = AX_ORDER)
            df <- df[order(df$Comparison), , drop = FALSE]
            df$Comparison <- as.character(df$Comparison)
            
            # add Section column (for RowGroup)
            df$Section <- AX_SECTION_FOR(df$Comparison)
            df$Section <- factor(df$Section, levels = names(AX_GROUPS))
            
          }
        }
      }
    }
    
    if (is.null(df)) return(NULL)
    
    datatable(
      df,
      rownames = FALSE,
      extensions = c("RowGroup"),
      options = list(
        dom = "t",
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "350px",
        rowGroup = list(
          dataSrc = 3,          # <-- "Section" column index (0-based in JS, but DT uses column position)
          startRender = JS(
            "function(rows, group) {
           return $('<tr/>')
             .append('<td colspan=\"3\" class=\"dt-group\">'+group+'</td>');
         }"
          )
        ),
        columnDefs = list(
          list(targets = 3, visible = FALSE)  # hide Section column
        ),
        rowCallback = JS(
          # thick border above first row of each section (optional; RowGroup header already separates)
          "function(row, data, index) {
         // data[3] is Section (hidden but still present in data array)
         // add a subtle style to numeric columns if you want:
       }"
        )
      ),
      class = "compact stripe hover",
      selection = "none",
      style = "bootstrap"
    ) %>%
      # highlight significant FDR cells
      formatStyle(
        "FDR",
        backgroundColor = styleInterval(0.05, c("#ffcccc", NA))
      ) %>%
      formatStyle(
        c("Comparison","FoldChange","FDR"),
        fontSize = "85%"
      )
  })
}

shinyApp(ui = ui, server = server)
