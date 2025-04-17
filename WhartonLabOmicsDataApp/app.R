library(shiny)
library(gplots)
library(DT)

# Load data from GitHub
git.dir <- "https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/GeneralDataFiles/"

GeneIDKey = read.csv(paste0(git.dir, "GeneIDKey.csv"), row.names = 1)

TKT.cpm = read.csv(paste0(git.dir,"TKT_cpmdata.csv"), row.names = 1)

TKT.groups  = read.csv(paste0(git.dir, "TKT.metadata.csv"), row.names = 1)

TKT.groups = TKT.groups[colnames(TKT.cpm),]

TKT.FC = read.csv(paste0(git.dir, "TKT.EdgeR.FCTable.csv"), row.names = 1)
TKT.FDR = read.csv(paste0(git.dir, "TKT.EdgeR.FDRTable.csv"), row.names = 1)

A4V.cpm = read.csv(paste0(git.dir, "A4V.cpmdata.csv"), row.names = 1)

A4V.FC = read.csv(paste0(git.dir, "A4V.FCdata.csv"), row.names = 1)

A4V.FDR = read.csv(paste0(git.dir, "A4V.FDRdata.csv"), row.names = 1)

Metab.data <- read.csv(paste0(git.dir, "RawMetabolomicData.csv"), row.names = 1)
Metab.Meta <- read.csv(paste0(git.dir, "MetabolomicMetadata.csv"), row.names = 1)
Metab.FC <- read.csv(paste0(git.dir, "MetaboliteFCs.csv"), row.names = 1)
Metab.FDR <- read.csv(paste0(git.dir, "MetaboliteFDRs.csv"), row.names = 1)

Metab.FC = Metab.FC[,-1]
Metab.FDR = Metab.FDR[,-1]

# KEGG mapping (optional future use)
KEGG.Key <- setNames(Metab.data[-c(1:2), 2], row.names(Metab.data)[-c(1:2)])

# Normalize
norm.func <- function(data) {
  nd <- (as.numeric(data) / Metab.Meta$TIC) * 1000
  return(nd)
}
norm.data <- t(apply(Metab.data[3:nrow(Metab.data), 3:ncol(Metab.data)], 1, norm.func))
colnames(norm.data) <- colnames(Metab.data)[3:ncol(Metab.data)]
colnames(norm.data) <- Metab.Meta[colnames(norm.data), "Genotypes"]


library(gplots)


##Function to plot all G85R groups
plot.G85R = function(ID){
  
  ID.row = GeneIDKey[c(grep(ID, GeneIDKey$FBgn),
                       grep(ID, GeneIDKey$Symbol)),]
  
  if(nrow(ID.row) < 50 & nrow(ID.row) > 0){
    if(nrow(ID.row) == 1
       | length(intersect(GeneIDKey$FBgn, ID)) == 1
       | length(intersect(GeneIDKey$Symbol, ID)) == 1){
      if(nrow(ID.row) == 1){
        FBgn = ID.row$FBgn
        name = ID.row$Symbol
      }
      if(length(intersect(GeneIDKey$FBgn, ID)) == 1){
        FBgn = ID.row$FBgn[ID.row$FBgn == ID]
        name = ID.row$Symbol[ID.row$FBgn == ID]
      }
      if(length(intersect(GeneIDKey$Symbol, ID)) == 1){
        FBgn = ID.row$FBgn[ID.row$Symbol == ID]
        name = ID.row$Symbol[ID.row$Symbol == ID]
      }
      
      
      
      gene = data.frame(cpm = as.numeric(TKT.cpm[FBgn,]), 
                        condition = TKT.groups$Group, 
                        group = paste0(TKT.groups$Genotype, TKT.groups$Sex))
      
      gene$group = factor(gene$group, levels = c("WTF", "GRF", "WTM", "GRM"))
      gene$condition = factor(gene$condition, levels = c('WT.F', "TktDfWT.F", "TktOEWT.F",
                                                         'GR.F', "TktDfGR.F", "TktOEGR.F",
                                                         'WT.M', "TktDfWT.M", "TktOEWT.M",
                                                         'GR.M', "TktDfGR.M", "TktOEGR.M"))
      
      if(nrow(na.omit(gene)) > 0){
        par(mar = c(5,7,5,2))
        plot(x = NA,
             y = NA,
             type = 'n',
             ylim = c(min(gene$cpm), max(gene$cpm)),
             xlim = c(.5, 12.5),
             xaxt = "n",
             ylab="",
             yaxt = "n",
             xlab = NA)
        
        rect(0, 0-(1.5*max(gene$cpm)), 13, 1.5*max(gene$cpm), col = 'gray95')
        
        boxplot(gene$cpm~gene$condition, 
                xlab="",
                las = 2, 
                cex.axis = .8,
                xaxt = 'none',
                yaxt = 'none',
                boxwex = .75,
                boxlwd = 2,
                lwd = 2,
                cex.ticks = 2,
                add = T,
                col = c(rep("honeydew", 3),
                        rep("thistle1", 3),
                        rep("lightcyan", 3),
                        rep("lightyellow", 3)))
        
        box(lwd = 2)
        
        points(x=as.numeric(factor(gene$condition)),
               y=gene$cpm,
               cex=1,
               pch=21,
               bg=c("darkgreen", 'firebrick',"dodgerblue","gold")[as.numeric(factor(gene$group))])
        
        axis(side = 1,
             at = 1:12,
             labels = rep(c("Control", "Df", "OE"), 4),
             line = 0,
             cex.axis = 1,
             lwd.ticks = 2)
        
        axis(side = 2,
             line = 0,
             cex.axis = 1,
             lwd.ticks = 2,
             las = 2)
        
        axis(side = 2,
             at = .5*(max(gene$cpm) + min(gene$cpm)),
             label = "Counts per million reads",
             tick = F,
             padj = -3,
             cex.axis = 1.25)
        
        axis(side = 1,
             at = c(2, 5, 8, 11),
             labels = c("Silent", "G85R","Silent", "G85R"),
             tick = F,
             padj = 2,
             cex.axis = 1.5)
        axis(side = 3,
             at = c(3.5, 9.5),
             labels = c("Female", "Male"),
             tick = F,
             padj = .5,
             cex.axis = 1.5)
        
        axis(side = 3,
             at = 6.5,
             labels = name,
             tick = F,
             padj = -1.5,
             cex.axis = 2)
        
        #axis(side = 3,
         #    at = c(3.5,9.5),
          #   tick = T,
           #  outer = T,
            # labels = F,
             #line = -2.25,
             #lwd.ticks = 0)
        
        xline1 = min(gene$cpm) - (max(gene$cpm)-min(gene$cpm))*.04
        xline2 = min(gene$cpm) - (max(gene$cpm)-min(gene$cpm))*.19
        lines(x=c(3.5,3.5), y = c(xline1, xline2), lty = 1, lwd = 2, xpd = T)
        lines(x=c(6.5,6.5), y = c(xline1, xline2), lty = 1, lwd = 2, xpd = T)
        lines(x=c(9.5,9.5), y = c(xline1, xline2), lty = 1, lwd = 2, xpd = T)
        abline(v = 6.5, lty = 2, lwd = 2)
      }
    }
  }
  
  if(nrow(ID.row) != 1
     | length(intersect(GeneIDKey$FBgn, ID)) > 1
     | length(intersect(GeneIDKey$Symbol, ID)) > 1){
    print(paste0("Multiple gene symbols match the ID ", ID))
    ID.row
  }
}

##Function to plot a specific gene expression levels in all A4V groups
plot.A4V = function(ID){
  ID.row = GeneIDKey[c(grep(ID, GeneIDKey$FBgn),
                       grep(ID, GeneIDKey$Symbol)),]
  
  if(nrow(ID.row) < 50 & nrow(ID.row) > 0){
    if(nrow(ID.row) == 1
       | length(intersect(GeneIDKey$FBgn, ID)) == 1
       | length(intersect(GeneIDKey$Symbol, ID)) == 1){
      if(nrow(ID.row) == 1){
        FBgn = ID.row$FBgn
        name = ID.row$Symbol
      }
      if(length(intersect(GeneIDKey$FBgn, ID)) == 1){
        FBgn = ID.row$FBgn[ID.row$FBgn == ID]
        name = ID.row$Symbol[ID.row$FBgn == ID]
      }
      if(length(intersect(GeneIDKey$Symbol, ID)) == 1){
        FBgn = ID.row$FBgn[ID.row$Symbol == ID]
        name = ID.row$Symbol[ID.row$Symbol == ID]
      }
      A4V.groups = as.vector(colnames(A4V.cpm))
      
      A4V.groups = substr(A4V.groups,1, nchar(A4V.groups)-1)
      
      
      if(nrow(na.omit(A4V.cpm[FBgn,])) > 0){
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        gene = data.frame(cpm = as.numeric(A4V.cpm[FBgn,]),
                          condition = A4V.groups,
                          group = paste0(substr(A4V.groups, 1,1), substr(A4V.groups, start = nchar(as.vector(A4V.groups))-1,nchar(as.vector(A4V.groups))-1)))
        
        gene$condition = factor(gene$condition, levels = c("S3FH", "S9FH", "S40FH",
                                                           "S3FT", "S9FT", "S40FT",
                                                           "S3FA", "S9FA", "S40FA",
                                                           "A3FH", "A9FH", "A40FH",
                                                           "A3FT", "A9FT", "A40FT",
                                                           "A3FA", "A9FA", "A40FA",
                                                           "S3MH", "S9MH",
                                                           "S3MT", "S9MT",
                                                           "S3MA", "S9MA",
                                                           "A3MH", "A9MH",
                                                           "A3MT", "A9MT",
                                                           "A3MA", "A9MA"))
        gene$group = factor(gene$group, levels = c("SF", "AF", "SM", "AM"))
        
        boxplot(gene$cpm~gene$condition, xlab="",     
                main= name,ylab="Counts per Million Reads", las = 2, cex.axis = 1, cex.main = 2)
        points(x=as.numeric(factor(gene$condition)),y=gene$cpm,cex=1,pch=21,bg=c('firebrick',"gold","darkgreen","dodgerblue")[as.numeric(factor(gene$group))])
        legend('topright',
               inset=c(-0.1,0),
               legend = c('Silent F', 'A4V F', 'Silent M', 'A4V M'), 
               fill = c('firebrick', "gold", "darkgreen", "dodgerblue"),
               cex = 1,
               bty = 'n',
               pt.cex = .5)
      }
    }
    
  }
}

##Plot levels of a specific metabolite
plot.metab = function(metabID){
  
  ID.row = grep(metabID, row.names(norm.data))
  
  
  if(length(ID.row) < 50 & length(ID.row) > 0){
    if(length(ID.row) == 1 | length(intersect(row.names(norm.data), metabID)) == 1){
      
      if(length(ID.row) == 1){
        metab = row.names(norm.data)[ID.row]
        name = row.names(norm.data)[ID.row]
      }
      
      if(length(intersect(row.names(norm.data), metabID)) == 1){
        metab = metabID
        name = metabID
      }
      
      
      
      metab = data.frame(data = as.numeric(norm.data[metab,]), 
                         group = colnames(norm.data), 
                         geno = "WT")
      
      
      metab$group = factor(metab$group, levels = c("WT", "TktDfWT", "TktOEWT",
                                                   "G85R", "TktDfGR", "TktOEGR"))
      
      metab$geno[grep("G", metab$group)] = "G85R"
      metab$geno = factor(metab$geno, 
                          levels = c("WT", "G85R"))
      
      par(mar = c(5,7,5,2))
      plot(x = NA,
           y = NA,
           type = 'n',
           ylim = c(min(na.omit(metab$data)), max(na.omit(metab$data))),
           xlim = c(.5, 6.5),
           xaxt = "n",
           ylab="",
           yaxt = "n",
           xlab = NA)
      
      rect(0, 0-(1.5*max(na.omit(metab$data))), 13, 1.5*max(na.omit(metab$data)), col = 'gray95')
      
      boxplot(metab$data~metab$group, 
              xlab="",
              las = 2, 
              cex.axis = .8,
              xaxt = 'none',
              yaxt = 'none',
              boxwex = .75,
              boxlwd = 2,
              lwd = 2,
              cex.ticks = 2,
              add = T,
              col = c(rep("lightyellow", 3),
                      rep("lightgreen", 3)))
      
      box(lwd = 2)
      
      points(x=as.numeric(factor(metab$group)),
             y=metab$data,
             cex=1.25,
             pch=21,
             bg=c('gold3',"darkgreen")[as.numeric(factor(metab$geno))])
      
      axis(side = 1,
           at = 1:12,
           labels = rep(c("Control", "Df", "OE"), 4),
           line = 0,
           cex.axis = .75,
           lwd.ticks = 2)
      
      axis(side = 2,
           line = 0,
           cex.axis = 1,
           lwd.ticks = 2,
           las = 2)
      
      axis(side = 1,
           at = c(2, 5, 8, 11),
           labels = c("Silent", "G85R","Silent", "G85R"),
           tick = F,
           padj = 2,
           cex.axis = 1)
      
      axis(side = 3,
           at = 3.5,
           labels = name,
           tick = F,
           padj = -.5,
           cex.axis = 2)
      
    }
    
    if(length(ID.row) != 1 | length(intersect(row.names(norm.data), metabID)) > 1){
      print(paste0("Multiple metabolites match the ID ", metabID))
      norm.data[ID.row,]
    }
    
  }
}

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
        DTOutput("metabTable", height = "300px"),
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
        scrollY = "300px",  # adjust height as needed
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
