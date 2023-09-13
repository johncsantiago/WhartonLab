##library(rsconnect)


library(shiny)
library(plotly)

git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/"
DfKey = read.csv(paste0(git.dir, "DeficiencyScreen/DataFiles/DfKey.csv"), row.names = 1)
row.names(DfKey) = DfKey$Deficiency
DfKey[DfKey$MLE == 0, "MLE"] = .05
GeneIDKey = read.csv(paste0(git.dir,"GeneralDataFiles/GeneIDKey.csv"),row.names = 1)

fbID = data.frame(FBgn = c(rep(GeneIDKey$FBgn, 2)),
                  Symbol = c(rep(GeneIDKey$Symbol, 2)),
                  Position = c(GeneIDKey$Start, GeneIDKey$End))

plot.pos = function(chr, start, end){
  goi = data.frame(Start = start,
                   End = end,
                   CHR = chr)
  range.min = start
  range.max = end
  
  ##genes in range
  gir = na.omit(GeneIDKey[(GeneIDKey$CHR == goi$CHR),])
  gir = na.omit(gir[(((gir$Start >= range.min) & (gir$Start <= range.max)) | 
                       ((gir$End >= range.min) & (gir$End <= range.max)) |
                       ((gir$End >= range.max) & (gir$Start <= range.min))), ])
  gir$x1 = gir$Start
  gir$x2 = gir$End
  
  Dfir = na.omit(DfKey[DfKey$Chr == goi$CHR,])
  Dfir = na.omit(Dfir[(((Dfir$Start >= range.min) & (Dfir$Start <= range.max)) | 
                         ((Dfir$End >= range.min) & (Dfir$End <= range.max)) |
                         ((Dfir$End >= range.max) & (Dfir$Start <= range.min))), ])
  Dfir$x1 = Dfir$Start
  Dfir$x2 = Dfir$End
  
  line.data = gir
  if(min(line.data$Start) < range.min){
    line.data$x2[line.data$Start < range.min] = range.min
  }
  
  if(max(line.data$End) > range.max){
    line.data$x2[line.data$End > range.max] = range.max
  }
  
  bar.data = Dfir
  if(min(bar.data$Start) < range.min){
    bar.data$x1[bar.data$Start < range.min] = range.min
  }
  
  if(max(bar.data$End) > range.max){
    bar.data$x2[bar.data$End > range.max] = range.max
  }
  
  line.data$Strand = line.data$Strand/10
  line.data = line.data[order(line.data$Start),]
  line.data = line.data[, c("Symbol", "Start", "Strand", "End", "x1", "x2")]
  line.data = unique(line.data)
  
  if(nrow(line.data)>1){
    line.data$Strand[seq(2, nrow(line.data), by = 2)] = 2 * line.data$Strand[seq(2, nrow(line.data), by = 2)]
  }
  
  Start.data = data.frame(Symbol = line.data$Symbol, 
                          x = line.data$x1,
                          y = line.data$Strand,
                          Start = line.data$Start,
                          End = line.data$End,
                          shape = '35')
  
  Start.data$shape[Start.data$y < 0] = "triangle-left"
  
  
  End.data = data.frame(Symbol = line.data$Symbol, 
                        x = line.data$x2,
                        y = line.data$Strand,
                        Start = line.data$Start,
                        End = line.data$End,
                        shape = '35')
  
  End.data$shape[End.data$y > 0] = "triangle-right"
  
  
  data = rbind(Start.data, End.data)
  data$color = "red"
  data$color[data$y < 0] = "blue"
  data = data[order(data$Symbol),]
  
  bar.data$color = bar.data$Modifier
  
  bar.data$color[bar.data$color == "Strong Suppressor"] = adjustcolor("brown", alpha.f = 0.3)
  bar.data$color[bar.data$color == "Suppressor"] = adjustcolor("gold", alpha.f = 0.3)
  bar.data$color[bar.data$color == "No Effect"] = adjustcolor("lightgrey", alpha.f = 0.3)
  bar.data$color[bar.data$color == "Enhancer"] = adjustcolor("steelblue", alpha.f = 0.3)
  
  p1.data = data.frame(dfName = bar.data$Deficiency, dfx = bar.data$x1, dfy = as.numeric(bar.data$MLE), dfcolor = bar.data$color, dfStart = bar.data$Start, dfEnd = bar.data$End)
  p2.data = data.frame(dfName = bar.data$Deficiency, dfx = bar.data$x2, dfy = as.numeric(bar.data$MLE), dfcolor = bar.data$color, dfStart = bar.data$Start, dfEnd = bar.data$End)
  
  df.line = data.frame(dfName = rep(NA, nrow(data)),
                       dfx = NA,
                       dfy = NA,
                       dfcolor = NA,
                       dfStart = NA,
                       dfEnd = NA)
  pdata = rbind(p1.data,p2.data)
  df.line[1:nrow(pdata),] = pdata
  data = cbind(data,df.line)
  data$Symbol = factor(data$Symbol, levels = unique(data$Symbol))
  
  
  fig = plot_ly(data,
                x = ~x,
                y = ~y,
                type = 'scatter',
                mode = "lines + markers",
                name = ~Symbol,
                marker = list(color = data$color, symbol = ~shape, size = 5),
                showlegend = F,
                line = list(
                  color = data$color
                ),
                hoverinfo = "text",
                hovertext = paste0(data$Symbol, 
                                   "\nStart: ", data$Start,
                                   "\nEnd: ", data$End))
  
  fig = fig %>% add_lines(x = ~dfx,
                          y = ~dfy,
                          type = 'scatter',
                          mode = "lines",
                          name = ~dfName,
                          showlegend = F,
                          fill = 'tozeroy',
                          fillcolor = ~dfcolor,
                          opacity = .1,
                          line = list(
                            color = 'black',
                            width = 3),
                          marker = list(color = 'black', size = 1),
                          hoverinfo = "text",
                          hovertext = paste0(data$dfName, 
                                             "\nStart: ", data$dfStart,
                                             "\nEnd: ", data$dfEnd))
  
  
  fig = layout(fig,
               xaxis = list(title = goi$CHR, range = c(start, end)),
               yaxis = list(title = "Mean Lifespan Extension", range = c(-4,6)))##,
  ##hovermode = "x unified")
  
  return(fig)
}

plot.gene = function(gene.symbol, start.range, end.range){
  
  goi = data.frame(Start = GeneIDKey[GeneIDKey$Symbol %in% gene.symbol, 'Start'],
                   End = GeneIDKey[GeneIDKey$Symbol %in% gene.symbol, 'End'],
                   CHR = GeneIDKey[GeneIDKey$Symbol %in% gene.symbol, 'CHR'])
  range.min = goi$Start - start.range
  range.max = goi$End + end.range 
  
  ##genes in range
  gir = na.omit(GeneIDKey[(GeneIDKey$CHR == goi$CHR),])
  gir = na.omit(gir[(((gir$Start >= range.min) & (gir$Start <= range.max)) | 
                       ((gir$End >= range.min) & (gir$End <= range.max)) |
                       ((gir$End >= range.max) & (gir$Start <= range.min))), ])
  gir$x1 = gir$Start
  gir$x2 = gir$End
  
  Dfir = na.omit(DfKey[DfKey$Chr == goi$CHR,])
  Dfir = na.omit(Dfir[(((Dfir$Start >= range.min) & (Dfir$Start <= range.max)) | 
                         ((Dfir$End >= range.min) & (Dfir$End <= range.max)) |
                         ((Dfir$End >= range.max) & (Dfir$Start <= range.min))), ])
  Dfir$x1 = Dfir$Start
  Dfir$x2 = Dfir$End
  
  line.data = gir
  if(min(line.data$Start) < range.min){
    line.data$x1[line.data$Start < range.min] = range.min
  }
  
  if(max(line.data$End) > range.max){
    line.data$x2[line.data$End > range.max] = range.max
  }
  
  bar.data = Dfir
  if(nrow(Dfir)>0){
    if(min(bar.data$Start) < range.min){
      bar.data$x1[bar.data$Start < range.min] = range.min
    }
    
    if(max(bar.data$End) > range.max){
      bar.data$x2[bar.data$End > range.max] = range.max
    }
  }
  
  line.data$Strand = line.data$Strand/10
  line.data = line.data[order(line.data$Start),]
  line.data = line.data[, c("Symbol", "Start", "Strand", "End", "x1", "x2")]
  line.data = unique(line.data)
  
  if(nrow(line.data)>1){
    line.data$Strand[seq(2, nrow(line.data), by = 2)] = 2 * line.data$Strand[seq(2, nrow(line.data), by = 2)]
  }
  
  Start.data = data.frame(Symbol = line.data$Symbol, 
                          x = line.data$x1,
                          y = line.data$Strand,
                          Start = line.data$Start,
                          End = line.data$End,
                          shape = '35')
  
  Start.data$shape[Start.data$y < 0] = "triangle-left"
  
  
  End.data = data.frame(Symbol = line.data$Symbol, 
                        x = line.data$x2,
                        y = line.data$Strand,
                        Start = line.data$Start,
                        End = line.data$End,
                        shape = '35')
  
  End.data$shape[End.data$y > 0] = "triangle-right"
  
  
  data = rbind(Start.data, End.data)
  data$color = "red"
  data$color[data$y < 0] = "blue"
  data = data[order(data$Symbol),]
  
  bar.data$color = bar.data$Modifier
  
  bar.data$color[bar.data$color == "Strong Suppressor"] = adjustcolor("brown", alpha.f = 0.3)
  bar.data$color[bar.data$color == "Suppressor"] = adjustcolor("gold", alpha.f = 0.3)
  bar.data$color[bar.data$color == "No Effect"] = adjustcolor("lightgrey", alpha.f = 0.3)
  bar.data$color[bar.data$color == "Enhancer"] = adjustcolor("steelblue", alpha.f = 0.3)
  
  p1.data = data.frame(dfName = bar.data$Deficiency, dfx = bar.data$x1, dfy = as.numeric(bar.data$MLE), dfcolor = bar.data$color, dfStart = bar.data$Start, dfEnd = bar.data$End)
  p2.data = data.frame(dfName = bar.data$Deficiency, dfx = bar.data$x2, dfy = as.numeric(bar.data$MLE), dfcolor = bar.data$color, dfStart = bar.data$Start, dfEnd = bar.data$End)
  
  df.line = data.frame(dfName = rep(NA, nrow(data)),
                       dfx = NA,
                       dfy = NA,
                       dfcolor = NA,
                       dfStart = NA,
                       dfEnd = NA)
  pdata = rbind(p1.data,p2.data)
  if(nrow(pdata)>0){
    df.line[1:nrow(pdata),] = pdata
  }
  data = cbind(data,df.line)
  data$Symbol = factor(data$Symbol, levels = unique(data$Symbol))
  
  
  fig = plot_ly(data,
                x = ~x,
                y = ~y,
                type = 'scatter',
                mode = "lines + markers",
                name = ~Symbol,
                marker = list(color = data$color, symbol = ~shape, size = 5),
                showlegend = F,
                line = list(
                  color = data$color
                ),
                hoverinfo = "text",
                hovertext = paste0(data$Symbol, 
                                   "\nStart: ", data$Start,
                                   "\nEnd: ", data$End))
  
  fig = fig %>% add_lines(x = ~dfx,
                          y = ~dfy,
                          type = 'scatter',
                          mode = "lines",
                          name = ~dfName,
                          showlegend = F,
                          fill = 'tozeroy',
                          fillcolor = ~dfcolor,
                          opacity = .1,
                          line = list(
                            color = 'black',
                            width = 3),
                          marker = list(color = 'black', size = 1),
                          hoverinfo = "text",
                          hovertext = paste0(data$dfName, 
                                             "\nStart: ", data$dfStart,
                                             "\nEnd: ", data$dfEnd))
  
  
  fig = layout(fig,
               xaxis = list(title = goi$CHR, range = c(range.min, range.max)),
               yaxis = list(title = "Mean Lifespan Extension", range = c(-4,6)),
               shapes = list(
                 type = "rect",
                 line = list(color = "green", dash = 'dash', width = 2),
                 x0 = data[data$Symbol == gene.symbol, "Start"][1], 
                 x1 = data[data$Symbol == gene.symbol, "End"][1],
                 y0 = -5, y1 = 7))##,
  ##hovermode = "x unified")
  
  return(fig)
}

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
server <- function(input, output) {
  
  output$geneControls =  renderUI({
    
    selectizeInput("gene", label = "Select Gene", 
                   choices = GeneIDKey$Symbol,
                   options = list(create = TRUE,
                                  ##maxOptions = 5,
                                  placeholder = 'select a gene name'),
                   selected="Arc1")
    
  })  
   
  output$plot <- renderPlotly({


    if(input$locfactor == "gene"){ 
      p = plot.gene(gene.symbol = input$gene,
                start.range = input$start.range,
                end.range = input$end.range)
    }
    else{
      p = plot.pos(chr = input$CHR,
                start = input$start,
                end = input$end)
    }
    
    p %>% 
      layout(dragmode = "select") %>%
      event_register("plotly_selecting")
    
  })
  
  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Flybase link will appear here after click (double-click to clear)" else{
      page = a(paste0("Flybase Page for ", fbID$Symbol[grep(d$x,fbID$Position)]), href = paste0("http://flybase.org/reports/",fbID$FBgn[grep(d$x,fbID$Position)]), target="_blank")
      tagList(page)
    }
  })
  
}





# Run the application 
shinyApp(ui = ui, server = server)
