##library(org.Dm.eg.db)

##gene.start = abs(unlist(as.list(org.Dm.egCHRLOC)))
##gene.end   = abs(unlist(as.list(org.Dm.egCHRLOCEND)))

##gene.dir = unlist(as.list(org.Dm.egCHRLOC))/abs(unlist(as.list(org.Dm.egCHRLOC)))

##gene.info = data.frame(Start = gene.start, End = gene.end, Strand = gene.dir, CHR = "", Name = "")

##i = 1
##while(i<=length(gene.start)){
  ##temp = strsplit(names(gene.start)[i], "\\.")
  ##gene.info[i, c("Name", "CHR")] = temp[[1]]
  ##i = i + 1
##}

##i = 1
##while(i<=length(unique(gene.info$Name))){
  ##gene.info[gene.info$Name %in% unique(gene.info$Name)[i],"Start"] =
  ##min(gene.info[gene.info$Name %in% unique(gene.info$Name)[i],"Start"])
  ##gene.info[gene.info$Name %in% unique(gene.info$Name)[i],"End"] =
  ##  max(gene.info[gene.info$Name %in% unique(gene.info$Name)[i],"End"])
  ##if(length(unique(gene.info[gene.info$Name %in% unique(gene.info$Name)[i],"Strand"])) == 2){
  ##  gene.info[gene.info$Name %in% unique(gene.info$Name)[i],"Strand"] = 0
  ##}
  ##i = i + 1
##}

##gene.info = unique(gene.info)

##use.genes = setdiff(unique(gene.info$Name), gene.info[duplicated(gene.info$Name),"Name"])
##gene.info = gene.info[gene.info$Name %in% use.genes, ]
##row.names(gene.info) = gene.info$Name
##git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/"
##GeneIDKey = read.csv(paste0(git.dir,"GeneralDataFiles/GeneIDKey.csv"),row.names = 1)
##gene.info = gene.info[GeneIDKey$ensembl,]
##row.names(gene.info) = row.names(GeneIDKey)
##GeneIDKey = cbind(GeneIDKey, gene.info)
##write.csv(GeneIDKey, "/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/GeneIDKey.csv")
    


Df2gene = vector("list", nrow(DfKey))
names(Df2gene) = DfKey$Deficiency


i=1
while(i <= length(Df2gene)){
  
  same.chr = setNames(GeneIDKey[grep(DfKey$Chr[i], GeneIDKey$CHR), "Start"], GeneIDKey[grep(DfKey$Chr[i], GeneIDKey$CHR), "FBgn"])
  start.inDf = same.chr[same.chr>DfKey$Start[i] & same.chr<DfKey$End[i]] 
  
  same.chr = setNames(GeneIDKey[grep(DfKey$Chr[i], GeneIDKey$CHR), "End"], GeneIDKey[grep(DfKey$Chr[i], GeneIDKey$CHR), "FBgn"])
  end.inDf = same.chr[same.chr>DfKey$Start[i] & same.chr<DfKey$End[i]] 
  
  Df2gene[[i]] = unique(c(names(start.inDf),names(end.inDf)))
  i = i + 1
}


library(plotly)
git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/"

DfKey = read.csv(paste0(git.dir, "DeficiencyScreen/DataFiles/DfKey.csv"), row.names = 1)
row.names(DfKey) = DfKey$Deficiency
DfKey[DfKey$MLE == 0, "MLE"] = .05
GeneIDKey = read.csv(paste0(git.dir,"GeneralDataFiles/GeneIDKey.csv"),row.names = 1)



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

plot.pos(chr = "2R",
         start = 5600000,
         end = 5800000)





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
               xaxis = list(title = goi$CHR, range = c(range.min, range.max)),
               yaxis = list(title = "Mean Lifespan Extension", range = c(-4,6)))##,
  ##hovermode = "x unified")
  
  return(fig)
}

plot.gene(gene.symbol = "Arc1",
          start.range = 57000,
          end.range = 58000)


fig = fig %>% add_trace(x = ~shapex,
                        y = ~shapey,
                        type = 'scatter',
                        mode = 'markers',
                        name = ~shape.name,
                        showlegend = F,
                        colors =  ~shape.color,
                        ##fillcolor = ~dfcolor,
                        ##opacity = .1,
                        marker = list(color = data$shape.color, symbol = data$shape, size = 5),
                        hoverinfo = "text",
                        hovertext = paste0(data$shape.name),
                        inherit = F)




fig = fig %>% add_lines(x = ~x,
                        y = ~y,
                        type = 'scatter',
                        mode = "lines",
                        name = ~Symbol,
                        showlegend = F,
                        fill = 'none',
                        opacity = 1,
                        colors = data$color,
                        line = list(
                          color = ~color
                        ),
                        hoverinfo = "text",
                        hovertext = paste0(data$Symbol,
                                           "\nStart: ", data$Start,
                                           "\nEnd: ", data$End))







gene.fb = goi[, "FBgn"]
gene.ens = goi[, "ensembl"]

gene.Dfs = names(Df2gene[grep(gene.fb, Df2gene)])

X = DfKey[grep("X", DfKey$Chr),]

X.xpos = matrix(0,
                nrow = 6,
                ncol = nrow(X))

colnames(X.xpos) = row.names(X)
X.xpos[2,] = X$Start
X.xpos[3,] = X$Start

X.xpos[4,] = X$End
X.xpos[5,] = X$End

X.xpos[6,] = max(X$End)

X.ypos = matrix(0,
                nrow = 6,
                ncol = nrow(X))

colnames(X.ypos) = row.names(X)
X.ypos[3,] = X$MLE
X.ypos[4,] = X$MLE

data = data.frame(x = X.xpos[,1], y = X.ypos[,1])
data$df = row.names(X)[1]

i=2
while(i<=nrow(X)){
  temp = data.frame(x = X.xpos[,i], y = X.ypos[,i])
  temp$df = row.names(X)[i]
  data = rbind(data,temp)
  i=i+1
}
data$x = as.numeric(data$x)
data$y = as.numeric(data$y)

data$color = (DfKey[data$df, "Modifier"])

data$color[data$color == "Strong Suppressor"] = adjustcolor("brown", alpha.f = 0.2)
data$color[data$color == "Suppressor"] = adjustcolor("gold", alpha.f = 0.2)
data$color[data$color == "No Effect"] = adjustcolor("lightgrey", alpha.f = 0.2)
data$color[data$color == "Enhancer"] = adjustcolor("steelblue", alpha.f = 0.2)

data = as.matrix(data2)
data = as.data.frame(data)
fig =plot_ly(data,
             x = ~x,
             y = ~y,
             type = 'scatter',
             mode = "lines",
             name = ~df,
             showlegend = F,
             fill = 'tozeroy',
             hoveron = 'points+fills',
             layer = "below",
             fillcolor = ~color,
             opacity = .1,
             line = list(
               color = 'black'
             ),
             hoverinfo = "text",
             hovertext = paste("Deficiency:", data$df))

fig


fig = layout(fig,
             shapes = list(
               type = "rect",
               fillcolor = "red", line = list(color = "red"), opacity = 0.3,
               x0 = X.xpos[i,2], x1 = X.xpos[i,4],
               y0 = data$y[2], y1 = data$y[4]
             ))


fig


fig = plot_ly(data,
              x = ~x,
              y = ~y,
              type = 'scatter',
              mode = "lines + markers",
              name = ~Symbol,
              hoveron = 'points+lines',
              showlegend = F,
              line = list(
                color = data$color
              ),
              marker = list(color = data$color, symbol = ~shape, size = 10),
              hoverinfo = "text",
              hovertext = paste0(data$Symbol, 
                                 "\nStart: ", data$Start,
                                 "\nEnd: ", data$End))




library(plotly)
library(crosstalk)
library(DT)


sd <- SharedData$new(iris)

a <- plot_ly(sd, x = ~Sepal.Width, y = ~Petal.Width) %>% 
  add_markers(alpha = 0.5) %>%
  highlight("plotly_selected", dynamic = TRUE)


options(persistent = TRUE)

p <- datatable(sd)

bscols(widths = c(6, 4), a, p)


plotly_example("shiny", "event_data")

library(shiny)
library(plotly)

ui <- fluidPage(
  radioButtons("plotType", "Plot Type:", choices = c("ggplotly", "plotly")),
  plotlyOutput("plot"),
  verbatimTextOutput("hover"),
  verbatimTextOutput("click"),
  verbatimTextOutput("brushing"),
  verbatimTextOutput("selecting"),
  verbatimTextOutput("brushed"),
  verbatimTextOutput("selected")
)

server <- function(input, output, session) {
  
  nms <- row.names(mtcars)
  
  output$plot <- renderPlotly({
    p <- if (identical(input$plotType, "ggplotly")) {
      ggplotly(ggplot(mtcars, aes(x = mpg, y = wt, customdata = nms)) + geom_point())
    } else {
      plot_ly(mtcars, x = ~mpg, y = ~wt, customdata = nms)
    }
    p %>% 
      layout(dragmode = "select") %>%
      event_register("plotly_selecting")
  })
  
  output$hover <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  })
  
  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Click events appear here (double-click to clear)" else d
  })
  
  output$brushing <- renderPrint({
    d <- event_data("plotly_brushing")
    if (is.null(d)) "Brush extents appear here (double-click to clear)" else d
  })
  
  output$selecting <- renderPrint({
    d <- event_data("plotly_selecting")
    if (is.null(d)) "Brush points appear here (double-click to clear)" else d
  })
  
  output$brushed <- renderPrint({
    d <- event_data("plotly_brushed")
    if (is.null(d)) "Brush extents appear here (double-click to clear)" else d
  })
  
  output$selected <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Brushed points appear here (double-click to clear)" else d
  })
  
}

shinyApp(ui, server, options = list(display.mode = "showcase"))







plotly_example("shiny", "event_data_parcoords")
