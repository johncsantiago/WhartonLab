list.of.packages <- c("shiny", "ggfortify")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages('shiny')
  install.packages("ggfortify")
}

library(shiny)
library(ggfortify)
library(heatmaply)
library(plotly)
library(colourpicker)

git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/TKT_RNAseq/CountTables/"
cpmdata = read.csv(paste0(git.dir,"TKT_cpmdata.csv"),row.names = 1)
groups  = read.csv(paste0(git.dir,"TKT.metadata.csv"),row.names = 1)
convert=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgnConversionTable.csv",row.names=1)
temp=convert[,"Symbol"]
temp2=setdiff(row.names(cpmdata),convert[,2])
names(temp)=convert[,2]
names(temp2)=temp2
convert=c(temp,temp2)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("TKT RNAseq Experiment: Custom PCA Figure"),
  sidebarLayout(position="left", sidebarPanel(
    
    fluidRow(
      column(6,
    checkboxGroupInput("genotype",
                       h5("Genotype"),
                       choices = list("G85R" = "GR",
                                      "WT" = "WT"),
                       selected = ""),
      ),
    
    column(6,
    checkboxGroupInput("sex",
                       h5("Sex"),
                       choices = list("Male"   = "M",
                                      "Female" = "F"),
                       selected = ""),
    ),),
    
    
    fluidRow(
      column(6,
    checkboxGroupInput("tkt",
                       h5("TKT"),
                       choices = list("Control" = "Control",
                                      "Deficient"    = "DF",
                                      "Over Expression" = "OE"),
                       selected = ""),
      ),
    
    column(6,
           radioButtons("pc",
                              h5("PC"),
                              choices = list("PC1/2" = 1,
                                             "PC1/3" = 2,
                                             "PC2/3" = 3),
                              selected = 1),
    ),),
    
    
    fluidRow(
      column(6,
             selectInput("color", h5("Color Variable"),
                         choices = list("Genotype"   = "Genotype",
                                        "Sex"       = "Sex",
                                        "TKT" = "TKT"),
                         selected="TKT"),
    
      ),
    
    
    column(6,
           selectInput("shape", h5("Shape Variable"),
                       choices = list("Genotype"   = "Genotype",
                                      "Sex"       = "Sex",
                                      "TKT" = "TKT"),
                       selected="Genotype"),
    ),),
    
    
    fluidRow(
      column(6,
    sliderInput("size", h5("Size"), min = 0, max = 50, value = 15, step=1),
      ),
    
    column(6,
           sliderInput("tophits", h5("Number of PC Genes"), min = 10, max = 50, value = 25, step=5)
    ),),
    
    selectInput("colselect", h5("Select Color"), 
                choices = list("Choose Condition" = 0,
                               "Color 1"            = 1, 
                               "Color 2"             = 2,
                               "Color 3"              = 3),
                selected = 0),
  
  conditionalPanel(
    condition = "input.colselect == 1",
    colourInput("color1", "Color 1",
                "brown",showColour = 'background')),
  
  
  conditionalPanel(
    condition = "input.colselect == 2",
    colourInput("color2", "Color 2",
                "lightskyblue",showColour = 'background')),
  
  conditionalPanel(
    condition = "input.colselect == 3",
    colourInput("color3", "Color 3",
                "darkgoldenrod1",showColour = 'background')),
  
  width = 3),
    
    
    mainPanel(
      plotlyOutput(outputId = "distPlot"),
      fluidRow(
        column(6,
      plotlyOutput(outputId = "topPC1hits", height = 550),
        ),
      column(6,
      plotlyOutput(outputId = "topPC2hits", height = 550))
      ),),
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  temp.groups=reactive({
    use.groups =  groups
    
    if(length(input$genotype) == 1){
      use.groups = use.groups[use.groups$Genotype == input$genotype,]
    }
    
    if(length(input$sex) == 1){
      use.groups = use.groups[use.groups$Sex == input$sex,]
    }
    
    if(length(input$tkt) == 1){
      use.groups = use.groups[use.groups$TKT == input$tkt,]
    }
    
    if(length(input$tkt) == 2){
      use.groups = use.groups[use.groups$TKT == input$tkt[1]|
                                use.groups$TKT == input$tkt[2],]
    }

    return(use.groups)
  })
  
  
  newdata=reactive({
    use.groups=temp.groups()
    use.data   = cpmdata[,row.names(use.groups)]
    use.data   = use.data[apply(use.data,1,sum)>0,]
    pca <- prcomp(t(use.data), scale.=TRUE) 
    gr <- as.data.frame(use.groups)
    PC1=scale(pca$x[,1])
    PC2=scale(pca$x[,2])
    eigs <- pca$sdev^2
    ve=signif(((eigs / sum(eigs))*100)[1:3],4)
    pcs=matrix(c(1,2,1,3,2,3),ncol =3)
    pcs=pcs[,as.numeric(input$pc)]
    
      
    pca.data=data.frame(Sample=row.names(gr),
                        ##PC1=scale(pca$x[,1]),
                        ##PC2=scale(pca$x[,2]),
                        PC1=scale(pca$x[,pcs[1]]),
                        PC2=scale(pca$x[,pcs[2]]),
                        sdev=pca$sdev,
                        Genotype=gr$Genotype, 
                        Sex=gr$Sex, 
                        TKT=gr$TKT,
                        Percent1=ve[pcs[1]],
                        Percent2=ve[pcs[2]],
                        Group=gr$Group)
    
    return(pca.data)
  })

output$distPlot <- renderPlotly({
  pca.data=newdata()
  pca.shape = factor(pca.data[,input$shape])
  pca.color = factor(pca.data[,input$color])
  percent1 = pca.data$Percent1[1]
  percent2 = pca.data$Percent2[1]
  pcs=matrix(c(1,2,1,3,2,3),ncol =3)
  pcs=pcs[,as.numeric(input$pc)]
  
  fig  = plot_ly(data = pca.data, 
                 x = ~PC1, 
                 y = ~PC2,
                 name = ~Group,
                 type = 'scatter',
                 mode = 'markers',
                 showlegend = T,
                 color = pca.color,
                 symbol = pca.shape,
                 ##symbols = as.character(seq_along(unique(point.data$Set))),
                 colors = c(input$color1, input$color2, input$color3)[c(1:length(levels(pca.color)))],
                 marker = list(size = input$size,
                               line = list(color = "black", width = 1.5)),
                 hoverinfo = "text",
                 hovertext = paste("Sample:", pca.data$Sample)) %>%
    layout(xaxis = list(title = paste0("PC", pcs[1], "(", percent1, "%)")),
           yaxis = list(title = paste0("PC", pcs[2], "(", percent2, "%)")))
  
  fig
  
    })
 




output$topPC1hits <- renderPlotly({
  use.groups = temp.groups()
  use.data   = cpmdata[,row.names(use.groups)]
  use.data   = use.data[apply(use.data,1,sum)>0,]
  pca.data <- prcomp(t(use.data), scale.=TRUE) 
  pcs=matrix(c(1,2,1,3,2,3),ncol =3)
  pcs=pcs[,as.numeric(input$pc)]
  PCcontribute=pca.data$rotation[,pcs[1]]
  PCcontribute=PCcontribute[order(PCcontribute,decreasing=T)]
  
  ##The column sum of squares of the loadings (pca$rotation) are the variances of PCs.
  #the squared factor loading is the percent of variance in that variable explained by the factor
  percents=(PCcontribute*PCcontribute)

  ##puts the genes in decreasing order
  percents=percents[order(percents,decreasing=T)]
  pca.genes=PCcontribute[names(percents)[1:input$tophits]]

  pca.genes=names(pca.genes[order(pca.genes)])
  PCorder=pca.data$x[,pcs[1]]
  PCorder=PCorder[order(PCorder)]
  
  top.pca=cpmdata[pca.genes,names(PCorder)]
  row.names(top.pca)=
    convert[row.names(top.pca)]
  color.groups=groups[names(PCorder),]

  
  col.text=1
  if(ncol(top.pca)>18){
    col.text=1-((ncol(top.pca)-18)*.011)
  }

  cpm.text = top.pca
  i=1
  while(i<ncol(cpm.text)){
    cpm.text[,i] = paste("CPM = ", round(cpm.text[,i], digits = 1))
    i=i+1
  }
  
  heatmaply(top.pca,trace="none",col=RdYlBu(100)[100:1], scale="row",
            dendrogram = "none",Rowv=F,Colv=F,
            cexRow = .75, na.color="grey",
            labRow = row.names(top.pca),
            cexCol = col.text,
            key = T,
            column_text_angle = 90,
            col_side_colors=color.groups[,1:4],
            col_side_palette=Spectral,
            main=paste0("PC",pcs[1]),
            custom_hovertext = cpm.text)
})

output$topPC2hits <- renderPlotly({
  use.groups = temp.groups()
  use.data   = cpmdata[,row.names(use.groups)]
  use.data   = use.data[apply(use.data,1,sum)>0,]
  pcs=matrix(c(1,2,1,3,2,3),ncol =3)
  pcs=pcs[,as.numeric(input$pc)]
  
  pca.data <- prcomp(t(use.data), scale.=TRUE) 
  PCcontribute=pca.data$rotation[,pcs[2]]
  
  PCcontribute=PCcontribute[order(PCcontribute,decreasing=T)]
  
  ##The column sum of squares of the loadings (pca$rotation) are the variances of PCs.
  #the squared factor loading is the percent of variance in that variable explained by the factor
  percents=(PCcontribute*PCcontribute)
  
  ##puts the genes in decreasing order
  percents=percents[order(percents,decreasing=T)]
  pca.genes=PCcontribute[names(percents)[1:input$tophits]]
  
  pca.genes=names(pca.genes[order(pca.genes)])
  PCorder=pca.data$x[,pcs[2]]
  PCorder=PCorder[order(PCorder)]
  
  top.pca=cpmdata[pca.genes,names(PCorder)]
  row.names(top.pca)=
    convert[row.names(top.pca)]
  color.groups=groups[names(PCorder),]

  col.text=1
  if(ncol(top.pca)>18){
    col.text=1-((ncol(top.pca)-18)*.011)
  }
  
  cpm.text = top.pca
  i=1
  while(i<ncol(cpm.text)){
    cpm.text[,i] = paste("CPM = ", round(cpm.text[,i], digits = 1))
    i=i+1
  }
  
  heatmaply(top.pca,trace="none",col=RdYlBu(100)[100:1], scale="row",
            dendrogram = "none",Rowv=F,Colv=F,
            cexRow = .75, na.color="grey",
            cexCol = col.text,
            labRow = row.names(top.pca),
            key = T,
            column_text_angle = 90,
            col_side_colors=color.groups[,1:4],
            col_side_palette=Spectral,
            main=paste0("PC",pcs[2]),
            custom_hovertext = cpm.text)
  
})

}


shinyApp(ui = ui, server = server, options)