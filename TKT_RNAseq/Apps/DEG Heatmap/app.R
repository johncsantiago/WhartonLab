list.of.packages <- c("shiny", "ggfortify")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages('shiny')
  install.packages("ggfortify")
}

library(shiny)
library(ggfortify)
library(heatmaply)

git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/TKT_RNAseq/CountTables/"
cpmdata = read.csv(paste0(git.dir,"TKT_cpmdata.csv"),row.names = 1)
groups  = read.csv(paste0(git.dir,"TKT.metadata.csv"),row.names = 1)
TKT.EdgeR = read.csv(paste0(git.dir,"TKT.EdgeR.Table.csv"),row.names = 1)

convert=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgnConversionTable.csv",row.names=1)
temp=convert[,"Symbol"]
temp2=setdiff(row.names(cpmdata),convert[,2])
names(temp)=convert[,2]
names(temp2)=temp2
convert=c(temp,temp2)

comparison.groups = data.frame(GRF.CxDf  = c("^GR.F","Tkt..GR.F"),
                               GRF.CxOE  = c("^GR.F", "Tkt..GR.F"),
                               WTF.CxDf  = c("^WT.F", "Tkt..WT.F"),
                               WTF.CxOE  = c("^WT.F", "Tkt..WT.F"),
                               GRxWT.FC  = c("^GR.F", "^WT.F"),
                               GRxWT.FDf = c("TktDfGR.F", "TktDfWT.F"),
                               GRxWT.FOE = c("TktOEGR.F", "TktOEWT.F"),
                               GRM.CxDf  = c("^GR.M","Tkt..GR.M"),
                               GRM.CxOE  = c("^GR.M","Tkt..GR.M"),
                               WTM.CxDf  = c("^WT.M", "Tkt..WT.M"),
                               WTM.CxOE  = c("^WT.M", "Tkt..WT.M"),
                               GRxWT.MC  = c("^GR.M", "^WT.M"),
                               GRxWT.MDf = c("TktDfGR.M", "TktDfWT.M"),
                               GRxWT.MOE = c("TktOEGR.M", "TktOEWT.M"))

Titles = c("GR.FC x GR.FDf","GR.FC x GR.FOE",
           "WT.FC x WT.FDf", "WT.FC x WT.FOE",
           "GR.FC x WT.FC", "GR.FDf x WT.FDf", "GR.FOE x WT.FOE",
           "GR.MC x GR.MDf", "GR.MC x GR.MOE",
           "WT.MC x WT.MDf", "WT.MC x WT.MOE",
           "GR.MC x WT.MC", "GR.MDf x WT.MDf", "GR.MOE x WT.MOE")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("TKT RNAseq Experiment: DEG Heatmaps"),
  sidebarLayout(position="left", sidebarPanel(
    
    fluidRow(
      column(8,
             selectInput("comparison",
                       h5("Comparison"),
                       choices = list("GR.FC x GR.FDf"  = 1,
                                      "GR.FC x GR.FOE"  = 2,
                                      "WT.FC x WT.FDf"  = 3,
                                      "WT.FC x WT.FOE"  = 4,
                                      "GR.FC x WT.FC"   = 5,
                                      "GR.FDf x WT.FDf" = 6,
                                      "GR.FOE x WT.FOE" = 7,
                                      "GR.MC x GR.MDf"  = 8,
                                      "GR.MC x GR.MOE"  = 9,
                                      "WT.MC x WT.MDf"  = 10,
                                      "WT.MC x WT.MOE"  = 11,
                                      "GR.MC x WT.MC"   = 12,
                                      "GR.MDf x WT.MDf" = 13,
                                      "GR.MOE x WT.MOE" = 14),
                       selected = 1),
      )),
    
    fluidRow(
    column(6,
           numericInput("tophits",
                        h5("Number of DE Genes"), 
                        min = 10, 
                        max = 13582,
                        value = 25)
    ),),
  
  width = 3),
    
    
    mainPanel(
      plotlyOutput(outputId = "DEG.Heatmap", , height = 1000),
      ),
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  temp.data=reactive({
    
    use.samples =  comparison.groups[,as.numeric(input$comparison)]
    use.samples = colnames(cpmdata)[c(grep(use.samples[1],colnames(cpmdata)),
                                      grep(use.samples[2],colnames(cpmdata)))]
    use.cpm = cpmdata[,use.samples]
    
    use.FDR = TKT.EdgeR[,as.numeric(input$comparison)]
    names(use.FDR) = row.names(TKT.EdgeR)
    use.FDR = signif(use.FDR[order(use.FDR)],3)
    
    use.genes = use.FDR[1:as.numeric(input$tophits)]
    
    use.cpm = use.cpm[names(use.genes),]
    
    use.data = cbind(use.genes, use.cpm)
    
    return(use.data)
  })


  
output$DEG.Heatmap <- renderPlotly({
  use.data = temp.data()
  hmdata = as.matrix(use.data[,2:ncol(use.data)])
  
  ##col.text=1
  ##if(ncol(top.pca)>18){
    ##col.text=1-((ncol(top.pca)-18)*.011)
  ##}

  color.groups = groups[colnames(hmdata),1:3]
  
  hover.text = round(hmdata, digits = 1)
  
  FDR = hmdata
  FDR[,1:ncol(FDR)] = use.data[,1]
  
  FBgn = hmdata
  FBgn[,1:ncol(FBgn)] = row.names(FBgn)
  
  i=1
  while(i<ncol(hover.text)){
    hover.text[,i] = paste("CPM = ", hover.text[,i],
                           "\nFDR = ", FDR[,i],
                           "\nFBgn = ", FBgn[,i])
    i=i+1
  }

  heatmaply(hmdata,
            trace="none",
            col=RdYlBu(100)[100:1],
            scale="row",
            dendrogram = "both",
            show_dendrogram = c(F,T),
            Rowv=T,
            Colv=T,
            cexRow = 1/(as.numeric(input$tophits)/25),
            na.color="grey",
            labRow = convert[row.names(hmdata)],
            ##cexCol = col.text,
            cexCol = 1,
            key = T,
            ##margins = c(50,350,NA,0),
            ##row_dend_left = T,
            column_text_angle = 90,
            col_side_colors=color.groups,
            col_side_palette=Spectral,
            main=paste0("Top ",input$tophits, " DEGs for ", Titles[as.numeric(input$comparison)]),
            custom_hovertext = hover.text)
})
}


shinyApp(ui = ui, server = server, options)