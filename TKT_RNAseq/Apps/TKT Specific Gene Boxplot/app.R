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
library(DT)

git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/TKT_RNAseq/CountTables/"
cpmdata = read.csv(paste0(git.dir,"TKT_cpmdata.csv"),row.names = 1)
groups  = read.csv(paste0(git.dir,"TKT.metadata.csv"),row.names = 1)
TKT.EdgeR = read.csv(paste0(git.dir,"TKT.EdgeR.Table.csv"),row.names = 1)
convert=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgnConversionTable.csv",row.names=1)

groups$GxS = paste0(groups$Genotype, groups$Sex)
groups$GxTKT = paste0(groups$Genotype, groups$TKT)

groups.use = matrix(0,nrow = nrow(groups), ncol = 8)
row.names(groups.use) = row.names(groups)
colnames(groups.use) = c("ALL","GR", "WT", "M", "F", "C", "DF", "OE")
groups.use[                       ,      "ALL"] = 0
groups.use[groups$Genotype == "GR",      "GR"]  = 1
groups.use[groups$Genotype == "WT",      "WT"]  = 1
groups.use[groups$Sex      == "M",       "M"]   = 1
groups.use[groups$Sex      == "F",       "F"]   = 1
groups.use[groups$TKT      == "Control", "C"]   = 1
groups.use[groups$TKT      == "DF",      "DF"]  = 1
groups.use[groups$TKT      == "OE",      "OE"]  = 1

temp                   = convert[,c("updated.FBgn","Symbol")]
temp2                  = setdiff(row.names(cpmdata),row.names(convert))
temp2                  = matrix(rep(temp2,2), ncol=2)
row.names(temp2)       = temp2[,1]
colnames(temp2)        = colnames(temp)
FBgn2Symbol            = rbind(temp,temp2)
colnames(FBgn2Symbol)  = c("FBgn", "Symbol")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("TKT RNAseq Experiment: Custom PCA Figure"),
  sidebarLayout(position="left", sidebarPanel(
    
    br(),
    uiOutput("geneControls"),
    p("Select Gene: Use gene symbol if available. autocompletes"),
    ##br(), 
    
    fluidRow(
      column(6,
             selectInput("choice.order", h5("DEG Order"),
                         choices = list("GR.FC~GR.FDf"    = "GRF.CxDf",
                                        "GR.FC~GR.FOE"    = "GRF.CxOE",
                                        
                                        "WT.FC~WT.FDf"    = "WTF.CxDf",
                                        "WT.FC~WT.FOE"    = "WTF.CxOE",
                                        
                                        "GR.FC~WT.FC"     = "GRxWT.FC",
                                        "GR.FDf~WT.FDf"   = "GRxWT.FDf",
                                        "GR.FOE~WT.FOE"   = "GRxWT.FOE",
                                        
                                        "GR.MC~GR.MDf"    = "GRM.CxDf",
                                        "GR.MC~GR.MOE"    = "GRM.CxOE",
                                        
                                        "WT.MC~WT.MDf"    = "WTM.CxDf",
                                        "WT.MC~WT.MOE"    = "WTM.CxOE",
                                        
                                        "GR.MC~WT.MC"     = "GRxWT.MC",
                                        "GR.MDf~WT.MDf"   = "GRxWT.MDf",
                                        "GR.MOE~WT.MOE"   = "GRxWT.MOE",
                                        "GR.FC~GR.MC"     = "GR.FxM",
                                        "WT.FC~WT.MC"     = "WT.FxM"),
                         selected="GRF.CxDf"),
             p("Choose Signifigance Order: Organize the 'Select Gene' drop down menu order by FDR observed for a comparison between specific conditions"),
             br(),
      ),
      column(6,
             
             checkboxGroupInput("cond.select", h5("Select Conditions to Exclude"), 
                                choices = list("G85R"             = 2, 
                                               "Wild Type"        = 3,
                                               "Male"             = 4,
                                               "Female"           = 5,
                                               "TKT Control"      = 6,
                                               "TKT Df"           = 7,
                                               "TKT OE"           = 8)),
             ##p("Select Any Conditions to Exclude"),
             ##br(), 
      ),),
    
    fluidRow(
      column(6,
             selectInput("colfactor", h5("Select Factor for Color"), 
                         choices = list("Genotype x Sex" = "GxS",
                                        "Group"     = "Group",
                                        "Genotype"  = "Genotype",
                                        "TKT"       = "TKT", 
                                        "Sex"       = "Sex",
                                        "Genotype x TKT" = "GxTKT"),
                         selected = "Genotype x Sex"),
      ),
      column(6,
             
             selectInput("colselect", h5("Select Color"), 
                         choices = list("Choose Condition" = 0,
                                        "Color 1"            = 1, 
                                        "Color 2"             = 2,
                                        "Color 3"              = 3,
                                        "Color 4"           = 4),
                         selected = 0),
      ),),
    
    conditionalPanel(
      condition = "input.colselect == 1",
      colourInput("color1", "Color 1",
                  "brown",showColour = 'background')),
    
    
    conditionalPanel(
      condition = "input.colselect == 2",
      colourInput("color2", "Color 2",
                  "cornflowerblue",showColour = 'background')),
    
    conditionalPanel(
      condition = "input.colselect == 3",
      colourInput("color3", "Color 3",
                  "chartreuse4",showColour = 'background')),
    
    conditionalPanel(
      condition = "input.colselect == 4",
      colourInput("color4", "Color 4",
                  "darkgoldenrod1",showColour = 'background')),
    
    width = 3),
    
    mainPanel(
      plotlyOutput(outputId = "plot",
                   height=600, width=1300),
      br(),
      DTOutput('data'),
      tableOutput("data1"),
      p("Table of FDR values generated in the indicated comparison")

    ),
),)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  output$geneControls =  renderUI({
    gene.choices<<-FBgn2Symbol[row.names(TKT.EdgeR)[order(TKT.EdgeR[,input$choice.order])],"Symbol"]
    selectizeInput("gene", label = "Select Gene", 
                   choices = gene.choices,
                   options = list(create = TRUE,
                                  ##maxOptions = 5,
                                  placeholder = 'select a gene name'),
                   selected=gene.choices[1])
  })  

  output$plot <- renderPlotly({
    
    gene.name = FBgn2Symbol[FBgn2Symbol[,"Symbol"]==input$gene, "FBgn"]
    minlim    = min(as.numeric(na.omit(cpmdata[gene.name,])))
    maxlim    = max(as.numeric(na.omit(cpmdata[gene.name,])))
      
    if(length(input$cond.select) == 0){
      use.samples = row.names(groups.use)
    }    

    
    if(length(input$cond.select)>0){
      use.samples=row.names(groups.use)[apply(groups.use[,c(1,as.numeric(input$cond.select))],1,sum)==0]
    }
      
    boxdata=data.frame(cpm   = as.numeric(cpmdata[gene.name,use.samples]),
                       color = groups[use.samples,input$colfactor],
                       group = groups[use.samples,"Group"],
                       ID    = groups[use.samples, "SampleID"])
    
    row.names(boxdata) = boxdata$ID
    
    boxcolors = c(input$color1,input$color2,
                  input$color3,input$color4)
    
    if(length(unique(boxdata$color))<5){
      boxcolors = c(input$color1,input$color2,
                    input$color3,input$color4)[1:length(unique(boxdata$color))]      
    }
    
    Table.names = groups[boxdata$ID,]
    Table.names = Table.names[order(Table.names$TKT),]
    Table.names = Table.names[order(Table.names$Genotype),]
    Table.names = Table.names[order(Table.names$Sex),]

    boxdata = boxdata[row.names(Table.names),]
    boxdata$group = factor(boxdata$group, levels=(unique(boxdata$group)))
        
    fig = plot_ly(data = boxdata,
                  x = ~group,
                  y = ~cpm,
                  name = ~group,
                  type = "box",
                  boxpoints = "all",
                  pointpos = 0,
                  hoveron = 'points',
                  jitter = 0,
                  hovertext = paste0("Sample = ", boxdata$ID, "\nCPM = ", round(boxdata$cpm, 1)),
                  hoverinfo = 'text',
                  color = ~color,
                  colors = boxcolors,
                  marker = list(size = 15,
                                line = list(color = "black", width = 1.5)))
    
  })  

  output$data <- DT::renderDataTable({
    
    gene.name = FBgn2Symbol[FBgn2Symbol[,"Symbol"]==input$gene, "FBgn"]
    
    if(length(input$cond.select) == 0){
      use.samples = row.names(groups.use)
    }    
    
    
    if(length(input$cond.select)>0){
      use.samples=row.names(groups.use)[apply(groups.use[,c(1,as.numeric(input$cond.select))],1,sum)==0]
    }
    
    FDR = as.numeric(TKT.EdgeR[gene.name,])
    names(FDR) = colnames(TKT.EdgeR)
    FDR = signif(FDR,2)
    
    Table.names = groups[use.samples,]
    Table.names = Table.names[order(Table.names$TKT),]
    Table.names = Table.names[order(Table.names$Genotype),]
    Table.names = Table.names[order(Table.names$Sex),]
    
    FDR.table = data.frame(GR.F      = FDR[c(NA, "GRF.CxDf", "GRF.CxOE", "GRxWT.FC", NA, NA)],
                           TktDfGR.F = FDR[c("GRF.CxDf", NA, NA, NA, "GRxWT.FDf", NA)],
                           TktOEGR.F = FDR[c("GRF.CxOE", NA, NA, NA, "GRxWT.FOE", NA)],
                           WT.F      = FDR[c("GRxWT.FC", NA, NA, NA, "WTF.CxDf", "WTF.CxOE")],
                           TktDfWT.F = FDR[c(NA, "GRxWT.FDf", NA, "WTF.CxDf", NA, NA)],
                           TktOEWT.F = FDR[c(NA, NA, "GRxWT.FOE", "WTF.CxOE", NA, NA)],
                           GR.M      = FDR[c(NA, "GRF.CxDf", "GRF.CxOE", "GRxWT.FC", NA, NA)],
                           TktDfGR.M = FDR[c("GRF.CxDf", NA, NA, NA, "GRxWT.FDf", NA)],
                           TktOEGR.M = FDR[c("GRF.CxOE", NA, NA, NA, "GRxWT.FOE", NA)],
                           WT.M      = FDR[c("GRxWT.MC", NA, NA, NA, "WTM.CxDf", "WTM.CxOE")],
                           TktDfWT.M = FDR[c(NA, "GRxWT.MDf", NA, "WTM.CxDf", NA, NA)],
                           TktOEWT.M = FDR[c(NA, NA, "GRxWT.MOE", "WTM.CxOE", NA, NA)])
    
    FDR.table = FDR.table[,unique(Table.names$Group)]
    
    Table.names = groups
    Table.names = Table.names[order(Table.names$TKT),]
    Table.names = Table.names[order(Table.names$Genotype),]
    Table.names = Table.names[order(Table.names$Sex),]
    
    row.names(FDR.table) = unique(substr(Table.names$Group, 1, (nchar(Table.names$Group)-2)))
    
    df = FDR.table
    
    datatable(df, options = list(dom = 't', autoHideNavigation=T)) %>%
      formatStyle(colnames(df),color='black', fontSize = '100%') %>%
      
      formatStyle(na.omit(colnames(df)[df[1,]<.05]), 
                  backgroundColor = styleRow(c(1),c('yellow'))) %>% 
      formatStyle(na.omit(colnames(df)[df[2,]<.05]), 
                  backgroundColor = styleRow(c(2),c('yellow'))) %>% 
      formatStyle(na.omit(colnames(df)[df[3,]<.05]), 
                  backgroundColor = styleRow(c(3),c('yellow'))) %>% 
      formatStyle(na.omit(colnames(df)[df[4,]<.05]), 
                  backgroundColor = styleRow(c(4),c('yellow'))) %>% 
      formatStyle(na.omit(colnames(df)[df[5,]<.05]), 
                  backgroundColor = styleRow(c(5),c('yellow'))) %>% 
      formatStyle(na.omit(colnames(df)[df[6,]<.05]), 
                  backgroundColor = styleRow(c(6),c('yellow'))) %>% 
      
      formatStyle(colnames(df)[is.na(df[1,])], 
                  backgroundColor = styleRow(c(1),c('darkgrey'))) %>%
      formatStyle(colnames(df)[is.na(df[2,])], 
                  backgroundColor = styleRow(c(2),c('darkgrey'))) %>%
      formatStyle(colnames(df)[is.na(df[3,])], 
                  backgroundColor = styleRow(c(3),c('darkgrey'))) %>%
      formatStyle(colnames(df)[is.na(df[4,])], 
                  backgroundColor = styleRow(c(4),c('darkgrey'))) %>%
      formatStyle(colnames(df)[is.na(df[5,])], 
                  backgroundColor = styleRow(c(5),c('darkgrey'))) %>%
      formatStyle(colnames(df)[is.na(df[6,])], 
                  backgroundColor = styleRow(c(6),c('darkgrey'))) %>%
      
      formatStyle(colnames(df)[1:ncol(df)], border = styleRow(c(1:6),'2px solid black'))
  })

}


shinyApp(ui = ui, server = server, options)