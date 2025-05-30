list.of.packages <- c("shiny", "ggfortify")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages('shiny')
  install.packages("ggfortify")
}

##Load necessary libraries
library(shiny)
library(ggfortify)
library(heatmaply)
library(plotly)
library(colourpicker)
library(DT)
library(rclipboard)

##load all files from github
git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/master/TKT_RNAseq/CountTables/"
cpmdata = read.csv(paste0(git.dir,"TKT_cpmdata.csv"),row.names = 1)
groups  = read.csv(paste0(git.dir,"TKT.metadata.csv"),row.names = 1)
TKT.EdgeR = read.csv(paste0(git.dir,"TKT.EdgeR.Table.csv"),row.names = 1)
convert=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgnConversionTable.csv",row.names=1)

##Add additional groups for color factors
groups$GxS = paste0(groups$Genotype, groups$Sex)
groups$GxTKT = paste0(groups$Genotype, groups$TKT)

##FBgn2Symbol; Variable for converting FBgn to symbol
temp                   = convert[,c("updated.FBgn","Symbol")]
temp2                  = setdiff(row.names(cpmdata),row.names(convert))
temp2                  = matrix(rep(temp2,2), ncol=2)
row.names(temp2)       = temp2[,1]
colnames(temp2)        = colnames(temp)
FBgn2Symbol            = rbind(temp,temp2)
colnames(FBgn2Symbol)  = c("FBgn", "Symbol")

##groups. use; Variable used by the "Filter gene set" options 
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

##all.genes; Variable used to copy gene set to clipboard so each gene is on a new line
temp=FBgn2Symbol[,"Symbol"]
i=2
temp2=paste0(temp[1],"\n")
while(i<=length(temp)){
  temp2=paste0(temp2,paste0(temp[i],"\n"))
  i=i+1
}
temp2
all.genes = temp2

##mean.cpm; Mean cpm for each gene in each condition (groups$Group). Used for filtering
mean.cpm = matrix(0, ncol = length(unique(groups$Group)), nrow = nrow(cpmdata))
colnames(mean.cpm) = unique(groups$Group)
row.names(mean.cpm) = row.names(cpmdata)
i=1
while(i<=ncol(mean.cpm)){
  mean.cpm[,i]= apply(cpmdata[,row.names(groups[groups$Group == colnames(mean.cpm)[i],])], 1, mean)
  i=i+1
}

##filter1-4.genes; Variable for filter conditions, set to no filter as default
filter1.genes = row.names(TKT.EdgeR)
filter2.genes = row.names(TKT.EdgeR)
filter3.genes = row.names(TKT.EdgeR)
filter4.genes = row.names(TKT.EdgeR)


##Starts the code for user interface. Sidebar controls and Main panel display
ui <- fluidPage(
  
##Sidebar control code 
  titlePanel("TKT RNAseq Experiment: Specific Gene Boxplot"),
  sidebarLayout(position="left", sidebarPanel(
    
##General plot customization tools 
    ##gene.Controls; Gene selection tool. Uses a reactive variable
    uiOutput("geneControls"),
    p("Select Gene: Use gene symbol if available. autocompletes"),
    
    ##url; Link to selected gene flybase summary page. Uses a reactive variable
    uiOutput("url"),
    br(),
    

    ##choice. order; A drop down tab used to select the EdgeR comparison of interest when sorting the gene list used in "geneControls" 
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
             p("Choose Signifigance Order: Organize the 'Select Gene' drop down menu order by FDR observed for a comparison between specific conditions"),),
      
    ##cond.select; A a checkbox selection with multiple selections allowed. Used to select conditions to exclude when creating the boxplot
      column(6,
             checkboxGroupInput("cond.select", h5("Select Conditions to Exclude"), 
                                choices = list("G85R"             = 2, 
                                               "Wild Type"        = 3,
                                               "Male"             = 4,
                                               "Female"           = 5,
                                               "TKT Control"      = 6,
                                               "TKT Df"           = 7,
                                               "TKT OE"           = 8)),
             ),
    ),
    
##This section is used to customize the colors of the box plots 
  fluidRow(
    ##colfactor; A drop down tab used to select the factor to differentiate by color
      column(6,
             selectInput("colfactor", h5("Select Factor for Color"), 
                         choices = list("Genotype x Sex" = "GxS",
                                        "Group"     = "Group",
                                        "Genotype"  = "Genotype",
                                        "TKT"       = "TKT", 
                                        "Sex"       = "Sex",
                                        "Genotype x TKT" = "GxTKT"),
                         selected = "Genotype x Sex"),),
      
    ##colselect; A drop down tab that lets you select a color to change. Can only customize up to 4 factor colors. Once selected, a new tab will open that can be selected which will open up a color pallette to choose from.
      column(6,
             
             selectInput("colselect", h5("Select Color"), 
                         choices = list("Choose Condition" = 0,
                                        "Color 1"          = 1, 
                                        "Color 2"          = 2,
                                        "Color 3"          = 3,
                                        "Color 4"          = 4),
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
    
br(),

##This section is used to filter the gene set used in geneControls 
  fluidRow(
    ##Header for the "Filter Genes Set" section
    column(12,
      h4("Filter Genes Set")
    ),),
    
    
  fluidRow(
  ##This section will be used to filter the gene set used in geneControls by EdgeR significance
    ##fdr.filter1; single selection buttons. Used to select the choice of significance used with deg.condition1 
      column(6,
           radioButtons("fdr.filter1",
                        h5(""),
                        choices = list("Significant in"      = 2,
                                       "Not Significant in"  = 3),
                        selected = 2),),

    ##deg.condtion1; drop down menu used to select the EdgeR comparison that will be use fdr.filter1 to filter genes by significance 
      column(6,
           selectInput("deg.condition1", h5(""),
                       choices = list("No Filter"       = "none",
                                      "GR.FC~GR.FDf"    = "GRF.CxDf",
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
                       selected="No Filter"),
           p(""),),),

  fluidRow(
  ##fdr.filter2; single selection buttons. Used to select the choice of significance used with deg.condition2 
    column(6,
           radioButtons("fdr.filter2",
                        h5(""),
                        choices = list("Significant in"      = 2,
                                       "Not Significant in"  = 3),
                        selected = 2),),
    
    ##deg.condtion2; drop down menu used to select the EdgeR comparison that will be use fdr.filter2 to filter genes by significance 
    column(6,
           selectInput("deg.condition2", h5(""),
                       choices = list("No Filter"       = "none",
                                      "GR.FC~GR.FDf"    = "GRF.CxDf",
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
                       selected="No Filter"),
           p(""),),),

  fluidRow(
  ##fdr.filter3; single selection buttons. Used to select the choice of significance used with deg.condition3 
    column(6,
         radioButtons("fdr.filter3",
                      h5(""),
                      choices = list("Significant in"      = 2,
                                     "Not Significant in"  = 3),
                      selected = 2),),
  
  ##deg.condtion3; drop down menu used to select the EdgeR comparison that will be use fdr.filter3 to filter genes by significance 
    column(6,
         selectInput("deg.condition3", h5(""),
                     choices = list("No Filter"       = "none",
                                    "GR.FC~GR.FDf"    = "GRF.CxDf",
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
                     selected="No Filter"),
         p(""),),),

  fluidRow(
  ##Filter gene set by a minimal expression level
  ##min.value; A numeric input used with min.cond to filter genes with a minimal mean cpm cutoff using mean.cpm
    column(6,
         numericInput("min.value",
                      h5("Minimum CPM level"),
                      min = 0,
                      max = 1000,
                      value = 0),),
  
  ##min.cond; Drop down box used to choose the experimental condition (groups$Group) that will be used with min.value to filter genes
    br(),
    column(6,
         selectInput("min.cond", h5(""),
                     choices = list("Any Condition" = "all",
                                    "GR.FC"         = "GR.F",
                                    "TktDfGR.F"     = "TktDfGR.F",
                                    "TktOEGR.F"     = "TktOEGR.F",
                                    
                                    "WT.FC"         = "WT.F",
                                    "TktDfWT.F"     = "TktDfWT.F",
                                    "TktOEWT.F"     = "TktOEWT.F",
                                    
                                    "GR.MC"         = "GR.M",
                                    "TktDfGR.M"     = "TktDfGR.M",
                                    "TktOEGR.M"     = "TktOEGR.M",
                                    
                                    "WT.MC"         = "WT.M",
                                    "TktDfWT.M"     = "TktDfWT.M",
                                    "TktOEWT.M"     = "TktOEWT.M"),
                     selected="Any Condition"),
         p(""),),),
  br(),

##This section is used to allow functional analysis of the filtered gene sets 
  fluidRow(
  ##Header for the "Gene Function Analysis" section
    column(12,
           h4("Gene Function Analysis")
    ),),
  
  ##GOrilla; Link to GOrilla home page
    uiOutput("GOrilla"),
  
  ##Pangea; Link to PANGEA home page
    uiOutput("Pangea"),

  ##Provides a button to copy the filtered gene set to the user clipboard
    rclipboardSetup(),
    ##clip; reactive function that stores the list of filtered genes
      uiOutput("clip"),

##Provides a button to copy the full gene set to the user clipboard  
    rclipboardSetup(),
    ##clip2; stores the full list of genes
      uiOutput("clip2"),
  
    width = 3),

##Figure display settings    
  mainPanel(
    ##Display the boxplot
    plotlyOutput(outputId = "plot",
                 height=600, width=1250),
    br(),
    
    br(),
    
    ##Display the FDR table
    DTOutput('data', width =1100),
    p("Table of FDR values generated in the indicated comparison.")
    ),),
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  output$geneControls =  renderUI({
    
    if(input$deg.condition1 == "none"){
      filter1.genes <<- row.names(TKT.EdgeR)
    } 
    
    if(input$fdr.filter1 == 2 & input$deg.condition1 != "none"){
      filter1.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition1] < .05,])
    } 
    if(input$fdr.filter1 == 3 & input$deg.condition1 != "none"){
      filter1.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition1] >= .05,])
    } 
      
      
    if(input$deg.condition2 == "none"){
      filter2.genes <<- row.names(TKT.EdgeR)
    }
    
    if(input$fdr.filter2 == 2 & input$deg.condition2 != "none"){
      filter2.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition2] < .05,])
    }
    if(input$fdr.filter2 == 3 & input$deg.condition2 != "none"){
      filter2.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition2] >= .05,])
    } 
    
    if(input$deg.condition3 == "none"){
      filter3.genes <<- row.names(TKT.EdgeR)
    }
    
    if(input$fdr.filter3 == 2 & input$deg.condition3 != "none"){
      filter3.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition3] < .05,])
    }
    if(input$fdr.filter3 == 3 & input$deg.condition3 != "none"){
      filter3.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition3] >= .05,])
    } 
    
    if(input$min.cond == "all"){
      filter4.genes = row.names(mean.cpm)[apply(mean.cpm,1,min)>input$min.value]
    }
    
    if(input$min.cond != "all"){
      filter4.genes = row.names(mean.cpm)[mean.cpm[,input$min.cond] > input$min.value]
    }
    
    gene.choices   <<- intersect(row.names(TKT.EdgeR)[order(TKT.EdgeR[,input$choice.order])],
                                 intersect(intersect(filter1.genes, filter2.genes),
                                           intersect(filter3.genes,filter4.genes)))
    
    gene.choices <<- as.character(FBgn2Symbol[gene.choices,"Symbol"])
    
    
    selectizeInput("gene", label = "Select Gene", 
                   choices = gene.choices,
                   options = list(create = TRUE,
                                  ##maxOptions = 5,
                                  placeholder = 'select a gene name'),
                   selected="CG8036")

  })  

  
  output$clip <- renderUI({
    output$clip <- renderUI({
     
      
      if(input$deg.condition1 == "none"){
        filter1.genes <<- row.names(TKT.EdgeR)
      } 
      
      if(input$fdr.filter1 == 2 & input$deg.condition1 != "none"){
        filter1.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition1] < .05,])
      } 
      if(input$fdr.filter1 == 3 & input$deg.condition1 != "none"){
        filter1.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition1] >= .05,])
      } 
      
      
      if(input$deg.condition2 == "none"){
        filter2.genes <<- row.names(TKT.EdgeR)
      }
      
      if(input$fdr.filter2 == 2 & input$deg.condition2 != "none"){
        filter2.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition2] < .05,])
      }
      if(input$fdr.filter2 == 3 & input$deg.condition2 != "none"){
        filter2.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition2] >= .05,])
      } 
      
      if(input$deg.condition3 == "none"){
        filter3.genes <<- row.names(TKT.EdgeR)
      }
      
      if(input$fdr.filter3 == 2 & input$deg.condition3 != "none"){
        filter3.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition3] < .05,])
      }
      if(input$fdr.filter3 == 3 & input$deg.condition3 != "none"){
        filter3.genes <<- row.names(TKT.EdgeR[TKT.EdgeR[,input$deg.condition3] >= .05,])
      } 
      
      if(input$min.cond == "all"){
        filter4.genes = row.names(mean.cpm)[apply(mean.cpm,1,min)>input$min.value]
      }
      
      if(input$min.cond != "all"){
        filter4.genes = row.names(mean.cpm)[mean.cpm[,input$min.cond] > input$min.value]
      }
      
      gene.choices   <<- intersect(row.names(TKT.EdgeR)[order(TKT.EdgeR[,input$choice.order])],
                                   intersect(intersect(filter1.genes, filter2.genes),
                                             intersect(filter3.genes,filter4.genes)))
      
      gene.choices <<- as.character(FBgn2Symbol[gene.choices,"Symbol"])
      
      temp=gene.choices
      i=2
      temp2=paste0(temp[1],"\n")
      while(i<=length(temp)){
        temp2=paste0(temp2,paste0(temp[i],"\n"))
        i=i+1
      }
      temp2
      
      gene.choices <<- temp2
      
      rclipButton(
        inputId = "clipbtn",
        label = "Copy Target Geneset to Clipboard",
        clipText = gene.choices, 
        icon = icon("clipboard")
      )
    })
  })
  
  output$clip2 <- renderUI({
    output$clip2 <- renderUI({
      rclipButton(
        inputId = "clipbtn",
        label = "Copy Background Geneset to Clipboard",
        clipText = all.genes, 
        icon = icon("clipboard")
      )
    })
  })
  
  
    output$url <- renderUI({
      page = a(paste0("Flybase Page for ", input$gene), href = paste0("http://flybase.org/reports/",FBgn2Symbol[FBgn2Symbol[,"Symbol"]==input$gene,"FBgn"]), target="_blank")
      tagList(page)
    })
    
    output$GOrilla <- renderUI({
      page = a("GOrilla Homepage", href = "https://cbl-gorilla.cs.technion.ac.il", target="_blank")
      tagList(page)
    })
    
    output$Pangea <- renderUI({
      page = a("PANGEA Homepage", href = " https://www.flyrnai.org/tools/pangea/web/home/7227", target="_blank")
      tagList(page)
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
    FDR = signif(FDR,3)
    
    Table.names = groups[use.samples,]
    Table.names = Table.names[order(Table.names$TKT),]
    Table.names = Table.names[order(Table.names$Genotype),]
    Table.names = Table.names[order(Table.names$Sex),]
    
    FDR.table = data.frame(GR.F      = FDR[c("GR.FxM", "GRF.CxDf", "GRF.CxOE", "GRxWT.FC", NA, NA)],
                           TktDfGR.F = FDR[c("GRF.CxDf", NA, NA, NA, "GRxWT.FDf", NA)],
                           TktOEGR.F = FDR[c("GRF.CxOE", NA, NA, NA, "GRxWT.FOE", NA)],
                           WT.F      = FDR[c("GRxWT.FC", NA, NA, "WT.FxM", "WTF.CxDf", "WTF.CxOE")],
                           TktDfWT.F = FDR[c(NA, "GRxWT.FDf", NA, "WTF.CxDf", NA, NA)],
                           TktOEWT.F = FDR[c(NA, NA, "GRxWT.FOE", "WTF.CxOE", NA, NA)],
                           GR.M      = FDR[c("GR.FxM", "GRM.CxDf", "GRM.CxOE", "GRxWT.MC", NA, NA)],
                           TktDfGR.M = FDR[c("GRM.CxDf", NA, NA, NA, "GRxWT.MDf", NA)],
                           TktOEGR.M = FDR[c("GRM.CxOE", NA, NA, NA, "GRxWT.MOE", NA)],
                           WT.M      = FDR[c("GRxWT.MC", NA, NA, "WT.FxM", "WTM.CxDf", "WTM.CxOE")],
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
      formatStyle(colnames(df),color='black', fontSize = '75%') %>%
      
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