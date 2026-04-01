if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(gplots, quietly = TRUE))
  install.packages("plotly")

if(!require(gplots, quietly = TRUE))
  BiocManager::install("goseq")

if(!require(org.Dm.eg.db, quietly = TRUE))
  BiocManager::install("org.Dm.eg.db")

if(!require(rrvgo, quietly = TRUE))
  BiocManager::install("rrvgo")

if(!require(pathview, quietly = TRUE))
  BiocManager::install("pathview")

if(!require(VennDiagram, quietly = TRUE))
  install.packages("VennDiagram")

library(plotly)
library(goseq)
library(org.Dm.eg.db)


git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/GeneralDataFiles/"

A4V.cpm = read.csv(paste0(git.dir, "A4V.cpmdata.csv"), row.names = 1)
A4V.meancpm = read.csv(paste0(git.dir, "A4V.meancpmdata.csv"), row.names = 1)

A4V.FDR = read.csv(paste0(git.dir, "A4V.FDRdata.csv"), row.names = 1)

A4V.FC = read.csv(paste0(git.dir, "A4V.FCdata.csv"), row.names = 1)

A4V.comparekey = setNames(colnames(A4V.FDR),
                          c('A3FHvS3FH', 'A3MHvS3MH', 'A9FHvS9FH', 'A9MHvS9MH', 'A40FHvS40FH',
                            'A3FHvA9FH', 'S3FHvS9FH', 'A3MHvA9MH', 'S3FHvS9MH',
                            'A9FHvA40FH', 'S9FHvS40FH', 'A3FHvA40FH', 'S3FHvS40FH',
                            
                            'A3FTvS3FT', 'A3MTvS3MT', 'A9FTvS9FT', 'A9MTvS9MT', 'A40FTvS40FT',
                            'A3FTvA9FT', 'S3FTvS9FT', 'A3MTvA9MT', 'S3MTvS9MT',
                            'A9FTvA40FT', 'S9FTvS40FT', 'A3FTvA40FT', 'S3FTvS40FT',
                            
                            'A3FAvS3FA', 'A3MAvS3MA', 'A9FAvS9FA', 'A9MAvS9MA', 'A40FAvS40FA',
                            'A3FAvA9FA', 'S3FAvS9FA', 'A3MAvA9MA', 'S3FAvS9MA',
                            'A9FAvA40FA', 'S9FAvS40FA', 'A3FAvA40FA', 'S3FAvS40FA',
                            
                            
                            'A3FHvA3MH', 'S3FHvS3MH', 'A9FHvA9MH', 'S9FHvS9MH',
                            'A3FTvA3MT', 'S3FTvS3MT', 'A9FTvA9MT', 'S9FTvS9MT',
                            'A3FAvA3MA', 'S3FAvS3MA', 'A9FAvA9MA', 'S9FAvS9MA'))

G85R.cpm = read.csv(paste0(git.dir, "TKT_cpmdata.csv"),row.names = 1)

G85R.cpm = G85R.cpm[grep('FBgn', row.names(G85R.cpm)),]

G85R.meancpm = read.csv(paste0(git.dir, "TKT_meancpmdata.csv"), row.names = 1)

colnames(G85R.meancpm) = c('GRFC', 'GRMC', 
                           'GRFDf', 'GRMDf', 'WTFDf', 'WTMDf',
                           'GRFOE', 'GRMOE', 'WTFOE', 'WTMOE',
                           'WTFC', 'WTMC')

G85R.FDR = read.csv(paste0(git.dir, "TKT.EdgeR.FDRTable.csv"),row.names = 1)

G85R.FC = read.csv(paste0(git.dir, "TKT.EdgeR.FCTable.csv"), row.names = 1)


G85R.comparekey = setNames(colnames(G85R.FDR),
                           c('GRFCvGRFDf', 'GRFCvGRFOE', 'WTFCvWTFDf', 'WTFCvWTFOE', 
                             'GRFCvWTFC', 'GRFDfvWTFDf', 'GRFOEvWTFOE',
                             'GRMCvGRMDf', 'GRMCvGRMOE', 'WTMCvWTMDf', 'WTMCvWTMOE', 
                             'GRMCvWTMC', 'GRMDfvWTMFDf', 'GRMOEvWTMOE',
                             'GRFCvGRMC', 'WTFCvWTMC'))


GeneIDKey = read.csv(paste0(git.dir, "GeneIDKey.csv"), row.names = 1)


genesingo = readRDS(gzcon(url(paste0(git.dir, "genesingo.RData"))))

GO.Names = read.csv(paste0(git.dir, "GO.names.csv"), row.names = 1)


plot.multipleGOterms = function(comparison,
                                keyword, 
                                sig.only = TRUE,
                                printallcategories = FALSE){
  
  all.specific.terms = GO.Names[grep(keyword, GO.Names$term),]

  
  if(comparison %in% colnames(G85R.FDR)){
    
    group1 = strsplit(names(G85R.comparekey)[1], split = "v")[[1]][1]
    group2 =  strsplit(names(G85R.comparekey)[1], split = "v")[[1]][2]
    
    GO.genes =  GeneIDKey[GeneIDKey$ensembl %in%  unique(unlist(genesingo[all.specific.terms$category])),
                          "FBgn"]
    GO.genes = intersect(unique(c(row.names(G85R.FDR), row.names(A4V.FDR))), GO.genes)
    GO.genes = GeneIDKey[GeneIDKey$FBgn %in% GO.genes, c("FBgn", "Symbol")]
    
    GO.genes$FDR = G85R.FDR[GO.genes$FBgn,comparison]
    GO.genes$FC = G85R.FC[GO.genes$FBgn,comparison]
    GO.genes$mean1 = G85R.meancpm[GO.genes$FBgn, group1]
    GO.genes$mean2 = G85R.meancpm[GO.genes$FBgn, group2]
    GO.genes = na.omit(GO.genes)
    
    ##insert cutoff filters here if wanted:
    if(sig.only == TRUE){
      GO.genes = GO.genes[GO.genes$FDR <= 0.05,]
    }
    
    volcano.data = GO.genes
    max.mean = max(na.omit(c(volcano.data$mean1, volcano.data$mean2)))
    volcano.data$size = apply(volcano.data[,c("mean1", "mean2")], MARGIN = 1, max)
    if(max.mean>200){
      volcano.data$size[volcano.data$size > 200] = 200
      max.mean = 200
    }  
    volcano.data$size = ((volcano.data$size/max.mean)+.5)*20
    volcano.data$FDR = -log2(volcano.data$FDR)
    volcano.data.color = volcano.data$FDR
    volcano.data.color[abs(volcano.data.color) < -log2(.05)] = 0
    volcano.data.color = 100*(volcano.data.color/max(na.omit(volcano.data.color)))
    volcano.data.color = 101+(volcano.data$FC/abs(volcano.data$FC))*volcano.data.color
    
    volcano.color <- colorRampPalette(c("royalblue", "deepskyblue", "lightblue",
                                        "white",
                                        "lightpink", "deeppink", "firebrick"))(201)
    volcano.color[101] = "lightgrey"
    volcano.data$Color = volcano.color[volcano.data.color]
    
    
    fig = plot_ly(data = volcano.data,
                  x = ~FC,
                  y = ~FDR,
                  type = 'scatter',
                  mode = 'markers',
                  marker = list(color = ~Color, colors = ~Color, size = volcano.data$size,
                                line = list(color = 'black', width = .5)),
                  hoverinfo = "text",
                  hovertext = paste("Gene:", volcano.data$FBgn,
                                    "\nSymbol:", volcano.data$Symbol,
                                    "\n-log2(FDR): ", round(volcano.data$FDR,2),
                                    "\nFC: ", round(volcano.data$FC,2),
                                    "\n", group1, "mean cpm: ",
                                    round(volcano.data$mean1, 1),
                                    "\n", group2, "mean cpm: ",
                                    round(volcano.data$mean2, 1)))
    fig = fig %>% layout(title = paste0(group1, " vs. ", group2))

  }
  
  
  if(comparison %in% colnames(A4V.FDR)){
    
    
    group1 = strsplit(names(A4V.comparekey)[1], split = "v")[[1]][1]
    group2 =  strsplit(names(A4V.comparekey)[1], split = "v")[[1]][2]
    
    GO.genes =  GeneIDKey[GeneIDKey$ensembl %in% 
                            unique(unlist(genesingo[all.specific.terms$category])),
                          "FBgn"]
    GO.genes = intersect(unique(c(row.names(G85R.FDR), row.names(A4V.FDR))), GO.genes)
    GO.genes = GeneIDKey[GeneIDKey$FBgn %in% GO.genes, c("FBgn", "Symbol")]
    
    GO.genes$FDR = A4V.FDR[GO.genes$FBgn,comparison]
    GO.genes$FC = A4V.FC[GO.genes$FBgn,comparison]
    GO.genes$mean1 = A4V.meancpm[GO.genes$FBgn, group1]
    GO.genes$mean2 = A4V.meancpm[GO.genes$FBgn, group2]
    GO.genes = na.omit(GO.genes)
    
    ##insert cutoff filters here if wanted:
    if(sig.only == TRUE){
      GO.genes = GO.genes[GO.genes$FDR <= 0.05,]
    }
    
    volcano.data = GO.genes
    max.mean = max(na.omit(c(volcano.data$mean1, volcano.data$mean2)))
    volcano.data$size = apply(volcano.data[,c("mean1", "mean2")], MARGIN = 1, max)
    if(max.mean>200){
      volcano.data$size[volcano.data$size > 200] = 200
      max.mean = 200
    }  
    volcano.data$size = ((volcano.data$size/max.mean)+.5)*20
    volcano.data$FDR = -log2(volcano.data$FDR)
    volcano.data.color = volcano.data$FDR
    volcano.data.color[abs(volcano.data.color) < -log2(.05)] = 0
    volcano.data.color = 100*(volcano.data.color/max(na.omit(volcano.data.color)))
    volcano.data.color = 101+(volcano.data$FC/abs(volcano.data$FC))*volcano.data.color
    
    volcano.color <- colorRampPalette(c("royalblue", "deepskyblue", "lightblue",
                                        "white",
                                        "lightpink", "deeppink", "firebrick"))(201)
    volcano.color[101] = "lightgrey"
    volcano.data$Color = volcano.color[volcano.data.color]
    
    
    fig = plot_ly(data = volcano.data,
                  x = ~FC,
                  y = ~FDR,
                  type = 'scatter',
                  mode = 'markers',
                  marker = list(color = ~Color, colors = ~Color, size = volcano.data$size,
                                line = list(color = 'black', width = .5)),
                  hoverinfo = "text",
                  hovertext = paste("Gene:", volcano.data$FBgn,
                                    "\nSymbol:", volcano.data$Symbol,
                                    "\n-log2(FDR): ", round(volcano.data$FDR,2),
                                    "\nFC: ", round(volcano.data$FC,2),
                                    "\n", group1, "mean cpm: ",
                                    round(volcano.data$mean1, 1),
                                    "\n", group2, "mean cpm: ",
                                    round(volcano.data$mean2, 1)))
    fig = fig %>% layout(title = paste0(group1, " vs. ", group2))

  }
  
  if(printallcategories = TRUE){
    all.specific.terms
  }
  
  fig
}

