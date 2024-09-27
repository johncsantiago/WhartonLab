
library(plotly)
library(goseq)
library(org.Dm.eg.db)
library(rrvgo)
library(pathview)
library(VennDiagram)

A4V.cpm = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/A4V_BodySections_RNAseq/A4V.cpmdata.csv", row.names = 1)
A4V.meancpm = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/A4V_BodySections_RNAseq/A4V.meancpmdata.csv", row.names = 1)
A4V.FDR = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/A4V_BodySections_RNAseq/A4V.FDRdata.csv", row.names = 1)
A4V.FC = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/A4V_BodySections_RNAseq/A4V.FCdata.csv", row.names = 1)
A4V.meta = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/A4V_BodySections_RNAseq/A4V.metadata.csv", row.names = 1)

A4V.comparekey = setNames(colnames(A4V.FDR),
                          c('A3FHvS3FH', 'A3MHvS3MH', 'A9FHvS9FH', 'A9MHvS9MH', 'A40FHvS40FH',
                            'A3FHvA9FH', 'S3FHvS9FH', 'A3MHvA9MH', 'S3FHvS9MH',
                            'A9FHvA40FH', 'S9FHvS40FH', 'A3FHvA40FH', 'S3FHvS40FH',
                            
                            'A3FTvS3FT', 'A3MTvS3MT', 'A9FTvS9FT', 'A9MTvS9MT', 'A40FTvS40FT',
                            'A3FTvA9FT', 'S3FTvS9FT', 'A3MTvA9MT', 'S3FTvS9MT',
                            'A9FTvA40FT', 'S9FTvS40FT', 'A3FTvA40FT', 'S3FTvS40FT',
                            
                            'A3FAvS3FA', 'A3MAvS3MA', 'A9FAvS9FA', 'A9MAvS9MA', 'A40FAvS40FA',
                            'A3FAvA9FA', 'S3FAvS9FA', 'A3MAvA9MA', 'S3FAvS9MA',
                            'A9FAvA40FA', 'S9FAvS40FA', 'A3FAvA40FA', 'S3FAvS40FA',
                            
                            
                            'A3FHvA3MH', 'S3FHvS3MH', 'A9FHvA9MH', 'S9FHvS9MH',
                            'A3FTvA3MT', 'S3FTvS3MT', 'A9FTvA9MT', 'S9FTvS9MT',
                            'A3FAvA3MA', 'S3FAvS3MA', 'A9FAvA9MA', 'S9FAvS9MA'))

G85R.cpm = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/TKT_RNAseq/CountTables/TKT_cpmdata.csv",row.names = 1)
G85R.cpm = G85R.cpm[grep('FBgn', row.names(G85R.cpm)),]
G85R.meancpm = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/TKT_RNAseq/CountTables/TKT_meancpmdata.csv", row.names = 1)
colnames(G85R.meancpm) = c('GRFC', 'GRMC', 
                           'GRFDf', 'GRMDf', 'WTFDf', 'WTMDf',
                           'GRFOE', 'GRMOE', 'WTFOE', 'WTMOE',
                           'WTFC', 'WTMC')
G85R.FDR = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/TKT_RNAseq/CountTables/TKT.EdgeR.FDRTable.csv",row.names = 1)
G85R.FC = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/TKT_RNAseq/CountTables/TKT.EdgeR.FCTable.csv", row.names = 1)
G85R.meta  = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/TKT_RNAseq/CountTables/TKT.metadata.csv",row.names = 1)

G85R.comparekey = setNames(colnames(G85R.FDR),
                           c('GRFCvGRFDf', 'GRFCvGRFOE', 'WTFCvWTFDf', 'WTFCvWTFOE', 
                             'GRFCvWTFC', 'GRFDfvWTFDf', 'GRFOEvWTFOE',
                             'GRMCvGRMDf', 'GRMCvGRMOE', 'WTMCvWTMDf', 'WTMCvWTMOE', 
                             'GRMCvWTMC', 'GRMDfvWTMFDf', 'GRMOEvWTMOE',
                             'GRFCvGRMC', 'WTFCvWTMC'))

G85R.metab.comparekey = setNames(c("GR.CxDf", "GR.CxOE", "WT.CxDf", "WT.CxOE",
                                   "GRxWT.C", "GRxWT.Df", "GRxWT.OE",
                                   "GR.CxDf", "GR.CxOE", "WT.CxDf", "WT.CxOE",
                                   "GRxWT.C", "GRxWT.Df", "GRxWT.OE"), names(G85R.comparekey)[1:14])

G85R.metabnorm = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/NormalizedMetabolomicData.csv", row.names = 1)
G85R.metabmean = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/MeanNormalizedMetabolomicData.csv", row.names = 1)
metab.id = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/MetaboliteMultiKEGGIDKey.csv", row.names = 1)
metab.id = setNames(metab.id$metab, metab.id$KEGG)
G85R.metabFC = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/MetaboliteFCs.csv", row.names = 1)
G85R.metabFDR = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/MetaboliteFDRs.csv", row.names = 1)
G85R.metabp = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/Metabolitepvals.csv", row.names = 1)

GeneIDKey = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/GeneIDKey.csv", row.names = 1)
KEGG.MetabKey = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Metabolomics/RawMetabolomicData.csv", row.names = 1)
KEGG.MetabKey = setNames(KEGG.MetabKey[-c(1:2), 2], row.names(KEGG.MetabKey)[-c(1:2)])
genesingo = readRDS("/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/genesingo.RData")
GenesInKegg = readRDS("/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/kegg.symbol2path.RData")
KEGG.Names = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/KEGG.names.csv", row.names = 1)
GO.Names = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/GeneralDataFiles/GO.names.csv", row.names = 1)





Enrichment = function(FDR){
  
  sigs = row.names(FDR)[FDR[,1] < .05]
  bg = row.names(FDR)
  
  
  sigsKEY = GeneIDKey[sigs,]
  sigs = na.omit(sigsKEY$Symbol)
  bgKEY = GeneIDKey[bg,]
  bg = na.omit(bgKEY$Symbol)
  
  genes = setNames(rep(0, length(bg)), bg)
  genes[sigs] = 1
  pwf=nullp(genes,"dm3","geneSymbol")
  GO.wall=goseq(pwf,"dm3","geneSymbol")
  GO.wall$adjp=p.adjust(GO.wall$over_represented_pvalue,method="BH")
  row.names(GO.wall) = GO.wall$category
  KEGG=goseq(pwf,gene2cat=GenesInKegg)
  KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
  row.names(KEGG)=substr(as.character(KEGG[,1]),6,13)
  KEGG$Name=KEGG.Names[row.names(KEGG),1]
  
  enrichment.data = list(GO.wall, KEGG)
  names(enrichment.data) = c("GO", "KEGG")
  return(enrichment.data)
}

Subset.Enrichment = function(bg, geneset){
  
  sigs = geneset
  
  
  sigsKEY = GeneIDKey[sigs,]
  sigs = na.omit(sigsKEY$Symbol)
  bgKEY = GeneIDKey[bg,]
  bg = na.omit(bgKEY$Symbol)
  
  genes = setNames(rep(0, length(bg)), bg)
  genes[sigs] = 1
  pwf=nullp(genes,"dm3","geneSymbol")
  GO.wall=goseq(pwf,"dm3","geneSymbol")
  GO.wall$adjp=p.adjust(GO.wall$over_represented_pvalue,method="BH")
  row.names(GO.wall) = GO.wall$category
  KEGG=goseq(pwf,gene2cat=GenesInKegg)
  KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
  row.names(KEGG)=substr(as.character(KEGG[,1]),6,13)
  KEGG$Name=KEGG.Names[row.names(KEGG),1]
  
  enrichment.data = list(GO.wall, KEGG)
  names(enrichment.data) = c("GO", "KEGG")
  return(enrichment.data)
}

tree = function(cat.data, GO.ontology){
  
  GO.data = cat.data$GO
  
  go_analysis = GO.data$category[GO.data$over_represented_pvalue <= .005 & 
                                   GO.data$ontology == GO.ontology]
  scores = setNames(-log10(GO.data$over_represented_pvalue), go_analysis)
  
  simMatrix <- calculateSimMatrix(go_analysis,
                                  orgdb="org.Dm.eg.db",
                                  ont=GO.ontology,
                                  method="Rel")
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Dm.eg.db")
  
  return(reducedTerms)
}


plotparentcat.genes = function(parent.category, FC, FDR, mean.data){
  
  FC.data = setNames(FC[,1], row.names(FC))
  FDR.data = setNames(FDR[,1], row.names(FDR))
  
  if(length(rt[rt$parentTerm == parent.category, "go"]) == 0){
    parent.category = rt$parentTerm[1]
  }
  
  ##select F if you want point size to scale with mean cpm of variable data
  set.size = F
  ##select T if you only want the genes in category to show
  catgenes.only = T
  
  GOterm.genes = vector('list', length = length(unique(rt$parent)))
  names(GOterm.genes) = unique(rt$parent)
  
  i=1
  while(i<=length(GOterm.genes)){
    goinparent = rt[grep(names(GOterm.genes)[i], rt$parent),'go']
    genesinparent = unlist(genesingo[goinparent])
    GOterm.genes[[i]] = GeneIDKey[GeneIDKey$ensembl %in% genesinparent, 'FBgn']
    i = i +1
  }
  
  volcano.data = data.frame(Symbol = GeneIDKey[names(FDR.data), "Symbol"], 
                            FDR = -log2(FDR.data),
                            FC = FC.data,
                            Color = 'grey',
                            size = 5*log(mean.data[,1]+2),
                            Variable.cpm = mean.data[,1],
                            Control.cpm = mean.data[,2],
                            GO.terms = "none")
  
  if(set.size == T){
    volcano.data$size = 15
    volcano.data$size[volcano.data$FDR <= -log2(.05)] = 5
  }
  
  parent.genes = GOterm.genes[[intersect(rt[rt$parentTerm == parent.category, "go"],
                                         names(GOterm.genes))]]
  parent.genes = intersect(parent.genes, row.names(volcano.data))
  volcano.data$Color = 'honeydew'
  volcano.data[parent.genes, 'Color'] = 'deeppink'
  volcano.data = volcano.data[order(volcano.data$Color, decreasing = T),]
  main.title = parent.category
  
  if(catgenes.only == T){
    volcano.data = volcano.data[volcano.data$Color == 'deeppink',]
    volcano.data$Color[volcano.data$FDR < -log2(.05)] = 'honeydew'
  }
  
  fig = plot_ly(data = volcano.data,
                x = ~FC,
                y = ~FDR,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color, colors = ~Color, size = volcano.data$size,
                              line = list(color = 'black', width = .5)),
                hoverinfo = "text",
                hovertext = paste("Gene:", volcano.data$Symbol,
                                  "\n-log2(FDR): ", round(volcano.data$FDR,2),
                                  "\nFC: ", round(volcano.data$FC,2),
                                  "\n ", colnames(mean.data)[1], "mean cpm: ",
                                  round(volcano.data$Variable.cpm, 1),
                                  "\n ", colnames(mean.data)[2], "mean cpm: ",
                                  round(volcano.data$Control.cpm, 1)))
  fig = fig %>% layout(title = main.title)
  
  return(fig)
}


GO20 = function(cat.data, GO.ontology, FC, FDR, usecats){

  GO.data = cat.data$GO
  
  FC.data = setNames(FC[,1], row.names(FC))
  FDR.data = setNames(FDR[,1], row.names(FDR))
  
  GO.data = GO.data[GO.data$ontology == GO.ontology,
                    c("category", "term","adjp", "numDEInCat", "numInCat")]
  
  data = GO.data[usecats[length(usecats):1],c("term","adjp","numDEInCat", "numInCat", "category")]
  if(usecats[1] == 'sigs.only'){
    data = GO.data[GO.data$adjp < .05,c("term","adjp","numDEInCat", "numInCat", "category")]
    if(nrow(data)>20){
      data = data[1:20,]
    }
    data = data[nrow(data):1,]
  }
  
  data$Score = -log10(data$adjp)
  
  temp = nchar(data$term)
  data$label = data$term
  ##data$label[temp > 40] = paste0(substr(data$term[temp>40], 1 ,35), "...")
  
  data$term = factor(data$term, levels = c(unique(data$term)))
  data$up = 0
  data$down = 0
  data$allup = 0
  data$alldown = 0
  
  sigs = names(FDR.data[FDR.data <= 0.05])
  
  i = 1
  while(i <= nrow(data)){
    go.id = row.names(data)[i]
    go.genes = genesingo[[go.id]]
    go.genes = GeneIDKey[GeneIDKey$ensembl %in% go.genes, "FBgn"]
    go.genes = intersect(go.genes, names(FC.data))
    sig.go.genes = intersect(go.genes, sigs)
    go.FC = FC.data[go.genes]
    sig.go.FC = FC.data[sig.go.genes]  
    
    data$allup[i] = length(go.FC[go.FC>0])
    data$alldown[i] = length(go.FC[go.FC<0])
    data$up[i] = length(sig.go.FC[sig.go.FC>0])
    data$down[i] = length(sig.go.FC[sig.go.FC<0])
    i = i +1
  }
  
  data$down = data$down * -1
  data$alldown = data$alldown * -1
  
  m <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 100,
    pad = 4
  )
  
  data$label = factor(data$label, levels = data$label)
  
  fig <- plot_ly(data, 
                 x = ~up, y = ~label, 
                 type = 'bar',
                 name = "Sig. Increased FC",
                 marker = list(color = 'brown',
                               line = list(color = "black", 
                                           width = 1.5)),
                 hoverinfo = "text",
                 hovertext = paste("Sample:", data$term,
                                   "\nGO ID: ", data$category,
                                   "\nTotal DE: ", data$numDEInCat,
                                   "\nTotal in Cat.: ", data$numInCat,
                                   "\nDE up: ", data$up,
                                   "\nDE down: ", data$down,
                                   "\nFDR: ", signif(data$adjp, 3))) %>%
    add_trace(x = ~down, 
              y ~term, 
              marker = list(color = 'steelblue', 
                            line = list(color = "black", width = 1.5)), 
              name = "Sig. Decreased FC")%>%
    
    layout(xaxis = list(title = "Number of Genes"),
           yaxis = list(title = ""),
           title = paste0("Top 20 Most Enriched ",GO.ontology, " GO Terms"),
           margin = m, barmode = 'overlay')
  
  return(fig)
}

reversal.GO20 = function(cat.data, GO.ontology, FC, FDR, usecats, sigs){
  
  GO.data = cat.data$GO
  
  FC.data = FC
  FDR.data = FDR
  
  GO.data = GO.data[GO.data$ontology == GO.ontology,
                    c("category", "term","adjp", "numDEInCat", "numInCat")]
  
  data = GO.data[usecats[length(usecats):1],c("term","adjp","numDEInCat", "numInCat", "category")]
  if(usecats[1] == 'sigs.only'){
    data = GO.data[GO.data$adjp < .05,c("term","adjp","numDEInCat", "numInCat", "category")]
    if(nrow(data)>20){
      data = data[1:20,]
    }
    data = data[nrow(data):1,]
  }
  
  data$Score = -log10(data$adjp)
  
  temp = nchar(data$term)
  data$label = data$term
  ##data$label[temp > 40] = paste0(substr(data$term[temp>40], 1 ,35), "...")
  
  data$term = factor(data$term, levels = c(unique(data$term)))
  data$up = 0
  data$down = 0
  data$allup = 0
  data$alldown = 0
  
  sigs = sigs
  
  i = 1
  while(i <= nrow(data)){
    go.id = row.names(data)[i]
    go.genes = genesingo[[go.id]]
    go.genes = GeneIDKey[GeneIDKey$ensembl %in% go.genes, "FBgn"]
    go.genes = intersect(go.genes, names(FC.data))
    sig.go.genes = intersect(go.genes, sigs)
    go.FC = FC.data[go.genes]
    sig.go.FC = FC.data[sig.go.genes]  
    
    data$allup[i] = length(go.FC[go.FC>0])
    data$alldown[i] = length(go.FC[go.FC<0])
    data$up[i] = length(sig.go.FC[sig.go.FC>0])
    data$down[i] = length(sig.go.FC[sig.go.FC<0])
    i = i +1
  }
  
  data$down = data$down * -1
  data$alldown = data$alldown * -1
  
  m <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 100,
    pad = 4
  )
  
  data$label = factor(data$label, levels = data$label)
  
  fig <- plot_ly(data, 
                 x = ~up, y = ~label, 
                 type = 'bar',
                 name = "Sig. Increased FC",
                 marker = list(color = 'brown',
                               line = list(color = "black", 
                                           width = 1.5)),
                 hoverinfo = "text",
                 hovertext = paste("Sample:", data$term,
                                   "\nGO ID: ", data$category,
                                   "\nTotal DE: ", data$numDEInCat,
                                   "\nTotal in Cat.: ", data$numInCat,
                                   "\nDE up: ", data$up,
                                   "\nDE down: ", data$down,
                                   "\nFDR: ", signif(data$adjp, 3))) %>%
    add_trace(x = ~down, 
              y ~term, 
              marker = list(color = 'steelblue', 
                            line = list(color = "black", width = 1.5)), 
              name = "Sig. Decreased FC")%>%
    
    layout(xaxis = list(title = "Number of Genes"),
           yaxis = list(title = ""),
           title = paste0("Top 20 Most Enriched ",GO.ontology, " GO Terms"),
           margin = m, barmode = 'overlay')
  
  return(fig)
}

KEGG20 = function(cat.data, FC, FDR, usecats){
  
  KEGG.data = cat.data$KEGG
  FC.data = setNames(FC[,1], row.names(FC))
  FDR.data = setNames(FDR[,1], row.names(FDR))
  
  data = KEGG.data[, c("category", "Name","adjp", "numDEInCat", "numInCat")]
  
  data$term = data$Name
  if(length(unique(data$term)) < length(data$term)){
    data$term[duplicated(data$term)] = paste(data$term[duplicated(data$term)], "2")
  }  
  data = na.omit(data[data$Name != 'Metabolic pathways',])
  
  if(usecats[1] != 'sigs.only'){
    data = data[usecats[length(usecats):1],]
  }
  
  if(usecats[1] == 'sigs.only'){
    data = data[data$adjp < .05,]
    data = data[ncol(data):1,]
    data = na.omit(data)
  }
  
  data$Score = -log10(data$adjp)
  data$term = factor(data$term, levels = c(unique(data$term)))
  data$up = 0
  data$down = 0
  data$allup = 0
  data$alldown = 0
  
  sigs = names(FDR.data[FDR.data <= 0.05])
  
  i = 1
  while(i <= nrow(data)){
    kegg.id = data$category[i]
    kegg.genes = names(GenesInKegg[grep(kegg.id, GenesInKegg)])
    kegg.genes = GeneIDKey[GeneIDKey$Symbol %in% kegg.genes, "FBgn"]
    kegg.genes = intersect(kegg.genes, names(FC.data))
    sig.kegg.genes = intersect(kegg.genes, sigs)
    kegg.FC = FC.data[kegg.genes]
    sig.kegg.FC = FC.data[sig.kegg.genes]
    data$allup[i] = length(kegg.FC[kegg.FC>0])
    data$alldown[i] = length(kegg.FC[kegg.FC<0])
    data$up[i] = length(sig.kegg.FC[sig.kegg.FC>0])
    data$down[i] = length(sig.kegg.FC[sig.kegg.FC<0])
    i = i +1
  }
  
  data$down = data$down * -1
  data$alldown = data$alldown * -1
  
  
  m <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 100,
    pad = 4
  )
  
  fig <- plot_ly(data,
                 x = ~up, y = ~term, 
                 type = 'bar',
                 name = "Sig. Increased FC",
                 marker = list(color = 'brown',
                               line = list(color = "black", 
                                           width = 1.5)),
                 hoverinfo = "text",
                 hovertext = paste("Sample:", data$term,
                                   "\nKEGG ID: ", data$category,
                                   "\nTotal DE: ", data$numDEInCat,
                                   "\nTotal in Cat.: ", data$numInCat,
                                   "\nDE up: ", data$up,
                                   "\nDE down: ", data$down,
                                   "\nFDR: ", signif(data$adjp, 3))) %>%
    add_trace(x = ~down, 
              y ~term, 
              marker = list(color = 'steelblue', 
                            line = list(color = "black", 
                                        width = 1.5)), 
              name = "Sig. Decreased FC")%>%
    layout(xaxis = list(title = "Number of Genes"),
           yaxis = list(title = ""),
           title = paste0("Top 20 Most Enriched KEGG Terms"),
           margin = m, barmode = 'overlay')
  
  return(fig)
  
}


KEGG.diagram = function(PathwayID, FC, use, metab.data){

  ##File name prefix
  out.id = paste(comparison, ' and ', secondary.comparison)
  
  ##table with the FC data you want to use for coloring with FBgn for row names 
  ##if only analyzing 1 condition, use setNames to make a named vector
  ##A4V.FC or G85R.FC
  
  FC.data =  FC
  
  ##all or a subset of row.names used in FC.data
  usegenes = row.names(FDR)
  
  #Any column name from G85R.metabFC or 'none'
  metab.data = metab.data
  
  
  if(use == 'primary'){
    
    ##File name prefix
    out.id = comparison
    
    ##table with the FC data you want to use for coloring with FBgn for row names 
    ##if only analyzing 1 condition, use setNames to make a named vector
    ##A4V.FC or G85R.FC
    
    FC.data =  FC[,1]
    
    #Any column name from G85R.metabFC or 'none'
    metab.data = metab.data[,1]
  }
  
  ##KEGG Mapping Function
  if(length(intersect(list.files(), "KEGG_Image_Files")) == 0){
    dir.create(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/KEGG_Image_Files/"))
  }
  Home = dirname(rstudioapi::getSourceEditorContext()$path)
  OutputDir = paste0(Home, "/KEGG_Image_Files/")
  setwd(OutputDir)
  sigs = usegenes
  ##sigs = names(FDR.data[FDR.data <= .05])
  
  if(class(FC.data)=="numeric"){
    deg.data=setNames(FC.data[sigs],GeneIDKey[sigs, "ensembl"])
  }
  
  if(class(FC.data)=="data.frame"){
    deg.data = as.matrix(FC.data)
    row.names(deg.data) = GeneIDKey[row.names(FC.data), "ensembl"]
  }
  
  
  if(length(G85R.metabFC[metab.id, metab.data]) > 1 ){
    metab.data = setNames(log2(G85R.metabFC[metab.id, metab.data]), names(metab.id))
  }
  
  
  
  pv.out <- pathview(gene.data =  deg.data,
                     cpd.data = metab.data,
                     pathway.id = PathwayID,
                     species = "dme",
                     kegg.native = T,
                     limit=list(gene=c(-2,2)),
                     node.sum="max.abs",
                     low="royalblue",mid = "grey",high="firebrick",
                     out.suffix = out.id,
                     bins=list(genes=30),
                     plot.col.key = F,
                     match.data = F)
  setwd(Home)
  
  
  if(class(FC.data)=="numeric"){
    img.file = paste0(OutputDir, "dme",PathwayID, ".", out.id, ".png")
  }
  
  if(class(FC.data)=="data.frame"){
    img.file = paste0(OutputDir, "dme",PathwayID, ".", out.id, ".multi.png")
  }
  knitr::include_graphics(img.file)
}

plotkegg.genes = function(specific.kegg, FC, FDR, mean.data){
  
  FC.data = setNames(FC[,1], row.names(FC))
  FDR.data = setNames(FDR[,1], row.names(FDR))
  ##table with mean cpm of variable and control conditions
  mean.data = mean.data[,1:2]
  
  ##select F if you want point size to scale with mean cpm of variable data
  set.size = F
  ##select T if you only want the genes in category to show
  catgenes.only = T
  
  volcano.data = data.frame(Symbol = GeneIDKey[names(FDR.data), "Symbol"],  
                            FBgn = names(FDR.data),
                            FDR = -log2(FDR.data),
                            FC = FC.data,
                            Color = 'grey',
                            size = 5*log(mean.data[,1]),
                            Variable.cpm = mean.data[,1],
                            Control.cpm = mean.data[,2])
  
  if(set.size == T){
    volcano.data$size = 15
    volcano.data$size[volcano.data$FDR <= -log2(.05)] = 5
  }
  
  
  kegg.id = row.names(KEGG.Names)[c(grep(specific.kegg, KEGG.Names), 
                                    grep(specific.kegg, row.names(KEGG.Names)))]
  kegg.genes = names(GenesInKegg[grep(kegg.id, GenesInKegg)])
  volcano.data$Color = 'honeydew'
  #row.names(volcano.data) = volcano.data$Symbol
  volcano.data[volcano.data$Symbol %in% kegg.genes, 'Color'] = 'deeppink'
  volcano.data = volcano.data[order(volcano.data$Color, decreasing = T),]
  main.title = KEGG.Names[c(grep(specific.kegg, KEGG.Names),
                            grep(specific.kegg, row.names(KEGG.Names))),]
  
  if(catgenes.only == T){
    volcano.data = volcano.data[volcano.data$Color == 'deeppink',]
    volcano.data$Color[volcano.data$FDR < -log2(.05)] = 'honeydew'
  }
  

  
  
  fig = plot_ly(data = volcano.data,
                x = ~FC,
                y = ~FDR,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color, colors = ~Color, size = volcano.data$size,
                              line = list(color = 'black', width = .5)),
                hoverinfo = "text",
                hovertext = paste("Gene:", volcano.data$Symbol,
                                  "\nFBgn:", volcano.data$FBgn,
                                  "\n-log2(FDR): ", round(volcano.data$FDR,2),
                                  "\nFC: ", round(volcano.data$FC,2),
                                  "\n ", colnames(mean.data)[1], "mean cpm: ",
                                  round(volcano.data$Variable.cpm, 1),
                                  "\n ", colnames(mean.data)[2], "mean cpm: ",
                                  round(volcano.data$Control.cpm, 1)))
  fig = fig %>% layout(title = main.title)
  
  
  return(fig)
}

volcano = function(FC, FDR, mean.data){
  
  FC.data = setNames(FC[,1], row.names(FC))
  FDR.data = setNames(FDR[,1], row.names(FDR))
  mean.data = mean.data[,1:2]
  
  ##select F if you want point size to scale with mean cpm of variable data
  set.size = F
  ##select T if you only want the genes in category to show
  catgenes.only = T
  
  volcano.data = data.frame(Symbol = GeneIDKey[names(FDR.data), "Symbol"], 
                            FDR = -log2(FDR.data),
                            FC = FC.data,
                            Color = 'grey',
                            size = 5*log(mean.data[,1]+2),
                            Variable.cpm = mean.data[,1],
                            Control.cpm = mean.data[,2])
  
  if(set.size == T){
    volcano.data$size = 15
    volcano.data$size[volcano.data$FDR <= -log2(.05)] = 5
  }
  
  volcano.data$Color[volcano.data$FDR >= -log2(.05) & volcano.data$FC < 0] = 'lightblue'
  volcano.data$Color[volcano.data$FDR >= -log2(.0001) & volcano.data$FC < 0] = 'steelblue'
  volcano.data$Color[volcano.data$FDR >= -log2(.000001) & volcano.data$FC < 0] = 'dodgerblue'
  
  volcano.data$Color[volcano.data$FDR >= -log2(.05) & volcano.data$FC > 0] = 'lightcoral'
  volcano.data$Color[volcano.data$FDR >= -log2(.0001) & volcano.data$FC > 0] = 'tomato'
  volcano.data$Color[volcano.data$FDR >= -log2(.000001) & volcano.data$FC > 0] = 'firebrick'
  
  volcano.data$size[volcano.data$FDR <= -log2(.05)] = 3

  main.title = paste(colnames(mean.data)[1], "vs. ", colnames(mean.data)[2])

if(set.size == T){
  volcano.data$size = 15
  volcano.data$size[volcano.data$FDR <= -log2(.05)] = 5
}


volcano.data = volcano.data[volcano.data$FDR >= -log2(.05),]

fig = plot_ly(data = volcano.data,
              x = ~FC,
              y = ~FDR,
              type = 'scatter',
              mode = 'markers',
              marker = list(color = ~Color, colors = ~Color, size = volcano.data$size,
                            line = list(color = 'black', width = .5)),
              hoverinfo = "text",
              hovertext = paste("Gene:", volcano.data$Symbol,
                                "\n-log2(FDR): ", round(volcano.data$FDR,2),
                                "\nFC: ", round(volcano.data$FC,2),
                                "\n ", colnames(mean.data)[1], "mean cpm: ",
                                round(volcano.data$Variable.cpm, 1),
                                "\n ", colnames(mean.data)[2], "mean cpm: ",
                                round(volcano.data$Control.cpm, 1)))
fig = fig %>% layout(title = main.title)

##color = ~Color,
##colors = 'Spectral')
##color = ~Color,
##colors = ~Color)
##marker = list(colorscale = list(c(0,.5,1), c("blue","yellow", "red")), color = ~Color))

fig

}

plotgo.genes = function(specific.go, FC.data, FDR.data, mean.data){
  
  FC.data = setNames(FC[,1], row.names(FC))
  FDR.data = setNames(FDR[,1], row.names(FDR))
  mean.data = mean.data[,c(1,2)]
  
  ##select F if you want point size to scale with mean cpm of variable data
  set.size = F
  ##select T if you only want the genes in category to show
  catgenes.only = T
  
  volcano.data = data.frame(Symbol = GeneIDKey[names(FDR.data), "Symbol"], 
                            FBgn = names(FDR.data),
                            FDR = -log2(FDR.data),
                            FC = FC.data,
                            Color = 'grey',
                            size = 5*log(mean.data[,1]),
                            Variable.cpm = mean.data[,1],
                            Control.cpm = mean.data[,2])
  
  if(set.size == T){
    volcano.data$size = 15
    volcano.data$size[volcano.data$FDR <= -log2(.05)] = 5
  }
  
  
  go.id = row.names(GO.Names)[c(grep(specific.go, GO.Names$category),
                                grep(specific.go, GO.Names$term))]
  go.genes = (genesingo[go.id])[[1]]
  go.genes = GeneIDKey[GeneIDKey$ensembl %in% go.genes, 'Symbol']
  volcano.data$Color = 'honeydew'
  #row.names(volcano.data) = volcano.data$Symbol
  volcano.data[volcano.data$Symbol %in% go.genes, 'Color'] = 'deeppink'
  volcano.data = volcano.data[order(volcano.data$Color, decreasing = T),]
  main.title = GO.Names[go.id, 'term']
  
  if(catgenes.only == T){
    volcano.data = volcano.data[volcano.data$Color == 'deeppink',]
    volcano.data$Color[volcano.data$FDR < -log2(.05)] = 'honeydew'
  }
  
  
  fig = plot_ly(data = volcano.data,
                x = ~FC,
                y = ~FDR,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color, colors = ~Color, size = volcano.data$size,
                              line = list(color = 'black', width = .5)),
                hoverinfo = "text",
                hovertext = paste("Gene:", volcano.data$Symbol,
                                  "\nFBgn:", volcano.data$FBgn,
                                  "\n-log2(FDR): ", round(volcano.data$FDR,2),
                                  "\nFC: ", round(volcano.data$FC,2),
                                  "\n ", colnames(mean.data)[1], "mean cpm: ",
                                  round(volcano.data$Variable.cpm, 1),
                                  "\n ", colnames(mean.data)[2], "mean cpm: ",
                                  round(volcano.data$Control.cpm, 1)))
  fig = fig %>% layout(title = main.title)
  
  
  return(fig)
}

comparego.genes = function(specific.go, FC, FDR, mean.data){
  
  FC.data = FC
  FDR.data = FDR
  
  go.id = row.names(GO.Names)[c(grep(specific.go, GO.Names$category),
                                grep(specific.go, GO.Names$term))]
  go.genes = (genesingo[go.id])[[1]]
  go.genes = GeneIDKey[GeneIDKey$ensembl %in% go.genes, 'FBgn']
  
  
data = data.frame(Symbol = GeneIDKey[go.genes, "Symbol"],
                  FBgn = go.genes,
                  FDR1 = -log2(FDR.data[go.genes, 1]),
                  FDR2 = -log2(FDR.data[go.genes, 2]),
                  FC1 = FC.data[go.genes, 1],
                  FC2 = FC.data[go.genes, 2],
                  Color = 'grey',
                  size = 1,
                  Variable1.cpm = mean.data[go.genes,1],
                  Control1.cpm = mean.data[go.genes,2],
                  Variable2.cpm = mean.data[go.genes,3],
                  Control2.cpm = mean.data[go.genes,4])

  main.title = GO.Names[go.id, 'term']
  
  data$Color[(data$FC1 - 1) > data$FC2] = 'royalblue'
  data$Color[(data$FC2 - 1) > data$FC1] = 'firebrick' 
  data$size = 2.5*apply(data[,c('FDR1', 'FDR2')], MARGIN = 1, max)
  data$size[data$size < -2.5*log2(.05)] = 5
  data = na.omit(data)
  
  fig = plot_ly(data = data,
                x = ~FC1,
                y = ~FC2,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color,
                              colors = ~Color,
                              line = list(color = "black", width = 1.5),
                              size = ~size),
                hoverinfo = "text",
                hovertext = paste("Symbol:", data$Symbol,
                                  "\nFBgn:", data$FBgn,
                                  "\n ", colnames(FDR.data)[1], "FDR: ", round(2^-data$FDR1, 2),
                                  "\n ", colnames(FDR.data)[2], "FDR: ", round(2^-data$FDR2,2),
                                  "\n ", colnames(FC.data)[1], "FC: ", round(data$FC1, 2),
                                  "\n ", colnames(FC.data)[2], "FC: ", round(data$FC2,2),
                                  "\n ", colnames(mean.data)[1], "mean CPM: ", signif(data$Variable1,2),
                                  "\n ", colnames(mean.data)[2], "mean CPM: ", signif(data$Control1,2),
                                  "\n ", colnames(mean.data)[3], "mean CPM: ", signif(data$Variable2,2),
                                  "\n ", colnames(mean.data)[4], "mean CPM: ", signif(data$Control2,2)))
  fig = fig %>% layout(title = main.title,
                       shapes = list(
                         type = "line",
                         line = list(color = "black", dash = 'dash', width = 2),
                         x0 = min(data$FC1, data$FC2) - 1,
                         x1 = max(data$FC1, data$FC2) + 1,
                         y0 = min(data$FC1, data$FC2) - 1,
                         y1 = max(data$FC1, data$FC2) + 1),
                       xaxis = list(title = paste(colnames(FC.data)[1], "FC")),
                       yaxis = list(title = paste(colnames(FC.data)[2], "FC")))
  
  
  fig
}

comparekegg.genes = function(specific.kegg, FC, FDR, mean.data){
  
  FC.data = FC
  FDR.data = FDR
  
  kegg.id = row.names(KEGG.Names)[c(grep(specific.kegg, KEGG.Names), 
                                    grep(specific.kegg, row.names(KEGG.Names)))]
  kegg.genes = names(GenesInKegg[grep(kegg.id, GenesInKegg)])
  
  kegg.fbgn = GeneIDKey[GeneIDKey$Symbol %in% kegg.genes,]
  kegg.fbgn = setNames(kegg.fbgn$FBgn, kegg.fbgn$Symbol)
  kegg.fbgn = na.omit(kegg.fbgn[kegg.genes])
  
  data = data.frame(Symbol = names(kegg.fbgn),
                    FBgn = kegg.fbgn,
                    FDR1 = -log2(FDR.data[kegg.fbgn, 1]),
                    FDR2 = -log2(FDR.data[kegg.fbgn, 2]),
                    FC1 = FC.data[kegg.fbgn, 1],
                    FC2 = FC.data[kegg.fbgn, 2],
                    Color = 'grey',
                    size = 1,
                    Variable1.cpm = mean.data[kegg.fbgn,1],
                    Control1.cpm = mean.data[kegg.fbgn,2],
                    Variable2.cpm = mean.data[kegg.fbgn,3],
                    Control2.cpm = mean.data[kegg.fbgn,4])
  
  main.title = KEGG.Names[kegg.id, 1]
  
  data$Color[(data$FC1 - 1) > data$FC2] = 'royalblue'
  data$Color[(data$FC2 - 1) > data$FC1] = 'firebrick' 
  data$size = 2.5*apply(data[,c('FDR1', 'FDR2')], MARGIN = 1, max)
  data$size[data$size < -2.5*log2(.05)] = 5
  data = na.omit(data)
  
  fig = plot_ly(data = data,
                x = ~FC1,
                y = ~FC2,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color,
                              colors = ~Color,
                              line = list(color = "black", width = 1.5),
                              size = ~size),
                hoverinfo = "text",
                hovertext = paste("Symbol:", data$Symbol,
                                  "\nFBgn:", data$FBgn,
                                  "\n ", colnames(FDR.data)[1], "FDR: ", round(2^-data$FDR1, 2),
                                  "\n ", colnames(FDR.data)[2], "FDR: ", round(2^-data$FDR2,2),
                                  "\n ", colnames(FC.data)[1], "FC: ", round(data$FC1, 2),
                                  "\n ", colnames(FC.data)[2], "FC: ", round(data$FC2,2),
                                  "\n ", colnames(mean.data)[1], "mean CPM: ", signif(data$Variable1,2),
                                  "\n ", colnames(mean.data)[2], "mean CPM: ", signif(data$Control1,2),
                                  "\n ", colnames(mean.data)[3], "mean CPM: ", signif(data$Variable2,2),
                                  "\n ", colnames(mean.data)[4], "mean CPM: ", signif(data$Control2,2)))
  fig = fig %>% layout(title = main.title,
                       shapes = list(
                         type = "line",
                         line = list(color = "black", dash = 'dash', width = 2),
                         x0 = min(data$FC1, data$FC2) - 1,
                         x1 = max(data$FC1, data$FC2) + 1,
                         y0 = min(data$FC1, data$FC2) - 1,
                         y1 = max(data$FC1, data$FC2) + 1),
                       xaxis = list(title = paste(colnames(FC.data)[1], "FC")),
                       yaxis = list(title = paste(colnames(FC.data)[2], "FC")))
  
  
  fig
}

compareparent.genes = function(parent.category, FC, FDR, mean.data){
  
  FC.data = FC
  FDR.data = FDR
  
  if(length(rt[rt$parentTerm == parent.category, "go"]) == 0){
    parent.category = rt$parentTerm[1]
  }
  
  ##select F if you want point size to scale with mean cpm of variable data
  set.size = F
  ##select T if you only want the genes in category to show
  catgenes.only = T
  
  GOterm.genes = vector('list', length = length(unique(rt$parent)))
  names(GOterm.genes) = unique(rt$parent)
  
  i=1
  while(i<=length(GOterm.genes)){
    goinparent = rt[grep(names(GOterm.genes)[i], rt$parent),'go']
    genesinparent = unlist(genesingo[goinparent])
    GOterm.genes[[i]] = GeneIDKey[GeneIDKey$ensembl %in% genesinparent, 'FBgn']
    i = i +1
  }
  
  parent.genes = GOterm.genes[[intersect(rt[rt$parentTerm == parent.category, "go"],
                                         names(GOterm.genes))]]
  parent.genes = intersect(parent.genes, row.names(FDR.data))
  
  
  data = data.frame(Symbol = GeneIDKey[parent.genes, "Symbol"],
                    FBgn = parent.genes,
                    FDR1 = -log2(FDR.data[parent.genes, 1]),
                    FDR2 = -log2(FDR.data[parent.genes, 2]),
                    FC1 = FC.data[parent.genes, 1],
                    FC2 = FC.data[parent.genes, 2],
                    Color = 'grey',
                    size = 1,
                    Variable1.cpm = mean.data[parent.genes,1],
                    Control1.cpm = mean.data[parent.genes,2],
                    Variable2.cpm = mean.data[parent.genes,3],
                    Control2.cpm = mean.data[parent.genes,4])
  
  main.title = parent.category
  
  data$Color[(data$FC1 - 1) > data$FC2] = 'royalblue'
  data$Color[(data$FC2 - 1) > data$FC1] = 'firebrick' 

  
  data$size = 2.5*apply(data[,c('FDR1', 'FDR2')], MARGIN = 1, max)
  data$size[data$size < -2.5*log2(.05)] = 5
  data = na.omit(data)
  
  fig = plot_ly(data = data,
                x = ~FC1,
                y = ~FC2,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color,
                              colors = ~Color,
                              line = list(color = "black", width = 1.5),
                              size = ~size),
                hoverinfo = "text",
                hovertext = paste("Symbol:", data$Symbol,
                                  "\nFBgn:", data$FBgn,
                                  "\n ", colnames(FDR.data)[1], "FDR: ", round(2^-data$FDR1, 2),
                                  "\n ", colnames(FDR.data)[2], "FDR: ", round(2^-data$FDR2,2),
                                  "\n ", colnames(FC.data)[1], "FC: ", round(data$FC1, 2),
                                  "\n ", colnames(FC.data)[2], "FC: ", round(data$FC2,2),
                                  "\n ", colnames(mean.data)[1], "mean CPM: ", signif(data$Variable1,2),
                                  "\n ", colnames(mean.data)[2], "mean CPM: ", signif(data$Control1,2),
                                  "\n ", colnames(mean.data)[3], "mean CPM: ", signif(data$Variable2,2),
                                  "\n ", colnames(mean.data)[4], "mean CPM: ", signif(data$Control2,2)))
  fig = fig %>% layout(title = main.title,
                       shapes = list(
                         type = "line",
                         line = list(color = "black", dash = 'dash', width = 2),
                         x0 = min(data$FC1, data$FC2) - 1,
                         x1 = max(data$FC1, data$FC2) + 1,
                         y0 = min(data$FC1, data$FC2) - 1,
                         y1 = max(data$FC1, data$FC2) + 1),
                       xaxis = list(title = paste(colnames(FC.data)[1], "FC")),
                       yaxis = list(title = paste(colnames(FC.data)[2], "FC")))
  
  
  fig
}



DEG.venn = function(FDR, FC){
  
  ##Sample titles as strings. Only fill in up to your number of selected categories
  set3=paste(colnames(FDR)[1], 'up')
  set1=paste(colnames(FDR)[1], 'down')
  set4=paste(colnames(FDR)[2], 'up')
  set2=paste(colnames(FDR)[2], 'down')
  
  ##items to be compared (ex: gene names) as a 1 dimensional vectors. Only give inputs up to your number of selected categories
  s3 = intersect(row.names(FDR[FDR[1] <.05,]), row.names(FC[FC[1] > 0,]))
  s1 = intersect(row.names(FDR[FDR[1] <.05,]), row.names(FC[FC[1] < 0,]))
  s4 = intersect(row.names(FDR[FDR[2] <.05,]), row.names(FC[FC[2] > 0,]))
  s2 = intersect(row.names(FDR[FDR[2] <.05,]), row.names(FC[FC[2] < 0,]))
  
  main = 
  
  grid.newpage()
  tempvenn=draw.quad.venn(area1=length(s1),area2=length(s2),area3=length(s3),area4=length(s4), n12=length(intersect(s1,s2)),n13 = length(intersect(s1,s3)),n14 = length(intersect(s1,s4)),n23 = length(intersect(s2,s3)),n24 = length(intersect(s2,s4)),n34 = length(intersect(s3,s4)),n123 =  length(intersect(s1,intersect(s2,s3))),n124 =  length(intersect(s1,intersect(s2,s4))), n134 = length(intersect(s1,intersect(s3,s4))),n234 = length(intersect(s2,intersect(s3,s4))),n1234 = length(intersect(s1,intersect(s2,intersect(s3,s4)))),category =c(set1,set2,set3,set4),fill = c("blue", "red","green","yellow"),cex=2,cat.cex = 1, title = "main")
}



compare.genes = function(FC, FDR, mean.data){
  
  FC.data = FC
  FDR.data = FDR
  
  data = data.frame(Symbol = GeneIDKey[row.names(FC.data), "Symbol"],
                    FBgn = row.names(FC.data),
                    FDR1 = -log2(FDR.data[row.names(FC.data), 1]),
                    FDR2 = -log2(FDR.data[row.names(FC.data), 2]),
                    FC1 = FC.data[, 1],
                    FC2 = FC.data[row.names(FC.data), 2],
                    Color = 'grey',
                    size = 1,
                    Variable1.cpm = mean.data[row.names(FC.data),1],
                    Control1.cpm = mean.data[row.names(FC.data),2],
                    Variable2.cpm = mean.data[row.names(FC.data),3],
                    Control2.cpm = mean.data[row.names(FC.data),4])
  
  main.title = ""
  
  data$Color[(data$FC1 - 1) > data$FC2] = 'royalblue'
  data$Color[(data$FC2 - 1) > data$FC1] = 'firebrick' 
  data$size = 2.5*apply(data[,c('FDR1', 'FDR2')], MARGIN = 1, max)
  data$size[data$size < -2.5*log2(.05)] = 5
  data = na.omit(data)
  
data = data[(data$FDR1 >= -log2(.05) | data$FDR2 >= -log2(.05)),]
  
  fig = plot_ly(data = data,
                x = ~FC1,
                y = ~FC2,
                type = 'scatter',
                mode = 'markers',
                marker = list(color = ~Color,
                              colors = ~Color,
                              line = list(color = "black", width = 1.5),
                              size = ~size),
                hoverinfo = "text",
                hovertext = paste("Symbol:", data$Symbol,
                                  "\nFBgn:", data$FBgn,
                                  "\n ", colnames(FDR.data)[1], "FDR: ", round(2^-data$FDR1, 2),
                                  "\n ", colnames(FDR.data)[2], "FDR: ", round(2^-data$FDR2,2),
                                  "\n ", colnames(FC.data)[1], "FC: ", round(data$FC1, 2),
                                  "\n ", colnames(FC.data)[2], "FC: ", round(data$FC2,2),
                                  "\n ", colnames(mean.data)[1], "mean CPM: ", signif(data$Variable1,2),
                                  "\n ", colnames(mean.data)[2], "mean CPM: ", signif(data$Control1,2),
                                  "\n ", colnames(mean.data)[3], "mean CPM: ", signif(data$Variable2,2),
                                  "\n ", colnames(mean.data)[4], "mean CPM: ", signif(data$Control2,2)))
  fig = fig %>% layout(title = main.title,
                       shapes = list(
                         type = "line",
                         line = list(color = "black", dash = 'dash', width = 2),
                         x0 = min(data$FC1, data$FC2) - 1,
                         x1 = max(data$FC1, data$FC2) + 1,
                         y0 = min(data$FC1, data$FC2) - 1,
                         y1 = max(data$FC1, data$FC2) + 1),
                       xaxis = list(title = paste(colnames(FC.data)[1], "FC")),
                       yaxis = list(title = paste(colnames(FC.data)[2], "FC")))
  
  
  fig
}


