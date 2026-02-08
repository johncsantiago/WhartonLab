# Load data from GitHub
git.dir <- "https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/GeneralDataFiles/"

G85R.meancpm = read.csv(paste0(git.dir, "TKT_meancpmdata.csv"), row.names = 1)
G85R.meancpm.meta = read.csv(paste0(git.dir, "TKT_meancpmdatameta.csv"), row.names = 1)

raw.nodes = read.csv(paste0(git.dir, "combined.nodes.csv"), row.names = 1)
raw.nodes$x = raw.nodes$x*1
raw.nodes$y = raw.nodes$y*1

raw.edges = read.csv(paste0(git.dir, "combined.edges.csv"))



GeneIDKey = read.csv(paste0(git.dir, "GeneIDKey.csv"), row.names = 1)

TKT.cpm = read.csv(paste0(git.dir,"TKT_cpmdata.csv"), row.names = 1)

TKT.groups  = read.csv(paste0(git.dir, "TKT.metadata.csv"), row.names = 1)

TKT.groups = TKT.groups[colnames(TKT.cpm),]

TKT.FC = read.csv(paste0(git.dir, "TKT.EdgeR.FCTable2.csv"), row.names = 1)
G85R.geneFC = TKT.FC

TKT.FDR = read.csv(paste0(git.dir, "TKT.EdgeR.FDRTable2.csv"), row.names = 1)
G85R.geneFDR = TKT.FDR

A4V.cpm = read.csv(paste0(git.dir, "A4V.cpmdata.csv"), row.names = 1)

A4V.FC = read.csv(paste0(git.dir, "A4V.FCdata.csv"), row.names = 1)

A4V.FDR = read.csv(paste0(git.dir, "A4V.FDRdata.csv"), row.names = 1)

Metab.data <- read.csv(paste0(git.dir, "RawMetabolomicData.csv"), row.names = 1)
Metab.Meta <- read.csv(paste0(git.dir, "MetabolomicMetadata.csv"), row.names = 1)
Metab.FC <- read.csv(paste0(git.dir, "MetaboliteFCs2.csv"), row.names = 1)
Metab.FDR <- read.csv(paste0(git.dir, "MetaboliteFDRs2.csv"), row.names = 1)
G85R.metabFDR = Metab.FDR
G85R.metabFC = Metab.FC 

Metab.FC = Metab.FC[,-1]
Metab.FDR = Metab.FDR[,-1]

temp = colorRampPalette(c("navy", "royalblue","white", "red", "firebrick"))(401)

temp2 = colorRampPalette(c("navy", "royalblue","grey80", "red", "firebrick"))(401)

hex = function(color, alpha){
  rgb2hex = color
  i = 1
  while(i <= length(color)){
    rgbvals = col2rgb(color[i])
    rgb2hex[i] = rgb(rgbvals[1]/255, rgbvals[2]/255, rgbvals[3]/255, alpha)
    i = i+1
  }
  return(rgb2hex)
}

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


METAB_CHOICES <- c(
  "GRxWT.C","GRxWT.Df","GRxWT.OE","GR.CxDf","WT.CxDf",
  "GR.CxOE","WT.CxOE","GR.DfxWT.C","GR.OExWT.C"
)

ENZYME_CHOICES <- c(
  "GRF.CxDf","GRF.CxOE","WTF.CxDf","WTF.CxOE","GRxWT.FC","GRxWT.FDf",
  "GRxWT.FOE","GRDfxWT.F","GROExWT.F","GRM.CxDf","GRM.CxOE","WTM.CxDf",
  "WTM.CxOE","GRxWT.MC","GRxWT.MDf","GRxWT.MOE","GRDfxWT.M","GROExWT.M",
  "GR.FxM","WT.FxM"
)

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
               cex=1.5,
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

fullnetwork = function(metab, 
                       enzyme,
                       fccutoff = 0, 
                       meansumcutoff = 0){
  nodes = raw.nodes
  edges = raw.edges
  mean.groups = row.names(G85R.meancpm.meta)[G85R.meancpm.meta[,enzyme]]
  
  nodes[nodes$shape == "dot", "label"] = NA
  
  nodes.color = G85R.metabFC[nodes[nodes$shape == "circle", "data.id"], metab]
  
  nodes$FC[nodes$shape == "circle"]= G85R.metabFC[nodes[nodes$shape == "circle", "data.id"], metab]
  nodes$FDR[nodes$shape == "circle"]= G85R.metabFDR[nodes[nodes$shape == "circle", "data.id"], metab]
  
  node.size = G85R.metabFC[nodes[nodes$shape == "circle", "data.id"], metab]
  node.size[na.action(na.exclude(node.size))] = 0
  node.size = ((abs(node.size)/11)*80)+10
  node.size[node.size>100] = 100
  nodes[nodes$shape == "circle", "size"] = node.size
  
  nodes$color.background[nodes$shape == "circle"] = "white"
  nodes$color.background[nodes$shape == "circle" &
                           nodes$FC > 0] = "lightpink"
  nodes$color.background[nodes$shape == "circle" &
                           nodes$FC < 0] = "lightblue"
  nodes$color.background[nodes$shape == "circle" &
                           nodes$FDR <= 0.05 &
                           nodes$FC > 0] = "firebrick"
  nodes$color.background[nodes$shape == "circle" &
                           nodes$FDR <= 0.05 &
                           nodes$FC < 0] = "royalblue"
  
  nodes$color.background[nodes$shape == "circle" &
                           abs(nodes$FC) <= fccutoff] = "grey"
  
  nodes$color.border = "black"
  nodes$borderWidth = 1
  nodes$shape = 'dot'
  nodes$font.background = "white"
  
  nodes$title = paste0(nodes$id,
                       "<br>-log2(FC): ", 
                       signif((G85R.metabFC[nodes$data.id, metab]), 3),
                       "<br>FDR: ", 
                       signif((G85R.metabFDR[nodes$data.id, metab]), 3))
  
  
  edges$title = paste0(GeneIDKey[edges$FBgn, "Symbol"],
                       "<br>-log2(FC): ", 
                       signif(G85R.geneFC[edges$FBgn,
                                          enzyme],3),
                       "<br>FDR: ", 
                       signif(G85R.geneFDR[edges$FBgn,
                                           enzyme],3),
                       "<br>", mean.groups[1]," mean cpm: ", 
                       signif(G85R.meancpm[edges$FBgn,
                                           mean.groups[1]],3),
                       "<br>", mean.groups[2]," mean cpm: ", 
                       signif(G85R.meancpm[edges$FBgn,
                                           mean.groups[2]],3))
  edges.meansum = G85R.meancpm[edges$FBgn,mean.groups]
  edges.meansum=apply(edges.meansum, 1, sum)
  
  edges.FC = G85R.geneFC[edges$FBgn, enzyme]
  
  edges.width = G85R.geneFC[edges$FBgn, enzyme]
  edges.width[na.action(na.exclude(edges.width))] = 0
  edges.width[apply(G85R.meancpm[edges$FBgn, row.names(G85R.meancpm.meta)[G85R.meancpm.meta[,enzyme]]], 1, sum) < 5] = 0
  edges.width = (2^abs(edges.width))-.5
  edges.width[edges.width>8] = 8
  edges$width = edges.width
  
  edges.sig = G85R.geneFDR[edges$FBgn, enzyme]
  edges.direction = G85R.geneFC[edges$FBgn, enzyme]
  edges.color = "black"
  edges.color[edges.sig <= .05 &
                edges.direction > 0] = "firebrick"
  edges.color[edges.sig > .05 &
                edges.direction > 0] = "lightpink"
  edges.color[edges.sig > .05 &
                edges.direction < 0] = "lightblue"
  edges.color[edges.sig <= .05 &
                edges.direction < 0] = "royalblue"
  edges$color = edges.color
  
  edges$color[abs(edges.FC) < fccutoff] = "grey"
  edges$color[edges.meansum < meansumcutoff] = "white"
  edges$color[na.action(na.exclude(edges.meansum))] = "black"
  edges = edges[edges$color != "white",]
  
  if(length(unique(G85R.meancpm.meta[mean.groups,"Sex"])) == "2"){
    sex.specific = G85R.meancpm[,mean.groups]
    sex.specific$sex.specific = F
    sex.specific[sex.specific[,1]<meansumcutoff/2 &
                   sex.specific[,2]>=meansumcutoff/2, "sex.specific"] = T
    sex.specific[sex.specific[,1]>=meansumcutoff/2 &
                   sex.specific[,2]<meansumcutoff/2, "sex.specific"] = T
    sex.specific = sex.specific[edges$FBgn,]
    edges$color[sex.specific$sex.specific] = "green"
  }
  
  
  
  
  visNetwork(nodes, edges, background = "white", main = enzyme)%>%
    visNodes(physics = F)%>%
    #visLayout(hierarchical=F, improvedLayout=T)%>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = list(enabled = T, values = nodes[nodes$size > 0, "label"]))%>%
    visNodes(font = list(color = "black", size = 10, background = "white"))%>%
    
    visEdges(arrows = edges$arrows, physics = T)
}
