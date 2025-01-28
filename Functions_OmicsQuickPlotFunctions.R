
##Load Libraries and Data

##library(plotly)

library(gplots)

##load all files from github
git.dir = "https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/GeneralDataFiles/"

GeneIDKey = read.csv(paste0(git.dir, "GeneIDKey.csv"), row.names = 1)

TKT.cpm = read.csv(paste0(git.dir,"TKT_cpmdata.csv"), row.names = 1)

TKT.groups  = read.csv(paste0(git.dir, "TKT.metadata.csv"), row.names = 1)

TKT.groups = TKT.groups[colnames(TKT.cpm),]


A4V.cpm = read.csv(paste0(git.dir, "A4V.cpmdata.csv"), row.names = 1)



gbbOE = read.csv(paste0(git.dir, "gbbOE_CountTable.csv"), row.names = 1)

gbbOE.groups = read.csv(paste0(git.dir, "gbbOE_Metadata.csv"), row.names = 1)

gbbOE = gbbOE[,row.names(gbbOE.groups)]
gbbOE.groups$group = paste0(gbbOE.groups$genotype, "_",gbbOE.groups$stage)

gbbKO.cpmdata = read.csv(paste0(git.dir, "gbbKO.cpmdata.csv"), row.names = 1)
gbbKO.groups = read.csv(paste0(git.dir, "gbbKO_Metadata.csv"), row.names = 1)


Metab.data = read.csv(paste0(git.dir, "RawMetabolomicData.csv"), row.names = 1)
Metab.Meta = read.csv(paste0(git.dir, "MetabolomicMetadata.csv"), row.names = 1)

KEGG.Key = setNames(Metab.data[-c(1:2), 2], row.names(Metab.data)[-c(1:2)])

norm.func = function(data){
  nd = (as.numeric(data)/Metab.Meta$TIC)*1000
  return(nd)
}

norm.data = t(apply(Metab.data[3:nrow(Metab.data),3:ncol(Metab.data)], 1, norm.func))

colnames(norm.data) = colnames(Metab.data)[3:ncol(Metab.data)]
colnames(norm.data) = Metab.Meta[colnames(norm.data), "Genotypes"]




##manually currated
serpe = read.csv("/Users/johncsantiago/Documents/GitHub/WhartonLab/Nguyen Serpe 2024 scRNA-seq VNC.csv", row.names = 1)

serpe.glia = read.csv(paste0(git.dir, "Nguyen%20Serpe%202024%20scRNA-seq%20VNC.csv"), row.names = 1)

serpe.mn = read.csv(paste0(git.dir, "Serpe%202024%20MNClusters.csv"), row.names = 1)

serpe.inter = read.csv(paste0(git.dir, "Serpe%202024%20InterneuronClusters.csv"), row.names = 1)

serpe.annotated.names = as.character(read.csv(paste0(git.dir, "Nguyen%20Serpe%202024%20scRNA-seq%20VNC.csv"), row.names = 1, header = F)[1,])

serpe.glia.names = as.character(read.csv(paste0(git.dir, "Serpe%202024%20GliaSubtypes.csv"), row.names = 1, header = F)[1,])



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
      
      
      
      gene = data.frame(cpm = as.numeric(TKT.cpm[FBgn,]), condition = TKT.groups$Group, group = paste0(TKT.groups$Genotype, TKT.groups$Sex))
      
      gene$group = factor(gene$group, levels = c("WT.f", "GR.F", "WT.M", "GR.M"))
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
               bg=c('firebrick',"gold","darkgreen","dodgerblue")[as.numeric(factor(gene$group))])
        
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
             cex.axis = 1)
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
        
        axis(side = 3,
             at = c(3.5,9.5),
             tick = T,
             outer = T,
             labels = F,
             line = -2.25,
             lwd.ticks = 0)
        
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
                main= name,ylab="Counts per Million Reads", las = 2, cex.axis = .8)
        points(x=as.numeric(factor(gene$condition)),y=gene$cpm,cex=1,pch=21,bg=c('firebrick',"gold","darkgreen","dodgerblue")[as.numeric(factor(gene$group))])
        legend('topright',
               inset=c(-0.275,0),
               legend = c('Silent F', 'A4V F', 'Silent M', 'A4V M', 'Silent M'), 
               fill = c('firebrick', "gold", "darkgreen", "dodgerblue"),
               cex = .65,
               bty = 'n',
               pt.cex = .5)
      }
    }
    
    if(nrow(ID.row) != 1
       | length(intersect(GeneIDKey$FBgn, ID)) > 1
       | length(intersect(GeneIDKey$Symbol, ID)) > 1){
      ID.row
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
                         geno = colnames(norm.data))
      
      
      metab$group = factor(metab$group, levels = c("WT", "TktDfWT", "TktOEWT",
                                                   "GR", "TktDfGR", "TktOEGR"))
      
      metab$geno = factor(substr(metab$group, (nchar(as.character(metab$group)) -1), nchar(as.character(metab$group))), 
                          levels = c("WT", "GR"))
      
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
      norm.data[ID.row,]
    }
    
  }
}

##Function to plot a specific gene expression levels in all cell types for the SC data
plot.SC = function(gene.name){
  ID.row = SCdata[grep(gene.name, SCdata$Gene.name),]
  
  if(nrow(ID.row) < 50 & nrow(ID.row) > 0){
    if(nrow(ID.row) == 1 | length(intersect(gene.name, SCdata$Gene.name)) == 1){
      
      if(length(intersect(gene.name, SCdata$Gene.name)) == 1){
        gene.levels = setNames(as.numeric(ID.row[ID.row$Gene.name == gene.name, colnames(ID.row) != "Gene.name"]), colnames(ID.row)[colnames(ID.row) != "Gene.name"]) 
      }
      
      if(nrow(ID.row) == 1){
        gene.name = ID.row$Gene.name
        gene.levels = setNames(as.numeric(ID.row[,colnames(ID.row) != "Gene.name"]), colnames(ID.row)[colnames(ID.row) != "Gene.name"]) 
      }
      
      gene.levels = gene.levels[18:1]
      highlight = rep(0, 18)
      highlight[(gene.levels > 2*(mean(gene.levels)))] = 1
      highlight = gene.levels * highlight
      
      nochange = gene.levels
      nochange[gene.levels > 1.5] = 0
      nochange = gene.levels * nochange
      
      label.colors = highlight
      
      label.colors[label.colors > 0] = "deeppink"
      label.colors[label.colors == 0] = "black"
      label.colors[nochange > 0] = 'steelblue'
      
      if(1.25*(max(gene.levels)) < 5){
        xmax = 5
      }
      
      if(1.25*(max(gene.levels)) >= 5){
        xmax = 1.25*(max(gene.levels))
      }
      
      
      par(mar=c(5, 12, 4, 6), xpd=TRUE)
      
      plot(x = NA,
           y = NA,
           type = 'n',
           ylim = c(.95, 21),
           xlim = c(0, xmax),
           xaxt = "n",
           ylab="",
           yaxt = "n",
           xlab = "Mean UMI Counts",
           bty = "n",
           main = gene.name)
      
      barplot(height = gene.levels, 
              horiz = T, 
              #las = 2,
              xlim = c(0, xmax),
              names.arg = NA,
              border = T,
              add = T)
      
      barplot(height = highlight, 
              horiz = T, 
              #las = 2, 
              names.arg = NA,
              border = T,
              col = "deeppink",
              add = T)
      
      barplot(height = nochange, 
              horiz = T, 
              #las = 2, 
              names.arg = NA,
              border = T,
              col = "lightsteelblue",
              add = T)
      
      mtext(side = 2,
            las = 2,
            at = .7+(1.2*(0:17)),
            text = names(gene.levels),
            col = label.colors,
            line = 0)
      
      abline(v = 1,
             lty = 2,
             lwd = 2,
             xpd = F)
      
      legend(x = .9*xmax,
             y = 23,
             legend = c("high", "mid", "low"),
             fill = c("deeppink", "grey", "lightsteelblue"),
             bty = "n",
             xpd = T)
      
    }
    
    if(nrow(ID.row) != 1 | length(intersect(gene.name, SCdata$Gene.name)) > 1){
      ID.row
    }
  }
}

##Function to plot a specific gene expression levels in the gbbOE transcriptomic data set
plot.gbbOE = function(ID){
  
  ID.row = GeneIDKey[grep(ID, GeneIDKey$FBgn),]
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
      
      gene = data.frame(cpm = as.numeric(gbbOE.cpmdata[FBgn,]), condition = gbbOE.groups$group, group = gbbOE.groups$genotype)
      
      gene$group = factor(gene$group, levels = c(unique(gene$group)))
      gene$condition = factor(gene$condition, levels = c('LoxP_ML3', "LoxP_LL3",
                                                         'G85R_ML3', "G85R_LL3",
                                                         'G85R_gbb_ML3', "G85R_gbb_LL3"))
      if(nrow(na.omit(gene)) > 0){
      boxplot(gene$cpm~gene$condition, xlab="",     
              main= name,ylab="Counts per Million Reads", las = 2, cex.axis = .8)
      points(x=as.numeric(factor(gene$condition)),y=gene$cpm,cex=1,pch=21,bg=c('firebrick',"dodgerblue","gold")[as.numeric(factor(gene$group))])
      }
    }
    
    if(nrow(ID.row) != 1){
      ID.row
    }
  }
}

##Function to plot a specific gene expression levels in the gbbKO transcriptomic data set
plot.gbbKO = function(ID){
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
      
      gene = data.frame(cpm = as.numeric(gbbKO.cpmdata[FBgn,]), condition = gbbKO.groups$Genotype)
      
      gene$condition = factor(gene$condition, levels = c('WT', "KO"))
      
      if(nrow(na.omit(gene)) > 0){
      boxplot(gene$cpm~gene$condition, xlab="",     
              main= name,ylab="Counts per Million Reads", las = 2, cex.axis = .8)
      points(x=as.numeric(factor(gene$condition)),y=gene$cpm,cex=1,pch=21,bg=c('firebrick',"dodgerblue","gold")[as.numeric(factor(gene$condition))])
      }
    }
    
    if(nrow(ID.row) != 1){
      ID.row
    }
  }
}

##Function to plot a specific gene expression levels in all cell types for the SC data
plot.serpe.allclusters = function(gene.name){
  
  ID.row = serpe[grep(gene.name, row.names(serpe)),]
  
  if(nrow(ID.row) < 50 & nrow(ID.row) > 0){
    if(nrow(ID.row) == 1 | length(intersect(gene.name, row.names(serpe))) == 1){
      
      
      serpe.annotated.gene = data.frame(level = as.numeric(serpe[gene.name,]),
                                        cluster = serpe.annotated.names,
                                        color = "darkgreen")
      
      
      serpe.mn.gene = data.frame(level = as.numeric(serpe.mn[gene.name,]),
                                 cluster = colnames(serpe.mn),
                                 color = "brown")
      
      
      serpe.inter.gene = data.frame(level = as.numeric(serpe.inter[gene.name,]),
                                    cluster = colnames(serpe.inter),
                                    color = "gold")
      
      
      serpe.glia.gene = data.frame(level = as.numeric(serpe.glia[gene.name,]),
                                   cluster = serpe.glia.names,
                                   color = "steelblue")
      
      serpe.all = rbind(serpe.annotated.gene,
                        rbind(serpe.mn.gene,
                              rbind(serpe.inter.gene,
                                    serpe.glia.gene)))
      
      serpe.all$color[serpe.all$level < 1.5] = "darkgrey"
      
      xmax = max(serpe.all$level)
      if(xmax < 10){
        xmax = 10
      }
      
      par(mar = c(4, 13, 4, 8), xpd=TRUE)
      
      plot(x = NA,
           y = NA,
           type = 'n',
           ylim = c(4.5, 120),
           xlim = c(0, xmax),
           xaxt = "n",
           ylab="",
           yaxt = "n",
           xlab = "Mean Expression Levels",
           bty = "n",
           main = gene.name)
      
      barplot(height = serpe.all$level, 
              horiz = T, 
              #las = 2,
              xlim = c(0, xmax),
              names.arg = NA,
              border = T,
              col = serpe.all$color,
              add = T)
      
      mtext(side = 2,
            las = 2,
            at = .7+(1.2*(0:104)),
            text = serpe.all$cluster,
            col = serpe.all$color,
            line = 0)
      
      legend(x = xmax,
             y = 125,
             pch = 21,
             legend = c("Glia", "Interneurons", "Motor Neurons", "Manually Curated"),
             pt.bg = c("steelblue", "gold", "brown", "darkgreen"),
             bty = "n",
             xpd = T,
             pt.cex = 1.5,
             y.intersp = 1)
      
      rect(0 - (xmax/40), 0, 0 - (xmax/100), 18*1.2, col = 'darkgreen')
      rect(0 - (xmax/40), 18*1.2, 0 - (xmax/100), 46*1.2, col = 'brown')
      rect(0 - (xmax/40), 46*1.2, 0 - (xmax/100), 95*1.2, col = 'gold')
      rect(0 - (xmax/40), 95*1.2, 0 - (xmax/100), 105*1.2, col = 'steelblue')
      
    }

  }
}

##Function to plot expression levels for a specified set of genes in all cell types for the SC data
serpe.heatmap = function(heat.genes){
  heat.data = na.omit(cbind(serpe[heat.genes,],
                            cbind(serpe.mn[heat.genes,],
                                  cbind(serpe.inter[heat.genes,], serpe.glia[heat.genes,]))))
  
  heat.data = apply(heat.data, MARGIN = 2, as.numeric)
  #heat.data = log2(heat.data)
  row.names(heat.data) = heat.genes
  
  #heat.data[heat.data == 0] = NA
  
  Cluster.id = c(serpe.annotated.names,
                 colnames(serpe.mn),
                 colnames(serpe.inter),
                 serpe.glia.names)
  
  Cluster.labelcolor = c(rep("darkgreen", length(serpe.annotated.names)),
                         rep("brown", ncol(serpe.mn)),
                         rep("orange", ncol(serpe.inter)),
                         rep("steelblue", length(serpe.glia.names)))
  
  
  heat.colors = colorRampPalette(c("grey95", "steelblue", "gold", "red3"))
  
  numcols = 100
  
  max.break = log2(max(heat.data))
  
  if(max.break > 7){
    max.break = 7
  }
  
  color.breaks = 2^(c(0:(numcols))*(max.break/(numcols + 1)))
  
  heatmap.2(heat.data,
            Colv = NA,
            #Rowv = NA,
            labCol = Cluster.id,
            scale = 'none',
            ColSideColors = Cluster.labelcolor,
            trace = "none",
            col = heat.colors(numcols),
            dendrogram = "none",
            na.color = "white",
            key.title = "Mean Expression",
            colCol = Cluster.labelcolor,
            margins = c(10, 5), 
            colsep = c(18, 46, 95),
            sepcolor = "black",
            breaks = color.breaks)
  legend("left",
         title = "Cell Types",
         legend = c("Manually Currated", "Motor Neurons", "Interneurons", "Glia"),
         fill = c("darkgreen", "brown", "orange", "steelblue"),
         xpd = T)
  
}

