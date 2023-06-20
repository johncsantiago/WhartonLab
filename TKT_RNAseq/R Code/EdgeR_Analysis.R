if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)

##Clean Trait Data and calculate CPMs
countdata=read.csv("/Users/johncsantiago/Documents/TKT RNAseq/CountTables/TKT.counttable.csv",row.names=1)
groups=read.csv("/Users/johncsantiago/Documents/TKT RNAseq/CountTables/TKT.metadata.csv",row.names=1)

##GR.F1.US = read.csv("/Users/johncsantiago/Downloads/GR-F1_CountTable_UnsortedBAM.txt", sep = "\t", row.names = 1, header = F)
##countdata$GR.F1.US = GR.F1.US [row.names(countdata),1]
##groups['GR.F1.US',] = c("GR.F1.US", "GR", "F", "4", "GR.F")

##normalize data
countdata=countdata[, row.names(groups)]
x <- countdata
gr <- factor(groups$Groups)
y <- DGEList(counts=x,group=gr)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM")
##normalized cpm
cpmdata=cpm(z)

plotMDS(z, col=as.numeric(factor(gr)))

##updated edger analysis pipeline
design<-model.matrix(~0 + gr)
colnames(design) = unique(gr)
z = estimateDisp(z, design)
fit <- glmQLFit(z, design)


##compare = makeContrasts(Lactate, levels=design)
##lrt = glmQLFTest(fit,contrast=as.vector(compare))

##lrt <- glmLRT(fit,contrast=as.vector(compare))		
##G_X_E<-topTags(lrt,adjust.method="BH",n = nrow(z$counts), sort.by="PValue")
##AFLF=G_X_E$table
##sigsAFLF=AFLF[AFLF$FDR<.05,]