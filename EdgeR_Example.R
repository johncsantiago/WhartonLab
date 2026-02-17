if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)

countdata=read.csv("https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/TKT_RNAseq/CountTables/TKT.counttable.csv",row.names=1)
groups=read.csv("https://raw.githubusercontent.com/johncsantiago/WhartonLab/refs/heads/master/TKT_RNAseq/CountTables/TKT.metadata.csv",row.names=1)

x = countdata[, row.names(groups)]
gr = factor(groups$Group)
y = DGEList(counts=x,group=gr)
keep = filterByExpr(y)
y = y[keep,,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM")

##normalized cpm
cpmdata=cpm(z)

design<-model.matrix(~0 + gr)
colnames(design) = unique(gr)
z = estimateDisp(z, design)
fit = glmQLFit(z,design)

compare = makeContrasts(GR.F-WT.F, levels=design)
lrt = glmQLFTest(fit,contrast=as.vector(compare))
G_X_E<-topTags(lrt,adjust.method="BH",n = nrow(z$counts), sort.by="PValue")
DEout=G_X_E$table
GRxWT.F2=DEout
sGRxWT.F2=DEout[DEout$FDR<.05,]



