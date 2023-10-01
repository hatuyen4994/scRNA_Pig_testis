markers = list(
  Spermatogonia=c("UCHL1", "CTCFL"),
  Spermatocytes=c("SYCP1", "PIWIL1"),
  Round.STids=c("ACRV1"),
  Elongating.STids=c("SPACA1"),
  Immature.Sperm=c("TNP1", "PRM2"),
  Sertoli=c("CLU", "AMH"),
  Leydig=c("INSL3", "STAR"),
  Myoid=c("MYH11","ACTA2"), 
  Endothelial=c("VWF","PECAM1"), 
  Macrophage=c("CD68","C1QA"), 
  T.cells=c("CD4", "IL7R")
)



setwd("D:/R/scRNA_Pig_testis/")
library(SingleCellExperiment)
library(DropletUtils)
set.seed(1000)
sce = read10xCounts("PND7/")


clusters = quickCluster(sce, assay.type="X")
sce = computeSumFactors(sce,cluster=clusters, assay.type="X")
sce = logNormCounts(sce, assay.type="X")
plot(librarySizeFactors(sce, assay.type="X"), sizeFactors(sce)
     ,pch=16,xlab="Library", ylab="Decon",log="xy")




###
dec = modelGeneVarByPoisson(sce, assay.type="logcounts")
top = getTopHVGs(dec,prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)


sce = scater::runPCA(sce,subset_row=top)
sce = scater::runTSNE(sce, dimred="PCA")
sce = scater::runUMAP(sce, dimred="PCA")


g = buildSNNGraph(sce, k=10, use.dimred="PCA")
clust = igraph::cluster_walktrap(g)$membership
colLabels(sce) = factor(clust)
table(colLabels(sce))
scater::plotTSNE(sce, colour_by="label")

###
sce = denoisePCA(sce,subset.row=top, technical=dec)
sce = scater::runTSNE(sce, dimred="PCA")
sce = scater::runUMAP(sce, dimred="PCA")


g = buildSNNGraph(sce, k=10, use.dimred="PCA")
clust = igraph::cluster_walktrap(g)$membership
colLabels(sce) = factor(clust)
table(colLabels(sce))
scater::plotTSNE(sce, colour_by="label")


###
dec = modelGeneVarByPoisson(sce, assay.type="X")
top = getTopHVGs(dec,prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)


sce = scater::runPCA(sce,subset_row=top)
sce = scater::runTSNE(sce, dimred="PCA")
sce = scater::runUMAP(sce, dimred="PCA")


g = buildSNNGraph(sce, k=10, use.dimred="PCA")
clust = igraph::cluster_walktrap(g)$membership
colLabels(sce) = factor(clust)
table(colLabels(sce))
scater::plotTSNE(sce, colour_by="label")
