library(scran)
library(scater)
sce = readRDS("Sertoli7to150.rds")
rownames(sce) = rowData(sce)$Symbol
sce = addPerCellQC(sce)

n.counts.idx = (colData(sce)$sum >2000) & (colData(sce)$sum < 100000)
sce.filtered.nc.1 = sce[,n.counts.idx]

n.detected.idx = (sce.filtered.nc.1$detected > 200) & (sce.filtered.nc.1$detected <4000)
sce.filtered.nd.2 = sce.filtered.nc.1[,n.detected.idx]

plot(sce.filtered.nd.2$detected, sce.filtered.nd.2$total, pch=16,xlab="ngenes", ylab="umis", log="xy")

plot(log10(sce$detected), log10(sce$total), pch=16,xlab="log10ngenes", ylab="log10umis", log="xy")

sce.filtered.nd.2$log10genesperumi = log10(sce.filtered.nd.2$detected)/log10(sce.filtered.nd.2$total)

plot(log10(sce.filtered.nd.2$detected), log10(sce.filtered.nd.2$total), pch=16,xlab="ngenes", ylab="umis", log="xy")

sce.filtered.nd.2 = addPerFeatureQC(sce.filtered.nd.2)
ncell.idx = rowData(sce.filtered.nd.2)$detected*3138/100 > 3
sce.filtered.ncell.3 = sce.filtered.nd.2[ncell.idx,]

scatter.smooth(sce.filtered.nd.2$log10genesperumi)





rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

#library(BSgenome.Sscrofa.UCSC.susScr3)
#location <- mapIds(BSgenome.Sscrofa.UCSC.susScr3, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")

sce = sce.filtered.ncell.3

set.seed(1000)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce<- logNormCounts(sce)

plot(librarySizeFactors(sce), sizeFactors(sce), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")

set.seed(1001)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)

plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)


set.seed(10000)
sce <- denoisePCA(sce, subset.row=top, technical=dec)

set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by="Sample")

sce = runPCA(sce)
sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by="Sample")

library(batchelor)
#d <- 32
sce$Study = to_vec(for (x in sce$Sample) (if(grepl("150",x)) "S2" else "S1"))
FMNN.out <-  fastMNN( sce  , batch=sce$Study , use.dimred="PCA" ) 
reducedDim (sce, "PCA.FMNN" ) <- FMNN.out$corrected 
sce <- runTSNE(sce, dimred="PCA.FMN")
plotTSNE(sce, colour_by="Sample")



clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)

reducedDim(z, "PCA" ) = rsvd::rpca(t( assay(z,"logcounts") ),k=32,retx=TRUE,center=TRUE,scale=FALSE)$x
mBN <- multiBatchNorm(S1 = z[,z$Study=="S1"],S2 = z[,z$Study=="S2"])
z <- cbind( mBN$S1, mBN$S2) 

mB.PCA <- multiBatchPCA( z, batch=z$Study, d=32, preserve.single = TRUE)
reducedDim(z , "PCA" )  <- mB.PCA[[1]]

z <- runTSNE(z, dimred="PCA")
plotTSNE(z, colour_by="Sample")

d <- 32
FMNN.out <- fastMNN( z  , batch=z$Sample , d=d ) 
reducedDim (z, "PCA.FMN" ) <- reducedDim(FMNN.out,"corrected")
z <- runTSNE(z, dimred="PCA.FMN")
plotTSNE(z, colour_by="Sample")


set.seed(1000000)
sce <- runUMAP(sce, dimred="PCA")


g <- buildSNNGraph(sce, k=50, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce) <- factor(clust)
table(colLabels(sce))
plotTSNE(sce, colour_by="Sample")


#markers <- findMarkers(sce, pval.type="some", direction="up")
plotExpression(sce, features=c("AMH", "SOX9",
                                    "CDH2", "CDH3","AR","FSHR"), x="Sample", colour_by="Sample")



mart.ss = useEnsembl("ensembl","sscrofa_gene_ensembl")

location = getBM(attributes = "chromosome_name", filters = "ensembl_gene_id", 
                 values = rowData(sce)$ID, mart = mart.ss)
