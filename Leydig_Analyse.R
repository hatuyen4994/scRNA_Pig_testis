library(scran)
library(scater)
library(comprehenr)
#sce = readRDS("Sertoli7to150.rds.rds")
sce = readRDS("Leydig7to150_2markers.rds")
snames = sce$Sample
snames = to_vec( for (x in snames) (if (x=="PND7") "PND007" else x))
snames = to_vec( for (x in snames) (if (x=="PND30") "PND030" else x))
snames = to_vec( for (x in snames) (if (x=="PND60") "PND060" else x))
snames = to_vec( for (x in snames) (if (x=="PND90") "PND090" else x))
snames = to_vec( for (x in snames) (if (x=="PND150/") "PND150" else x))
sce$Sample = snames
rownames(sce) = rowData(sce)$Symbol

##QC
sce = addPerFeatureQC(sce)
ncells.idx = rowData(sce)$detected*dim(sce)[2]/100 > 3
sce = sce[ncells.idx,]

ID_MT = read.table("ID_Chromosome.tsv")
DF = DataFrame(
  chromosome=ID_MT$V3,
  row.names=ID_MT$V1
)
rowData(sce)$location = DF[rowData(sce)$ID,]
is.mito = grepl("MT", rowData(sce)$location)
sce=addPerCellQCMetrics(sce, subsets=list(Mito=is.mito))

mito.idx <- sce$subsets_Mito_percent < 20
sce=sce[,mito.idx]
umis.idx = (colData(sce)$sum >2000) & (colData(sce)$sum < 100000)
sce = sce[,umis.idx]
ngenes.idx = (sce$detected > 200) & (sce$detected <4000)
sce = sce[,ngenes.idx]
sce$log10genesperumi = log10(sce$detected)/log10(sce$total)

#ensembl_ID to orthologue genes
df.in=read.table("Ley.ENS2Name.tsv", header = FALSE, sep='\t')
dim(df.in)
to.convert = df.in$V1
gene.names = df.in$V2
rowData(sce[to.convert])$Symbol = gene.names
rownames(sce) = rowData(sce)$Symbol

#Vis
gridExtra::grid.arrange(
  plotColData(sce, x="Sample", y="sum", colour_by = "Sample") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="Sample", y="detected", colour_by = "Sample") +
    ggtitle("Detected genes"),
  plotColData(sce, x="Sample", y="subsets_Mito_percent", colour_by = "Sample") 
  + ggtitle("Mito percent"),
  ncol=1
)


sce.corrected = sce
plot(log10(sce$detected), log10(sce$total), pch=16,xlab="log10ngenes", ylab="log10umis", log="xy")

######Uncorrected######
d=NULL
set.seed(1000)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce<- logNormCounts(sce)

#normal PCA
set.seed(1001)
sce <- runPCA(sce)
plotPCA(sce, colour_by="Sample")

set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by="Sample")

sce <- runUMAP(sce, dimred="PCA")
plotUMAP(sce, colour_by="Sample")

#denoise
set.seed(1001)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)
set.seed(10000)
sce <- denoisePCA(sce, subset.row=top, technical=dec)
plotPCA(sce, colour_by="Sample")

set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by="Sample")

sce <- runUMAP(sce, dimred="PCA")
plotUMAP(sce, colour_by="Sample")

####Corrected#####
library(batchelor)
sce.corrected$Study = to_vec(for (x in sce.corrected$Sample) (if(grepl("150",x)) "S2" else "S1"))

d=50
clusters <- quickCluster(sce.corrected)
sce.corrected <- computeSumFactors(sce.corrected, cluster=clusters)
mBN <- multiBatchNorm(S1 = sce.corrected[,sce.corrected$Study=="S1"],S2 = sce.corrected[,sce.corrected$Study=="S2"])
sce.corrected <- cbind( mBN$S1, mBN$S2) 
mB.PCA <- multiBatchPCA( sce.corrected, batch=sce.corrected$Study, d=d, preserve.single = TRUE)
reducedDim(sce.corrected , "PCA" )  <- mB.PCA[[1]]
FMNN.out <- fastMNN( sce.corrected  , batch=sce.corrected$Study, d=d ) 
reducedDim (sce.corrected, "PCA.FMNN" ) <- reducedDim(FMNN.out,"corrected")

test = sce.corrected
reducedDim(test, "PCA") = reducedDim (test, "PCA.FMNN" )
plotPCA(test, colour_by="Sample")


sce.corrected <- runTSNE(sce.corrected, dimred="PCA.FMNN")
plotTSNE(sce.corrected, colour_by="Sample")

sce.corrected <- runUMAP(sce.corrected, dimred="PCA.FMNN")
plotUMAP(sce.corrected, colour_by="Sample")


###Visual###
plotExpression(sce.corrected, features=c("CYP11A1","CYP17A1","HSD3B1","HSD17B3", 
                                         "INSL3","STAR","NR5A1","LHCGR"), x="Sample", colour_by="Sample")

plotExpression(sce.corrected, features=c("CDH1", "CDH2","CDH3", "GJA1",
                                         "CLDN11", "CTNNB1","TJP1","OCLN"), x="Sample", colour_by="Sample")

cdh = rownames(sce)[grep("^CDH+[0-9]",rownames(sce))]
cdh = cdh[rowData(sce.corrected[cdh])$detected > 2]
plotExpression(sce.corrected, features=cdh, x="Sample", colour_by="Sample")

cdh = rownames(sce)[grep("^CDH+[0-9]",rownames(sce))]
cdh = cdh[(rowData(sce.corrected[cdh])$detected > 1)&(rowData(sce.corrected[cdh])$detected < 2)]
plotExpression(sce.corrected, features=cdh, x="Sample", colour_by="Sample")

gja = rownames(sce.corrected)[grepl("GJA",rownames(sce))]
plotExpression(sce.corrected, features=gja, x="Sample", colour_by="Sample")

sox = rownames(sce.corrected)[grepl("^SOX",rownames(sce))]
sox = sox[(rowData(sce.corrected[sox])$detected > 1)]
plotExpression(sce.corrected, features=sox, x="Sample", colour_by="Sample")

plotExpression(sce.corrected, 
               features=c("JUN","JUNB","ENSSSCG00000031657" ,"FOS", "FOSB", 
                          "FOSL1","FOSL2","CREB1" ), x="Sample", colour_by="Sample")

plotExpression(sce.corrected, 
               features=c("NR4A1","NR2F2","CEBPB","SP1","INSL3",
                          "PCNA","MKI67","CCND1"), x="Sample", colour_by="Sample")


####Markers####
colLabels(sce.corrected) = factor(sce.corrected$Sample)
marker.info <- scoreMarkers(sce.corrected, colLabels(sce.corrected))

chosen <- marker.info[[1]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4]) # showing basic stats only, for brevity.
plotExpression(sce.corrected, features=head(rownames(ordered)), 
               x="label", colour_by="label")


###Heatmap
chosen <- marker.info[[1]]
ordered <- chosen[order(chosen$rank.logFC.cohen),]
top.ranked <- ordered[ordered$rank.logFC.cohen <= 5,]
rownames(top.ranked)
plotGroupedHeatmap(sce.corrected, features=rownames(top.ranked), group="label", 
                   center=TRUE, zlim=c(-3, 3))


###FULL EFFECT###
marker.info <- scoreMarkers(sce.corrected, colLabels(sce.corrected), full.stats=TRUE)
chosen <- marker.info[[1]]
chosen$full.AUC
lyz.high <- c("PND007", "PND030", "PND060", "PND090", "PND150") # based on inspection of the previous Figure.
subset <- chosen$full.AUC[,colnames(chosen$full.AUC) %in% lyz.high]
to.show <- subset[computeMinRank(subset) <= 5,]
chosen.genes = rownames(to.show)
chosen.genes = c(chosen.genes, rownames(to.show))

for (x in c(2,3,4,5)){
  chosen <- marker.info[[x]]
  subset <- chosen$full.AUC[,colnames(chosen$full.AUC) %in% lyz.high]
  to.show <- subset[computeMinRank(subset) <= 5,]
  chosen.genes = c(chosen.genes, rownames(to.show))
}
chosen.genes = unique(chosen.genes)

#plotGroupedHeatmap(sce.corrected[,colLabels(sce.corrected) %in% lyz.high],features=rownames(to.show), group="label", center=TRUE, zlim=c(-3, 3))
fontsize = 20
scale=FALSE
plotGroupedHeatmap(sce.corrected[,colLabels(sce.corrected) %in% lyz.high],
                   features=chosen.genes, group="label", center=TRUE, scale=scale,
                   zlim=c(-3, 3), fontsize=fontsize, fontsize_row=fontsize/2)
