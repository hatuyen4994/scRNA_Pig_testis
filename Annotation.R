setwd("D:/R/scRNA_Pig_testis/")
sample = "PND7"
raw.path = file.path(sample)

markers.2 = list(
  "Spermatogonia_names"=c("UCHL1","CTCFL"), 
  "Spermatocytes_names"=c("SYCP1", "PIWIL1"), 
  "Round.spermatids_names"=c("ACRV1"), 
  "Elongating.STids_names"=c("SPACA1"), 
  "Immature.Sperm_names"=c("TNP1","PRM2"),
  "Sertoli_names"=c("CLU","AMH"), 
  "Leydig_names"=c("INSL3","STAR"), 
  "Myoid_names"=c("MYH11","ACTA2"), 
  "Endothelial_names"=c("VWF","PECAM1"), 
  "Macrophage_names"=c("CD68", "C1QA"),
  "T_names"=c("CD4","IL7R")
)
library(GSEABase)
data <- readxl::read_excel("markers.xlsx")
col.idx = grepl("_names", colnames(data))
markers.150 = data[1:150, col.idx]


z = as.list(markers.2)
cell.names = names(z)
cell.names  = lapply(cell.names, function(x) {
  paste(unlist(strsplit(x,"_"))[1],"_cell",sep="")
})
names(z) = unlist(cell.names)
all.sets <- lapply(names(z), function(x) {
  GeneSet(z[[x]], setName=x)        
})
all.sets <- GeneSetCollection(all.sets)


library(SingleCellExperiment)
library(DropletUtils)
genes.idx = read.table("common.genes.tsv")
genes.idx= genes.idx$V1
sce = read10xCounts(raw.path)
rownames(sce)= rowData(sce)$ID
sce = sce[genes.idx,]
gc()

rownames(sce)= rowData(sce)$Symbol
library(AUCell)
rankings <- AUCell_buildRankings(counts(sce),
                                 plotStats=FALSE, verbose=FALSE, splitByBlocks=TRUE)
cell.aucs <- AUCell_calcAUC(all.sets, rankings)
results <- t(assay(cell.aucs))
head(results)
new.labels <- colnames(results)[max.col(results)]
write.table(new.labels, file='new.labels.tsv', quote=FALSE, sep='\t', col.names = F, row.names=F)

cell.type = "Sertoli"
idx = new.labels==paste(cell.type,"_cell",sep="")
out.path = file.path(cell.type, raw.path, paste(sample,cell.type,".rds",sep=""))
saveRDS(sce[,idx],out.path)

cell.type = "Leydig"
idx = new.labels==paste(cell.type,"_cell",sep="")
out.path = file.path(cell.type, raw.path, paste(sample,cell.type,".rds",sep=""))
saveRDS(sce[,idx],out.path)