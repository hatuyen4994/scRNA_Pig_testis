#Things to remember
"' SCE data have slots:   assay (primary data),
                          colData(cell data), 
                          rowData and rowranges (feature data),
                          metada,
                          reducedDims (like cell data)
                          altExps (SE in SCE)

'"

#Create DF about location
ID_MT = read.table("ID_Chromosome.tsv")
DF = DataFrame(
  chromosome=ID_MT$V3,
  row.names=ID_MT$V1
)
rowData(sce)$location = DF[rowData(sce)$ID,]

#Add element to a vector
vector_new = c(vector_old, new_element)

#Download something
mm10.gtf <- bfcrpath(bfc, file.path("http://ftp.ensembl.org/pub/release-82",
                                    "gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz"))

#add assay
sce = SingleCellExperiment(assay = list(counts=matrix1))
assay(sce, "new_counts") = matrix2

#choose assay to keep
assays(sce) = assays(sce)[1]

#add coldata
colData(sce) = coldata #en masse
colData$phenotype = coldata$phenotype #1 by 1
#check colname of assay slot is the same as row name in colData 
stopifnot(identical(rownames(coldata), colnames(mat)))


#combines 2 sce by columns assuming all objects have the same column annotation 
#values and compatible row annotation fields.
sce_cbind = cbind(sce1,sce2)

#combines 2 sce by row assuming all objects have the same column annotation 
#values and compatible row annotation fields.
sce_rbind = rbind(sce1, sce2)

##
rowRanges(sce) <- gene.data[rownames(sce)]


#scran::computerSumFactors is similar to scater::librarySizeFactors
#Both put a sizeFactor colData in the SCE


#Labeling cells 
"'colLabels(sce) = factors with level (containing categorical values)'"


#Markers for annotation
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