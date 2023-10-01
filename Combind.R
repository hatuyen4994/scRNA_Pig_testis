sce.7     =  readRDS("Sertoli/PND7/PND7Sertoli.rds")
sce.30   =   readRDS("Sertoli/PND30/PND30Sertoli.rds")
sce.60   =   readRDS("Sertoli/PND60/PND60Sertoli.rds")
sce.90   =   readRDS("Sertoli/PND90/PND90Sertoli.rds")
sce.150  =   readRDS("Sertoli/PND150/PND150Sertoli.rds")

sce_cbind = cbind(sce.7,sce.30,sce.60,sce.90)
rownames(sce_cbind) = rowData(sce_cbind)$ID
rownames(sce.150) = rowData(sce.150)$ID
common = intersect(rownames(sce_cbind), rownames(sce.150))
test = sce_cbind[common]
test2 = sce.150[common]
rowData(test) = NULL
rowData(test2) = NULL
sce =  cbind(test,test2)
rowData(sce) = rowData(sce_cbind[common])
colnames(sce) = sce$Barcode
#saveRDS(sce,"Sertoli7to150.rds")
saveRDS(sce,"Sertoli7to150_2markers.rds")
