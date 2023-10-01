library(SingleCellExperiment)
library(scran)
filter_h5ad = "PND7_filtered2.h5ad"
output = "output.h5ad"
sce = zellkonverter::readH5AD(filter_h5ad)
mean.unorm = mean(assay(sce))
clusters = quickCluster(sce, assay.type="X")
sce = computeSumFactors(sce,cluster=clusters, assay.type="X")
sce = logNormCounts(sce, assay.type="X")
plot(librarySizeFactors(sce, assay.type="X"), sizeFactors(sce)
     ,pch=16,xlab="Library", ylab="Decon",log="xy")

rnames = rownames(sce)
normalized_mat = assay(sce,"logcounts")
gc()
sce = zellkonverter::readH5AD(filter_h5ad)
gc()
assay(sce) = normalized_mat
rownames(sce) = rnames
mean.norm = mean(assay(sce))
sce = 
mean.unorm
mean.norm
zellkonverter::writeH5AD(sce, output)
