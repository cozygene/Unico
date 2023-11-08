library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("../Data/RNA/Simulation-PBMC/stephenson.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc <- LoadH5Seurat("../Data/RNA/Simulation-PBMC/stephenson.h5seurat")
saveRDS(pbmc, "../Data/RNA/Simulation-PBMC/stephenson.seu")
