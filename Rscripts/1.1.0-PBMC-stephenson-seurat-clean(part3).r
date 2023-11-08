# original paper from Stephenson el al
# https://www.nature.com/articles/s41591-021-01329-2

# Data available from here 
# https://www.covid19cellatlas.org/index.patient.html
# under "COVID-19 PBMC Ncl-Cambridge-UCL"

# This notebook takes
# source.file: the source R file of the single cell experiement 
# res.dir: the result directory to save the cleaned version of the object

# This notebook loads the single cell object, kept only the top major celltypes: CD4,8, NK, Mono, B, under col "decon.L1"
# it also pushes these major celltype cells throught the pipeline again just as a sanity check
# new low dims are saved under "decon.pca" and "decon.umap"
library(Seurat)
library(ggplot2)
library(ggpubr)
set.seed(2023)

source.file = "../Data//RNA//Simulation-PBMC/stephenson.seu"
res.dir = "../Data/RNA/Simulation-PBMC"

stephenson = readRDS(source.file)
DefaultAssay(stephenson) <- "raw"

stephenson

stephenson@meta.data

colnames(stephenson@meta.data)

options(repr.plot.width = 12, repr.plot.height = 6, repr.plot.res = 250)
p1 = DimPlot(stephenson, reduction = "umap", group.by = "initial_clustering",    label = T, label.size = 3, repel = T) + NoLegend()
p2 = DimPlot(stephenson, reduction = "umap", group.by = "full_clustering",       label = T, label.size = 3, repel = T) + NoLegend()
p1 + p2 

table(stephenson@meta.data$full_clustering)

#directly group by full_clustering
print(sort(unique(stephenson@meta.data$full_clustering)))

NK.filter = stephenson@meta.data[,"full_clustering"] %in% 
            c("NKT", "NK_16hi", "NK_56hi", "NK_prolif")

CD4.filter = stephenson@meta.data[,"full_clustering"] %in% 
             c("CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Naive", "CD4.Prolif", "CD4.Tfh", "CD4.Th1", "CD4.Th17", "CD4.Th2")

CD8.filter = stephenson@meta.data[,"full_clustering"] %in% 
             c("CD8.EM", "CD8.Naive", "CD8.Prolif", "CD8.TE")

Mono.CD14.filter = stephenson@meta.data[,"full_clustering"] %in% 
              c("CD14_mono", "CD83_CD14_mono", "Mono_prolif")

Mono.CD16.filter = stephenson@meta.data[,"full_clustering"] %in% 
              c("C1_CD16_mono", "CD16_mono")


B.filter    = stephenson@meta.data[,"full_clustering"] %in% 
              c("B_exhausted", "B_immature", "B_malignant", "B_naive", "B_non-switched_memory", "B_switched_memory")

Plasma.filter = stephenson@meta.data[,"full_clustering"] %in% 
                c("Plasma_cell_IgA", "Plasma_cell_IgG", "Plasma_cell_IgM", "Plasmablast")

# L2
stephenson@meta.data[NK.filter,        "decon.L2"] = "NK"
stephenson@meta.data[CD4.filter,       "decon.L2"] = "CD4"
stephenson@meta.data[CD8.filter,       "decon.L2"] = "CD8"
stephenson@meta.data[Mono.CD14.filter, "decon.L2"] = "Mono.CD14"
stephenson@meta.data[Mono.CD16.filter, "decon.L2"] = "Mono.CD16"
stephenson@meta.data[B.filter,         "decon.L2"] = "B"
stephenson@meta.data[Plasma.filter,    "decon.L2"] = "Plasma"

# L1
stephenson@meta.data[NK.filter,   "decon.L1"] = "NK"
stephenson@meta.data[CD4.filter,  "decon.L1"] = "CD4"
stephenson@meta.data[CD8.filter,  "decon.L1"] = "CD8"
stephenson@meta.data[Mono.CD14.filter | Mono.CD16.filter, "decon.L1"] = "Mono"
stephenson@meta.data[B.filter | Plasma.filter,    "decon.L1"] = "B"

#check is missing any thing major
keep.cells = NK.filter | CD4.filter| CD8.filter| Mono.CD14.filter | Mono.CD16.filter | B.filter | Plasma.filter 
filter.cells.table = table(stephenson@meta.data[!keep.cells, "full_clustering"])
filter.cells.table[filter.cells.table>0]

stephenson = stephenson[, keep.cells]
stephenson

options(repr.plot.width = 18, repr.plot.height = 6, repr.plot.res = 250)
p1 = DimPlot(stephenson, reduction = "umap", group.by = "decon.L1",  label = T, label.size = 3, repel = T) + NoLegend()
p2 = DimPlot(stephenson, reduction = "umap", group.by = "decon.L2",  label = T, label.size = 3, repel = T) + NoLegend()
p3 = DimPlot(stephenson, reduction = "umap", group.by = "full_clustering", label = T, label.size = 3, repel = T) + NoLegend()
p1 + p2 + p3

stephenson@meta.data["nCount_raw"]   = as.numeric(unlist(stephenson@meta.data["nCount_raw"]))
stephenson@meta.data["nFeature_raw"] = as.numeric(unlist(stephenson@meta.data["nFeature_raw"]))

stephenson <- PercentageFeatureSet(stephenson, "^MT-",      col.name = "percent_mito", assay = "raw")
stephenson <- PercentageFeatureSet(stephenson, "^RP[SL]",   col.name = "percent_ribo", assay = "raw")
stephenson <- PercentageFeatureSet(stephenson, "^HB[^(P)]", col.name = "percent_hb",   assay = "raw")

feats <- c("nFeature_raw", "nCount_raw", "percent_mito", "percent_ribo", "percent_hb")
options(repr.plot.width = 20, repr.plot.height = 5, repr.plot.res = 250)
VlnPlot(stephenson, features = feats, pt.size = 0.01, ncol = 5) + NoLegend()

quantile(stephenson@meta.data$nFeature_raw, c(0.05, 0.95))

quantile(stephenson@meta.data$nCount_raw, c(0.05, 0.95))

keep_cells <- (stephenson@meta.data$nFeature_raw > 500) & (stephenson@meta.data$nFeature_raw < 2500) &
              (stephenson@meta.data$nCount_raw   > 2000) & (stephenson@meta.data$nCount_raw < 15000) &
              (stephenson@meta.data$percent_mito < 5) &
              (stephenson@meta.data$percent_ribo > 1) & 
              (stephenson@meta.data$percent_hb   < 1)

selected_f <- rownames(stephenson)[Matrix::rowSums(stephenson@assays[["raw"]]@counts>=1) >= 5]
stephenson <- subset(stephenson, features = selected_f, cells = colnames(stephenson)[keep_cells])
stephenson

feats <- c("nCount_raw", "nFeature_raw", "percent_mito", "percent_ribo", "percent_hb")
options(repr.plot.width = 20, repr.plot.height = 5, repr.plot.res = 250)
VlnPlot(stephenson, features = feats, pt.size = 0.01, ncol = 5) + NoLegend()

stephenson <- NormalizeData(stephenson, normalization.method = "LogNormalize", scale.factor = 10000)
stephenson <- FindVariableFeatures(stephenson, selection.method = "vst", nfeatures = 2000)

options(repr.plot.width = 15, repr.plot.height = 4, repr.plot.res = 250)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(stephenson), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(stephenson)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

stephenson <- ScaleData(stephenson)
stephenson <- RunPCA(stephenson, reduction.name = "decon_pca",)

#seu.obj: a seurat object with meaningful meta.data field
#embedding: a low dim matrix of cells by low dime loading
#reduction.name: usually choose from "pca", "cca" for x and y axis name only 
#from: starting index of the low dimension
#to: ending index of the low dimension
#group.by: the column in the meta.data that we use to color the dots in the scatter plot
low_dim_plots = function(seu.obj, embedding, reduction.name, from, to, group.by = NULL){
    plts <- vector(mode = "list", length = 6)
    counter = 1
    for (i in seq(from = from, to = to, by = 2)){
        g = ggplot(data.frame(x = embedding[,i], 
                              y = embedding[,(i+1)],
                              color = if(is.null(group.by)) 1 else as.numeric(seu.obj@meta.data[,group.by])), 
               aes(x = x, y = y, color = color)) +
        geom_point(alpha = 0.5, size = 0.1) +  
        xlab(paste0(reduction.name, " ", i)) + 
        ylab(paste0(reduction.name, " ", (i+1))) 

        plts[[counter]] = g
        counter = counter + 1
    }
    g <- ggarrange(plotlist = plts, ncol = 4, nrow = floor(counter/4) + 1)
    return(g)
}

options(repr.plot.width = 12, repr.plot.height = 10, repr.plot.res = 150)
g = low_dim_plots(seu.obj = stephenson, 
                  embedding = stephenson@ reductions[["decon_pca"]]@ cell.embeddings, 
                  reduction.name = "decon_pca", from = 1, to = 30, group.by = "nCount_raw")
g

# using reduction = "decon_pca" is not good. some strong batch effect. so here just using the pca provided 
stephenson <- FindNeighbors(stephenson,  dims = 1:20)
stephenson <- RunUMAP(stephenson, dims = 1:20, reduction.name = "decon_umap")

options(repr.plot.width = 12, repr.plot.height = 4, repr.plot.res = 250)
p1 = DimPlot(stephenson, reduction = "decon_umap", group.by = "full_clustering", label = T, label.size = 3, repel = T) + NoLegend()
p2 = DimPlot(stephenson, reduction = "decon_umap", group.by = "decon.L1",        label = T, label.size = 3, repel = T) + NoLegend()
p3 = DimPlot(stephenson, reduction = "decon_umap", group.by = "decon.L2",        label = T, label.size = 3, repel = T) + NoLegend()
p1 + p2 + p3


stephenson

if (!dir.exists(res.dir)){
  dir.create(res.dir, recursive = T)
}
saveRDS(stephenson, file.path(res.dir, "stephenson.clean.seu"))
