# original paper from Sikkema el al
# https://www.biorxiv.org/content/10.1101/2022.03.10.483747v1

# Data available from here 
# https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
# under "An integrated cell atlas of the human lung in health and disease (core)"

# This notebook takes
# source.file: the source rds file of the single cell experiement 
# media1.file: the first supplementary file in the paper, exported from xlsx to txt format
# media2.file: the second supplementary file in the paper, exported from xlsx to txt format
# res.dir: the result directory to save the cleaned version of the object

# This notebook loads the single cell object, 
# kept only the top major celltypes: 
# under col "decon.L1"  "Endothelial" "Stroma" "Myeloid" "Lymphoid" "Airway Epithelium" "Alveolar Epithelium"
# under col "decon.L2": "Endothelial" "Stroma" "Immune" "Epithelial"

# keep only sample from location: parenchyma

# it also pushes these filtered cells throught the pipeline again just as a sanity check
# new low dims are saved under "decon.pca" and "decon.umap"
library(Seurat)
library(data.table)

library(ggplot2)
library(ggpubr)
set.seed(2023)

source.file = "../Data/RNA/Simulation-Lung/HLCA.rds"
media1.file = "../Data/RNA/Simulation-Lung/media-1.txt"
media2.file = "../Data/RNA/Simulation-Lung/media-2.txt"

res.dir = "../Data/RNA/Simulation-Lung"

HLCA   = readRDS(source.file)
table1 = as.matrix(data.frame(fread(media1.file)))
table2 = as.matrix(data.frame(fread(media2.file, fill=TRUE ), row.names=1))

HLCA

HLCA@meta.data

colnames(HLCA@meta.data)

options(repr.plot.width = 15, repr.plot.height = 5, repr.plot.res = 250)
p1 = DimPlot(HLCA, group.by = "ann_level_1", label = T) + NoLegend()
p2 = DimPlot(HLCA, group.by = "ann_level_2", label = T) + NoLegend()
p3 = DimPlot(HLCA, group.by = "cell_type", label = T, label.size = 2, repel = T) + NoLegend()
p1 + p2 + p3

Endothelial.filter = HLCA@meta.data[,"ann_level_1"] %in% c("Endothelial")
Stroma.filter      = HLCA@meta.data[,"ann_level_1"] %in% c("Stroma")
Myeloid.filter              = HLCA@meta.data[,"ann_level_2"] %in% c("Myeloid")
Lymphoid.filter             = HLCA@meta.data[,"ann_level_2"] %in% c("Lymphoid")
AirwayEpithelium.filter     = HLCA@meta.data[,"ann_level_2"] %in% c("Airway epithelium")
AlveolarEpithelium.filter   = HLCA@meta.data[,"ann_level_2"] %in% c("Alveolar epithelium")

HLCA@meta.data[Endothelial.filter,        "decon.L1"] = "Endothelial"
HLCA@meta.data[Stroma.filter,             "decon.L1"] = "Stroma"
HLCA@meta.data[Myeloid.filter,            "decon.L1"] = "Immune"
HLCA@meta.data[Lymphoid.filter,           "decon.L1"] = "Immune"
HLCA@meta.data[AirwayEpithelium.filter,   "decon.L1"] = "Epithelial"
HLCA@meta.data[AlveolarEpithelium.filter, "decon.L1"] = "Epithelial"

HLCA@meta.data[Endothelial.filter,        "decon.L2"] = "Endothelial"
HLCA@meta.data[Stroma.filter,             "decon.L2"] = "Stroma"
HLCA@meta.data[Myeloid.filter,            "decon.L2"] = "Myeloid"
HLCA@meta.data[Lymphoid.filter,           "decon.L2"] = "Lymphoid"
HLCA@meta.data[AirwayEpithelium.filter,   "decon.L2"] = "Airway Epithelium"
HLCA@meta.data[AlveolarEpithelium.filter, "decon.L2"] = "Alveolar Epithelium"

#check is missing any thing major
keep.cells = Endothelial.filter | Stroma.filter| Myeloid.filter| Lymphoid.filter | AirwayEpithelium.filter | AlveolarEpithelium.filter
filter.cells.table = table(HLCA@meta.data[!keep.cells, "ann_level_2"])
filter.cells.table[filter.cells.table>0]

HLCA = HLCA[, keep.cells]

options(repr.plot.width = 15, repr.plot.height = 5, repr.plot.res = 250)
p1 = DimPlot(HLCA, group.by = "ann_level_2", label = T) + NoLegend()
p2 = DimPlot(HLCA, group.by = "decon.L1", label = T) + NoLegend()
p3 = DimPlot(HLCA, group.by = "decon.L2", label = T) + NoLegend()
p1 + p2 + p3

# the seurat object only contains "core", which comes from 11 studies, a total of 14 datasets. 
# some study contains multiple datasets
# 2 studies are labeled as "not primary" dataset. this concept only shows up in the meta data part of the seurat object. not mentioned in the paper
# for core: a total of 107 unique individuals, but over 166 unique samples 

print("in the full samples file media-2.txt")
print(paste0("number of samples: ",        length(unique(rownames(table2)))))
print(paste0("number of core samples: ",   length(unique(rownames(table2[table2[,"HLCA_core_or_extension"] == "core",])))))
print(paste0("number of core individual: ",length(unique(table2[table2[,"HLCA_core_or_extension"] == "core","subject_ID"]))))


print("seu object oly contains cells from core part of the study. so 11 unique study, 14 datasets, 107 unique individuals, 166 samples")
print("in the actual seu object")
print(paste0("number of samples: ",     length(unique((HLCA@meta.data$sample)))))
print(paste0("number of individuals: ", length(unique((HLCA@meta.data$donor_id)))))

print("Primary study: ")
print(unique((HLCA@meta.data[HLCA@meta.data$is_primary_data, "study"])))

print("Not primary study: ")
print(unique((HLCA@meta.data[!HLCA@meta.data$is_primary_data, "study"])))

strict = F # if to exclude unhealthy
if(strict){
    # 82 samples
    keep.sample.ids = rownames(table2[(table2[,"HLCA_core_or_extension"] == "core") & 
                                      (table2[,"lung_condition"] == "healthy") &      
                                      (table2[,"anatomical_region_level_1"] =="parenchyma"), ])

}else{
    # 90 samples
    keep.sample.ids = rownames(table2[(table2[,"HLCA_core_or_extension"] == "core") & 
                                      (table2[,"anatomical_region_level_1"] =="parenchyma"), ])
}
keep.sample.ids = as.character(unique(keep.sample.ids))
print(paste0("keeping: ", length(keep.sample.ids), " samples"))


HLCA = HLCA[, (HLCA@meta.data$"sample" %in% keep.sample.ids)]
print(paste0("kept porportion of cells: ", dim(HLCA)[2]/dim(HLCA)[2]))
print(paste0("raw counts dimension: ", dim(HLCA@assays$RNA@counts)[1], " by ", dim(HLCA@assays$RNA@counts)[2]))

HLCA

# ################## redo scvi way ##################
# scanvi_check = F
# if (scanvi_check){
# 	HLCA.sub <- FindNeighbors(HLCA.sub, reduction = "scanvi_emb", dims = 1:10)
# 	HLCA.sub <- FindClusters(HLCA.sub,  reduction = "scanvi_emb", resolution = 0.1)
# 	HLCA.sub <- RunUMAP(HLCA.sub,  reduction = "scanvi_emb", dims = 1:10, reduction.name = "umap_scvi", reduction.key = "UMAPSCVI_")
	
# 	Idents(HLCA.sub) <- 'ann_level_1'
# 	p1 = DimPlot(HLCA.sub, reduction = "umap_scvi")
	
# 	Idents(HLCA.sub) <- 'ann_level_2'
# 	p2 = DimPlot(HLCA.sub, reduction = "umap_scvi")

# 	p1 + p2
# }

HLCA <- PercentageFeatureSet(HLCA, "^MT-",      col.name = "percent_mito", assay = "RNA")
HLCA <- PercentageFeatureSet(HLCA, "^RP[SL]",   col.name = "percent_ribo", assay = "RNA")
HLCA <- PercentageFeatureSet(HLCA, "^HB[^(P)]", col.name = "percent_hb",   assay = "RNA")

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
options(repr.plot.width = 20, repr.plot.height = 5, repr.plot.res = 250)
VlnPlot(HLCA, features = feats, pt.size = 0.01, ncol = 5) + NoLegend()

quantile(HLCA@meta.data$nFeature_RNA, c(0.05, 0.95))

quantile(HLCA@meta.data$nCount_RNA, c(0.05, 0.95))

keep_cells <- (HLCA@meta.data$nFeature_RNA > 500) & (HLCA@meta.data$nFeature_RNA < 4000) &
              (HLCA@meta.data$nCount_RNA > 1000) & (HLCA@meta.data$nCount_RNA < 20000) 
             
selected_f <- rownames(HLCA)[Matrix::rowSums(HLCA@assays[["RNA"]]@counts>=1) >= 5]
HLCA <- subset(HLCA, features = selected_f, cells = colnames(HLCA)[keep_cells])
HLCA

options(repr.plot.width = 20, repr.plot.height = 5, repr.plot.res = 250)
VlnPlot(HLCA, features = feats, pt.size = 0.01, ncol = 5) + NoLegend()



HLCA <- NormalizeData(HLCA, normalization.method = "LogNormalize", scale.factor = 10000)
HLCA <- FindVariableFeatures(HLCA, selection.method = "vst", nfeatures = 2000)
HLCA <- ScaleData(HLCA) 
HLCA <- RunPCA(HLCA, features = VariableFeatures(object = HLCA), reduction.name = "decon_pca")

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
g = low_dim_plots(seu.obj = HLCA, 
                  embedding = HLCA@ reductions[["decon_pca"]]@ cell.embeddings, 
                  reduction.name = "decon_pca", from = 1, to = 30, group.by = "nCount_RNA")
g

HLCA <- FindNeighbors(HLCA, dims = 1:10, reduction = "decon_pca")
HLCA <- RunUMAP(HLCA, reduction = "decon_pca", dims = 1:10, reduction.name = "decon_umap")

options(repr.plot.width = 15, repr.plot.height = 5, repr.plot.res = 250)
p1 = DimPlot(HLCA, reduction = "decon_umap", group.by = "ann_level_2", label = T) + NoLegend()
p2 = DimPlot(HLCA, reduction = "decon_umap", group.by = "decon.L1", label = T) + NoLegend()
p3 = DimPlot(HLCA, reduction = "decon_umap", group.by = "decon.L2", label = T) + NoLegend()
p1 + p2 + p3

# #decide to delay this step to simulation data generation
# HLCA = FindVariableFeatures(HLCA, selection.method = "vst", nfeatures = 10000)
# HLCA = HLCA[VariableFeatures(HLCA), ]

HLCA

if (!dir.exists(res.dir)){
  dir.create(res.dir, recursive = T)
}
saveRDS(HLCA, file.path(res.dir, "HLCA.clean.seu"))
