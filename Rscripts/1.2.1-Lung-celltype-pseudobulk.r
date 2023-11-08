# This notebook takes 
# res.dir: the path to the directory that have the seu.rds of the pbmc single cell data
# and construct the raw celltype level psudobulk material needed for simluation 

# if partial.W is T
# it filter out samples with low cell counts or imbalance cell type prop
# and construct the raw celltype level psudobulk material needed for simluation 
# W: remove samples with lower than 1000 cells (around 10 percent)
#    remove samples with any celltype less than 1%

# if partial.W is F
# W: use all samples' W 
# Z: calculated as the AverageExpression on the CPM. (average is taken in the unloged space)
# params: $mus, $sigmas, $corrs: calculated based on Z with samples outside max_sds removed 
# entropy: calculated based on params$corrs

#save the final object with W, Z, params, entorpy as list element as pbmc.pseudobulk.rds in the res.dir

library("Seurat")
library("pracma")
library("ggplot2")
library("ggpubr")

source("analysis.utils.r")
source("simulate.expression.utils.r")
set.seed(2023)

partial.W = F

source.col  = "decon.L1"
#source.col  = "decon.L2"
sample.col = "sample"

res.dir = "../Data/RNA/Simulation-Lung"

seu.obj  = readRDS(file.path(res.dir, "HLCA.clean.seu"))

DefaultAssay(seu.obj) <- "RNA"
seu.obj <- NormalizeData(seu.obj, normalization.method = "LogNormalize", scale.factor = 1e6)

seu.obj

options(repr.plot.width = 15, repr.plot.height = 5, repr.plot.res = 250)
p1 = DimPlot(seu.obj, reduction = "decon_umap", group.by = "ann_level_2", label = T) + NoLegend()
p2 = DimPlot(seu.obj, reduction = "decon_umap", group.by = "decon.L1", label = T) + NoLegend()
p3 = DimPlot(seu.obj, reduction = "decon_umap", group.by = "decon.L2", label = T) + NoLegend()
p1 + p2 + p3

colnames(seu.obj@meta.data)

sample.ids = unique(seu.obj@meta.data[,sample.col])
source.ids = unique(seu.obj@meta.data[,source.col])

print(paste0("total number of sample: ", length(sample.ids)))
print(paste0("total number of source: ", length(source.ids)))

W.counts = matrix(0, length(sample.ids), length(source.ids))
rownames(W.counts) = sample.ids
colnames(W.counts) = source.ids

for (sample.id in sample.ids){
    for (source.id in source.ids){
        W.counts[sample.id,source.id] = sum((seu.obj@meta.data[, sample.col] == sample.id) & (seu.obj@meta.data[,source.col] == source.id))
    }
}

W.prop = W.counts/repmat(as.matrix(rowSums(W.counts)), 1, ncol(W.counts))

hist(log10(rowSums(W.counts)))

sum(rowSums(W.counts) > 1000)

sum(rowSums(W.prop < 0.01) == 0 & rowSums(W.counts) > 1000)

if (partial.W){
    #samples with at least 1000 cells (will remove 10% of the sample)
    keep.samples.counts = rowSums(W.counts) > 1000
    #samples with at least 1% cell proportion in all celltypes 
    keep.samples.prop = rowSums(W.prop < 0.01) == 0 

    #filter W.prop
    keep.samples = keep.samples.counts & keep.samples.prop
    W.counts = W.counts[keep.samples, ]
    W.prop   = W.prop[keep.samples, ]

    #update sample.ids, source.ids
    W.prop = W.prop[, order(-colMeans(W.prop))]
    sample.ids  = rownames(W.prop)
    source.ids  = colnames(W.prop)
}
#set W to filtered W.prop
W = W.prop

W.counts

W

avg_counts = AverageExpression(object = seu.obj, 
                        group.by  = c(source.col,sample.col), 
                        assays = "RNA",
                        slot = "data") # if slot is "data", the average is computed in the original "un-log" space

feature.ids = rownames(seu.obj)
Z = array(0,c(length(source.ids), length(feature.ids), length(sample.ids)))
dimnames(Z)[[1]] = source.ids
dimnames(Z)[[2]] = feature.ids
dimnames(Z)[[3]] = sample.ids

for (source.id in source.ids){
    for (sample.id in sample.ids){
        Z[source.id,,sample.id] = tryCatch(
            {
                avg_counts$RNA[, paste(source.id, sample.id, sep = "_")]
            },
            error=function(cond) {
                print(paste("no cells for", source.id, sample.id))
                0
            })
    }
}

Z[1, 1:10, 1:10]

# #manually do that 
# gen.pseudobulk = function(seu.meta, seu.counts, sample.group.by, source.group.by){
    
#     source.ids  = as.character(unique(seu.meta[[source.group.by]]))
#     sample.ids  = as.character(unique(seu.meta[[sample.group.by]]))
#     feature.ids = as.character(rownames(seu.counts))

#     k = length(source.ids)
#     m = length(feature.ids)
#     n = length(sample.ids)

#     W = matrix(0, n, k)
#     Z = array(0,c(k, m, n))
#     for (i in 1:n){
#         sample.id = sample.ids[i]        
#         for (h in 1:k){
#             source.id  = source.ids[h]
#             cells.use = rownames(seu.meta[(seu.meta[sample.group.by] == sample.id & seu.meta[source.group.by] == source.id),])
            
#             # if more than one cell, take the average 
#             if(length(cells.use) >= 2){
#                 Z[h,,i] = as.matrix(rowMeans(seu.counts[feature.ids,cells.use]))
#                 W[i,h]  = length(cells.use)
#             # if only one cell, use it 
#             }else if (length(cells.use) == 1) {
#                 Z[h,,i] = as.matrix(seu.counts[feature.ids,cells.use])
#                 W[i,h]  = 1
#             # else just put all 0s
#             }else{
#                 Z[h,,i] = 0
#                 W[i,h]  = 0
#             }
#         }# end of source
#     }# end of sample
    
#     # process W and reorder things
#     W = W/repmat(as.matrix(rowSums(W)), 1, dim(W)[2])
    
#     reorder.index = order(-colMeans(W))
#     source.ids = source.ids[reorder.index]
#     W = W[,reorder.index]
#     rownames(W) = sample.ids
#     colnames(W) = source.ids
    
#     Z = Z[reorder.index,,]
#     dimnames(Z)[[1]] = source.ids
#     dimnames(Z)[[2]] = feature.ids
#     dimnames(Z)[[3]] = sample.ids
    
#     # generate X
#     X = matrix(0, m, n)
#     rownames(X) = feature.ids
#     colnames(X) = sample.ids
#     for (h in 1:k){
#         X = X + Z[h,,] * repmat(t(W[,h]), m ,1)
#     }
    
#     res = list()
#     res$Z = Z
#     res$W = W
#     res$X = X
#     res$source.ids  = source.ids
#     res$feature.ids = feature.ids
#     res$sample.ids  = sample.ids
#     return(res)
# }

# seu.obj  = NormalizeData(object = seu.obj, 
#                             normalization.method = "RC",
#                             scale.factor = 1e6)
# pseudobulk.cpm = gen.pseudobulk(seu.meta    = seu.obj@meta.data,
#                                 seu.counts  = seu.obj@assays$raw@data,
#                                 sample.group.by = "sample_id",
#                                 source.group.by = "decon.L1")

# pseudobulk.cpm$Z["CD4", 1:10, 1:10]

# pseudobulk.cpm$X[1:10, 1:10]

# DefaultAssay(seu.obj) <- "raw"
# seu.obj <- NormalizeData(seu.obj, normalization.method = "LogNormalize", scale.factor = 1e6)

# # Generate X

# X = matrix(0, length(feature.ids), length(sample.ids))
# rownames(X) = feature.ids
# colnames(X) = sample.ids
# for (h in 1:length(source.ids)){
#     X = X + Z[h,,] * repmat(t(W[,h]), length(feature.ids), 1)
# }

# X[1:10, 1:10] 

params    = calc_params_from_Z(Z, max_sds = 2)
entropy   = calc_entropy(params$corrs)

hist(entropy, main = "All genes entropy")

seu.obj = FindVariableFeatures(seu.obj, nfeatures = 10000)

hist(entropy[VariableFeatures(seu.obj), ], main = "Top 10k HVF entropy")

HLCA.pseudobulk = list(Z = Z, W = W, 
                       params = params, entropy = entropy,
                       HVF = VariableFeatures(seu.obj),
                       HEF = highly_expressed_features(Z = Z, expression.qtl = 0)[1:10000])

saveRDS(HLCA.pseudobulk,  file.path(res.dir, paste0("HLCA.pseudobulk.", 
                                                    source.col, ".",
                                                    if (partial.W) "partial.W." else "all.W.",
                                                    "rds")))

str(HLCA.pseudobulk)

low.entropy.thr  = quantile(entropy, 0.1)
high.entropy.thr = quantile(entropy, 0.9)

low_entropy_gene  = rownames(HLCA.pseudobulk$entropy)[sample(which (HLCA.pseudobulk$entropy < low.entropy.thr), 1)]
high_entropy_gene = rownames(HLCA.pseudobulk$entropy)[sample(which (HLCA.pseudobulk$entropy > high.entropy.thr), 1)]

HLCA.pseudobulk$params$mus[low_entropy_gene,]

HLCA.pseudobulk$params$sigmas[low_entropy_gene,,]

HLCA.pseudobulk$params$corrs[low_entropy_gene,,]

HLCA.pseudobulk$params$mus[high_entropy_gene,]

HLCA.pseudobulk$params$sigmas[high_entropy_gene,,]

HLCA.pseudobulk$params$corrs[high_entropy_gene,,]

plot_example_genes = function(pseudobulk, gene.names){
    plots = list()
     
    counter = 1
    for (gene in gene.names){
        cor_mat = pseudobulk$params$corrs[gene,,]
        
        #most correlated 2 sources
        idx = which(cor_mat == max(cor_mat[lower.tri(cor_mat,diag = F)]), arr.ind = TRUE)[1,]
        title = paste(gene, "source", idx[1], "vs", idx[2])
        
        plot_df = data.frame(t(pseudobulk$Z[idx,gene,]))
        colnames(plot_df) = c("source.1" , "source.2")
        
        
        plots[[counter]] = ggplot(plot_df, aes(x = log10(1 + source.1), y = log10(1 + source.2))) + 
                           geom_point(alpha = 0.5, size = 1) + 
                           ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio = 1)
        
        plots[[counter + 1]] = ggplot(plot_df, aes(x = source.1, y = source.2)) + 
                               geom_point(alpha = 0.5, size = 1) + 
                               ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio = 1)
        counter = counter + 2
        
    }
    
    g <- egg::ggarrange(plots = plots, ncol = 2)  
    return(g)
}

g = plot_example_genes(HLCA.pseudobulk, c(low_entropy_gene, high_entropy_gene))


