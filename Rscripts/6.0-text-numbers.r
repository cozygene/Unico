library(Seurat)
library(xtable)
library(matrixStats)
set.seed(2023)

sim.pbmc.stats = readRDS("../Figure//Simulation-PBMC/decon.L1.max_stds.2/RNA_Simulation-PBMC_decon.L1_500_stats.list.rds")
sim.lung.stats = readRDS("../Figure//Simulation-Lung/decon.L1.max_stds.2/RNA_Simulation-Lung_decon.L1_500_stats.list.rds")

str(sim.pbmc.stats)

Unico.pbmc.covar = sim.pbmc.stats$covar.barplot.df[sim.pbmc.stats$covar.barplot.df$method == "Unico",]
bMIND.pbmc.covar = sim.pbmc.stats$covar.barplot.df[sim.pbmc.stats$covar.barplot.df$method == "bMIND",]

Unico.lung.covar = sim.lung.stats$covar.barplot.df[sim.lung.stats$covar.barplot.df$method == "Unico",]
bMIND.lung.covar = sim.lung.stats$covar.barplot.df[sim.lung.stats$covar.barplot.df$method == "bMIND",]

Unico.pbmc.covar

bMIND.pbmc.covar

Unico.lung.covar

bMIND.lung.covar

pbmc.improve = (mean(Unico.pbmc.covar$mean) - mean(bMIND.pbmc.covar$mean))/(mean(bMIND.pbmc.covar$mean))
lung.improve = (mean(Unico.lung.covar$mean) - mean(bMIND.lung.covar$mean))/(mean(bMIND.lung.covar$mean))

pbmc.improve
lung.improve

(pbmc.improve + lung.improve)/2



Unico.pbmc.Z.corr = sim.pbmc.stats$ Z.corrs.list$ Unico
TCA.pbmc.Z.corr   = sim.pbmc.stats$ Z.corrs.list$ TCA

Unico.lung.Z.corr = sim.lung.stats$ Z.corrs.list$ Unico
TCA.lung.Z.corr   = sim.lung.stats$ Z.corrs.list$ TCA

colMedians(Unico.pbmc.Z.corr, na.rm = T)

colMedians(TCA.pbmc.Z.corr, na.rm = T)

colMedians(Unico.lung.Z.corr, na.rm = T)

colMedians(TCA.lung.Z.corr, na.rm = T)

pbmc.improve = (mean(colMedians(Unico.pbmc.Z.corr, na.rm = T)) - mean(colMedians(TCA.pbmc.Z.corr, na.rm = T)))/mean(colMedians(TCA.pbmc.Z.corr, na.rm = T))

lung.improve = (mean(colMedians(Unico.lung.Z.corr, na.rm = T)) - mean(colMedians(TCA.lung.Z.corr, na.rm = T)))/mean(colMedians(TCA.lung.Z.corr, na.rm = T))

pbmc.improve

lung.improve

(pbmc.improve + lung.improve)/2

Unico.pbmc.low.corr = mean(c(sim.pbmc.stats$boxplot.meta$CD4 $median.mat["Low Entropy", "Unico"],
                             sim.pbmc.stats$boxplot.meta$NK  $median.mat["Low Entropy", "Unico"],
                             sim.pbmc.stats$boxplot.meta$CD8 $median.mat["Low Entropy", "Unico"],
                             sim.pbmc.stats$boxplot.meta$Mono$median.mat["Low Entropy", "Unico"],
                             sim.pbmc.stats$boxplot.meta$B   $median.mat["Low Entropy", "Unico"]))

TCA.pbmc.low.corr = mean(c(sim.pbmc.stats$boxplot.meta$CD4 $median.mat["Low Entropy", "TCA"],
                           sim.pbmc.stats$boxplot.meta$NK  $median.mat["Low Entropy", "TCA"],
                           sim.pbmc.stats$boxplot.meta$CD8 $median.mat["Low Entropy", "TCA"],
                           sim.pbmc.stats$boxplot.meta$Mono$median.mat["Low Entropy", "TCA"],
                           sim.pbmc.stats$boxplot.meta$B   $median.mat["Low Entropy", "TCA"]))


Unico.pbmc.high.corr = mean(c(sim.pbmc.stats$boxplot.meta$CD4 $median.mat["High Entropy", "Unico"],
                              sim.pbmc.stats$boxplot.meta$NK  $median.mat["High Entropy", "Unico"],
                              sim.pbmc.stats$boxplot.meta$CD8 $median.mat["High Entropy", "Unico"],
                              sim.pbmc.stats$boxplot.meta$Mono$median.mat["High Entropy", "Unico"],
                              sim.pbmc.stats$boxplot.meta$B   $median.mat["High Entropy", "Unico"]))

TCA.pbmc.high.corr = mean(c(sim.pbmc.stats$boxplot.meta$CD4 $median.mat["High Entropy", "TCA"],
                            sim.pbmc.stats$boxplot.meta$NK  $median.mat["High Entropy", "TCA"],
                            sim.pbmc.stats$boxplot.meta$CD8 $median.mat["High Entropy", "TCA"],
                            sim.pbmc.stats$boxplot.meta$Mono$median.mat["High Entropy", "TCA"],
                            sim.pbmc.stats$boxplot.meta$B   $median.mat["High Entropy", "TCA"]))

pbmc.low.improve = (Unico.pbmc.low.corr - TCA.pbmc.low.corr)/TCA.pbmc.low.corr
pbmc.low.improve

pbmc.high.improve = (Unico.pbmc.high.corr - TCA.pbmc.high.corr)/TCA.pbmc.high.corr
pbmc.high.improve

Unico.lung.low.corr = mean(c(sim.lung.stats$boxplot.meta$Immune     $median.mat["Low Entropy", "Unico"],
                             sim.lung.stats$boxplot.meta$Epithelial $median.mat["Low Entropy", "Unico"],
                             sim.lung.stats$boxplot.meta$Endothelial$median.mat["Low Entropy", "Unico"],
                             sim.lung.stats$boxplot.meta$Stroma     $median.mat["Low Entropy", "Unico"]))

TCA.lung.low.corr = mean(c(sim.lung.stats$boxplot.meta$Immune     $median.mat["Low Entropy", "TCA"],
                           sim.lung.stats$boxplot.meta$Epithelial $median.mat["Low Entropy", "TCA"],
                           sim.lung.stats$boxplot.meta$Endothelial$median.mat["Low Entropy", "TCA"],
                           sim.lung.stats$boxplot.meta$Stroma     $median.mat["Low Entropy", "TCA"]))


Unico.lung.high.corr = mean(c(sim.lung.stats$boxplot.meta$Immune     $median.mat["High Entropy", "Unico"],
                              sim.lung.stats$boxplot.meta$Epithelial $median.mat["High Entropy", "Unico"],
                              sim.lung.stats$boxplot.meta$Endothelial$median.mat["High Entropy", "Unico"],
                              sim.lung.stats$boxplot.meta$Stroma     $median.mat["High Entropy", "Unico"]))

TCA.lung.high.corr = mean(c(sim.lung.stats$boxplot.meta$Immune     $median.mat["High Entropy", "TCA"],
                            sim.lung.stats$boxplot.meta$Epithelial $median.mat["High Entropy", "TCA"],
                            sim.lung.stats$boxplot.meta$Endothelial$median.mat["High Entropy", "TCA"],
                            sim.lung.stats$boxplot.meta$Stroma     $median.mat["High Entropy", "TCA"]))

lung.low.improve = (Unico.lung.low.corr - TCA.lung.low.corr)/TCA.lung.low.corr
lung.low.improve

lung.high.improve = (Unico.lung.high.corr - TCA.lung.high.corr)/TCA.lung.high.corr
lung.high.improve

mean(c(lung.low.improve, pbmc.low.improve))

mean(c(lung.high.improve, pbmc.high.improve))

reinius.meth.stats = readRDS("../Figure/Purified-Reinius/hvf.10k/Reinius_eval_hvf.10k_stats.list.rds")

Unico.meth.RMedS = reinius.meth.stats$ RMedS.bar.dfs$all[reinius.meth.stats$ RMedS.bar.dfs$all$method == "Unico", ] 
bMIND.meth.RMedS = reinius.meth.stats$ RMedS.bar.dfs$all[reinius.meth.stats$ RMedS.bar.dfs$all$method == "bMIND", ] 

Unico.meth.RMedS

bMIND.meth.RMedS

(mean(bMIND.meth.RMedS$mean) - mean(Unico.meth.RMedS$mean))/mean(bMIND.meth.RMedS$mean)

Unico.meth.corr = reinius.meth.stats$cor.bar.dfs$all[reinius.meth.stats$ cor.bar.dfs$all$method == "Unico", ] 
bMIND.meth.corr = reinius.meth.stats$cor.bar.dfs$all[reinius.meth.stats$ cor.bar.dfs$all$method == "bMIND", ] 

Unico.meth.corr

bMIND.meth.corr

(mean(Unico.meth.corr$mean) - mean(bMIND.meth.corr$mean))/mean(bMIND.meth.corr$mean)

load("../Data//Methylation//Consistency//hannon1.processed.RData")
hannon1.feature.ids = rownames(hannon1$X)
rm(hannon1)
load("../Data//Methylation//Consistency//hannon2.processed.RData")
hannon2.feature.ids = rownames(hannon2$X)
rm(hannon2)
load("../Data//Methylation//Consistency//liu.processed.RData")
liu.feature.ids = rownames(liu$X)
rm(liu)
load("../Data//Methylation//Consistency//hannum.processed.RData")
hannum.feature.ids = rownames(hannum$X)
rm(hannum)

union.features.ids = union(union(hannon1.feature.ids, hannon2.feature.ids), 
                           union(liu.feature.ids, hannum.feature.ids))

print(length(union.features.ids))

# readRDS("../Result//Methylation/Consistency/XY/gender.meta.mats.F1.rds")
# readRDS("../Result//Methylation/Consistency/XY/gsex.meta.mats.F1.rds")

gender.meta.mats.MCC = readRDS("../Result//Methylation/Consistency/XY/gender.meta.mats.MCC.rds")
age.meta.mats.MCC    = readRDS("../Result//Methylation/Consistency/XY/age.meta.mats.MCC.rds")

gender.meta.mats.MCC$Unico.parametric.X2.marginal

gender.meta.mats.MCC$Baseline.parametric.full.joint

d = 4
Unico.gender.mcc = (sum(gender.meta.mats.MCC$Unico.parametric.X2.marginal) - d)/(d*(d-1))
bulk.gender.mcc  = (sum(gender.meta.mats.MCC$Baseline.parametric.full.joint) - d)/(d*(d-1))

Unico.age.mcc    = (sum(age.meta.mats.MCC$Unico.parametric.X2.marginal) - d)/(d*(d-1))
bulk.age.mcc     = (sum(age.meta.mats.MCC$Baseline.parametric.full.joint) - d)/(d*(d-1))

Unico.gender.mcc

bulk.gender.mcc

gender.meta.mats.MCC$Baseline.parametric.full.joint

(Unico.gender.mcc - bulk.gender.mcc)/bulk.gender.mcc

Unico.age.mcc 

bulk.age.mcc

(Unico.age.mcc - bulk.age.mcc)/bulk.age.mcc

reinius = readRDS("../Data/Methylation/Purified-Reinius/reinius.rds")

#reinius = readRDS("../Data/Methylation/Purified-Reinius/")

str(reinius)

length(liu.feature.ids)

length(hannum.feature.ids)

length(hannon1.feature.ids)

length(hannon2.feature.ids)

pbmc.pb = readRDS("../Data/RNA/Simulation-PBMC/pbmc.pseudobulk.decon.L1.all.W.rds")

str(pbmc.pb)

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "Mono", "B"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "Mono", "CD4"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "Mono", "CD8"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "Mono", "NK"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "B", "CD4"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "B", "CD8"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "B", "NK"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "CD4", "CD8"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "CD4", "NK"])

quantile(pbmc.pb$params$corrs[pbmc.pb$HEF, "CD8", "NK"])

data.names = c("liu", "hannon1", "hannon2", "hannum")

marg.method.names.full = c("Unico.parametric.marginal", "TCA.parametric.marginal", 
                           "CellDMC.parametric.marginal", "Baseline.parametric.marginal", 
                           "bMIND.YX.marginal") 



marg.method.names = c("Unico", "TCA", "CellDMC", "Baseline", "bMIND")

joint.method.names.full = c("Unico.parametric.joint", "TCA.parametric.joint", 
                            "CellDMC.parametric.joint", "Baseline.parametric.joint", 
                            "bMIND.XY.joint") 
joint.method.names = c("Unico", "TCA", "CellDMC", "Bulk", "bMIND")




partial.W  = FALSE
source.col = "decon.L1"

pbmc.pb = readRDS(file.path("../Data/RNA/Simulation-PBMC/", 
                            paste0("pbmc.pseudobulk.",
                                    source.col, ".",
                                    if (partial.W) "partial.W." else "all.W.",
                                    "rds")))

round(colMeans(pbmc.pb$W)[order (-colMeans(pbmc.pb$W))], 3) * 100

partial.W  = FALSE
source.col = "decon.L2"

pbmc.pb = readRDS(file.path("../Data/RNA/Simulation-PBMC/", 
                            paste0("pbmc.pseudobulk.",
                                    source.col, ".",
                                    if (partial.W) "partial.W." else "all.W.",
                                    "rds")))

round(colMeans(pbmc.pb$W)[order (-colMeans(pbmc.pb$W))], 3) * 100

partial.W  = FALSE
source.col = "decon.L1"

HLCA.pb = readRDS(file.path("../Data/RNA/Simulation-Lung/", 
                            paste0("HLCA.pseudobulk.",
                                    source.col, ".",
                                    if (partial.W) "partial.W." else "all.W.",
                                    "rds")))

round(colMeans(HLCA.pb$W)[order (-colMeans(HLCA.pb$W))], 3) * 100

partial.W  = FALSE
source.col = "decon.L2"

HLCA.pb = readRDS(file.path("../Data/RNA/Simulation-Lung/", 
                            paste0("HLCA.pseudobulk.",
                                    source.col, ".",
                                    if (partial.W) "partial.W." else "all.W.",
                                    "rds")))

round(colMeans(HLCA.pb$W)[order (-colMeans(HLCA.pb$W))], 3) * 100
