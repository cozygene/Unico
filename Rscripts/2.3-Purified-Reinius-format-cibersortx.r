library(data.table)
set.seed(2023)

feature.set = "hvf.10k"
#feature.set = "random.10k"

data.dir = paste0("../Data//Methylation/Purified-Reinius/", feature.set)
res.dir  = paste0("../Result//Methylation/Purified-Reinius/", feature.set)

if (!file.exists(res.dir)){dir.create(file.path(res.dir),recursive = T)}
print(data.dir)
print(res.dir)

hannum  = readRDS(file.path(data.dir, paste0("hannum.", feature.set, ".rds")))
reinius = readRDS(file.path(data.dir, paste0("reinius.", feature.set, ".rds")))

colnames(reinius$W)

source.ids  = colnames(reinius$W)
sample.ids  = colnames(reinius$X)
feature.ids = rownames(reinius$X)

cibersortx.mdl = list()

print("read in cibersortx's estimated Z")
cibersortx.mdl$Z.hat = array(0, c(length(source.ids), length(feature.ids), length(sample.ids)))
dimnames(cibersortx.mdl$Z.hat)[[1]] = source.ids
dimnames(cibersortx.mdl$Z.hat)[[2]] = feature.ids 
dimnames(cibersortx.mdl$Z.hat)[[3]] = sample.ids

for (source.id in source.ids){
    Z.hat.file = file.path(res.dir, 
                           paste0("CIBERSORTxHiRes_NA_", 
                                  gsub("\\ ", "", source.id), 
                                  "_Window",round(4 * length(source.ids)),".",
                                  feature.set, ".txt"))
    #just take reinius and scale back
    cibersortx.mdl$Z.hat[source.id,,] = as.matrix(data.frame(fread(Z.hat.file), row.names=1))[feature.ids, sample.ids]/10000
}

length(feature.ids)

saveRDS(cibersortx.mdl,  file.path(res.dir, paste0("cibersortx.mdl.rds")))

str(cibersortx.mdl)

cibersortx.mdl$Z.hat[1,,]

cibersortx.mdl$Z.hat[5,,]




