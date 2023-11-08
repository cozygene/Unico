# This notebook read the cibersortx estiamted tensor from local laptop upload under "Result/RNA/CREBBP/Cibersortx/"
# will iterate thru all the configurations/folder, currently just the "genes"
library("data.table")
set.seed(2023)

data.dir   = file.path("../Data/RNA/CREBBP/")
res.dir    = file.path("../Result/RNA/CREBBP/")

CREBBP.dat = readRDS(file.path(data.dir, "CREBBP.dat.rds"))
#all X have the same order of samples
#all DE genes or probes are unique and present in the corresponding X

str(CREBBP.dat)

source.ids   = colnames(CREBBP.dat$expm.dat$genes$W)
feature.ids = rownames(CREBBP.dat$expm.dat$genes$X)
sample.ids   = colnames(CREBBP.dat$expm.dat$genes$X)

cibersortx.dir = file.path(res.dir, "Cibersortx")

version = "genes"
key = version
print(paste0("working on: ",key))
    
cibersortx.mdl = list()
print("read in cibersortx's estimated Z")
Z.hat = array(0, c(length(source.ids), length(feature.ids), length(sample.ids)))
dimnames(Z.hat)[[1]] = source.ids
dimnames(Z.hat)[[2]] = feature.ids
dimnames(Z.hat)[[3]] = sample.ids

for (source.id in source.ids){
    Z.hat.file = file.path(cibersortx.dir, key, "res", paste0("CIBERSORTxHiRes_NA_", gsub("\\ ", "", source.id) , "_Window",round(4 * length(source.ids)), ".txt"))
    Z.hat[source.id,,] = as.matrix(data.frame(fread(Z.hat.file), row.names=1))
}
cibersortx.mdl$Z.hat = Z.hat
saveRDS(cibersortx.mdl, file.path(res.dir, paste0("cibersortx.mdl.", key, ".rds"))       


