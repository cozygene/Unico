library("MIND")
library("TCA")
source("Unico.r")
source("analysis.utils.r")

set.seed(2023)

#feature.set = "hvf.10k"
feature.set = "random.10k"

data.dir = paste0("../Data//Methylation/Purified-Reinius/", feature.set)
res.dir  = paste0("../Result//Methylation/Purified-Reinius/", feature.set)

if (!file.exists(res.dir)){dir.create(file.path(res.dir),recursive = T)}
print(data.dir)
print(res.dir)

hannum  = readRDS(file.path(data.dir, paste0("hannum.", feature.set, ".rds")))
reinius = readRDS(file.path(data.dir, paste0("reinius.", feature.set, ".rds")))

Unico.mdl = list()

Unico.mdl$params.hat <- Unico(hannum$X, hannum$W, C1 = NULL, C2 = NULL, 
                              parallel = TRUE, num_cores = NULL, 
                              log_file = file.path(res.dir, paste0("Unico.log")))

#capping 
Unico.mdl$params.hat$sigmas_hat  = cap_values(Unico.mdl$params.hat$sigmas_hat,  max.val = 10 ** (4), min.val = 10**(-4))

#tensor 
Unico.mdl$Z.hat = tensor(X = reinius$X, W = reinius$W, C1 = NULL, C2 = NULL, 
                        Unico.mdl$params.hat, parallel = FALSE)


#save
saveRDS(Unico.mdl,  file.path(res.dir, paste0("Unico.mdl.rds")))

str(Unico.mdl)

tca.mdl = list()
tca.mdl$params.hat <- tca(hannum$X, hannum$W,  
                          constrain_mu = TRUE, log_fil=file.path(res.dir, paste0("TCA.log")))

#capping 
tca.mdl$params.hat$sigmas_hat = cap_values(tca.mdl$params.hat$sigmas_hat, 
                                           max.val = 10**(4), min.val = 10**(-4))

#mimic the structure to allow TCA to treat reinius as if trained on it
tca.mdl.params.hat    = copy(tca.mdl$params.hat)
tca.mdl.params.hat$W  = reinius$W
tca.mdl.params.hat$C1 = matrix(0, ncol(reinius$X), 0)
tca.mdl.params.hat$C2 = matrix(0, ncol(reinius$X), 0)

#tensor
tca.mdl$Z.hat = TCA::tensor(X = reinius$X, tca.mdl.params.hat, log_fil=file.path(res.dir, paste0("TCA.log")))
tca.mdl$Z.hat = list_2_array(tca.mdl$Z.hat, colnames(reinius$W))  

#save
saveRDS(tca.mdl,  file.path(res.dir, paste0("tca.mdl.rds")))

str(tca.mdl)

base.mdl = list()

base.mdl$Z.hat = copy(tca.mdl$Z.hat) 
base.mdl$Z.hat[,,] = 0

print("constructing simple Z by distribute X by W")
for (h in 1:dim(base.mdl$Z.hat)[1]){
    base.mdl$Z.hat[h,,] = reinius$X * repmat(t(as.matrix(reinius$W[,h])), dim(base.mdl$Z.hat)[2], 1)
} 

saveRDS(base.mdl,  file.path(res.dir, paste0("base.mdl.rds")))

str(base.mdl)

# concat 
bMIND.X = cbind(reinius$X, hannum$X) 
bMIND.W = rbind(reinius$W, hannum$W) 

# remove . in sample name 
colnames(bMIND.X) = paste0("sample", 1:ncol(bMIND.X))
rownames(bMIND.W) = paste0("sample", 1:nrow(bMIND.W))


b = bMIND(bMIND.X, frac = bMIND.W)

bMIND.X

bMIND.W 

#format result
bMIND.mdl = list()

#params at original scale
bMIND.mdl$params.hat.orig = list(mus_hat = b$mu,
                                 sigmas_hat = b$Sigma_c)

#tensor at original scale
bMIND.mdl$Z.hat = array(0, c(ncol(reinius$W), nrow(reinius$X), ncol(reinius$X))) #source by feature by sample 
for(h in 1:ncol(reinius$W)){
    bMIND.mdl$Z.hat[h,,] = b$A[, h, paste0("sample", 1:ncol(reinius$X))]
}
dimnames(bMIND.mdl$Z.hat)[[1]] =  colnames(reinius$W)
dimnames(bMIND.mdl$Z.hat)[[2]] =  rownames(reinius$X)
dimnames(bMIND.mdl$Z.hat)[[3]] =  paste0("sample.", 1:ncol(reinius$X))

str(bMIND.mdl)

bMIND.mdl$Z.hat[1,,]

bMIND.mdl$Z.hat[5,,]

saveRDS(bMIND.mdl,  file.path(res.dir, paste0("bMIND.mdl.rds")))

 


