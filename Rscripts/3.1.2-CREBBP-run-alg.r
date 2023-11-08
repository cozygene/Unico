# This notebook will iteratate thru all configures and  
# run baseline, bulk, Unico, TCA under the mixture X and W.hat estiamted by CIBERSORTx
# save the result to Result/RNA/CREBBP/ with configuration name "key" added to the file name
library("MIND")
library("TCA")
source("Unico.r")
source("analysis.utils.r")

set.seed(2023)  

data.dir = file.path("../Data/RNA/CREBBP/")
res.dir  = file.path("../Result/RNA/CREBBP/")
CREBBP.dat = readRDS(file.path(data.dir, "CREBBP.dat.rds"))

#Unico algorithm parameters
max_stds = 2

configures = list(c("genes",  F, F))

str(CREBBP.dat$expm.dat[["genes"]])


version = "genes"
key = version
    
W  = CREBBP.dat$expm.dat[[key]]$W
X  = CREBBP.dat$expm.dat[[key]]$X
C1 = CREBBP.dat$expm.dat[[key]]$C1
C2 = CREBBP.dat$expm.dat[[key]]$C2
source.ids  = colnames(W)
feature.ids = rownames(X)
sample.ids  = rownames(W)
k = length(source.ids)
m = length(feature.ids)
n = length(sample.ids)
     

############ Unico ###########
Unico.log.file = file.path(res.dir, paste0("Unico.mdl.", key, ".log"))
Unico.res.file = file.path(res.dir,paste0("Unico.mdl.", key,  ".rds"))
Unico.mdl = list()
Unico.mdl$params.hat <- Unico(X  = X, W  = W, C1 = C1, C2 = C2,
                              parallel = TRUE, num_cores = NULL, 
                              log_file = Unico.log.file)

#capping at interal scale
Unico.mdl$params.hat$sigmas_hat= cap_values(Unico.mdl$params.hat$sigmas_hat,  max.val = 10 ** (4), min.val = 10**(-4))
#hat internal scaled version
Unico.mdl$params.hat$corrs     = calc_corrs_from_sigmas(Unico.mdl$params.hat$sigmas_hat)
Unico.mdl$params.hat$entropies = calc_entropy(Unico.mdl$params.hat$corrs)


#tensor at original scale
Unico.mdl$Z.hat = tensor(X  = X,  
                        W  = W, 
                        C1 = C1, 
                        C2 = C2,
                        Unico.mdl$params.hat, 
                        parallel = FALSE)
saveRDS(Unico.mdl, Unico.res.file)
             
############ TCA ###########            
tca.log.file = file.path(res.dir, paste0("tca.mdl.", key, ".log"))
tca.res.file = file.path(res.dir,paste0("tca.mdl.", key, ".rds"))
tca.mdl = list()
tca.mdl$params.hat <- tca(X  = X, 
                          W  = W, 
                          C1 = C1,  
                          C2 = C2, 
                          verbose = TRUE, parallel = F, num_cores = 1, log_file = tca.log.file)
             
             
tca.mdl$Z.hat  <- TCA::tensor(X = as.matrix(CREBBP.dat$expm.dat[[key]]$X), 
                              tca.mdl$params.hat,
                              parallel = FALSE, num_cores = 1, log_file=tca.log.file)
             
tca.mdl$Z.hat = list_2_array(tca.mdl$Z.hat, source.ids)
dimnames(tca.mdl$Z.hat)[[1]] = colnames(CREBBP.dat$expm.dat[[key]]$W)
dimnames(tca.mdl$Z.hat)[[2]] = rownames(CREBBP.dat$expm.dat[[key]]$X)
dimnames(tca.mdl$Z.hat)[[3]] = colnames(CREBBP.dat$expm.dat[[key]]$X)
saveRDS(tca.mdl, tca.res.file)


############ Base ###########      
base.res.file = file.path(res.dir,paste0("base.mdl.", key, ".rds"))       
base.mdl                 = list()
base.mdl$Z.hat = array(0, c(k, m, n))
dimnames(base.mdl$Z.hat)[[1]] = source.ids
dimnames(base.mdl$Z.hat)[[2]] = feature.ids
dimnames(base.mdl$Z.hat)[[3]] = sample.ids
for(l in 1:k){
    base.mdl$Z.hat[l,,] = X[feature.ids, sample.ids] * repmat(t(as.matrix(W[sample.ids,l])), m,1)
}
saveRDS(base.mdl, base.res.file)
             
############ Bulk ###########  
bulk.res.file = file.path(res.dir,paste0("bulk.mdl.", key, ".rds"))    
bulk.mdl = list()
bulk.mdl$Z.hat = array(0, c(k, m, n))
dimnames(bulk.mdl$Z.hat)[[1]] = source.ids
dimnames(bulk.mdl$Z.hat)[[2]] = feature.ids
dimnames(bulk.mdl$Z.hat)[[3]] = sample.ids
for(l in 1:k){
    bulk.mdl$Z.hat[l,,] = X[feature.ids, sample.ids]
}
saveRDS(bulk.mdl, bulk.res.file)
 
             
############ bMIND.log ###########    
b = bMIND(log2(X + 1), # as recommended by the aurthor 
       W)        
#format result
bMIND.mdl = list()

#params at original scale "logged"
bMIND.mdl$params.hat = list(mus_hat = b$mu,
                         sigmas_hat = b$Sigma_c)
#tensor at original scale "logged"
bMIND.mdl$Z.hat = array(0, c(ncol(W), nrow(X), ncol(X)))
dimnames(bMIND.mdl$Z.hat)[[1]] = dimnames(b$A)[[2]]
dimnames(bMIND.mdl$Z.hat)[[2]] = dimnames(b$A)[[1]]
dimnames(bMIND.mdl$Z.hat)[[3]] = dimnames(b$A)[[3]]

for(h in 1:ncol(W)){
    bMIND.mdl$Z.hat[h,,] = 2 ** (b$A[,h,]) # make it back to original space
}     
bMIND.res.file = file.path(res.dir,paste0("bMIND.log.mdl.", key, ".rds"))       
saveRDS(bMIND.mdl, bMIND.res.file)

         
         
             
############ bMIND.scale ###########          
b = bMIND(X/(max(X)/50), # hack to make it less than 50
       W)
         #format result
bMIND.mdl = list()
#params at original scale 
bMIND.mdl$params.hat = list(mus_hat = b$mu,
                         sigmas_hat = b$Sigma_c)
#tensor at original scale
bMIND.mdl$Z.hat = array(0, c(ncol(W), nrow(X), ncol(X)))
dimnames(bMIND.mdl$Z.hat)[[1]] = dimnames(b$A)[[2]]
dimnames(bMIND.mdl$Z.hat)[[2]] = dimnames(b$A)[[1]]
dimnames(bMIND.mdl$Z.hat)[[3]] = dimnames(b$A)[[3]]

for(h in 1:ncol(W)){
 bMIND.mdl$Z.hat[h,,] = b$A[,h,] * (max(X)/50)
}
         
bMIND.res.file = file.path(res.dir,paste0("bMIND.scale.mdl.", key, ".rds"))       
saveRDS(bMIND.mdl, bMIND.res.file)
                        


bMIND.mdl

# shrink and expand by (max(X)/50)
bMIND.mdl$Z.hat[1,,]




