library(data.table)
source("analysis.utils.r")

set.seed(2023)

feature.set = "hvf.10k"
#feature.set = "random.10k"

robust = T
qtl = 0.95

data.dir = paste0("../Data//Methylation/Purified-Reinius/", feature.set)
res.dir  = paste0("../Result//Methylation/Purified-Reinius/", feature.set)

if (!file.exists(res.dir)){print("no result yet")}
print(data.dir)
print(res.dir)

hannum  = readRDS(file.path(data.dir, paste0("hannum.", feature.set, ".rds")))
reinius = readRDS(file.path(data.dir, paste0("reinius.", feature.set, ".rds")))



max_stds = 2
#max_stds = 3


base.mdl      = readRDS(file.path(res.dir, paste0("base.mdl.rds")))
cibersortx.mdl= readRDS(file.path(res.dir, paste0("cibersortx.mdl.rds")))
tca.mdl       = readRDS(file.path(res.dir, paste0("tca.mdl.rds")))
Unico.mdl      = readRDS(file.path(res.dir, paste0("Unico.mdl.rds")))
bMIND.mdl      = readRDS(file.path(res.dir, paste0("bMIND.mdl.rds")))

#after adjust for source, feature specific mean, 
#collapse all features to one meta feature
#return a list with two key, ref and hat. 
#each key is just a long meta vector, essentially one meta feature for this source
calc_centered_Z_h = function(Z.ref, Z.hat, h, mask){
    n = dim(Z.ref)[3]
    mask = as.vector(mask)
    
    #source, feature specific mean
    mus.ref = as.matrix(rowMeans(Z.ref[h,mask,]))
    mus.hat = as.matrix(rowMeans(Z.hat[h,mask,]))
    
    #collapse
    ref = as.vector(Z.ref[h,mask,] - repmat(mus.ref, 1, n))
    hat = as.vector(Z.hat[h,mask,] - repmat(mus.hat, 1, n))
    return (data.frame(ref = ref, hat = hat))  
}

#Z.ref: ground truth tensor with propor dim names: #sources by #features by #samples
#Z.hat: estimated tensor with propor dim names: #sources by #features by #samples
#group.mat:  matrix. number of features by number of different stratification.
#            each column should be booleans indicating if we should keep/True or exclude/false under that column name
#qtl: a number between 0 and 1. indicating the good portion of the data 
#returns a list with 5 elements: robcorrelation, mean and median absolute error, root mean and median square error
#each element is number of sources by number of stratification. essetially we have one meta feature per source now
calc_center_Z_hat_metrics = function(Z.ref, Z.hat, group.mat, qtl){
    
    center.Z.corrs           = matrix(0, dim(Z.ref)[1], ncol(group.mat))
    rownames(center.Z.corrs) = dimnames(Z.ref)[[1]]
    colnames(center.Z.corrs) = colnames(group.mat)  
   
    
    center.Z.MAE   = copy(center.Z.corrs)
    center.Z.MedAE = copy(center.Z.corrs)
    center.Z.RMS   = copy(center.Z.corrs)
    center.Z.RMedS = copy(center.Z.corrs)
    
    
    for (h in 1:dim(Z.ref)[1]){
        for (mask.i in 1:ncol(group.mat)){
           
            res = calc_centered_Z_h(Z.ref, Z.hat, h, mask = group.mat[,mask.i])
            ref = res$ref
            hat = res$hat
           
            center.Z.corrs[h,mask.i] = safe_cor(ref, hat, robust = T, qtl = qtl)
           
            center.Z.MAE[h,mask.i]   = mean(abs(ref - hat))
            center.Z.MedAE[h,mask.i] = median(abs(ref - hat))
            
            center.Z.RMS[h,mask.i]   = sqrt(mean((ref - hat)**2))
            center.Z.RMedS[h,mask.i] = sqrt(median((ref - hat)**2))
        }
    }
    
    metrics = list(center.Z.corrs = center.Z.corrs,
                   center.Z.MAE  = center.Z.MAE, center.Z.MedAE = center.Z.MedAE,
                   center.Z.RMS  = center.Z.RMS, center.Z.RMedS = center.Z.RMedS)
   
    return(metrics)
}

#Z.ref: ground truth tensor with propor dim names: #sources by #features by #samples
#Z.hat: estimated tensor with propor dim names: #sources by #features by #samples
#group.mat:  matrix. number of features by number of different stratification.
#            each column should be booleans indicating if we should keep/True or exclude/false under that column name
#qtl: a number between 0 and 1. indicating the good portion of the data 
#sample.itr: number of sampling to do 
#sample.size: number of samples to draw at each iteration 
#seed.number: default to 2023
#return a list, each entry is a random sample's 5 metrics caluclated using calc_center_Z_hat_metrics
calc_center_Z_hat_metrics_CI = function(Z.ref, Z.hat, group.mat, qtl, sample.itr, sample.size, seed.number = 2023){
    set.seed(seed.number)
    
    feature.ids = dimnames(Z.ref)[[2]]
    if (sample.itr*sample.size > length(feature.ids)){
        message("sample with replacement")
        feature.ids = sample(feature.ids, sample.itr * sample.size, replace = T)
    }else{
        message("sample without replacement")
        feature.ids = sample(feature.ids, sample.itr * sample.size, replace = F)
    }

    features.list = list()
    metrics.list = list()
    for (t in 1:sample.itr){
        features.list[[t]] = feature.ids[((t-1)*sample.size + 1):(t*sample.size)]
        metrics.list[[t]]  = calc_center_Z_hat_metrics(Z.ref[,features.list[[t]],], 
                                                       Z.hat[,features.list[[t]],], 
                                                       group.mat[features.list[[t]], ], qtl)
    }
    res = list(features.list = features.list, 
               metrics.list  = metrics.list)
    
    return(res)
}

str(reinius)



feature.ids = dimnames(reinius$Z)[[2]]
source.ids  = dimnames(reinius$Z)[[1]]

group.mat = matrix(FALSE, length(feature.ids), 3)
rownames(group.mat) = feature.ids
colnames(group.mat) = c("low entropy", "high entropy", "all")
group.mat[,"low entropy"]  = reinius$params$entropies <= quantile(reinius$params$entropies, 0.5)
group.mat[,"high entropy"] = reinius$params$entropies >  quantile(reinius$params$entropies, 0.5)
group.mat[,"all"] = TRUE

group.mat

colSums(group.mat)

calc_centered_Z_h(Z.ref = reinius$Z, Z.hat = base.mdl$Z.hat, h = 1, mask  = group.mat[,3])

base.mdl$evals = calc_center_Z_hat_metrics_CI(Z.ref = reinius$Z, 
                                              Z.hat = base.mdl$Z.hat, 
                                              group.mat = group.mat, qtl = qtl,
                                              sample.itr = 20, sample.size = 1000)

cibersortx.mdl$evals = calc_center_Z_hat_metrics_CI(Z.ref = reinius$Z, 
                                                    Z.hat = cibersortx.mdl$Z.hat, 
                                                    group.mat = group.mat, qtl = qtl,
                                                    sample.itr = 20, sample.size = 1000)

tca.mdl$evals = calc_center_Z_hat_metrics_CI(Z.ref = reinius$Z, 
                                              Z.hat = tca.mdl$Z.hat, 
                                              group.mat = group.mat, qtl = qtl,
                                              sample.itr = 20, sample.size = 1000)

Unico.mdl$evals = calc_center_Z_hat_metrics_CI(Z.ref = reinius$Z, 
                                              Z.hat = Unico.mdl$Z.hat, 
                                              group.mat = group.mat, qtl = qtl,
                                              sample.itr = 20, sample.size = 1000)






bMIND.mdl$evals = calc_center_Z_hat_metrics_CI(Z.ref = reinius$Z, 
                                              Z.hat = bMIND.mdl$Z.hat, 
                                              group.mat = group.mat, qtl = qtl,
                                              sample.itr = 20, sample.size = 1000)


base.mdl$evals$metrics.list[[1]]$center.Z.corrs

cibersortx.mdl$evals$metrics.list[[1]]$center.Z.corrs

tca.mdl$evals$metrics.list[[1]]$center.Z.corrs

Unico.mdl$evals$metrics.list[[1]]$center.Z.corrs

bMIND.mdl$evals$metrics.list[[1]]$center.Z.corrs

base.mdl$evals$metrics.list[[1]]$center.Z.RMedS

cibersortx.mdl$evals$metrics.list[[1]]$center.Z.RMedS

tca.mdl$evals$metrics.list[[1]]$center.Z.RMedS

Unico.mdl$evals$metrics.list[[1]]$center.Z.RMedS

bMIND.mdl$evals$metrics.list[[1]]$center.Z.RMedS

saveRDS(base.mdl,       file.path(res.dir, paste0("base.mdl.rds")))
saveRDS(cibersortx.mdl, file.path(res.dir, paste0("cibersortx.mdl.rds")))
saveRDS(tca.mdl,        file.path(res.dir, paste0("tca.mdl.rds")))
saveRDS(Unico.mdl,       file.path(res.dir, paste0("Unico.mdl.rds")))
saveRDS(bMIND.mdl,      file.path(res.dir, paste0("bMIND.mdl.rds")))




