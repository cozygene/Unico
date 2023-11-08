library("MASS")   # for robust correlation
library("pracma") # for eigen

#sigmas: an array of m features by k sources by k sources, indicating the k by k covariance matrix for each feature
#return the corresponding correlation matrices for all features: an array of m features by k sources by k sources
calc_corrs_from_sigmas = function(sigmas){
    
    m = dim(sigmas)[1]
    k = dim(sigmas)[2]
    
    corrs = array(0, c(m, k, k))
    for (j in 1:m){
        corrs[j,,]  = sigmas[j,,]/sqrt(diag(sigmas[j,,]) %*% t(diag(sigmas[j,,])))
    }
    dimnames(corrs)[[1]] = dimnames(sigmas)[[1]]
    dimnames(corrs)[[2]] = dimnames(sigmas)[[2]]
    dimnames(corrs)[[3]] = dimnames(sigmas)[[3]]
    
    #if there are no variation in one celltype. corresponding sigmas will be 0, related cor will be NA
    #here set those to 0
    corrs[is.na(corrs)] = 0
    return(corrs)
}

#ltri_val: numbers in the lower triangle of a matrix
#k: the dimension (number of rows or cols) of the full symmetric matrix 
#return a k by k symmetric matrix constructed from the ltri_val
fill_lower_tri = function(ltri_val, k){
    full_mat <- matrix(0,k,k)
    full_mat[lower.tri(full_mat,diag = T)] <- ltri_val
    full_mat <- t(full_mat)
    full_mat[lower.tri(full_mat,diag = T)] <- ltri_val
    return (full_mat)
}

#Z: an array of k sources by m features by n samples
#max_sds: numeric. Per (source,feature), samples whose value fall outside of max_sds will be ignored
#returns a list with the following entires: mus (m by k), sigmas (m by k by k), corrs (m by k by k), entropies (m by 1), 
#for features that are not expressed in a source. The varaince term on the diagnol is set to 0
#the correlation term is undefined, also set to 0
calc_params_from_Z <- function(Z, max_sds){
    
    k <- dim(Z)[1]
    m <- dim(Z)[2]
    n <- dim(Z)[3] 
    source.ids  = dimnames(Z)[[1]]
    feature.ids = dimnames(Z)[[2]]
    
    means  <- matrix(0,m,k)
    rownames(means) = feature.ids 
    colnames(means) = source.ids
    
    sigmas <- matrix(0,m,(choose(k,2) + k))
    for (j in 1:m){
        counter <- 1
        for (h in 1:k){    
            sd_jh <- sd(Z[h,j,])
            indices_jh <- abs(Z[h,j,] - mean(Z[h,j,])) <= max_sds*sd_jh
            #no such sample left
            if (sum(indices_jh) < 1){
                means[j,h] <- 0
            }else{
                means[j,h] <- mean(Z[h,j,indices_jh])
            }
            for (l in h:k){
                sd_jl <- sd(Z[l,j,])
                indices_jl <- abs(Z[l,j,]-mean(Z[l,j,])) <= max_sds*sd_jl
                indices_jhl <- indices_jh & indices_jl
                #only 1 or no such sample left
                if ((sum(indices_jhl) <=1)|(sd(Z[h,j,indices_jhl]) == 0) | (sd(Z[l,j,indices_jhl]) == 0)){
                    sigmas[j, counter] = 0
                }else{
                    sigmas[j, counter] <- cov(Z[h,j,indices_jhl],Z[l,j,indices_jhl])
                }
                counter <- counter + 1
            }
        }  
    } 
    
    sigmas.tensor = array(0,c(m,k,k))
    dimnames(sigmas.tensor)[[1]] = feature.ids 
    dimnames(sigmas.tensor)[[2]] = source.ids
    dimnames(sigmas.tensor)[[3]] = source.ids
   
    for (j in 1:m){
        sigmas.tensor[j,,] = fill_lower_tri(sigmas[j,], k)
    }
    corrs.tensor  = calc_corrs_from_sigmas(sigmas.tensor)
    
    
    params = list()
    params$sigmas    = sigmas.tensor
    params$corrs     = corrs.tensor
    params$mus       = means
    params$entropies = calc_entropy(params$corrs)
    return(params)
}

# corrs: an array of shape m features by k sources by k sources
# return a matrix of m features by 1
# per feature, calculate entropy given correlations (k by k) matrix
calc_entropy <- function(corrs){
    m <- dim(corrs)[1]
    k <- dim(corrs)[2]
    entropies <- numeric(m)
    for (j in 1:m){
        sigma_j <- corrs[j,,]
        # von Neumann entropy for Pearson correlation matrix (https://arxiv.org/pdf/2106.05379.pdf)
        # make sure it is a valid full rank correlation matrix
        diag(sigma_j) = 1
        sigma_j_eig <- tryCatch({
            eig(sigma_j/k)
        },error=function(cond){
            print(j)
        }) 
    
    #this is trying to catch the case when the covariance matrix is not PSD
    #which might happen because of the way we are constructing the covariance matrix (by exluding different outlier samples)
    sigma_j_eig[sigma_j_eig<=0] <- 1e-5
    entropies[j] <- -sum(sigma_j_eig*log(sigma_j_eig))  
    } 
    
    entropies = as.matrix(entropies)
    rownames(entropies) = dimnames(corrs)[[1]]
    colnames(entropies) = "entropy"
    return(entropies)
}

##############################################################################################################################

# data: any data structure, could be 2d matrix or 3d array 
# min.val: numeric. minimal absolute value to be capped 
# max.val: numeric. maximal absolute value to be capped 
# NA in the data structure is turned into + min.val
# absolute values larger than max.val or smaller than min.val (including Inf, -Inf)are replaced by max.val and min.val
# sign information is also preserved
cap_values = function(data, min.val = 10**(-4), max.val = 10**(4)){
    if(sum(is.na(data)| is.infinite(data)) != 0){
        print(paste0("there are NAN: ", sum(is.na(data))))
        data[is.na(data)] <- min.val
    }
    if(sum(abs(data) < min.val) != 0){
        print(paste0("there are extrmemely close to 0 values: ", sum(abs(data) < min.val)))
        sign.info = sign(data[which(abs(data) < min.val)])
        data[which(abs(data) < min.val)] <- sign.info * min.val
    }
    if(sum(abs(data) > max.val) != 0){
        print(paste0("there are extrmemely large values: " , sum(abs(data) < max.val)))
        sign.info = sign(data[which(abs(data) > max.val)])
        data[which(abs(data) > max.val)] <- sign.info * max.val
    }
    return(data)
} 

# Z: TCA's tensor in the form of a list of matrices. each list element is a celltype. each matrix is features by samples
# source.ids: a list of charactors, indicating the sources/celltypes
# turn TCA's tensor Z to an 3d array format with all the dimnames 
list_2_array = function(Z, source.ids){
    
    k = length(Z)
    m = dim(Z[[1]])[1]
    n = dim(Z[[1]])[2]
    Z.array = array(0, c(k, m, n))
    
    dimnames(Z.array)[[1]] = source.ids
    dimnames(Z.array)[[2]] = rownames(Z[[1]])
    dimnames(Z.array)[[3]] = colnames(Z[[1]])
    
    for(h in 1:k){
        Z.array[h,,] = Z[[h]]
    }
    return(Z.array)
}

# mat: a matrix that is number of features by number of sources
# return a tensor that is number of features by number of sources by number of sources
# per feature, it turns the original vector (of length number of sources) to entries in a diagnol matrix
matrix_to_diag_tensor = function(mat){
    m = nrow(mat)
    k = ncol(mat)
    feature.ids = rownames(mat)
    source.ids  = colnames(mat)
    
    dtensor = array(0, c(m,k,k))
    dimnames(dtensor)[[1]] = feature.ids
    dimnames(dtensor)[[2]] = source.ids
    dimnames(dtensor)[[3]] = source.ids
    
    for (j in 1:m){
        dtensor[j,,] = diag(mat[j,])
    }
    return(dtensor) 
}

##############################################################################################################################

# x: list of numbers  
# y: list of numbers 
# robust: boolean. indicate if using robust correlation
# qtl: numeric between (0,1), only meaningful when robust is turned on. indicate the fraction of the data that is considered to be "good portion" and will participate in the correlation calculation
# check the correlation between x and y after excluding outliers defined by the qtl
# note that if in IQR is 0 in either x or y, OR x and y are completely colinear, robust correlation will be NA and thus filled by 0
safe_cor = function(x, y, robust = FALSE, qtl = 0.95){
    res = tryCatch({
        if(robust){
            cov.rob(cbind(x,y), cor = TRUE, quantile.used = round(qtl*length(x)) ,method = "mve")$cor[1,2] 
        }else{
            cor(x,y)
        }
        
    }, error = function(cond){
        0
    })
    
    # in regular correlation if all constant, then it is going to be NA
    if (is.na(res)){
        return(0)
    }else{
        return(res)
    }
    
}



# Z.true is an array of k by m by n, the true tensor 
# Z.hat is an array of k by m by n, the estimated tensor 
# eval.feature.source is a matrix of logical values that is number of features by number of sources. only calculate correlation on feature-source that has TRUE in this matrix 
# robust: boolean. indicate if using robust correlation
# qtl: numeric between (0,1), only meaningful when robust is turned on. indicate the fraction of the data that is considered to be "good portion" and will participate in the correlation calculation
# return a m by k correltion matrix, each column is a different source and each row is a different feature, 
# each entry is a correlation score, indicating that feautre, that source's tensor estimated across samples
# those entry with FALSE in the eval.feature.source is going to be set to NA 
calc_Z_corrs <- function(Z.true, Z.hat, eval.feature.source = NULL, robust = TRUE, qtl = 0.95){
    
    m = dim(Z.true)[2]
    k = dim(Z.true)[1]
    Z.corrs = matrix(NA, m, k)
    rownames(Z.corrs) = dimnames(Z.true)[[2]]
    colnames(Z.corrs) = dimnames(Z.true)[[1]]
    
    # if eval.feature.source is not present, return all results
    if(is.null(eval.feature.source)){
        eval.feature.source = matrix(TRUE, m, k)
    }
    
    for (h in 1:k){
        for (j in 1:m){
            if(eval.feature.source[j,h]){
                Z.corrs[j,h] <- safe_cor(Z.true[h,j,], Z.hat[h,j,], robust, qtl)
            }else{
                Z.corrs[j,h] <- NA
            }
                
        }
    }
    return(Z.corrs)
}

#source.ids: list of sources
#generate all the unique combination of 2 difference sources  
get_covar_ids = function(source.ids){
    k = length(source.ids)
    covar.ids = list(choose( k,2))
    counter = 1
    for (h1 in 1:( k -1)){
        for (h2 in (h1 + 1): k){
            covar.ids[counter] = paste(source.ids[h1], source.ids[h2], sep = "-")
            counter = counter +1
        }
    }
    return(unlist(covar.ids))
}

# params.true: a list of 2 keys: mus, sigmas. mus is number of features by number of sources. sigmas is number of features by source by source. these are true parameters
# params.hat: similar to params.true but are the estimated 
# robust: boolean if to use robust correlation. default is TRUE
# qtl: numeric, the quantile of robust correlation, default is 0.95
get_moment_corrs = function (params.true, params.hat, robust = TRUE, qtl = 0.95){
            
    feature.ids = rownames(params.true$mus)
    source.ids  = colnames(params.true$mus)
    covar.ids   = get_covar_ids(source.ids)
    
    mus.rob.corrs   = matrix(0, 1, length(source.ids))
    colnames(mus.rob.corrs)   = source.ids
    rownames(mus.rob.corrs)   = c("mus.rob.corrs")
    
    var.rob.corrs   = matrix(0, 1, length(source.ids))
    colnames(var.rob.corrs)   = source.ids
    rownames(var.rob.corrs)   = c("var.rob.corrs")
    
    covar.rob.corrs = matrix(0, 1, choose(length(source.ids), 2))
    colnames(covar.rob.corrs) = covar.ids
    rownames(covar.rob.corrs) = c("covar.rob.corrs")

    for (h in 1:length(source.ids)){
        mus.rob.corrs[1,h]  = safe_cor(params.true$mus[,h], 
                                       params.hat$mus[,h], robust, qtl)
        
        var.rob.corrs[1,h]  = safe_cor(params.true$sigmas[,h,h], 
                                       params.hat$sigmas[,h,h], robust, qtl)
    }
    for (covar.id in covar.ids){
        covar.pair = strsplit(covar.id,"-")[[1]]
        covar.rob.corrs[1,covar.id] = safe_cor(params.true$sigmas[,covar.pair[1],covar.pair[2]], 
                                               params.hat$sigmas[,covar.pair[1],covar.pair[2]], robust, qtl)
    }
                      
    return(list(mus.rob.corrs   = mus.rob.corrs,         
                var.rob.corrs   = var.rob.corrs,      
                covar.rob.corrs = covar.rob.corrs))                    
}



#across runs
##############################################################################################################################

# mdl.list: a list of models, each represent one run
# key: character, the name of the variable to extract. 
# return concated version of that varibale, stacked vertically. columns remain the same
concat_key = function(mdl.list, key){
    res = list(length(mdl.list))
    for (t in 1:length(mdl.list)){
        res[[t]] = mdl.list[[t]][[key]]
    }
    res = Reduce(rbind, res)
    return(res)
}

# mdl.list: a list of models, each represent one run
# key1: character, the first name of the variable to extract. 
# key2: character, the second name of the variable to extract. 
# return concated version of that varibale, stacked vertically. columns remain the same
concat_2_keys = function(mdl.list, key1, key2){
    res = list(length(mdl.list))
    for (t in 1:length(mdl.list)){
        res[[t]] = mdl.list[[t]][[key1]][[key2]]
    }
    res = Reduce(rbind, res)
    return(res)
}


##############################################################################################################################
#param: matrix with rows as features 
#feature.scale.factor: a matrix of number of features by 1. 
#return param multiplied on the feature level by the feature.scale.factor
scale_feature_matrix  = function(param, feature.scale.factor){
    return(param * repmat(feature.scale.factor, 1, ncol(param)))
}

#param: array of number of features by number of sources by number of sources
#feature.scale.factor: a matrix of number of features by 1. 
#return param multiplied on the feature level by the feature.scale.factor
scale_feature_source_source  = function(param, feature.scale.factor){
    param.scale = copy(param)
    for(j in 1:dim(param)[1]){
        param.scale[j,,] = param.scale[j,,]*(feature.scale.factor[j,])**2
    }
    return(param.scale)
}

#param: array of number of sources by number of features by number of samples
#feature.scale.factor: a matrix of number of features by 1. 
#return param multiplied on the feature level by the feature.scale.factor
scale_source_feature_sample  = function(param, feature.scale.factor){
    param.scale = copy(param)
    for(l in 1:dim(param)[1]){
        param.scale[l,,] = param.scale[l,,]*repmat(feature.scale.factor, 1, dim(param)[3])
    }
    return(param.scale)
}

# mat: a matrix that is number of features by number of sources
# return a tensor that is number of features by number of sources by number of sources
# per feature, it turns the original vector (of length number of sources) to entries in a diagnol matrix
matrix_to_diag_tensor = function(mat){
    m = nrow(mat)
    k = ncol(mat)
    feature.ids = rownames(mat)
    source.ids  = colnames(mat)
    
    dtensor = array(0, c(m,k,k))
    dimnames(dtensor)[[1]] = feature.ids
    dimnames(dtensor)[[2]] = source.ids
    dimnames(dtensor)[[3]] = source.ids
    
    for (j in 1:m){
        dtensor[j,,] = diag(mat[j,])
    }
    return(dtensor) 
}





# a slightly more reasonable way to make sure Z is non-negative while preserving the between sample relationship
####################################################################################

#Z: a 3D array of sources by features by samples
#return a new version with no negative values, per source, feature, if the minimal value is below zero, 
#shift the entire distribution across all samples by the - of the minimal value  
none_neg_Z = function(Z){
    k = dim(Z)[1]
    m = dim(Z)[2]
    n = dim(Z)[3]
    counter = 0
    
    for (l in 1:k){
        for (j in 1:m){
            sample.min = min(Z[l,j,])
            if(sample.min < 0){
                counter = counter + 1
                #move entire distribution to non-neg range 
                Z[l,j,] = Z[l,j,] + (-sample.min)
            }
        }
    }
    
    message(paste0(round((counter/(m * k)) * 100, 2), " percent of the feature-source are shifted to be non negative"))
    return(Z)
}


