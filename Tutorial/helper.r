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

# Z: an array of c(sources, features, samples). The raw source tensor 
# max_sds: numeric. the number of standard deviations for which if outside is deemed as outliers and removed before consider if a feature is variable
# return a matrix of booleans: number of feautures by number of sources. 
# features-source that have mean and variance larger than active_thr are set to TRUE. 
calc_variable_feature_source = function(Z, variable_thr = 10**(-4), max_sds = 2){
    params = calc_params_from_Z(Z, max_sds = max_sds)
    variable.feature.source = matrix(FALSE, dim(Z)[2], dim(Z)[1])
    rownames(variable.feature.source) = dimnames(Z)[[2]]
    colnames(variable.feature.source) = dimnames(Z)[[1]]
    
    for (j in 1:dim(Z)[2]){
        for (l in 1:dim(Z)[1]){
            variable.feature.source[j,l] = (params$mus[j,l] > variable_thr) & (params$sigmas[j,l,l] > variable_thr) 
        }
    }
    return(variable.feature.source)
}

# x: list of numbers  
# y: list of numbers 
# robust: boolean. indicate if using robust correlation
# qtl: numeric between (0,1), only meaningful when robust is turned on. indicate the fraction of the data that is considered to be "good portion" and will participate in the correlation calculation
# check the correlation between x and y after excluding outliers defined by the qtl
# note that if in IQR is 0 in either x or y, OR x and y are completely colinear, robust correlation will be NA and thus filled by 0
safe_cor = function(x, y, robust = TRUE, qtl = 0.95){
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