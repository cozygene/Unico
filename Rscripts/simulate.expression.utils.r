library("testit")
library("data.table") # also has copy
library("compositions")

# W.src: a matrix with rows as samples and columns as source names. source names needs to be meaningful
# n: numeric.  number of sample
# k: numeric. number of sources. 
# fit.dirichlet: logical.  if T, it will fit a dirichlet distribution and draw new W from it. 
# seed.num: numeric. 

# This function pioritize randomly choosing without replcaement from sources with > 5% 
# if the request k is larger than the number of sources with > 5% aboundance, 
# it will pick the top k most aboundant sources to work on
# if sample n bigger than what's in W then check if fit.dicichlet or sample with replacement 
# esle sample without replacement 
# return simulated W, row scaled to sum up to 1, col names coems from W.src, rownames are generated automatically 
simulate_W = function(W.src, n, k, fit.dirichlet = FALSE, seed.num = 2023){
    set.seed(seed.num)
    assert(k <= ncol(W.src))
    
    # working on k
    if (k > sum(colMeans(W.src)>0.05)){
        W <- W.src[,order(-colMeans(W.src))[1:k]]
        print("taking the top most aboundant celltypes (certain celltypes will be <5%)")
    }else{
        W <- W.src[,colMeans(W.src)>0.05] # exclude lowly abundant cell types  
        W <- W[,sample(ncol(W))[1:k]]
        print("sampling in the more aboundant celltypes (>5%)")
    }
    # rescale and reorder
    W <- W/repmat(as.matrix(rowSums(W)), 1, k) 
    W <- W[,order(-colMeans(W))]
    source.ids = colnames(W) # for dirichlet
    

    # working on samples in W
    if (n > nrow(W)){
        print("n is larger than the number of individuals in W")
        if(fit.dirichlet){
            print("fitting dirichlet")
            dirichlet.alpha = fitDirichlet(W)$alpha
            W = rDirichlet.acomp(n,dirichlet.alpha)
            W = matrix(W, n, k)    
        } else{
            print("draw with replacement")
            W = W[sample(1:nrow(W), n, replace = T),]  
        }
    }else{
        print("n is smaller than the number of individuals in W")
        print("draw without replacement")
        W = W[sample(1:nrow(W), n, replace = F),]
    }
    
    # regardless use new sample.ids (Z with also be randomly drawn)
    sample.ids  = as.vector(vapply(1:n, function(i) paste0("sample.",i), character(1)))
    rownames(W) = sample.ids
    # use real celltype names
    colnames(W) = source.ids
    return(as.matrix(W))
}



# Z: an array of c(sources, features, samples). The raw source tensor 
# expression.qtl: numeric between (0, 1). the top quantiles deemed as highly expressed
# return a ranked list of charactors: features that have are in the top quantile of (sum of expression across sources, samples)
highly_expressed_features = function(Z, expression.qtl){
    m = dim(Z)[2]
    expression.level = matrix(0, m, 1) 
    rownames(expression.level) = dimnames(Z)[[2]]
    
    for (j in 1:m){
        expression.level[j,] = sum(Z[,j,])
    }
    print(paste0("extracting top quantile :", expression.qtl, 
                 " features, based on sum of expression in all sources, all samples"))
    he.mask = expression.level >= quantile(expression.level, expression.qtl)
    expression.level = expression.level[he.mask, , drop = F]
    expression.level = expression.level[order(expression.level, decreasing = T), ,drop = F]
    return(rownames(expression.level))
}

# Z: an array of c(sources, features, samples). The raw source tensor 
# max_sds: numeric. the number of standard deviations for which if outside is deemed as outliers and removed before consider if a feature is variable
# return a matrix of booleans: number of feautures by number of sources. 
# features-source that have mean and variance larger than active_thr are set to TRUE. 
calc_variable_feature_source = function(Z, variable_thr = 10**(-4), max_sds = 3){
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


# variable.feature.source: a matrix of booleans: number of feautures by number of sources. entries that are set to FALSE indicate no expression 
# max_sds: numeric. the number of standard deviations for which if outside is deemed as outliers and removed before consider if a feature is variable
# return a list of charactors: features that have non-zero variance in all sources 
variable_features = function(variable.feature.source){
    vf.mask = logical(dim(variable.feature.source)[1])
    for (j in 1:length(vf.mask)){
        vf.mask[j] = all(variable.feature.source[j,])
    }
    return(rownames(variable.feature.source)[vf.mask])
}




# Z: an array of c(sources, features, samples). The raw source tensor 
# variable.feature.source: a matrix of booleans: number of feautures by number of sources. entries that are set to FALSE indicate no expression 
# fill_thr: numeric. the mean and variance of the noise to be filled
# seed.num: numeric. 
# return a modified version of Z: an array of c(sources, features, samples). 
# (source, feature) related entries in Z with no variation or no expression (as indicated by variable.feature.source)
# are filled with gaussian noise with fill_thr magnitude.

add_noise_Z = function(Z, variable.feature.source, fill_thr = 10**(-4), seed.num = 2023){
    set.seed(seed.num)
    
    for (h in 1:dim(Z)[1]){
        for (j in 1:dim(Z)[2]){
            if(!variable.feature.source[j,h]){
                # add gaussian noise
                Z[h,j,] = Z[h,j,] + abs(rnorm(dim(Z)[3], 0, fill_thr))
            }
        }
    }
    return(Z)
}

# Z: an array of c(sources, features, samples). The raw source tensor 
# n: number of samples to simulate 
# m: number of features to simulate 
# enrich.low: logical. set to T if we want to enrich the low entropy set 
# entropy.thr.ratio: numerical (0,1). only relavent if enrich.low is TRUE. set the threshold for what is considered as a low entropy feature. actual threshold is calculated as the ratio of the max possible entropy under that setup
# enrich.ratio: numerical (0,1) only relavent if enrich.low is TRUE. set the ratio of features out of m to be low entropy features
# max_sds: numerical. specify when learning the intial parameters (including entropy) in Z, what max_stds to use
# seed.num: numerical

# return Z
simulate_Z = function(Z, n, m, enrich.low, entropy.thr.ratio, enrich.ratio, max_sds = 3, seed.num = 2023){
    set.seed(seed.num)
    assert(m <= dim(Z)[2])
    assert(enrich.ratio >= 0 & enrich.ratio <= 1)
    assert(entropy.thr.ratio >= 0 & entropy.thr.ratio <= 1)
    
    params = calc_params_from_Z(Z, max_sds = 3)
    entropies = params$entropies
    k = dim(Z)[1]
    
    # low entropy enrich
    if (enrich.low){ 
        entropy.thr = as.numeric(entropy.thr.ratio * calc_entropy(array(eye(k), c(1,k,k)))) # needs to be adjusted by the number of celltype
        print(paste0("low entropy thr on ", k," celltypes: ", entropy.thr))

        low.entropy.genes.idx  = which(entropies <= entropy.thr)
        print(paste0("low entropy genes: ", length(low.entropy.genes.idx)))
        high.entropy.genes.idx  = which(entropies > entropy.thr)
        print(paste0("high entropy genes: ", length(high.entropy.genes.idx)))
        
        if (length(low.entropy.genes.idx) <= enrich.ratio * m){
            print(paste0("sampling all low entropy genes: ", length(low.entropy.genes.idx), " ",m - length(low.entropy.genes.idx)))
            
            feature.idx = c(low.entropy.genes.idx, sample(high.entropy.genes.idx, m - length(low.entropy.genes.idx), replace = F))
        }else{
            print(paste0("sampling ", enrich.ratio * 100, "% genes to be low entropy, rest to be high entropy: ",round(enrich.ratio * m), " ", m-round(enrich.ratio * m)))
            feature.idx = c(sample(low.entropy.genes.idx, round(enrich.ratio * m), replace = F), sample(high.entropy.genes.idx, (m - round(enrich.ratio * m)), replace = F))
        }
    }else{
        print(paste0("sampling ", m ," random genes no enrichment on low entropy"))
        feature.idx = sample(1: dim(Z)[2], m, replace = F)
    }

    Z = Z[,feature.idx,]
    
    
    
    # sample individuals in Z
    if (n > dim(Z)[3]){
        print("n is larger than the number of individuals in Z")
        print("draw with replacement")
        Z = Z[,,sample(1:dim(Z)[3], n, replace = T)]
    }else{
        print("n is smaller than the number of individuals in Z")
        print("draw without replacement")
        Z = Z[,,sample(1:dim(Z)[3], n, replace = F)]
    }
    dimnames(Z)[[3]] = as.vector(vapply(1:n, function(i) paste0("sample.",i), character(1)))
    
    return(Z)
}


#sim.data: a list that should include at least key "X" and "Z"
#scale.max_sds: numeric. scale factor is calculated on the subset of data that is within scale.max_stds of the mean in "X"
#scale.factor.th: numeric. the minimal bound on the scale.factor. mostly used for numerical stability, preventing almost 0 scale factor
#max_sds: numeric, specify when learning the new parameters in the scaled version, what max_stds to use

#This function takes the sim.data. per feature, calculated the scale.factor based on samples in side the scale.max_stds
#for those features with too small variance, scale.factor.th is set to the corresponding entry in the scale factor
#scale.factor is saved as a number of features by 1 matrix in side the sim.data list with corresponding rownames set to feature name
#generate X.scale, Z.scale by dividing each feature with scale.factor's corresponding entry
#relearn the parameter in this scaled space: params.scale with good portion of the data within "max_stds"
#return a list with the following key; "X","Z","feature.scale.factor", "X.scale", "Z.scale", "params.scale"
normalize_data = function (sim.data, scale.max_sds = 3, scale.factor.thr = 10**(-4), max_sds = 3){ 
    
    scale.factor = matrix(0, nrow(sim.data$X), 1)
    rownames(scale.factor) = rownames(sim.data$X)
    colnames(scale.factor) = "scale.factor"

    for (j in 1:nrow(sim.data$X)){
        sd_j = sd(sim.data$X[j,])
        
        if (sd_j < scale.factor.thr){
            #need a check here because if sd_j = 0 and scale.max_sds is INF, then scale.max_sds*sd_j is NA
            scale.factor[j,] = scale.factor.thr
        }else{
            indices_j <- abs(sim.data$X[j,] - mean(sim.data$X[j,])) <= scale.max_sds*sd_j
            if (sum(indices_j) <= 1){
                scale.factor[j,] = 0
            }else{
                scale.factor[j,] = sd(sim.data$X[j,indices_j])
            }
        }
    }

    if(any(scale.factor < scale.factor.thr)){
        print("certain feature has extremely close to 0 varinace in X")
        scale.factor[scale.factor < scale.factor.thr,] = scale.factor.thr 
    } 
    
    sim.data$feature.scale.factor = scale.factor
    sim.data$X.scale = sim.data$X/repmat(scale.factor, 1, ncol(sim.data$X))
    sim.data$Z.scale = copy(sim.data$Z)
    for (h in 1:dim(sim.data$Z)[1]){
        sim.data$Z.scale[h,,] = sim.data$Z[h,,]/repmat(scale.factor, 1, dim(sim.data$Z)[3])
    }
    sim.data$params.scale = calc_params_from_Z(sim.data$Z.scale, max_sds = max_sds)
    return(sim.data)
}
                                        
# Z.src: an array of c(sources, features, samples). The raw source tensor. order of the sources will be linked to correpsonding W.src. so first source in Z.src will be mixed with first source in W.src. source names dont have to be meaningful # W.src: a matrix with rows as samples and columns as source names. source names needs to be meaningful
# k: numerical. Number of sources
# m: numerical. Number of features
# n: numerical. Number of samples
# fit.dirichlet: logical.  if T, it will fit a dirichlet distribution and draw new W from it. 
# add.noise: logical.  if T, it will add small Gaussian noise abs(N(0, fill_thr)) to feature-source that show no meaningful expression
# variable_thr: numerical, minimal variance and mean across samples within max_sds required to be considered as a variable feature-source. 
# fill_thr: numerical, variance of the small Gaussian noise
# expression.qtl: numeric between (0, 1). the top quantiles deemed as highly expressed
# enrich.low: logical. set to T if we want to enrich the low entropy set 
# entropy.thr.ratio: numerical (0,1). only relavent if enrich.low is TRUE. set the threshold for what is considered as a low entropy feature. actual threshold is calculated as the ratio of the max possible entropy under that setup
# enrich.ratio: numerical (0,1) only relavent if enrich.low is TRUE. set the ratio of features out of m to be low entropy features
# max_sds: numerical. specify when learning the intial parameters (including entropy) in Z, what max_stds to use
# scale.max_sds: numerical. specify when scale to unit variance, what max_stds to use (samples inside will have unit variance)
# scale.factor.thr: numerical. specify the minimal scaling factor. this is mostly for numerical stability, set to 10 ** -4 by default
# seed.num: numeric. 
simulate_expression_mixture = function(Z.src, W.src, k, m, n, fit.dirichlet, 
                                       add.noise, variable_thr, fill_thr, expression.qtl, 
                                       enrich.low, entropy.thr.ratio, enrich.ratio, 
                                       max_sds, scale.max_sds, scale.factor.thr, seed.num = 2023){
    
    set.seed(seed.num)
    
    print("rename feature.ids to replace - or / with . ")
    feature.ids = dimnames(Z.src)[[2]]
    feature.ids = vapply(1: length(feature.ids), function(j) gsub("\\-", ".", feature.ids[j]), character(1))
    feature.ids = vapply(1: length(feature.ids), function(j) gsub("\\/", ".", feature.ids[j]), character(1))

    print("simulate W")
    W = simulate_W(W.src, n, k, fit.dirichlet, seed.num = seed.num)
    source.ids = colnames(W)
    sample.ids = rownames(W)
    
    print("extract sources selected in W, only kept those in Z")
    source.idx = vapply(1:k, function(h) which(colnames(W.src) == colnames(W)[h]), numeric(1))
    Z = Z.src[source.idx,,]
    # forces the sourc.ids to be the celltypes from Z (W is always estimated from blood)
    source.ids = dimnames(Z)[[1]]
    colnames(W)= source.ids
    # update Z's feature names 
    dimnames(Z)[[2]] = feature.ids
                              
                 
    print("simulate Z")
    variable.feature.source = calc_variable_feature_source(Z = Z, variable_thr = variable_thr, max_sds = max_sds)
    if(!add.noise){
        Z = Z[, variable_features(Z, variable.feature.source), ]
    }  
                        
    hef =  highly_expressed_features(Z, expression.qtl)                   
    Z = simulate_Z(Z[,hef,], n, m, enrich.low, entropy.thr.ratio, enrich.ratio, max_sds, seed.num)
    
    if(add.noise){
        Z = add_noise_Z(Z, variable.feature.source, fill_thr, seed.num)
    }
    #update feature.ids to the subset
    feature.ids = dimnames(Z)[[2]]                   
                        
                        
    print("generate X")
    X = matrix(0, m, n)
    rownames(X) = feature.ids
    colnames(X) = sample.ids   
    for (h in 1:k){
        X = X + Z[h,,] * repmat(t(W[,h]), m ,1)
    }
                        
    print("recording parameters used for reproducibility")
    entropy.thr = as.numeric(entropy.thr.ratio * calc_entropy(array(eye(k), c(1,k,k)))) # needs to be adjusted by the number of celltype   
    data.gen.params = list(fit.dirichlet = fit.dirichlet, add.noise = add.noise, fill_thr = fill_thr, 
                           expression.qtl = expression.qtl, enrich.low  = enrich.low, 
                           entropy.thr = entropy.thr, entropy.thr.ratio = entropy.thr.ratio, enrich.ratio = enrich.ratio, max_sds = max_sds,
                           seed.num = seed.num)
   
    print("put everything together")
    sim.data = list(source.ids = source.ids, 
                    feature.ids= feature.ids,
                    sample.ids = sample.ids,
                    W = W,
                    Z = Z,
                    X = X,
                    params = calc_params_from_Z(Z, max_sds), 
                    variable.feature.source = if(add.noise) variable.feature.source[feature.ids, ] else NULL,
                    data.gen.params = data.gen.params)
               
    
    print("normalization")
    sim.data = normalize_data(sim.data, scale.max_sds, scale.factor.thr, max_sds)
    return(sim.data)                  
}