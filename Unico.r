library(pbapply)
library(config)
library(pracma)
library(testit)
library(data.table)
library(matrixStats)
library(matrixcalc)

library(mgcv)
library(nloptr)

source("utils.r")


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


# related to mean updates
format_mean_param_j = function(W, C1, C2, mus_j, gammas_j, betas_j){
    
    k  = ncol(W)
    p1 = ncol(C1)
    p2 = ncol(C2)

    #https://stackoverflow.com/questions/15143829/repeat-the-values-of-a-row-in-matrix-n-number-of-times
    cols = rep(1:ncol(W), rep(p1,k))
               #source patter, repeat pattern of each element in source
    C1_tilde   = W[,cols] * repmat(C1,1,k)
    design_mat = cbind(W,C1_tilde,C2)
    #https://stackoverflow.com/questions/3823211/convert-a-matrix-to-a-1-dimensional-array
    # row1, row2, row3,...
    variables  = c(mus_j, gammas_j, betas_j)
    
    return(list("design_mat" = design_mat, "variables" = variables))
}

#design_mat = design_mat_j; variables = variables_j;
mean_updates_j = function (X_j, W, C1, C2, design_mat, variables, Us, config){
  
    k  = dim(W)[2]
    p1 = dim(C1)[2]
    p2 = dim(C2)[2]
    # penalty should be just a number 
    # penalty should be proportional to weighted total sample size and average/median value
    mean_penalty = config[["mean_penalty"]] * sum((Us)**2) * median(X_j)
    
    G <- list(y  = X_j*Us,
              X  = as.matrix(design_mat)*repmat(Us,1,ncol(design_mat)),    # individual by coeffs
              w  = numeric(length(X_j))+1,            # weight of each individuals
              S  = list(diag(k + p1 * k + p2)),     # a list of just 1 penalty matrix (identity)
              sp = mean_penalty,             # lambda: penalize weight/parameter (no penalty on the variable)
              p  = variables,     # a feasible starting point 
              C   = matrix(0,0,0),   # no equality constraint 
              Ain = rbind(diag(k + p1 * k + p2), -diag(k + p1 * k + p2)),              # inequality constrain mat 
              bin = c(rep(0,k),rep(config[["effect_size_negative_lb"]], p1 * k + p2), rep(-config[["max_mean"]],k+p1*k+p2)), # inequality constraint vec
              off = 0)

    res <- pcls(G);
    res[c(res[1:k]<config[["min_mean"]],rep(FALSE,p1*k+p2))] <- config[["min_mean"]]
    return(res)
}

get_X_j_ept = function (mus_j, W, C1, gammas_j, C2, betas_j){
    res = format_mean_param_j(W, C1, C2, mus_j, gammas_j, betas_j)
    X_j_ept = res$"design_mat" %*% res$"variables"
    return(X_j_ept)
}


# related to variance update
vec_2_tril = function(v,k){
    assert("wrong length or wrong dimension",length(v) == (k+1) * k /2)
    tril = matrix(0, nrow = k, ncol = k)
    tril[lower.tri(tril, diag = TRUE)] = v  
    return(tril)
}

tril_2_vec = function(tril){
    v = tril[which(lower.tri(tril, diag = TRUE), arr.ind = TRUE)]
    return(v)
}

get_WW_tensor = function(W){
    n = dim(W)[1]
    k = dim(W)[2]
    
    WW_tensor = array(rep(0,n * k * k), c(n, k, k))
    for (i in 1:n){
        WW_tensor[i,,] = W[i,] %*% t(as.matrix(W[i,]))
    }
    return (WW_tensor)
} 

get_var_j_sample = function (X_j, mus_j, W, C1, gammas_j, C2, betas_j){
    X_j_ept = get_X_j_ept (mus_j, W, C1, gammas_j, C2, betas_j)
    var_j_ept = (X_j - X_j_ept)**2
    return(var_j_ept)
}

# variables = c(Ls_j, taus_j)
var_objective = function(WW_tensor, var_j_sample, Vs, variables, config, return_grads = TRUE){
    
    Ls_j   = matrix(variables[1:(length(variables)-1)])
    taus_j = variables[length(variables)]
    
    n = dim(WW_tensor)[1]
    k = dim(WW_tensor)[2]
    
    # proportional to weighted sample size and proportional to variance size
    covar_penalty = config[["covar_penalty"]] * sum(Vs**2) * median(var_j_sample) 
    var_penalty   = config[["var_penalty"]] * sum(Vs**2) * median(var_j_sample)
    
    L_j_mat = vec_2_tril(Ls_j,k)
    sigmas_j = L_j_mat %*% t(L_j_mat) 
    var_j_ept = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_j) + taus_j, numeric(1)))
    
    ## TODO put this V outside
    objective = sum((Vs**2)*((var_j_sample - var_j_ept)**2)) 
    # add entire sigmas but then remove the extra penalty on diagnal
    objective = objective + sum(covar_penalty*(sigmas_j**2)) - sum((covar_penalty - var_penalty)*(diag(sigmas_j)**2))
    #print(c(sum((Vs**2)*((var_j_sample - var_j_ept)**2)), sum(adjusted_penalty*(Ls_j**2)), sum(adjusted_penalty*(Ls_j**2)) / objective))
 
    
#     grads = matrix(0,n, (k+1) * k /2)
#     for (i in 1:n){
#         grads[i,] = tril_2_vec(WW_tensor[i,,] %*% L_j_mat)
#     }
#     grads_coeff = Vs**2 * 2 * (var_j_sample - var_j_ept) * (-1) 
#     grads = repmat(grads_coeff, 1, (k+1) * k /2) * 2 * grads 
#     grads = colSums(grads)
#     grads = c(grads, sum(grads_coeff))
#     return(list("objective" = objective, "gradient"  = grads))
    
    # only calculate the gradient when using a derivative based optimizer
    if (return_grads){
        ## new implementation
        grad_coeff = Vs**2 * 2 * (var_j_sample - var_j_ept) * (-1)
        # braodcast
        grad_coeff_tensor = array(0,c(n,k,k))
        grad_coeff_tensor[1:n,,] = grad_coeff
        # multiply by 2 here for the gradients of Ls and collapse across individuals
        grads = colSums(2 * grad_coeff_tensor * WW_tensor, dim = 1)

        grads = tril_2_vec(grads %*% L_j_mat)
        # TODO: finish the gradiant calculation with penalty terms
        # grads = grads + ......
        # adding gradients of taus
        grads = c(grads, sum(grad_coeff))
        return(list("objective" = objective, 
                    "gradient"  = grads))
    }else{
        #cannot return a list, LN algorithm only want one number 
        #if provide a list with just $objective and no $gradient it will complain 
        return(objective)
    }                       
    
}

get_bounds = function(k, fit_tau = TRUE){
    # L bounds
    col_num = 1
    pos = 1
    idx = c()
    while (col_num <= k) {
        idx <- c(idx, pos)
        pos = pos + k - (col_num -1)
        col_num = col_num + 1
    }
    
    num_var = (k + 1) * k / 2
    lb <- rep(-Inf, num_var)
    ub <- rep(Inf,  num_var)
    lb[idx] = 0
    
    #tau bounds 
    if (fit_tau){
        lb = c(lb, 0)
        ub = c(ub, Inf) 
    }else{
        lb = c(lb, 0)
        ub = c(ub, 0) 
    }
    
    return(list("lb" = lb, "ub" = ub))
}

# get the opts for variance upate
get_opts = function(nloptr_opts_algorithm, scope_type, grad_type, config){    
    if (scope_type == "L"){
        opts = list("algorithm"= nloptr_opts_algorithm, 
                    "xtol_rel"=  config[["nloptr_local_xtol_rel"]], 
                    "maxeval" =  config[["nloptr_local_maxeval"]])
    
    }else if (scope_type == "G"){
        opts = list ("algorithm"= nloptr_opts_algorithm, 
                     "xtol_rel"=  config[["nloptr_global_xtol_rel"]], 
                     "maxeval" =  config[["nloptr_global_maxeval"]])  
    }else{
        return(NULL)
    }
    
    if(grad_type == "D"){
        opts = c(opts, list("check_derivatives" = config[["nloptr_grad_check_derivatives"]]))
    }
       
    return(opts)
}
# the core of the updates with parallel                              
#X_j = X[,j]; mus_j = mus_hat[j,]; gammas_j = gammas_hat[j,]; betas_j = betas_hat[j,]; Ls_j = Ls_hat[j,]; taus_j = taus_hat[j];
fit_mean_var = function(X_j, W, C1, C2, mus_j, gammas_j, betas_j, Ls_j, taus_j, mean_max_iterations, var_max_iterations, fit_tau, nloptr_opts_algorithm, config){
    
    k  = dim(W)[2]
    n  = dim(C1)[1]
    p1 = dim(C1)[2]
    p2 = dim(C2)[2]

    #mean updates
    #Us         = matrix(1,n,1)
    if (config[["init_weight"]] == "W_norm"){
        Us         = as.matrix(rowSums(W**2))
    }else if(config[["init_weight"]] == "inv_W_norm"){
        Us         = 1/as.matrix(rowSums(W**2)) 
    }else{
        Us         = matrix(1,n,1)
    }
    
    #masking outliers for mean updates
    max_stds = config[["max_stds"]]
    if(!is.null(max_stds)){
        sd_j <- sd(X_j)
        outlier_mask = abs(X_j - mean(X_j)) >= max_stds*sd_j
        #flog.debug(paste0(sum(outlier_mask)," outliers are masked for mean updates"))
        Us[outlier_mask] = 0
    }
    

    Us_array   = c()
    for (itr in 1:mean_max_iterations){
        
        # TODO; separate the format_j function, design mat is constant 
        format_j      = format_mean_param_j(W, C1, C2, mus_j, gammas_j, betas_j)
        design_mat_j  = as.matrix(format_j$"design_mat")
        variables_j   = format_j$"variables"
    
        res     = mean_updates_j(X_j, W, C1, C2, design_mat_j, variables_j, Us, config)
        mus_j   = res[1:k]
        gammas_j= res[seq(k+1, k + k * p1, length=k * p1)]
        betas_j = res[seq(k + k * p1 + 1, k + k*p1 + p2, length=p2)] 
        
        ##check cost with current Us but best means and effect size 
        #print("mean cost")
        #print(sum(Us**2 * (X_j - design_mat_j %*% res)**2))
        
        # update Us 
        Us_array = cbind(Us_array, Us)
        
        Us = abs(1/(X_j - design_mat_j %*% res))
        Us = as.matrix(replace(Us, which (Us > config[["max_u"]]), config[["max_u"]]))
        Us = as.matrix(Us)
        Us = Us/max(Us)
        
        # masking those outliers 
        if(!is.null(max_stds)){
            Us[outlier_mask] = 0
        }
        
    }

    #variance updates
    WW_tensor = get_WW_tensor(W) 
    #Vs         = matrix(1,n,1)
    if (config[["init_weight"]] == "W_norm"){
        Vs         = as.matrix(rowSums(W**2))
    }else if(config[["init_weight"]] == "inv_W_norm"){
        Vs         = 1/as.matrix(rowSums(W**2)) 
    }else{
        Vs         = matrix(1,n,1)
    }
        
    #masking outliers for variance updates
    max_stds = config[["max_stds"]]
      if(!is.null(max_stds)){
        sd_j <- sd(X_j)
        outlier_mask = abs(X_j - mean(X_j)) >= max_stds*sd_j
        #print(paste0(sum(outlier_mask)," outliers are masked for variance updates"))
        Vs[outlier_mask] = 0
    }
                           
    #prepare input to nlpotr
    bounds = get_bounds(k, fit_tau) 
    #get scope_type and grad_type 
    scope_type = strsplit(strsplit(nloptr_opts_algorithm, split = "_")[[1]][2], "")[[1]][1] 
    grad_type  = strsplit(strsplit(nloptr_opts_algorithm, split = "_")[[1]][2], "")[[1]][2] 
    return_grads = grad_type == "D"
    opts = get_opts(nloptr_opts_algorithm, scope_type, grad_type, config) 
    var_j_sample = get_var_j_sample(X_j, mus_j, W, C1, gammas_j, C2, betas_j)
    
    
    Vs_array     = list(var_max_iterations)
    Ls_array     = list(var_max_iterations)
    sigmas_array = list(var_max_iterations)
    
    for (itr in 1:var_max_iterations){
        res <- nloptr(x0 = c(Ls_j, taus_j),
                      eval_f = function(x, WW_tensor, var_j_sample, Vs, return_grads) var_objective(WW_tensor, var_j_sample, Vs, x, config, return_grads),
                      lb = bounds$"lb",
                      ub = bounds$"ub", 
                      WW_tensor = WW_tensor, 
                      var_j_sample = var_j_sample, 
                      Vs = Vs,
                      return_grads = return_grads,
                      opts=opts)
        
        # update variables 
        Ls_j  = res$solution[1:(length(res$solution)-1)]
        taus_j = res$solution[length(res$solution)]  
        sigmas_j = vec_2_tril(Ls_j,k) %*% t(vec_2_tril(Ls_j,k))
        
        ##print("varaince")
        ##print(var_objective(WW_tensor, var_j_sample, Vs, c(Ls_j, taus_j)))
        
        # track
        Vs_array[[itr]]     = Vs
        Ls_array[[itr]]     = Ls_j
        sigmas_array[[itr]] = as.vector(sigmas_j)
        
        # update Vs
        var_j_ept = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_j) + taus_j, numeric(1)))
        Vs = abs(1/(var_j_sample - var_j_ept))
        Vs = as.matrix(replace(Vs, which (Vs > config[["max_v"]]), config[["max_v"]]))
        Vs = Vs/max(Vs)
                                  
                                  
        # masking those outliers 
        if(!is.null(max_stds)){
            Vs[outlier_mask] = 0
        }
        
    }
                                  
    #return
    params = list("mus_j"    = mus_j,
                  "gammas_j" = gammas_j,
                  "betas_j"  = betas_j,
                  "Ls_j"     = Ls_j,
                  "taus_j"   = taus_j,
                  "Us"       = Us,
                  "Us_array" = Us_array,
                  "Vs"       = Vs,
                  "Vs_array" = Vs_array,
                  "Ls_array" = Ls_array,
                  "sigmas_array" = sigmas_array)

    
    return(params)
}
  
# Unico main function   
# input X should be the bulk level matrix of m features by n samples
# input W should be the celltype composition matrix of n samples by k celltypes
# input C1 should be the celltype specific covariates matrix of n sample by p1 features
# input C2 should be the global covariates matrix of n sample by p2 features     
# align_var will force the tensor to match the estimated variance. however, it will no longer respect the reconstructed X is the same as the real X      
Unico = function(X, W, C1, C2, 
                mean_max_iterations = 2, var_max_iterations = 3, fit_tau = FALSE,nloptr_opts_algorithm = "NLOPT_LN_COBYLA", 
                mean_penalty = 0, var_penalty = 0.01, covar_penalty = 0.01, 
                max_u = 1, max_v =1, init_weight = "default", max_stds = 2, align_var = FALSE, Ls.init = NULL, 
                parallel = TRUE, num_cores= 12, 
                config_file = NULL, log_file = "Unico.log", debug = FALSE, verbose = TRUE, scale.factor = NULL){
    
    
    if (is.null(config_file)){
        config = config::get(file = file.path("inst","extdata", "config.yml"))
    }else{
        config::get(file = config_file)
    } 
    
    start_logger(log_file, debug, verbose)
    flog.info("Starting Unico...")
    op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")
        
    sample.ids = colnames(W)
    feature.ids= rownames(X)
    source.ids = colnames(W)
    C1.ids = colnames(C1)
    C2.ids = colnames(C2)
   
    # here calculate the scale.factor and later save this information in mdl
    # scale and data to unit variance and move on, all downstream calculation like parameter estimatation are on the scaled data
    # only use scale.factor at the tensor stage to scale it back to original size
    # TODO: modified here. double check 
    if (is.null(scale.factor)){
        scale.max_stds = max_stds
        scale.min.thr = 10**(-4)
        
        scale.factor = matrix(0, nrow(X), 1)
        rownames(scale.factor) = feature.ids
        colnames(scale.factor) = c("scale.factor")

        for (j in 1:nrow(X)){
            sd_j = sd(X[j,])
            indices_j <- abs(X[j,] - mean(X[j,])) <= scale.max_stds*sd_j
            if (sum(indices_j) <= 1){
                scale.factor[j,] = 0
            }else{
                scale.factor[j,] = sd(X[j,indices_j])
            }
        }

        if(any(scale.factor < scale.min.thr)){
            print("certain feature has extremely close to 0 varinace in X")
            scale.factor[scale.factor < scale.min.thr,] = scale.min.thr
        }  
    }
    
    
    
    X = scale_feature_matrix(X, 1/scale.factor) 
    #make it samples by features
    X = t(X)
    #TODO: add validate input 
    if (is.null(C1)) C1 <- matrix(0, nrow=nrow(X), ncol=0)
    if (is.null(C2)) C2 <- matrix(0, nrow=nrow(X), ncol=0)
    
    n = dim(X)[1]
    m = dim(X)[2]
    k = dim(W)[2]
    p1 = dim(C1)[2]
    p2 = dim(C2)[2]
    

     # assumes
    assert("covar_penalty should be larger than var_penalty",covar_penalty >= var_penalty)
    config[["mean_penalty"]] <- mean_penalty
    config[["var_penalty"]] <- var_penalty
    config[["covar_penalty"]] <- covar_penalty
    config[["max_u"]] <- max_u 
    config[["max_v"]] <- max_v
    config[["init_weight"]] <- init_weight
    config[["max_stds"]]    <- max_stds
    config[["align_var"]] <- align_var

    #TODO: add initialization
    mus_hat    = matrix(1, m, k)
    gammas_hat = matrix(0, m, k * p1) # already flatten to [row(ct)1, row(ct)2, ...]    
    betas_hat  = matrix(0, m, p2)
    if (is.null(Ls.init)){
        Ls_hat     = matrix(1, m, ((k + 1) * k / 2) )
    }else{
        Ls_hat = Ls.init
    }
    taus_hat   = matrix(0, m)
    
    # tracking
    Vs_hat_list     = list(var_max_iterations)
    Ls_hat_list     = list(var_max_iterations)
    sigmas_hat_list = list(var_max_iterations)
    
    for (itr in 1:var_max_iterations){
        Vs_hat_list[[itr]]     = matrix(0, m, n)
        Ls_hat_list[[itr]]     = matrix(0, m, (k + 1) * k / 2)
        sigmas_hat_list[[itr]] = matrix(0, m, k*k)
    }
    
    flog.info("Starting parameter learning ...")
    
    cl <- if (parallel) init_cluster(num_cores) else NULL
    if (parallel) clusterExport(cl, c("X", "W", "C1", "C2", 
                                      "mus_hat", "gammas_hat", "betas_hat", "Ls_hat", "taus_hat",
                                      "mean_max_iterations", "var_max_iterations", "config", "fit_tau", "nloptr_opts_algorithm",
                                      "format_mean_param_j","mean_updates_j", "get_X_j_ept",
                                      "get_opts", "get_bounds", "var_objective", "get_var_j_sample", "get_WW_tensor", "tril_2_vec", "vec_2_tril", "fit_mean_var"), envir=environment())


    res <- pblapply(1:m,function(j) fit_mean_var(X[,j], W, C1, C2, 
                                                 mus_hat[j,], gammas_hat[j,], betas_hat[j,], Ls_hat[j,],taus_hat[j], 
                                                 mean_max_iterations, var_max_iterations, fit_tau, nloptr_opts_algorithm, config), 
                    cl = cl)

    if (parallel) stop_cluster(cl)

    #update parameters
    for (j in 1:m){
        mus_hat[j,]    = res[[j]][["mus_j"]]
        gammas_hat[j,] = res[[j]][["gammas_j"]]
        betas_hat[j,]  = res[[j]][["betas_j"]]
        Ls_hat[j,]     = res[[j]][["Ls_j"]]
        taus_hat[j]    = res[[j]][["taus_j"]]
        for (itr in 1:var_max_iterations){
            Vs_hat_list[[itr]][j,] = res[[j]][["Vs_array"]][[itr]]
            Ls_hat_list[[itr]][j,] = res[[j]][["Ls_array"]][[itr]]
            sigmas_hat_list[[itr]][j,] = res[[j]][["sigmas_array"]][[itr]]
            
        }
        
    }

    flog.info("Formating results ...")
    
    Unico.mdl = c()
    Unico.mdl$W = W
    Unico.mdl$C1 = C1
    Unico.mdl$C2 = C2
    
    rownames(mus_hat) = feature.ids
    colnames(mus_hat) = source.ids
    Unico.mdl$mus_hat  = mus_hat
    
    C1.ids_ <- matrix("0",k*p1,1)
    if (p1){
        for (j in 1:k){
            C1.ids_[((j-1)*p1+1):(p1*j)] <- unlist(lapply(C1.ids,function(i) paste(source.ids[j],".",i,sep="")))
        }
    }
    
    rownames(gammas_hat) = feature.ids
    colnames(gammas_hat) = C1.ids_
    Unico.mdl$gammas_hat = gammas_hat

    rownames(betas_hat) = feature.ids
    colnames(betas_hat) = C2.ids
    Unico.mdl$betas_hat  = betas_hat 

    Unico.mdl$Ls_hat     = Ls_hat
    
    Unico.mdl$sigmas_hat = array(0, c(m,k,k))
    for (j in 1:m){
        Unico.mdl$sigmas_hat[j,,] = vec_2_tril(Ls_hat[j,], k) %*%  t(vec_2_tril(Ls_hat[j,], k))
    }
    dimnames(Unico.mdl$sigmas_hat)[[1]] = feature.ids
    dimnames(Unico.mdl$sigmas_hat)[[2]] = source.ids
    dimnames(Unico.mdl$sigmas_hat)[[3]] = source.ids
    
    rownames(taus_hat) = feature.ids
    colnames(taus_hat) = c("taus_hat")
    Unico.mdl$taus_hat   = taus_hat

    Unico.mdl$Vs_hat_list = Vs_hat_list 
    Unico.mdl$Ls_hat_list = Ls_hat_list 
    Unico.mdl$sigmas_hat_list = sigmas_hat_list
    Unico.mdl$config   = config
    Unico.mdl$scale.factor   = scale.factor
                                                          
    flog.info("Finished parameter learning")
    return(Unico.mdl)
}                            
   

# related to estimate Z
#X_j=X[,j]; mus_j=Unico.mdl$"mus_hat"[j,];gammas_j=Unico.mdl$"gammas_hat"[j,];betas_j=Unico.mdl$"betas_hat"[j,]; sigmas_hat_j=Unico.mdl$"sigmas_hat"[j,,]; taus_j =Unico.mdl$"taus_hat"[j]
get_Z_j = function(X_j, W, C1, C2, mus_j, gammas_j, betas_j, sigmas_hat_j, taus_j, config){
    
    n = dim(W)[1]
    k = dim(W)[2]
    
    ## TODO provide WW_tensor as argument 
    WW_tensor = get_WW_tensor(W)
    X_j_ept = get_X_j_ept (mus_j, W, C1, gammas_j, C2, betas_j)
 
    var_j_ept    = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_hat_j) + taus_j, numeric(1)))
    #if no variance in any entry (var_j_ept is n by 1)
    #(rare but consider celltype1 is estimated to have 0 variance and sample's only has celltype1 and no other celltypes)       
    #setting pure 0 variance people to have pseudo variance 1. exact number does not matter 
    #the final Z_j takes Z_ept         
    if (sum(var_j_ept <= 10**(-4)) !=0){var_j_ept[var_j_ept <= 10**(-4)] = 1}
    
    covar_z_x = W %*% sigmas_hat_j
    var_j_ept_inv = -(1/var_j_ept)
    B = repmat(var_j_ept_inv, 1, k) * covar_z_x
    Z_ept = t(repmat(as.matrix(mus_j), 1,n)) + C1  %*% t(matrix(gammas_j, nrow = dim(W)[2], byrow = TRUE))
    X_diff = X_j_ept - X_j
    delta = B * repmat(X_diff, 1, k)
                                 
    if(config[["align_var"]]){
                                    
        #set the ratio between real sigmas_hat and empirical covaraince from the perturbation as scale.factor
        var.scale.factor = sqrt(diag(sigmas_hat_j)/diag(cov(delta)))
        Z_j = Z_ept  +  delta * repmat(var.scale.factor, n, 1) 
    }
    else{
        Z_j = Z_ept  +  delta
    }

    return(Z_j)
}
                                                     
    
# # input Z should be a array of size (m,n,k)
# # output Z as a new array of size (k,m,n)                  
# format_Z = function (Z){
#     m = dim(Z)[1]
#     n = dim(Z)[2]
#     k = dim(Z)[3]
#     # format Z from (m,n,k) to (k,m,n) 
#     Z_shifted = array(0, c(k,m,n))
#     for (h in 1:k){
#         Z_shifted[h,,] = Z[,,h]
#     }
#     return(Z_shifted)
# }

# tensor main function
# Unico.mdl = Unico.mdl
# input X should be the bulk level matrix of m features by n samples
# input W should be the celltype composition matrix of n samples by k celltypes
# input C1 should be the celltype specific covariates matrix of n sample by p1 features
# input C2 should be the global covariates matrix of n sample by p2 features                                 
tensor = function(X, W, C1, C2, Unico.mdl, parallel = TRUE, num_cores = 12, 
                  config_file = NULL, log_file = "Unico.log", debug = FALSE, verbose = TRUE){
    if (is.null(config_file)){
        config = config::get(file = file.path("inst","extdata", "config.yml"))
    }else{
        config::get(file = config_file)
    } 
    start_logger(log_file, debug, verbose)
    flog.info("Starting tensor ...")
    op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")
        
    #transfer it to scaled version
    X = X/repmat(Unico.mdl$scale.factor[rownames(X), , drop = F], 1, ncol(X))

    source.ids  = colnames(W)
    feature.ids = rownames(X)
    sample.ids  = colnames(X)
    
    X = t(X)

    if (is.null(C1)) C1 <- matrix(0, nrow=nrow(X), ncol=0)
    if (is.null(C2)) C2 <- matrix(0, nrow=nrow(X), ncol=0)

    n = dim(X)[1]
    m = dim(X)[2]
    k = dim(W)[2]
    p1 = dim(C1)[2]
    p2 = dim(C2)[2]
    
    Z_hat   = array(rep(0, m * n * k), c(m,n,k))
    # parallel
    cl <- if (parallel) init_cluster(num_cores) else NULL
    if (parallel) clusterExport(cl, c("X","W","C1", "C2", "Unico.mdl", 
                                      "format_mean_param_j","get_X_j_ept", "get_var_j_sample", "get_WW_tensor", 
                                      "tril_2_vec", "vec_2_tril",
                                      "get_Z_j"), envir=environment())
    res <- pblapply(1:m,function(j) get_Z_j(X[,j], W, C1, C2, 
                                            Unico.mdl$"mus_hat"[j,], Unico.mdl$"gammas_hat"[j,], Unico.mdl$"betas_hat"[j,], 
                                            Unico.mdl$"sigmas_hat"[j,,],  Unico.mdl$"taus_hat"[j], Unico.mdl$"config"), 
                    cl = cl)
    if (parallel) stop_cluster(cl)

    flog.info("Formating tensor result...")
    
    #update 
    for (j in 1:m){
        Z_hat[j,,] = res[[j]]
    }
    Z_hat = aperm(Z_hat, c(3,1,2))
    dimnames(Z_hat)[[1]] = source.ids
    dimnames(Z_hat)[[2]] = feature.ids
    dimnames(Z_hat)[[3]] = sample.ids
    
    # scale the data back to its original scale
    for (h in 1:k){
        # updated the scale.factor
        Z_hat[h,,] = Z_hat[h,,] * repmat(Unico.mdl$scale.factor[feature.ids, , drop = F], 1, length(sample.ids))
    }
    flog.info("Finished tensor estimation")
    return(Z_hat)
}
      
# Association 
                                 
calc_C1_W_interactions <- function(W,C1){
    n <- nrow(W)
    k <- ncol(W)
    p1 <- ncol(C1)
    if (p1){
        return( hadamard.prod(Reshape(Reshape(apply(W, 2, function(v) repmat(v,1,p1)), n*p1*k,1), n,p1*k), repmat(C1, 1, k)))
    }else{
        return(matrix(0,n,0))
    }
}
                                                

mask_outliers = function(x, std.thr, min.thr = NULL, max.thr = NULL, filter.direction = "both"){
    if(is.null(min.thr)){min.thr = min(x)}
    if(is.null(max.thr)){max.thr = max(x)}
    
    mu = mean(x)
    sigmas = std(x)
    min.std.thr = mu - std.thr * sigmas 
    max.std.thr = mu + std.thr * sigmas 
    
    if (filter.direction == "both"){
        mask = (x <= max.thr) &  (x >= min.thr) & (x <= max.std.thr) & (x >= min.std.thr) 
    }else if (filter.direction == "greater"){
        mask = (x <= max.thr) &  (x >= min.thr) & (x <= max.std.thr) 
    }else{
        mask = (x <= max.thr) &  (x >= min.thr) & (x >= min.std.thr) 
    }
    
    return(mask)
}

add_C1_C2_pvals_parametric = function(X, Unico.mdl, slot_name = "parametric", 
                                      diag_only = FALSE, intercept = TRUE, 
                                      X_max_stds  = 2, 
                                      Q_max_stds  = Inf, 
                                      XQ_max_stds = Inf,
                                      parallel = FALSE, num_cores = NULL, 
                                      config_file = NULL, log_file = "Unico.log", debug = FALSE, verbose = TRUE){
    if (is.null(config_file)){
        config = config::get(file = file.path("inst","extdata", "config.yml"))
    }else{
        config::get(file = config_file)
    } 
    start_logger(log_file, debug, verbose)
    
    #input format
    X          = X/Unico.mdl$scale.factor[rownames(X), ] #scale features so that ept calculation make sense 
    X          = t(X) #X should be samples by features
    
    W          = Unico.mdl$W
    C1         = Unico.mdl$C1
    C2         = Unico.mdl$C2
  
    n = nrow(X)
    m = ncol(X)
    k = ncol(W)
    p1 = ncol(C1)
    p2 = ncol(C2)

    source.ids  = colnames(W)
    feature.ids = colnames(X)
    sample.ids  = rownames(X)
    C1.ids      = colnames(C1)
    C2.ids      = colnames(C2)

    taus_hat   = Unico.mdl$taus_hat[feature.ids, ,drop = F]
    sigmas_hat = Unico.mdl$sigmas_hat[feature.ids, , ]

    flog.info("Preparing weights for parametric pvals calculation ...")
    op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")
      
    # ct1 related interations, ct2 related related interations, ....
    C1.ids_ <- matrix("0",k*p1,1)
    if (p1){
        for (j in 1:k){
            C1.ids_[((j-1)*p1+1):(p1*j)] <- unlist(lapply(C1.ids,function(i) paste(source.ids[j],".",i,sep="")))
        }
    }
      
    C1_ <- calc_C1_W_interactions(W,C1)
    rownames(C1_) = sample.ids
    colnames(C1_) = C1.ids_
      
    s_len = k + p1*k + p2 
      
    # result matrix
    phi_hat = matrix(1, m, s_len)
    rownames(phi_hat) = feature.ids
    colnames(phi_hat) = c(source.ids, C1.ids_, C2.ids)
      
    phi_se = matrix(1, m, s_len)
    rownames(phi_se) = feature.ids
    colnames(phi_se) = c(source.ids, C1.ids_, C2.ids)

    phi_pvals = matrix(1, m, s_len)
    rownames(phi_pvals) = feature.ids
    colnames(phi_pvals) = c(source.ids, C1.ids_, C2.ids)

    gammas_hat_pvals.joint = matrix(1, m, p1)
    rownames(gammas_hat_pvals.joint) = feature.ids
    colnames(gammas_hat_pvals.joint) = c(C1.ids)

    masks = matrix(FALSE, n, m)
    rownames(masks) = sample.ids
    colnames(masks) = feature.ids

    # weighting scheme
    W_norms = matrix(0, n, m) 
    rownames(W_norms) = sample.ids
    colnames(W_norms) = feature.ids
      
    WW_tensor = get_WW_tensor(W)
    for(j in 1:m){
        if(diag_only){
            # just use variance
            sigmas_hat_j = diag(diag(sigmas_hat[j,,] ))
        }else{
            sigmas_hat_j = sigmas_hat[j,,] 
        }
        W_norms[,j] =  vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_hat_j) + taus_hat[j,], numeric(1))
    }
  
  
    # expected variance used for the inverse variance weighting scheme
    Q = 1/W_norms


    flog.info(paste0("Starting parametric pvals calculation: ", slot_name))
    cl <- NULL
    if (parallel){
        cl <- init_cluster(num_cores)
        clusterExport(cl, varlist = c("X","W","C1_","C2","Q",
                                      "k","p1","p2",
                                      "lm",
                                      "X_max_stds","Q_max_stds","mask_outliers","intercept"), envir=environment())
    }
  
    res <- pblapply(1:m,function(j) {
        X_j = as.matrix(X[,j, drop = F])
        Q_j = as.matrix(Q[,j, drop = F])
        df = data.frame(y = as.vector(X_j), cbind(W, C1_, C2))
        df = df * repmat(Q_j**0.5, 1 ,ncol(df)) #take the sqaure root as the OLS has a square in objective function

        ####### masking #######
        #filter 1: if X_j is outlier 
        mask1 = mask_outliers(as.vector(X_j),  std.thr = X_max_stds)
        #filter 2: if Q_j is outlier 
        mask2 = mask_outliers(as.vector(Q_j),  std.thr = Q_max_stds)
        #filter 3: if XQ_j is outlier 
        mask3 = mask_outliers(as.vector(df$y), std.thr = XQ_max_stds)
        
        mask =  mask1 & mask2 & mask3
        df = df[mask, ]
        
        ####### fitting #######
        if (intercept){
            mdl1.fit <- lm(y ~ ., data = df)
        }else{
            mdl1.fit <- lm(y ~ 0 + ., data = df)
        }

        mdl1.coef <- summary(mdl1.fit)$coefficients
        mdl1.cov.names <- colnames(df)[colnames(df) != 'y']
        
        phi_hat_j   <- mdl1.coef[mdl1.cov.names,"Estimate"] 
        phi_se_j    <- mdl1.coef[mdl1.cov.names,"Std. Error"]  
        phi_pvals_j <- mdl1.coef[mdl1.cov.names,"Pr(>|t|)"] 
        
        gammas_hat_pvals_joint_j <- numeric(p1)+1
        if (p1){
            #johnsonc                       
            for (d in 1:p1){
                C1_null = C1_[,setdiff(1:(p1*k),seq(d,k*p1,p1))]
                df = data.frame(y = as.vector(X_j), cbind(W, C1_null, C2))
                df = df * repmat(Q_j**0.5, 1 ,ncol(df))
                df = df[mask, ]
        
                if (intercept){
                    mdl0.fit <- lm(y ~ ., data = df)
                }else{
                    mdl0.fit <- lm(y ~ 0 + ., data = df)
                }
                
                anova.fit <- anova(mdl0.fit, mdl1.fit)
                gammas_hat_pvals_joint_j[d] <- anova.fit$"Pr(>F)"[2]
            }
        }
        return(c(phi_hat_j, phi_se_j, phi_pvals_j, mask, gammas_hat_pvals_joint_j));
        }, cl = cl )

    if (parallel) stop_cluster(cl)
    for (j in 1:m){
        assert(length(res[[j]]) == 3* s_len + n + p1 )
        phi_hat[j,]   <- res[[j]][(0*s_len + 1):(1*s_len)]
        phi_se[j,]    <- res[[j]][(1*s_len + 1):(2*s_len)]
        phi_pvals[j,] <- res[[j]][(2*s_len + 1):(3*s_len)]
        masks[,j]     <- res[[j]][(3*s_len + 1):(3*s_len + n)]
        gammas_hat_pvals.joint[j, ]     <- res[[j]][(3*s_len + n + 1):(3*s_len + n + p1)]
    }
  
    param.res = list(Q      = t(Q),     # now features by samples
                     masks  = t(masks), # now features by samples
                   
                     phi_hat   = phi_hat, 
                     phi_se    = phi_se,
                     phi_hat_pvals = phi_pvals,
                   
                     gammas_hat       = phi_hat[, (k+1):(k + p1*k)],
                     betas_hat        = phi_hat[, (k + p1*k+1):(k + p1*k + p2)],
                     gammas_hat_pvals = phi_pvals[, (k+1):(k + p1*k)],
                     betas_hat_pvals  = phi_pvals[, (k + p1*k+1):(k + p1*k + p2)], 
                   
                     gammas_hat_pvals.joint = gammas_hat_pvals.joint)

    Unico.mdl[[slot_name]] = param.res
    flog.info(paste0("Finished parametric pvals calculation: ", slot_name))
    return(Unico.mdl)
}
         
add_C1_C2_pvals_asymptotic = function(X, Unico.mdl, slot_name = "asymptotic", 
                                      diag_only = FALSE, intercept = TRUE,
                                      X_max_stds  = 2, #bulk
                                      Q_max_stds  = Inf, #expect variance 
                                      V_min_qlt   = 0.05,   #moment condition variance 
                                      parallel = FALSE, num_cores = NULL, 
                                      config_file = NULL, log_file = "Unico.log", debug = FALSE, verbose = TRUE){

    if (is.null(config_file)){
        config = config::get(file = file.path("inst","extdata", "config.yml"))
    }else{
        config::get(file = config_file)
    } 
    start_logger(log_file, debug, verbose)
    
    #input format
    X          = X/Unico.mdl$scale.factor[rownames(X), ] #scale features so that ept calculation make sense 
    X          = t(X) #X should be samples by features
    
    W          = Unico.mdl$W
    C1         = Unico.mdl$C1
    C2         = Unico.mdl$C2
  
    n = nrow(X)
    m = ncol(X)
    k = ncol(W)
    p1 = ncol(C1)
    p2 = ncol(C2)

    source.ids  = colnames(W)
    feature.ids = colnames(X)
    sample.ids  = rownames(X)
    C1.ids      = colnames(C1)
    C2.ids      = colnames(C2)

    mus_hat     = Unico.mdl$mus_hat[feature.ids, ]
    gammas_hat  = Unico.mdl$gammas_hat[feature.ids,]
    betas_hat   = Unico.mdl$betas_hat[feature.ids,]

    taus_hat   = Unico.mdl$taus_hat[feature.ids, ,drop = F]
    sigmas_hat = Unico.mdl$sigmas_hat[feature.ids, , ]

    flog.info("Preparing weights for asymptotic pvals calculation ...")
    op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")
      
    # ct1 related interations, ct2 related related interations, ....
    C1.ids_ <- matrix("0",k*p1,1)
    if (p1){
        for (j in 1:k){
            C1.ids_[((j-1)*p1+1):(p1*j)] <- unlist(lapply(C1.ids,function(i) paste(source.ids[j],".",i,sep="")))
        }
    }
      
    C1_ <- calc_C1_W_interactions(W,C1)
    rownames(C1_) = sample.ids
    colnames(C1_) = C1.ids_
    
    S = cbind(W, C1_, C2)                                                       
    s_len = k + p1*k + p2 
      
    # result matrix
    phi_hat = matrix(1, m, s_len)
    rownames(phi_hat) = feature.ids
    colnames(phi_hat) = c(source.ids, C1.ids_, C2.ids)
      
    phi_se = matrix(1, m, s_len)
    rownames(phi_se) = feature.ids
    colnames(phi_se) = c(source.ids, C1.ids_, C2.ids)

    phi_pvals = matrix(1, m, s_len)
    rownames(phi_pvals) = feature.ids
    colnames(phi_pvals) = c(source.ids, C1.ids_, C2.ids)
    
    #estimated moment condition variance                                                   
    fphi_var = matrix(-1, n, m)
    rownames(fphi_var) = sample.ids
    colnames(fphi_var) = feature.ids
                                                          
    masks = matrix(FALSE, n, m)
    rownames(masks) = sample.ids
    colnames(masks) = feature.ids

    # weighting scheme
    W_norms = matrix(0, n, m) 
    rownames(W_norms) = sample.ids
    colnames(W_norms) = feature.ids
      
    WW_tensor = get_WW_tensor(W)
    for(j in 1:m){
        if(diag_only){
            # just use variance
            sigmas_hat_j = diag(diag(sigmas_hat[j,,] ))
        }else{
            sigmas_hat_j = sigmas_hat[j,,] 
        }
        W_norms[,j] =  vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_hat_j) + taus_hat[j,], numeric(1))
    }
  
  
    # expected variance used for the inverse variance weighting scheme
    Q = 1/W_norms

    flog.info(paste0("Starting asymptotic pvals calculation: ", slot_name))
    cl <- NULL
    if (parallel){
        cl <- init_cluster(num_cores)
        clusterExport(cl, varlist = c("X","W","C1_","C2","Q","S",
                                      "k","p1","p2","s_len",
                                      "mus_hat","gammas_hat","betas_hat",
                                      "lm", "intercept",
                                      "get_X_j_ept", "format_mean_param_j",
                                      "X_max_stds","Q_max_stds","V_min_qlt","mask_outliers"), envir=environment())
    }
  
    res <- pblapply(1:m,function(j) {
        X_j = as.matrix(X[,j, drop = F])
        Q_j = as.matrix(Q[,j, drop = F])
        X_j_ept    = get_X_j_ept (mus_hat[j,], W, C1, gammas_hat[j,], C2, betas_hat[j,])
        fphi_var_j = (X_j - X_j_ept)**2 
        
        ####### masking #######
        #filter 1: if X_j is outlier 
        mask1 = mask_outliers(as.vector(X_j),  std.thr = X_max_stds)
        #filter 2: if Q_j is outlier 
        mask2 = mask_outliers(as.vector(Q_j),  std.thr = Q_max_stds)
        #filter 3: if V_j is too small. inverse will target those with super small value which suggest overfitting
        mask3 = as.vector(fphi_var_j) >= quantile(fphi_var_j, V_min_qlt)
        
        mask =  mask1 & mask2 & mask3
        S    =   S[mask, ]
        X_j  = X_j[mask]
        Q_j  = Q_j[mask]
        fphi_var_j = fphi_var_j[mask]
        
        ####### fitting #######
        if (intercept){S = cbind(S, 1/Q_j**0.5)} #after weighting by Q_j**0.5, will be all 1 intercept 
        
        Q_star_j = fphi_var_j * Q_j**2
        
        inv_res        = inv(t(S) %*% diag(as.vector(Q_j)) %*% S)
        #phi_var_j      = inv_res %*% t(S) %*% diag(as.vector(Q_star_j)) %*% S %*% t(inv_res)
        phi_var_j      = inv_res %*% (t(S) * repmat(t(as.matrix(Q_star_j)),  ncol(S),  1)) %*% S %*% t(inv_res)     
  
        phi_var_j_diag = diag(phi_var_j)[1:s_len] #ignore the last intercept if fitted
        phi_se_j       = sqrt(phi_var_j_diag)
        
        #phi_hat_j   = inv_res %*% t(S) %*% diag(as.vector(Q_j)) %*% X_j
        phi_hat_j   = inv_res %*% (t(S) * repmat(t(as.matrix(Q_j)),  ncol(S),  1)) %*% X_j
        phi_hat_j   = phi_hat_j[1:s_len]#ignore the last intercept if fitted
        phi_test_j  = phi_hat_j/phi_se_j
        phi_pvals_j = 2 * pnorm(abs(phi_test_j), lower.tail=F)
        return(c(phi_hat_j, phi_se_j, phi_pvals_j, mask, fphi_var_j));
        
    }, cl = cl )

    if (parallel) stop_cluster(cl)
                              
    for (j in 1:m){
         
        phi_hat[j,]   <- res[[j]][(0*s_len + 1):(1*s_len)]
        phi_se[j,]    <- res[[j]][(1*s_len + 1):(2*s_len)]
        phi_pvals[j,] <- res[[j]][(2*s_len + 1):(3*s_len)]
        masks[,j]     <- res[[j]][(3*s_len + 1):(3*s_len + n)]
        fphi_var[which(masks[,j] == 1),j]  <- res[[j]][(3*s_len + n + 1):length(res[[j]])]
        assert(length(res[[j]]) == 3*s_len + n + sum(masks[,j]))
    }
                              
    asym.res = list(Q      = t(Q),     # now features by samples
                    masks  = t(masks), # now features by samples
                    fphi_var = fphi_var, 
    
                    phi_hat       = phi_hat,
                    phi_se        = phi_se,
                    phi_hat_pvals = phi_pvals,
                  
                    gammas_hat       = phi_hat[, (k+1):(k + p1*k)],
                    betas_hat        = phi_hat[, (k + p1*k+1):(k + p1*k + p2)],
                    gammas_hat_pvals = phi_pvals[, (k+1):(k + p1*k)],
                    betas_hat_pvals  = phi_pvals[, (k + p1*k+1):(k + p1*k + p2)])
                     
    Unico.mdl[[slot_name]] = asym.res
    flog.info(paste0("Finished asymptotic pvals calculation: ", slot_name))
    return(Unico.mdl)
}
                              
 