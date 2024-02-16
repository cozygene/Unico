#' @importFrom futile.logger appender.tee appender.console flog.appender flog.threshold flog.info flog.debug
#' @importFrom parallel detectCores makeCluster clusterEvalQ stopCluster
#' @importFrom testit assert
#' @importFrom data.table copy
#' @importFrom pracma repmat
#' @importFrom matrixStats rowSds
#' @importFrom utils str
#' @importFrom stats median sd

start_logger <- function(log_file, debug, verbose = TRUE){
  config_level <- if (debug) "debug" else "default"
  Sys.setenv(R_CONFIG_ACTIVE = config_level)
  if (is.null(log_file)){
    invisible(flog.appender(appender.console()))
  }else{
    invisible(flog.appender(appender.tee(log_file)))
  }
  invisible(flog.threshold(if(debug) "DEBUG" else "INFO"))
  if (!verbose) (flog.threshold("ERROR"))
}


get_num_cores <- function(num_cores){
  if (is.null(num_cores)){
    return (max(1,detectCores() - 1))
  }else{
    return (num_cores)
  }
}


init_cluster <- function(num_cores = NULL){
  flog.info("Initiate cluster ...")
  # check if user asked too many cores
	max_cores = get_num_cores(NULL)
  if (is.null(num_cores)){
  	num_cores = max_cores
  }else{
  	if(num_cores > max_cores){
  		num_cores = max_cores
  		flog.info("num_cores is set to maximum number of cores detected: %s .", max_cores)
  	}
  }

  flog.info(paste0("DELETE THIS max core: ", max_cores))
  flog.info(paste0("DELETE THIS number core: ", num_cores))
  flog.info("DELETE THIS right before makeCluster")
  cl <- makeCluster(num_cores)
  flog.info("DELETE THIS right after makeCluster")


  flog.info("Parallel is on with %s nodes.", get_num_cores(num_cores))
  invisible(clusterEvalQ(cl, c(library("pracma"), library("mgcv"), library("nloptr"), library("testit"))))
  flog.info("Packages were loaded into the cluster nodes.")
  return(cl)
}
stop_cluster <- function(cl){
  flog.info("Stop cluster")
  stopCluster(cl)
}


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


validate_inputs_unico = function(X, W, C1, C2, fit_tau,
                                 mean_penalty, var_penalty, covar_penalty,
                                 mean_max_iterations, var_max_iterations, nloptr_opts_algorithm,
                                 max_stds, init_weight, max_u, max_v,
                                 parallel, num_cores, log_file, verbose, debug){

	config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)

  flog.debug("Validating Unico input types ...")
  assert("X must be of class 'matrix'",  is.matrix(X))
  assert("W must be of class 'matrix'",  is.matrix(W))
  assert("C1 must be of class 'matrix'", is.matrix(C1))
  assert("C2 must be of class 'matrix'", is.matrix(C2))
  assert("fit_tau must be of class 'logical'", is.logical(fit_tau))

  assert("mean_penalty must be of class 'numeric'",  is.numeric(mean_penalty))
  assert("var_penalty must be of class 'numeric'",   is.numeric(var_penalty))
  assert("covar_penalty must be of class 'numeric'", is.numeric(covar_penalty))

  assert("mean_max_iterations must be of class 'numeric'",     is.numeric(mean_max_iterations))
  assert("var_max_iterations must be of class 'numeric'",      is.numeric(var_max_iterations))
  assert("nloptr_opts_algorithm must be of class 'character'", is.character(nloptr_opts_algorithm))

  assert("max_stds must be of class 'numeric'",      is.numeric(max_stds))
  assert("init_weight must be of class 'character'", is.character(init_weight))
  assert("max_u must be of class 'numeric'",         is.numeric(max_u))
  assert("max_v must be of class 'numeric'",         is.numeric(max_v))

  assert("parallel must be of class 'logical'",   is.logical(parallel))
  assert("num_cores must be of class 'numeric' or NULL",  is.numeric(num_cores) | is.null(num_cores))
  assert("log_file must be of class 'character' or NULL", is.character(log_file) | is.null(log_file))
  assert("verbose must be of class 'logical'",    is.logical(verbose))
  assert("debug must be of class 'logical'",      is.logical(debug))

  flog.debug("Validating Unico matrix structure ...")
  assert("X must have row names and column names", !is.null(rownames(X)) & !is.null(colnames(X)))
  assert("W must have row names and column names", !is.null(rownames(W)) & !is.null(colnames(W)))
  assert("C1 must have row names and column names (or 0 column)", !is.null(rownames(C1)) & (!is.null(colnames(C1)) | ncol(C1)==0))
  assert("C2 must have row names and column names (or 0 column)", !is.null(rownames(C2)) & (!is.null(colnames(C2)) | ncol(C2)==0))


  flog.debug("Validating Unico matrix dimensions ...")
  assert("The number of columns in X is inconsistent with the number of rows in W", dim(X)[2] == dim(W)[1])
  assert("The number of columns in X is inconsistent with the number of rows in C1", dim(X)[2] == dim(C1)[1])
  assert("The number of columns in X is inconsistent with the number of rows in C2", dim(X)[2] == dim(C2)[1])

  flog.debug("Validating the order of observations across matrices ...")
  assert("The order of observations in W (in the rows) must match the order of the observations in X (in the columns).", all(colnames(X) == rownames(W)))
  assert("The order of observations in C1 (in the rows) must match the order of the observations in X (in the columns).", all(colnames(X) == rownames(C1)))
  assert("The order of observations in C2 (in the rows) must match the order of the observations in X (in the columns).", all(colnames(X) == rownames(C2)))

  flog.debug("Validating Unico input values ...")
  th <- config[["epsilon"]]
  assert(paste0("X must not include features with standard deviation less than ", as.character(th), sep=""), all(rowSds(X) >= th))
  assert("The entries of W must be non-negative and sum up to 1", all(W >= 0 & abs(rowSums(W) - 1) <= th))

  assert("The penalty terms must be non-negative", mean_penalty >=0 & var_penalty >=0 & covar_penalty >=0)
  assert("The iteration terms must be positive integers", mean_max_iterations > 0 & mean_max_iterations%%1 == 0 & var_max_iterations > 0 & var_max_iterations%%1 == 0)
  nloptr_candidate_algs = strsplit(config[["nloptr_candidate_algs"]], ";")[[1]]
  assert(paste0("nloptr_opts_algorithm must be one of the following: ", paste(nloptr_candidate_algs, collapse = " ")),  nloptr_opts_algorithm %in% nloptr_candidate_algs)

  assert("max_stds must be positive", max_stds > 0)
  init_weight_list = c("default", "W_norm", "inv_W_norm")
  assert(paste0("init_weight must be one of the following: ", paste(init_weight_list, collapse = " ")), init_weight %in% c(init_weight_list))
  assert("max_u must be positive", max_u > 0)
  assert("max_v must be positie", max_v > 0)
  if (!is.null(num_cores)) assert("num_cores must be positive integers",  num_cores > 0 & num_cores %%1 ==0)
}


validate_inputs_mdl = function(mdl){
	config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
  flog.debug("Validating W, C1, C2 ...")
  assert("W must be of class 'matrix'",  is.matrix(mdl$W))
  assert("C1 must be of class 'matrix' (can have 0 column)", is.matrix(mdl$C1))
  assert("C2 must be of class 'matrix' (can have 0 column)", is.matrix(mdl$C2))
  assert("W must have row names and column names",  !is.null(rownames(mdl$W)) & !is.null(colnames(mdl$W)))
  assert("C1 must have row names and column names (or 0 column)", !is.null(rownames(mdl$C1)) & (!is.null(colnames(mdl$C1)) | ncol(mdl$C1)==0))
  assert("C2 must have row names and column names (or 0 column)", !is.null(rownames(mdl$C2)) & (!is.null(colnames(mdl$C2)) | ncol(mdl$C2)==0))
  assert("The number of rows in W is inconsistent with the number of rows in C1", dim(mdl$W)[1] == dim(mdl$C1)[1])
  assert("The number of rows in W is inconsistent with the number of rows in C2", dim(mdl$W)[1] == dim(mdl$C2)[1])
  assert("The order of observations in C1 must match the order of the observations in W.", all(rownames(mdl$W) == rownames(mdl$C1)))
  assert("The order of observations in C2 must match the order of the observations in W.", all(rownames(mdl$W) == rownames(mdl$C2)))
  assert("The entries of W must be non-negative and sum up to 1", all(mdl$W >= 0 & abs(rowSums(mdl$W) - 1) <= config[["epsilon"]]))

  flog.debug("Validating model parameters types ...")
  assert("mus_hat must be of class 'matrix'",      is.matrix(mdl$mus_hat))
  assert("gammas_hat must be of class 'matrix'",   is.matrix(mdl$gammas_hat))
  assert("betas_hat must be of class 'matrix'",    is.matrix(mdl$betas_hat))
  assert("sigmas_hat must be of class 'array'",    is.array(mdl$sigmas_hat))
  assert("taus_hat must be of class 'matrix'",     is.matrix(mdl$taus_hat))
  assert("scale.factor must be of class 'matrix'", is.matrix(mdl$scale.factor))

  flog.debug("Validating model parameters structure ...")
  assert("mus_hat must have row names and column names", !is.null(rownames(mdl$mus_hat)) & !is.null(colnames(mdl$mus_hat)))
  assert("gammas_hat must have row names", !is.null(rownames(mdl$gammas_hat)))
  assert("betas_hat must have row names", !is.null(rownames(mdl$betas_hat)))
  assert("sigmas_hat must have all dimension names", !is.null(dimnames(mdl$sigmas_hat)[[1]]) & !is.null(dimnames(mdl$sigmas_hat)[[2]]) & !is.null(dimnames(mdl$sigmas_hat)[[3]]))
  assert("taus_hat must have row names", !is.null(rownames(mdl$taus_hat)))
  assert("scale.factor must have row names", !is.null(rownames(mdl$scale.factor)))

  flog.debug("Validating model parameters dimensions ...")
  feature.ids = rownames(mdl$mus_hat)
  source.ids  = colnames(mdl$mus_hat)

  assert("mus_hat must be features by sources", all(rownames(mdl$mus_hat) == feature.ids) & all(colnames(mdl$mus_hat) == source.ids))
  assert("gammas_hat must have features on the row", all(rownames(mdl$gammas_hat) == feature.ids))
  assert("betas_hat must have features on the row", all(rownames(mdl$betas_hat) == feature.ids))
  assert("sigmas_hat must have features by sources by sources", all(dimnames(mdl$sigmas_hat)[[1]] == feature.ids) & all(dimnames(mdl$sigmas_hat)[[2]] == source.ids) & all(dimnames(mdl$sigmas_hat)[[3]] == source.ids))
  assert("taus_hat must have features on the row", all(rownames(mdl$taus_hat) == feature.ids))
  assert("scale.factor must have features on the row", all(rownames(mdl$scale.factor) == feature.ids))

  flog.debug("Validating model parameters values ...")
  assert("The entries of taus_hat must be non-negative", all(mdl$taus_hat >= 0))
  assert("The entries of scale.factor must be non-negative", all(mdl$scale.factor >= 0))
}


validate_inputs_tensor = function(X, W, C1, C2, Unico.mdl, parallel, num_cores,
                                  log_file, verbose, debug){

	config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
  flog.debug("Validating Unico inputs ...")
  validate_inputs_unico(X = X, W = W, C1 = C1, C2 = C2, parallel = parallel, num_cores = num_cores,
                        log_file = log_file, debug = debug, verbose = verbose,

                        #rest are set to valid inputs
                        fit_tau = T,
                        mean_penalty = 0, var_penalty = 0, covar_penalty = 0,
                        mean_max_iterations = 1, var_max_iterations = 1, nloptr_opts_algorithm = "NLOPT_LN_COBYLA",
                        max_stds = 2, init_weight = "default", max_u = 1, max_v = 1)

  flog.debug("Validating model parameters ...")
  validate_inputs_mdl(Unico.mdl)

  flog.debug("Validating inputs are consistent ...") # interally already consistent, so checking just one pair
  # there is no restriction on matching sample.ids between tensor input and learned model parameters,
  # allowing a model to be applied on samples not used in the parameter estimation.
  assert("features must align between inputs and model parameters", all(rownames(X) == rownames(Unico.mdl$mus_hat)))
  assert("sources must align between inputs and model parameters", all(colnames(W) == colnames(Unico.mdl$mus_hat)))
}


validate_inputs_parametric = function(X, Unico.mdl, slot_name,
                                      diag_only, intercept,
                                      X_max_stds, Q_max_stds, XQ_max_stds,
                                      parallel, num_cores, log_file, debug, verbose){

	config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
  flog.debug("Validating Unico inputs ...")
  validate_inputs_unico(X = X, W = Unico.mdl$W, C1 = Unico.mdl$C1, C2 = Unico.mdl$C2,
                        #checks that X has the sample.ids in the correct order as learned model parameters
                        parallel = parallel, num_cores = num_cores,
                        log_file = log_file, debug = debug, verbose = verbose,
                        #rest are set to valid inputs
                        fit_tau = T,
                        mean_penalty = 0, var_penalty = 0, covar_penalty = 0,
                        mean_max_iterations = 1, var_max_iterations = 1, nloptr_opts_algorithm = "NLOPT_LN_COBYLA",
                        max_stds = 2, init_weight = "default", max_u = 1, max_v = 1)

  flog.debug("Validating model parameters ...")
  validate_inputs_mdl(Unico.mdl)

  flog.debug("Validating parametric inputs ...")
  assert("slot_name must be of type 'character'", is.character(slot_name))
  assert("diag_only must be of type 'logical'", is.logical(diag_only))
  assert("intercept must be of type 'logical'", is.logical(intercept))
  assert("X_max_stds must be of type 'numeric' and positive", is.numeric(X_max_stds) & X_max_stds > 0)
  assert("Q_max_stds must be of type 'numeric' and positive", is.numeric(Q_max_stds) & X_max_stds > 0)
  assert("XQ_max_stds must be of type 'numeric' and positive", is.numeric(XQ_max_stds) & X_max_stds > 0)
}


validate_inputs_asymptotic = function(X, Unico.mdl, slot_name,
                                      diag_only, intercept,
                                      X_max_stds, Q_max_stds, V_min_qlt,
                                      parallel, num_cores, log_file, debug, verbose){

	config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
  flog.debug("Validating Unico inputs ...")
  validate_inputs_unico(X = X, W = Unico.mdl$W, C1 = Unico.mdl$C1, C2 = Unico.mdl$C2,
                        #checks that X has the sample.ids in the correct order as learned model parameters
                        parallel = parallel, num_cores = num_cores,
                        log_file = log_file, debug = debug, verbose = verbose,
                        #rest are set to valid inputs
                        fit_tau = T,
                        mean_penalty = 0, var_penalty = 0, covar_penalty = 0,
                        mean_max_iterations = 1, var_max_iterations = 1, nloptr_opts_algorithm = "NLOPT_LN_COBYLA",
                        max_stds = 2, init_weight = "default", max_u = 1, max_v = 1)

  flog.debug("Validating model parameters ...")
  validate_inputs_mdl(Unico.mdl)

  flog.debug("Validating parametric inputs ...")
  assert("slot_name must be of type 'character'", is.character(slot_name))
  assert("diag_only must be of type 'logical'", is.logical(diag_only))
  assert("intercept must be of type 'logical'", is.logical(intercept))
  assert("X_max_stds must be of type 'numeric' and positive", is.numeric(X_max_stds) & X_max_stds > 0)
  assert("Q_max_stds must be of type 'numeric' and positive", is.numeric(Q_max_stds) & X_max_stds > 0)
  assert("V_min_qlt must be of type 'numeric' and between 0 and 1", is.numeric(V_min_qlt) & V_min_qlt > 0 & V_min_qlt < 1)
}


# data: any data structure, could be 2d matrix or 3d array
# min.val: numeric. minimal absolute value to be capped
# max.val: numeric. maximal absolute value to be capped
# NA in the data structure is turned into + min.val
# absolute values larger than max.val or smaller than min.val (including Inf, -Inf)are replaced by max.val and min.val
# sign information is also preserved
cap_values = function(data, min.val = 10**(-4), max.val = 10**(4)){
    flog.info("Capping extreme values in estimated parameters")
    if(sum(is.na(data)| is.infinite(data)) != 0){
        flog.debug(paste0(sum(is.na(data)), " NAN"))
        data[is.na(data)] <- min.val
    }
    if(sum(abs(data) < min.val) != 0){
        flog.debug(paste0(sum(abs(data) < min.val), " extremely close to 0 values"))
        sign.info = sign(data[which(abs(data) < min.val)])
        data[which(abs(data) < min.val)] <- sign.info * min.val
    }
    if(sum(abs(data) > max.val) != 0){
        flog.debug(paste0(sum(abs(data) < max.val), " extremely large values"))
        sign.info = sign(data[which(abs(data) > max.val)])
        data[which(abs(data) > max.val)] <- sign.info * max.val
    }
    return(data)
}


assert_cov = function(C){
	assert("matrix has to be square", ncol(C) == nrow(C))
	assert("matrix has to be symmetric", all((C - t(C)) < 0.0001))
	assert("Not all diag elements in a covariance matrix is non-negative", all(diag(C) >= - 10**(-8)))
	assert("matrix in not PD: not all eigen values are greater than 0",all(eigen(C)$values > - 10**(-8)))
}

#' @importFrom mgcv pcls
#' @importFrom nloptr nloptr
#' @importFrom matrixcalc hadamard.prod
#' @importFrom pracma Reshape

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


mean_updates_j = function (X_j, W, C1, C2, design_mat, variables, Us, config){

	k  = dim(W)[2]
	p1 = dim(C1)[2]
	p2 = dim(C2)[2]

	# penalty is proportional to weighted total sample size and average/median value
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


get_X_j_var_emp = function (X_j, mus_j, W, C1, gammas_j, C2, betas_j){
	X_j_ept = get_X_j_ept (mus_j, W, C1, gammas_j, C2, betas_j)
	X_j_var_emp = (X_j - X_j_ept)**2
	return(X_j_var_emp)
}


var_objective = function(WW_tensor, X_j_var_emp, Vs, variables, config, return_grads = TRUE){

	Ls_j   = matrix(variables[1:(length(variables)-1)])
	taus_j = variables[length(variables)]

	n = dim(WW_tensor)[1]
	k = dim(WW_tensor)[2]

	# proportional to weighted sample size and proportional to feature j's variance size
	covar_penalty = config[["covar_penalty"]] * sum(Vs**2) * median(X_j_var_emp)
	var_penalty   = config[["var_penalty"]]   * sum(Vs**2) * median(X_j_var_emp)

	L_j_mat = vec_2_tril(Ls_j,k)
	sigmas_j = L_j_mat %*% t(L_j_mat)
	X_j_var_ana = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_j) + taus_j, numeric(1)))

	objective = sum((Vs**2)*((X_j_var_emp - X_j_var_ana)**2))
	# add entire sigmas but then remove the extra penalty on diagnal
	objective = objective + sum(covar_penalty*(sigmas_j**2)) - sum((covar_penalty - var_penalty)*(diag(sigmas_j)**2))

	# TODO: finish the gradiant calculation with penalty terms
	if (return_grads){
		## new implementation
		grad_coeff = Vs**2 * 2 * (X_j_var_emp - X_j_var_ana) * (-1)
		# braodcast
		grad_coeff_tensor = array(0,c(n,k,k))
		grad_coeff_tensor[1:n,,] = grad_coeff
		# multiply by 2 here for the gradients of Ls and collapse across individuals
		#grads = base::colSums(2 * grad_coeff_tensor * WW_tensor, dims = 1)

		grads = tril_2_vec(grads %*% L_j_mat)
		# TODO: adding penalty
		# grads = grads + ......
		# adding gradients of taus
		grads = c(grads, sum(grad_coeff))
		return(list("objective" = objective,
								"gradient"  = grads))
	}else{
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


# parameter fitting
fit_mean_var = function(X_j, W, C1, C2, mus_j, gammas_j, betas_j, Ls_j, taus_j, fit_tau,
												mean_max_iterations, var_max_iterations, nloptr_opts_algorithm, config, debug){

	k  = dim(W)[2]
	n  = dim(C1)[1]
	p1 = dim(C1)[2]
	p2 = dim(C2)[2]

	#mean updates
	if (config[["init_weight"]] == "W_norm"){
		Us = as.matrix(rowSums(W**2))
	}else if(config[["init_weight"]] == "inv_W_norm"){
		Us = 1/as.matrix(rowSums(W**2))
	}else{
		Us = matrix(1,n,1)
	}

	#masking outliers for mean updates
	max_stds = config[["max_stds"]]
	if(!is.null(max_stds)){
		sd_j <- sd(X_j)
		outlier_mask = abs(X_j - mean(X_j)) >= max_stds*sd_j
		Us[outlier_mask] = 0
	}

	if(debug) Us_j_list = c()
	for (itr in 1:mean_max_iterations){

		# TODO; separate the format_j function, design mat is constant
		format_j      = format_mean_param_j(W, C1, C2, mus_j, gammas_j, betas_j)
		design_mat_j  = as.matrix(format_j$"design_mat")
		variables_j   = format_j$"variables"

		res     = mean_updates_j(X_j, W, C1, C2, design_mat_j, variables_j, Us, config)
		mus_j   = res[1:k]
		gammas_j= res[seq(k+1, k+k*p1, length=k*p1)]
		betas_j = res[seq(k+k*p1+1, k+k*p1+p2, length=p2)]

		# track
		if(debug) Us_j_list = cbind(Us_j_list, Us)

		# update Us
		Us = abs(1/(X_j - design_mat_j %*% res))
		Us = as.matrix(replace(Us, which (Us > config[["max_u"]]), config[["max_u"]]))
		Us = as.matrix(Us)
		Us = Us/max(Us)

		# update outliers mask
		if(!is.null(max_stds)){
			Us[outlier_mask] = 0
		}
	}

	#variance updates
	WW_tensor = get_WW_tensor(W)
	if (config[["init_weight"]] == "W_norm"){
		Vs = as.matrix(rowSums(W**2))
	}else if(config[["init_weight"]] == "inv_W_norm"){
		Vs = 1/as.matrix(rowSums(W**2))
	}else{
		Vs = matrix(1,n,1)
	}

	#masking outliers for variance updates
	max_stds = config[["max_stds"]]
	if(!is.null(max_stds)){
		sd_j <- sd(X_j)
		outlier_mask = abs(X_j - mean(X_j)) >= max_stds*sd_j
		Vs[outlier_mask] = 0
	}

	#prepare input to nlpotr
	bounds = get_bounds(k, fit_tau)
	#get scope_type and grad_type
	scope_type = strsplit(strsplit(nloptr_opts_algorithm, split = "_")[[1]][2], "")[[1]][1]
	grad_type  = strsplit(strsplit(nloptr_opts_algorithm, split = "_")[[1]][2], "")[[1]][2]
	return_grads = grad_type == "D"
	opts = get_opts(nloptr_opts_algorithm, scope_type, grad_type, config)
	X_j_var_emp = get_X_j_var_emp(X_j, mus_j, W, C1, gammas_j, C2, betas_j)

	if(debug){
		Vs_j_list     = list(var_max_iterations)
		Ls_j_list     = list(var_max_iterations)
		sigmas_j_list = list(var_max_iterations)
	}

	for (itr in 1:var_max_iterations){
		res <- nloptr(x0 = c(Ls_j, taus_j),
									eval_f = function(x, WW_tensor, X_j_var_emp, Vs, return_grads) var_objective(WW_tensor, X_j_var_emp, Vs, x, config, return_grads),
									lb = bounds$"lb",
									ub = bounds$"ub",
									WW_tensor = WW_tensor,
									X_j_var_emp = X_j_var_emp,
									Vs = Vs,
									return_grads = return_grads,
									opts=opts)

		Ls_j  = res$solution[1:(length(res$solution)-1)]
		taus_j = res$solution[length(res$solution)]
		sigmas_j = vec_2_tril(Ls_j,k) %*% t(vec_2_tril(Ls_j,k))

		# track
		if(debug){
			Vs_j_list[[itr]]     = as.vector(Vs)
			Ls_j_list[[itr]]     = as.vector(Ls_j)
			sigmas_j_list[[itr]] = as.vector(sigmas_j)
		}

		# update Vs
		X_j_var_ana = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_j) + taus_j, numeric(1)))
		Vs = abs(1/(X_j_var_emp - X_j_var_ana))
		Vs = as.matrix(replace(Vs, which (Vs > config[["max_v"]]), config[["max_v"]]))
		Vs = Vs/max(Vs)

		# update outliers mask
		if(!is.null(max_stds)){
			Vs[outlier_mask] = 0
		}
	}
	params = list("mus_j"     = mus_j,
								"gammas_j"  = gammas_j,
								"betas_j"   = betas_j,
								"Ls_j"      = Ls_j,
								"taus_j"    = taus_j,
								"Us"        = Us,
								"Vs"        = Vs)
	if(debug){
		params = c(params, list("Us_j_list" = Us_j_list, "Vs_j_list" = Vs_j_list,
														"Ls_j_list" = Ls_j_list, "sigmas_j_list" = sigmas_j_list))
	}
	return(params)
}


# related to estimate Z
get_Z_j = function(X_j, W, C1, C2, mus_j, gammas_j, betas_j, sigmas_hat_j, taus_j, WW_tensor){

	config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
	n = dim(W)[1]
	k = dim(W)[2]

	X_j_ept     = get_X_j_ept (mus_j, W, C1, gammas_j, C2, betas_j)
	X_j_var_ana = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * sigmas_hat_j) + taus_j, numeric(1)))
	#if no variance in any entry (X_j_var_ana is n by 1)
	#(rare but consider celltype1 is estimated to have 0 variance and sample's only has celltype1 and no other celltypes)
	#setting pure 0 variance people to have pseudo variance 1. exact number does not matter
	#the final Z_j takes Z_ept
	if (sum(X_j_var_ana <= config[["epsilon"]]) !=0){X_j_var_ana[X_j_var_ana <= config[["epsilon"]]] = 1}

	covar_z_x = W %*% sigmas_hat_j
	X_j_var_ana_inv = -(1/X_j_var_ana)
	B = repmat(X_j_var_ana_inv, 1, k) * covar_z_x
	Z_ept = t(repmat(as.matrix(mus_j), 1,n)) + C1  %*% t(matrix(gammas_j, nrow = dim(W)[2], byrow = TRUE))
	X_diff = X_j_ept - X_j
	delta = B * repmat(X_diff, 1, k)
	Z_j = Z_ept  +  delta
	return(Z_j)
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
	sigmas = sd(x)
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
