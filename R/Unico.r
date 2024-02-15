#' @title Fitting the Unico model
#'
#' @description Fits the Unico model for an input matrix of features by observations that are coming from a mixture of \code{k} sources, under the assumption that each observation is a mixture of unique (unobserved) source-specific values (in each feature in the data). Specifically, for each feature, it standardizes the data and learns the source-specific mean and full \code{k} by \code{k} variance-covariance matrix.
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations.
#' @param W An \code{n} by \code{k} matrix of weights - the weights of \code{k} sources for each of the \code{n} mixtures (observations). All the weights must be positive and each row - corresponding to the weights of a single observation - must sum up to 1. Note that \code{W} must include row names and column names and that NA values are currently not supported.
#' @param C1 An \code{n} by \code{p1} design matrix of covariates that may affect the hidden source-specific values (possibly a different effect size in each source). Note that \code{C1} must include row names and column names and should not include an intercept term. NA values are currently not supported. Note that each covariate in \code{C1} results in \code{k} additional parameters in the model of each feature, therefore, in order to alleviate the possibility of model overfitting, it is advised to be mindful of the balance between the size of \code{C1} and the sample size in \code{X}.
#' @param C2 An \code{n} by \code{p2} design matrix of covariates that may affect the mixture (i.e. rather than directly the sources of the mixture; for example, variables that capture biases in the collection of the measurements). Note that \code{C2} must include row names and column names and should not include an intercept term. NA values are currently not supported.
#' @param fit_tau A logical value indicating whether to fit the standard deviation of the measurement noise (i.e. the i.i.d. component of variation in the model denoted as \eqn{\tau}).
#' @param mean_penalty A non-negative numeric value indicating the regularization strength on the source-specific mean estimates.
#' @param var_penalty A non-negative numeric value indicating the regularization strength on the diagonal entries of the full \code{k} by \code{k} variance-covariance matrix.
#' @param covar_penalty A non-negative numeric value indicating the regularization strength on the off diagonal entries of the full \code{k} by \code{k} variance-covariance matrix.
#' @param mean_max_iterations A non-negative numeric value indicating the number of iterative updates performed on the mean estimates.
#' @param var_max_iterations A non-negative numeric value indicating the number of iterative updates performed on the variance-covariance matrix.
#' @param nloptr_opts_algorithm A string indicating the optimization algorithm to use.
#' @param max_stds A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers. Only samples within \code{max_stds} standard deviations from the mean will be used for the moments estimation of a given feature.
#' @param init_weight A string indicating the initial weights on the samples to start the iterative optimization.
#' @param max_u A non-negative numeric value indicating the maximum weights/influence a sample can have on mean estimates.
#' @param max_v A non-negative numeric value indicating the maximum weights/influence a sample can have on variance-covariance estimates.
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; set \code{debug} to \code{TRUE} before reporting issues.
#'
#' @importFrom futile.logger flog.info
#' @importFrom pbapply pboptions pblapply
#' @importFrom parallel clusterExport
#'
#' @details Unico assumes the following model: \deqn{X_{ij} =  w_{i}^T Z_{ij} +(c_i^{(2)})^T \beta_j+ e_{ij}}
#' The mixture value at sample \eqn{i} feature \eqn{j}: \eqn{X_{ij}} is modeled as a weighted linear combination, specified by weights \eqn{w_i = (w_{i1},...,w_{ik})}, of a total of \eqn{k} source-specific levels, specified by \eqn{Z_{ij}=(Z_{ij1},...,Z_{ijk})}.
#' In addition, we also consider global-level covariates \eqn{c_i^{(2)}} that systematically affect the observed mixture values and their effect sizes \eqn{\beta_j}. \eqn{e_{ij}} denotes the i.i.d measurement noise with variance \eqn{\tau} across all samples.
#' Weights have be to non-negative and sum up to \eqn{1} across all sources for each sample. In practice, we assume that the weights are fixed and estimated by external methods.
#'
#' Source specific profiles are further modeled as: \deqn{Z_{ijh} = \mu_{jh} + (c_i^{(1)})^T \gamma_{jh} + \epsilon_{ijh}}
#' \eqn{\mu_{jh}} denotes the population level mean of feature \eqn{j} at source \eqn{h}.
#' We also consider covariates \eqn{c_i^{(1)}} that systematically affect the source-specific values and their effect sizes \eqn{\gamma_{jh}} on each source.
#' Finally, we actively model the \eqn{k} by \eqn{k} covariance structure of a given feature \eqn{j} across all \eqn{k} sources \eqn{Var[\vec{\epsilon_{ij}}] = \Sigma_{j} \in \mathbf{R}^{k \times k}}.
#'
#'
#' @return A list with the estimated parameters of the model. This list can be then used as the input to other functions such as \link{tensor}.
#' \item{W}{An \code{n} by \code{k} matrix of weights. This is the same as \code{W} from input.}
#' \item{C1}{An \code{n} by \code{p1} design matrix of source-specific covariates. This is the same as \code{C1} from input.}
#' \item{C2}{An \code{n} by \code{p2} design matrix of not source-specific covariates. This is the same as \code{C2} from input.}
#' \item{mus_hat}{An \code{m} by \code{k} matrix of estimates for the mean of each source in each feature.}
#' \item{gammas_hat}{An \code{m} by \code{k*p1} matrix of the estimated effects of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} covariates on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
#' \item{betas_hat}{An \code{m} by \code{p2} matrix of the estimated effects of the \code{p2} covariates in \code{C2} on the mixture values of each of the \code{m} features in \code{X}.}
#' \item{sigmas_hat}{An \code{m} by \code{k} by \code{k} tensor of estimates for the cross source \code{k} by \code{k} variance-covariance matrix in each feature.}
#' \item{taus_hat}{An \code{m} by \code{1} matrix of estimates for the variance of the measurement noise.}
#' \item{scale.factor}{An \code{m} by \code{1} matrix of scaling factors for standardizing each feature.}
#' \item{config}{A list with hyper-parameters used for fitting the model and configurations for in the optimization algorithm.}
#' \item{Us_hat_list}{A list tracking, for each feature, the sample weights used for each iteration of the mean optimization  (activated only if \code{debug == TRUE}).}
#' \item{Vs_hat_list}{A list tracking, for each feature, the sample weights used for each iteration of the variance-covariance optimization  (activated only if \code{debug == TRUE}).}
#' \item{Ls_hat_list}{A list tracking, for each feature, the computed estimates of the upper triangular cholesky decomposition of variance-covariance matrix at each iteration of the variance-covariance optimization (activated only if \code{debug == TRUE}).}
#' \item{sigmas_hat_list}{A list tracking, for each feature, the computed estimates of the variance-covariance matrix at each iteration of the variance-covariance optimization (activated only if \code{debug == TRUE}).}
#'
#' @examples
#' data = simulate_data(n=100, m=2, k=3, p1=1, p2=1, taus_std=0, log_file=NULL)
#' res = list()
#' res$params.hat = Unico(data$X, data$W, data$C1, data$C2, parallel=FALSE, log_file=NULL)
#'
#' @export Unico
Unico = function(X, W, C1, C2, fit_tau = FALSE,
                 mean_penalty = 0, var_penalty = 0.01, covar_penalty = 0.01,
                 mean_max_iterations = 2, var_max_iterations = 3, nloptr_opts_algorithm = "NLOPT_LN_COBYLA",
                 max_stds = 2, init_weight = "default", max_u = 1, max_v = 1,
                 parallel = TRUE, num_cores = NULL, log_file = "Unico.log", verbose = FALSE, debug = FALSE){

    start_logger(log_file, debug, verbose)
    flog.info("Starting Unico ...")
    config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
    op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")

    n = dim(X)[2]
    m = dim(X)[1]
    k = dim(W)[2]
    sample.ids = rownames(W)
    feature.ids= rownames(X)
    source.ids = colnames(W)

    if (is.null(C1)) C1 <- matrix(0, nrow=nrow(W), ncol=0)
    if (is.null(C2)) C2 <- matrix(0, nrow=nrow(W), ncol=0)
    rownames(C1) = sample.ids
    rownames(C2) = sample.ids
    p1 = dim(C1)[2]
    p2 = dim(C2)[2]
    C1.ids = colnames(C1)
    C2.ids = colnames(C2)
    C1.ids_ <- matrix("",k*p1,1)
    if (p1){
        for (j in 1:k){
            C1.ids_[((j-1)*p1+1):(p1*j)] <- unlist(lapply(C1.ids,function(i) paste(source.ids[j],".",i,sep="")))
        }
    }

    flog.info("Validate Unico inputs")
    validate_inputs_unico(X, W, C1, C2, fit_tau,
                          mean_penalty, var_penalty, covar_penalty,
                          mean_max_iterations, var_max_iterations, nloptr_opts_algorithm,
                          max_stds, init_weight, max_u, max_v,
                          parallel, num_cores, log_file, debug, verbose)

    config[["mean_max_iterations"]] <- mean_max_iterations
    config[["var_max_iterations"]]  <- var_max_iterations
    config[["mean_penalty"]]  <- mean_penalty
    config[["var_penalty"]]   <- var_penalty
    config[["covar_penalty"]] <- covar_penalty
    config[["max_u"]]         <- max_u
    config[["max_v"]]         <- max_v
    config[["max_stds"]]      <- max_stds
    config[["init_weight"]]   <- init_weight

    # here calculate the scale.factor and later save this information in mdl
    # scale and data to unit variance and move on, all downstream calculation like parameter estimatation are on the scaled data
    # only use scale.factor at the tensor stage to scale it back to original size
    scale.max_stds = max_stds
    scale.min.thr = config[["epsilon"]]
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
        flog.info(paste0("Certain features have extremely close to 0 varinace in X, set their stds to ", scale.min.thr))
        scale.factor[scale.factor < scale.min.thr,] = scale.min.thr
    }

    X = scale_feature_matrix(X, 1/scale.factor)
    #make it samples by features
    X = t(X)

    #initialization
    mus_hat    = matrix(1, m, k)
    rownames(mus_hat) = feature.ids
    colnames(mus_hat) = source.ids

    gammas_hat = matrix(0, m, k * p1) # already flatten to [row(ct)1, row(ct)2, ...]
    rownames(gammas_hat) = feature.ids
    colnames(gammas_hat) = C1.ids_

    betas_hat  = matrix(0, m, p2)
    rownames(betas_hat) = feature.ids
    colnames(betas_hat) = C2.ids

    Ls_hat     = matrix(1, m, ((k + 1) * k / 2))
    sigmas_hat = array(0, c(m,k,k))
    dimnames(sigmas_hat)[[1]] = feature.ids
    dimnames(sigmas_hat)[[2]] = source.ids
    dimnames(sigmas_hat)[[3]] = source.ids

    taus_hat   = matrix(0, m, 1)
    rownames(taus_hat) = feature.ids
    colnames(taus_hat) = c("taus_hat")

    # tracking
    if(debug){
        Us_hat_list     = list(mean_max_iterations)
        Vs_hat_list     = list(var_max_iterations)
        Ls_hat_list     = list(var_max_iterations)
        sigmas_hat_list = list(var_max_iterations)
        for (itr in 1:mean_max_iterations){
            Us_hat_list[[itr]]     = matrix(0, m, n)
            rownames(Us_hat_list[[itr]]) = feature.ids
            colnames(Us_hat_list[[itr]]) = sample.ids
        }
        for (itr in 1:var_max_iterations){
            Vs_hat_list[[itr]]     = matrix(0, m, n)
            rownames(Vs_hat_list[[itr]]) = feature.ids
            colnames(Vs_hat_list[[itr]]) = sample.ids
            Ls_hat_list[[itr]]     = matrix(0, m, (k + 1) * k / 2)
            sigmas_hat_list[[itr]] = matrix(0, m, k*k)
        }
    }

    flog.info("Starting parameter learning ...")
    cl <- if (parallel) init_cluster(num_cores) else NULL
    if (parallel) clusterExport(cl, c("X", "W", "C1", "C2",
                                      "mus_hat", "gammas_hat", "betas_hat", "Ls_hat", "taus_hat", "fit_tau",
                                      "mean_max_iterations", "var_max_iterations", "nloptr_opts_algorithm", "config", "debug",
                                      "format_mean_param_j","mean_updates_j", "get_X_j_ept",
                                      "get_opts", "get_bounds", "var_objective", "get_X_j_var_emp",
                                      "get_WW_tensor", "tril_2_vec", "vec_2_tril", "fit_mean_var"), envir=environment())


    res <- pblapply(1:m,function(j) fit_mean_var(X[,j], W, C1, C2,
                                                 mus_hat[j,], gammas_hat[j,], betas_hat[j,], Ls_hat[j,],taus_hat[j], fit_tau,
                                                 mean_max_iterations, var_max_iterations, nloptr_opts_algorithm, config, debug),
                    cl = cl)

    if (parallel) stop_cluster(cl)

    flog.info("Formating results ...")
    for (j in 1:m){
        mus_hat[j,]     = res[[j]][["mus_j"]]
        gammas_hat[j,]  = res[[j]][["gammas_j"]]
        betas_hat[j,]   = res[[j]][["betas_j"]]
        Ls_hat[j,]      = res[[j]][["Ls_j"]]
        sigmas_hat[j,,] = vec_2_tril(Ls_hat[j,], k) %*%  t(vec_2_tril(Ls_hat[j,], k))
        taus_hat[j]     = res[[j]][["taus_j"]]

        if(debug){
            for (itr in 1:mean_max_iterations){
                Us_hat_list[[itr]][j,]     = res[[j]][["Us_j_list"]][[itr]]
            }
            for (itr in 1:var_max_iterations){
                Vs_hat_list[[itr]][j,]     = res[[j]][["Vs_j_list"]][[itr]]
                Ls_hat_list[[itr]][j,]     = res[[j]][["Ls_j_list"]][[itr]]
                sigmas_hat_list[[itr]][j,] = res[[j]][["sigmas_j_list"]][[itr]]
            }
        }
    }
    sigmas_hat = cap_values(sigmas_hat,  max.val = config[["max_var"]], min.val = config[["min_var"]])

    Unico.mdl = list(W = W, C1 = C1, C2 = C2,
                     mus_hat = mus_hat, gammas_hat = gammas_hat, betas_hat = betas_hat,
                     sigmas_hat= sigmas_hat, taus_hat = taus_hat,
                     scale.factor = scale.factor, config = config)
    if(debug){
        Unico.mdl = c(Unico.mdl, list(Us_hat_list = Us_hat_list, Vs_hat_list = Vs_hat_list,
                                      Ls_hat_list = Ls_hat_list, sigmas_hat_list = sigmas_hat_list))
    }

    flog.info("Finished parameter learning")
    return(Unico.mdl)
}


#' @title Inferring the underlying source-specific 3D tensor
#'
#' @description Infers the underlying (sources by features by observations) 3D tensor from the observed (features by observations) 2D mixture, under the assumption of the Unico model that each observation is a mixture of unique source-specific values (in each feature in the data). In the context of bulk genomics containing a mixture of cell types (i.e. the input could be CpG sites by individuals for DNA methylation and genes by individuals for RNA expression), \code{tensor} allows to estimate the cell-type-specific levels for each individual in each CpG site/gene (i.e. a tensor of CpG sites/genes by individuals by cell types).
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations. Note that \code{X} could potentially be different from the \code{X} used to learn \code{Unico.mdl} (i.e. the original observed 2D mixture used to fit the model).
#' @param W An \code{n} by \code{k} matrix of weights - the weights of \code{k} sources for each of the \code{n} mixtures (observations). All the weights must be positive and each row - corresponding to the weights of a single observation - must sum up to 1. Note that \code{W} must include row names and column names and that NA values are currently not supported.
#' @param C1 An \code{n} by \code{p1} design matrix of covariates that may affect the hidden source-specific values (possibly a different effect size in each source). Note that \code{C1} must include row names and column names and should not include an intercept term. NA values are currently not supported. Note that all covariates in \code{C1} must be present and match the order of the set of covariates in \code{C1} stored in \code{Unico.mdl} (i.e. the original set of source-specific covariates available when initially fitting the model).
#' @param C2 An \code{n} by \code{p2} design matrix of covariates that may affect the mixture (i.e. rather than directly the sources of the mixture; for example, variables that capture biases in the collection of the measurements). Note that \code{C2} must include row names and column names and should not include an intercept term. NA values are currently not supported. Note that all covariates in \code{C2} must be present and match the order of the set of covariates in \code{C2} stored in \code{Unico.mdl} (i.e. the original set of not source-specific covariates available when initially fitting the model).
#' @param Unico.mdl The entire set of model parameters estimated by Unico on the 2D mixture matrix (i.e. the list returned by applying function \code{Unico} to \code{X}).
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; set \code{debug} to \code{TRUE} before reporting issues.
#'
#' @importFrom futile.logger flog.info
#' @importFrom pbapply pboptions pblapply
#' @importFrom parallel clusterExport
#'
#' @details After obtaining all the estimated parameters in the Unico model (by calling \link{Unico}), \code{tensor} uses the conditional distribution \eqn{Z_{jh}^i|X_{ij}=x_{ij}} for estimating the \eqn{k} source-specific levels of each sample \eqn{i} at each feature \eqn{j}.
#
#' @return A \code{k} by \code{m} by \code{n} array with the estimated source-specific values. The first axis/dimension in the array corresponds to the different sources.
#'
#' @examples
#' data = simulate_data(n=100, m=2, k=3, p1=1, p2=1, taus_std=0, log_file=NULL)
#' res = list()
#' res$params.hat = Unico(data$X, data$W, data$C1, data$C2, parallel=FALSE, log_file=NULL)
#' res$Z = tensor(data$X, data$W, data$C1, data$C2, res$params.hat, parallel=FALSE, log_file=NULL)
#'
#' @export tensor
tensor = function(X, W, C1, C2, Unico.mdl, parallel = TRUE, num_cores = NULL,
                  log_file = "Unico.log", verbose = FALSE, debug = FALSE){

    start_logger(log_file, debug, verbose)
	  config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)
    op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")

    n  = dim(W)[1]
    m  = dim(X)[1]
    k  = dim(W)[2]
    source.ids  = colnames(W)
    feature.ids = rownames(X)
    sample.ids  = colnames(X)

    if (is.null(C1)) C1 <- matrix(0, nrow=nrow(W), ncol=0)
    if (is.null(C2)) C2 <- matrix(0, nrow=nrow(W), ncol=0)
    rownames(C1) = sample.ids
    rownames(C2) = sample.ids
    p1 = dim(C1)[2]
    p2 = dim(C2)[2]

    flog.info("Validate tensor inputs ...")
    validate_inputs_tensor(X, W, C1, C2, Unico.mdl, parallel, num_cores,
                           log_file, verbose, debug)

    flog.info("Starting tensor ...")
    #transfer it to scaled version
    X = X/repmat(Unico.mdl$scale.factor[rownames(X), , drop = F], 1, ncol(X))
    X = t(X)

    #initialize
    Z_hat   = array(rep(0, m * n * k), c(m,n,k))
    #shared calculation across features
    WW_tensor = get_WW_tensor(W)

    # parallel
    cl <- if (parallel) init_cluster(num_cores) else NULL
    if (parallel) clusterExport(cl, c("X","W","C1", "C2", "Unico.mdl",
                                      "format_mean_param_j","get_X_j_ept", "get_X_j_var_emp", "get_WW_tensor",
                                      "tril_2_vec", "vec_2_tril",
                                      "get_Z_j", "WW_tensor"), envir=environment())
    res <- pblapply(1:m,function(j) get_Z_j(X[,j], W, C1, C2,
                                            Unico.mdl$"mus_hat"[j,], Unico.mdl$"gammas_hat"[j,], Unico.mdl$"betas_hat"[j,],
                                            Unico.mdl$"sigmas_hat"[j,,],  Unico.mdl$"taus_hat"[j], WW_tensor),
                    cl = cl)
    if (parallel) stop_cluster(cl)

    flog.info("Formating tensor result ...")
    for (j in 1:m){
        Z_hat[j,,] = res[[j]]
    }
    Z_hat = aperm(Z_hat, c(3,1,2))
    dimnames(Z_hat)[[1]] = source.ids
    dimnames(Z_hat)[[2]] = feature.ids
    dimnames(Z_hat)[[3]] = sample.ids

    # scale the data back to its original scale
    for (h in 1:k){
        Z_hat[h,,] = Z_hat[h,,] * repmat(Unico.mdl$scale.factor[feature.ids, , drop = F], 1, length(sample.ids))
    }
    flog.info("Finished tensor estimation")
    return(Z_hat)
}


#' @title Performs parametric statistical testing
#'
#' @description Performs parametric statistical testing (T-test) on (1) the marginal effect of each covariate in \code{C1} at source-specific level (2) the joint effect across all sources for each covariate in \code{C1} (3) non-source-specific effect for each covariate in \code{C2}. In the context of bulk genomic data containing a mixture of cell types, these correspond to the marginal effect of each covariate in \code{C1} (potentially including the phenotype of interest) at each cell type, joint tissue-level effect for each covariate in \code{C1}, and tissue-level effect for each covariate in \code{C2}.
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations. Note that \code{X} must be the same \code{X} used to learn \code{Unico.mdl} (i.e. the original observed 2D mixture used to fit the model).
#' @param Unico.mdl The entire set of model parameters estimated by Unico on the 2D mixture matrix (i.e. the list returned by applying function \code{Unico} to \code{X}).
#' @param slot_name A string indicating the key for storing the results under \code{Unico.mdl}
#' @param diag_only A logical value indicating whether to only use the estimated source-level variances (and thus ignoring the estimate covariance) for controlling the heterogeneity in the observed mixture. if set to FALSE, Unico instead estimates the observation- and feature-specific variance in the mixture by leveraging the entire \code{k} by \code{k} variance-covariance matrix.
#' @param intercept A logical value indicating whether to fit the intercept term when performing the statistical testing.
#' @param X_max_stds A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the observed mixture value. Only samples whose observed mixture value fall within \code{X_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.
#' @param Q_max_stds A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the estimated mixture variance. Only samples whose estimated mixture variance fall within \code{Q_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.
#' @param XQ_max_stds A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the weighted mixture value. Only samples whose weighted mixture value fall within \code{XQ_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; set \code{debug} to \code{TRUE} before reporting issues.
#'
#' @importFrom futile.logger flog.info
#' @importFrom pbapply pboptions pblapply
#' @importFrom parallel clusterExport
#' @importFrom stats anova lm
#' @importFrom pracma repmat
#'
#' @details If we assume that source-specific values \eqn{Z_{ijh}} are normally distributed, under the Unico model, we have the following:
#' \deqn{Z_{ij} \sim \mathcal{N}\left(\mu_{j} + (c_i^{(1)})^T \gamma_{jh}, \sigma_{jh}^2 \right)}
#' \deqn{X_{ij} \sim \mathcal{N}\left(w_{i}^T (\mu_{j} + (c_i^{(1)})^T \gamma_{jh}) + (c_i^{(2)})^T \beta_j,  \text{Sum}\left((w_i w_i^T ) \odot \Sigma_j\right) + \tau_j^2\right)}
#' For a given feature \eqn{j} under test, the above equation corresponds to a heteroskedastic regression problem with \eqn{X_{ij}} as the dependent variable and \eqn{\{\{w_i\}, \{w_i c_i^{(1)}\}, \{c_i^{(2)}\}\}} as the set of independent variables.
#' This view allows us to perform parametric statistical testing (T-test for marginal effects and partial F-test for joint effects) by solving a generalized least squares problem with sample \eqn{i} scaled by the inverse of its estimated standard deviation.
#'
#' @return An updated \code{Unico.mdl} object with the the following list of effect size and p-value estimates stored in an additional key specified by \code{slot_name}
#' \item{gammas_hat}{An \code{m} by \code{k*p1} matrix of the estimated effects of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} covariates on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
#' \item{betas_hat}{An \code{m} by \code{p2} matrix of the estimated effects of the \code{p2} covariates in \code{C2} on the mixture values of each of the \code{m} features in \code{X}.}
#' \item{gammas_hat_pvals}{An \code{m} by \code{k*p1} matrix of p-values for the estimates in \code{gammas_hat} (based on a T-test).}
#' \item{betas_hat_pvals}{An \code{m} by \code{p2} matrix of p-values for the estimates in \code{betas_hat} (based on a T-test).}
#' \item{gammas_hat_pvals.joint}{An \code{m} by \code{p1} matrix of p-values for the joint effects (i.e. across all \code{k} sources) of each of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X} (based on a partial F-test). In other words, these are p-values for the combined statistical effects (across all sources) of each one of the \code{p1} covariates on each of the \code{m} features under the Unico model.}
#' \item{Q}{An \code{m} by \code{n} matrix of weights used for controlling the heterogeneity of each observation at each feature (activated only if \code{debug == TRUE}).}
#' \item{masks}{An \code{m} by \code{n} matrix of logical values indicating whether observation participated in statistical testing at each feature  (activated only if \code{debug == TRUE}).}
#' \item{phi_hat}{An \code{m} by \code{k+p1*k+p2} matrix containing the entire estimated effect sizes (including those on source weights) for each feature (activated only if \code{debug == TRUE}).}
#' \item{phi_se}{An \code{m} by \code{k+p1*k+p2} matrix containing the estimated standard errors associated with \code{phi_hat} for each feature (activated only if \code{debug == TRUE}).}
#' \item{phi_hat_pvals}{An \code{m} by \code{k+p1*k+p2} matrix containing the p-values associated with \code{phi_hat} for each feature (activated only if \code{debug == TRUE}).}
#'
#' @examples
#' data = simulate_data(n=100, m=2, k=3, p1=1, p2=1, taus_std=0, log_file=NULL)
#' res = list()
#' res$params.hat = Unico(data$X, data$W, data$C1, data$C2, parallel=FALSE, log_file=NULL)
#' res$params.hat = association_parametric(data$X, res$params.hat, parallel=FALSE, log_file=NULL)
#'
#' @export association_parametric
association_parametric = function(X, Unico.mdl, slot_name = "parametric",
                                  diag_only = FALSE, intercept = TRUE,
                                  X_max_stds  = 2,
                                  Q_max_stds  = Inf,
                                  XQ_max_stds = Inf,
                                  parallel = TRUE, num_cores = NULL,
                                  log_file = "Unico.log", verbose = FALSE, debug = FALSE){

    start_logger(log_file, debug, verbose)
	  config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)

    flog.info("Validate parametric inputs ...")
    validate_inputs_parametric(X, Unico.mdl, slot_name,
                               diag_only, intercept,
                               X_max_stds, Q_max_stds, XQ_max_stds,
                               parallel, num_cores, log_file, debug, verbose)

    W = Unico.mdl$W
    n = nrow(W)
    m = nrow(X)
    k = ncol(W)
    source.ids  = colnames(W)
    feature.ids = rownames(X)
    sample.ids  = rownames(W)

    C1 = Unico.mdl$C1
    C2 = Unico.mdl$C2
    p1 = ncol(C1)
    p2 = ncol(C2)
    C1.ids = colnames(C1)
    C2.ids = colnames(C2)

    taus_hat   = Unico.mdl$taus_hat  [feature.ids,,drop = F]
    sigmas_hat = Unico.mdl$sigmas_hat[feature.ids,,]

    #transfer it to scaled version
    X          = X/Unico.mdl$scale.factor[rownames(X), ] #scale features so that ana calculation make sense
    X          = t(X)

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

    flog.info(paste0("Starting parametric pvals calculation: ", slot_name, " ..."))
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
        assert(length(res[[j]]) == 3*s_len + n + p1 )
        phi_hat[j,]   <- res[[j]][seq(0*s_len+1, 1*s_len,   length=s_len)]
        phi_se[j,]    <- res[[j]][seq(1*s_len+1, 2*s_len,   length=s_len)]
        phi_pvals[j,] <- res[[j]][seq(2*s_len+1, 3*s_len,   length=s_len)]
        masks[,j]     <- res[[j]][seq(3*s_len+1, 3*s_len+n, length=n)]
        gammas_hat_pvals.joint[j, ] <- res[[j]][seq(3*s_len+n+1, 3*s_len+n+p1, length=p1)]
    }

    param.res = list(gammas_hat       = phi_hat  [, seq(k+1,      k+p1*k,    length=p1*k)],
                     betas_hat        = phi_hat  [, seq(k+p1*k+1, k+p1*k+p2, length=p2)],
                     gammas_hat_pvals = phi_pvals[, seq(k+1,      k+p1*k,    length=p1*k)],
                     betas_hat_pvals  = phi_pvals[, seq(k+p1*k+1, k+p1*k+p2, length=p2)],

                     gammas_hat_pvals.joint = gammas_hat_pvals.joint)

    if(debug){
    	param.res = c(param.res, list(Q = t(Q),  masks = t(masks), # now features by samples
    																phi_hat = phi_hat, phi_se = phi_se, phi_hat_pvals = phi_pvals))
    }

    Unico.mdl[[slot_name]] = param.res
    flog.info(paste0("Finished parametric pvals calculation: ", slot_name))
    return(Unico.mdl)
}


#' @title Performs asymptotic statistical testing under no distribution assumption
#'
#' @description Performs asymptotic statistical testing on (1) the marginal effect of each covariate in \code{C1} at source-specific level (2) non-source-specific effect for each covariate in \code{C2}. In the context of bulk genomic data containing a mixture of cell types, these correspond to the marginal effect of each covariate in \code{C1} (potentially including the phenotype of interest) at each cell type and tissue-level effect for each covariate in \code{C2}.
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations. Note that \code{X} must be the same \code{X} used to learn \code{Unico.mdl} (i.e. the original observed 2D mixture used to fit the model).
#' @param Unico.mdl The entire set of model parameters estimated by Unico on the 2D mixture matrix (i.e. the list returned by applying function \code{Unico} to \code{X}).
#' @param slot_name A string indicating the key for storing the results under \code{Unico.mdl}
#' @param diag_only A logical value indicating whether to only use the estimated source-level variances (and thus ignoring the estimate covariance) for controlling the heterogeneity in the observed mixture. if set to FALSE, Unico instead estimates the observation- and feature-specific variance in the mixture by leveraging the entire \code{k} by \code{k} variance-covariance matrix.
#' @param intercept A logical value indicating whether to fit the intercept term when performing the statistical testing.
#' @param X_max_stds A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the observed mixture value. Only samples whose observed mixture value fall within \code{X_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.
#' @param Q_max_stds A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the estimated mixture variance. Only samples whose estimated mixture variance fall within \code{Q_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.
#' @param V_min_qlt A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the estimated moment condition variance. This value should be between 0 and 1. Only samples whose estimated moment condition variance fall outside the bottom \code{V_min_qlt} quantile will be used for the statistical testing of a given feature.
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; set \code{debug} to \code{TRUE} before reporting issues.
#'
#' @importFrom futile.logger flog.info
#' @importFrom pbapply pboptions pblapply
#' @importFrom parallel clusterExport
#' @importFrom stats pnorm quantile
#' @importFrom pracma inv repmat
#'
#' @details Under no distribution assumption, we can solve for the following weighted least square problem, which is similar to the heteroskedastic regression view described in \link{association_parametric}.
#' \deqn{\hat{\phi_j}^{\text{asym}} = \text{argmin}_{\phi_j} (x_j - S\phi_j) ^T Q_j (x_j - S\phi_j)}
#'
#' \eqn{S} denotes the design matrix formed by stacking samples in the rows and dependent variables \eqn{\{\{w_i\}, \{w_i c_i^{(1)}\}, \{c_i^{(2)}\}\}} on the columns.
#' \eqn{\phi_j} denotes the corresponding effect sizes on the dependent variables.
#' \eqn{Q_j} denotes the feature-specific weighting scheme. Similar to the parametric counterpart, \eqn{Q_j=\text{diag}(q_{1j}^2,...,q_{nj}^2)}, where for each sample \eqn{i}, its corresponding weight will be the inverse of the estimated variance in the mixture: \eqn{q_{ij}^2 = \frac{1}{sum(w_i w_i^T \odot \hat{\Sigma}_j)}}.
#' Marginal testing can thus be carried out on each dependent variable via the asymptotic distribution of the estimator \eqn{\hat{\phi_j}^{\text{asym}}}.
#'
#' @return An updated \code{Unico.mdl} object with the the following list of effect size and p-value estimates stored in an additional key specified by \code{slot_name}
#' \item{gammas_hat}{An \code{m} by \code{k*p1} matrix of the estimated effects of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} covariates on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
#' \item{betas_hat}{An \code{m} by \code{p2} matrix of the estimated effects of the \code{p2} covariates in \code{C2} on the mixture values of each of the \code{m} features in \code{X}.}
#' \item{gammas_hat_pvals}{An \code{m} by \code{k*p1} matrix of p-values for the estimates in \code{gammas_hat} (based on a T-test).}
#' \item{betas_hat_pvals}{An \code{m} by \code{p2} matrix of p-values for the estimates in \code{betas_hat} (based on a T-test).}
#' \item{Q}{An \code{m} by \code{n} matrix of weights used for controlling the heterogeneity of each observation at each feature (activated only if \code{debug == TRUE}).}
#' \item{masks}{An \code{m} by \code{n} matrix of logical values indicating whether observation participated in statistical testing at each feature  (activated only if \code{debug == TRUE}).}
#' \item{fphi_hat}{An \code{m} by \code{n} matrix containing the entire estimated moment condition variance for each feature. Note that observations who are considered as outliers due to any of the criteria will be marked as -1 in the estimated moment condition variance (activated only if \code{debug == TRUE}).}
#' \item{phi_hat}{An \code{m} by \code{k+p1*k+p2} matrix containing the entire estimated effect sizes (including those on source weights) for each feature (activated only if \code{debug == TRUE}).}
#' \item{phi_se}{An \code{m} by \code{k+p1*k+p2} matrix containing the estimated standard errors associated with \code{phi_hat} for each feature (activated only if \code{debug == TRUE}).}
#' \item{phi_hat_pvals}{An \code{m} by \code{k+p1*k+p2} matrix containing the p-values associated with \code{phi_hat} for each feature (activated only if \code{debug == TRUE}).}
#'
#' @examples
#' data = simulate_data(n=100, m=2, k=3, p1=1, p2=1, taus_std=0, log_file=NULL)
#' res = list()
#' res$params.hat = Unico(data$X, data$W, data$C1, data$C2, parallel=FALSE, log_file=NULL)
#' res$params.hat = association_asymptotic(data$X, res$params.hat, parallel=FALSE, log_file=NULL)
#'
#' @export association_asymptotic
association_asymptotic = function(X, Unico.mdl, slot_name = "asymptotic",
                                  diag_only = FALSE, intercept = TRUE,
                                  X_max_stds  = 2,    #bulk profile
                                  Q_max_stds  = Inf,  #analytical variance
                                  V_min_qlt   = 0.05, #moment condition variance
                                  parallel = TRUE, num_cores = NULL,
                                  log_file = "Unico.log", verbose = FALSE, debug = FALSE){

    start_logger(log_file, debug, verbose)
	  config = config::get(file = system.file("extdata", "config.yml", package = "Unico"), use_parent = FALSE)

    flog.info("Validate asymptotic inputs ...")
    validate_inputs_asymptotic(X, Unico.mdl, slot_name,
                               diag_only, intercept,
                               X_max_stds, Q_max_stds, V_min_qlt,
                               parallel, num_cores, log_file, debug, verbose)

    #input format
    W = Unico.mdl$W
    n = nrow(W)
    m = nrow(X)
    k = ncol(W)
    source.ids  = colnames(W)
    feature.ids = rownames(X)
    sample.ids  = rownames(W)

    C1 = Unico.mdl$C1
    C2 = Unico.mdl$C2
    p1 = ncol(C1)
    p2 = ncol(C2)
    C1.ids = colnames(C1)
    C2.ids = colnames(C2)

    mus_hat    = Unico.mdl$mus_hat   [feature.ids,]
    gammas_hat = Unico.mdl$gammas_hat[feature.ids,,drop = F]
    betas_hat  = Unico.mdl$betas_hat [feature.ids,,drop = F]
    taus_hat   = Unico.mdl$taus_hat  [feature.ids,,drop = F]
    sigmas_hat = Unico.mdl$sigmas_hat[feature.ids,,]

    #transfer it to scaled version
    X = X/Unico.mdl$scale.factor[rownames(X), ] #scale features so that ana calculation make sense
    X = t(X) #X should be samples by features

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

    flog.info(paste0("Starting asymptotic pvals calculation: ", slot_name, " ..."))
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
        S    =  S[mask, ]
        X_j  =  X_j[mask]
        Q_j  =  Q_j[mask]
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

    }, cl = cl)

    if (parallel) stop_cluster(cl)
    for (j in 1:m){
        phi_hat[j,]   <- res[[j]][seq(0*s_len+1, 1*s_len,  length=s_len)]
        phi_se[j,]    <- res[[j]][seq(1*s_len+1, 2*s_len,  length=s_len)]
        phi_pvals[j,] <- res[[j]][seq(2*s_len+1, 3*s_len,  length=s_len)]
        masks[,j]     <- res[[j]][seq(3*s_len+1, 3*s_len+n,length=n)]
        fphi_var[which(masks[,j] == 1),j] <- res[[j]][seq(3*s_len+n+1, length(res[[j]]), length=sum(masks[,j]))]
        assert(length(res[[j]]) == 3*s_len+n+sum(masks[,j]))
    }

    asym.res = list(gammas_hat       = phi_hat  [, seq(k+1,      k+p1*k,    length=p1*k)],
                    betas_hat        = phi_hat  [, seq(k+p1*k+1, k+p1*k+p2, length=p2)],
                    gammas_hat_pvals = phi_pvals[, seq(k+1,      k+p1*k,    length=p1*k)],
                    betas_hat_pvals  = phi_pvals[, seq(k+p1*k+1, k+p1*k+p2, length=p2)])

    if(debug){
    	asym.res = c(asym.res, list(Q = t(Q), masks  = t(masks), fphi_var = t(fphi_var), # now features by samples,
    															phi_hat = phi_hat, phi_se = phi_se, phi_hat_pvals = phi_pvals))
    }

    Unico.mdl[[slot_name]] = asym.res
    flog.info(paste0("Finished asymptotic pvals calculation: ", slot_name))
    return(Unico.mdl)
}



#' @title Simulate data under Unico model assumption
#'
#' @description Simulate all model parameters and sample source specific data from multivariate gaussian with full covariance structure.
#'
#' @param n A positive integer indicating the number of observations to simulate.
#' @param m A positive integer indicating the number of features to simulate.
#' @param k A positive integer indicating the number of sources to simulate.
#' @param p1 A non-negative integer indicating the number of source-specific covariates to simulate.
#' @param p2 A non-negative indicating the number of non-source-specific covariates to simulate.
#' @param mus_mean A numerical value indicating the average of the source specific means.
#' @param mus_std A positive value indicating the variation of the source specific means across difference sources.
#' @param gammas_mean A numerical value indicating the average effect sizes of the source-specific covariates.
#' @param gammas_std A non-negative numerical value indicating the variation of the effect sizes of the source-specific covariates.
#' @param betas_mean A numerical value indicating the average effect sizes of the non-source-specific covariates.
#' @param betas_std A non-negative numerical value indicating the variation of the effect sizes of the non-source-specific covariates.
#' @param sigmas_lb A numerical value indicating the lower bound of a uniform distribution from which we sample entries of matrix \code{A} used to construct the feature specific \code{k} by \code{k} variance-covariance matrix.
#' @param sigmas_ub A numerical value indicating the upper bound of a uniform distribution from which we sample entries of matrix \code{A} used to construct the feature specific \code{k} by \code{k} variance-covariance matrix.
#' @param taus_std non-negative numerical value indicating the variation of the measurement noise across difference features.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @importFrom futile.logger flog.info
#' @importFrom compositions rDirichlet.acomp
#' @importFrom stats rnorm runif
#' @importFrom pracma repmat
#' @importFrom MASS mvrnorm
#'
#' @details Simulate data based on the generative model described in function \link{Unico}.
#'
#' @return A list of simulated model parameters, covariates, observed mixture, and source-specific data.
#' \item{X}{An \code{m} by \code{n} matrix of the simulated mixture for \code{m} features and \code{n} observations.}
#' \item{W}{An \code{n} by \code{k} matrix of the weights/proportions of \code{k} source for each of the \code{n} observations.}
#' \item{C1}{An \code{n} by \code{p1} matrix of the simulated covariates that affect the source-specific values.}
#' \item{C2}{An \code{n} by \code{p2} matrix of the simulated covariates that affect the mixture values.}
#' \item{Z}{A \code{k} by \code{m} by \code{n} tensor of the source specific values for each of the \code{k} sources}
#' \item{mus}{An \code{m} by \code{k} matrix of the mean of each of the \code{m} features for each of the \code{k} sources.}
#' \item{gammas}{An \code{m} by \code{k*p1} matrix of the effect sizes of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} covariates on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
#' \item{betas}{An \code{m} by \code{p2} matrix of the effect sizes of the \code{p2} covariates in \code{C2} on the mixture values of each of the \code{m} features.}
#' \item{sigmas}{An \code{m} by \code{k} by \code{k} tensor of the variance-covariance matrix of each of the \code{m} features.}
#' \item{taus}{An \code{m} by \code{1} matrix of the feature specific variance of the measurement noise for all m features.}

#'
#' @examples
#' data = simulate_data(n=100, m=2, k=3, p1=1, p2=1, taus_std=0, log_file=NULL)
#'
#' @export simulate_data

simulate_data = function(n, m, k, p1, p2,
												 mus_mean = 10, mus_std = 2,
												 gammas_mean = 1, gammas_std = 0.1,
												 betas_mean = 1, betas_std = 0.1,
												 sigmas_lb = 0, sigmas_ub = 1, taus_std = 0.1,
												 log_file = "Unico.log", verbose = FALSE){
	start_logger(log_file, FALSE, verbose)
	flog.info("Start simulation ...")

	sample_ids  <- vapply(1:n, function(i) paste0("sample.", i),  character(1))
	feature_ids <- vapply(1:m, function(j) paste0("feature.", j), character(1))
	source_ids  <- vapply(1:k, function(h) paste0("source.", h),  character(1))


	W.alpha = sort(runif(k, min = 0, max = 20), decreasing = T)
	W = matrix(rDirichlet.acomp(n,W.alpha), n, k)
	rownames(W) <- sample_ids
	colnames(W) <- source_ids

	# mean related
	mus    = abs(matrix(rnorm(m * k, mus_mean, mus_std), m, k))
	rownames(mus) <- feature_ids
	colnames(mus) <- source_ids

	C1 <- matrix(rnorm(n*p1, mean = 0, sd = 1), n, p1)
	gammas = array(rnorm(m * k * p1, gammas_mean, gammas_std), c(m, k, p1))
	if (p1){
		C1_names  <- vapply(1:p1, function(p) paste0("C1.", p), character(1))
		C1_names_ <- c()
		for (j in 1:k){
			C1_names_ = c(C1_names_, unlist(lapply(C1_names,function(i) paste(source_ids[j],"-",i,sep=""))))
		}
		rownames(C1) <- sample_ids
		colnames(C1) <- C1_names
	}

	C2 <- matrix(rnorm(n*p2, mean = 0, sd = 1), n, p2)
	betas  = matrix(rnorm(m * p2, betas_mean, betas_std), m, p2)
	if (p2){
		C2_names <- vapply(1:p2, function(p) paste0("C2.", p),  character(1))
		rownames(betas) <- feature_ids
		colnames(betas) <- C2_names
		rownames(C2) <- sample_ids
		colnames(C2) <- C2_names
	}

	# var related
	sigmas = array(0, c(m, k, k))
	for (j in 1:m){
		A = matrix(runif(k**2, sigmas_lb, sigmas_ub), k, k)
		sigmas[j,,] = A %*% t(A)
		assert_cov(sigmas[j,,])
	}
	dimnames(sigmas)[[1]] <- feature_ids
	dimnames(sigmas)[[2]] <- source_ids
	dimnames(sigmas)[[3]] <- source_ids

	taus = matrix(abs(rnorm(m, 0, taus_std)), m, 1)
	rownames(taus) <- feature_ids
	colnames(taus) <- c("tau.square")

	Z = array(0, c(k, m, n))
	for (j in 1: m){
		sigma = sigmas[j,,]
		mean_j = pracma::repmat(mus[j,], n, 1) + C1 %*% t(gammas[j,,])
		Z[,j,] = t(mean_j + MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = sigma))
	}
	dimnames(Z)[[1]] = source_ids
	dimnames(Z)[[2]] = feature_ids
	dimnames(Z)[[3]] = sample_ids

	gammas = t(vapply(1:m, function(j) as.vector(t(gammas[j,,])), numeric(k * p1))) #flatten by rows -> m by (k*p1)
	if (p1){
		colnames(gammas) <- C1_names_
		rownames(gammas) <- feature_ids
	}

	# generative model on global signals
	X = matrix(0, n, m)
	for (h in 1:k){
		X = X + repmat(as.matrix(W[,h]),1, m) * t(Z[h,,])
	}
	X = X + C2 %*% t(betas)

	# measurement noise
	errors = matrix(0, n, m)
	for (j in 1: m){
		errors[,j] = rnorm(n, mean = 0, sd = sqrt(taus[j]))
	}
	X = X + errors
	rownames(X) <- sample_ids
	colnames(X) <- feature_ids
	X = t(X)

	sim.data = list(X = X, W = W, C1 = C1, C2 = C2, Z = Z,
									mus = mus, gammas = gammas, betas = betas,
									sigmas = sigmas, taus = taus)
	flog.info(paste0("Finished simulation"))
	return(sim.data)
}
