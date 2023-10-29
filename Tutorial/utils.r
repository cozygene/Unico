library(futile.logger)
library(parallel)
library(pbapply)

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
  flog.info("Initiate cluster...")
  # check if user asked too many cores
  max_cores = get_num_cores(NULL)
  if (!is.null(num_cores)){
      if(num_cores > max_cores){
          num_cores = max_cores
          flog.info("num_cores is set to maximum number of cores detected: %s .", max_cores)
      }   
  }
    
  cl <- makeCluster(get_num_cores(num_cores))
  flog.info("Parallel is on with %s nodes.",get_num_cores(num_cores))
  invisible(clusterEvalQ(cl, c(library("pracma"), library("mgcv"), library("nloptr"), library("testit"))))
  flog.info("Packages were loaded into the cluster nodes.")
  return(cl)
}
stop_cluster <- function(cl){
  flog.info("Stop cluster")
  stopCluster(cl)
}