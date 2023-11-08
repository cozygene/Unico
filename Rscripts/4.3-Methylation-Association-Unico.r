# Depends on the input this notebook fits a Unico model with appropriate C1, C2 and performs association testing 
# it takes the following arguments
# dir.version: choose from "XY" or "YX" only relavent if fit_assoc is turned on
#              "XY" indicates that methylation is modeled as the response variable, phenotype as explanatory variable 
#              "YX" indicates that phenotype is modeled as the response variable, methylation as explanatory variable

# data.version:  choose from "liu" "hannum" "hannon1" "hannon2"
# study.version: choose from covariates namesa and "NA"
#                when fitting the "XY" direction, explicitly set it to NA so that all covars are fitted together
#                if set it to a specfic covar, we will test the pval calibration by shuffling that covar's sample assignment
#                when fitting the "YX" direction, explcitly set it to the covar name under test. 
# model_tau:     model_tau (boolean) used in the Unico model fitting, only relavent if fit_model is turned on

# fit_model:    boolean, whether to fit the model or not 
# fit_assoc     boolean, whether to fit association or not
# fit_tensor    boolean, whether to perform tensor operation or not

# res.dir:      character, which directory to save the result 
# additional parameters: specify the version of association test to perform 

# example 
# Rscript 4.3-Methylation-Association-Unico.r XY liu NA FALSE TRUE TRUE FALSE /u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/
# Rscript 4.3-Methylation-Association-Unico.r YX hannon1 disease FALSE TRUE TRUE TRUE /u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/

source("Unico.r")
source("analysis.utils.r")

set.seed(2023)

args = commandArgs(trailingOnly=TRUE)

dir.version   = as.character(args[1]) 
data.version  = as.character(args[2])
study.version = as.character(args[3])
model_tau     = as.logical(args[4])
fit_model     = as.logical(args[5])
fit_assoc     = as.logical(args[6])
fit_tensor    = as.logical(args[7])
res.dir       = as.character(args[8])

if (length(args) >=9){
    assoc_versions = as.list(args[9:length(args)])
}else{
    assoc_versions = list("parametric")
}


print(paste0("dir.version: ", dir.version))
print(paste0("data.version: ", data.version))
print(paste0("study.version: ", study.version))
print(paste0("model_tau: ", model_tau))
print(paste0("fit_model: ", fit_model))
print(paste0("fit_assoc: ", fit_assoc))
print(paste0("fit_tensor: ", fit_tensor))
print(paste0("res.dir: ", res.dir))
print(paste0("assoc_versions: ", assoc_versions))

data.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/Data/Methylation/Consistency/"

file.paths = list()
file.paths[["liu"]]        = file.path(data.dir, "liu.processed.RData")
file.paths[["hannum"]]     = file.path(data.dir, "hannum.processed.RData")
file.paths[["hannon1"]]    = file.path(data.dir, "hannon1.processed.RData")
file.paths[["hannon2"]]    = file.path(data.dir, "hannon2.processed.RData")

if (dir.version == "XY" & study.version == "NA"){
    res.dir = file.path(res.dir, dir.version, data.version)
}else if(dir.version == "XY" & study.version != "NA"){
    print("################################### WARNING #################################")
    print(paste0("in XY direction with study version provided, running the shuffled version on: ", study.version))
    res.dir = file.path(res.dir, dir.version, paste0(data.version, "-", study.version,  "-shuffled"))
    print("################################### WARNING #################################")
}else if ((dir.version == "YX")){
    res.dir = file.path(res.dir, dir.version, paste0(data.version, "-", study.version))
}else{
    print("wrong dir.version")
}

res.file = file.path(res.dir, paste0("Unico.mdl.rds"))
log.file = file.path(res.dir, paste0("Unico.mdl.log"))

if (!file.exists(res.dir)){
    dir.create(res.dir, recursive = T)
}

load(file.paths[[data.version]])

# by default: for XY direction all biological covars (including study.version y) are included as C1
if (data.version == "liu"){
    X  = liu$X; 
    W  = liu$W; 
    C1 = liu$cov[, c("age", "gender", "disease","smoking")];
    C2 = liu$ctrl_pcs;
}else if(data.version == "hannum"){
    X  = hannum$X; 
    W  = hannum$W; 
    C1 = hannum$cov[, c("age", "gender", "ethnicity", "smoking")];
    C2 = cbind(hannum$ctrl_pcs, hannum$cov[,"plate", drop = F])

}else if(data.version == "hannon1"){
    X  = hannon1$X; 
    W  = hannon1$W; 
    C1 = hannon1$cov[, c("age", "gender", "disease")];
    C2 = hannon1$ctrl_pcs

}else if(data.version == "hannon2"){
    X  = hannon2$X; 
    W  = hannon2$W; 
    C1 = hannon2$cov[, c("age", "gender", "disease")];
    C2 = hannon2$ctrl_pcs
}else{
    print("check your input")    
}


if (dir.version == "YX"){
    # C1 excludes y so each get a different directory
    y = C1[,study.version, drop = F]
    # TODO add conversion to binary if needed 
    C1 = C1[, setdiff(colnames(C1), c(study.version))]
    
    if(study.version == "disease"){
        message("make the control 0 and case 1")
        y = y - 1
        table(y)
    }
    if(study.version == "smoking"){
        message("converting never smoke 0 and ever smoked 1")
        y[y>0,] = 1
        table(y)
    }
    if(study.version == "gender"){
        message("converting gender to 0 and 1")
        y = y - 1
        table(y)
    }
}

if (dir.version == "XY" & study.version != "NA"){
    print(paste0("shuffeling the pehnotype of interest:", study.version))
    print(paste0("before shuffeling "))
    print(C1[1:10,])
    C1[,study.version] = sample(as.numeric(C1[,study.version]))
    print(paste0("after shuffeling "))
    print(C1[1:10,])
}

#reorder W by aboundance 
W = W[,order(-colMeans(W))]


if (fit_model){
    Unico.mdl = list()
    Unico.mdl$params.hat <- Unico(X = as.matrix(X), W = as.matrix(W), C1 = as.matrix(C1), C2 = as.matrix(C2), 
                                  parallel = TRUE, num_cores = 12, 
                                  log_file = log.file, verbose = TRUE)
                                
    #capping at interal scale
    Unico.mdl$params.hat$sigmas_hat  = cap_values(Unico.mdl$params.hat$sigmas_hat,  max.val = 10 ** (4), min.val = 10**(-4))
    
    #hat original input scale 
    Unico.mdl$params.hat.orig = list(mus_hat = scale_feature_matrix(Unico.mdl$params.hat$mus_hat, 
                                                               Unico.mdl$params.hat$scale.factor),
                                     sigmas_hat = scale_feature_source_source(Unico.mdl$params.hat$sigmas_hat, 
                                                                         Unico.mdl$params.hat$scale.factor),
                                     taus  = scale_feature_matrix(Unico.mdl$params.hat$taus_hat, 
                                                                 Unico.mdl$params.hat$scale.factor))

    #save the object with estiamted parameters
    saveRDS(Unico.mdl, res.file)
}


if (dir.version == "XY" & fit_assoc){
    message("running association")
    Unico.mdl = readRDS(res.file)    
                              
    if ("parametric" %in% assoc_versions){                   
        Unico.mdl$params.hat = add_C1_C2_pvals_parametric(X = X, Unico.mdl = Unico.mdl$params.hat, slot_name = "parametric", 
                                                          parallel = F, num_cores = NULL, 
                                                          config_file = NULL, log_file = NULL, debug = FALSE, verbose = TRUE)

        saveRDS(Unico.mdl, res.file)}
    
    if ("asymptotic" %in% assoc_versions){   
    Unico.mdl$params.hat = add_C1_C2_pvals_asymptotic(X = X, Unico.mdl = Unico.mdl$params.hat, slot_name = "asymptotic", 
                                                      parallel = F, num_cores = NULL, 
                                                      config_file = NULL, log_file = NULL, debug = FALSE, verbose = TRUE)

    saveRDS(Unico.mdl, res.file)} 
    
}