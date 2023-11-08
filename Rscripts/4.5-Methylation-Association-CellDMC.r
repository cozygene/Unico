# Depends on the input this notebook fits a CellDMC model with appropriate C1, C2 and performs association testing on both "age" and "gender" 
# it takes the following arguments
# data.version:  choose from "liu" "hannum" "hannon1" "hannon2"
# res.dir:      character, which directory to save the result 
# additional parameters: specify the version of association test to perform 

# example 
# Rscript 4.5-Methylation-Association-CellDMC.r liu /u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/


require("data.table")
library("stringr")
library("testit")

source("Association.CellDMC.utils.r")
set.seed(2023)

args = commandArgs(trailingOnly=TRUE)
data.version  = as.character(args[1])
res.dir       = as.character(args[2])

dir.version = "XY"
study.versions = c("age", "gender")
assoc_versions = list("parametric")

print(paste0("dir.version: ", dir.version))
print(paste0("data.version: ", data.version))
print(paste0("study.versions: ", study.versions))
print(paste0("res.dir: ", res.dir))
print(paste0("assoc_versions: ", assoc_versions))

data.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/Data/Methylation/Consistency/"

file.paths = list()
file.paths[["liu"]]        = file.path(data.dir, "liu.processed.RData")
file.paths[["hannum"]]     = file.path(data.dir, "hannum.processed.RData")
file.paths[["hannon1"]]    = file.path(data.dir, "hannon1.processed.RData")
file.paths[["hannon2"]]    = file.path(data.dir, "hannon2.processed.RData")

res.dir  = file.path(res.dir, dir.version, data.version)
res.file = file.path(res.dir, paste0("CellDMC.mdl.rds"))
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

#reorder W by aboundance 
W = W[,order(-colMeans(W))]


start.t = Sys.time()
CellDMC.mdl = list()
for (study.version in study.versions){# age and gender
    
    print(study.version)
    study.res = list()
    y       = C1[, study.version]
    cov.mod = cbind(C1[, colnames(C1)!=study.version], C2)
    
    print(head(y))
    print(cov.mod[1:min(nrow(cov.mod), 5), 1:min(nrow(cov.mod), 5)])
    
    
    if ("parametric" %in% assoc_versions){                     
        print(paste0("fitting ", study.version, " parametric")) 
        study.res[["parametric"]]  = celldmc_w_joint(X = X, y = y, W = W, cov.mod = cov.mod)
    } 
    
    CellDMC.mdl[[study.version]] = study.res
} 
end.t = Sys.time()
print(start.t - end.t)

saveRDS(CellDMC.mdl, res.file)
