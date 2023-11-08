# Depends on the input this notebook fits a CellDMC model with appropriate C1, C2 and performs association testing on both "age" and "gender" 
# it takes the following arguments
# data.version:  choose from "liu" "hannum" "hannon1" "hannon2"
# res.dir:      character, which directory to save the result 
# additional parameters: specify the version of association test to perform 

# example 
# Rscript 4.6-Methylation-Association-Base.r liu /u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/
source("Association.Baseline.utils.r") 
set.seed(2023)

args = commandArgs(trailingOnly=TRUE)
data.version  = as.character(args[1])
res.dir       = as.character(args[2])

dir.version = "XY"
assoc_versions = list("parametric")

print(paste0("dir.version: ", dir.version))
print(paste0("data.version: ", data.version))
print(paste0("res.dir: ", res.dir))
print(paste0("assoc_versions: ", assoc_versions))

data.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/Data/Methylation/Consistency/"

file.paths = list()
file.paths[["liu"]]        = file.path(data.dir, "liu.processed.RData")
file.paths[["hannum"]]     = file.path(data.dir, "hannum.processed.RData")
file.paths[["hannon1"]]    = file.path(data.dir, "hannon1.processed.RData")
file.paths[["hannon2"]]    = file.path(data.dir, "hannon2.processed.RData")

res.dir  = file.path(res.dir, dir.version, data.version)
res.file = file.path(res.dir, paste0("base.mdl.rds"))
if (!file.exists(res.dir)){
    dir.create(res.dir, recursive = T)
}

res.dir

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
base.mdl = list()
   
C = cbind(C1, C2)

         
if ("parametric" %in% assoc_versions){                     
    print(paste0("fitting parametric")) 
    base.mdl[["parametric"]]  = base_assos(X = X, W, C)
} 

end.t = Sys.time()
print(end.t - start.t)

saveRDS(base.mdl, res.file)

