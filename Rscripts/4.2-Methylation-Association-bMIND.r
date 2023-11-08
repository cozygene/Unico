#TODO
library(MIND)
set.seed(2023)

args = commandArgs(trailingOnly=TRUE)
data.version  = as.character(args[1])
res.dir       = as.character(args[2])

# data.version = "liu" # choose from "liu" "hannum" "hannon1" "hannon2"
# res.dir = "/u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/Debug"

dir.version = "XY"
pheno = "gender"

print(paste0("dir.version: ", dir.version))
print(paste0("data.version: ", data.version))
print(paste0("res.dir: ", res.dir))

data.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/Data/Methylation/Consistency/"

file.paths = list()
file.paths[["liu"]]        = file.path(data.dir, "liu.processed.RData")
file.paths[["hannum"]]     = file.path(data.dir, "hannum.processed.RData")
file.paths[["hannon1"]]    = file.path(data.dir, "hannon1.processed.RData")
file.paths[["hannon2"]]    = file.path(data.dir, "hannon2.processed.RData")

res.dir  = file.path(res.dir, dir.version, data.version)
res.file = file.path(res.dir, paste0("bMIND.mdl.rds"))
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




# pval for one gene in case some cell types have missing pval
get_pval = function(pval, cell_type, K) {
  pval0 = rep(NA, K)
  names(pval0) = cell_type
  names = intersect(names(pval), cell_type)
  pval0[names] = pval[names]
  return(pval0)
}


test = function(A, y, covariate = NULL) {
  
  if(dim(A)[3] != length(y)) print('CTS estimates and y have different length')
  if(!is.null(covariate)) if(dim(A)[3] != nrow(covariate)) print('CTS estimates and covariate have different number of samples/subjects') else {
    if(!is.null(rownames(covariate)) & any(rownames(covariate) != dimnames(A)[[3]])) covariate = covariate[dimnames(A)[[3]],]
  }
  
  K = ncol(A)
  cell_type = colnames(A)
  if(is.null(covariate)) pval = apply(A, 1, function(x) {
    pval = coef(summary(glm(y ~ ., data = data.frame(t(x)), family = 'binomial')))[,4]
    return(get_pval(pval, cell_type, K))
  }) else
    pval = apply(A, 1, function(x) {
      pval = coef(summary(glm(y ~ ., data = data.frame(t(x), covariate), family = 'binomial')))[,4]
      return(get_pval(pval, cell_type, K))
    })
  
  qval = pval2qval(pval, A, y, covariate)
  # rownames(qval) = rownames(pval) = substring(rownames(pval), 5)
  return(list(manova.pvals = qval, glm.pvals = pval))
}


# MANOVA; pval: K x ngene
pval2qval = function(pval, A, y, covariate = NULL) {
    ng = nrow(A)
    # pval for each gene
    if(is.null(covariate)){
        pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,,]) ~ y))$stats[1, "Pr(>F)"], silent = T))   
    }else{
        pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,,]) ~ y + covariate))$stats[1, "Pr(>F)"], silent = T))
    }
    
    #pval1 is m by 1 vector from the MANOVA
    #this step is filtering pval on features that has valid result from the MANOVA
    pval = pval[,!is.na(as.numeric(pval1))]

    # qval1 is still working on MANOVA stats
    pval1 = na.omit(as.numeric(pval1))
    qval1 = p.adjust(pval1, 'fdr')

    # the final qval result. first just copy the shape on the pval: K by M
    qval.adjust   = pval
    qval.unadjust = pval
    K = ncol(A)

    #per feature
    for(i in 1:ncol(pval)) {
        #initialize all celltypes to 1
        qval.adjust[,i] = 1
        qval.unadjust[,i] = 1
        #if the pval is decently small by 0.05/K. copy over the on MANOVA stat for the gene to the final qval. put it in the celltype with smallest p
        if(min(pval[,i], na.rm = T) < .05/K) qval.adjust[,i][which.min(pval[,i])] = qval1[i]
        if(min(pval[,i], na.rm = T) < .05/K) qval.unadjust[,i][which.min(pval[,i])] = pval1[i]
    }
                       
    joint.pvals = as.matrix(pval1)
    rownames(joint.pvals) = colnames(qval.adjust)
    colnames(joint.pvals) = c("joint")
                       
    return(list(joint.unadjust.pvals    = joint.pvals, 
                marginal.adjust.pvals   = qval.adjust, 
                marginal.unadjust.pvals = qval.unadjust))
}
                    


start.t = Sys.time()
# CTS-DE (np = TRUE: use non-informative prior)
bMIND.mdl = bMIND(X, W)
end.t = Sys.time()
print(end.t - start.t)

y =  as.numeric(as.factor(C1[,pheno]))-1 # used to be 1,2 coded, subtract 1 to be 0,1 coded
covariate = cbind(C1[,-which(colnames(C1) == pheno)], C2)
bMIND.mdl =  c(bMIND.mdl, test(bMIND.mdl$A, y, covariate))

saveRDS(bMIND.mdl, res.file)
