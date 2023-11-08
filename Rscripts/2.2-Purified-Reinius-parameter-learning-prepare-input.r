library(data.table)
library(matrixStats)
source("analysis.utils.r")
set.seed(2023)

data.dir = "../Data/Methylation/Purified-Reinius/"

#read in the original input 
Reinius = readRDS(file.path(data.dir, "reinius.rds"))
Hannum  = readRDS(file.path(data.dir, "hannum.rds"))

str(Hannum)

str(Reinius)

feature.ids = intersect(rownames(Hannum$X), rownames(Reinius$X))
source.ids  = colnames(Hannum$W)

Reinius$X = Reinius$X [feature.ids,]
Reinius$W = Reinius$W [,source.ids]
Reinius$Z = Reinius$Z [source.ids,feature.ids,]
Reinius$params = calc_params_from_Z(Reinius$Z, max_sds = Inf) # only 6 samples dont remove any sample

#keep Caucasian samples only
keep.samples    = rownames(Hannum$cov[Hannum$cov[,"ethnicity"] == "Caucasian - European",])
Hannum$X        = Hannum$X[feature.ids,keep.samples]
Hannum$cov      = Hannum$cov[keep.samples, ]
Hannum$W        = Hannum$W[keep.samples, source.ids]
Hannum$ctrl_pcs = Hannum$ctrl_pcs[keep.samples, ]

str(Hannum)

str(Reinius)

source.ids = colnames(Reinius$W)
sum.var = matrix(0, length(feature.ids), 1)
for (source.id in source.ids){
    sum.var = sum.var + rowVars(Reinius$Z[source.id,,])
}
rownames(sum.var)  = feature.ids

sum.var[order(-sum.var[,1]),,drop = F]

feature.ids[order(-sum.var)[1:10000]]

for (version in c("hvf.10k", "random.10k")){
    print(paste0("working on ", version))
    
    if (version == "hvf.10k"){
        keep.cpgs = feature.ids[order(-sum.var)[1:10000]]
    }else if(version == "random.10k"){
        keep.cpgs = sample(feature.ids, 10000)        
    }else{
        print("wrong version")
    }
    
 
    #Reinius
    Reinius.sub   = copy(Reinius)
    Reinius.sub$X = Reinius.sub$X [keep.cpgs,]
    Reinius.sub$Z = Reinius.sub$Z [source.ids,keep.cpgs,]
    Reinius.sub$params = calc_params_from_Z(Reinius.sub$Z, max_sds = Inf) # only 6 samples dont remove any sample
    
    #Hannum
    Hannum.sub   = copy(Hannum)
    Hannum.sub$X = Hannum.sub$X[keep.cpgs, ]
    
    #save subset of data
    Hannum.sub.file   = file.path(data.dir, version, paste0("hannum.", version, ".rds"))
    Reinius.sub.file  = file.path(data.dir, version, paste0("reinius.", version, ".rds"))
    
    saveRDS(Hannum.sub,   Hannum.sub.file) 
    saveRDS(Reinius.sub,  Reinius.sub.file)

    #save cibersortx related data 
    cibersortx.X.file = file.path(data.dir, version, paste0("cibersortx.", version, ".X.txt"))
    cibersortx.W.file = file.path(data.dir, version, paste0("cibersortx.", version, ".W.txt"))
    #so that cibersortx will not consider this as log transformed version of the expression count
    cibersortx.X = cbind(Reinius.sub$X, Hannum.sub$X) * 10000 
    cibersortx.W = rbind(Reinius.sub$W, Hannum.sub$W) 
    fwrite(as.data.frame(cibersortx.X),   
                 file = cibersortx.X.file,  
                 sep = "\t", quote=FALSE, row.names = T, col.names = T)

    fwrite(as.data.frame(cibersortx.W),   
                 file = cibersortx.W.file,  
                 sep = "\t", quote=FALSE, row.names = T, col.names = T)
   
}

str(Hannum.sub)

str(Reinius.sub)
