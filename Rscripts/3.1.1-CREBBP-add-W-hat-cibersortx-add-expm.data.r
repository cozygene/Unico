# This notebook add the CIBERSORTx estimated celltype proportion (run on local computer) to the CREBBP.dat.rds
# it then sets up the (X,W,C1) pair needed for the experiment, saved as entries in CREBBP.dat$expm.dat
# txt version of each (X,W,C1) pair is also exported so that CIBERSORTx can read and compute tensor in local computer
library("data.table")
set.seed(2023)

data.dir   = file.path("../Data/RNA/CREBBP/")
res.dir    = file.path("../Result/RNA/CREBBP/")

CREBBP.dat = readRDS(file.path(data.dir, "CREBBP.dat.rds"))
#all X have the same order of samples
#all DE genes or probes are unique and present in the corresponding X

W.hat = as.matrix(data.frame(fread(file.path(res.dir, "Cibersortx-W-hat", "CIBERSORTxGEP_Fractions.txt"), header = T), row.names = 1))
W.hat = W.hat[, 1:22]
# colapse to 4 celltypes
merge.class = colnames(fread(file.path(data.dir,"Cibersortx-fig3/Fig3b-f-LM4_merged_classes.txt")))
colnames(W.hat) = merge.class
W.hat = t(rowsum(t(W.hat), group = colnames(W.hat), na.rm = T))
# remove space in celltype names
colnames(W.hat) = gsub(" ", "", colnames(W.hat))
# reorder to be most abundant to least abundant
W.hat = W.hat[,order(-colMeans(W.hat))]
# reorder sample.ids
W.hat = W.hat[colnames(CREBBP.dat$X.genes),]

W.hat

CREBBP.dat$W.hat = W.hat 

# logistic stuff to save cibsersortx need data
# key: character: used to identify the experiment set up
# expm: a list containing all necessary (W,X,C1) information
# ciberosrtx.dir: parent directory that the txt files of this experiment will be saved
dumpCibersortx = function(key, expm, cibersortx.dir){
  
  fp.src = file.path(cibersortx.dir, key, "src")
  fp.res = file.path(cibersortx.dir, key, "res")
  
  if (!file.exists(fp.src)){dir.create(fp.src, recursive = T)}
  if (!file.exists(fp.res)){dir.create(fp.res, recursive = T)}
  print(paste0("cibersortx data saved here: ",fp.src))
  print(paste0("cibersortx res dir: ",fp.res))
  
  #cibersortx data dump 
  fwrite(as.data.frame(expm$X),   
         file = file.path(fp.src,"X.txt"),  
         sep = "\t", quote=FALSE, row.names = T, col.names = T)
  
  fwrite(as.data.frame(expm$W),   
         file = file.path(fp.res,"W.txt"),  
         sep = "\t", quote=FALSE, row.names = T, col.names = T)
}


CREBBP.dat[["expm.dat"]] = list()
version = "genes"
key = paste0(version)
print(paste0("working on: ",key))

expm = list()

feature.ids = c(CREBBP.dat$inc.genes,  CREBBP.dat$dec.genes)
expm$X = CREBBP.dat$X.genes[feature.ids,]
sample.ids = colnames(expm$X)
expm$W  = CREBBP.dat$W.hat[sample.ids, ]
expm$C1 = NULL

CREBBP.dat[["expm.dat"]][[key]] = expm
dumpCibersortx(key, expm, file.path(res.dir,"Cibersortx"))


saveRDS(CREBBP.dat, file.path(data.dir, "CREBBP.dat.rds"))
