# This notebook takes
# n.copies: number of copies to generate 
# res.dir: path to the directory where the celltype level pseudobulk of Lung HLCA is
# all the parameters needed to simulate mixtures from the pseudobulk RDS file
# This notebooks calls the functions in simulate.expression.utils.r (checked in the 0.0.2 notebook)
# and save n.copies of simulated datasets (generated under different seed)
# in a subfolder of res.dir. the name of the subfolder contains all parameters used to generate the data
source("analysis.utils.r")
source("simulate.expression.utils.r")
set.seed(2023)

# partial.W = F # T or F
# source.col = "decon.L1" # "decon.L1" or "decon.L2"
# k = 4
# feature.set = "HEF.10k" # "HVF.10k" or "HEF.10k"

# partial.W = F # T or F
# source.col = "decon.L2" # "decon.L1" or "decon.L2"
# k = 6
# feature.set = "HEF.10k" # "HVF.10k" or "HEF.10k"
# n = 500

# partial.W = F # T or F
# source.col = "decon.L1" # "decon.L1" or "decon.L2"
# k = 4
# feature.set = "HEF.10k" # "HVF.10k" or "HEF.10k"
# n = 250

partial.W = F # T or F
source.col = "decon.L1" # "decon.L1" or "decon.L2"
k = 4
feature.set = "HEF.10k" # "HVF.10k" or "HEF.10k"
n = 100

n.copies = 20
res.dir = "../Data/RNA/Simulation-Lung/"

max_sds = 2
m = 600;

fit.dirichlet = F;
add.noise = T; variable_thr = 10^(-4); fill_thr = 10^(-4); expression.qtl = 0; # dont do additional filtering
enrich.low = F; entropy.thr.ratio = 0; enrich.ratio = 0;
scale.max_sds = Inf; scale.factor.thr = 10**(-4)

res.folder  = paste("sc-HLCA", if (partial.W) "partial.W" else "all.W", source.col, feature.set, 
                    "k",k,"m",m,"n",n,
                    "dirichlet",substr(as.character(fit.dirichlet), 1,1),
                    "noiseZ",substr(as.character(add.noise), 1,1), "varThr",variable_thr, "filThr",fill_thr, "expQtl",expression.qtl,
                    "enrich",substr(as.character(enrich.low), 1,1), "etpRat",entropy.thr.ratio, "enrichRat",enrich.ratio, 
                    "maxSds",max_sds, "scale.maxSds",scale.max_sds, "scale.factor.thr",scale.factor.thr, sep =  "_")

res.folder  = file.path(res.dir,res.folder)
if (!file.exists(res.folder)){dir.create(file.path(res.folder),recursive = T)}

res.folder

HLCA.pb = readRDS(file.path(res.dir, 
                            paste0("HLCA.pseudobulk.",
                                    source.col, ".",
                                    if (partial.W) "partial.W." else "all.W.",
                                    "rds")))
                                   
Z.src = HLCA.pb$Z
W.src = HLCA.pb$W

length(intersect(HLCA.pb$HEF, HLCA.pb$HVF))

rowMeans(HLCA.pb$params$mus[head(HLCA.pb$HVF),])

rowMeans(HLCA.pb$params$mus[head(HLCA.pb$HEF),])

W.src

colMeans(W.src > 0.05)

Z.src[1,,]

if (feature.set == "HVF.10k"){
    Z.src = Z.src[,HLCA.pb$HVF,]
}else{
    Z.src = Z.src[,HLCA.pb$HEF,]
}

head(dimnames(Z.src)[[2]])

sim.data.list = list()
for (t in 1:n.copies){
    sim.data = simulate_expression_mixture(Z.src = Z.src, W.src = W.src, 
                                           k = k, m = m, n = n,
                                           fit.dirichlet = fit.dirichlet, 
                                           add.noise = add.noise, variable_thr = variable_thr, fill_thr = fill_thr, expression.qtl = expression.qtl, 
                                           enrich.low = enrich.low, entropy.thr.ratio = entropy.thr.ratio, enrich.ratio = enrich.ratio, 
                                           max_sds = max_sds, scale.max_sds = scale.max_sds, scale.factor.thr = scale.factor.thr, seed.num = t)
    
    # save files for cibersortx to read
    fwrite(as.data.frame(sim.data$X),   
           file = file.path(res.folder,paste0("X.txt.", t)), 
           sep = "\t", quote=FALSE, row.names = T, col.names = T) 

    fwrite(as.data.frame(sim.data$W),   
           file = file.path(res.folder,paste0("W.txt.", t)), 
           sep = "\t", quote=FALSE, row.names = T, col.names = T) 

    # add to list
    sim.data.list[[t]] = sim.data
}
saveRDS(sim.data.list, file.path(res.folder, "sim.data.list.rds"))

sim.data$variable.feature.source


