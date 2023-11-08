library(matrixStats)# colMeans

library("ggplot2")
library("ggpubr")

source("analysis.utils.r")
source("simulate.expression.utils.r")

set.seed(2023)

# source.col = "decon.L1"
# n = 500

# source.col = "decon.L2"
# n = 500

# source.col = "decon.L1"
# n = 250

source.col = "decon.L1"
n = 100

max_stds = 2
#which copy
ts = 1:20

project.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/"

#relative path of the data dir to project.dir
if(source.col == "decon.L1"){
    data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L1_HEF.10k_k_4_m_600_n_", n ,"_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}else{
    data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L2_HEF.10k_k_6_m_600_n_", n ,"_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}


robust = T
qtl = 0.95

data.dir = file.path(project.dir, data.dir)
data.name = strsplit(data.dir, "/")[[1]]
res.dir = file.path(project.dir, 
                    paste(c("Result", data.name[(length(data.name)-2) : length(data.name)]), collapse = "/"))

if (!file.exists(res.dir)){print("no result in the result directory")}
print(data.dir)
print(res.dir)

sim.data.list = readRDS(file.path(data.dir, "sim.data.list.rds"))

# load TCA, baseline
tca.mdl.list  = readRDS(file.path(res.dir, paste0("tca.mdl.list.rds")))
base.mdl.list = readRDS(file.path(res.dir, paste0("base.mdl.list.rds")))
Unico.mdl.list = readRDS(file.path(res.dir, paste0("Unico.mdl.list.rds")))

str(Unico.mdl.list)

for (t in ts){
    #hat
    tca.mdl.list[[t]]$moment.hat.corrs    = get_moment_corrs(params.true = sim.data.list[[t]]$ params.scale,
                                                             params.hat = tca.mdl.list[[t]]$params.hat.eval, robust,qtl)
    Unico.mdl.list[[t]]$moment.hat.corrs   = get_moment_corrs(params.true = sim.data.list[[t]]$ params.scale,
                                                             params.hat = Unico.mdl.list[[t]]$params.hat.eval, robust,qtl)
    
    #recon
    base.mdl.list[[t]]$moment.recon.corrs = get_moment_corrs(params.true = sim.data.list[[t]]$ params.scale,
                                                             params.hat = base.mdl.list[[t]]$params.recon.eval, robust,qtl)
    tca.mdl.list[[t]]$moment.recon.corrs  = get_moment_corrs(params.true = sim.data.list[[t]]$ params.scale,
                                                             params.hat = tca.mdl.list[[t]]$params.recon.eval, robust,qtl)
    Unico.mdl.list[[t]]$moment.recon.corrs = get_moment_corrs(params.true = sim.data.list[[t]]$ params.scale,
                                                             params.hat = Unico.mdl.list[[t]]$params.recon.eval, robust,qtl)
}

tca.mdl.list[[1]]$moment.hat.corrs

Unico.mdl.list[[1]]$moment.hat.corrs

base.mdl.list[[1]]$moment.recon.corrs

start = Sys.time()

# first adding eval.features.source, 
# we dont want to compare correlation against gene-celltype that is not really active and we injected noise
#for (t in ts){
for (t in ts){
    Z.scale = sim.data.list[[t]]$Z.scale
    eval.feature.source = calc_variable_feature_source(Z = Z.scale, variable_thr = 0.1, max_sds = max_stds)
    sim.data.list[[t]]$eval.feature.source = eval.feature.source
    
    #all the gene-celltype that has inject noise are not entering the evaluation
    assert(!any(sim.data.list[[t]]$eval.feature.source[!sim.data.list[[t]]$variable.feature.source]))
    
    base.mdl.list[[t]]$Z.corrs = calc_Z_corrs(Z.true = Z.scale, 
                                              Z.hat = base.mdl.list[[t]]$Z.hat.eval, 
                                              eval.feature.source = eval.feature.source,
                                              robust = robust, qtl = qtl)
    
    tca.mdl.list[[t]]$Z.corrs  = calc_Z_corrs(Z.true = Z.scale,
                                              Z.hat = tca.mdl.list[[t]]$Z.hat.eval, 
                                              eval.feature.source = eval.feature.source,
                                              robust = robust, qtl = qtl)
    
    Unico.mdl.list[[t]]$Z.corrs = calc_Z_corrs(Z.true = Z.scale,
                                              Z.hat = Unico.mdl.list[[t]]$Z.hat.eval, 
                                              eval.feature.source = eval.feature.source,
                                              robust = robust, qtl = qtl)
}

end = Sys.time()
print(end - start)


str(Unico.mdl.list)

Unico.mdl.list[[1]]$ Z.corrs

saveRDS(base.mdl.list,       file.path(res.dir, paste0("base.mdl.list.rds"))) 
saveRDS(tca.mdl.list,        file.path(res.dir, paste0("tca.mdl.list.rds")))
saveRDS(Unico.mdl.list,       file.path(res.dir, paste0("Unico.mdl.list.rds")))

##placeholder
#cibersortx.mdl.list = copy(base.mdl.list)
#saveRDS(cibersortx.mdl.list, file.path(res.dir, paste0("cibersortx.mdl.list.rds")))


