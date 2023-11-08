library("data.table")
library("matrixStats")# colMeans

source("analysis.utils.r")
source("simulate.expression.utils.r")

set.seed(2023)

#source.col = "decon.L1"
#n = 500

#source.col = "decon.L2"
#n = 500

# source.col = "decon.L1"
# n = 250

source.col = "decon.L1"
n = 100

max_stds = 2
project.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/"
#relative path of the data dir to project.dir
if(source.col == "decon.L1"){
    data.dir = paste0("Data/RNA/Simulation-PBMC/sc-Stephenson_all.W_decon.L1_HEF.10k_k_5_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}else{
    data.dir = paste0("Data/RNA/Simulation-PBMC/sc-Stephenson_all.W_decon.L2_HEF.10k_k_7_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}

ts = 1:20
robust = T
qtl = 0.95

data.dir = file.path(project.dir, data.dir)
data.name = strsplit(data.dir, "/")[[1]]
res.dir = file.path(project.dir, 
                    paste(c("Result", data.name[(length(data.name)-2) : length(data.name)]), collapse = "/"))


if (!file.exists(res.dir)){dir.create(file.path(res.dir),recursive = T)}
print(project.dir)
print(data.dir)
print(res.dir)

sim.data.list = readRDS(file.path(data.dir, "sim.data.list.rds"))

str(sim.data.list)

cibersortx.mdl.list = list()
for (t in ts){
    sim.data = sim.data.list[[t]]
    cibersortx.mdl = list()
    
    print("read in cibersortx's estimated Z")
    cibersortx.mdl$Z.hat.orig = copy(sim.data.list[[t]]$Z) 
    cibersortx.mdl$Z.hat.orig[,,] = 0
    
    for (source.id in dimnames(cibersortx.mdl$Z.hat.orig)[[1]]){
        Z.hat.file = file.path(res.dir, paste0("CIBERSORTxHiRes_NA_", gsub("\\ ", "", source.id) , "_Window",round(4 * ncol(sim.data$W)), ".txt.", t))
        cibersortx.mdl$Z.hat.orig[source.id,,] = as.matrix(data.frame(fread(Z.hat.file), row.names=1))
    }

    #recon evalulation scale 
    cibersortx.mdl$Z.hat.eval = scale_source_feature_sample(cibersortx.mdl$Z.hat.orig, 1/sim.data$feature.scale.factor)
    cibersortx.mdl$params.recon.eval = calc_params_from_Z(cibersortx.mdl$Z.hat.eval, max_sds = max_stds)

    cibersortx.mdl.list[[t]] = cibersortx.mdl
}

saveRDS(cibersortx.mdl.list,  file.path(res.dir, paste0("cibersortx.mdl.list.rds")))

for (t in ts){
    #recon
    cibersortx.mdl.list[[t]]$moment.recon.corrs = get_moment_corrs(params.true = sim.data.list[[t]]$ params.scale,
                                                                   params.hat = cibersortx.mdl.list[[t]]$params.recon.eval, robust, qtl)
}

cibersortx.mdl.list[[1]]$moment.recon.corrs

start = Sys.time()

# first adding eval.features.source, 
# we dont want to compare correlation against gene-celltype that is not really active and we injected noise
#for (t in ts){
for (t in ts){
    Z.scale = sim.data.list[[t]]$Z.scale
    eval.feature.source = calc_variable_feature_source(Z = Z.scale, variable_thr = 0.1, max_sds = max_stds)
    sim.data.list[[t]]$eval.feature.source = eval.feature.source
    
    #all the gene-celltype that has inject noise are not entering the evaluation
    #first explicitly set this
    #edge case: the original gene is not expressed in any celltypes (after outliers removed), so all celltypes are noise
    #after scaling. the injeted noise passed the eval.feature.source. 
    sim.data.list[[t]]$eval.feature.source[!sim.data.list[[t]]$variable.feature.source] = FALSE
    assert(!any(sim.data.list[[t]]$eval.feature.source[!sim.data.list[[t]]$variable.feature.source]))
    
    
    cibersortx.mdl.list[[t]]$Z.corrs = calc_Z_corrs(Z.true = Z.scale, 
                                                    Z.hat = cibersortx.mdl.list[[t]]$Z.hat.eval, 
                                                    eval.feature.source = eval.feature.source,
                                                    robust = robust, qtl = qtl)
    
}

end = Sys.time()
print(end - start)


str(cibersortx.mdl.list[[1]])

cibersortx.mdl.list[[1]]$ Z.corrs

saveRDS(cibersortx.mdl.list, file.path(res.dir, paste0("cibersortx.mdl.list.rds")))


