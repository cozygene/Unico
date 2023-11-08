library("TCA")
source("Unico.r")
source("analysis.utils.r")

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
ts = 1:20

project.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/"
#relative path of the data dir to project.dir
if(source.col == "decon.L1"){
    data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L1_HEF.10k_k_4_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}else{
    data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L2_HEF.10k_k_6_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}

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

Unico.mdl.list = list()
for (t in ts){
    sim.data = sim.data.list[[t]]
    Unico.mdl = list()
    
    Unico.mdl$params.hat <- Unico(sim.data$X, sim.data$W, C1 = NULL, C2 = NULL, 
                                 parallel = TRUE, num_cores = NULL, 
                                 log_file = file.path(res.dir, paste0("Unico.t_",t,".log")))

    #capping at interal scale
    Unico.mdl$params.hat$sigmas_hat  = cap_values(Unico.mdl$params.hat$sigmas_hat,  max.val = 10 ** (4), min.val = 10**(-4))
    #hat internal scaled version
    Unico.mdl$params.hat$corrs     = calc_corrs_from_sigmas(Unico.mdl$params.hat$sigmas_hat)
    Unico.mdl$params.hat$entropies = calc_entropy(Unico.mdl$params.hat$corrs)

    #hat original input scale 
    Unico.mdl$params.hat.orig = list(mus_hat = scale_feature_matrix(Unico.mdl$params.hat$mus_hat, 
                                                               Unico.mdl$params.hat$scale.factor),
                                    sigmas_hat = scale_feature_source_source(Unico.mdl$params.hat$sigmas_hat, 
                                                                         Unico.mdl$params.hat$scale.factor),
                                    taus  = scale_feature_matrix(Unico.mdl$params.hat$taus_hat, 
                                                                 Unico.mdl$params.hat$scale.factor),
                                    corrs = Unico.mdl$params.hat$corrs,
                                    entropies = Unico.mdl$params.hat$entropies)

    #tensor at original scale
    Unico.mdl$Z.hat.orig = tensor(sim.data$X, W = sim.data$W, C1 = NULL, C2 = NULL, 
                                 Unico.mdl$params.hat, parallel = FALSE)


    #hat evalulation scale 
    Unico.mdl$params.hat.eval = list(mus = scale_feature_matrix(Unico.mdl$params.hat.orig$mus_hat, 
                                                               1/sim.data$feature.scale.factor),
                                    sigmas = scale_feature_source_source(Unico.mdl$params.hat.orig$sigmas_hat, 
                                                                         1/sim.data$feature.scale.factor),
                                    taus  = scale_feature_matrix(Unico.mdl$params.hat.orig$taus, 
                                                                 1/sim.data$feature.scale.factor),
                                    corrs = Unico.mdl$params.hat.orig$corrs,
                                    entropies = Unico.mdl$params.hat.orig$entropies)

    #recon evalulation scale 
    Unico.mdl$Z.hat.eval = scale_source_feature_sample(Unico.mdl$Z.hat.orig, 1/sim.data$feature.scale.factor)
    Unico.mdl$params.recon.eval = calc_params_from_Z(Unico.mdl$Z.hat.eval, max_sds = max_stds)
    
    Unico.mdl.list[[t]] = Unico.mdl
}
saveRDS(Unico.mdl.list,  file.path(res.dir, paste0("Unico.mdl.list.rds")))

str(Unico.mdl.list)

Unico.mdl.list[[1]]$params.hat.eval$sigmas[1,,]

Unico.mdl.list[[1]]$params.recon.eval$sigmas[1,,]

tca.mdl.list = list()

for (t in ts){
    sim.data = sim.data.list[[t]]
    tca.mdl = list()
    
    #directly operate in eval scale
    tca.mdl$params.hat.eval <- tca(sim.data$X.scale, sim.data$W, 
                                   constrain_mu = TRUE, log_fil=NULL)

    #capping at eval scale
    tca.mdl$params.hat.eval$sigmas_hat = cap_values(tca.mdl$params.hat.eval$sigmas_hat, max.val = 10**(4), min.val = 10**(-4))
    tca.mdl$Z.hat.eval = TCA::tensor(sim.data$X.scale, tca.mdl$params.hat.eval, log_fil=NULL)
    
    #format similar to Unico
    tca.mdl$params.hat.eval$mus    = tca.mdl$params.hat.eval$mus_hat
    tca.mdl$params.hat.eval$sigmas = matrix_to_diag_tensor(tca.mdl$params.hat.eval$sigmas_hat)
    tca.mdl$Z.hat.eval = list_2_array(tca.mdl$Z.hat.eval, colnames(sim.data$W))      
    
    #recon    
    tca.mdl$params.recon.eval = calc_params_from_Z(tca.mdl$Z.hat.eval, max_sds = max_stds)
    
    #tensor original scale (for orthognal signal check)
    tca.mdl$Z.hat.orig = scale_source_feature_sample(tca.mdl$Z.hat.eval, sim.data$feature.scale.factor)
                     
    tca.mdl.list[[t]] = tca.mdl
}

saveRDS(tca.mdl.list,  file.path(res.dir, "tca.mdl.list.rds"))

str(tca.mdl.list)

round(tca.mdl.list[[1]]$params.hat.eval$sigmas[3,,], 2)

round(tca.mdl.list[[1]]$params.recon.eval$sigmas[3,,], 2)

base.mdl.list = list()
for (t in ts){
    sim.data = sim.data.list[[t]]
    base.mdl = list()
    
    base.mdl$Z.hat.orig = copy(sim.data$Z.scale) 
    print("constructing simple Z by distribute X.scale by W")
    for (h in 1:dim(base.mdl$Z.hat.orig)[1]){
        base.mdl$Z.hat.orig[h,,] = sim.data$X * repmat(t(as.matrix(sim.data$W[,h])), dim(base.mdl$Z.hat)[2], 1)
    } 


    #recon evalulation scale 
    base.mdl$Z.hat.eval = scale_source_feature_sample(base.mdl$Z.hat.orig, 1/sim.data$feature.scale.factor)
    base.mdl$params.recon.eval = calc_params_from_Z(base.mdl$Z.hat.eval, max_sds = max_stds)

    base.mdl.list[[t]] = base.mdl
}
saveRDS(base.mdl.list,  file.path(res.dir, paste0("base.mdl.list.rds")))

str(base.mdl.list)


