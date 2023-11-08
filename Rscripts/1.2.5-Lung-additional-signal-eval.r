library(matrixStats)
source("analysis.utils.r")
source("simulate.expression.utils.r")

# source.col = "decon.L1" 
# n= 500

# source.col = "decon.L2" 
# n= 500

# source.col = "decon.L1" 
# n= 250

source.col = "decon.L1" 
n= 100

max_stds = 2
project.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/"

#relative path of the data dir to project.dir
if(source.col == "decon.L1"){
    data.dir   = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L1_HEF.10k_k_4_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds,"_scale.maxSds_Inf_scale.factor.thr_1e-04/")
    res.dir    = paste0("Result/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L1_HEF.10k_k_4_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds,"_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}else{
    data.dir   = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L2_HEF.10k_k_6_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds,"_scale.maxSds_Inf_scale.factor.thr_1e-04/")
    res.dir    = paste0("Result/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L2_HEF.10k_k_6_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds,"_scale.maxSds_Inf_scale.factor.thr_1e-04/")
}

ts = 1:20

data.dir   = file.path(project.dir, data.dir)
res.dir   = file.path(project.dir, res.dir)

if (!file.exists(res.dir)){print("no result in the result directory")}

sim.data.list = readRDS(file.path(data.dir, "sim.data.list.rds"))

# load TCA, baseline
tca.mdl.list  = readRDS(file.path(res.dir, paste0("tca.mdl.list.rds")))
base.mdl.list = readRDS(file.path(res.dir, paste0("base.mdl.list.rds")))
cibersortx.mdl.list = readRDS(file.path(res.dir, paste0("cibersortx.mdl.list.rds")))
Unico.mdl.list = readRDS(file.path(res.dir, paste0("Unico.mdl.list.rds")))



#X: matrix of features by sample, bulk expression or methylation 
#W: matric of samples by sources, proportions
#Z: array of sources by features by samples, ground truth Z
#Z.hat: array of sources by features by samples, estiamted Z
#eval.feature.source: matrix of logical features by sources, 
#if True, evaluate this feature-source, if False, skip it. all components will have pval1
#num.stds: numeric, number of standard deviations away to be considered as outliers and removed from regression 
#constant_thr: numeric, miminal amount of variation in Z.hat to be considered worth fitting the regression 
#if fall below this number, Z.hat automoatically get pval 1. regression only fit on the rest of the components

#per feature-source fit the following regression
#log(1 + Z[feature-source]) ~ 1 + W[drop one celltype] + log(1 + bulk[feature]) + log(1 + Z.hat[feature-source])
#add the - (min Z.hat) to all Z.hat if there are any negative values in that feature,source
#remove 2sd outliers based on bulk

#return an array of pvalues that is source.ids by features by compents participates in the regression 

calc_joint_bulk_p_vals_array = function(X, W, Z, Z.hat, eval.feature.source, num.stds = 2, constant_thr = 10**(-4)){
    # if there are negative value, shift entire distribution 
    Z.hat = none_neg_Z(Z.hat)

    m = nrow(X)
    k = ncol(W)
    feature.ids = rownames(X)
    source.ids  = colnames(W)
    
    # init to all 1 pval
    p.vals.array = array(1, c(k, m, (1 + 1 + 1 + (k-1)))) # intercept, bulk, Z[,jl], W related,
    dimnames(p.vals.array)[[1]] = source.ids  
    dimnames(p.vals.array)[[2]] = feature.ids
    dimnames(p.vals.array)[[3]] = c("(Intercept)", "bulk", "Z.hat", colnames(W)[-1])

    for (source.id in source.ids){
        for (feature.id in feature.ids){
            if (!eval.feature.source[feature.id, source.id]){
                p.vals.array[source.id, feature.id, ] = NA
                next # skip this feature source due to no meaningful variable in the ground truth Z
                # everything is set to 1
            }else{
                df = data.frame(Z     = log(1 + Z[source.id,feature.id,]),     
                                bulk  = log(1 + X[feature.id, ]), 
                                Z.hat = log(1 + Z.hat[source.id,feature.id,]))
                df = cbind(df, data.frame(W[,-1])) # drop one celltype
                
                #filter out 2SD
                x =  df$bulk
                mask = abs(x - mean(x)) <= num.stds * sd(x)
                df = df[mask, ]

                if(sd(df$Z.hat) <= constant_thr){
                    #if no variance detected in celltype h, gene j, 
                    df = df[, colnames(df)!= "Z.hat"]
                    fit = lm(formula = Z ~ 1 + ., data = df) 
                    coef = summary(fit)$coefficients
                    #prevent lm replacing the special characters in the colnames so hard code this
                    p.vals.array[source.id, feature.id, c("(Intercept)", "bulk", colnames(W[,-1]))] = t(coef[, "Pr(>|t|)", drop = F])
                }else{
                    fit = lm(formula = Z ~ 1 + ., data = df) 
                    coef = summary(fit)$coefficients
                    p.vals.array[source.id, feature.id, ] = t(coef[, "Pr(>|t|)", drop = F])
                }
             }
        }#end of feature
    }#end of source
    return(p.vals.array)
}

for (t in ts){
    sim.data.list[[t]]$eval.feature.source = calc_variable_feature_source(sim.data.list[[t]]$Z.scale, 
                                                                          variable_thr = 0.1, max_sds = max_stds)

    base.mdl.list[[t]]$joint.bulk.p = calc_joint_bulk_p_vals_array(X = sim.data.list[[t]]$X,
                                                                   W = sim.data.list[[t]]$W,
                                                                   Z = sim.data.list[[t]]$Z,
                                                                   Z.hat = base.mdl.list[[t]]$Z.hat.orig,
                                                                   eval.feature.source = sim.data.list[[t]]$eval.feature.source)

    
    cibersortx.mdl.list[[t]]$joint.bulk.p = calc_joint_bulk_p_vals_array(X = sim.data.list[[t]]$X,
                                                                         W = sim.data.list[[t]]$W,
                                                                         Z = sim.data.list[[t]]$Z,
                                                                         Z.hat = cibersortx.mdl.list[[t]]$Z.hat.orig,
                                                                         eval.feature.source = sim.data.list[[t]]$eval.feature.source)

    tca.mdl.list[[t]]$joint.bulk.p = calc_joint_bulk_p_vals_array(X = sim.data.list[[t]]$X,
                                                                  W = sim.data.list[[t]]$W,
                                                                  Z = sim.data.list[[t]]$Z,
                                                                  Z.hat = tca.mdl.list[[t]]$Z.hat.orig,
                                                                  eval.feature.source = sim.data.list[[t]]$eval.feature.source)
    
    Unico.mdl.list[[t]]$joint.bulk.p = calc_joint_bulk_p_vals_array(X = sim.data.list[[t]]$X,
                                                                   W = sim.data.list[[t]]$W,
                                                                   Z = sim.data.list[[t]]$Z,
                                                                   Z.hat = Unico.mdl.list[[t]]$Z.hat.orig,
                                                                   eval.feature.source = sim.data.list[[t]]$eval.feature.source)
     
    
    
    #celltype by features
    base.mdl.list[[t]]$joint.bulk.Z.hat.log10p.diff = t(-log10(base.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) - t(-log10(base.mdl.list[[t]]$joint.bulk.p[,,"bulk"]))
    base.mdl.list[[t]]$joint.bulk.Z.hat.log10p      = t(-log10(base.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) 
    
    cibersortx.mdl.list[[t]]$joint.bulk.Z.hat.log10p.diff = t(-log10(cibersortx.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) - t(-log10(cibersortx.mdl.list[[t]]$joint.bulk.p[,,"bulk"]))
    cibersortx.mdl.list[[t]]$joint.bulk.Z.hat.log10p      = t(-log10(cibersortx.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) 
   
    tca.mdl.list[[t]]$joint.bulk.Z.hat.log10p.diff = t(-log10(tca.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) - t(-log10(tca.mdl.list[[t]]$joint.bulk.p[,,"bulk"]))
    tca.mdl.list[[t]]$joint.bulk.Z.hat.log10p      = t(-log10(tca.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"]))
    
    Unico.mdl.list[[t]]$joint.bulk.Z.hat.log10p.diff = t(-log10(Unico.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) - t(-log10(Unico.mdl.list[[t]]$joint.bulk.p[,,"bulk"])) 
    Unico.mdl.list[[t]]$joint.bulk.Z.hat.log10p      = t(-log10(Unico.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"]))

}

# just on the first round
colMedians(base.mdl.list[[1]]$joint.bulk.Z.hat.log10p.diff, na.rm = T)

colMedians(cibersortx.mdl.list[[1]]$joint.bulk.Z.hat.log10p.diff, na.rm = T)

colMedians(tca.mdl.list[[1]]$joint.bulk.Z.hat.log10p.diff, na.rm = T)

colMedians(Unico.mdl.list[[1]]$joint.bulk.Z.hat.log10p.diff, na.rm = T)

# #Time permit
# options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)

# plot(-log10(Unico.mdl.list[[1]]$joint.bulk.p[1,,"bulk"]), 
#      -log10(Unico.mdl.list[[1]]$joint.bulk.p[1,,"Z.hat"]), 
#      xlim=range(0,40), ylim=range(0,40), col = alpha("red", 0.1))

# options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)
# plot(-log10(tca.mdl.list[[1]]$joint.bulk.p[1,,"bulk"]), 
#      -log10(tca.mdl.list[[1]]$joint.bulk.p[1,,"Z.hat"]),
#      xlim=range(0,40), ylim=range(0,40), col = alpha("red", 0.1))

# options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)

# plot(-log10(tca.mdl.list[[1]]$joint.bulk.p[1,,"Z.hat"]), 
#      -log10(Unico.mdl.list[[1]]$joint.bulk.p[1,,"Z.hat"]), 
#      xlim=range(0,40), ylim=range(0,40), col = alpha("red", 0.1))

saveRDS(base.mdl.list,       file.path(res.dir, paste0("base.mdl.list.rds"))) 
saveRDS(cibersortx.mdl.list, file.path(res.dir, paste0("cibersortx.mdl.list.rds")))
saveRDS(tca.mdl.list,        file.path(res.dir, paste0("tca.mdl.list.rds")))
saveRDS(Unico.mdl.list,       file.path(res.dir, paste0("Unico.mdl.list.rds")))




