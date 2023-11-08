library("data.table")
library("matrixStats")# colMeans
library("matrixcalc") #is.positive.definite
library("MIND")

source("analysis.utils.r")
source("simulate.expression.utils.r")

set.seed(2023)

# args = commandArgs(trailingOnly=TRUE)

# data.name   = as.character(args[1]) 
# source.col  = as.character(args[2]) 
# n           = as.numeric(args[3])
# sc.prior    = as.logical(args[4])
# prior.sample.ratio  = as.numeric(args[5])


data.name = "Lung" # "PBMC" or "Lung"
source.col = "decon.L1"
n = 100
sc.prior = F
prior.sample.ratio = 0 #from (0-1]

print(paste0("data.name: ", data.name))
print(paste0("source.col: ", source.col))
print(paste0("n: ", n))
print(paste0("sc.prior: ", sc.prior))
print(paste0("prior.sample.ratio: ", prior.sample.ratio))

#which copy
ts = 1:20
max_stds = 2
robust = T
qtl = 0.95

project.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/"
if(data.name == "PBMC"){
    # for prior 
    seu.file = "../Data/RNA/Simulation-PBMC/stephenson.clean.seu"
    sample.col = "sample_id"
    count.assay = "raw"
    
    # for data 
    if(source.col == "decon.L1"){
        data.dir = paste0("Data/RNA/Simulation-PBMC/sc-Stephenson_all.W_decon.L1_HEF.10k_k_5_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
    }else{
        data.dir = paste0("Data/RNA/Simulation-PBMC/sc-Stephenson_all.W_decon.L2_HEF.10k_k_7_m_600_n_", n, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
    }
}else{ 
    # for prior 
    seu.file = "../Data/RNA/Simulation-Lung/HLCA.clean.seu"
    sample.col = "sample"
    count.assay = "RNA"
    
    # for data 
    if(source.col == "decon.L1"){
        data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L1_HEF.10k_k_4_m_600_n_", n ,"_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
    }else{
        data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L2_HEF.10k_k_6_m_600_n_", n ,"_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
    }
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

cpm = function(x){
    x/sum(x+1)* 10**6
}

get_prior = function(sc, sample = NULL, cell_type = NULL, meta_sc = NULL, filter_pd = T) {
  
  if(is.null(meta_sc)) meta_sc = data.frame(sample = sample, cell_type = cell_type)
  meta_sc$sample = as.character(meta_sc$sample)
  sample = unique(meta_sc[, c('sample')])
  
  cell_type = sort(unique(meta_sc$cell_type))
  K = length(cell_type)
  
  cts = array(NA, dim = c(nrow(sc), length(sample), K))
  rownames(cts) = rownames(sc)
  colnames(cts) = sample
  dimnames(cts)[[3]] = cell_type
  for(j in sample) {
    for(k in dimnames(cts)[[3]]) {
      id = which(meta_sc$sample == j & meta_sc$cell_type == k)
      if(length(id) > 0) cts[,j,k] = rowMeans(sc[, id, drop = F])
    }
    id = which(colMeans(is.na(cts[,j,])) == 0)
    #cts[,j,id] = log2(cpm(cts[,j,id]) + 1) # make it log2 CPM + 1
    cts[,j,id] = cpm(cts[,j,id])  # make it log2 CPM + 1
  }
  
  cov = array(NA, dim = c(nrow(cts), K, K))
  rownames(cov) = rownames(sc)
  colnames(cov) = dimnames(cov)[[3]] = cell_type
  for(i in 1:nrow(cov)) {
    cov[i,,] = cov(cts[i,,], use = 'pairwise')
  }
  
  if(filter_pd) {
    gene_pd = apply(cov, 1, is.positive.definite)
    print(paste(sum(1-gene_pd), 'genes are filtered out because cell-type covariance matrix is not positive-definite (PD);'))
    print('The filtering can be disabled by setting filter_pd = FALSE. Note the prior cell-type covariance matrix for each gene is required to be PD.')
  } else gene_pd = 1:nrow(cov)
  cov <- cov[gene_pd,,]
  
  profile = matrix(NA, nrow(sc), K)
  rownames(profile) = rownames(sc)
  colnames(profile) = cell_type
  for(i in cell_type) {
      #profile[,i] = log2(cpm(rowMeans(sc[, meta_sc$cell_type == i])) + 1)
      profile[,i] = cpm(rowMeans(sc[, meta_sc$cell_type == i])) 
  }
  
  return(list(profile = profile[gene_pd,], covariance = cov)) # ctsExp = cts, 
}

if (sc.prior){
    seu.obj = readRDS(seu.file)
    if (prior.sample.ratio != 1){
        prior.sample.size = round(length(unique(seu.obj@meta.data[,sample.col])) * prior.sample.ratio)
        prior.sample.ids = sample (unique(seu.obj@meta.data[,sample.col]), prior.sample.size)
        seu.obj = seu.obj[,(seu.obj@meta.data[,sample.col] %in% prior.sample.ids)]
        print(paste0("total sample used for learning prior: ", length(unique(seu.obj@meta.data[,sample.col]))))
    }else{
        print("using all samples for learning prior")
    } 
   
    counts = seu.obj@assays[[count.assay]]@counts

    # format the gene names to be consistent with how we simulated X
    feature.ids = rownames(counts)
    feature.ids = vapply(1: length(feature.ids), function(j) gsub("\\-", ".", feature.ids[j]), character(1))
    feature.ids = vapply(1: length(feature.ids), function(j) gsub("\\/", ".", feature.ids[j]), character(1))
    rownames(counts) = feature.ids


    prior.list = list()
    for (t in ts){
        sim.data = sim.data.list[[t]]

        # prior 
        prior = get_prior(sc = as.matrix(counts[rownames(sim.data$X),]), 
                          sample = seu.obj@meta.data[,sample.col], 
                          cell_type = seu.obj@meta.data[[source.col]], 
                          meta_sc = NULL, filter_pd = F)

        #log 2 based CPM to 10^6, same code as in bMIND
        #adding a pseudo count 1 to the sample pseuodbulk. certain sample is completely 0
        #this can happen when you sample a W with certain celltype completely missing 
        #and a Z with only that celltype present 
        X.CPM = sim.data$X/repmat(t(as.matrix(colSums(sim.data$X)+1)), nrow(sim.data$X), 1) * 10**6
        W = sim.data$W
        K = ncol(W)

        # for those not pd, replace with rough prior from bulk, same code as in bMIND
        gene_pd = apply(prior$covariance, 1, is.positive.definite)
        print("number of pd covar genes: ")
        print(sum(gene_pd))
        for (j in which(!gene_pd)){
            prior$covariance[j,,] = diag(K) * var(X.CPM[j,]) / sum(colMeans(W)^2)
        }
        
        # variance is too close to 0 will also trigger some numerical instability. 
        # check it here and correct it 
        for (j in 1:dim(prior$covariance)[[1]]){
            if(any(abs(diag(prior$covariance[j,,])) <= 10**(-4))){
                prior$covariance[j,,] = diag(K) * var(X.CPM[j,]) / sum(colMeans(W)^2)
            }
        }
        prior.list[[t]] = prior
        
    }
    #free up memory 
    rm(seu.obj)
    rm(counts)
}

bMIND.mdl.list = list()

for (t in ts){
    
    sim.data = sim.data.list[[t]]

    # format name
    sample.ids = colnames(sim.data$X)
    source.ids = colnames(sim.data$W)
    new.sample.ids = paste0("sample", 1:ncol(sim.data$X))
    new.source.ids = vapply(1: length(source.ids), function(h) gsub(" ", "_", source.ids[h]), character(1))
                    
    colnames(sim.data$X) = new.sample.ids
    rownames(sim.data$W) = new.sample.ids
    colnames(sim.data$W) = new.source.ids

    if (sc.prior){
        # format name on prior as well. first order to simulation source order then replace with new source id
        prior = list(profile = prior.list[[t]]$profile[rownames(sim.data$X), source.ids],
                     covariance = prior.list[[t]]$covariance[rownames(sim.data$X), source.ids, source.ids])

        colnames(prior$profile) = new.source.ids
        dimnames(prior$covariance)[[2]] = new.source.ids
        dimnames(prior$covariance)[[3]] = new.source.ids
    }
                            
    # check for samples with pure 0 mixutre expression, remove that sample
    keep.samples = (colSums(sim.data$X) != 0)
    sim.data$X = sim.data$X[,keep.samples]
    sim.data$W = sim.data$W[keep.samples,]
    
    # run bMIND
    scale.constant = max(sim.data$X)
    if (sc.prior){
        b = bMIND(sim.data$X/scale.constant,
                  frac = sim.data$W,
                  profile = prior$profile[rownames(sim.data$X), colnames(sim.data$W)]/scale.constant,
                  covariance = prior$covariance[rownames(sim.data$X), 
                                                          colnames(sim.data$W), 
                                                          colnames(sim.data$W)]/scale.constant**2)
    }else{
        b = bMIND(sim.data$X/scale.constant, frac = sim.data$W)
    }
    # change back the source and sample name 
    dimnames(b$A)[[3]] = sample.ids[keep.samples]
    dimnames(b$A)[[2]] = source.ids
    colnames(b$mu)     = source.ids
    dimnames(b$Sigma_c)[[2]] = source.ids
    dimnames(b$Sigma_c)[[3]] = source.ids
    if (sc.prior){
        #change prior                    
        colnames(prior$profile) = source.ids
        dimnames(prior$covariance)[[2]] = source.ids
        dimnames(prior$covariance)[[3]] = source.ids
    }
    
    #format result
    bMIND.mdl = list()
   
    #params at original scale 
    bMIND.mdl$params.hat.orig = list(mus_hat = b$mu * scale.constant,
                                     sigmas_hat = b$Sigma_c * scale.constant**2)
    #tensor at original scale "logged"
    bMIND.mdl$Z.hat.orig = array(0, dim(sim.data$Z))
    dimnames(bMIND.mdl$Z.hat.orig) = dimnames(sim.data$Z)
       
    for(source.id in source.ids){
        #samples got removed due to pure 0 in X will automatically left as default as all 0 in Z
        bMIND.mdl$Z.hat.orig[source.id, dimnames(b$A)[[1]], dimnames(b$A)[[3]]] = b$A[, source.id,]*scale.constant
    }

    #hat evalulation scale 
    bMIND.mdl$params.hat.eval = list(mus = scale_feature_matrix(bMIND.mdl$params.hat.orig$mus_hat, 
                                                                1/sim.data$feature.scale.factor),
                                     sigmas = scale_feature_source_source(bMIND.mdl$params.hat.orig$sigmas_hat, 
                                                                          1/sim.data$feature.scale.factor))

    #recon evalulation scale 
    bMIND.mdl$Z.hat.eval = scale_source_feature_sample(bMIND.mdl$Z.hat.orig, 1/sim.data$feature.scale.factor)
    bMIND.mdl$params.recon.eval = calc_params_from_Z(bMIND.mdl$Z.hat.eval, max_sds = max_stds)
    
    #prior
    if (sc.prior){bMIND.mdl$prior = prior}
    bMIND.mdl.list[[t]] = bMIND.mdl

}
                            

str(bMIND.mdl.list)

for (t in ts){
    #hat
    bMIND.mdl.list[[t]]$moment.hat.corrs = get_moment_corrs(params.true = sim.data.list[[t]]$params.scale,
                                                            params.hat  = bMIND.mdl.list[[t]]$params.hat.eval, robust, qtl)

    #recon
    bMIND.mdl.list[[t]]$moment.recon.corrs = get_moment_corrs(params.true = sim.data.list[[t]]$params.scale,
                                                              params.hat  = bMIND.mdl.list[[t]]$params.recon.eval, robust, qtl)
} 

bMIND.mdl.list[[1]]$moment.hat.corrs

bMIND.mdl.list[[1]]$moment.recon.corrs

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
    
    
    bMIND.mdl.list[[t]]$Z.corrs = calc_Z_corrs(Z.true = Z.scale, 
                                               Z.hat = bMIND.mdl.list[[t]]$Z.hat.eval, 
                                               eval.feature.source = eval.feature.source,
                                               robust = robust, qtl = qtl)
    
}

end = Sys.time()
print(end - start)


colMedians(bMIND.mdl.list[[1]]$Z.corrs, na.rm = T)

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
                # everything is set to NA
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

#bMIND.mdl.list = readRDS(file.path(res.dir, paste0("bMIND.mdl.rough.prior.list.rds")))
for (t in ts){
    sim.data.list[[t]]$eval.feature.source = calc_variable_feature_source(sim.data.list[[t]]$Z.scale, 
                                                                          variable_thr = 0.1, max_sds = max_stds)
    bMIND.mdl.list[[t]]$joint.bulk.p = calc_joint_bulk_p_vals_array(X = sim.data.list[[t]]$X,
                                                                    W = sim.data.list[[t]]$W,
                                                                    Z = sim.data.list[[t]]$Z,
                                                                    
                                                                    Z.hat = none_neg_Z(bMIND.mdl.list[[t]]$Z.hat.orig), 
                                                                    eval.feature.source = sim.data.list[[t]]$eval.feature.source)
    
                                                                   #celltype by features
    bMIND.mdl.list[[t]]$joint.bulk.Z.hat.log10p.diff = t(-log10(bMIND.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) - t(-log10(bMIND.mdl.list[[t]]$joint.bulk.p[,,"bulk"]))
    bMIND.mdl.list[[t]]$joint.bulk.Z.hat.log10p      = t(-log10(bMIND.mdl.list[[t]]$joint.bulk.p[,,"Z.hat"])) 
}

bMIND.mdl.list[[t]]$joint.bulk.Z.hat.log10p.diff

colMedians(bMIND.mdl.list[[1]]$joint.bulk.Z.hat.log10p.diff, na.rm = T)

colMedians(bMIND.mdl.list[[2]]$joint.bulk.Z.hat.log10p.diff, na.rm = T)

saveRDS(bMIND.mdl.list,  file.path(res.dir, 
                                   paste0("bMIND.mdl.",
                                          if(sc.prior) paste0("sc.prior.", prior.sample.ratio) else "rough.prior",
                                          ".list.rds")))








