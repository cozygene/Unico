# simulation 
library(ggplot2)
library(ggpubr)
library(pracma)
source("Unico.r")
source("analysis.utils.r")
set.seed(2023)

#countour with know distribution 
get_contour_df = function(mu, sigma){
    data.grid <- expand.grid(s.1 = seq(mu[1] - 4*sqrt(sigma[1,1]), mu[1] + 4*sqrt(sigma[1,1]), length.out=200), 
                             s.2 = seq(mu[2] - 4*sqrt(sigma[2,2]), mu[2] + 4*sqrt(sigma[2,2]), length.out=200))

    q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu, sigma = sigma))
    return(q.samp)
}

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 150)
mu <- c(.5, -.5)
sigma <- matrix(c(1,   -0.5,
                  -0.5,   1), nrow=2)

ggplot(get_contour_df(mu, sigma), aes(x=s.1, y=s.2, z=prob)) + 
       geom_contour()



#real Z
PBMC.pb = readRDS("../Data/RNA/Simulation-PBMC/pbmc.pseudobulk.decon.L1.all.W.rds")

PBMC.entropy = as.matrix(PBMC.pb$entropy[PBMC.pb$HEF,])
PBMC.low.gs  = rownames(PBMC.entropy[PBMC.entropy < quantile(PBMC.entropy, 0.25),,drop = F])
PBMC.high.gs = rownames(PBMC.entropy[PBMC.entropy > quantile(PBMC.entropy, 0.5),,drop = F])

# three relatively distinct population
source.ids = c("CD4", "B", "Mono")

# low entropy 
tmp = PBMC.pb$params$mus[PBMC.low.gs, source.ids]
tmp.low = rownames(tmp[(rowSums(tmp> 150) == 3) & (rowSums(tmp< 300) == 3), ])
tmp.low = tmp.low[order(PBMC.pb$entropy[tmp.low,])]
tmp[tmp.low,]

PBMC.pb$params$mus["CHCHD2",]

PBMC.pb$params$corrs["CHCHD2",source.ids,source.ids]

# high entropy express in all three cell types
tmp = PBMC.pb$params$mus[PBMC.high.gs, source.ids]
tmp.high = rownames(tmp[(rowSums(tmp> 100) == 3) & (rowSums(tmp< 300) == 3), ])
tmp.high = tmp.high[order(PBMC.pb$entropy[tmp.high,])]
tmp[tmp.high,]

PBMC.pb$params$corrs["SLC2A3",source.ids,source.ids]

PBMC.pb$params$mus["SLC2A3",]

# first low, second high 
PBMC.pb$entropy[c("CHCHD2", "SLC2A3"), ]

Z = PBMC.pb$Z[source.ids, c("CHCHD2", "SLC2A3", "NME2", "NDUFB11", "TRIR", "PRR13", "COX6C", "CDC42", "CDC42", "ATP5F1D"),]

source.ids  = dimnames(Z)[[1]]
feature.ids = dimnames(Z)[[2]]
sample.ids  = dimnames(Z)[[3]]

k = length(source.ids)
m = length(feature.ids)
n = length(sample.ids)

# n = 1000
# k = 3
# m = 4

# source.ids  = paste0("soure.", 1:3) 
# feature.ids = paste0("feature.",1:4)
# sample.ids  = paste0("sample.", 1:n)


# Z = array(0, c(k,m,n))
# dimnames(Z)[[1]] = source.ids
# dimnames(Z)[[2]] = feature.ids
# dimnames(Z)[[3]] = sample.ids
# mu <- c(X = 10, Y = 10, Z = 10)
# R1 <- matrix(c(1,              sqrt(2) * 0.75, sqrt(3) * 0.95,
#                sqrt(2) * 0.75, 2,              sqrt(6) * 0.85, 
#                sqrt(3) * 0.95, sqrt(6) * 0.85, 3 ), 
#              nrow = 3, ncol = 3, byrow = TRUE)


# R2 <- matrix(c(1,              sqrt(2) * 0.25, sqrt(3) * 0.05,
#                sqrt(2) * 0.25, 2,              sqrt(6) * 0.15, 
#                sqrt(3) * 0.05, sqrt(6) * 0.15, 3 ), 
#              nrow = 3, ncol = 3, byrow = TRUE)
# Z[,1,] = t(MASS::mvrnorm(n, mu = mu, Sigma = R1))
# Z[,2,] = t(MASS::mvrnorm(n, mu = mu, Sigma = R2))

# R1
# R2

W = matrix(1/k, n, k)
rownames(W) = sample.ids
colnames(W) = source.ids


X = matrix(0, m, n)
rownames(X) = feature.ids
colnames(X) = sample.ids   
for (h in 1:k){
    X = X + Z[h,,] * repmat(t(W[,h]), m ,1)
}

head(W)

X 

# fit a fake model to get the structure
# original parameter learning is garbage, but has the right structure,
# under this simulation there is no variation in W, so it is impossible to learn



Unico.mdl = list()
Unico.mdl$params.hat <- Unico(X, W, C1 = NULL, C2 = NULL, 
                              parallel = TRUE, num_cores = NULL, log_file = NULL)


params.true = calc_params_from_Z(Z, max_sds = Inf)
str(params.true)

params.true$mus

params.true$sigmas[1,,]

params.true$sigmas[2,,]

# plug in true parameters
Unico.mdl$params.hat$mus_hat    = scale_feature_matrix(params.true$mus, 1/Unico.mdl$params.hat$scale.factor)
Unico.mdl$params.hat$sigmas_hat = scale_feature_source_source(params.true$sigmas, 1/Unico.mdl$params.hat$scale.factor)

Unico.mdl$params.hat$mus_hat

Unico.mdl$params.hat$sigmas_hat[1,,]

Unico.mdl$params.hat$sigmas_hat[2,,]

C1 = Unico.mdl$params.hat$C1
C2 = Unico.mdl$params.hat$C2

betas_hat  = Unico.mdl$params.hat$betas_hat
gammas_hat = Unico.mdl$params.hat$gammas_hat
taus_hat   = Unico.mdl$params.hat$taus_hat

Unico.mdl$Z.hat <- tensor(X, W = W, C1 = NULL, C2 = NULL, 
                         Unico.mdl$params.hat, parallel = FALSE)

genes = list("Z" = Z, "X" = X, "W" = W, "C1"= C1, "C2" = C2, 
             "betas_hat" = betas_hat, "gammas_hat" = gammas_hat, "taus_hat" = taus_hat,
             "Unico.mdl" = Unico.mdl)
for (j in 1:length(feature.ids)){
    feature.id = feature.ids[j]
    Z_ept = t(repmat(as.matrix(params.true$mus[j,]), 1, n)) + C1  %*% t(matrix(gammas_hat[j,], nrow = dim(W)[2], byrow = TRUE))

    covar_z_x = W %*% params.true$sigmas[j,,]

    WW_tensor = get_WW_tensor(W)
    var_X_j = matrix(vapply(1:n, function(i) sum(WW_tensor[i,,] * params.true$sigmas[j,,]) + taus_hat[j], numeric(1)))

    X_j_ept = get_X_j_ept (params.true$mus[j,], W, C1, gammas_hat, C2, betas_hat)
    X_diff = X[j,] - X_j_ept                    

    # external implementation to get conditional mean: sample by source 
    mu_cond = Z_ept + (covar_z_x /repmat(var_X_j, 1, k)) * repmat(X_diff, 1, k)

    # conditional covaraince: sample by source by source
    sigmas_cond = array(0, c(n,k,k))
    dimnames(sigmas_cond)[[1]] = sample.ids      
    dimnames(sigmas_cond)[[2]] = source.ids      
    dimnames(sigmas_cond)[[3]] = source.ids            
    for (i in 1:n){
        simgas_j = params.true$sigmas[j,,]
        W_i = W[i,]
        sigmas_cond[i,,] = simgas_j - t(t(W_i) %*% simgas_j) %*% ((var_X_j[i,])**(-1)) %*% (t(W_i) %*% simgas_j)
    } 
                            
    # conditional correlation: sample by source by source
    corr_cond = calc_corrs_from_sigmas(sigmas_cond)
    genes[[feature.id]] = list(Z_ept = Z_ept, covar_z_x = covar_z_x, var_X_j = var_X_j, 
                               mu_cond = mu_cond, sigmas_cond = sigmas_cond, corr_cond = corr_cond)
}

str(genes[[feature.ids[1]]])

j = 1
for(h in 1:k){
    print(plot(Z[h,j,], Unico.mdl$Z.hat[h,j,]) + 
          title (paste0("source ",h, " rob cor : " ,round(safe_cor(Unico.mdl$Z.hat[h,j,], Z[h,j,], robust = T, qtl = 0.95), 3))))
}

j = 2
for(h in 1:k){
    print(plot(Z[h,j,], Unico.mdl$Z.hat[h,j,]) + 
          title (paste0("source ",h, " rob cor : " ,round(safe_cor(Unico.mdl$Z.hat[h,j,], Z[h,j,], robust = T, qtl = 0.95), 3))))
}

Unico.mdl$Z.hat[,feature.ids[1],1]

genes[[feature.ids[1]]]$mu_cond[1,]

# on the conditional corr Z|X

# check on the low entropy gene first sample
genes[[feature.ids[1]]]$corr_cond[1,,]

# check on the high entropy gene
genes[[feature.ids[2]]]$corr_cond[1,,]



visual.source.ids = c("CD4", "Mono")

# low entropy
# check on first sample but all samples share the same  sigmas_cond due to the same W
genes[[feature.ids[1]]]$corr_cond[1, visual.source.ids,visual.source.ids]

# high entropy
genes[[feature.ids[2]]]$corr_cond[1, visual.source.ids,visual.source.ids]

options(repr.plot.width = 6, repr.plot.height =6, repr.plot.res = 100)
feature.id = feature.ids[1]

for (visual.sample.id in 1:10){
    mu    = genes[[feature.id]]$mu_cond[visual.sample.id, visual.source.ids]
    sigma = matrix(genes[[feature.id]]$sigmas_cond[visual.sample.id, visual.source.ids,visual.source.ids], nrow=2)
    plot.df = get_contour_df(mu, sigma)
    p= ggplot(plot.df, aes(x=s.1, y=s.2, z=prob)) +  
       geom_contour() + 
       geom_point(aes(x=mu[1], y=mu[2]), colour="blue") + 
       geom_point(aes(x=Z[,feature.id,visual.sample.id][visual.source.ids[1]], 
                      y=Z[,feature.id,visual.sample.id][visual.source.ids[2]]), colour="red")
    print(p)
}

options(repr.plot.width = 6, repr.plot.height =6, repr.plot.res = 100)
feature.id = feature.ids[2]

for (visual.sample.id  in 1:10){
    mu    = genes[[feature.id]]$mu_cond[visual.sample.id, visual.source.ids]
    sigma = matrix(genes[[feature.id]]$sigmas_cond[visual.sample.id, visual.source.ids,visual.source.ids], nrow=2)
    plot.df = get_contour_df(mu, sigma)
    p= ggplot(plot.df, aes(x=s.1, y=s.2, z=prob)) +  
       geom_contour() + 
       geom_point(aes(x=mu[1], y=mu[2]), colour="blue") + 
       geom_point(aes(x=Z[,feature.id,visual.sample.id][visual.source.ids[1]], 
                      y=Z[,feature.id,visual.sample.id][visual.source.ids[2]]), colour="red")
    print(p)
}



if (!file.exists("../Figure/Illustrative/")){dir.create(file.path("../Figure/Illustrative/"))}

saveRDS(genes, "../Figure/Illustrative/toy.genes.rds")




