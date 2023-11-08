set.seed(2023)

permute.size = 5
res.dir    = "../Result/Runtime/mdl/"

time.df = c()
for (t in 1:permute.size){
    bMIND.mdl = readRDS(file.path(res.dir, paste0("bMIND.mcmck_6.m_1000.n_500.t_", t, ".rds")))
    time.df = rbind(time.df, c(as.numeric(bMIND.mdl$mcmc.assoc.time), sum(1/bMIND.mdl$mcmc.deconv$pval)))
}

time.df = as.data.frame(time.df)
rownames(time.df) = paste0("iteration.", 1:permute.size)
colnames(time.df) = c("time", "nsamples")

time.df$per.sample = time.df$time/ time.df$nsamples

time.df

bon.cor.p = 0.05/(150000 * 6)
sig.nsamples = round(1/bon.cor.p)

min(time.df$per.sample) * sig.nsamples/3600

max(time.df$per.sample) * sig.nsamples/3600

median(time.df$per.sample) * sig.nsamples/3600
