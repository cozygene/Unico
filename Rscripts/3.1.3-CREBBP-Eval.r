# This notebook reads in all the fitted model and tensor from different methods
# calculates the DE signal in log transformed estiamted celltype profile on those known down and up regulated genes
# calculates the "effect size" when we log transform the estiamted celltype profile and fit it against proportion and mutation status 
library("testit")
library("matrixStats")

source("analysis.utils.r")
set.seed(2023)

data.dir = file.path("../Data/RNA/CREBBP/")
res.dir  = file.path("../Result/RNA/CREBBP/")
CREBBP.dat = readRDS(file.path(data.dir, "CREBBP.dat.rds"))

inc.genes  = CREBBP.dat$inc.genes
dec.genes  = CREBBP.dat$dec.genes
samples.mt = CREBBP.dat$samples.mt
samples.wt = CREBBP.dat$samples.wt

print(paste0("length inc.genes: ", length(inc.genes)))
print(paste0("length dec.genes: ", length(dec.genes)))
print(paste0("length samples.mt: ", length(samples.mt)))
print(paste0("length samples.wt: ", length(samples.wt)))

key     = "genes"

W = CREBBP.dat$expm.dat[[key]]$W
X = CREBBP.dat$expm.dat[[key]]$X
source.ids  = colnames(W)
feature.ids = rownames(X)
sample.ids  = rownames(W)

k = length(source.ids)
m = length(feature.ids)
n = length(sample.ids)

print(source.ids)

max_stds = 2

base.mdl       = readRDS(file.path(res.dir, paste0("base.mdl.",       key, ".rds")))
bulk.mdl       = readRDS(file.path(res.dir, paste0("bulk.mdl.",       key, ".rds")))
cibersortx.mdl = readRDS(file.path(res.dir, paste0("cibersortx.mdl.", key, ".rds")))
tca.mdl        = readRDS(file.path(res.dir, paste0("tca.mdl.",        key, ".rds")))
Unico.mdl       = readRDS(file.path(res.dir, paste0("Unico.mdl.",       key, ".rds")))
bMIND.scale.mdl= readRDS(file.path(res.dir, paste0("bMIND.scale.mdl.",key,  ".rds")))
bMIND.log.mdl  = readRDS(file.path(res.dir, paste0("bMIND.log.mdl.",  key,  ".rds")))

# now force none negativity for the log transformation later
base.mdl$Z.hat = none_neg_Z(base.mdl$Z.hat)
bulk.mdl$Z.hat = none_neg_Z(bulk.mdl$Z.hat)
cibersortx.mdl$Z.hat = none_neg_Z(cibersortx.mdl$Z.hat)
tca.mdl$Z.hat  = none_neg_Z(tca.mdl$Z.hat)
Unico.mdl$Z.hat = none_neg_Z(Unico.mdl$Z.hat)
bMIND.scale.mdl$Z.hat = none_neg_Z(bMIND.scale.mdl$Z.hat)
bMIND.log.mdl$Z.hat  = none_neg_Z(bMIND.log.mdl$Z.hat)

bMIND.scale.mdl$Z.hat[1,,]

bMIND.log.mdl$Z.hat[1,,]

test.source.id    = "Bcells"
main.method.names = c("CIBERSORTx", "TCA", "Unico", "bMIND.scale", "bMIND.log")
method.names      = c("Bulk", "Baseline", "CIBERSORTx", "TCA", "Unico", "bMIND.scale", "bMIND.log")

method.mdls       = list("Bulk"       = bulk.mdl, 
                         "Baseline"   = base.mdl, 
                         "CIBERSORTx" = cibersortx.mdl, 
                         "TCA"        = tca.mdl, 
                         "Unico"       = Unico.mdl,
                         "bMIND.scale"= bMIND.scale.mdl,
                         "bMIND.log"  = bMIND.log.mdl)

logDiff.df = matrix(NA, m, length(method.names))
rownames(logDiff.df) = feature.ids
colnames(logDiff.df) = method.names

for(method.name in method.names){
    Z.hat = method.mdls[[method.name]]$Z.hat

    log.mt.mean = log2(1 + rowMeans(Z.hat[test.source.id, feature.ids, samples.mt]))
    log.wt.mean = log2(1 + rowMeans(Z.hat[test.source.id, feature.ids, samples.wt])) 

    logDiff.df[, method.name] = log.mt.mean - log.wt.mean
}




# rowMedians(log)
colMedians(logDiff.df[inc.genes,])

colMedians(logDiff.df[dec.genes,])




length(logDiff.df[c(inc.genes, dec.genes), "Bulk"])

df = data.frame(diff = c(logDiff.df[c(inc.genes, dec.genes), "Bulk"],
                         logDiff.df[c(inc.genes, dec.genes), "Baseline"],
                         logDiff.df[c(inc.genes, dec.genes), "CIBERSORTx"],
                         logDiff.df[c(inc.genes, dec.genes), "TCA"],
                         logDiff.df[c(inc.genes, dec.genes), "Unico"], 
                         logDiff.df[c(inc.genes, dec.genes), "bMIND.scale"],
                         logDiff.df[c(inc.genes, dec.genes), "bMIND.log"]),

                models = c(rep("Bulk", m), rep("Baseline", m), rep("CIBERSORTx", m),
                           rep("TCA", m), rep("Unico", m), 
                           rep("bMIND.scale", m), rep("bMIND.log", m)),
                
                DE.type = rep(c(rep("Up Regulated", length(inc.genes)), rep("Down Regulated", length(dec.genes))), 
                              length(method.names)))

df$models <- factor(df$models, # Change ordering manually
                    levels = c("Bulk", "Baseline", "CIBERSORTx", "TCA", "Unico", "bMIND.scale", "bMIND.log"))

pvals.inc = matrix(1, length(method.names), length(method.names))
rownames(pvals.inc) = method.names
colnames(pvals.inc) = method.names
pvals.dec = matrix(1, length(method.names), length(method.names))
rownames(pvals.dec) = method.names
colnames(pvals.dec) = method.names

for (method.name1 in method.names){
    for (method.name2 in method.names){
        # does method2/col outperforms method1/row
        pvals.inc[method.name1, method.name2] = wilcox.test(logDiff.df[inc.genes, method.name1],   
                                                            logDiff.df[inc.genes, method.name2], paired = T, "greater")$p.value
        pvals.dec[method.name1, method.name2] = wilcox.test(logDiff.df[dec.genes, method.name1],   
                                                            logDiff.df[dec.genes, method.name2], paired = T, "less")$p.value
    }
}

pvals.inc 
pvals.dec



test.sample.ids = c(samples.mt, samples.wt)

mutation = matrix(0, length(test.sample.ids),1)
rownames(mutation) = test.sample.ids
colnames(mutation) = c("mutation.status")
mutation[samples.mt, ] = 1
mutation[samples.wt, ] = 0

method.beta.list = list()
method.log10p.list = list()

for (method.name in method.names){
    Z.hat = method.mdls[[method.name]]$Z.hat

    beta.mat = matrix(0, length(feature.ids), (3)) # 1 intercept + 1 celltype + 1 mutation status
    rownames(beta.mat) = feature.ids
    p.vals.mat = matrix(1, length(feature.ids), (3)) 
    rownames(p.vals.mat) = feature.ids

    for (feature.id in feature.ids){
        regression.df = data.frame(Z.hat         = log2(1 + Z.hat[test.source.id, feature.id, test.sample.ids]),
                                   W.test.source = W[test.sample.ids, test.source.id] ,
                                   mutation      = mutation)
        tryCatch({
            #normalize to unit variance 
            regression.df$Z.hat = regression.df$Z.hat/sd(regression.df$Z.hat)
            
            fit <- lm(Z.hat ~ ., data = regression.df)
            p.vals.mat[feature.id, ] = t(summary(fit)$coefficients)["Pr(>|t|)", , drop = F]
            colnames(p.vals.mat)     = rownames(summary(fit)$coefficients)
            
            beta.mat[feature.id, ] = t(summary(fit)$coefficients)["Estimate", , drop = F]
            colnames(beta.mat)     = rownames(summary(fit)$coefficients)
            
        }, error=function(cond){})
    }
   
    method.log10p.list[[method.name]] = -log10(p.vals.mat)
    method.beta.list[[method.name]] = beta.mat
}

regression.df

method.log10p.list$Unico

method.beta.list$Unico

reg.pvals.inc = matrix(1, length(method.names), length(method.names))
rownames(reg.pvals.inc) = method.names
colnames(reg.pvals.inc) = method.names
reg.pvals.dec = matrix(1, length(method.names), length(method.names))
rownames(reg.pvals.dec) = method.names
colnames(reg.pvals.dec) = method.names

for (method.name1 in method.names){
    for (method.name2 in method.names){
        # does method2/col outperforms method1/row
        reg.pvals.inc[method.name1, method.name2] = wilcox.test(method.beta.list[[method.name1]][inc.genes, "mutation.status"],   
                                                                method.beta.list[[method.name2]][inc.genes, "mutation.status"], paired = T, "greater")$p.value
        reg.pvals.dec[method.name1, method.name2] = wilcox.test(method.beta.list[[method.name1]][dec.genes, "mutation.status"],   
                                                                method.beta.list[[method.name2]][dec.genes, "mutation.status"], paired = T, "less")$p.value

    }
}

#showing bMIND after scaling s.t max is 50
#if directly divided by max(X). somehow the result look better on the up regulated set worse on the down regulated for bMIND

reg.pvals.inc
reg.pvals.dec

reg.pvals = matrix(1, length(method.names), length(method.names))
rownames(reg.pvals) = method.names
colnames(reg.pvals) = method.names

for (method.name1 in method.names){
    for (method.name2 in method.names){
        # does method2/col outperforms method1/row
        reg.pvals[method.name1, method.name2] = wilcox.test(c(method.beta.list[[method.name1]][inc.genes, "mutation.status"], -method.beta.list[[method.name1]][dec.genes, "mutation.status"]), 
                                                            c(method.beta.list[[method.name2]][inc.genes, "mutation.status"], -method.beta.list[[method.name2]][dec.genes, "mutation.status"]),
                                                            paired = T, "greater")$p.value
       
    }
}



CREBBP.eval = list(tensor = list(logDiff.df = logDiff.df,
                                 pvals.inc = pvals.inc,
                                 pvals.dec = pvals.dec),
                   regression = list(method.log10p.list = method.log10p.list,
                                     method.beta.list   = method.beta.list,
                                     reg.pvals.inc = reg.pvals.inc,
                                     reg.pvals.dec = reg.pvals.dec))

#constant and failed the regression
method.beta.list$CIBERSORTx[method.beta.list$CIBERSORTx[, "mutation.status"] == 0 , ]

saveRDS(CREBBP.eval, file.path(res.dir, "CREBBP.eval.rds"))


