library("GEOquery")
set.seed(2023)

data.dir = "../Data/Methylation/Purified-Reinius/"

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
gse <- GEOquery::getGEO("GSE40279", destdir = data.dir, GSEMatrix = TRUE)

str(gse)

#general bad cpgs
p <- 10000 # num of "control probes" for calculating control probe PCs
nonspecific_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/nonspecific_probes.txt")[,1]
XY_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/HumanMethylationSites_X_Y.txt")[,1]
polymorphic_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/polymorphic_cpgs.txt")[,1]
exclude <- union(nonspecific_probes,union(XY_probes,polymorphic_probes))


# covariates
cov.hannum <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["age (y):ch1"]])), 
                    as.numeric(as.factor(Biobase::pData(gse[[1]])[["gender:ch1"]])),
                    as.numeric(as.factor(Biobase::pData(gse[[1]])[["plate:ch1"]])), 
                    as.character(as.factor(Biobase::pData(gse[[1]])[["ethnicity:ch1"]])))

rownames(cov.hannum) <- Biobase::pData(gse[[1]])[["geo_accession"]]
colnames(cov.hannum) <- c("age", "gender", "plate","ethnicity")

# methylaion data
X.hannum <- Biobase::exprs(gse[[1]])
# remove sites with missing values 
X.hannum <- X.hannum[rowSums(is.na(X.hannum)) == 0,]
X.hannum = X.hannum[, rownames(cov.hannum)] #make sure they are the same order

# estimate cell-type proportions using a reference-based approach 
ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
W.hannum <- EpiDISH::epidish(X.hannum, ref)$estF
hannum <- list(X=X.hannum,
               cov=cov.hannum,
               W=W.hannum)

# calculate PCs from low variance probes, to be treated as control probes (as in Lenhe et al. 2015); here, since we don't work with IDAT files and therefore don't have actual control probes, we use sites with the least variation in the data as control probes (as in Rahmani et al. 2019).
site.variances <- matrixStats::rowVars(hannum$X)
names(site.variances) <- rownames(hannum$X)

low.var.sites <- names(head(sort(site.variances), p))
low.var.pca <- prcomp(t(hannum$X[low.var.sites,]), center=TRUE, scale=TRUE, rank=20)
hannum$ctrl_pcs <- low.var.pca$x

# remove polymorphic or cross-reactive sites and non-autosomal sites; exclude low variance cpgs
low_var_sites <- names(which(site.variances<0.001))
keep <- setdiff(rownames(hannum$X),union(low_var_sites, exclude))

hannum$X <- hannum$X[keep,]
W_colnames <- c("Gran",setdiff(colnames(hannum$W), c("Neutro","Eosino")))
hannum$W <- cbind(rowSums(hannum$W[,c("Neutro","Eosino")]), hannum$W[,setdiff(colnames(hannum$W), c("Neutro","Eosino"))])
colnames(hannum$W) <- W_colnames
hannum$W = hannum$W[, order(colMeans(hannum$W), decreasing = T)]

str(hannum)
#length(unique(hannum$cov[,"plate"])) returns 9 

hannum$W

saveRDS(hannum, file.path(data.dir, "hannum.rds"))

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
gse <- GEOquery::getGEO("GSE35069", destdir = data.dir, GSEMatrix = TRUE)

# covars
cov.reinius <- cbind(as.character(Biobase::pData(gse[[1]])[["title"]]),
                     as.character(Biobase::pData(gse[[1]])[["source_name_ch1"]]),
                     as.character(Biobase::pData(gse[[1]])[["gender:ch1"]]))
rownames(cov.reinius) <- Biobase::pData(gse[[1]])[["geo_accession"]]
colnames(cov.reinius) <- c("sample.id", "tissue", "gender")
cov.reinius[,"sample.id"] = vapply(1:nrow(cov.reinius), 
                                   function(i) paste0("sample.",strsplit(cov.reinius[i,"sample.id"], "_")[[1]][2]), character(1))

#same name used in EpiDISH
cov.reinius[cov.reinius[,"tissue"] == "Granulocytes",    "tissue"] = "Gran"
cov.reinius[cov.reinius[,"tissue"] == "CD4+ T cells",    "tissue"] = "CD4T"
cov.reinius[cov.reinius[,"tissue"] == "CD8+ T cells",    "tissue"] = "CD8T"
cov.reinius[cov.reinius[,"tissue"] == "CD14+ Monocytes", "tissue"] = "Mono"
cov.reinius[cov.reinius[,"tissue"] == "CD19+ B cells",   "tissue"] = "B"
cov.reinius[cov.reinius[,"tissue"] == "CD56+ NK cells",  "tissue"] = "NK"

# methylaion data
X.reinius <- Biobase::exprs(gse[[1]])
X.reinius = X.reinius[, rownames(cov.reinius)] #same order as covars

# remove sites with missing values 
X.reinius <- X.reinius[rowSums(is.na(X.reinius)) == 0,]

source.ids  = c("Gran", "CD4T", "CD8T", "Mono", "B", "NK")
feature.ids = rownames(X.reinius)
sample.ids  = sort(unique(cov.reinius[,"sample.id"]))

k = length(source.ids)
m = length(feature.ids)
n = length(sample.ids)

Z = array(0, c(k,m,n))
dimnames(Z)[[1]] = source.ids
dimnames(Z)[[2]] = feature.ids
dimnames(Z)[[3]] = sample.ids

for (source.id in source.ids){
    for(sample.id in sample.ids){
        gse.sample.id = rownames(cov.reinius[cov.reinius[,"tissue"] == source.id & cov.reinius[,"sample.id"] == sample.id, ,drop = F])
        Z[source.id,,sample.id] = X.reinius[feature.ids, gse.sample.id] 
    }
}

# just the mixture X 
X = matrix(0,m,n)
rownames(X) = feature.ids
colnames(X) = sample.ids

gse.sample.ids = c()
for(sample.id in sample.ids){
    gse.sample.id = rownames(cov.reinius[cov.reinius[,"tissue"] == "Whole blood" & cov.reinius[,"sample.id"] == sample.id, ,drop = F])
    gse.sample.ids = c(gse.sample.ids, gse.sample.id)
    X[feature.ids,sample.id] = X.reinius[feature.ids, gse.sample.id] 
}

# estimate cell-type proportions using a reference-based approach 
ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
W.reinius <- EpiDISH::epidish(X, ref)$estF
W_colnames <- c("Gran",setdiff(colnames(W.reinius), c("Neutro","Eosino")))
W.reinius <- cbind(rowSums(W.reinius[,c("Neutro","Eosino")]), W.reinius[,setdiff(colnames(W.reinius), c("Neutro","Eosino"))])
colnames(W.reinius) <- W_colnames
W.reinius = W.reinius[, order(colMeans(W.reinius), decreasing = T)]

reinius <- list(X   = X,
                cov = cov.reinius[gse.sample.ids, "gender", drop = F],
                W   = W.reinius,
                Z   = Z)

str(reinius)

reinius$W

saveRDS(reinius, file.path(data.dir, "reinius.rds"))
