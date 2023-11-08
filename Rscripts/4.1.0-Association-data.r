# The script downloads and runs a processing and filtering piepline of the Hannum and Liu DNA methylation datasets
# This script is based on the function prep_data in the following script: https://github.com/cozygene/CellTypeSpecificMethylationAnalysis/blob/master/R/diff_meth_age_sex.R
# data_dir - destination directory for RData files with the data
# hannum_smk_status_path - Path to hannum_smoking_status.txt that must be downloaded separately.
suppressPackageStartupMessages({
  library("TCA")
  library("EpiDISH")
  library("matrixStats")
  require("GEOquery")
  library("R.utils")
  require("data.table")
  library("pracma")
  library("stats")
  library("stringr")
  library("ggplot2")
  library("ggpubr")
  library("grid")
  library("hrbrthemes")
  library("extrafont")
  extrafont::loadfonts()
})

hannum_smk_status_path = "/u/project/halperin/johnsonc/Unico/Unico2023/Data/Methylation/Consistency/hannum_smoking_status.txt"
data_dir <- "/u/project/halperin/johnsonc/Unico/Unico2023/Data/Methylation/Consistency"

prep_data <- function(data_dir, hannum_smk_status_path){
  options(timeout=10000)
  Sys.setenv(VROOM_CONNECTION_SIZE=5000720)
  p <- 10000 # num of "control probes" for calculating control probe PCs
  nonspecific_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/nonspecific_probes.txt")[,1]
  XY_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/HumanMethylationSites_X_Y.txt")[,1]
  polymorphic_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/polymorphic_cpgs.txt")[,1]
  exclude <- union(nonspecific_probes,union(XY_probes,polymorphic_probes))

  
  file_name1 <- paste(data_dir,"liu.processed.RData",sep="/")
  if (!file.exists(file_name1)){
    # Download the Liu et al. data
    # note that the GEO series matrix of the Liu et al. data was updates on 2021-09-20; the new version now includes only ~140k CpGs. We therefore download the full processed version of the data.
    gse <- GEOquery::getGEO("liu", destdir = data_dir, GSEMatrix = TRUE)
    download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/suppl/GSE42861_processed_methylation_matrix.txt.gz", file.path(data_dir,"GSE42861_processed_methylation_matrix.txt.gz"))
    gunzip(file.path(data_dir,"GSE42861_processed_methylation_matrix.txt.gz"))
    X <- fread(file.path(data_dir,"GSE42861_processed_methylation_matrix.txt"))
    cpg_names <- X$ID_REF
    sample_names <- colnames(X)[2:(ncol(X)-1)] # the last column is pvals
    X <- as.matrix(X[,2:(ncol(X)-1)])
    colnames(X) <- sample_names
    rownames(X) <- cpg_names
    
    #str_split(Biobase::pData(gse[[1]])[["supplementary_file"]],"_")
    ids_map <- rownames(Biobase::pData(gse[[1]]))
    l <- str_split(Biobase::pData(gse[[1]])[["supplementary_file"]],"_")
    ids_map_names <- character(length(l))
    for (i in 1:length(l)){
      ids_map_names[i] <- paste(l[[i]][2],l[[i]][3],sep="_")
    }
    names(ids_map) <- ids_map_names
    X <- X[,names(ids_map)]
    colnames(X) <- as.character(ids_map)
    
    # smoking status; consider never-smokers [0], ex-smokers [1] and current smokers [2] (refer to occasional smokers as smokers)
    smk.liu <- Biobase::pData(gse[[1]])[["smoking status:ch1"]]
    # two samples have NA values for smoking; remove them from the analysis
    keep <- setdiff(1:length(smk.liu), which(smk.liu == "na"))
    smk.liu <- smk.liu[keep]
    smk.liu[which(smk.liu == "occasional")] <- "current"
    smk.liu <- 3-as.matrix(as.numeric(as.factor(smk.liu)))
    colnames(smk.liu) <- "smoking"
    # covariates
    cov.liu <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["disease state:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["gender:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["age:ch1"]])))
    colnames(cov.liu) <- c("disease", "gender", "age")
    cov.liu <- cov.liu[keep,]
    rownames(cov.liu) <- rownames(Biobase::pData(gse[[1]]))[keep]
    rownames(smk.liu) <- rownames(cov.liu)
    # methylaion data
    X.liu <- X[,rownames(smk.liu)]
    
    # estimate cell-type proportions using a reference-based approach 
    # methylation reference
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.liu <- EpiDISH::epidish(X.liu, ref)$estF
    
    liu <- list(X=X.liu,
                smk=smk.liu,
                cov=cov.liu,
                W=W.liu)

    # calculate PCs from low variance probes , to be treated as control probes (as in Lenhe et al. 2015); here, since we don't work with IDAT files and therefore don't have actual control probes, we use sites with the least variation in the data as control probes (as in Rahmani et al. 2019).
    site.variances <- matrixStats::rowVars(liu$X)
    names(site.variances) <- rownames(liu$X)
    low.var.sites <- names(head(sort(site.variances), p))
    low.var.pca <- prcomp(t(liu$X[low.var.sites,]), center=TRUE, scale=TRUE, rank=20)
    liu$ctrl_pcs <- low.var.pca$x

    # remove polymorphic or cross-reactive sites and non-autosomal sites; exclude low variance cpgs
    low_var_sites <- names(which(site.variances<0.001))
    keep <- setdiff(rownames(liu$X),union(low_var_sites, exclude))

    liu$X <- liu$X[keep,]
    liu$cov <- cbind(liu$cov,liu$smk)
    liu$smk <- NULL
    W_colnames <- c("Gran",setdiff(colnames(liu$W), c("Neutro","Eosino")))
    liu$W <- cbind(rowSums(liu$W[,c("Neutro","Eosino")]), liu$W[,setdiff(colnames(liu$W), c("Neutro","Eosino"))])
    colnames(liu$W) <- W_colnames
    save(liu, file=file_name1)
    rm(liu)
  }

   
  file_name2 <- paste(data_dir,"hannum.processed.RData",sep="/")
  if (!file.exists(file_name2)){
    # Download the Hannum et al. data 
    gse <- GEOquery::getGEO("GSE40279", destdir = data_dir, GSEMatrix = TRUE)
    # smoking status; consider never-smokers [0], ex-smokers [1] and current smokers [2] (refer to occasional smokers as smokers)
    smk.hannum <- read.table(hannum_smk_status_path, header = T, row.names = 1, sep=",")
    # remove samples with NA values for smoking
    keep <- setdiff(1:nrow(smk.hannum), which(is.na(smk.hannum)))
    keep.ids <- rownames(smk.hannum)[keep]
    smk.hannum <- smk.hannum[keep.ids,,drop=F]
    smk.hannum[smk.hannum=="Never"] <- 0
    smk.hannum[smk.hannum=="Past Smoker"] <- 1
    smk.hannum[smk.hannum=="Current Smoker"] <- 2
    smk.hannum <- as.numeric(smk.hannum$smk)
    names(smk.hannum) <- keep.ids
    # covariates
    cov.hannum <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["age (y):ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["gender:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["plate:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["ethnicity:ch1"]])))
    rownames(cov.hannum) <- Biobase::pData(gse[[1]])[["geo_accession"]]
    colnames(cov.hannum) <- c("age", "gender", "plate","ethnicity")
    cov.hannum <- cov.hannum[keep.ids,]
    # methylaion data
    X.hannum <- Biobase::exprs(gse[[1]])[,keep.ids]
    # remove sites with missing values 
    X.hannum <- X.hannum[rowSums(is.na(X.hannum)) == 0,]
    
    names(smk.hannum) <- colnames(X.hannum)
    rownames(cov.hannum) <- colnames(X.hannum)
    
    # estimate cell-type proportions using a reference-based approach 
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.hannum <- EpiDISH::epidish(X.hannum, ref)$estF
    
    hannum <- list(X=X.hannum,
                   smk=smk.hannum,
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
    cov_colnames <- c(colnames(hannum$cov),"smoking") 
    hannum$cov <- cbind(hannum$cov,hannum$smk)
    colnames(hannum$cov) <- cov_colnames
    hannum$smk <- NULL
    W_colnames <- c("Gran",setdiff(colnames(hannum$W), c("Neutro","Eosino")))
    hannum$W <- cbind(rowSums(hannum$W[,c("Neutro","Eosino")]), hannum$W[,setdiff(colnames(hannum$W), c("Neutro","Eosino"))])
    colnames(hannum$W) <- W_colnames
    save(hannum, file=file_name2)
    rm(hannum)
  }

  
  file_name3 <- paste(data_dir,"hannon1.processed.RData",sep="/")
  if (!file.exists(file_name3)){
    gse <- GEOquery::getGEO("GSE80417", destdir = data_dir, GSEMatrix = TRUE)

    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80417&format=file&file=GSE80417%5FnormalizedBetas%2Ecsv%2Egz", file.path(data_dir,"GSE80417_normalizedBetas.csv.gz"))
    gunzip(file.path(data_dir,"GSE80417_normalizedBetas.csv.gz"))
    X <- fread(file.path(data_dir,"GSE80417_normalizedBetas.csv"))
    cpg_names <- X$V1
    sample_names <- colnames(X)[2:ncol(X)]
    X <- as.matrix(X[,2:ncol(X)])
    colnames(X) <- sample_names
    rownames(X) <- cpg_names
    
    ids_map <- rownames(Biobase::pData(gse[[1]]))
    names(ids_map) <- Biobase::pData(gse[[1]])[["description"]]
    X <- X[,names(ids_map)]
    colnames(X) <- as.character(ids_map)

    # covariates
    cov.hannon1 <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["disease status:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["Sex:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["age:ch1"]])))
    colnames(cov.hannon1) <- c("disease", "gender", "age")
    rownames(cov.hannon1) <- rownames(Biobase::pData(gse[[1]]))
    X.hannon1 <- X
    
    # estimate cell-type proportions using a reference-based approach 
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.hannon1 <- EpiDISH::epidish(X.hannon1, ref)$estF
    
    hannon1 <- list(X=X.hannon1,
                cov=cov.hannon1,
                W=W.hannon1)

    # calculate PCs from low variance probes , to be treated as control probes (as in Lenhe et al. 2015); here, since we don't work with IDAT files and therefore don't have actual control probes, we use sites with the least variation in the data as control probes (as in Rahmani et al. 2019).
    site.variances <- matrixStats::rowVars(hannon1$X)
    names(site.variances) <- rownames(hannon1$X)
    low.var.sites <- names(head(sort(site.variances), p))
    low.var.pca <- prcomp(t(hannon1$X[low.var.sites,]), center=TRUE, scale=TRUE, rank=20)
    hannon1$ctrl_pcs <- low.var.pca$x

    # remove polymorphic or cross-reactive sites and non-autosomal sites; exclude low variance cpgs
    low_var_sites <- names(which(site.variances<0.001))
    keep <- setdiff(rownames(hannon1$X),union(low_var_sites, exclude))
    hannon1$X <- hannon1$X[keep,]
    
    W_colnames <- c("Gran",setdiff(colnames(hannon1$W), c("Neutro","Eosino")))
    hannon1$W <- cbind(rowSums(hannon1$W[,c("Neutro","Eosino")]), hannon1$W[,setdiff(colnames(hannon1$W), c("Neutro","Eosino"))])
    colnames(hannon1$W) <- W_colnames
    save(hannon1, file=file_name3)
    rm(hannon1)
  }


  file_name4 <- paste(data_dir,"hannon2.processed.RData",sep="/")
  if (!file.exists(file_name3)){
    gse <- GEOquery::getGEO("GSE84727", destdir = data_dir, GSEMatrix = TRUE)
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84727&format=file&file=GSE84727%5FnormalisedBetas%2Ecsv%2Egz", file.path(data_dir,"GSE84727_normalisedBetas.csv.gz"))
    gunzip(file.path(data_dir,"GSE84727_normalisedBetas.csv.gz"))
    X <- fread(file.path(data_dir,"GSE84727_normalisedBetas.csv"))

    cpg_names <- X$V1
    sample_names <- colnames(X)[2:ncol(X)]
    X <- as.matrix(X[,2:ncol(X)])
    colnames(X) <- sample_names
    rownames(X) <- cpg_names

    # covariates; remove samples with missing values
    cov.hannon2 <- cbind(Biobase::pData(gse[[1]])[["disease_status:ch1"]], Biobase::pData(gse[[1]])[["Sex:ch1"]], Biobase::pData(gse[[1]])[["age:ch1"]])
    rownames(cov.hannon2) <- rownames(Biobase::pData(gse[[1]]))
    cov.hannon2 <- cov.hannon2[rowSums(cov.hannon2 == "NA")==0,]
    sample_ids <- rownames(cov.hannon2)
    cov.hannon2 <- cbind(as.numeric(as.factor(cov.hannon2[,1])), as.numeric(as.factor(cov.hannon2[,2])), as.numeric(cov.hannon2[,3]))
    rownames(cov.hannon2) <- sample_ids
    colnames(cov.hannon2) <- c("disease", "gender", "age")
    
    ids_map_names <- c()
    l <- str_split(Biobase::pData(gse[[1]])[["characteristics_ch1"]],": ")
    for (i in 1:length(l)){
      ids_map_names[i] <- l[[i]][[2]]
    }
    ids_map <- ids_map_names
    names(ids_map) <- rownames(Biobase::pData(gse[[1]]))

    X.hannon2 <- X[,ids_map[rownames(cov.hannon2)]]
    #colnames(X.hannon2) <- rownames(cov.hannon2) #johnsonc

    # estimate cell-type proportions using a reference-based approach 
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.hannon2 <- EpiDISH::epidish(X.hannon2, ref)$estF
    
    hannon2 <- list(X=X.hannon2,
                cov=cov.hannon2,
                W=W.hannon2)

    # calculate PCs from low variance probes , to be treated as control probes (as in Lenhe et al. 2015); here, since we don't work with IDAT files and therefore don't have actual control probes, we use sites with the least variation in the data as control probes (as in Rahmani et al. 2019).
    site.variances <- matrixStats::rowVars(hannon2$X)
    names(site.variances) <- rownames(hannon2$X)
    low.var.sites <- names(head(sort(site.variances), p))
    low.var.pca <- prcomp(t(hannon2$X[low.var.sites,]), center=TRUE, scale=TRUE, rank=20)
    hannon2$ctrl_pcs <- low.var.pca$x

    # remove polymorphic or cross-reactive sites and non-autosomal sites; exclude low variance cpgs
    low_var_sites <- names(which(site.variances<0.001))
    keep <- setdiff(rownames(hannon2$X),union(low_var_sites, exclude))
    hannon2$X <- hannon2$X[keep,]
    
    W_colnames <- c("Gran",setdiff(colnames(hannon2$W), c("Neutro","Eosino")))
    hannon2$W <- cbind(rowSums(hannon2$W[,c("Neutro","Eosino")]), hannon2$W[,setdiff(colnames(hannon2$W), c("Neutro","Eosino"))])
    colnames(hannon2$W) <- W_colnames
    save(hannon2, file=file_name4)
    rm(hannon2)
  }
    
}

prep_data(data_dir)
