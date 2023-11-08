library("rlang") #duplicate
set.seed(2023)

data.names  = c("liu", "hannum", "hannon1", "hannon2")
source.ids = c("Gran", "CD4T", "CD8T", "Mono", "B", "NK")

results_dir = "/u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/XY/"
if (!file.exists(results_dir)){
     print("check results_dir")
}

pheno = "age"

pval.list.file  = file.path(results_dir, paste0(pheno,".pval.list"))
coeff.list.file = file.path(results_dir, paste0(pheno,".coeff.list"))

pval.list  = readRDS(pval.list.file)
coeff.list = readRDS(coeff.list.file)

str(pval.list$liu)

get_confusion_df = function(dis_pvals_mat, val_pvals_mat,
                            dis_coeffs_mat, val_coeffs_mat,
                            sig_thr,  dis_sig_mode,  val_sig_mode,
                            null_thr, dis_null_mode, val_null_mode, verbose = F){
    
    message(paste("discovery shape: ", dim(dis_pvals_mat)[1], dim(dis_pvals_mat)[2]))
    message(paste("validation shape: ", dim(val_pvals_mat)[1], dim(val_pvals_mat)[2])) 
     
    keep.probs     = intersect(rownames(dis_pvals_mat), rownames(val_pvals_mat))
    dis_pvals_mat  = duplicate(dis_pvals_mat)[keep.probs, , drop = F]
    dis_coeffs_mat = duplicate(dis_coeffs_mat)[keep.probs, , drop = F]
    val_pvals_mat  = duplicate(val_pvals_mat)[keep.probs, , drop = F]
    val_coeffs_mat = duplicate(val_coeffs_mat)[keep.probs, , drop = F]
    
    message("after")
    message(paste("discovery shape: ", dim(dis_pvals_mat)[1], dim(dis_pvals_mat)[2]))
    message(paste("validation shape: ", dim(val_pvals_mat)[1], dim(val_pvals_mat)[2])) 
    message("###################")
    
    #Positive
    if (dis_sig_mode == "strict"){
        dis_sig_thr = sig_thr/(dim(dis_pvals_mat)[1] * dim(dis_pvals_mat)[2])
    }else if (dis_sig_mode == "unadjust"){
        dis_sig_thr = sig_thr
    }else{
        print("error")
    }
    if(verbose){print(paste0("dis_sig_thr: ",  dis_sig_thr))}

    dis_sig_mask = dis_pvals_mat < dis_sig_thr
    n_ds = sum(dis_sig_mask)
    if(verbose){print(paste0("n_ds: ",  n_ds))}

    if (val_sig_mode == "strict"){
        val_sig_thr = sig_thr/(dim(val_pvals_mat)[1] * dim(val_pvals_mat)[2])
    }else if (val_sig_mode == "loose"){
        val_sig_thr = sig_thr/n_ds
    }else if (val_sig_mode == "unadjust"){
        val_sig_thr = sig_thr
    }else{
        print("error")
    }
    if(verbose){print(paste0("val_sig_thr: ",  val_sig_thr))}
    
    val_sig_mask_pvals = val_pvals_mat[dis_sig_mask] < val_sig_thr
    print("val_sig_mask_pvals")
    print(sum(val_sig_mask_pvals))
    
    val_sig_mask_ceoffs = sign(val_coeffs_mat[dis_sig_mask]) == sign(dis_coeffs_mat[dis_sig_mask]) 
    print("val_sig_mask_ceoffs")
    print(sum(val_sig_mask_ceoffs))
    val_sig_mask = val_sig_mask_pvals & val_sig_mask_ceoffs
    TP = sum(val_sig_mask)
    FN = sum(!val_sig_mask)


    # Negative
    if (dis_null_mode == "strict"){
        print("error")
    }else if (dis_null_mode == "unadjust"){
        dis_null_thr = null_thr
    }else{
        print("error")
    }
    if(verbose){print(paste0("dis_null_thr: ",  dis_null_thr))}

    dis_null_mask = dis_pvals_mat > dis_null_thr

    n_dn = sum(dis_null_mask)
    if(verbose){print(paste0("n_dn: ",  n_dn))}

    if (val_null_mode == "strict"){
        print("error")
    }else if (val_null_mode == "loose"){
        val_null_thr = sig_thr/n_dn
    }else if (val_null_mode == "unadjust"){
        val_null_thr = sig_thr
    }else{
        print("error")
    }
    if(verbose){print(paste0("val_null_thr: ",  val_null_thr))}

    val_null_mask = val_pvals_mat[dis_null_mask] > val_null_thr 

    FP = sum(!val_null_mask)
    TN = sum(val_null_mask)   
    
    recall    = TP/(TP + FN)
    precision = TP/(TP + FP)
    specificity = TN/(TN + FP)
    F1        = 1/(0.5 * (1/recall + 1/precision))
    accuracy = (TP + TN)/(TP + FP + TN + FN)
    TP = as.numeric(TP)
    TN = as.numeric(TN)
    FP = as.numeric(FP)
    FN = as.numeric(FN)
    MCC = as.numeric((TP*TN) - (FP*FN))/ sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    
    confusion_df = data.frame(TP = TP, 
                              FN = FN, 
                              FP = FP, 
                              TN = TN, 
                              recall = recall, 
                              precision = precision, 
                              F1 = F1,
                              specificity = specificity,
                              accuracy = accuracy,
                              MCC = MCC)
    return(confusion_df)
}

#loop through each method
get_confusion_dfs =function(methods, 
                            dis_pvals_mats, val_pvals_mats,
                            dis_coeffs_mats, val_coeffs_mats,
                            sig_thr,  dis_sig_mode,  val_sig_mode,
                            null_thr, dis_null_mode, val_null_mode, verbose = F){
    confusion_dfs = c()
    for (method in methods){
        print(method)
        confusion_df = get_confusion_df(dis_pvals_mats[[method]], val_pvals_mats[[method]],
                                        dis_coeffs_mats[[method]], val_coeffs_mats[[method]],
                                        sig_thr,  dis_sig_mode,  val_sig_mode,
                            null_thr, dis_null_mode, val_null_mode, verbose = verbose)
        rownames(confusion_df) = c(method)
        confusion_dfs = rbind(confusion_dfs, confusion_df)
    }
    return(confusion_dfs)
}

#loop through both direction 
get_bidir_confusion_dfs = function(methods, 
                                   pvals_mats, coeffs_mats, 
                                   sig_thr,  dis_sig_mode,  val_sig_mode,
                                   null_thr, dis_null_mode, val_null_mode, verbose = F){
    
    
    
    bidir_confusion_dfs = list()

    bidir_confusion_dfs[["forward"]] = get_confusion_dfs(methods,
                                                         #dataset 1.      # dataset 2
                                                         pvals_mats[[1]], pvals_mats[[2]],
                                                         coeffs_mats[[1]], coeffs_mats[[2]],
                                                         sig_thr,  dis_sig_mode,  val_sig_mode,
                                                         null_thr, dis_null_mode, val_null_mode, verbose = verbose)
    bidir_confusion_dfs[["backward"]] = get_confusion_dfs(methods,
                                                          pvals_mats[[2]], pvals_mats[[1]],
                                                          coeffs_mats[[2]], coeffs_mats[[1]],
                                                          sig_thr,  dis_sig_mode,  val_sig_mode,
                                                          null_thr, dis_null_mode, val_null_mode, verbose = verbose)
    return(bidir_confusion_dfs)
}

methods = c()
for (data.name in names(pval.list)){
    methods = c(methods, list(names(pval.list[[data.name]])))
}
methods = Reduce(intersect, methods)
print("methods with result on all 4 datasets")
print(methods)

bidir_confusion_dfs_list = list()
for (i in 1:length(data.names)){
    for(j in i:length(data.names)){
        print(paste0(data.names[i],".", data.names[j]))
        bidir_confusion_dfs = get_bidir_confusion_dfs(methods = methods,
                                                      pvals_mats = list(pval.list[[i]], pval.list[[j]]),
                                                      coeffs_mats = list(coeff.list[[i]], coeff.list[[j]]),
                                                      sig_thr  = 0.05, dis_sig_mode  = "strict", val_sig_mode = "loose",
                                                      null_thr = 0.95, dis_null_mode = "unadjust", val_null_mode = "loose", verbose = T)
        bidir_confusion_dfs_list[[paste0(data.names[i],".", data.names[j])]] = bidir_confusion_dfs
    }
}

str(pval.list)



#biggest sample
bidir_confusion_dfs_list$'liu.liu'$forward

#smallest sample
bidir_confusion_dfs_list$'hannum.hannum'$forward

bidir_confusion_dfs_list$'liu.hannum'

get_meta_mats = function(bidir_confusion_dfs, key){
    
    meta.mats = list()
    for (method in methods){
        meta.mats[[method]] = matrix(0,length(data.names),length(data.names))
        rownames(meta.mats[[method]]) = data.names
        colnames(meta.mats[[method]]) = data.names
    }

    for (i in 1:length(data.names)){
        for(j in i:length(data.names)){
            #print(paste0(data.names[i],".", data.names[j]))
            bidir_confusion_dfs = bidir_confusion_dfs_list[[paste0(data.names[i],".", data.names[j])]]
            for (method in methods){
                meta.mats[[method]][i,j] = bidir_confusion_dfs$forward[method, key]
                meta.mats[[method]][j,i] = bidir_confusion_dfs$backward[method, key]
            }
        }
    }
    return(meta.mats)
}

meta.mats.F1  = get_meta_mats(bidir_confusion_dfs, "F1")
meta.mats.MCC = get_meta_mats(bidir_confusion_dfs, "MCC")

saveRDS(bidir_confusion_dfs_list, file.path(results_dir, paste0(pheno,".bidir_confusion_dfs_list",".rds")))
saveRDS(meta.mats.F1,  file.path(results_dir, paste0(pheno,".meta.mats.F1",".rds")))
saveRDS(meta.mats.MCC, file.path(results_dir, paste0(pheno,".meta.mats.MCC",".rds")))

