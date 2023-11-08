set.seed(2023)

pheno = "age"
#pheno = "gender"

model_tau     = FALSE


data.names = c("liu", "hannum", "hannon1", "hannon2")
source.ids = c("Gran", "CD4T", "CD8T", "Mono", "B", "NK")

results_dir = "/u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/XY/"
if (!file.exists(results_dir)){
     print("check results_dir")
}

pval.list  = list()
coeff.list = list()

for (data.name in data.names){
    
    
    Unico.mdl      = readRDS(file.path(results_dir, data.name, paste0("Unico.mdl.rds")))
    tca.mdl       = readRDS(file.path(results_dir, data.name,'tca.mdl.rds'))
    celldmc.mdl  = readRDS(file.path(results_dir,  data.name,'CellDMC.mdl.rds'))
    base.mdl = readRDS(file.path(results_dir,  data.name,'base.mdl.rds'))
    
    pval  = list()
    coeff = list()

    #Unico
    for (slot.name in c("parametric")){

        pval [[paste0("Unico.",slot.name,".marginal")]] = Unico.mdl$params.hat[[slot.name]]$gammas_hat_pvals[,paste(source.ids, pheno, sep=".")]
        coeff[[paste0("Unico.",slot.name,".marginal")]] = Unico.mdl$params.hat[[slot.name]]$gammas_hat[,paste(source.ids, pheno, sep=".")]

        pval [[paste0("Unico.",slot.name,".joint")]] = Unico.mdl$params.hat[[slot.name]]$gammas_hat_pvals.joint[, pheno, drop = F]
        #all set to 1 for the coef
        coeff[[paste0("Unico.",slot.name,".joint")]] = (pval [[paste0("Unico.",slot.name,".joint")]] > -1) * 1
    }

    #TCA
    for (slot.name in c("parametric")){

        pval [[paste0("TCA.",slot.name,".marginal")]] = tca.mdl$params.hat[[slot.name]]$gammas_hat_pvals[,paste(source.ids, pheno, sep=".")]
        coeff[[paste0("TCA.",slot.name,".marginal")]] = tca.mdl$params.hat[[slot.name]]$gammas_hat[,paste(source.ids, pheno, sep=".")]

        pval [[paste0("TCA.",slot.name,".joint")]] = tca.mdl$params.hat[[slot.name]]$gammas_hat_pvals.joint[, pheno, drop = F]
        #all set to 1 for the coef
        coeff[[paste0("TCA.",slot.name,".joint")]] = (pval [[paste0("TCA.",slot.name,".joint")]] > -1) * 1
    } 

    #CellDMC
    for (slot.name in c("parametric")){
        slot.res = celldmc.mdl[[pheno]][[slot.name]]
        pval [[paste0("CellDMC.",slot.name,".marginal")]] = cbind(slot.res$marg.res$Gran[,"p"], 
                                                                  slot.res$marg.res$CD4T[,"p"],
                                                                  slot.res$marg.res$CD8T[,"p"],
                                                                  slot.res$marg.res$Mono[,"p"],
                                                                  slot.res$marg.res$B[,"p"],
                                                                  slot.res$marg.res$NK[,"p"])
        colnames(pval [[paste0("CellDMC.",slot.name,".marginal")]]) = paste0(source.ids, ".", pheno)

        coeff[[paste0("CellDMC.",slot.name,".marginal")]] = cbind(slot.res$marg.res$Gran[,"Estimate"], 
                                                                  slot.res$marg.res$CD4T[,"Estimate"],
                                                                  slot.res$marg.res$CD8T[,"Estimate"],
                                                                  slot.res$marg.res$Mono[,"Estimate"],
                                                                  slot.res$marg.res$B[,"Estimate"],
                                                                  slot.res$marg.res$NK[,"Estimate"])
        colnames(coeff [[paste0("CellDMC.",slot.name,".marginal")]]) = paste0(source.ids, ".", pheno)


        pval [[paste0("CellDMC.",slot.name,".joint")]] = slot.res$joint.res
        colnames(pval [[paste0("CellDMC.",slot.name,".joint")]]) = c(pheno)
        #all set to 1 for the coef
        coeff[[paste0("CellDMC.",slot.name,".joint")]] = (pval [[paste0("CellDMC.",slot.name,".joint")]] > -1) * 1
    }
    
    for (slot.name in c("parametric")){
        slot.res = base.mdl[[slot.name]]
        pval [[paste0("Baseline.",slot.name,".marginal")]] = t(slot.res$marg.pvals[source.ids,,pheno])
        colnames(pval [[paste0("Baseline.",slot.name,".marginal")]]) = paste0(source.ids, ".", pheno)

        coeff [[paste0("Baseline.",slot.name,".marginal")]] = t(slot.res$marg.coefs[source.ids,,pheno])
        colnames(coeff [[paste0("Baseline.",slot.name,".marginal")]]) = paste0(source.ids, ".", pheno)

        

        pval [[paste0("Baseline.",slot.name,".joint")]] = slot.res$joint.pvals[,pheno, drop = F]
        coeff[[paste0("Baseline.",slot.name,".joint")]] = slot.res$joint.coefs[,pheno, drop = F]
        
    } 

    pval.list[[data.name]]  = pval
    coeff.list[[data.name]] = coeff
}

str(pval.list)

str(pval.list$liu)

pval.list.file  = file.path(results_dir, paste0(pheno,".pval.list"))
coeff.list.file = file.path(results_dir, paste0(pheno,".coeff.list"))

saveRDS(pval.list,   pval.list.file)
saveRDS(coeff.list,  coeff.list.file)


