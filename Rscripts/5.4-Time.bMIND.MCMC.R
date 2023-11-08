library(MIND)
set.seed(2023)
#################################################################################
m = 1000
permute.size = 10

#################################################################################
# path related
rds.dir  =  "../Data/Runtime/" 
res.dir  = "../Result/Runtime/mdl/"
#################################################################################

cofigures = list(c(6,  m, 500), c(5,  m, 500), c(4,  m, 500), c(3,  m, 500),
								 c(6,  m, 100), c(6,  m, 250)) #automatically also get the 500  sample size result from previous

#################################################################################
# pval for one gene in case some cell types have missing pval
get_pval = function(pval, cell_type, K) {
	pval0 = rep(NA, K)
	names(pval0) = cell_type
	names = intersect(names(pval), cell_type)
	pval0[names] = pval[names]
	return(pval0)
}


test = function(A, y, covariate = NULL) {
	
	if(dim(A)[3] != length(y)) print('CTS estimates and y have different length')
	if(!is.null(covariate)) if(dim(A)[3] != nrow(covariate)) print('CTS estimates and covariate have different number of samples/subjects') else {
		if(!is.null(rownames(covariate)) & any(rownames(covariate) != dimnames(A)[[3]])) covariate = covariate[dimnames(A)[[3]],]
	}
	
	K = ncol(A)
	cell_type = colnames(A)
	if(is.null(covariate)) pval = apply(A, 1, function(x) {
		pval = coef(summary(glm(y ~ ., data = data.frame(t(x)), family = 'binomial')))[,4]
		return(get_pval(pval, cell_type, K))
	}) else
		pval = apply(A, 1, function(x) {
			pval = coef(summary(glm(y ~ ., data = data.frame(t(x), covariate), family = 'binomial')))[,4]
			return(get_pval(pval, cell_type, K))
		})
	
	qval = pval2qval(pval, A, y, covariate)
	# rownames(qval) = rownames(pval) = substring(rownames(pval), 5)
	return(list(manova.pvals = qval, glm.pvals = pval))
}


# MANOVA; pval: K x ngene
pval2qval = function(pval, A, y, covariate = NULL) {
	ng = nrow(A)
	# pval for each gene
	if(is.null(covariate)){
		pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,,]) ~ y))$stats[1, "Pr(>F)"], silent = T))   
	}else{
		pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,,]) ~ y + covariate))$stats[1, "Pr(>F)"], silent = T))
	}
	
	#pval1 is m by 1 vector from the MANOVA
	#this step is filtering pval on features that has valid result from the MANOVA
	pval = pval[,!is.na(as.numeric(pval1))]
	
	# qval1 is still working on MANOVA stats
	pval1 = na.omit(as.numeric(pval1))
	qval1 = p.adjust(pval1, 'fdr')
	
	# the final qval result. first just copy the shape on the pval: K by M
	qval.adjust   = pval
	qval.unadjust = pval
	K = ncol(A)
	
	#per feature
	for(i in 1:ncol(pval)) {
		#initialize all celltypes to 1
		qval.adjust[,i] = 1
		qval.unadjust[,i] = 1
		#if the pval is decently small by 0.05/K. copy over the on MANOVA stat for the gene to the final qval. put it in the celltype with smallest p
		if(min(pval[,i], na.rm = T) < .05/K) qval.adjust[,i][which.min(pval[,i])] = qval1[i]
		if(min(pval[,i], na.rm = T) < .05/K) qval.unadjust[,i][which.min(pval[,i])] = pval1[i]
	}
	
	joint.pvals = as.matrix(pval1)
	rownames(joint.pvals) = colnames(qval.adjust)
	colnames(joint.pvals) = c("joint")
	
	return(list(joint.unadjust.pvals    = joint.pvals, 
							marginal.adjust.pvals   = qval.adjust, 
							marginal.unadjust.pvals = qval.unadjust))
}
#################################################################################
#for(t in 1:1){
for(t in 1:permute.size){
	set.seed(t)
	print(paste0("working on iteration: ", t))
	
	for (cofigure in cofigures){
		k  = cofigure[1]
		m  = cofigure[2]
		n  = cofigure[3]

		config = paste(paste0("k_",k),
									 paste0("m_",m),
									 paste0("n_",n),
									 paste0("t_",t), sep = ".")
		print(paste0("working on simulated data: ", config))
		
		sim.data = readRDS(file.path(rds.dir, paste0("sim.", config, ".rds")))
		X = sim.data$X
		W = sim.data$W
		C1 = sim.data$C1
		C2 = sim.data$C2
	
	
		######################### bMIND ############################
		res.file = file.path(res.dir, paste0("bMIND.mcmc", config, ".rds"))
		#mcmc assoc
		start.t = Sys.time()
		bMIND.mdl = list()
		y =  as.numeric(as.factor(C1[,"gender"]))-1 # used to be 1,2 coded, subtract 1 to be 0,1 coded
		covariate = cbind(C1[,-which(colnames(C1) == "gender")], C2)
		deconv = bmind_de(X, W, y = y, # used to be 1,2 coded, subtract 1 to be 0,1 coded
											covariate      = covariate, 
											covariate_bulk = colnames(C1[,-which(colnames(C1) == "gender")]), 
											covariate_cts  = colnames(C2), 
											np = T, ncore = 30, max_samp = 10**7 * 2)
		end.t = Sys.time()
		bMIND.mdl$mcmc.deconv = deconv
		bMIND.mdl$mcmc.assoc.time = difftime(end.t, start.t, units = "secs")
		print(bMIND.mdl$mcmc.assoc.time)
		saveRDS(bMIND.mdl, res.file)
	}
}

