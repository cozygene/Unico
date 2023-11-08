library(TCA)
library(MIND)
library(EpiDISH)


setwd("/Users/johnsonchen/Documents/Research/Unico/Unico2023/Notebook/")
source("Unico.r")
source("analysis.utils.r")
set.seed(2023)
#################################################################################
m = 1000
permute.size = 10

#################################################################################
# path related
rds.dir  =  "../Data/Runtime/" 
res.dir  = "../Result/Runtime/mdl/"
#################################################################################



####################
cofigures = list(c(6,  m, 500), c(5,  m, 500), c(4,  m, 500), c(3,  m, 500),
								 c(6,  m, 100), c(6,  m, 250)) #automatically also get the 500  sample size result from previous

#cofigures = list(c(6,  m, 250)) #automatically also get the 500  sample size result from previous


for(t in 1:permute.size){
	for (cofigure in cofigures){
	#for (cofigure in cofigures[1:1]){
	k  = cofigure[1]
	m  = cofigure[2]
	n  = cofigure[3]
	
	
	set.seed(t)
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
	
		
	######################### Unico ############################
	res.file = file.path(res.dir, paste0("Unico.", config, ".rds"))
	log.file = file.path(res.dir, paste0("Unico.", config, ".log"))
	

	Unico.mdl = list()
	Unico.mdl$params.hat <- Unico(X = as.matrix(X), W = as.matrix(W), C1 = as.matrix(C1), C2 = as.matrix(C2), 
															parallel = T, num_cores = 8, 
															log_file = log.file, verbose = TRUE)
	
	#capping at interal scale
	Unico.mdl$params.hat$sigmas_hat  = cap_values(Unico.mdl$params.hat$sigmas_hat,  max.val = 10 ** (4), min.val = 10**(-4))
	#tensor 
	Unico.mdl$Z.hat  <- tensor(X = X,  
														W = W, 
														C1 = C1, 
														C2 = C2, 
														Unico.mdl$params.hat, parallel = F, 
														num_cores = 1, 
														log_file = log.file)
	#parametric
	Unico.mdl$params.hat = add_C1_C2_pvals_parametric(X = X, 
													 Unico.mdl = Unico.mdl$params.hat, slot_name = "parametric", 
													 parallel = T, num_cores = 8, 
													 config_file = NULL, log_file = log.file, debug = FALSE, verbose = TRUE)
		
	
		
	#Asymptotic
	Unico.mdl$params.hat = add_C1_C2_pvals_asymptotic(X = X, 
													 Unico.mdl = Unico.mdl$params.hat, slot_name = "asymptotic", 
													 parallel = T, num_cores = 8, 
													 config_file = NULL, log_file = log.file, debug = FALSE, verbose = TRUE)
	saveRDS(Unico.mdl, res.file)
	
		
	######################## TCA ############################
	res.file = file.path(res.dir, paste0("TCA.", config, ".rds"))
	log.file = file.path(res.dir, paste0("TCA.", config, ".log"))
	
	tca.mdl = list()
	tca.mdl$params.hat <- tca(X  = as.matrix(X), 
							W  = as.matrix(W), 
							C1 = as.matrix(C1),  
							C2 = as.matrix(C2), verbose = TRUE, parallel = TRUE, num_cores = 8, log_file=log.file)
	
	tca.mdl$params.hat$sigmas_hat  = cap_values(tca.mdl$params.hat$sigmas_hat,  max.val = 10 ** (4), min.val = 10**(-4))
	tca.mdl$Z.hat  <- TCA::tensor(X = as.matrix(X), 
								tca.mdl$params.hat,
								parallel = FALSE, 
								num_cores = 1, 
								log_file=log.file)
	saveRDS(tca.mdl, res.file)
	
	
	######################## bMIND.tensor ############################
	res.file = file.path(res.dir, paste0("bMIND.tensor.", config, ".rds"))

	#tensor
	start.t = Sys.time()
	bMIND.mdl = bMIND(X, W, ncore = 8) 
	end.t = Sys.time()
	bMIND.mdl$tensor.time = difftime(end.t, start.t, units = "secs")
	
	#tensor assoc doesnt have parallel
	start.t = Sys.time()
	y =  as.numeric(as.factor(C1[,"gender"]))-1 # used to be 1,2 coded, subtract 1 to be 0,1 coded
	covariate = cbind(C1[,-which(colnames(C1) == "gender")], C2)
	bMIND.mdl =  c(bMIND.mdl, test(bMIND.mdl$A, y, covariate))
	end.t = Sys.time()
	bMIND.mdl$tensor.assoc.time = difftime(end.t, start.t, units = "secs")
	
	saveRDS(bMIND.mdl, res.file)
	print('bMIND')
	print(bMIND.mdl$tensor.time)
	print(bMIND.mdl$tensor.assoc.time)
	
	########################## cellDMC ############################
	res.file = file.path(res.dir, paste0("CellDMC.", config, ".rds"))
	CellDMC.mdl = list()
	y =  as.numeric(as.factor(C1[,"gender"]))-1 # used to be 1,2 coded, subtract 1 to be 0,1 coded
	covariate = cbind(C1[,-which(colnames(C1) == "gender")], C2)
	
	start.t = Sys.time()
	CellDMC.mdl$res <- CellDMC(X, y, W, cov.mod = covariate, mc.cores = 8)
	end.t = Sys.time()
	
	CellDMC.mdl$assoc.time = difftime(end.t, start.t, units = "secs")
	saveRDS(CellDMC.mdl, res.file)
	print('CellDMC')
	print(CellDMC.mdl$assoc.time)
	}
}





#confirm that new implementation and the old has the same p val manually and load the qqcompare 
#conform runtime is behaving as expected 
#gammas_hat_res_old = Unico.mdl$params.hat$asymptotic$gammas_hat
#gammas_hat_pvals_res_old = Unico.mdl$params.hat$asymptotic$gammas_hat_pvals

#all(gammas_hat_res_old == Unico.mdl$params.hat$asymptotic$gammas_hat)
#all(abs(gammas_hat_pvals_res_old - Unico.mdl$params.hat$asymptotic$gammas_hat_pvals)< 0.00000001)
#all(round(log(gammas_hat_pvals_res_old), 1)== round(log(Unico.mdl$params.hat$asymptotic$gammas_hat_pvals),1))


# there is an issue with this operation:
# inv_res %*% t(S) %*% diag(as.vector(Q_star_j)) %*% S %*% t(inv_res)
# the diag matrix in the middle is n by n 

#n = 1000
#m = 20
#s_len = 20 

#S = matrix(rnorm(n * s_len), n, s_len)
#Q_vec  = rnorm(n)
#Q_diag = diag(Q_vec)

#T = 1000
#start.time <- Sys.time()
#for (t in 1:T){
#	tmp = t(S) %*% Q_diag %*% S
#}
#end.time <- Sys.time()
#end.time - start.time
 
#tmp[1:4,1:4]

#T = 1000
#start.time <- Sys.time()
#for (t in 1:T){
#	tmp = (t(S) * repmat(t(as.matrix(Q_vec)),  s_len,  1)) %*% S
#}
#end.time <- Sys.time()
#end.time - start.time
#tmp[1:4,1:4]
