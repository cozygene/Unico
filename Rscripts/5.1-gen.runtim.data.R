library("testit")
library("data.table") 
library("matrixStats")
set.seed(2023)
#################################################################################
# path related
project.dir = "/Users/johnsonchen/Documents/Research/Unico/Unico2023"
rds.dir     = file.path(project.dir, paste0("Data/Runtime/"))
txt.dir     = file.path(project.dir, paste0("Result/Runtime/cibersortx"))
#################################################################################


load(file.path(project.dir, paste0("Data/Methylation/Consistency/hannum.processed.RData")))
X  = hannum$X
C1 = hannum$cov
C2 = hannum$ctrl_pcs[,c(1:10)]

#################################################################################


#loop through random shuffle fo data 
#permute.size = 10
#m = 1000

#loop through celltype numbers 3, 4, 5, 6
#loop through n 100, 250, 500
#cofigures = list(c(6,  m, 500), c(5,  m, 500), c(4,  m, 500), c(3,  m, 500),
#								 c(6,  m, 100), c(6,  m, 250)) #automatically also get the 500  sample size result from previous

permute.size = 1
m = 500
cofigures = list(c(6,  m, 500))

for (cofigure in cofigures){
	k  = cofigure[1]
	m  = cofigure[2]
	n  = cofigure[3]
	
	# get appropriate W
	if (k == 6){
		W = hannum$W[,c("Gran","CD4T","CD8T","Mono","B","NK")]
	}else if(k == 5){
		W = cbind(hannum$W[,c("Gran","Mono","B","NK")], rowSums(hannum$W[,c("CD4T","CD8T")]))
		colnames(W) = c("Gran","Mono","B","NK", "T")
	}else if(k == 4){
		W = cbind(hannum$W[,c("Gran","Mono","B")], rowSums(hannum$W[,c("CD4T","CD8T", "NK")]))
		colnames(W) = c("Gran","Mono","B", "TNK")
	}else if(k == 3){
		W = cbind(hannum$W[,c("Gran","Mono")], rowSums(hannum$W[,c("CD4T","CD8T", "NK", "B")]))
		colnames(W) = c("Gran","Myeloid", "Lymphoid")
	}else{
		print("error")
	}
	assert(all(abs(rowSums(W) - 1) < 0.001))
	
	
	for(t in 1:permute.size){
		set.seed(t)
		config = paste(paste0("k_",k),
									 paste0("m_",m),
									 paste0("n_",n),
									 paste0("t_",t), sep = ".")
		print(paste0("generating simulated data: ", config))
	

		
		sim.X  = X[sample(nrow(X), replace = F, m), 
					     sample(ncol(X), replace = F, n)]
		
		sim.sample.ids  = colnames(sim.X)
		sim.feature.ids = rownames(sim.X)

		
		sim.W  = W[sim.sample.ids, ]
		sim.C1 = C1[sim.sample.ids, ] 
		sim.C2 = C2[sim.sample.ids, ] 
		
		sim.data = list(X = sim.X,
										W = sim.W,
										C1 = sim.C1,
										C2 = sim.C2)
		
		saveRDS(sim.data, file.path(rds.dir, paste0("sim.", config, ".rds")))
		fwrite(as.data.frame(sim.W),   
					 file = file.path(txt.dir, "res", paste0("W.", config, ".txt")),
					 sep = "\t", quote=FALSE, row.names = T, col.names = T)
		
		fwrite(as.data.frame(sim.X),   
					 file = file.path(txt.dir, "src", paste0("X.", config, ".txt")), 
					 sep = "\t", quote=FALSE, row.names = T, col.names = T)
	
	}
	
}
