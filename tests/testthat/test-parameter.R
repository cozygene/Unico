library(Unico)

test_that("estimated params and tensors are consisitent with frozen result", {
	skip_on_cran()


	#basedir <- "../assets"
	download.file("https://github.com/cozygene/Unico/raw/main/tests/assets/simulation.rds","simulation.rds")
	download.file("https://github.com/cozygene/Unico/raw/main/tests/assets/frozen.mdl.rds","frozen.mdl.rds")

	sim.data = readRDS("simulation.rds")
	frozen.mdl = readRDS("frozen.mdl.rds")

	file.remove("simulation.rds")
	file.remove("frozen.mdl.rds")


	#model fit
	set.seed(2023)
	Unico.mdl = list()
	Unico.mdl$params.hat <- Unico(sim.data$X, sim.data$W, C1 = sim.data$C1, C2 = sim.data$C2, parallel = F, log_file=NULL)
	n = nrow(sim.data$W)
	m = nrow(sim.data$X)
	k = ncol(sim.data$W)
	p1 = ncol(sim.data$C1)
	p2 = ncol(sim.data$C2)

	# mean check
	for (h in 1:k){
		expect_equal(cor(frozen.mdl$mus_hat[,h], Unico.mdl$params.hat$mus_hat[,h]) > 0.99, TRUE)
	}

	for (p in 1:(p1*k)){
		expect_equal(cor(frozen.mdl$gammas_hat[,p], Unico.mdl$params.hat$gammas_hat[,p]) > 0.99, TRUE)
	}

	for (p in 1:p2){
		expect_equal(cor(frozen.mdl$betas_hat[,p], Unico.mdl$params.hat$betas_hat[,p]) > 0.99, TRUE)
	}

	# variance check
	for (h in 1:k){
		if(colMeans(sim.data$W)[h] > 0.2){
			thr = 0.8
		}else if (colMeans(sim.data$W)[h] > 0.1){
			thr = 0.5
		}else{
			thr = -1
		}
		expect_equal(cor(frozen.mdl$sigmas_hat[,h,h], Unico.mdl$params.hat$sigmas_hat[,h,h]) > thr, TRUE)
	}

	# tensor check
	Unico.mdl$params.hat$mus_hat = frozen.mdl$mus_hat
	Unico.mdl$params.hat$gammas_hat = frozen.mdl$gammas_hat
	Unico.mdl$params.hat$betas_hat = frozen.mdl$betas_hat
	Unico.mdl$params.hat$sigmas_hat = frozen.mdl$sigmas_hat
	Unico.mdl$Z.hat <- tensor(sim.data$X, sim.data$W, C1 = sim.data$C1, C2 = sim.data$C2, Unico.mdl = Unico.mdl$params.hat, parallel = F, log_file=NULL)
	thr = 10**(-4)
	for (h in 1:k){
		expect_equal(mean((Unico.mdl$Z.hat[h,,] - frozen.mdl$Z_hat[h,,])**2) < thr, TRUE)
	}

})


