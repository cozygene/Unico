library(Unico)

test_that("parallel works on model fitting, tensor and association", {
	skip_on_cran()
	expect_equal(detectCores() > 0, TRUE)
	sim.data = simulate_data(n = 100, m = 5, k = 3, p1 = 1, p2 = 2, taus_std = 0)

	pass <- FALSE
	result = tryCatch({
		unico.res = list()
		unico.res$params.hat <- Unico(sim.data$X, sim.data$W, C1 = sim.data$C1, C2 = sim.data$C2, parallel = T, num_cores = 2, log_file=NULL)
		unico.res$Z.hat <- tensor(sim.data$X, sim.data$W, C1 = sim.data$C1, C2 = sim.data$C2, Unico.mdl = unico.res$params.hat, parallel = T, num_cores = 2, log_file=NULL)
		unico.res$params.hat <- association_parametric(sim.data$X, unico.res$params.hat, parallel = T, num_cores = 2, log_file=NULL)
		unico.res$params.hat <- association_asymptotic(sim.data$X, unico.res$params.hat, parallel = T, num_cores = 2, log_file=NULL)
		pass <- TRUE
	})
	expect_equal(pass, TRUE)
})

