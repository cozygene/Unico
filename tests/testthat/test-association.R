library(Unico)

test_that("association testing on cell-type and tissue level covariates has good power and FP control", {

	skip_on_cran()

	#source-specific association

	download.file("https://github.com/cozygene/Unico/raw/main/tests/assets/simulation.gammas.10.rds", "simulation.gammas.10.rds")
	sim.data = readRDS("simulation.gammas.10.rds")
	file.remove("simulation.gammas.10.rds")

	Unico.mdl = list()
	Unico.mdl$params.hat <- Unico(sim.data$X, sim.data$W, C1 = sim.data$C1, C2 = sim.data$C2, parallel = F, log_file=NULL)
	Unico.mdl$params.hat = association_parametric(X = sim.data$X, Unico.mdl$params.hat, parallel = F, log_file=NULL)
	#marginal power
	marg.pvals = Unico.mdl$params.hat$parametric$gammas_hat_pvals
	expect_equal(sum(marg.pvals < 0.05/(nrow(marg.pvals) * ncol(marg.pvals)))/(nrow(marg.pvals) * ncol(marg.pvals)) > 0.9, TRUE)
	#joint power
	joint.pvals = Unico.mdl$params.hat$parametric$gammas_hat_pvals.joint
	expect_equal(sum(joint.pvals < 0.05/(length(joint.pvals)))/(length(joint.pvals)) > 0.95, TRUE)
	#FP control on non-source-specific
	global.marg.pvals = Unico.mdl$params.hat$parametric$betas_hat_pvals
	expect_equal(sum(global.marg.pvals < 0.05/(nrow(global.marg.pvals) * ncol(global.marg.pvals)))/(nrow(global.marg.pvals) * ncol(global.marg.pvals)) < 0.05, TRUE)

	#non-source-specific association
	download.file("https://github.com/cozygene/Unico/raw/main/tests/assets/simulation.betas.10.rds","simulation.betas.10.rds")
	sim.data = readRDS("simulation.betas.10.rds")
	file.remove("simulation.betas.10.rds")

	Unico.mdl = list()
	Unico.mdl$params.hat <- Unico(sim.data$X, sim.data$W, C1 = sim.data$C1, C2 = sim.data$C2, parallel = F, log_file=NULL)
	Unico.mdl$params.hat = association_parametric(X = sim.data$X, Unico.mdl$params.hat, parallel = F, log_file=NULL)
	#marginal FP control
	marg.pvals = Unico.mdl$params.hat$parametric$gammas_hat_pvals
	expect_equal(sum(marg.pvals < 0.05/(nrow(marg.pvals) * ncol(marg.pvals)))/(nrow(marg.pvals) * ncol(marg.pvals)) < 0.05, TRUE)
	#joint FP control
	joint.pvals = Unico.mdl$params.hat$parametric$gammas_hat_pvals.joint
	expect_equal(sum(joint.pvals < 0.05/(length(joint.pvals)))/(length(joint.pvals)) < 0.05, TRUE)
	#non-source-specific power
	global.marg.pvals = Unico.mdl$params.hat$parametric$betas_hat_pvals
	expect_equal(sum(global.marg.pvals < 0.05/(nrow(global.marg.pvals) * ncol(global.marg.pvals)))/(nrow(global.marg.pvals) * ncol(global.marg.pvals)) > 0.95, TRUE)

})
