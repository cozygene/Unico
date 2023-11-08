#X: matrix of features by samples, methylation beta values
#y: vector of length number of samples, phenotype of interest to test, 
#W: matrix of number of samples by number of celltypes. 
#cov.mod: matrix of number of samples by number of "global" covariates, to be added to regression 
#this functions essentially fit x ~ 0 + W + yW + cov.mod

#returns a list of two keys
#first key: marg.res, a list of number of celltypes of matrix, each matrix is number of features by coef, SE, t, and p (from T test)
#second key: joint.res, a matrix of number of features by 1. joint p value from a F test
celldmc_w_joint = function(X, y, W, cov.mod = NULL){
    
    # check pheno
    if (nlevels(factor(y)) == 2) {
        message("Binary phenotype detected. Predicted change will be 1 - 0.")
        y <- factor(y)
    }
    if (!is.factor(y) & !is.character(y)) 
        message("y is not factor or character. Treating as continuous variables.")

    # design matrix
    design <- model.matrix(~ W + y:W)[, -1] # -1 removes the intercept term
    if (!is.null(cov.mod)) design <- cbind(design, cov.mod)
    IntNames.v <- str_c(colnames(W), "Pheno")
    colnames(design)[(1 + ncol(W)):(2*ncol(W))] <- IntNames.v 


    #testing 
    marg.res = array(1,c(nrow(X), ncol(W), 4))
    dimnames(marg.res)[[1]] = rownames(X)
    dimnames(marg.res)[[2]] = IntNames.v
    dimnames(marg.res)[[3]] = c("Estimate", "SE", "t", "p")

    joint.res = matrix(1, nrow(X), 1)
    rownames(joint.res) = rownames(X)
    colnames(joint.res) = c("p")
    
        
    for (feature.id in rownames(X)){
        x.df = as.data.frame(X[feature.id, ])
        colnames(x.df) = c("x")


        d1.df = cbind(x.df, design)
          
        # full model
        mdl1.fit = lm(x ~ 0 + ., data = d1.df) 
        marg.res[feature.id,,] = summary(mdl1.fit)$coe[IntNames.v, ]

        # restricted model
        d0.df = d1.df[, !(colnames(d1.df) %in% IntNames.v)]
        mdl0.fit = lm(x ~ 0 + ., data = d0.df) 
        anova.fit <- anova(mdl0.fit, mdl1.fit)
        joint.res[feature.id,] = anova.fit$"Pr(>F)"[2]
    }
    marg.res.list = list()
    for (h in 1:ncol(W)){
        source.id  = colnames(W)[h]
        marg.res.list[[source.id]] = marg.res[,h,]
    }
  
    return(list(marg.res = marg.res.list, joint.res = joint.res))
}