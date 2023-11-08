library(pracma)
library(data.table)


#X: matrix of features by samples, methylation beta values
#W: matrix of number of samples by number of celltypes. 
#C: matrix of number of samples by number of covariates, both global and celltype level covars. these will be the phenotypes as well
#this functions essentially fit Z.hat(XW) ~ 1 + W + C for marginal test
#this functions essentially fit X ~ 1 + W + C for joint test

#returns a list of 6 keys
#first two keys:  marg.pvals (k by m by number of covars), marg.coefsf (k by m by number of covars), 
#second two keys: joint.pvals(m by number of covars), joint.coefsf(m by number of covars). 

base_assos = function(X, W, C){

    source.ids  = colnames(W)
    feature.ids = rownames(X)
    sample.ids  = colnames(X)
    covar.ids   = colnames(C)


    # marginal related
    marg.pvals = array(1, c(length(source.ids), length(feature.ids), length(covar.ids)))
    dimnames(marg.pvals)[[1]] = source.ids
    dimnames(marg.pvals)[[2]] = feature.ids
    dimnames(marg.pvals)[[3]] = covar.ids
    marg.coefs = copy(marg.pvals)



    # joint related
    joint.pvals = matrix(1, length(feature.ids), length(covar.ids))
    rownames(joint.pvals) = feature.ids
    colnames(joint.pvals) = covar.ids
    joint.coefs = copy(joint.pvals)



    for (feature.id in feature.ids){    
        # marginal model
        for (source.id in source.ids){
            Z.hat = X[feature.id, sample.ids] * W[sample.ids, source.id]
            df   = data.frame(y = Z.hat, 
                              cbind(W[sample.ids,-1], C[sample.ids,])) # drop W to prevent co-linear

            mdl.fit = lm(y ~ 1 + ., data = df) 
            marg.pvals[source.id, feature.id, covar.ids] = summary(mdl.fit)$coef[covar.ids, "Pr(>|t|)"]
            marg.coefs[source.id, feature.id, covar.ids] = summary(mdl.fit)$coef[covar.ids, "Estimate"]   
        }

        df   = data.frame(y = as.vector(X[feature.id, sample.ids, drop = F]), 
                          cbind(W[sample.ids, -1], C[sample.ids,])) # drop W to prevent co-linear
        
        # joint model
        mdl.fit = lm(y ~ 1 + ., data = df) 
        joint.pvals[feature.id, covar.ids] = summary(mdl.fit)$coef[covar.ids, "Pr(>|t|)"]
        joint.coefs[feature.id, covar.ids] = summary(mdl.fit)$coef[covar.ids, "Estimate"] 


    }


    res = list(marg.pvals = marg.pvals, marg.coefs = marg.coefs, 
               joint.pvals = joint.pvals, joint.coefs = joint.coefs)
    
    return(res)
}