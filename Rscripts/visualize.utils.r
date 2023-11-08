################################################################################################
########unchecked
########################################################unchecked

library("ggplot2")
library("ggpubr")
library("grid")
library("hrbrthemes")
library("extrafont")
library("hash")
library("rlang")
library("scales")

plot_qq <- function(pvals, labels, ggarrange.nrow = 1, ggarrange.ncol = 1, alpha = 0.05, experiment_wide_line = TRUE){
  significance_th <- list(alpha/length(pvals[[1]]))
  if(length(pvals)-1) significance_th[[2]] <- alpha/(length(pvals)*length(pvals[[1]]))
  qqplots <- lapply(1:length(pvals), function(p){
    df <- data.frame(pvals.obs = -log10(sort(pvals[[p]])), pvals.exp = -log10(sort((1:length(pvals[[1]]))/length(pvals[[1]]))));
    qqplot <- ggplot(df, aes(x = pvals.exp, y = pvals.obs)) +
      stat_binhex(geom = "point", bins=1000, size=1) +
      geom_abline() +
      ggtitle(labels[p]) +
      xlab(expression(Expected~-log[10](P))) + ylab(expression(Observed~-log[10](P))) +
      theme_bw() +
      guides(fill="none") +
      geom_hline(yintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)
    if (length(significance_th)-1 & experiment_wide_line) qqplot <- qqplot + geom_hline(yintercept=-log10(significance_th[[2]]), linetype="dashed", color = "red", size=0.5)
    return(qqplot)
  })
  ggarrange(plotlist = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
}

plot_qq_prime = function(pvals_mat, experiment, cts){
  pvals_list = list()
  for (l in 1:length(cts)){
    pvals_list[[l]] = pvals_mat[,l]
  }
  plot_qq(pvals_list,
          labels = paste(experiment, cts),
          ggarrange.nrow = 1,
          ggarrange.ncol = length(cts),
          experiment_wide_line = TRUE)
}

plot_qq_compare <- function(pvals1, pvals2, labels, ggarrange.nrow = 1, ggarrange.ncol = 1, alpha = 0.05, experiment_wide_line = TRUE){
  significance_th <- list(alpha/length(pvals1[[1]]))
  if(length(pvals1)-1) significance_th[[2]] <- alpha/(length(pvals1)*length(pvals1[[1]]))
  qqplots <- lapply(1:length(pvals1), function(p){
    df <- data.frame(pvals.1 = -log10(pvals1[[p]]), 
                     pvals.2= -log10(pvals2[[p]]));
    qqplot <- ggplot(df, aes(x = pvals.1, y = pvals.2)) +
      stat_binhex(geom = "point", bins=1000, size=1) +
      geom_abline() +
      ggtitle(labels[p]) +
      xlab(expression(pval1~-log[10](P))) + ylab(expression(pval2~-log[10](P))) +
      theme_bw() +
      guides(fill="none") +
      geom_hline(yintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)+
      geom_vline(xintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)
    if (length(significance_th)-1 & experiment_wide_line) qqplot <- qqplot + geom_hline(yintercept=-log10(significance_th[[2]]), linetype="dashed", color = "red", size=0.5)
    return(qqplot)
  })
  ggarrange(plotlist = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
}

plot_qq_compare_prime = function(pvals_mat1, pvals_mat2, experiment, cts){
  pvals_list1 = list()
  pvals_list2 = list()
  for (l in 1:length(cts)){
    pvals_list1[[l]] = pvals_mat1[,l]
    pvals_list2[[l]] = pvals_mat2[,l]
  }
  plot_qq_compare(pvals_list1,
                  pvals_list2,
                  labels = paste(experiment, cts),
                  ggarrange.nrow = 1,
                  ggarrange.ncol = length(cts),
                  experiment_wide_line = FALSE)
}

get_recon_X = function(Z, W, C2 = NULL, betas = NULL, scale.factor = NULL){
    #C2, and betas need to be tall matrix
    n = dim(Z)[3]
    m = dim(Z)[2]
    k = dim(Z)[1]
    
    feature.ids =  dimnames(Z)[[2]]
    sample.ids =  dimnames(Z)[[3]]
    
    if(is.null(C2)){
        C2 = matrix(0, n,0)
    }
    if (is.null(betas)){
        betas = matrix(0, 0, 1)
    }
    
    X = matrix(0,m,n)
    rownames(X) = dimnames(Z)[[2]]
    colnames(X) = dimnames(Z)[[3]]
    
    for (l in 1:k){
        X = X + Z[l,,] * repmat(t(as.matrix(W[,l])), m, 1)
    }
    #adding C2 effect
    X = X + t(C2 %*% t(betas)) * repmat(scale.factor[feature.ids, , drop = F], 1, length(sample.ids))
    return(X)
}

