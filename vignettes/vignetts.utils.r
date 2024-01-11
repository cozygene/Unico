library(MASS)

# x: list of numbers
# y: list of numbers
# robust: boolean. indicate if using robust correlation
# qtl: numeric between (0,1), only meaningful when robust is turned on. indicate the fraction of the data that is considered to be "good portion" and will participate in the correlation calculation
# check the correlation between x and y after excluding outliers defined by the qtl
# note that if in IQR is 0 in either x or y, OR x and y are completely colinear, robust correlation will be NA and thus filled by 0
safe_cor = function(x, y, robust = TRUE, qtl = 0.95){
	res = tryCatch({
		if(robust){
			cov.rob(cbind(x,y), cor = TRUE, quantile.used = round(qtl*length(x)) ,method = "mve")$cor[1,2]
		}else{
			cor(x,y)
		}

	}, error = function(cond){
		0
	})

	# in regular correlation if all constant, then it is going to be NA
	if (is.na(res)){
		return(0)
	}else{
		return(res)
	}

}


# Z.true is an array of k by m by n, the true tensor
# Z.hat is an array of k by m by n, the estimated tensor
# eval.feature.source is a matrix of logical values that is number of features by number of sources. only calculate correlation on feature-source that has TRUE in this matrix
# robust: boolean. indicate if using robust correlation
# qtl: numeric between (0,1), only meaningful when robust is turned on. indicate the fraction of the data that is considered to be "good portion" and will participate in the correlation calculation
# return a m by k correltion matrix, each column is a different source and each row is a different feature,
# each entry is a correlation score, indicating that feautre, that source's tensor estimated across samples
# those entry with FALSE in the eval.feature.source is going to be set to NA
calc_Z_corrs <- function(Z.true, Z.hat, eval.feature.source = NULL, robust = TRUE, qtl = 0.95){

	m = dim(Z.true)[2]
	k = dim(Z.true)[1]
	Z.corrs = matrix(NA, m, k)
	rownames(Z.corrs) = dimnames(Z.true)[[2]]
	colnames(Z.corrs) = dimnames(Z.true)[[1]]

	# if eval.feature.source is not present, return all results
	if(is.null(eval.feature.source)){
		eval.feature.source = matrix(TRUE, m, k)
	}

	for (h in 1:k){
		for (j in 1:m){
			if(eval.feature.source[j,h]){
				Z.corrs[j,h] <- safe_cor(Z.true[h,j,], Z.hat[h,j,], robust, qtl)
			}else{
				Z.corrs[j,h] <- NA
			}

		}
	}
	return(Z.corrs)
}


library("ggplot2")
library("ggpubr")
library("hexbin")

plot_qq_compare <- function(pvals_mat1, pvals_mat2, labels, ggarrange.nrow = 1, ggarrange.ncol = 1,
														alpha = 0.05, text.size = 10,
														xlab = "pval1", ylab = "pval2", title = ""){

	pvals1 = list()
	pvals2 = list()

	for (l in 1:ncol(pvals_mat1)){
		pvals1[[l]] = pvals_mat1[,l]
		pvals2[[l]] = pvals_mat2[,l]
	}



	qqplots <- lapply(1:length(pvals1), function(p){
		df <- data.frame(pvals.1 = -log10(pvals1[[p]]),
										 pvals.2 = -log10(pvals2[[p]]));
		qqplot <- ggplot(df, aes(x = pvals.1, y = pvals.2)) +
			stat_binhex(geom = "point", bins=1000, size=1, alpha = alpha) +
			geom_abline() +
			coord_cartesian(xlim = c(min(df), max(df)),
											ylim = c(min(df), max(df))) +

			theme_bw() +
			guides(fill="none") +
			ggtitle(paste0(labels[p],
										 " (r = ",
										 round(safe_cor(df$pvals.1,  df$pvals.2, robust = T, qtl = 0.95),2),
										 ")")) +
			xlab(parse(text = paste0("\'",xlab,"\'~", expression(-log[10](P))))) +
			ylab(parse(text = paste0("\'",ylab,"\'~", expression(-log[10](P))))) +


			theme(plot.title = element_text(hjust = 0.5, size = text.size, face="bold"))+
			theme(axis.title.x = element_text(size = text.size)) +
			theme(axis.text.x = element_text(size = text.size * 0.75)) +
			theme(axis.title.y = element_text(size = text.size)) +
			theme(axis.text.y = element_text(size = text.size * 0.75))
		return(qqplot)
	})
	g = egg::ggarrange(plots = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
	g = annotate_figure(g,top = text_grob(title,
																				color = "black", face = "bold", size = text.size + 5))
	return(g)
}

plot_qq<- function(pvals_mat, labels, ggarrange.nrow = 1, ggarrange.ncol = 1,
									 alpha = 0.5, text.size = 10,
									 title = ""){

	pvals = list()

	for (l in 1:ncol(pvals_mat)){
		pvals[[l]] = pvals_mat[,l]
	}

	significance_th <- list(alpha/length(pvals[[1]]))
	if(length(pvals)-1) significance_th[[2]] <- alpha/(length(pvals)*length(pvals[[1]]))

	qqplots <- lapply(1:length(pvals), function(p){
		df <- data.frame(pvals.obs = -log10(sort(pvals[[p]])),
										 pvals.exp = -log10(sort((1:length(pvals[[p]]))/length(pvals[[p]]))));
		qqplot <- ggplot(df, aes(x = pvals.exp, y = pvals.obs)) +
			stat_binhex(geom = "point", bins=1000, size=1, alpha = alpha) +
			geom_abline() +
			coord_cartesian(#xlim = c(min(df), max(-log10(df$pvals.exp))),
				ylim = c(min(df), max(max(df), -log10(significance_th[[2]])))
				#ylim = c(min(df),  -log10(significance_th[[2]]))
			) +

			theme_bw() +
			guides(fill="none") +
			ggtitle(labels[p]) +
			xlab(expression(Expected~-log[10](P))) +
			ylab(expression(Observed~-log[10](P))) +


			theme(plot.title = element_text(hjust = 0.5, size = text.size, face="bold"))+
			theme(axis.title.x = element_text(size = text.size)) +
			theme(axis.text.x = element_text(size = text.size * 0.75)) +
			theme(axis.title.y = element_text(size = text.size)) +
			theme(axis.text.y = element_text(size = text.size * 0.75))

		#qqplot = qqplot + geom_hline(yintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1  ) +
		#                  geom_hline(yintercept=-log10(significance_th[[2]]), linetype="dashed", color = "red", size=0.5)
		# only show the celltpye cpg wide
		qqplot = qqplot + geom_hline(yintercept=-log10(significance_th[[2]]), linetype="dashed", color = "red", size=1)
		return(qqplot)
	})
	g = egg::ggarrange(plots = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
	g = annotate_figure(g,top = text_grob(title,
																				color = "black", face = "bold", size = text.size + 5))
	return(g)
}
