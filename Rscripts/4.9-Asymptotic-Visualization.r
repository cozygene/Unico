source("visualize.utils.r")
source("analysis.utils.r")
set.seed(2023)

pheno = "age"
#pheno = "gender"

figure.dir  = "../Figure/Asymptotics"
result.dir = "/u/project/halperin/johnsonc/Unico/Unico2023/Result/Methylation/Consistency/XY/"
if (!file.exists(figure.dir)){dir.create(file.path(figure.dir),recursive = T)}
if (!file.exists(result.dir)){print("check results_dir")}

data.names = c("liu", "hannon1", "hannon2", "hannum")
full.data.names = list("liu" = "Liu et al. (n= 687)",
                       "hannon1" = "Hannon et al. I (n= 675)", 
                       "hannon2" = "Hannon et al. II (n= 665)",
                       "hannum" = "Hannum et al. (n= 590)")

#make sure every dataset pvals columns are in the same order 
source.ids = c("Gran", "CD4T", "CD8T", "Mono", "B", "NK")
k = length(source.ids)

pval.list  = list()
for (data.name in data.names){
    pval  = list()
    
    # signal
    Unico.mdl = readRDS(file.path(result.dir, data.name, paste0("Unico.mdl.rds")))
    pval[[paste0("Unico.parametric")]] = Unico.mdl$params.hat[["parametric"]]$gammas_hat_pvals[,paste(source.ids, pheno, sep=".")]
    pval[[paste0("Unico.asymptotic")]] = Unico.mdl$params.hat[["asymptotic"]]$gammas_hat_pvals[,paste(source.ids, pheno, sep=".")]
    
    # null 
    Unico.mdl = readRDS(file.path(result.dir, paste0(data.name, "-", pheno, "-shuffled"), paste0("Unico.mdl.rds")))
    pval[[paste0("Unico.asymptotic.null")]] = Unico.mdl$params.hat[["asymptotic"]]$gammas_hat_pvals[,paste(source.ids, pheno, sep=".")]
    
    # collect result
    pval.list[[data.name]]  = pval
    
}

str(pval.list)

#source("visualize.utils.r")
text.size = 20
labels = paste0(if(pheno == "gender") "Sex: " else "Age: ", source.ids)
labels

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


options(repr.plot.width = 22, repr.plot.height = 4.5, repr.plot.res = 200)
plts = list()
for (t in 1:length(data.names)){
    plts[[t]] = plot_qq_compare(pvals_mat1 = pval.list[[t]][["Unico.parametric"]],
                                pvals_mat2 = pval.list[[t]][["Unico.asymptotic"]],
                               labels = labels, 
                               ggarrange.nrow = 1, ggarrange.ncol = k, 
                               alpha = 0.05, text.size = text.size,
                               xlab = "Parametric", ylab = "Asymptotic", 
                               title = full.data.names[[data.names[t]]])
}
compare = ggarrange(plts[[1]], NULL, plts[[2]], NULL, plts[[3]], NULL, plts[[4]], 
                    nrow = 2*length(data.names)-1, 
                    ncol = 1, 
                    heights = c(1, 0.01, 1, 0.05, 1, 0.05, 1)) 


options(repr.plot.width = 22, repr.plot.height = 18, repr.plot.res = 200)
compare



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


options(repr.plot.width = 22, repr.plot.height = 4.5, repr.plot.res = 200)
plts.null = list()
for (t in 1:length(data.names)){
    plts.null[[t]] = plot_qq(pvals_mat = pval.list[[t]][["Unico.asymptotic.null"]],
                                  labels = labels, 
                                  ggarrange.nrow = 1, ggarrange.ncol = k, 
                                  alpha = 0.5, text.size = text.size,
                                  title = full.data.names[[data.names[t]]])
}

calibrate = ggarrange(plts.null[[1]], NULL, plts.null[[2]], NULL, 
                       plts.null[[3]], NULL, plts.null[[4]], 
                       nrow = 2*length(data.names)-1, 
                       ncol = 1, 
                       heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1)) 

options(repr.plot.width = 22, repr.plot.height = 18, repr.plot.res = 200)
calibrate

file.path(figure.dir, "Asymptotic_fig.list.rds")

fig.list.file = file.path(figure.dir, "Asymptotic_fig.list.rds")
if(file.exists(fig.list.file)){
    fig.list = readRDS(fig.list.file)
}else{
    fig.list = list()
}


fig.list[[paste0(pheno,".plts.null")]]   = plts.null
fig.list[[paste0(pheno,".plts")]]        = plts

saveRDS(fig.list, fig.list.file)

ggsave(file.path(figure.dir, paste0("Parametric.vs.Asymptotic.", pheno, ".png")), 
       compare,
       device = "png", width = 22, height = 18, bg = 'white')

ggsave(file.path(figure.dir, paste0("Asymptotic.", pheno, ".null.png")), 
       calibrate,
       device = "png", width = 22, height = 18, bg = 'white')

