library(ggplot2)
library(ggExtra)
library(extrafont)
library(ggpubr)
source("analysis.utils.r")
set.seed(2023)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

PBMC.color = brewer.pal(n = 12, name = "Paired")[2]
HLCA.color = brewer.pal(n = 12, name = "Paired")[8]


hat.color = brewer.pal(n = 12, name = "Paired")[1]
true.color = brewer.pal(n = 12, name = "Paired")[6]

source.1.color = brewer.pal(n = 8, name = "Dark2")[1]
source.2.color = brewer.pal(n = 8, name = "Dark2")[2]
source.3.color = brewer.pal(n = 8, name = "Dark2")[3]

color.df = data.frame(x = c("source.1", "source.2", "source.3"), 
                      y = c(1,2,3))
colar.bar.g <- ggplot(color.df, aes(x=y, fill = x, color=x)) + 
               geom_bar(colour = "black") + 
               scale_fill_manual( values = c(source.1.color, source.2.color, source.3.color)) + 
               theme_classic() 
             

colar.bar.g

PBMC.pb = readRDS("../Data//RNA//Simulation-PBMC/pbmc.pseudobulk.decon.L1.all.W.rds")
HLCA.pb = readRDS("../Data//RNA//Simulation-Lung//HLCA.pseudobulk.decon.L1.all.W.rds")

PBMC.entropy = PBMC.pb$entropy[PBMC.pb$HEF, ]
PBMC.entropy = PBMC.entropy/max(PBMC.entropy)

HLCA.entropy = HLCA.pb$entropy[HLCA.pb$HEF, ]
HLCA.entropy = HLCA.entropy/max(HLCA.entropy)

# entropy.df = data.frame(entropy = c(PBMC.entropy, HLCA.entropy),
#                         dataset = c(rep("PBMC", length(PBMC.entropy)), rep("Lung", length(HLCA.entropy))))

# entropy.df$dataset = factor(entropy.df$dataset, levels = c("PBMC", "Lung"))

entropy.df = data.frame(entropy = c(PBMC.entropy),
                        dataset = c(rep("PBMC", length(PBMC.entropy))))

entropy.df$dataset = factor(entropy.df$dataset, levels = c("PBMC"))

entropy.df

lab.size = 20
entropy.p <- ggplot(entropy.df, aes(x=entropy, fill = dataset, color=dataset)) + 
             geom_histogram(position="identity", alpha=0.5, binwidth=0.02, colour='black') + 
             scale_fill_manual( values = c(PBMC.color, HLCA.color)) + 
             xlab("Normalized entropy") + 
             ylab("Number of genes") + 
              
             theme_classic() + 
             theme(text=element_text(size=lab.size)) +   
             coord_cartesian(xlim = c(0, 1), ylim = c(0, 1000)) + 
             theme(legend.title = element_blank())  + 
             theme(legend.position = "top") + 
             theme(legend.direction = "horizontal",
                   #legend.position = c(0.6, 0.99),
                   legend.position = "none",
                   legend.text = element_text(size = lab.size,
                                              margin = margin(r = 100, unit = "pt")))

options(repr.plot.width = 7, repr.plot.height = 7, repr.plot.res = 100)
entropy.p

lab.size = 18
title.size = 18

plot_density = function(df, title = "", title.size = 20,  lab.size = 20,
                        xlab.name = "Expression level (cell type 1)", 
                        ylab.name = "Expression level (cell type 2)",
                        size.dot = 1, size.scatter = 4, size.marg.line = 0.8,
                        alpha.desity = 0.15, bins.densigram = 30,
                        extend.lb.iqr = 0.75, extend.ub.iqr = 1.5, 
                        mark.x = NULL, mark.y = NULL, mark.color = "red", mark.size = 1){
    
    p = ggplot(df, aes(x = x, y = y, )) +
        geom_point(data=df[(df$x!=0) & (df$y!=0), ],
                   size=size.dot) + 
    
        stat_density2d(data=df[(df$x!=0) & (df$y!=0), ], geom="polygon", size = 10,  alpha=alpha.desity) + 

        xlim(c(quantile(df$x, 0.05) - extend.lb.iqr * IQR(df$x), 
                 quantile(df$x, 0.95) + extend.ub.iqr * IQR(df$x)))+
        ylim(c(quantile(df$y, 0.05) - extend.lb.iqr * IQR(df$y), 
                 quantile(df$y, 0.95) + extend.ub.iqr * IQR(df$y))) + 
        #ggtitle(title) + 
        xlab(xlab.name) + 
        ylab(ylab.name) + 
        theme_classic() + 
        theme(plot.title = element_text(hjust = 0.5, size = title.size),
              legend.position = "none",
              text=element_text(size=lab.size)) 


    if(!is.null(mark.x) & !is.null(mark.y)){
        mark.df = data.frame(x = mark.x, y = mark.y, color = mark.color)
        p = p + geom_point(data = mark.df, aes(x=x, y=y, color = factor(color)), size = mark.size) + 
                scale_color_manual(values = c(mark.color))
    }
    

    p = ggMarginal(p,type="densigram",  
                   size = size.scatter,
                   bins = bins.densigram,
                   xparams = list(fill = source.1.color),
                   yparams = list(fill = source.2.color, size=size.marg.line))
    return(p)
}

str(PBMC.pb)

options(repr.plot.width = 8, repr.plot.height =6, repr.plot.res = 100)
hist(PBMC.pb$entropy[PBMC.pb$HEF[1:1000],])

# celltypes are not sorted in the PBMC.pb$W
sort(colMeans(PBMC.pb$W))

PBMC.entropy = as.matrix(PBMC.pb$entropy[PBMC.pb$HEF[1:1000],])

plot.low.gs  = sample(rownames(PBMC.entropy[PBMC.entropy < quantile(PBMC.entropy, 0.25),,drop = F]), 20)
plot.high.gs = sample(rownames(PBMC.entropy[PBMC.entropy > quantile(PBMC.entropy, 0.95),,drop = F]), 20)

low.plts = list()
for (gene in plot.low.gs){
    df = as.data.frame(t(PBMC.pb$Z[c("CD4","Mono"), gene ,]))
    colnames(df) = c("x", "y")
    low.plts[[gene]] = plot_density(df, title = gene, 
                                    xlab.name = "Expression (CD4)", ylab.name = "Expression (Mono)")
}

high.plts = list()
for (gene in plot.high.gs){
    df = as.data.frame(t(PBMC.pb$Z[c("CD4","Mono"), gene ,]))
    colnames(df) = c("x", "y")
    high.plts[[gene]] = tryCatch({
        plot_density(df, title = gene, 
                     xlab.name = "Expression (CD4)", ylab.name = "Expression (Mono)")
    }, error = function(e) {
        print("ERROR")
    })                  
}

options(repr.plot.width = 24, repr.plot.height =20, repr.plot.res = 100)
ggarrange(plotlist = low.plts)

options(repr.plot.width = 24, repr.plot.height =20, repr.plot.res = 100)
ggarrange(plotlist = high.plts)

round(PBMC.pb$params$corrs["CHCHD2",
                           c("CD4", "Mono", "B", "CD8", "NK"), 
                           c("CD4", "Mono", "B", "CD8", "NK")], 2)

PBMC.entropy["CHCHD2"]

round(PBMC.pb$params$corrs["SLC2A3",
                           c("CD4", "Mono", "B", "CD8", "NK"), 
                           c("CD4", "Mono", "B", "CD8", "NK")], 2)

PBMC.entropy["SLC2A3"]

plts = list()
for (gene in c("CHCHD2", "SLC2A3",   #other candidate listed below
               "NME2", "NDUFB11", "TRIR", "PRR13", "COX6C", "CDC42", "CDC42", "ATP5F1D")){
    df = as.data.frame(t(PBMC.pb$Z[c("CD4","Mono"), gene ,]))
    #df = as.data.frame(t(PBMC.pb$Z[c("CD4","CD8"), gene ,]))
    colnames(df) = c("x", "y")
    plts[[gene]] = tryCatch({
        plot_density(df, title = gene,
                     xlab.name = "Expression (CD4)", ylab.name = "Expression (Mono)",
                     title.size = title.size,  lab.size = lab.size)
                     
    }, error = function(e) {
        print("ERROR")
    })                  
}

PBMC.pb$entropy[c("CHCHD2", "SLC2A3"),]

options(repr.plot.width = 4, repr.plot.height = 4, repr.plot.res = 150)
plts[["CHCHD2"]]

plts[["SLC2A3"]]

#countour with know distribution 
get_contour_df = function(mu, sigma){
    data.grid <- expand.grid(s.1 = seq(mu[1] - 4*sqrt(sigma[1,1]), mu[1] + 4*sqrt(sigma[1,1]), length.out=200), 
                             s.2 = seq(mu[2] - 4*sqrt(sigma[2,2]), mu[2] + 4*sqrt(sigma[2,2]), length.out=200))

    q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu, sigma = sigma))
    return(q.samp)
}

cond_contour = function(plot.df, mu, true, 
                        extend.lb.iqr = 0.5, extend.ub.iqr = 0.75, 
                        xlab.name  = "Deconvolution (cell type 1)", 
                        ylab.name = "Deconvolution (cell type 2)",
                        title =  expression(paste(E, "[", Z, "|", X, "]")),
                        hat.color = "black", true.color = "red", alpha.desity = 0.15, 
                        title.size = 20, lab.size = 20, point.size = 1){
    
    mark.df = data.frame(x = c(mu[1], true[1]),
                         y = c(mu[2], true[2]),
                         color = c("Estimate", "True"))
    mark.df$color = factor(mark.df$color, level = c("Estimate", "True"))

    g = ggplot(plot.df, aes(x = x, y = y, )) +
        stat_density2d(data=plot.df, geom="polygon", size = 10,  alpha=alpha.desity)+
        geom_point(data = mark.df, aes(x=x, y=y, color = factor(color)), size = point.size)+
        #scale_shape_manual(values = c('Estimate' = 16, 'True' = 8)) + 
        scale_color_manual(values = c(hat.color, true.color))+
    
        xlim(c(quantile(plot.df$x, 0.05) - extend.lb.iqr * IQR(plot.df$x), 
                 quantile(plot.df$x, 0.95) + extend.ub.iqr * IQR(plot.df$x)))+
        ylim(c(quantile(plot.df$y, 0.05) - extend.lb.iqr * IQR(plot.df$y), 
                 quantile(plot.df$y, 0.95) + extend.ub.iqr * IQR(plot.df$y))) + 
    
        #ggtitle(title) +
        xlab(xlab.name) + 
        ylab(ylab.name) + 
        theme_classic() + 
        theme(plot.title = element_text(hjust = 0.5, size = title.size),
              legend.key = element_rect(fill = "white", linewidth =0, colour = "white", color ="white", inherit.blank = T),
           
              #legend.background = element_rect(fill = "darkgray"), 
              legend.position = "none",
              legend.title = element_blank(),
              text=element_text(size=lab.size)) 
    return(g)
}

genes = readRDS("../Figure/Illustrative/toy.genes.rds")

str(genes[["CHCHD2"]])

mark.point.size = 4
visual.sample.id = 2
visual.source.ids = c("CD4", "Mono")

feature.id = "CHCHD2"
mu    = genes[[feature.id]]$mu_cond[visual.sample.id, visual.source.ids]
sigma = matrix(genes[[feature.id]]$sigmas_cond[visual.sample.id, visual.source.ids,visual.source.ids], nrow=2)
low.cond.df = as.data.frame(mvtnorm::rmvnorm(n = 10000, mean = mu, sigma = sigma))
colnames(low.cond.df) = c("x", "y")

true = c(genes$Z[,feature.id,visual.sample.id][visual.source.ids[1]], 
         genes$Z[,feature.id,visual.sample.id][visual.source.ids[2]])

low.cond.g = cond_contour(plot.df = low.cond.df, 
                           mu = mu, true = true,
                           extend.lb.iqr = 0.75, extend.ub.iqr = 0.75, 
                           xlab.name  = "CD4 deconvolution", 
                           ylab.name  = "Monocyte deconvolution",
                           hat.color = hat.color, true.color = true.color,
                           title.size = title.size, lab.size = lab.size, point.size = mark.point.size)




df = as.data.frame(t(PBMC.pb$Z[c("CD4","Mono"), feature.id,]))
colnames(df) = c("x", "y")
low.Z.g = plot_density(df, 
                        xlab.name = "CD4 expression", ylab.name = "Monocyte expression",
                        title.size = title.size,  lab.size = lab.size,
                        extend.lb.iqr = 0.5, extend.ub.iqr = 0.75, 
                        mark.x = true[1], mark.y = true[2], mark.color = true.color, mark.size = mark.point.size)

options(repr.plot.width = 4, repr.plot.height = 4, repr.plot.res = 150)
low.Z.g

options(repr.plot.width = 3.5, repr.plot.height = 3.25, repr.plot.res = 150)
low.cond.g

feature.id = "SLC2A3"
mu    = genes[[feature.id]]$mu_cond[visual.sample.id, visual.source.ids]
sigma = matrix(genes[[feature.id]]$sigmas_cond[visual.sample.id, visual.source.ids,visual.source.ids], nrow=2)
high.cond.df = as.data.frame(mvtnorm::rmvnorm(n = 10000, mean = mu, sigma = sigma))
colnames(high.cond.df) = c("x", "y")

true = c(genes$Z[,feature.id,visual.sample.id][visual.source.ids[1]], 
         genes$Z[,feature.id,visual.sample.id][visual.source.ids[2]])

high.cond.g = cond_contour(plot.df = high.cond.df, 
                           mu = mu, true = true,
                           extend.lb.iqr = 0.75, extend.ub.iqr = 0.75, 
                           xlab.name  = "CD4 deconvolution", 
                           ylab.name  = "Monocyte deconvolution",
                           hat.color = hat.color, true.color = true.color,
                           title.size = title.size, lab.size = lab.size, point.size = mark.point.size)




df = as.data.frame(t(PBMC.pb$Z[c("CD4","Mono"), feature.id,]))
colnames(df) = c("x", "y")
high.Z.g = plot_density(df, 
                        #extend.lb.iqr = 0.5, extend.ub.iqr = 0.75, 
                        xlab.name = "CD4 expression", ylab.name = "Monocyte expression",
                        title.size = title.size,  lab.size = lab.size,
                        mark.x = true[1], mark.y = true[2], mark.color = "red", mark.size = mark.point.size)

options(repr.plot.width = 4, repr.plot.height = 4, repr.plot.res = 150)
high.Z.g

options(repr.plot.width = 3.5, repr.plot.height = 3.25, repr.plot.res = 150)
high.cond.g

mark.df = data.frame(x = c(1, 3),
                     y = c(2, 2),
                     color = c("Estimate", "True"))
mark.df$color = factor(mark.df$color, level = c("Estimate", "True"))

legend.g = ggplot(mark.df, aes(x=x, y=y, color=color)) +
           geom_point(size = 7)+
           scale_color_manual(values = c(hat.color, true.color))+
           theme_classic() + 
           theme(legend.title = element_blank(),
                 legend.position = "top", 
                 text=element_text(size=lab.size))  



legend.g = as_ggplot(get_legend(legend.g + 
                                 theme(legend.direction = "horizontal",
                                      legend.position = c(0.575, 0.5),
                                      legend.text = element_text(size=lab.size, 
                                                                 margin = margin(r = 30, unit = "pt"))))) 

options(repr.plot.width = 4, repr.plot.height = 0.5, repr.plot.res = 100)
legend.g

#options(repr.plot.width = 16, repr.plot.height = 5, repr.plot.res = 100)
top.g = ggarrange(NULL, NULL, nrow = 1, 
                  widths = c(1.5, 1.25),
                  labels = c("a", "b"),
                  font.label = list(size = 30, color = "black", face = "bold", family = NULL))

bottom.g = ggarrange(NULL, NULL, nrow = 1, 
                     widths = c(1, 1),
                     labels = c("c", "d"),
                     font.label = list(size = 30, color = "black", face = "bold", family = NULL))

layout.g = ggarrange(top.g, bottom.g, nrow = 2, 
                     heights = c(1.1, 1), widths  = c(1,1))

options(repr.plot.width = 16, repr.plot.height = 12, repr.plot.res = 100)
layout.g

ggsave("../Figure/Illustrative/colar.bar.g.pdf", colar.bar.g,
       bg = 'white',
       device = "pdf", width = 5, height = 5, dpi = 600)

ggsave("../Figure/Illustrative/layout.g.pdf", layout.g,
       bg = 'white',
       device = "pdf", width = 16, height = 12, dpi = 600)

ggsave("../Figure/Illustrative/entropy.g.pdf", entropy.p,
       bg = 'white',
       device = "pdf", width = 8, height = 6.6, dpi = 600) # later shrink in post to height 6

ggsave("../Figure/Illustrative/low.Z.g.pdf", low.Z.g,
       bg = 'white',
       device = "pdf", width = 4, height = 4, dpi = 600)

ggsave("../Figure/Illustrative/low.cond.g.pdf", low.cond.g,
       bg = 'white',
       device = "pdf", width = 3.5, height = 3.25, dpi = 600)

ggsave("../Figure/Illustrative/high.Z.g.pdf", high.Z.g,
       bg = 'white',
       device = "pdf", width = 4, height = 4, dpi = 600)

ggsave("../Figure/Illustrative/high.cond.g.pdf", high.cond.g,
       bg = 'white',
       device = "pdf", width = 3.5, height = 3.25, dpi = 600)

ggsave("../Figure/Illustrative/cond.legend.g.pdf", legend.g,
       bg = 'white',
       device = "pdf", width = 4, height = 0.5, dpi = 600)

fig.list = list(entropy.g = entropy.p,
                low.corrs = PBMC.pb$params$corrs["CHCHD2",,],
                low.Z.g = low.Z.g, 
                low.cond.g = low.cond.g,
                high.corrs = PBMC.pb$params$corrs["SLC2A3",,],
                high.Z.g = high.Z.g, 
                high.cond.g = high.cond.g,
                cond.legend.g = legend.g)

saveRDS(fig.list, "../Figure/Illustrative/illustrative.fig.list.rds")








