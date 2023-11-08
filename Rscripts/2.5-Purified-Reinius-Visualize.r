# save each individaul evaluation bar plots as ggplots
# save scatter plots as pdf

library("matrixStats")
library("data.table")
library("ggplot2")
library("ggpubr")
source("analysis.utils.r")
set.seed(2023)

feature.set = "hvf.10k"
#feature.set = "random.10k"


robust = T
qtl = 0.95

data.dir    = paste0("../Data/Methylation/Purified-Reinius/", feature.set)
res.dir     = paste0("../Result/Methylation/Purified-Reinius/", feature.set)
figure.dir  = paste0("../Figure/Purified-Reinius/", feature.set)

if (!file.exists(res.dir)){print("no result yet")}
if (!file.exists(figure.dir)){dir.create(file.path(figure.dir),recursive = T)}
print(data.dir)
print(res.dir)
print(figure.dir)

hannum  = readRDS(file.path(data.dir, paste0("hannum.", feature.set, ".rds")))
reinius = readRDS(file.path(data.dir, paste0("reinius.", feature.set, ".rds")))


max_stds = 2


base.mdl      = readRDS(file.path(res.dir, paste0("base.mdl.rds")))
cibersortx.mdl= readRDS(file.path(res.dir, paste0("cibersortx.mdl.rds")))
tca.mdl       = readRDS(file.path(res.dir, paste0("tca.mdl.rds")))
Unico.mdl      = readRDS(file.path(res.dir, paste0("Unico.mdl.rds")))
bMIND.mdl      = readRDS(file.path(res.dir, paste0("bMIND.mdl.rds")))

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

base.color       = brewer.pal(n = 8, name = "Set2")[8]
cibersortx.color = brewer.pal(n = 8, name = "Set2")[5]
Unico.color       = brewer.pal(n = 8, name = "Set2")[4]
tca.color        = brewer.pal(n = 8, name = "Set2")[6]
bMIND.color      = brewer.pal(n = 8, name = "Set2")[3]

model.list   = list("Baseline" = base.mdl, 
                    "CIBERSORTx" = cibersortx.mdl, 
                    'TCA'= tca.mdl, 
                    'bMIND' = bMIND.mdl,
                    'Unico' = Unico.mdl)
model.names  = c("Baseline", "CIBERSORTx", 'TCA', "bMIND", "Unico")
model.colors = c(base.color, cibersortx.color, tca.color,  bMIND.color, Unico.color)

str(Unico.mdl$eval$ metrics.list[[1]])

model.list[["Unico"]]$evals$metrics.list[[1]][["center.Z.corrs"]]

model.list[["bMIND"]]$evals$metrics.list[[1]][["center.Z.corrs"]]

model.list[["Unico"]]$evals$metrics.list[[1]][["center.Z.RMedS"]]

model.list[["bMIND"]]$evals$metrics.list[[1]][["center.Z.RMedS"]]

lab.size    = 20
title.size  = 20
legend.size = 20  

#model.list: list of mdl
#model.names: list of character
#plot.metric: character choose from c("center.Z.corrs", "center.Z.MAE", "center.Z.MedAE", "center.Z.RMS", "center.Z.RMedS")
#plot.group: character, choose from c("low entropy", "high entropy", "all")
#sample.itr: number of iterations to pool
get_bar_df = function(model.list, model.names, plot.metric, plot.group, sample.itr = 20){
    
    sample.itr = min(length(model.list[[1]]$evals$metrics.list), sample.itr)
    source.ids = rownames(model.list[[model.names[1]]]$evals$metrics.list[[1]][[1]])# first iter first eval metric
    
    
    plot.df = list()
    for (model.name in model.names){

        #construt a celltype by number of iteration matrix, each col is the desired metric in that desired group at at iteration 
        metrics.mat = vapply(1:sample.itr, function(itr) model.list[[model.name]]$evals$metrics.list[[itr]][[plot.metric]][,plot.group], numeric(length(source.ids)))

        #construt ploting related matrix for one method
        metrics.meta = matrix(0, length(source.ids), 5) 
        colnames(metrics.meta) = c("source", "mean", "lb_sd", "ub_sd", "method")
                                                    
        metrics.meta[,1] = source.ids
        metrics.meta[,2] = rowMeans(metrics.mat)
        metrics.meta[,3] = rowMeans(metrics.mat) - rowSds(metrics.mat)
        metrics.meta[,4] = rowMeans(metrics.mat) + rowSds(metrics.mat)
        metrics.meta[,5] = rep(model.name, nrow(metrics.meta))
        plot.df[[model.name]] = metrics.meta
    }
    plot.df = as.data.frame(Reduce(rbind, plot.df))
    plot.df$mean   =  as.numeric(plot.df$mean)   
    plot.df$lb_sd  =  as.numeric(plot.df$lb_sd)      
    plot.df$ub_sd  =  as.numeric(plot.df$ub_sd)   
                             
    plot.df$method = factor(plot.df$method, levels = model.names)  
    plot.df$source = factor(plot.df$source, levels = source.ids)
    levels(plot.df$method)[match("Unico",levels(plot.df$method))] <- "Unico"
    return (plot.df)
}

source.ids = colnames(reinius$W)
sub.titles = c("Granulocytes", "CD8 T Cells", "CD4 T Cells", "Monocytes", "NK Cells", "B Cells")
message("making sure the source.ids align with the full name one line below")
print(source.ids)
print(sub.titles)

cor.bar.dfs = list()
cor.bar.plots = list()
for (group.name in c("all", "low entropy", "high entropy")){
    if(group.name == "low entropy"){
        title = "Low Entropy"
    }else if(group.name == "high entropy"){
        title = "High Entropy"
    }else{
        title = "All CpGs"
    }

    cor.bar.df = get_bar_df(model.list, model.names, plot.metric = "center.Z.corrs", 
                         plot.group = group.name, sample.itr = 20)


    cor.bar.g = ggplot(cor.bar.df, aes(x=as.factor(source), y=mean, fill=method)) +
                geom_bar(position=position_dodge(.9), stat="identity", colour='black') + 
                geom_errorbar(aes(ymin = lb_sd, ymax= ub_sd), width=.5, position=position_dodge(.9)) +
                coord_cartesian(ylim = c(0, max(cor.bar.df$ub_sd))) + 
                scale_fill_manual(values = model.colors) + 

                theme_classic() +
                ggtitle(title) +
                scale_x_discrete(labels=sub.titles) +

                xlab(paste0("Cell Types")) + 
                ylab(paste0("Correlation")) + 

           
                theme(plot.title = element_text(hjust = 0.5, size = title.size),
                      text=element_text(size=lab.size),
                      axis.title.x = element_blank(),
                      axis.text.x  = element_text(angle = 25, vjust = 0.5, hjust = 0.5),
                      legend.title = element_blank())                                 
              

    cor.bar.plots[[group.name]] = cor.bar.g
    cor.bar.dfs[[group.name]]   = cor.bar.df
}

options(repr.plot.width = 6, repr.plot.height = 5, repr.plot.res = 150)
cor.bar.plots[["low entropy"]]

cor.bar.plots[["high entropy"]]

RMedS.bar.dfs = list()
RMedS.bar.plots = list()
for (group.name in c("all", "low entropy", "high entropy")){
    if(group.name == "low entropy"){
        title = "Low Entropy"
    }else if(group.name == "high entropy"){
        title = "High Entropy"
    }else{
        title = "All CpGs"
    }

    RMedS.bar.df = get_bar_df(model.list, model.names, plot.metric = "center.Z.RMedS", 
                         plot.group = group.name, sample.itr = 20)


    RMedS.bar.g = ggplot(RMedS.bar.df, aes(x=as.factor(source), y=mean, fill=method)) +
                  geom_bar(position=position_dodge(.9), stat="identity", colour='black') + 
                  geom_errorbar(aes(ymin = lb_sd, ymax= ub_sd), width=.5, position=position_dodge(.9)) +
                  scale_fill_manual(values = model.colors) + 

                  theme_classic() +
                  ggtitle(title) +
                  scale_x_discrete(labels=sub.titles) +

                  xlab(paste0("Cell Types")) + 
                  ylab(paste0("RMSE")) + 

                  theme(plot.title = element_text(hjust = 0.5, size = title.size),
                      text=element_text(size=lab.size),
                      axis.title.x = element_blank(),
                      axis.text.x  = element_text(angle = 25, vjust = 0.5, hjust = 0.5),
                      legend.title = element_blank())                                 
              
    RMedS.bar.plots[[group.name]] = RMedS.bar.g
    RMedS.bar.dfs[[group.name]]   = RMedS.bar.df
}

RMedS.bar.plots[["low entropy"]]

RMedS.bar.plots[["high entropy"]]

RMS.bar.dfs = list()
RMS.bar.plots = list()
for (group.name in c("all", "low entropy", "high entropy")){
    if(group.name == "low entropy"){
        title = "Low Entropy"
    }else if(group.name == "high entropy"){
        title = "High Entropy"
    }else{
        title = "All CpGs"
    }

    RMS.bar.df = get_bar_df(model.list, model.names, plot.metric = "center.Z.RMS", 
                         plot.group = group.name, sample.itr = 20)


    RMS.bar.g = ggplot(RMS.bar.df, aes(x=as.factor(source), y=mean, fill=method)) +
                  geom_bar(position=position_dodge(.9), stat="identity", colour='black') + 
                  geom_errorbar(aes(ymin = lb_sd, ymax= ub_sd), width=.5, position=position_dodge(.9)) +
                  scale_fill_manual(values = model.colors) + 

                  theme_classic() +
                  ggtitle(title) + 
                  scale_x_discrete(labels=sub.titles) +

                  xlab(paste0("Cell Types")) + 
                  ylab(paste0("RMSE")) + 

                  theme(plot.title = element_text(hjust = 0.5, size = title.size),
                      text=element_text(size=lab.size),
                      axis.title.x = element_blank(),
                      axis.text.x  = element_text(angle = 25, vjust = 0.5, hjust = 0.5),
                      legend.title = element_blank())                                 
         
    RMS.bar.plots[[group.name]] = RMS.bar.g
    RMS.bar.dfs[[group.name]]   = RMS.bar.df
}

RMS.bar.plots[["low entropy"]] 

RMS.bar.plots[["high entropy"]] 

legend.g = as_ggplot(get_legend(cor.bar.plots[["low entropy"]] + 
                                theme(legend.direction = "horizontal",
                                      legend.position = c(0.575, 0.5),
                                      legend.text = element_text(size=legend.size, 
                                                                 margin = margin(r = 100, unit = "pt"))))) 


options(repr.plot.width = 16, repr.plot.height = 0.5, repr.plot.res = 100)
legend.g 

cor.bar.g = ggarrange(cor.bar.plots[["low entropy"]]   + theme(legend.position = "none"),
                      cor.bar.plots[["high entropy"]]  + theme(legend.position = "none", axis.title.y = element_blank()) , 
                      nrow = 1,
                      widths = c(1, 1))

rmse.bar.g = ggarrange(RMedS.bar.plots[["low entropy"]]   + theme(legend.position = "none"),
                       RMedS.bar.plots[["high entropy"]]  + theme(legend.position = "none", axis.title.y = element_blank()),
                       nrow = 1,
                       widths = c(1, 1))
                      
                       

multi.g = ggarrange(rmse.bar.g, cor.bar.g,
                    nrow = 1,
                    widths = c(1, 1),
                    labels = c("a","b"),
                    font.label = list(size = 30, color = "black", face = "bold", family = NULL))

multi.g = ggarrange(legend.g, multi.g,
                    nrow = 2,
                    heights = c(0.15, 1))
                    

options(repr.plot.width = 16, repr.plot.height = 6, repr.plot.res = 200)
multi.g

ggsave(file.path(figure.dir, paste0("Reinius_eval_multi_", feature.set, ".pdf")), 
       multi.g,
       device = "pdf", width = 16, height = 6, dpi = 600)

fig.list = list(cor.bar.plots = cor.bar.plots, 
                RMedS.bar.plots = RMedS.bar.plots,
                RMS.bar.plots = RMS.bar.plots)

stats.list = list(cor.bar.dfs = cor.bar.dfs,
                  RMS.bar.dfs = RMS.bar.dfs, 
                  RMedS.bar.dfs = RMedS.bar.dfs)


saveRDS(fig.list,   file.path(figure.dir, paste0("Reinius_eval_", feature.set, "_fig.list.rds")))
saveRDS(stats.list, file.path(figure.dir, paste0("Reinius_eval_", feature.set, "_stats.list.rds")))

# this is seperate from the barplots
lab.size = 16
title.size = 15

#after adjust for source, feature specific mean, 
#collapse all features to one meta feature
#return a list with two key, ref and hat. 
#each key is just a long meta vector, essentially one meta feature for this source
calc_centered_Z_h = function(Z.ref, Z.hat, h, mask){
    n = dim(Z.ref)[3]
    mask = as.vector(mask)
    
    #source, feature specific mean
    mus.ref = as.matrix(rowMeans(Z.ref[h,mask,]))
    mus.hat = as.matrix(rowMeans(Z.hat[h,mask,]))
    
    #collapse
    ref = as.vector(Z.ref[h,mask,] - repmat(mus.ref, 1, n))
    hat = as.vector(Z.hat[h,mask,] - repmat(mus.hat, 1, n))
    return (data.frame(ref = ref, hat = hat))  
}

#for density
calc_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, 
                      h = c(ifelse(bandwidth.nrd(x) == 0, 0.001, bandwidth.nrd(x)),
                            ifelse(bandwidth.nrd(y) == 0, 0.001, bandwidth.nrd(y))), ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

feature.ids = dimnames(reinius$Z)[[2]]
source.ids  = dimnames(reinius$Z)[[1]]
k = length(source.ids)   

group.mat = matrix(FALSE, length(feature.ids), 3)
rownames(group.mat) = feature.ids
colnames(group.mat) = c("low entropy", "high entropy", "all")
group.mat[,"low entropy"]  = reinius$params$entropies <= median(reinius$params$entropies)
group.mat[,"high entropy"] = reinius$params$entropies >  median(reinius$params$entropies)
group.mat[,"all"] = TRUE

#group.name = "all"
for (group.name in c("low entropy", "high entropy", "all")){ 

    plts <- vector(mode = "list", length = length(model.names) * k)
    counter  = 1
    #loop over source.ids
    for (h in 1:k){ 
        #loop over models
        for (t in 1:length(model.names)){
            
            model.name = model.names[t]
            Z.hat =  model.list[[model.name]]$Z.hat
            df = calc_centered_Z_h(reinius$Z, Z.hat, h, group.mat[,group.name])

            # remove outliers
            df = as.matrix(df)
            df = as.data.frame(df[(df[,"ref"] > 1.25 * quantile(df[,"ref"], 0.01)) & 
                                  (df[,"ref"] < 1.25 * quantile(df[,"ref"], 0.99)) & 
                                  (df[,"hat"] > 1.25 * quantile(df[,"hat"], 0.01)) &
                                  (df[,"hat"] < 1.25 * quantile(df[,"hat"], 0.99)),])

            # for space of annotation on top
            corr = round(safe_cor(df$ref, df$hat, robust, qtl=qtl), 2)
            
            
            # subsample to plot
            df = df[sample(nrow(df),500), ] # subsample 500 dots
            df$density <- calc_density(df$ref, df$hat, n = 100) # larger -> more smooth 

            coord.min = 1.25 * min(quantile(df$ref, 0.01), quantile(df$hat, 0.01)) 
            coord.max = 1.25 * max(quantile(df$ref, 0.99), quantile(df$hat, 0.99)) 
            
            #top row needs to have more space for title/method name
            if (counter < length(model.list)){
                ylim = c(coord.min, coord.max * 1.2)
            }

            g = ggplot(df, aes(x = ref, y = hat, color = density)) +
                geom_point(alpha = 1, size = 1.5) + 

                scale_color_gradient() + 
                geom_abline(slope=1, colour= "black", linetype = "dashed", size = 1) + 
                scale_x_continuous(n.breaks = 4) + 
                scale_y_continuous(n.breaks = 4) + 

                coord_cartesian(xlim = c(coord.min, coord.max), 
                                ylim = c(coord.min, coord.max))+
            
                #coord_cartesian(xlim = c(1.25 * quantile(df$ref, 0.01), 1.25 * quantile(df$ref, 0.99)), 
                #                ylim = c(1.25 * quantile(df$hat, 0.01), 1.25 * quantile(df$hat, 0.99))) + 
                
            
                #clean background 
                theme_classic()+

                ggtitle(paste0("r = ",corr)) +
                theme(plot.title = element_text(hjust = 0.5, size = title.size))+

                ylab(paste0("Deconvolution")) + 
                xlab(paste0("Ground truth")) +
                theme(text=element_text(size=lab.size))+ 
                theme(axis.title.x = element_text(size = lab.size)) +
                theme(axis.text.x  = element_text(size = lab.size)) + 
                theme(axis.title.y = element_text(size = lab.size)) +
                theme(axis.text.y  = element_text(size = lab.size)) 



            #if not the first method, turn off y label 
            if(t != 1){
                g = g + theme(axis.title.y  = element_blank())
            }
            #if not the last source, turn off x axis 
            if(h != length(source.ids)){
                g = g + theme(axis.title.x  = element_blank())
            }

            plts[[counter]] = g
            counter = counter + 1

        }
    }


    #add title/method names on the column
    for (counter in 1:length(model.list)){
        model.name = model.names[counter]
        if(model.name == "Unico"){model.name = "Unico"}
        plts[[counter]] = annotate_figure(plts[[counter]],
                                          top = text_grob(model.name, #vjust = 0.5,  
                                                          color = "black", face = "bold", size = title.size + 5))
    
    }

    #add source.ids on the row
    counter = 1
    for (h in 1:length(source.ids)){

        plts[[counter]] = annotate_figure(plts[[counter]],
                                          left = text_grob(paste0(source.ids[h], "(", round(colMeans(reinius$W)[h], 2) * 100, "%)"), 
                                                           rot = 90, #vjust = 0.5,  
                                                           color = "black", face = "bold", size = title.size + 5))
        counter = counter + length(model.names)
    }



    g <- ggarrange(plotlist = plts, ncol = length(model.names), nrow = length(source.ids), 
                   widths  = c(1.1,rep (1, (length(model.names) - 1))),
                   heights = c(rep (1, (length(source.ids) - 1), 1.1)))

    
    ggsave(file.path(figure.dir,paste0("Reinius_", gsub(" ",".",group.name), "_", feature.set, "_scatter.pdf")), 
           g,
           device = "pdf", width = 18, height = 17, dpi = 600)    
    
} 

options(repr.plot.width = 18, repr.plot.height = 17, repr.plot.res = 200)
g 
