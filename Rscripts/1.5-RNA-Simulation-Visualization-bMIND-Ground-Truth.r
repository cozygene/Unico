library("ggplot2")
library("ggpubr")

library("matrixStats")# colMeans
source("analysis.utils.r")
source("simulate.expression.utils.r")

set.seed(2023)

data.name = "PBMC"
#data.name = "Lung"

source.col = "decon.L1" 
N = 500

max_stds = 2
ts = 1:20

project.dir = "/u/home/j/johnsonc/project-halperin/Unico/Unico2023/"

#relative path of the data dir to project.dir
if(data.name == "PBMC"){
    if(source.col == "decon.L1"){
        data.dir = paste0("Data/RNA/Simulation-PBMC/sc-Stephenson_all.W_decon.L1_HEF.10k_k_5_m_600_n_", N, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
        figure.dir  = paste0("Figure/Simulation-PBMC/decon.L1.max_stds.",max_stds)
    }else{
        data.dir = paste0("Data/RNA/Simulation-PBMC/sc-Stephenson_all.W_decon.L2_HEF.10k_k_7_m_600_n_", N, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
        figure.dir  = paste0("Figure/Simulation-PBMC/decon.L2.max_stds.",max_stds)
    }
}else{
    if(source.col == "decon.L1"){
        data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L1_HEF.10k_k_4_m_600_n_", N, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
        figure.dir  = paste0("Figure/Simulation-Lung/decon.L1.max_stds.",max_stds)
    }else{
        data.dir = paste0("Data/RNA/Simulation-Lung/sc-HLCA_all.W_decon.L2_HEF.10k_k_6_m_600_n_", N, "_dirichlet_F_noiseZ_T_varThr_1e-04_filThr_1e-04_expQtl_0_enrich_F_etpRat_0_enrichRat_0_maxSds_", max_stds, "_scale.maxSds_Inf_scale.factor.thr_1e-04/")
        figure.dir  = paste0("Figure/Simulation-Lung/decon.L2.max_stds.",max_stds)
    }

}

data.dir   = file.path(project.dir, data.dir)
figure.dir = file.path(project.dir, figure.dir)
data.dir.split  = strsplit(data.dir, "/")[[1]]
res.dir    = file.path(project.dir, paste(c("Result", data.dir.split[(length(data.dir.split)-2) : length(data.dir.split)]), collapse = "/"))

if (!file.exists(res.dir)){print("no result in the result directory")}
if (!file.exists(figure.dir)) {dir.create(figure.dir, recursive = T)}
print(data.dir)
print(res.dir)
print(figure.dir)

sim.data.list = readRDS(file.path(data.dir, "sim.data.list.rds"))

# load models
bMIND_rp.mdl.list = readRDS(file.path(res.dir, paste0("bMIND.mdl.rough.prior.list.rds")))
bMIND_sc.mdl.list = readRDS(file.path(res.dir, paste0("bMIND.mdl.sc.prior.1.list.rds")))
Unico.mdl.list = readRDS(file.path(res.dir, paste0("Unico.mdl.list.rds")))

title.size = 20
lab.size = 20
legend.size = 20  

source.ids = colnames(sim.data.list[[1]]$W)
print(source.ids)
k = length(source.ids)
m = nrow(sim.data.list[[1]]$X)
n = ncol(sim.data.list[[1]]$X)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

Unico.color       = brewer.pal(n = 8, name = "Set2")[4]
bMIND_sc.color   = brewer.pal(n = 11, name = "RdYlBu")[10]
bMIND_rp.color   = brewer.pal(n = 8, name = "Set2")[3]

#keys are the method
#per key a matrix, rows are the runs
mus.corrs.list = list(
    Unico       = concat_2_keys(Unico.mdl.list,       key1 ="moment.hat.corrs",   key2 = "mus.rob.corrs"),
    bMIND_rp   = concat_2_keys(bMIND_rp.mdl.list,   key1 ="moment.hat.corrs",   key2 = "mus.rob.corrs"),
    bMIND_sc   = concat_2_keys(bMIND_sc.mdl.list,   key1 ="moment.hat.corrs",   key2 = "mus.rob.corrs")
)
var.corrs.list = list(
    Unico       = concat_2_keys(Unico.mdl.list,       key1 ="moment.hat.corrs",   key2 = "var.rob.corrs"),
    bMIND_rp   = concat_2_keys(bMIND_rp.mdl.list,   key1 ="moment.hat.corrs",   key2 = "var.rob.corrs"),
    bMIND_sc   = concat_2_keys(bMIND_sc.mdl.list,   key1 ="moment.hat.corrs",   key2 = "var.rob.corrs")
)
covar.corrs.list = list(
    Unico       = concat_2_keys(Unico.mdl.list,       key1 ="moment.hat.corrs",   key2 = "covar.rob.corrs"),
    bMIND_rp   = concat_2_keys(bMIND_rp.mdl.list,   key1 ="moment.hat.corrs",   key2 = "covar.rob.corrs"),
    bMIND_sc   = concat_2_keys(bMIND_sc.mdl.list,   key1 ="moment.hat.corrs",   key2 = "covar.rob.corrs")
)
covar.ids = colnames(covar.corrs.list[[1]])

methods = c("bMIND_rp", "bMIND_sc", "Unico")
method.colors = c(bMIND_rp.color, bMIND_sc.color, Unico.color)

mean.barplot.df = data.frame(
    method = as.vector(vapply(1:length(methods), function(a) rep(methods[a], k), character(k))),
    source = as.vector(rep(source.ids, length(methods))),
    mean   = as.vector(vapply(1:length(methods), function(a) colMeans(mus.corrs.list[[methods[a]]]), numeric(k))),
    lb_sd  = as.vector(vapply(1:length(methods), function(a) colMeans(mus.corrs.list[[methods[a]]]) - colSds(mus.corrs.list[[methods[a]]]),  numeric(k))),
    ub_sd  = as.vector(vapply(1:length(methods), function(a) colMeans(mus.corrs.list[[methods[a]]]) + colSds(mus.corrs.list[[methods[a]]]),  numeric(k)))
)

# set some to factors so that the figures respect the order 
mean.barplot.df$source = factor(mean.barplot.df$source, levels = source.ids)
mean.barplot.df$method = factor(mean.barplot.df$method, levels = methods)
levels(mean.barplot.df$method)[match("Unico",levels(mean.barplot.df$method))] <- "Unico"

str(mus.corrs.list)

head(mean.barplot.df)

mean.bar.g = ggplot(mean.barplot.df, aes(x=as.factor(source), y=mean, fill=method)) +
             geom_bar(position=position_dodge(.9), stat="identity", colour='black') + 
             geom_errorbar(aes(ymin = lb_sd, ymax= ub_sd), width=.5, position=position_dodge(.9)) +
             scale_fill_manual(values = method.colors) + 
             coord_cartesian(ylim = c(0, max(mean.barplot.df$ub_sd))) +  
             ggtitle("Means") + 
             xlab(paste0("Cell Types")) + 
             ylab(paste0("Correlation")) + 
             theme_classic() +  
             theme(plot.title = element_text(hjust = 0.5, size = title.size),
                   text=element_text(size=lab.size), #text size
                   axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1), #rotate
                   axis.title.x=element_blank(), # no x tilte
                   legend.title = element_blank())
 
options(repr.plot.width = 6, repr.plot.height = 4, repr.plot.res = 100)
mean.bar.g

var.barplot.df = data.frame(
    method = as.vector(vapply(1:length(methods), function(a) rep(methods[a], k), character(k))),
    source = as.vector(rep(source.ids, length(methods))),
    mean   = as.vector(vapply(1:length(methods), function(a) colMeans(var.corrs.list[[methods[a]]]), numeric(k))),
    lb_sd  = as.vector(vapply(1:length(methods), function(a) colMeans(var.corrs.list[[methods[a]]]) - colSds(var.corrs.list[[methods[a]]]),    numeric(k))),
    ub_sd  = as.vector(vapply(1:length(methods), function(a) colMeans(var.corrs.list[[methods[a]]]) + colSds(var.corrs.list[[methods[a]]]),    numeric(k)))
)

# set some to factors so that the figures respect the order 
var.barplot.df$source = factor(var.barplot.df$source, levels = source.ids)
var.barplot.df$method = factor(var.barplot.df$method, levels = methods)  
levels(var.barplot.df$method)[match("Unico",levels(var.barplot.df$method))] <- "Unico"

head(var.barplot.df)

var.bar.g = ggplot(var.barplot.df, aes(x=as.factor(source), y=mean, fill=method)) +
            geom_bar(position=position_dodge(.9), stat="identity", colour='black') + 
            geom_errorbar(aes(ymin = lb_sd, ymax= ub_sd), width=.5, position=position_dodge(.9)) +  
            coord_cartesian(ylim = c(0, max(var.barplot.df$ub_sd))) + 
            scale_fill_manual(values = method.colors) + 
            ggtitle("Variances") + 
            xlab(paste0("Cell Types")) + 
            ylab(paste0("Correlation")) + 
            theme_classic() + 
            theme(plot.title = element_text(hjust = 0.5, size = title.size),
                   text=element_text(size=lab.size), #text size
                   axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1), #rotate
                   axis.title.x=element_blank(), # no x tilte
                   legend.title = element_blank())
 
options(repr.plot.width = 6, repr.plot.height = 6, repr.plot.res = 100)
var.bar.g

covar.barplot.df = data.frame(
    method = as.vector(vapply(1:length(methods), function(a) rep(methods[a], length(covar.ids)), character(length(covar.ids)))),
    source = as.vector(rep(covar.ids, length(methods))),
    mean   = as.vector(vapply(1:length(methods), function(a) colMeans(covar.corrs.list[[methods[a]]]), numeric(length(covar.ids)))),
    lb_sd  = as.vector(vapply(1:length(methods), function(a) colMeans(covar.corrs.list[[methods[a]]]) - colSds(covar.corrs.list[[methods[a]]]),    numeric(length(covar.ids)))),
    ub_sd  = as.vector(vapply(1:length(methods), function(a) colMeans(covar.corrs.list[[methods[a]]]) + colSds(covar.corrs.list[[methods[a]]]),    numeric(length(covar.ids))))
)

# extract top TCA perform
covar.ids = covar.ids[order(-covar.barplot.df[covar.barplot.df["method"] == "bMIND_rp","mean"])]
                       
# set some to factors so that the figures respect the order 
covar.barplot.df$source = factor(covar.barplot.df$source, levels = covar.ids)
covar.barplot.df$method = factor(covar.barplot.df$method, levels = methods)     
levels(covar.barplot.df$method)[match("Unico",levels(covar.barplot.df$method))] <- "Unico"

head(covar.barplot.df)

covar.bar.g = ggplot(covar.barplot.df, aes(x=as.factor(source), y=mean, fill=method)) +
              geom_bar(position=position_dodge(.9), stat="identity", colour='black') + 
              
              geom_errorbar(aes(ymin = lb_sd, ymax= ub_sd), width=.5, position=position_dodge(.9)) +                                                
              coord_cartesian(ylim = c(0, max(covar.barplot.df$ub_sd))) + 
              scale_fill_manual(values = method.colors) + 
              ggtitle("Covariances") + 
              xlab(paste0("Cell Type Pairs")) + 
              ylab(paste0("Correlation")) + 
              theme_classic() + 
              theme(plot.title = element_text(hjust = 0.5, size = title.size),
                    text=element_text(size=lab.size), #text size
                    axis.text.x  = element_text(angle = 25, vjust = 1, hjust = 1), #rotate
                    axis.title.x = element_blank(), # no x tilte
                    legend.title = element_blank())
 
options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 100)
covar.bar.g

plts = list(3)
plts[[1]] = mean.bar.g  + theme(legend.position="none") 
plts[[2]] = var.bar.g   + theme(legend.position="none") + theme(axis.title.y = element_blank()) 
plts[[3]] = covar.bar.g + theme(legend.position="none") + theme(axis.title.y = element_blank())

#too many celltypes. shrink the size
if(source.col == "decon.L2"){
    plts[[1]] = plts[[1]] +  theme(text=element_text(size=lab.size - 2.5))
    plts[[2]] = plts[[2]] +  theme(text=element_text(size=lab.size - 2.5))
    plts[[3]] = plts[[3]] +  theme(text=element_text(size=lab.size - 2.5))
}

covar.width   = round(2 * length(covar.ids)/length(source.ids))
#figure.width  = (4 + covar.width) * 3
#figure.height = 7 

options(repr.plot.width = 16, repr.plot.height = 3, repr.plot.res = 200)
params.bar.g = egg::ggarrange(plots=list(plts[[1]], plts[[2]], plts[[3]]),  
                    #labels = c("a", "", ""), 
                    align = "h",
                    widths = c(2,2,covar.width),  
                    #label.args = list(gp = grid::gpar(font = 20, cex =3)), 
                    debug=F)

for (t in ts){
    sim.data.list[[t]]$eval.feature.source = calc_variable_feature_source(sim.data.list[[t]]$Z.scale, 
                                                                          variable_thr = 0.1, max_sds = max_stds)
    
    #sim.data.list[[t]]$eval.entropies = sim.data.list[[t]]$params$entropies <= sim.data.list[[t]]$data.gen.params$entropy.thr
    # arbitray cut to low and high entropy group
    sim.data.list[[t]]$eval.entropies = sim.data.list[[t]]$params$entropies <= median(sim.data.list[[t]]$params$entropies)
    sim.data.list[[t]]$eval.entropies = as.matrix(vapply(1:m, function(j) if (sim.data.list[[t]]$eval.entropies[j]) "Low Entropy" else "High Entropy", character(1)))
}

Z.corrs.list = list(
    #m*t by k 
    Unico           = concat_key(Unico.mdl.list,       key = "Z.corrs"),
    bMIND_rp       = concat_key(bMIND_rp.mdl.list,   key = "Z.corrs"),
    bMIND_sc       = concat_key(bMIND_sc.mdl.list,   key = "Z.corrs"),
    
    mask           = concat_key(sim.data.list,       key = "eval.feature.source"),
    #m*t by 1
    eval.entropies = concat_key(sim.data.list,       key = "eval.entropies")
)

median(sim.data.list[[t]]$params$entropies)

hist(sim.data.list[[t]]$params$entropies)

table(sim.data.list[[t]]$eval.entropies)

boxplot.meta = list()
for (source.id in source.ids){
    
    median.mat       = matrix(0, 2, length(methods))
    rownames(median.mat) = c("High Entropy", "Low Entropy")
    colnames(median.mat) = methods
    
    whitney.pval.mat = copy(median.mat)
    binom.pval.mat   = copy(median.mat)
    n.success.mat    = copy(median.mat)
    p.success.mat    = copy(median.mat)
    
    for (entropy in c("High Entropy" ,"Low Entropy")){
        
        # selected gene in the ideal entropy group and deemed to be evaluated 
        keep.mask = (Z.corrs.list$eval.entropies == entropy) & Z.corrs.list$mask[,source.id]
        
        for (method in methods){
            method.corrs = Z.corrs.list[[method]][keep.mask, source.id]
            Unico.corrs   = Z.corrs.list[["Unico"]][keep.mask, source.id]
            
            #paired one side non-parametric whitney test
            median.mat[entropy, method] = median(method.corrs)
            res = wilcox.test(x = Unico.corrs, 
                              y = method.corrs, 
                              alternative = "g", paired = T) # x > y one side test 
            whitney.pval.mat[entropy,method] = res$p.value
            
            #binomial test
            n_s = sum(Unico.corrs > method.corrs)                          
            n   = sum(keep.mask)
            res = binom.test(x = n_s, n = n, p = 0.5, alternative = c("greater"))
            binom.pval.mat[entropy,method] = res$p.value
            n.success.mat[entropy,method]  = n_s
            p.success.mat[entropy,method]  = n_s/n
            
        }
        boxplot.meta[[source.id]] = list(median.mat       = median.mat,
                                         whitney.pval.mat = whitney.pval.mat,
                                         binom.pval.mat   = binom.pval.mat,
                                         n.success.mat    = n.success.mat,
                                         p.success.mat    = p.success.mat)
    }
}

method

source.id

head(Z.corrs.list[["Unico"]])

boxplot.meta

boxplot.meta[[1]]$whitney.pval.mat

if(data.name == "PBMC"){
    if (source.col == "decon.L1"){
        sub.titles = c("CD4 T Cells", "NK Cells", "CD8 T Cells", "Monocytes", "B Cells")
    }else{
        sub.titles = c("CD4 T Cells", "NK Cells", "CD8 T Cells", "Monocytes (CD14)", "B Cells", "Monocytes (CD16)", "Plasma" )
    }
}else{
    if (source.col == "decon.L1"){
        sub.titles = source.ids
    }else{
        sub.titles = source.ids
    }
}

print(source.ids)
print(sub.titles)

h.loc = linspace(0.7, 1.3, 3)
v.loc = linspace(1.0, 1.1, 2)
v.size = 0.02
p.size = 4

under.flow = "< 1e-300"
#under.flow = "***"

W = concat_key(sim.data.list, "W")
plts <- vector(mode = "list", length = k)
m.m = nrow(Z.corrs.list$eval.entropies)
for (h in 1:k){
    Z.boxplot.df <- data.frame(cor     = as.vector(vapply(1:length(methods), function(a) Z.corrs.list[[methods[a]]][,h],   numeric(m.m))),
                               mask    = as.vector(vapply(1:length(methods), function(a) Z.corrs.list[["mask"]][,h],       logical(m.m))),
                               entropy = as.vector(vapply(1:length(methods), function(a) Z.corrs.list[["eval.entropies"]], character(m.m))),
                               method  = as.vector(vapply(1:length(methods), function(a) rep(methods[a], m.m),             character(m.m))))
    # set some to factors so that the figures respect the order                                                     
    Z.boxplot.df$method  = factor(Z.boxplot.df$method,  levels = methods) 
    Z.boxplot.df$entropy = factor(Z.boxplot.df$entropy, levels = c("High Entropy", "Low Entropy"))
    levels(Z.boxplot.df$method)[match("Unico",levels(Z.boxplot.df$method))] <- "Unico"
                                                          
    #keep only those gene-celltype that has evaluation turned on 
    Z.boxplot.df = Z.boxplot.df[Z.boxplot.df$mask, ]
    print(dim(Z.boxplot.df))                                                 
                                                     
    plts[[h]] <- ggplot(Z.boxplot.df, aes(x = entropy, y=cor)) +
                 geom_boxplot(aes(fill = method), position = position_dodge(0.9))+    
                 scale_fill_manual (values = method.colors) +
                 scale_color_manual(values = method.colors) +
                                                                
                 coord_cartesian(xlim = NULL, ylim = c(0,1 + (length(methods)-1)/10)) +  
                 xlab(paste0(" "))+                                    
                 ylab(paste0("Correlation"))+
                 scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
                                                          
                 ggtitle(paste0(sub.titles[h], " (", round(mean(W[,h]) * 100), "%)")) + 
                 theme_classic() + 
                 theme(plot.title = element_text(hjust = 0.5, size = title.size),
                       text=element_text(size=lab.size),
                       legend.title = element_blank())
                                                     
    #adding p-vales
    whitney.pval.mat = formatC(boxplot.meta[[h]]$whitney.pval.mat, format = "e", digits = 1) 
    whitney.pval.mat[which(whitney.pval.mat == "0.0e+00")] = under.flow
          
    df1 <- data.frame(a = c(h.loc[1],h.loc[1], h.loc[3], h.loc[3]), b = c(v.loc[1]-v.size, v.loc[1], v.loc[1], v.loc[1]-v.size))
    df2 <- data.frame(a = c(h.loc[2],h.loc[2], h.loc[3], h.loc[3]), b = c(v.loc[2]-v.size, v.loc[2], v.loc[2], v.loc[2]-v.size))
    
    df3 <- data.frame(a = 1 + c(h.loc[1],h.loc[1], h.loc[3], h.loc[3]), b = c(v.loc[1]-v.size, v.loc[1], v.loc[1], v.loc[1]-v.size))
    df4 <- data.frame(a = 1 + c(h.loc[2],h.loc[2], h.loc[3], h.loc[3]), b = c(v.loc[2]-v.size, v.loc[2], v.loc[2], v.loc[2]-v.size))
    
    plts[[h]] <- plts[[h]] + geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = (h.loc[1]+ h.loc[3])/2,   y = v.loc[1] + 2*v.size, label = whitney.pval.mat[1,1], size = p.size) +
                             geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = (h.loc[2]+ h.loc[3])/2,   y = v.loc[2] + 2*v.size, label = whitney.pval.mat[1,2], size = p.size) +
                                                          
                             geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 1+(h.loc[1]+ h.loc[3])/2, y = v.loc[1] + 2*v.size, label = whitney.pval.mat[2,1], size = p.size) +
                             geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 1+(h.loc[2]+ h.loc[3])/2, y = v.loc[2] + 2*v.size, label = whitney.pval.mat[2,2], size = p.size) 
 
    plts[[h]] = plts[[h]] + theme(legend.position="none") 
}

Z.boxplot.df

options(repr.plot.width = 16, repr.plot.height = 10, repr.plot.res = 200)
cor.box.g = egg::ggarrange(plots = plts,
                           #labels = c("b", rep("", (length(source.ids) ))),
                           align = "h",
                           widths = rep(2, ceil(k/2)),  
                           #label.args = list(gp = grid::gpar(font = 20, cex =3)), 
                           debug=F)


# might see underflow if Z.hat has pval 1 due to not fitted, and bulk has p val whose -log10 value smaller than 300 
min.log10p.diff = -300

#list of matricies
#each entry is a method
joint.bulk.Z.hat.log10p.diff.list = list(
    #m*t by k 
    Unico           = concat_key(Unico.mdl.list,       key = "joint.bulk.Z.hat.log10p.diff"),
    bMIND_rp       = concat_key(bMIND_rp.mdl.list,   key = "joint.bulk.Z.hat.log10p.diff"),
    bMIND_sc       = concat_key(bMIND_sc.mdl.list,   key = "joint.bulk.Z.hat.log10p.diff"),
    mask           = concat_key(sim.data.list,       key = "eval.feature.source"),
    #m*t by 1
    eval.entropies = concat_key(sim.data.list,       key = "eval.entropies")
)

names(joint.bulk.Z.hat.log10p.diff.list )

#those that have too small variance are ignored
joint.bulk.Z.hat.log10p.diff.list$Unico

plot.list = joint.bulk.Z.hat.log10p.diff.list

plot.list["Unico"]

addition.box.dfs = data.frame()
m.m = nrow(plot.list$eval.entropies)
for (h in 1:k){
    boxplot.df <- data.frame(source.id = as.vector(vapply(1:length(methods), function(a) rep(source.ids[h], m.m),       character(m.m))),
                             y         = as.vector(vapply(1:length(methods), function(a) plot.list[[methods[a]]][,h],   numeric(m.m))),
                             mask      = as.vector(vapply(1:length(methods), function(a) plot.list[["mask"]][,h],       logical(m.m))),
                             entropy   = as.vector(vapply(1:length(methods), function(a) plot.list[["eval.entropies"]], character(m.m))),
                             method    = as.vector(vapply(1:length(methods), function(a) rep(methods[a], m.m),          character(m.m))))
    #underflow issue
    boxplot.df[which(boxplot.df$y == -Inf), "y"] = min.log10p.diff
                                                          
    # set some to factors so that the figures respect the order                                                     
    boxplot.df$method  = factor(boxplot.df$method,  levels = methods) 
    boxplot.df$entropy = factor(boxplot.df$entropy, levels = c("High Entropy", "Low Entropy"))
      
    #keep only those gene-celltype that has evaluation turned on 
    boxplot.df = boxplot.df[boxplot.df$mask, ]   
    addition.box.dfs = rbind(addition.box.dfs, boxplot.df)                                                                             
}

# update source to be more informative, change it to full names in sub.titles and add proportion
sub.w.titles = c()
W   = concat_key(sim.data.list, "W")
for (l in 1:k){   
    sub.w.title = paste0(sub.titles[l], " (", round(mean(W[,l]) * 100), "%)")
    addition.box.dfs[addition.box.dfs$source.id == source.ids[l], "source.id"] = sub.w.title
    sub.w.titles = c(sub.w.titles, sub.w.title)
}
addition.box.dfs$source.id = factor(addition.box.dfs$source.id, levels = sub.w.titles)

# name swap to Unico
levels(addition.box.dfs$method)[match("Unico",levels(addition.box.dfs$method))] <- "Unico"

addition.box.dfs

add.box.g = ggplot(addition.box.dfs, aes(x = source.id, y=y)) +
           geom_boxplot(aes(fill = method), position = position_dodge(0.65), width=2/k)+   
           geom_hline(yintercept=0, linetype='dotted', col = 'black', linewidth = 1) + 
           scale_fill_manual (values = method.colors) +
           scale_color_manual(values = method.colors) +

           theme_classic() + 
           ggtitle("") + 
           ylab(expression(Delta * - Log[10] * " " * P * " " * value))+                                    
           xlab(" ")+    
           scale_y_continuous(breaks = c(-100, -10, 0, 10, 100),
                              trans = scales::pseudo_log_trans(sigma = 0.25)) + 
 
           theme(plot.title = element_text(hjust = 0.5, size = title.size),
                 text=element_text(size=lab.size), 
                 axis.text.y = element_text(size=lab.size),  
                 axis.text.x = element_text(size=lab.size,
                                            #face = "bold", 
                                            vjust = 0.5, 
                                            #hjust = 0, 
                                            angle =  if (source.col == "decon.L1") 0 else 15),
                 legend.title = element_blank(),
                 legend.position="none") 
                
           

#figure.width  = (k) * 4
#figure.height = if (source.col == "decon.L1") 4 + 1 else 4 + 2
options(repr.plot.width = 16, repr.plot.height = 3, repr.plot.res = 200)

add.box.p  = egg::ggarrange(plots = list(add.box.g), 
                           #labels = c("c"),  
                           #label.args = list(gp = grid::gpar(font = 20, cex =3)), 
                           debug=F)

methods.full = c("bMIND", "bMIND (oracle prior)", " Unico")
legend.g = as_ggplot(get_legend(covar.bar.g + 
                                scale_fill_manual(labels = methods.full, values = method.colors) + 
                                theme(legend.direction = "horizontal",
                                      legend.position = c(0.575, 0.5),
                                      legend.text = element_text(size=legend.size, 
                                                                 margin = margin(r = 100, unit = "pt"))))) 


options(repr.plot.width = 16, repr.plot.height = 0.5, repr.plot.res = 100)
legend.g 

multi.g = ggarrange(legend.g,
                    params.bar.g, NULL,
                    cor.box.g, NULL,
                    add.box.p, 
                    nrow = 6,
                    heights = c(0.15, 1, 0.025, 2.25, 0.025, 1),
                    labels = c("", "a","", "b","", "c"),
                    font.label = list(size = 30, color = "black", face = "bold", family = NULL))

options(repr.plot.width = 16, repr.plot.height = 20, repr.plot.res = 200)
multi.g

ggsave(file.path(figure.dir, paste0("RNA_Simulation-", data.name, "_bMINDsc_multi.pdf")), multi.g,
       device = "pdf", width = 16, height = 20, dpi = 600)

fig.list = list(mean.bar.g = mean.bar.g, 
                var.bar.g = var.bar.g, 
                covar.bar.g = covar.bar.g, 
                Z.corrs.g = plts,
                add.box.g = add.box.g)

stats.list = list(covar.barplot.df = covar.barplot.df,
                  mean.barplot.df = mean.barplot.df,
                  var.barplot.df = var.barplot.df,
                  boxplot.meta = boxplot.meta, 
                  Z.corrs.list = Z.corrs.list,
                  addition.box.dfs = addition.box.dfs)


saveRDS(fig.list,   file.path(figure.dir, paste0("RNA_Simulation-", data.name, "_bMINDsc_fig.list.rds")))
saveRDS(stats.list, file.path(figure.dir, paste0("RNA_Simulation-", data.name, "_bMINDsc_stats.list.rds")))


