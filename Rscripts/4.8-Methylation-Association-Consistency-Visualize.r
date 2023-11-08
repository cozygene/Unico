library("ggplot2")
library("ggpubr")
library("scales")
library("reshape2")
library("ggpattern")

# mutation plot 
library("pracma")
library("ggplot2")
library("gghalves")
library("ggbeeswarm")

set.seed(2023)

metric = "MCC" # choose from "MCC" "F1"

res.dir    = "../Result/Methylation/Consistency/XY"
figure.dir = "../Figure/EWAS-Consistency/"

if (!file.exists(figure.dir)){dir.create(file.path(figure.dir),recursive = T)}
if (!file.exists(res.dir)){print("missing result")}

print(figure.dir)
print(res.dir)

age.meta.mats    = readRDS(file.path(res.dir, paste0("age.meta.mats.", metric, ".rds")))
gender.meta.mats = readRDS(file.path(res.dir, paste0("gender.meta.mats.", metric, ".rds")))

# add in bMIND for gender
gender.bMIND.meta.mats = readRDS(file.path(res.dir, "gender.bMIND.meta.mats.MCC.rds"))
gender.meta.mats = c(gender.meta.mats, gender.bMIND.meta.mats)
# for age just fill with NA
for (bmind.name in names(gender.bMIND.meta.mats)){   
    age.meta.mats[[bmind.name]] = matrix(NA, 4, 4)
    rownames(age.meta.mats[[bmind.name]]) = c("liu", "hannon1", "hannon2", "hannum")
    colnames(age.meta.mats[[bmind.name]]) = c("liu", "hannon1", "hannon2", "hannum")
}

gender.meta.mats

score1 = as.numeric(age.meta.mats$"Unico.parametric.marginal")
score1 = score1[score1!=1]

score0 = as.numeric(age.meta.mats$"TCA.parametric.marginal")
score0 = score0[score0!=1]

res = wilcox.test(x = score1, 
                  y = score0, 
                  alternative = "g", paired = T) # x > y one side test 

res

score1 = as.numeric(gender.meta.mats$"Unico.parametric.marginal")
score1 = score1[score1!=1]

score0 = as.numeric(gender.meta.mats$"TCA.parametric.marginal")
score0 = score0[score0!=1]

res = wilcox.test(x = score1, 
                  y = score0, 
                  alternative = "g", paired = T) # x > y one side test 

res

score1 = gender.meta.mats$"Unico.parametric.marginal" 
diag(score1) = NA
score1 = as.numeric(score1)
score1 = as.numeric(score1)[!is.na(as.numeric(score1))]

score0 = gender.meta.mats$"TCA.parametric.marginal" 
diag(score0) = NA
score0 = as.numeric(score0)
score0 = as.numeric(score0)[!is.na(as.numeric(score0))]

res = wilcox.test(x = score1, 
                  y = score0, 
                  alternative = "g", paired = T) # x > y one side test 

res

label.size = 20

F1_heatmap = function(res.meta.mat, phenotype, method.names, 
                      data.names, full.data.names,
                      legend.name = "F1 score", sub.title.name = "Mean F1", 
                      min.quantile = 0.1, max.quantile = 1, 
                      title.size = 15, label.size = 20, number.size = 5, 
                      barheight = 5, barwidth = 1, barres = 0.1, hide.x = F){


    f1_vals = c()
    mean_f1_vals = c()
    for (h in 1:length(res.meta.mat)){
        
        df = res.meta.mat[[h]][data.names, data.names]
        rownames(df) = full.data.names
        colnames(df) = full.data.names
        
        if(sum(!is.na(df)) == 0){
           
            df <- melt(df)
            colnames(df) = c("Discovery", "Validation", "Effective.F1")
            df$ignore = "Yes"
            res.meta.mat[[h]] = df
            f1_vals = c(f1_vals, df$"Effective.F1")
            mean_f1_vals[[h]] = NA
           
        }else{
           
            df <- melt(df)
            df$value = round(df$value, 2)

            colnames(df) = c("Discovery", "Validation", "Effective.F1")

            # to not show the diagnol, first make it NA, then the plotting wont show this
            df[df["Discovery"] == df["Validation"], "Effective.F1"] = NA

            # also make another column: "ignore" to explicitly indicate if we need to add stripes to show ignore 
            df$ignore = "No"
            df[df["Discovery"] == df["Validation"], "ignore"] = "Yes"

            res.meta.mat[[h]] = df
            f1_vals = c(f1_vals, df$"Effective.F1")
            mean_f1_vals[[h]] = mean(df$"Effective.F1", na.rm = T)
            
        }
    }


    f1_vals = f1_vals[!is.na(f1_vals)]
    min.f1 = quantile(f1_vals, min.quantile)
    max.f1 = quantile(f1_vals, max.quantile)
    
    plts <- vector(mode = "list", length = length(res.meta.mat))
    for (h in 1:length(method.names)){
        df = res.meta.mat[[h]]
        
        df$Effective.F1 = as.numeric(df$Effective.F1)
        # make a factor variable
        df$Discovery  <- as.factor(df$Discovery)
        df$Discovery  <- factor(df$Discovery , levels = levels(df$Discovery)[c(4,3,2,1)])

        g = ggplot(data = df, aes(x = Validation, y = Discovery, pattern = ignore, fill = Effective.F1)) +
            geom_tile(color = "black") +

            geom_tile_pattern(pattern_color = NA, pattern_fill = "black",
                              pattern_angle = 325, pattern_density = 0.25, 
                              pattern_spacing = 0.025, pattern_key_scale_factor = 1,
                              show.legend = FALSE) +
            scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +

            geom_text(aes(label = Effective.F1), 
                      color = "white", size = number.size) +

            scale_fill_gradientn(colors = hcl.colors(20, "cividis"), na.value="white",
                                 limits=c(min.f1, max.f1), oob=squish, 
                                 breaks = seq(round(min.f1, 1), round(max.f1,1), barres),
                                 name = legend.name) +
            coord_fixed() + 

            ggtitle(paste0(method.names[h], "-", phenotype, 
                           "\n (", sub.title.name, ": ", round(mean_f1_vals[[h]], 2), ")"))+
            
            theme(plot.title = element_text(hjust = 0.5, size = title.size),
                  text = element_text(size=label.size, face = "bold"),
                  axis.title.x = element_blank(),#remove x axis labels
                  axis.text.x  = element_text(angle = 25, vjust = 1, hjust = 1))+
            guides(fill = guide_colourbar(barwidth = barwidth, barheight = barheight)) 



        if (h != length(method.names)){g = g + theme(legend.position="none")}
        # dont show the y axis label names
        if (h != 1){
            g = g + theme(axis.title.y  = element_blank(),
                          axis.text.y  = element_blank(),  #remove y axis labels
                          axis.ticks.y = element_blank()   #remove y axis ticks
            )}
        if(hide.x){
            g = g + theme(axis.text.x  = element_blank(), 
                          axis.ticks.x = element_blank()  #remove x axis ticks
        )}

        plts[[h]] = g +  labs(tag = " ")
    }
    return(plts)
}



############################# Aggranging them with just ggarrange(not egg::ggarrange)  ############################
consistency.gender.marg.plts = F1_heatmap(res.meta.mat = gender.meta.mats[c("bMIND.YX.marginal",
                                                                            "Baseline.parametric.marginal", 
                                                                            "CellDMC.parametric.marginal", 
                                                                            "TCA.parametric.marginal", 
                                                                            "Unico.parametric.marginal")],
                                          phenotype = "Sex",
                                          method.names = c("bMIND","Baseline", "CellDMC", "TCA", "Unico"),
                                          data.names = c("liu", "hannon1", "hannon2", "hannum"),
                                          full.data.names = c("Liu et al. (n= 687)",
                                                              "Hannon et al. I (n= 675)", 
                                                              "Hannon et al. II (n= 665)",
                                                              "Hannum et al. (n= 590)"),

                                          legend.name = metric, sub.title.name = paste0("Mean ", metric), 
                                          min.quantile = 0.4, max.quantile = 1, 
                                          title.size = 15, label.size = label.size, number.size = 5, 
                                          barheight = 5, barwidth = 1, 
                                          barres = 0.1, hide.x = T)

consistency.age.marg.plts = F1_heatmap(res.meta.mat = age.meta.mats[c("bMIND.YX.marginal",
                                                                      "Baseline.parametric.marginal", 
                                                                      "CellDMC.parametric.marginal", 
                                                                      "TCA.parametric.marginal", 
                                                                      "Unico.parametric.marginal")],
                                       phenotype = "Age",
                                       method.names = c("bMIND", "Baseline", "CellDMC", "TCA", "Unico"),
                                       data.names = c("liu", "hannon1", "hannon2", "hannum"),
                                       full.data.names = c("Liu et al. (n= 687)",
                                                           "Hannon et al. I (n= 675)", 
                                                           "Hannon et al. II (n= 665)",
                                                           "Hannum et al. (n= 590)"),

                                       legend.name = metric, sub.title.name = paste0("Mean ", metric), 
                                       min.quantile = 0.1, max.quantile = 0.95, 
                                       title.size = 15, label.size = label.size, number.size = 5, 
                                       barheight = 5, barwidth = 1, 
                                       barres = 0.1, hide.x = F)


consistency.gender.marg.g = egg::ggarrange(plots=consistency.gender.marg.plts,
                                           align = "h",
                                           widths = rep(1, length(consistency.gender.marg.plts)),  
                                           debug=F)

consistency.age.marg.g = egg::ggarrange(plots= consistency.age.marg.plts,
                                        align = "h",
                                        widths = rep(1, length(consistency.age.marg.plts)),  
                                        debug=F)


marg.g  = ggarrange(consistency.gender.marg.g,
                    consistency.age.marg.g,
                    nrow = 2, 
                    heights = c(0.90,  1.275),
                    widths =  c(1,  1))
                    #font.label = list(size = 25, color = "black", face = "bold", family = NULL)) 

marg.g =  marg.g + annotate("text", x = 0.555, y = -0.01, label = 'bold("Validation")', 
                            parse = TRUE, size = 7.5) +
                   theme(plot.margin = unit(c(0,0,3,0), "lines"))


options(repr.plot.width = 16, repr.plot.height = 8, repr.plot.res = 200)
marg.g 

ggsave(file.path(figure.dir, paste0("Consistency.marg.", metric, ".pdf")), marg.g,
       bg = 'white',
       device = "pdf", width = 16, height = 8, dpi = 600)

############################ Aggranging them with just ggarrange(not egg::ggarrange)  ############################
consistency.gender.joint.plts = F1_heatmap(res.meta.mat = gender.meta.mats[c("bMIND.XY.joint",
                                                                             "Baseline.parametric.joint", 
                                                                             "CellDMC.parametric.joint", 
                                                                             "TCA.parametric.joint", 
                                                                             "Unico.parametric.joint")],
                                           phenotype = "Sex",
                                           method.names = c("bMIND", "Bulk", "CellDMC", "TCA", "Unico"),
                                           data.names = c("liu", "hannon1", "hannon2", "hannum"),
                                           full.data.names = c("Liu et al. (n= 687)",
                                                               "Hannon et al. I (n= 675)", 
                                                               "Hannon et al. II (n= 665)",
                                                               "Hannum et al. (n= 590)"),
                                           

                                           legend.name = metric, sub.title.name = paste0("Mean ", metric), 
                                           min.quantile = 0.15, max.quantile = 0.9, 
                                           title.size = 15, label.size = label.size, number.size = 5, 
                                           barheight = 5, barwidth = 1, 
                                           barres = 0.15, hide.x = T)

consistency.age.joint.plts = F1_heatmap(res.meta.mat = age.meta.mats[c("bMIND.XY.joint",
                                                                       "Baseline.parametric.joint", 
                                                                       "CellDMC.parametric.joint", 
                                                                       "TCA.parametric.joint", 
                                                                       "Unico.parametric.joint")],
                                        phenotype = "Age",
                                        method.names = c("bMIND", "Bulk", "CellDMC", "TCA", "Unico"),
                                        data.names = c("liu", "hannon1", "hannon2", "hannum"),
                                        full.data.names = c("Liu et al. (n= 687)",
                                                            "Hannon et al. I (n= 675)", 
                                                            "Hannon et al. II (n= 665)",
                                                            "Hannum et al. (n= 590)"),

                                        legend.name = metric, sub.title.name = paste0("Mean ", metric), 
                                        min.quantile = 0.05, max.quantile = 0.9, 
                                        title.size = 15, label.size = label.size, number.size = 5, 
                                        barheight = 5, barwidth = 1, 
                                        barres = 0.15, hide.x = F)

consistency.gender.joint.g = egg::ggarrange(plots=consistency.gender.joint.plts,
                                            align = "h",
                                            widths = rep(1, length(consistency.gender.joint.plts)),  
                                            debug=F)

consistency.age.joint.g = egg::ggarrange(plots=consistency.age.joint.plts,
                                         align = "h",
                                         widths = rep(1, length(consistency.age.joint.plts)),  
                                         debug=F)

joint.g  = ggarrange(consistency.gender.joint.g,
                     consistency.age.joint.g,
                     nrow = 2, 
                     heights = c(0.90,  1.275),
                     widths =  c(1,  1))         

joint.g = joint.g +  annotate("text", x = 0.555, y = -0.01, label = 'bold("Validation")', 
                              parse = TRUE, size = 7.5) +
                    theme(plot.margin = unit(c(0,0,3,0), "lines"))

options(repr.plot.width = 16, repr.plot.height = 8, repr.plot.res = 200)
joint.g

ggsave(file.path(figure.dir, paste0("Consistency.joint.", metric, ".pdf")), joint.g,
       bg = 'white',
       device = "pdf", width = 16, height = 8, dpi = 600)

fig.list = list(consistency.gender.marg.plts  = consistency.gender.marg.plts, 
                consistency.age.marg.plts     = consistency.age.marg.plts, 
                consistency.gender.joint.plts = consistency.gender.joint.plts, 
                consistency.age.joint.plts    = consistency.age.joint.plts)

saveRDS(fig.list,   file.path(figure.dir, paste0("Consistency.", metric, ".fig.list.rds")))




