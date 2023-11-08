library("data.table")

library("ggplot2")
library("ggpubr")
library("dplyr")

set.seed(2023)

res.dir    = "../Result/Runtime-local/"
figure.dir = "../Figure/Runtime/"

if (!file.exists(res.dir)){print("no result yet")}
if (!file.exists(figure.dir)){dir.create(file.path(figure.dir),recursive = T)}



runtime.list = readRDS(file.path(res.dir, "runtime.list.rds"))

tensor.k.df = runtime.list$tensor.k.df
asso.k.df   = runtime.list$asso.k.df
tensor.n.df = runtime.list$tensor.n.df
asso.n.df   = runtime.list$asso.n.df


tensor.k.df$method = factor(tensor.k.df$method, levels = c("CIBERSORTx","TCA", "bMIND", "Unico"))
levels(tensor.k.df$method)[match("Unico",levels(tensor.k.df$method))] <- "Unico"
tensor.k.summary.df <- tensor.k.df %>% group_by(k, method) %>% dplyr::summarise(sd = sd(time), time = mean(time))


asso.k.df $method  = factor(asso.k.df $method, levels = c("CellDMC", "TCA","bMIND", "Unico-parametric", "Unico-asymptotic"))
levels(asso.k.df $method)[match("Unico-parametric",levels(asso.k.df$method))] <- "Unico (Parametric)"
levels(asso.k.df $method)[match("Unico-asymptotic",levels(asso.k.df$method))] <- "Unico (Asymptotic)"
asso.k.summary.df <- asso.k.df %>% group_by(k, method) %>% dplyr::summarise(sd = sd(time), time = mean(time))



tensor.n.df$method = factor(tensor.n.df$method, levels = c("CIBERSORTx","TCA", "bMIND", "Unico"))
levels(tensor.n.df$method)[match("Unico",levels(tensor.n.df$method))] <- "Unico"
tensor.n.summary.df <- tensor.n.df %>% group_by(n, method) %>% dplyr::summarise(sd = sd(time), time = mean(time))

asso.n.df $method  = factor(asso.n.df $method, levels = c("CellDMC", "TCA","bMIND", "Unico-parametric", "Unico-asymptotic"))
levels(asso.n.df $method)[match("Unico-parametric",levels(asso.n.df$method))] <- "Unico-parametric"
levels(asso.n.df $method)[match("Unico-asymptotic",levels(asso.n.df$method))] <- "Unico-asymptotic"
asso.n.summary.df   <- asso.n.df %>% group_by(n, method) %>% dplyr::summarise(sd = sd(time), time = mean(time))


tensor.k.df[(tensor.k.df$k == 4) & (tensor.k.df$method == "CIBERSORTx"), ]

tensor.k.df[(tensor.k.df$k == 6) & (tensor.k.df$method == "CIBERSORTx"), ]



tensor.k.summary.df

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

base.color       = brewer.pal(n = 8, name = "Set2")[8]
cibersortx.color = brewer.pal(n = 8, name = "Set2")[5]
CellDMC.color    = brewer.pal(n = 12, name = "Paired")[4]
Unico.color       = brewer.pal(n = 8, name = "Set2")[4]
Unico.asymp.color = brewer.pal(n = 12, name = "Paired")[6]

tca.color        = brewer.pal(n = 8, name = "Set2")[6]
bMIND_sc.color   = brewer.pal(n = 11, name = "RdYlBu")[10]
bMIND_rp.color   = brewer.pal(n = 8, name = "Set2")[3]

method.tensor.colors = c(cibersortx.color, tca.color, bMIND_rp.color, Unico.color)
method.asso.colors   = c(CellDMC.color, tca.color, bMIND_rp.color, Unico.color, Unico.asymp.color)

title.size = 20
lab.size = 15

options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 100)

####################################### plotting tensor (k is changing) ##########################################
pk = ggplot(tensor.k.df, aes(k, time, color = method)) +
     
     #geom_jitter(position = position_jitter(0.2)) + 
     geom_line(aes(group = method), data = tensor.k.summary.df) +
     geom_errorbar(data = tensor.k.summary.df, 
                   aes(ymin = time-sd, ymax = time+sd),  width = 0.2) +
     scale_colour_manual(values = method.tensor.colors, name=NULL) +
     scale_fill_manual(values = method.tensor.colors, name=NULL) + 

     scale_x_continuous(breaks = c(3, 4, 5, 6),
                        labels=c("k = 3", "k = 4","k = 5", "k = 6")) + 
     scale_y_continuous(breaks = c(50, 100, 250, 1000),
                        trans='log10') +
     #ylab(expression(atop("Tensor Deconvolution", "Seconds (log scale)"))) + 
     ylab(expression(atop("Deconvolution runtime", "(seconds)"))) + 
   

     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     # label size
     theme(axis.title.x = element_blank()) + 
     theme(axis.title.y = element_text(size = title.size)) + 
     #ticks size
     theme(axis.text.y = element_text(size = lab.size + 5 , face="bold")) +
     theme(axis.text.x = element_text(size = lab.size + 5, face="bold", 
                                      #angle = 20, vjust = 0.5
                                     ))+
     theme(legend.text = element_text(size = lab.size + 5)) + 
     theme(legend.position = "none")

pk


####################################### plotting tensor (n is changing) ##########################################
pn = ggplot(tensor.n.df, aes(n, time, color = method)) +
     #geom_jitter(position = position_jitter(0.04)) + 
     geom_line(aes(group = method), data = tensor.n.summary.df) +
     geom_errorbar(data = tensor.n.summary.df, aes(ymin = time-sd, ymax = time+sd),  width = 0.04) +
     scale_colour_manual(values = method.tensor.colors, name=NULL) +
     scale_fill_manual(values = method.tensor.colors, name=NULL) + 

     scale_x_continuous(breaks = c(100, 250, 500),
                                         labels=c("n = 100", "n = 250","n = 500"),
                                         trans='log10') + 
     scale_y_continuous(breaks = c(50,  250, 1000),
                        trans='log10') +
     ylab("") + 

     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     # label size
     theme(axis.title.x = element_blank()) + 
     theme(axis.text.y = element_blank())+
     #ticks size
     theme(axis.text.y = element_text(size = lab.size + 5 , face="bold")) +
     theme(axis.text.x = element_text(size = lab.size + 5, face="bold", 
                                      #angle = 20, vjust = 0.5
                                     ))+
     theme(legend.text = element_text(size = lab.size + 5)) 


pn

ak = ggplot(asso.k.df, aes(k, time, color = method)) +

     #geom_jitter(position = position_jitter(0.2)) + 
     geom_line(aes(group = method), data = asso.k.summary.df) +
     geom_errorbar(data = asso.k.summary.df, aes(ymin = time-sd, ymax = time+sd),  width = 0.2) +
     scale_colour_manual(values = method.asso.colors, name=NULL) +
     scale_fill_manual(values = method.asso.colors, name=NULL) + 

     scale_x_continuous(breaks = c(3, 4, 5, 6),
                        labels=c("k = 3", "k = 4","k = 5", "k = 6")) + 
     scale_y_continuous(breaks = c(0.5, 2.5, 10),
                        trans='log10') +
     #ylab(expression(atop("Cell Type Specific Association", "Seconds (log scale)"))) + 
     ylab(expression(atop("Cell-type association testing", "(seconds)"))) + 
     

     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     # label size
     theme(axis.title.x = element_blank()) + 
     theme(axis.title.y = element_text(size = title.size)) + 
     #ticks size
     theme(axis.text.y = element_text(size = lab.size + 5 , face="bold")) +
     theme(axis.text.x = element_text(size = lab.size + 5, face="bold", 
                                      #angle = 20, vjust = 0.5
                                     ))+
     theme(legend.text = element_text(size = lab.size + 5)) + 
     theme(legend.position = "none")



ak

an = ggplot(asso.n.df, aes(x = n, y = time, color = method)) +
     #ggplot(asso.n.df, aes(x = factor(n), y = time)) +
     #geom_violin(aes(fill =  method), position = position_dodge(0.2), 
     #            scale="count", width=5) + 
     
     #geom_jitter(position = position_jitter(0.04)) + 
     geom_line(aes(group = method), data = asso.n.summary.df) +
     geom_errorbar(data = asso.n.summary.df, aes(ymin = time-sd, ymax = time+sd),  width = 0.04) +
     scale_colour_manual(values = method.asso.colors, name=NULL) +
     scale_fill_manual(values = method.asso.colors, name=NULL) + 

     scale_x_continuous(breaks = c(100, 250, 500),
                        labels=c("n = 100", "n = 250","n = 500"), #) + 
                        trans='log10') + 
     scale_y_continuous(breaks = c(0.5, 2.5, 10),
                        trans='log10') +
     ylab("") + 

     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     # label size
     theme(axis.title.x = element_blank()) + 
     theme(axis.text.y = element_blank())+
     #ticks size
     theme(axis.text.y = element_text(size = lab.size + 5, face="bold")) +
     theme(axis.text.x = element_text(size = lab.size + 5, face="bold", 
                                      #angle = 20, vjust = 0.5
                                     ))+
     theme(legend.text = element_text(size = lab.size + 5)) 

an

options(repr.plot.width = 15, repr.plot.height = 10, repr.plot.res = 100)
multi.g = egg::ggarrange(plots=list(pk + theme(legend.position = "none"), 
                                    pn + theme(legend.position = "right"),
                                    ak + theme(legend.position = "none"),
                                    an + theme(legend.position = "right")),  
               labels = c("a", "", "b", ""), 
               align = "h",
               widths = c(3,3.5),  
               heights = c(4,4),
               label.args = list(gp = grid::gpar(font = 20, cex =3)), 
               debug=F)

ggsave(file.path(figure.dir, paste0("Runtime_multi.pdf")), multi.g,
       device = "pdf", width = 15, height = 10, dpi = 600)

fig.list = list(tensor.k.g = pk, 
                tensor.n.g = pn, 
                asso.k.g = ak, 
                asso.n.g = an)

stats.list = list(tensor.k.df = tensor.k.df,
                  tensor.n.df = tensor.n.df,
                  asso.k.df   = asso.k.df,
                  asso.n.df   = asso.n.df)


saveRDS(fig.list,   file.path(figure.dir, paste0("Runtime_fig.list.rds")))
saveRDS(stats.list, file.path(figure.dir, paste0("Runtime_stats.list.rds")))


