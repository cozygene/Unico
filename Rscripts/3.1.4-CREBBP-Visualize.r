library("plyr")
library("pracma")
library("ggpubr")
library("ggplot2")
library("gghalves")
library("ggbeeswarm")
set.seed(2023)

data.dir    = file.path("../Data/RNA/CREBBP/")
res.dir     = file.path("../Result/RNA/CREBBP/")
figure.dir  = file.path("../Figure/CREBBP/")

if (!file.exists(figure.dir)){
    dir.create(figure.dir, recursive = T)
}

CREBBP.dat  = readRDS(file.path(data.dir, "CREBBP.dat.rds"))
CREBBP.eval = readRDS(file.path(res.dir, "CREBBP.eval.rds"))

inc.genes  = CREBBP.dat$inc.genes
dec.genes  = CREBBP.dat$dec.genes
samples.mt = CREBBP.dat$samples.mt
samples.wt = CREBBP.dat$samples.wt

print(paste0("length inc.genes: ", length(inc.genes)))
print(paste0("length dec.genes: ", length(dec.genes)))
print(paste0("length samples.mt: ", length(samples.mt)))
print(paste0("length samples.wt: ", length(samples.wt)))

key     = "genes"

W = CREBBP.dat$expm.dat[[key]]$W
X = CREBBP.dat$expm.dat[[key]]$X
source.ids  = colnames(W)
feature.ids = rownames(X)
sample.ids  = rownames(W)

k = length(source.ids)
m = length(feature.ids)
n = length(sample.ids)

print(source.ids)

########################################### Color code ################################################
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

method.colors = list("Baseline"   = brewer.pal(n = 8, name = "Set2")[8],
                     "Bulk"       = brewer.pal(n = 8, name = "Dark2")[7],
                     "CIBERSORTx" = brewer.pal(n = 8, name = "Set2")[5],
                     "TCA"        = brewer.pal(n = 8, name = "Set2")[6],
                     "bMIND"      = brewer.pal(n = 8, name = "Set2")[3],
                     "Unico"       = brewer.pal(n = 8, name = "Set2")[4],
                     "Unico"      = brewer.pal(n = 8, name = "Set2")[4])



mt.color = brewer.pal(n = 11, name = "RdYlBu")[1]
wt.color = brewer.pal(n = 11, name = "RdBu")[6]

scientific_pvals_10x <- function(values, digits = 1) {
    if(!is.numeric(values)){
        stop("values must be numbers")
    }
    if(grepl("^\\d{2}$", digits)){
        stop("digits must a one or two digit whole number")
    }

    x <- sprintf(paste0("%.", digits, "e"), values)
    x <- gsub("^(.*)e", "'\\1'e", x)
    longestExponent <- max(sapply(gregexpr("\\d{1,}$", x), attr, 'match.length'))
    zeroTrimmed <- ifelse(longestExponent > 2,
                          paste0("\\1", paste(rep("~", times = longestExponent-1), collapse = "")),
                          "\\1")
    x <- gsub("(e[+|-])[0]", zeroTrimmed, x)

    x <- gsub("e", "~x~10^", x)

    if(any(grepl("\\^\\-", x))){
        x <- gsub("\\^\\+", "\\^~~", x)
    } else {
        x <- gsub("\\^\\+", "\\^", x)
    }
    # return this as an expression
    #parse(text="P = " ~ x)
    return(x)
} 

methods = c("Baseline", "CIBERSORTx", "TCA", "Bulk", "bMIND", "Unico")
logDiff.df = CREBBP.eval$tensor$logDiff.df

colnames(logDiff.df)[which(colnames(logDiff.df) == "bMIND.scale")] = "bMIND"

logDiff.df

de.plot.df = data.frame()
for (method in methods){
    method.df = data.frame(diff     = logDiff.df[c(inc.genes, dec.genes), method],
                           method   = rep(method, m),
                           DE.type  = c(rep("Up Regulated", length(inc.genes)), 
                                        rep("Down Regulated", length(dec.genes))))
    de.plot.df = rbind(de.plot.df, method.df)
}
de.plot.df$method = factor(de.plot.df$method, levels = methods)
levels(de.plot.df$method)[match("Unico",levels(de.plot.df$method))] <- "Unico"

de.plot.df

pvals.inc = CREBBP.eval$tensor$pvals.inc
pvals.dec = CREBBP.eval$tensor$pvals.dec

colnames(pvals.inc)[which(colnames(pvals.inc) == "bMIND.scale")] = "bMIND"
colnames(pvals.dec)[which(colnames(pvals.dec) == "bMIND.scale")] = "bMIND"

pvals.inc

title.size = 20
lab.size = 18


tp = ggplot(de.plot.df, aes(x=DE.type, y=diff, fill = method)) + 
     geom_beeswarm(aes(DE.type, y = diff, fill = method, colour = method),
                   priority = 'density',
                   side = -1, dodge.width = 0.8, cex = 0.15, size = 1.75, colour = "black", #dot border 
                   pch = 21, shape = 22, stroke = 0, show.legend = F) +

     geom_half_boxplot(aes(fill = method), side = 'r', color = 'black', 
                       outlier.size = NA, outlier.color = NA, notchwidth = 0.5, size = 0.7, #thickness of the line 
                       width = 0.6, position = position_dodge(0.8), alpha = 1) +  
     geom_hline(yintercept=0, linetype='dotted', col = 'grey', size=1) + 
    
    
     scale_fill_manual(values = as.character(method.colors[levels(de.plot.df$method)]), name=NULL) +  
     coord_cartesian(ylim = c(-65,2)) + # for annotation to have space
     scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                        trans = scales::pseudo_log_trans(sigma = 0.25)) + 
     scale_x_discrete(labels=c("Down Regulated" = "CREBBP down-regulated genes \n (lower is better)", 
                               "Up Regulated"   = "CREBBP up-regulated genes \n (higher is better)")) +
 
     #ylab(expression(atop("Differences in "*Log[2]*" Expression", "(Deconvolution B cells)"))) + 
     ylab("Log fold change") + 

     #ggtitle("CREBBP Mutation")  + 
     ggtitle("")  + 

     theme_classic() + # later operation will overide this
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     # label size
     theme(axis.title.x = element_blank()) + 
     theme(axis.title.y = element_text(size = lab.size + 2, face="bold")) + 
     #ticks size
     theme(axis.text.y = element_text(size = lab.size, face="bold")) +
     theme(axis.text.x = element_text(size = lab.size, face="bold"))+
     theme(legend.text = element_text(size = lab.size, face="bold"))

options(repr.plot.width = 16, repr.plot.height = 6, repr.plot.res = 150)
tp

h.loc  = linspace(0.665, 1.33, 6)
v.loc  = -exp(linspace(0.9, 3.75, 5))
v.size = -(exp(linspace(0.4, 1.5, 5))-1)
p.size = 4.75

df1 <- data.frame(a = c(h.loc[1],           h.loc[1], h.loc[6], h.loc[6]), 
                  b = c(v.loc[1]-v.size[1], v.loc[1], v.loc[1], v.loc[1]-v.size[1]))

df2 <- data.frame(a = c(h.loc[2],           h.loc[2], h.loc[6], h.loc[6]), 
                  b = c(v.loc[2]-v.size[2], v.loc[2], v.loc[2], v.loc[2]-v.size[2]))

df3 <- data.frame(a = c(h.loc[3],           h.loc[3], h.loc[6], h.loc[6]), 
                  b = c(v.loc[3]-v.size[3], v.loc[3], v.loc[3], v.loc[3]-v.size[3]))

df4 <- data.frame(a = c(h.loc[4],           h.loc[4], h.loc[6], h.loc[6]), 
                  b = c(v.loc[4]-1.5 * v.size[4], v.loc[4], v.loc[4], v.loc[4]-1.5 * v.size[4]))

df5 <- data.frame(a = c(h.loc[5],           h.loc[5], h.loc[6], h.loc[6]), 
                  b = c(v.loc[5]-2 * v.size[5], v.loc[5], v.loc[5], v.loc[5]-2 * v.size[5]))


df6  <- data.frame(a = df1$a + 1, b = df1$b)
df7  <- data.frame(a = df2$a + 1, b = df2$b)
df8  <- data.frame(a = df3$a + 1, b = df3$b)
df9  <- data.frame(a = df4$a + 1, b = df4$b)
df10 <- data.frame(a = df5$a + 1, b = df5$b)


tp = tp + geom_line(data = df1, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df2, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df3, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df4, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df5, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df6, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df7, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df8, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df9, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df10,aes(x = a, y = b), inherit.aes = FALSE) +

          annotate("text", x = (h.loc[1]+ h.loc[6])/2,       y = v.loc[1] + 1.5*v.size[1],      label = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "Baseline"],   digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = (h.loc[2]+ h.loc[6])/2,       y = v.loc[2] + 1.5*v.size[2],      label = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "CIBERSORTx"], digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = (h.loc[3]+ h.loc[6])/2,       y = v.loc[3] + 2.0*v.size[3],      label = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "TCA"],        digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = (h.loc[4]+ h.loc[6])/2,       y = v.loc[4] + 2.5*v.size[4],      label = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "Bulk"],       digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = (h.loc[5]+ h.loc[6])/2,       y = v.loc[5] + 2.5*v.size[5],      label = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "bMIND"],      digits = 1)), parse=TRUE, size = p.size) +

          annotate("text", x = 1 + (h.loc[1]+ h.loc[6])/2,   y = (v.loc[1] + 1.5*v.size[1]),    label = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "Baseline"],   digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = 1 + (h.loc[2]+ h.loc[6])/2,   y = (v.loc[2] + 1.5*v.size[2]),    label = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "CIBERSORTx"], digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = 1 + (h.loc[3]+ h.loc[6])/2,   y = (v.loc[3] + 2.0*v.size[3]),    label = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "TCA"],        digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = 1 + (h.loc[4]+ h.loc[6])/2,   y = (v.loc[4] + 2.5*v.size[4]),    label = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "Bulk"],       digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = 1 + (h.loc[5]+ h.loc[6])/2,   y = (v.loc[5] + 2.5*v.size[5]),    label = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "bMIND"],      digits = 1)), parse=TRUE, size = p.size) #+

options(repr.plot.width = 16, repr.plot.height = 6, repr.plot.res = 150)
tp

############################ figure boxplot ############################
W.df  = data.frame(W        = c(W[samples.mt, "Bcells"], W[samples.wt, "Bcells"]),
                   mutation = c(rep("MT", length(samples.mt)), rep("WT", length(samples.wt))))
W.df$mutation = factor(W.df$mutation, levels = c("WT","MT"))


wp = ggplot(W.df, aes(x=mutation, y=W, fill = mutation)) + #geom_boxplot()
     geom_dotplot(binaxis = "y", stackdir = "down", dotsize = 2.5) + 
     geom_half_boxplot(aes(fill = mutation), side = 'r', color = 'black', 
                       outlier.size = NA, outlier.color = NA, notchwidth = 0.5, size = 1, #thickness of the line 
                       width = 0.8, position = position_dodge(0.9), alpha = 1) + 
     #scale_y_continuous(position = "right") + 
     scale_fill_manual(values = c(wt.color, mt.color), name=NULL) + 
     coord_cartesian(ylim = c(0,0.8))  + 
     ylab(expression(atop("CIBERSORTx Estimated","B cells Proportion"))) + 


     theme_classic() + # later operation will overide this
     theme(legend.position = "none") + 
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     theme(axis.title.x = element_blank()) + 
     theme(axis.title.y = element_text(size = title.size)) + 
     #ticks size
     theme(axis.text.y = element_text(size = lab.size, face="bold")) +
     theme(axis.text.x = element_text(size = lab.size + 5, face="bold"))

options(repr.plot.width = 3, repr.plot.height = 6, repr.plot.res = 150)
wp

str(CREBBP.eval$regression)

method.beta.list = CREBBP.eval$regression$ method.beta.list 
reg.pvals.inc    = CREBBP.eval$regression$ reg.pvals.inc
reg.pvals.dec    = CREBBP.eval$regression$ reg.pvals.dec

names(method.beta.list)[which(names(method.beta.list) == "bMIND.scale")] = "bMIND"
colnames(reg.pvals.inc)[which(colnames(reg.pvals.inc) == "bMIND.scale")] = "bMIND"
rownames(reg.pvals.inc)[which(rownames(reg.pvals.inc) == "bMIND.scale")] = "bMIND"
colnames(reg.pvals.dec)[which(colnames(reg.pvals.dec) == "bMIND.scale")] = "bMIND"
rownames(reg.pvals.dec)[which(rownames(reg.pvals.dec) == "bMIND.scale")] = "bMIND"

method.beta.list

reg.plot.df = data.frame()
for (method in methods){
    method.df = data.frame(beta     = method.beta.list[[method]][c(inc.genes, dec.genes), "mutation.status"],
                           method   = rep(method, m),
                           DE.type  = c(rep("Up Regulated", length(inc.genes)), 
                                        rep("Down Regulated", length(dec.genes))))
    reg.plot.df = rbind(reg.plot.df, method.df)
}
reg.plot.df$method = factor(reg.plot.df$method, levels = c('Baseline','CIBERSORTx','TCA', 'Bulk', 'bMIND', 'Unico'))
levels(reg.plot.df$method)[match("Unico",levels(reg.plot.df$method))] <- "Unico"

reg.plot.df

rp = ggplot(reg.plot.df, aes(x=DE.type, y=beta, fill = method)) + 
     geom_beeswarm(aes(DE.type, y = beta, fill = method, colour = method),
                   priority = 'density', cex = 0.3, #point space
                   side = -1, dodge.width = 0.8, size = 1.5, colour = "black", #dot border 
                   pch = 21, shape = 22, stroke = 0, show.legend = F) +

     geom_half_boxplot(aes(fill = method), side = 'r', color = 'black', 
                       outlier.size = NA, outlier.color = NA, notchwidth = 0.5, size = 0.7, #thickness of the line 
                       width = 0.6, position = position_dodge(0.8), alpha = 1) +  
     geom_hline(yintercept=0, linetype='dotted', col = "black", size=1) + 
    
    
     scale_fill_manual(values = as.character(method.colors[levels(reg.plot.df$method)]), name=NULL) +  
     coord_cartesian(ylim = c(-60,1.25)) + # for annotation to have space
     #scale_y_continuous(trans = "exp") + 
     scale_y_continuous(breaks = c(-1, -0.25, 0, 0.25, 1),
                        trans = scales::pseudo_log_trans(sigma = 0.025)) + 
     scale_x_discrete(labels=c("Down Regulated" = "CREBBP down-regulated genes \n (lower is better)", 
                               "Up Regulated"   = "CREBBP up-regulated genes \n (higher is better)")) +

     ylab("Estimated effect size") + 
     ggtitle("")  + 
     #ggtitle("CREBBP Mutation")  + 

     theme_classic() + # later operation will overide this
     theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
     # label size
     theme(axis.title.x = element_blank()) + 
     theme(axis.title.y = element_text(size = lab.size + 3, face="bold")) + 
     #ticks size
     theme(axis.text.y = element_text(size = lab.size, face="bold")) +
     theme(axis.text.x = element_text(size = lab.size, face="bold"))+
     theme(legend.text = element_text(size = lab.size, face="bold"))



options(repr.plot.width = 16, repr.plot.height = 6, repr.plot.res = 150)
rp

h.loc  = linspace(0.665, 1.33, 6)
v.loc  = -exp(linspace(0.8, 4.00, 5))
v.size = -(exp(linspace(0.45, 1.8, 5))-1)
p.size = 4.5

df1 <- data.frame(a = c(h.loc[1],           h.loc[1], h.loc[6], h.loc[6]), 
                  b = c(v.loc[1]-v.size[1], v.loc[1], v.loc[1], v.loc[1]-v.size[1]))

df2 <- data.frame(a = c(h.loc[2],           h.loc[2], h.loc[6], h.loc[6]), 
                  b = c(v.loc[2]-v.size[2], v.loc[2], v.loc[2], v.loc[2]-v.size[2]))

df3 <- data.frame(a = c(h.loc[3],           h.loc[3], h.loc[6], h.loc[6]), 
                  b = c(v.loc[3]-v.size[3], v.loc[3], v.loc[3], v.loc[3]-v.size[3]))

df4 <- data.frame(a = c(h.loc[4],           h.loc[4], h.loc[6], h.loc[6]), 
                  b = c(v.loc[4]-1.5 * v.size[4], v.loc[4], v.loc[4], v.loc[4]-1.5 * v.size[4]))

df5 <- data.frame(a = c(h.loc[5],           h.loc[5], h.loc[6], h.loc[6]), 
                  b = c(v.loc[5]-2 * v.size[5], v.loc[5], v.loc[5], v.loc[5]-2 * v.size[5]))


df6  <- data.frame(a = df1$a + 1, b = df1$b)
df7  <- data.frame(a = df2$a + 1, b = df2$b)
df8  <- data.frame(a = df3$a + 1, b = df3$b)
df9  <- data.frame(a = df4$a + 1, b = df4$b)
df10 <- data.frame(a = df5$a + 1, b = df5$b)


rp = rp + geom_line(data = df1, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df2, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df3, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df4, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df5, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df6, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df7, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df8, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df9, aes(x = a, y = b), inherit.aes = FALSE) +
          geom_line(data = df10,aes(x = a, y = b), inherit.aes = FALSE) +

          annotate("text", x = (h.loc[1]+ h.loc[6])/2,       y = v.loc[1] + 1.5*v.size[1],      label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.dec["Unico", "Baseline"],   digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = (h.loc[2]+ h.loc[6])/2,       y = v.loc[2] + 1.5*v.size[2],      label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.dec["Unico", "CIBERSORTx"], digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = (h.loc[3]+ h.loc[6])/2,       y = v.loc[3] + 2.0*v.size[3],      label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.dec["Unico", "TCA"],        digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = (h.loc[4]+ h.loc[6])/2,       y = v.loc[4] + 2.5*v.size[4],      label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.dec["Unico", "Bulk"],       digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = (h.loc[5]+ h.loc[6])/2,       y = v.loc[5] + 2.5*v.size[5],      label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.dec["Unico", "bMIND"],      digits = 1)), parse=TRUE, size = p.size) +

          annotate("text", x = 1 + (h.loc[1]+ h.loc[6])/2,   y = (v.loc[1] + 1.5*v.size[1]),    label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.inc["Unico", "Baseline"],   digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = 1 + (h.loc[2]+ h.loc[6])/2,   y = (v.loc[2] + 1.5*v.size[2]),    label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.inc["Unico", "CIBERSORTx"], digits = 1)), parse=TRUE, size = p.size) + 
          annotate("text", x = 1 + (h.loc[3]+ h.loc[6])/2,   y = (v.loc[3] + 2.0*v.size[3]),    label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.inc["Unico", "TCA"],        digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = 1 + (h.loc[4]+ h.loc[6])/2,   y = (v.loc[4] + 2.5*v.size[4]),    label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.inc["Unico", "Bulk"],       digits = 1)), parse=TRUE, size = p.size) +
          annotate("text", x = 1 + (h.loc[5]+ h.loc[6])/2,   y = (v.loc[5] + 2.5*v.size[5]),    label = paste0("\'P = \'~", scientific_pvals_10x(reg.pvals.inc["Unico", "bMIND"],      digits = 1)), parse=TRUE, size = p.size) #+

reg.pvals.dec

options(repr.plot.width = 16, repr.plot.height = 6, repr.plot.res = 150)
rp







CREBBP.supp.g = egg::ggarrange(plots=list(wp,tp),  
                               labels = c("a", "b"), 
                               align = "h",
                               widths = c(2.5,9),  
                               label.args = list(gp = grid::gpar(font = 30, cex =3)), 
                               debug=F)

# #doesnt align that well the y axis
# CREBBP.supp.g = ggarrange(wp,tp,
#                           nrow = 1, 
#                           widths = c(2.5,9), 
#                           labels = c("a", "b"), 
#                           font.label = list(size = 30, color = "black", face = "bold", family = NULL))

options(repr.plot.width = 16, repr.plot.height = 7, repr.plot.res = 150)
CREBBP.supp.g

ggsave(file.path(figure.dir, "CREBBP.supp.multi.pdf"), CREBBP.supp.g,
       bg = 'white',
       device = "pdf", width = 16, height = 7, dpi = 600)

fig.list = list(w.g = wp, 
                de.g = tp, 
                regress.g = rp)

stats.list = list(w.df = W.df,
                  de.plot.df = de.plot.df, 
                  reg.plot.df = reg.plot.df)


saveRDS(fig.list,   file.path(figure.dir, paste0("CREBBP_fig.list.rds")))
saveRDS(stats.list, file.path(figure.dir, paste0("CREBBP_stats.list.rds")))



# # DE Main figure

# main.methods.plot = c("CIBERSORTx","TCA", "Unico")
# main.plot.df = de.plot.df[de.plot.df$method %in% main.methods.plot, ]
# main.plot.df$method = factor(main.plot.df$method, levels = main.methods.plot)

# main.plot.df

# title.size = 20
# lab.size = 15

# mp = ggplot(main.plot.df, aes(x=DE.type, y=diff, fill = method)) + 
#      geom_beeswarm(aes(DE.type, y = diff, fill = method, colour = method),
#                    priority = 'density',
#                    side = -1, dodge.width = 0.99, cex = 0.25,  size = 1.75, colour = "black", #dot border 
#                    pch = 21, shape = 22, stroke = 0, show.legend = F) +

#      geom_half_boxplot(aes(fill = method), side = 'r', color = 'black', 
#                        outlier.size = NA, outlier.color = NA, notchwidth = 0.5, size = 0.7, #thickness of the line 
#                        width = 0.6, position = position_dodge(0.99), alpha = 1) +  
#      geom_hline(yintercept=0, linetype='dotted', col = 'grey', size=1) + 
    
    
#      scale_fill_manual(values = method.colors[main.methods.plot], name=NULL) +  
#      coord_cartesian(ylim = c(-8,8)) + # for annotation to have space
#      scale_y_continuous(breaks = c(-4, -2, -1, 0, 1, 2, 4),
#                         trans = scales::pseudo_log_trans(sigma = 0.25)) + 
#      scale_x_discrete(labels=c("Down Regulated" = "MT down genes \n (Green et al., 2015)", 
#                                "Up Regulated"= "MT up genes \n (Green et al., 2015)")) + 
#      ylab(expression(atop("Differences in "*Log[2]*" Expression", "(Deconvolution B cells)"))) + 
#      #ggtitle("CREBBP Mutation")  + 
#      ggtitle("")  + 

#      theme_classic() + # later operation will overide this
#      theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
#      # label size
#      theme(axis.title.x = element_blank()) + 
#      theme(axis.title.y = element_text(size = title.size - 2)) + 
#      #ticks size
#      theme(axis.text.y = element_text(size = lab.size  , face="bold")) +
#      theme(axis.text.x = element_text(size = lab.size + 5, face="bold"))+
#      theme(legend.text = element_text(size=lab.size ))

# options(repr.plot.width = 9, repr.plot.height = 6, repr.plot.res = 150)
# mp



# h.loc  = linspace(0.68, 1.32, 3)
# v.loc  = exp(linspace(0.25, 1.25, 2))  
# v.size = linspace(0.5, 1, 2)
# p.size = 5

# df1 <- data.frame(a = c(h.loc[1],           h.loc[1], h.loc[3], h.loc[3]), 
#                   b = c(v.loc[1]-v.size[1], v.loc[1], v.loc[1], v.loc[1]-v.size[1]))

# df2 <- data.frame(a = c(h.loc[2],           h.loc[2], h.loc[3], h.loc[3]), 
#                   b = c(v.loc[2]-v.size[2], v.loc[2], v.loc[2], v.loc[2]-v.size[2]))



# df3 <- data.frame(a = 1 + (c(h.loc[1],           h.loc[1], h.loc[3], h.loc[3])), 
#                   b = c(-v.loc[1]+v.size[1], -v.loc[1], -v.loc[1], -v.loc[1]+v.size[1]))

# df4 <- data.frame(a = 1 + (c(h.loc[2],           h.loc[2], h.loc[3], h.loc[3])), 
#                   b = c(-v.loc[2]+v.size[2], -v.loc[2], -v.loc[2], -v.loc[2]+v.size[2]))


# mp =  mp + geom_line(data = df1, aes(x = a, y = b), inherit.aes = FALSE) +
#           geom_line(data = df2, aes(x = a, y = b), inherit.aes = FALSE) +
     
       
#           geom_line(data = df3, aes(x = a, y = b), inherit.aes = FALSE) +
#           geom_line(data = df4, aes(x = a, y = b), inherit.aes = FALSE) +
        
#           annotate("text", x = (h.loc[1]+ h.loc[3])/2,       y = v.loc[1] + v.size[1],      label = parse(text = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "CIBERSORTx"],   digits = 1))), size = p.size) + 
#           annotate("text", x = (h.loc[2]+ h.loc[3])/2,       y = v.loc[2] + v.size[2],      label = parse(text = paste0("\'P = \'~", scientific_pvals_10x(pvals.dec["Unico", "TCA"], digits = 1))), size = p.size) + 
       
#           annotate("text", x = 1 + (h.loc[1]+ h.loc[3])/2,   y = -(v.loc[1] + v.size[1]),    label = parse(text = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "CIBERSORTx"],   digits = 1))), size = p.size) + 
#           annotate("text", x = 1 + (h.loc[2]+ h.loc[3])/2,   y = -(v.loc[2] + v.size[2]),    label = parse(text = paste0("\'P = \'~", scientific_pvals_10x(pvals.inc["Unico", "TCA"], digits = 1))), size = p.size) 
        

# options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 150)
# mp
# # ggsave(file.path(figure.dir, "CREBBP.main.DE.box.png"), device = "png", width = 10, height = 5, dpi = 600)

# # save related ggplots and stats

# fig.list = list(main.g = mp, 
#                 w.g = wp, 
#                 de.g = p, 
#                 regress.g = rp)

# stats.list = list(main.df = main.plot.df,
#                   w.df = W.df,
#                   de.plot.df = de.plot.df, 
#                   reg.plot.df = reg.plot.df)


# saveRDS(fig.list,   file.path(figure.dir, paste0("CREBBP_fig.list.rds")))
# saveRDS(stats.list, file.path(figure.dir, paste0("CREBBP_stats.list.rds")))





# # Histogram 

# library(see)

# ggplot(reg.plot.df, aes(x=DE.type, y=beta, fill = method)) + 
#      geom_boxplot(aes(fill = method), side = 'r', color = 'black', 
#                        outlier.size = NA, outlier.color = NA, notchwidth = 0.5, size = 0.7, #thickness of the line 
#                        width = 0.6, position = position_dodge(0.8), alpha = 1) +  
#      geom_hline(yintercept=0, linetype='dotted', col = "black", size=1) + 
    
    
#      scale_fill_manual(values = as.character(method.colors[levels(reg.plot.df$method)]), name=NULL) +  
#      #coord_cartesian(ylim = c(-60,1.25)) + # for annotation to have space
#      #scale_y_continuous(trans = "exp") + 
#      #scale_y_continuous(breaks = c(-1, -0.25, 0, 0.25, 1),
#      #                   trans = scales::pseudo_log_trans(sigma = 0.01)) + 
#      scale_x_discrete(labels=c("Down Regulated" = "CREBBP down-regulated genes \n (lower is better)", 
#                                "Up Regulated"   = "CREBBP up-regulated genes \n (higher is better)")) +

#      ylab(expression(atop("Estimated effect size", "on mutation status"))) + 
#      ggtitle("")  + 
#      #ggtitle("CREBBP Mutation")  + 

#      theme_classic() + # later operation will overide this
#      theme(plot.title = element_text(hjust = 0.5, size = title.size, face="bold")) + 
#      # label size
#      theme(axis.title.x = element_blank()) + 
#      theme(axis.title.y = element_text(size = lab.size + 2, face="bold")) + 
#      #ticks size
#      theme(axis.text.y = element_text(size = lab.size, face="bold")) +
#      theme(axis.text.x = element_text(size = lab.size, face="bold"))+
#      theme(legend.text = element_text(size = lab.size, face="bold"))

# bMIND.betas.dec = reg.plot.df[reg.plot.df$method == "bMIND" & reg.plot.df$"DE.type" == "Down Regulated", "beta"]
# Unico.betas.dec  = reg.plot.df[reg.plot.df$method == "Unico" & reg.plot.df$"DE.type" == "Down Regulated", "beta"]

# boxplot.stats(Unico.betas.dec)

# min(max(Unico.betas.dec), quantile(Unico.betas.dec, 0.75) + 1.5 * IQR(Unico.betas.dec))

# max(min(Unico.betas.dec), quantile(Unico.betas.dec, 0.25) - 1.5 * IQR(Unico.betas.dec))

# median(Unico.betas.dec)

# median(bMIND.betas.dec)

# tmp.df = reg.plot.df[reg.plot.df$"method" %in% c("Unico", "bMIND") & reg.plot.df$"DE.type" == "Down Regulated",]

# tmp.df

# mu <- ddply(tmp.df, "method", summarise, grp.mean=median(beta))
# tmp.g = ggplot(tmp.df, aes(x=beta, color=method, fill=method)) +
#         geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=100)+
#         geom_density(alpha=0.2)+
#         geom_vline(data=mu, aes(xintercept=grp.mean, color=method),linetype="dashed")+
#         scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#         scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#         scale_x_continuous(breaks = c(-1, -0.25, 0, 0.25, 1),
#                            trans = scales::pseudo_log_trans(sigma = 0.01)) + 
     
#         theme_classic()

# options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 150)
# tmp.g




