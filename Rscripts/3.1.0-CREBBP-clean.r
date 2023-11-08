# This notebook organize all the source files and data needed for the CREBBP mutation experiment
# all related data are saved under ../Data/RNA/CREBBP/CREBBP.dat.rds
# Data/RNA/CREBBP/Green
    #original array data downloaded from GEO: GSE127462.rds

# Data/RNA/CREBBP/Cibersortx-fig3
    #mutation status: cibersortx.supp.table.3.csv (downloaded from cibersortx paper) 41587_2019_114_MOESM5_ESM.xlsx
    #downloaded from cibersortx website -> Downloads -> Group Level GEPs - FL (Fig. 3b-f) (zip)
        #reference celltype specific mean: Fig3b-f-FL-arrays-groundtruth.RMA.txt  
        #bulk level expression: Fig3b-f-FL-arrays-mixture.txt  
        #major celltype group: Fig3b-f-LM4_merged_classes.txt  
        #reference pannel: LM22.txt
        
#known up and down regulated gene: Green.supp.CREBBP.gene.set.csv(page 37-50 supp converted to csv file)
    #straitified to genes and probes
    #make sure that the up and down regulated set: 
        #unique within each group,
        #no overlap
        #all present in corresponding bulk 
        
#final mixture files of genes and probes are save with the same sample order 
require(GEOquery)
require(Biobase)
library(testit)
library(data.table)
set.seed(2023)

data.dir = file.path("../Data/RNA/CREBBP/")

if(! file.exists(file.path(data.dir,"Green","GSE127462.rds"))){
    gset <- getGEO("GSE127462", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    saveRDS(gset, file.path(data.dir,"Green","GSE127462.rds"))
}else{
    gset = readRDS(file.path(data.dir,"Green","GSE127462.rds"))
}

X.probes = Biobase::exprs(gset)
X.probes

#converting the GSM ids to sample ids in CREBBP.dat$X.probes
get_geo2sample = function(gset){
  geo.ids = gset@phenoData@data$geo_accession
  sample.ids = gset@phenoData@data$source_name_ch1
  sample.ids = vapply(1: length(sample.ids), function(i) strsplit(sample.ids[i], ", ")[[1]][3], character(1))
  geo2sample = as.matrix(sample.ids)
  rownames(geo2sample) = geo.ids
  return(geo2sample)
}

geo2sample = get_geo2sample(gset)
colnames(X.probes) = paste0("FL_", as.character(geo2sample[colnames(X.probes), ]))
X.probes

X.genes = as.matrix(data.frame(fread(file.path(data.dir, "Cibersortx-fig3/Fig3b-f-FL-arrays-mixture.txt")), row.names = 1))
X.genes

LM22     = as.matrix(data.frame(fread(file.path(data.dir,"Cibersortx-fig3/LM22.txt")), row.names = 1))

mus.ref  = as.matrix(data.frame(fread(file.path(data.dir, "Cibersortx-fig3/Fig3b-f-FL-arrays-groundtruth.RMA.txt")), row.names = 1))

data.frame(fread(file.path(data.dir,"Cibersortx-fig3/cibersortx.supp.table.3.csv")))

mutation = as.matrix(data.frame(fread(file.path(data.dir,"Cibersortx-fig3/cibersortx.supp.table.3.csv"))))[,c(2,3)]
#convert mutation entries to same samples.ids format
mutation[,"PID"] = paste0("FL_", as.numeric(mutation[,"PID"]))
rownames(mutation) = mutation[,"PID"]
mutation = mutation[,-1, drop = F]

mutation

samples.mt = unique(rownames(mutation[mutation[,"CREBBP.mutation.status"] == "MT", , drop = F]))
samples.wt = unique(rownames(mutation[mutation[,"CREBBP.mutation.status"] == "WT", , drop = F]))

gene.set = as.matrix(data.frame(fread(file.path(data.dir, "Green","Green.supp.CREBBP.gene.set.csv"))))[,1:9]
gene.set[, "Probe.ID"]  = gsub("‐","-", gene.set[, "Probe.ID"])
gene.set[, "Gene.Symbol"]  = gsub("‐","-", gene.set[, "Gene.Symbol"])
gene.set = gene.set[gene.set[, "Gene.Symbol"] != '---', ]
gene.set = gene.set[gene.set[, "Gene.Symbol"] != '', ]
gene.set

inc.genes = unique(gene.set[gene.set[, "Mutant.Mean"] > gene.set[, "Wild.type.Mean"], ][,"Gene.Symbol"])
dec.genes = unique(gene.set[gene.set[, "Mutant.Mean"] < gene.set[, "Wild.type.Mean"], ][,"Gene.Symbol"])



length(inc.genes)

length(dec.genes)

gene.set = as.matrix(data.frame(fread(file.path(data.dir, "Green","Green.supp.CREBBP.gene.set.csv"))))[,1:9]
gene.set[, "Probe.ID"]  = gsub("‐","-", gene.set[, "Probe.ID"])
gene.set[, "Gene.Symbol"]  = gsub("‐","-", gene.set[, "Gene.Symbol"])
gene.set = gene.set[gene.set[, "Gene.Symbol"] != '---', ]
gene.set = gene.set[gene.set[, "Gene.Symbol"] != '', ]
gene.set

inc.genes = unique(gene.set[gene.set[, "Mutant.Mean"] > gene.set[, "Wild.type.Mean"], ][,"Gene.Symbol"])
dec.genes = unique(gene.set[gene.set[, "Mutant.Mean"] < gene.set[, "Wild.type.Mean"], ][,"Gene.Symbol"])

inc.genes = intersect(inc.genes, rownames(X.genes))
dec.genes = intersect(dec.genes, rownames(X.genes))

print(paste(length(inc.genes), length(dec.genes)))

umb.genes = intersect(inc.genes, dec.genes)
if(length(umb.genes)!= 0){
    message("notice there are shared genes between up and down regulated genes")
    inc.genes = inc.genes [ - which(inc.genes == umb.genes)]
    dec.genes = dec.genes [ - which(dec.genes == umb.genes)]
}
print(paste(length(inc.genes), length(dec.genes)))

# probes
inc.probes = unique(gene.set[gene.set[, "Mutant.Mean"] > gene.set[, "Wild.type.Mean"], ][,"Probe.ID"])
dec.probes = unique(gene.set[gene.set[, "Mutant.Mean"] < gene.set[, "Wild.type.Mean"], ][,"Probe.ID"])
inc.probes = intersect(inc.probes, rownames(X.probes))
dec.probes = intersect(dec.probes, rownames(X.probes))
assert(length(intersect(inc.probes, dec.probes)) ==0)
print(paste(length(inc.probes), length(dec.probes)))

#main text page 4 
MHC.C2.green =  c("HLA-DRA", "HLA-DRB1", 
                  "HLA-DMA", "HLA-DMB", 
                  "HLA-DPA1", "HLA-DQA1", "HLA-DQB1") 
#SI Appendix, Fig. S10 pannel F (some are not present in their provided list)                  
# MHC.C2.green = c("HLA-DMA",  "HLA-DMB", 
#                  "HLA-DOA",  "HLA-DOB",
#                  "HLA-DPA1", "HLA-DPB1", 
#                  "HLA-DQA1", "HLA-DQB1", 
#                  "HLA-DRA" , "HLA-DRB1", "HLA-DRB4")

#Fig. 5 pannel B           
MHC.C2.cibersortx =  c("HLA-DMA",  "HLA-DMB", 
                       "HLA-DOA",  "HLA-DOB",
                       "HLA-DPA1", "HLA-DPB1", 
                       "HLA-DQA1", "HLA-DQB1", "HLA-DQB2", 
                       "HLA-DRA" , "HLA-DRB1", "HLA-DRB6")

CREBBP.dat = list("X.probes"   = X.probes[,colnames(X.genes)], # make sure X.genes and X.probes has the same column sample ids order
                  "X.genes"    = X.genes, 
                  
                  "LM22"       = LM22,
                  "mus.ref"    = mus.ref,
                  "mutation"   = mutation,
                  "samples.mt" = samples.mt,
                  "samples.wt" = samples.wt,
                 
                  "inc.probes" = inc.probes, "inc.genes" = inc.genes,
                  "dec.probes" = dec.probes, "dec.genes" = dec.genes, 
                  
                  "MHC.C2.green" = MHC.C2.green, 
                  "MHC.C2.cibersortx" = MHC.C2.cibersortx)

saveRDS(CREBBP.dat, file.path(data.dir, "CREBBP.dat.rds"))
