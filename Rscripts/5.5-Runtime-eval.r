# m = 1000
# k = 3, 4, 5, 6,  n = 500
#     (gran, CD4T CD8T NK B  Mono)
#     (gran, CDT NK B  Mono)
#     (gran, CDT+NK, B  Mono)
#     (gran, Lym  Myeloid )

# n = 100, 250, 500  k = 6 (note that 500 case is covered in previous)


# parameter learning and tensor: 
#     CIBERSORTx
#     Unico
#     TCA
#     bMIND rough prior 

# association (after parameters are learned): 
#     TCA-parametric       (Joint and marginal)
#     Unico-parametric     (Joint and marginal)
#     Unico-asymptotics. (marginal only)

#     cellDMC (Joint and marginal)
#     bMIND Tensor  (Tensor, Joint MANOVA X|Y and marginal Y|X)
#     bMIND MCMC  (MCMC based pval) ask if want to run on server with 30 nodes or something. 
library("lubridate")
library("dplyr")
library(data.table)
set.seed(2023)

res.dir = "../Result/Runtime-local/"

m = 1000
permute.size = 10
configures = list(c(6,  m, 500), c(5,  m, 500), c(4,  m, 500), c(3,  m, 500),
                  c(6,  m, 100), c(6,  m, 250)) #automatically also get the 500  sample size result from previous

extract_time_R_log = function(line){
    brac.s = unlist(gregexpr('\\[', line))
    brac.e = unlist(gregexpr('\\]', line))
    t = substr(line, (brac.s + 1), (brac.e - 1))
    return(parse_date_time(t, 'ymd HMS!'))
}


parse_Unico_log = function(logfile){
    param.s = extract_time_R_log(logfile[grepl("Starting parameter learning*", logfile)])
    param.e = extract_time_R_log(logfile[grepl("Finished parameter learning*", logfile)])

    tensor.s = extract_time_R_log(logfile[grepl("Starting tensor*", logfile)])
    tensor.e = extract_time_R_log(logfile[grepl("Finished tensor estimation*", logfile)])

    ewas.asymp.s = extract_time_R_log(logfile[grepl("Starting asymptotic pvals calculation*", logfile)])
    ewas.asymp.e = extract_time_R_log(logfile[grepl("Finished asymptotic pvals calculation*", logfile)])

    ewas.param.s = extract_time_R_log(logfile[grepl("Starting parametric pvals calculation*", logfile)])
    ewas.param.e = extract_time_R_log(logfile[grepl("Finished parametric pvals calculation*", logfile)])

    param.t      = as.numeric(difftime(param.e,  param.s, units = c("secs")))
    tensor.t     = as.numeric(difftime(tensor.e,  tensor.s, units = c("secs")))
    ewas.asymp.t = as.numeric(difftime(ewas.asymp.e,  ewas.asymp.s, units = c("secs")))
    ewas.param.t = as.numeric(difftime(ewas.param.e,  ewas.param.s, units = c("secs")))

    return(list(param.t  = param.t, tensor.t = tensor.t,
                ewas.asymp.t = ewas.asymp.t, ewas.param.t = ewas.param.t))
}

parse_tca_log = function(logfile){
    param.s = extract_time_R_log(logfile[grepl("Fitting the TCA model*", logfile)])
    if(any(grepl("Internal loop converged*", logfile))){
        param.e = extract_time_R_log(logfile[grepl("Internal loop converged*", logfile)])
    }else{
        param.e = extract_time_R_log(logfile[grepl("Iteration 10 out of 10 internal iterations*", logfile)])
    }
    tensor.s = extract_time_R_log(logfile[grepl("Starting tensor for estimating Z*", logfile)])
    tensor.e = extract_time_R_log(logfile[grepl("Finished estimating tensor.*", logfile)])

    ewas.param.s = extract_time_R_log(logfile[grepl("Calculate p-values for deltas and gammas.*", logfile)])
    ewas.param.e = extract_time_R_log(logfile[grepl("Finished tca.*", logfile)])

    param.t      = as.numeric(difftime(param.e,  param.s, units = c("secs")))
    tensor.t     = as.numeric(difftime(tensor.e,  tensor.s, units = c("secs")))
    ewas.param.t = as.numeric(difftime(ewas.param.e,  ewas.param.s, units = c("secs")))

    return(list(param.t  = param.t, 
                tensor.t = tensor.t,
                ewas.param.t = ewas.param.t))
}


Unico.time = list()
tca.time = list()
for(configure in configures){
    k  = configure[1]
    m  = configure[2]
    n  = configure[3]
    for(t in 1:permute.size){
        config = paste(paste0("k_",k),
                       paste0("m_",m),
                       paste0("n_",n),
                       paste0("t_",t), sep = ".")
        print(paste0("working on simulated data: ", config))

        log.file = file.path(res.dir,"mdl", paste0("Unico.", config, ".log"))

        tca.logfile  = readLines(file.path(res.dir, "mdl", paste0("TCA.", config, ".log")))
        Unico.logfile = readLines(file.path(res.dir, "mdl", paste0("Unico.", config, ".log")))

        tca.time[[config]]   = parse_tca_log(tca.logfile)
        Unico.time [[config]] = parse_Unico_log(Unico.logfile)
    }
}

str(Unico.time)

str(tca.time)

cibersortx.logfile =  readLines(file.path(res.dir, "cibersortx", "CIBERSORTx.log"))
cibersortx.logfile.list = list()

counter = 1
record = F
#read in all the log, parse by each configuration
while(counter < length(cibersortx.logfile)){
    if(!record & grepl("^start", cibersortx.logfile[counter])){#match on start
         record = T
         logfile = list()
         line = cibersortx.logfile[counter]
         config = substr(line , (unlist(gregexpr(" ",line ))[1] +1), (unlist(gregexpr(":",line ))[1] -2))
    }

    if(record){
        logfile = c(logfile, cibersortx.logfile[counter])
    }
    if(record & grepl("^>Running time", cibersortx.logfile[[counter]])){ #match on finished 
        cibersortx.logfile.list[[config]] = logfile
        record = F
    }
    counter = counter + 1
}

# extract time 
cibersortx.time = list()
for (config in names(cibersortx.logfile.list)){
    #print(config)
    logfile = cibersortx.logfile.list[[config]]
    runtime = logfile[[length(logfile)]]
    runtime = as.numeric(substr(runtime, (unlist(gregexpr(":",runtime))[1] + 2), nchar(runtime)))
    cibersortx.time[[config]] = runtime 
}

str(cibersortx.time)

CellDMC.time = list()
bMIND.time = list()
for(configure in configures){
    k  = configure[1]
    m  = configure[2]
    n  = configure[3]
    for(t in 1:permute.size){
        config = paste(paste0("k_",k),
                       paste0("m_",m),
                       paste0("n_",n),
                       paste0("t_",t), sep = ".")
        print(paste0("working on simulated data: ", config))
        bMIND.mdl     =  readRDS(file.path(res.dir, "mdl", paste0("bMIND.tensor.", config, ".rds")))
        CellDMC.mdl   =  readRDS(file.path(res.dir, "mdl", paste0("CellDMC.", config, ".rds")))

        CellDMC.time[[config]] = list(assoc.time = as.numeric(CellDMC.mdl$assoc.time))
        bMIND.time [[config]] = list(tensor.time = as.numeric(bMIND.mdl$tensor.time),
                                     tensor.assoc.time = as.numeric(bMIND.mdl$tensor.assoc.time))

    }
}

str(CellDMC.time)

str(bMIND.time)

tensor.k.df = matrix(0, 0, 3)
colnames(tensor.k.df) = c("time", "k", "method")
asso.k.df = matrix(0, 0, 3)
colnames(asso.k.df) = c("time", "k", "method")

k.configures = list(c(3,  m, 500), c(4,  m, 500), c(5,  m, 500), c(6,  m, 500))
for(configure in k.configures){
    k  = configure[1]
    m  = configure[2]
    n  = configure[3]
    for(t in 1:permute.size){
        config  = paste(paste0("k_",k), paste0("m_",m), paste0("n_",n), paste0("t_",t), sep = ".")

        tensor.k.df = rbind(tensor.k.df, c(cibersortx.time[config ],                                     k, "CIBERSORTx"))
        tensor.k.df = rbind(tensor.k.df, c(bMIND.time[[config]]$tensor.time,                             k, "bMIND"))
        tensor.k.df = rbind(tensor.k.df, c((tca.time[[config]]$param.t  + tca.time[[config]]$tensor.t),  k, "TCA"))
        tensor.k.df = rbind(tensor.k.df, c((Unico.time[[config]]$param.t + Unico.time[[config]]$tensor.t), k, "Unico"))

        asso.k.df = rbind(asso.k.df, c((CellDMC.time[[config]]$assoc.time ),     k, "CellDMC"))
        asso.k.df = rbind(asso.k.df, c((bMIND.time[[config]]$tensor.assoc.time), k, "bMIND"))
        asso.k.df = rbind(asso.k.df, c((tca.time[[config]]$ewas.param.t ),       k, "TCA"))
        asso.k.df = rbind(asso.k.df, c((Unico.time[[config]]$ewas.param.t ),      k, "Unico-parametric"))
        asso.k.df = rbind(asso.k.df, c((Unico.time[[config]]$ewas.asymp.t ),      k, "Unico-asymptotic"))
    }
}

tensor.k.df = data.frame("time"   = as.numeric(tensor.k.df[,"time"]),
                         "k"      = as.numeric(tensor.k.df[,"k"]),
                         "method" = as.character(tensor.k.df[,"method"]))

asso.k.df = data.frame("time"    = as.numeric(asso.k.df[,"time"]),
                       "k"      = as.numeric(asso.k.df[,"k"]),
                       "method" = as.character(asso.k.df[,"method"]))



head(tensor.k.df)

head(asso.k.df)

####################################### constructing dataframe (n is changing) ##########################################
tensor.n.df = matrix(0, 0, 3)
colnames(tensor.n.df) = c("time", "n", "method")
asso.n.df = matrix(0, 0, 3)
colnames(asso.n.df) = c("time", "n", "method")

n.configures = list(c(6,  m, 100), c(6,  m, 250), c(6,  m, 500))
for(configure in n.configures){
    k  = configure[1]
    m  = configure[2]
    n  = configure[3]
    for(t in 1:permute.size){
        config  = paste(paste0("k_",k), paste0("m_",m), paste0("n_",n), paste0("t_",t), sep = ".")
        print(config)
        tensor.n.df = rbind(tensor.n.df, c(cibersortx.time[config ],                                     n, "CIBERSORTx"))
        tensor.n.df = rbind(tensor.n.df, c(bMIND.time[[config]]$tensor.time,                             n, "bMIND"))
        tensor.n.df = rbind(tensor.n.df, c((tca.time[[config]]$param.t  + tca.time[[config]]$tensor.t),  n, "TCA"))
        tensor.n.df = rbind(tensor.n.df, c((Unico.time[[config]]$param.t + Unico.time[[config]]$tensor.t), n, "Unico"))

        asso.n.df = rbind(asso.n.df, c((CellDMC.time[[config]]$assoc.time ),     n, "CellDMC"))
        asso.n.df = rbind(asso.n.df, c((bMIND.time[[config]]$tensor.assoc.time), n, "bMIND"))
        asso.n.df = rbind(asso.n.df, c((tca.time[[config]]$ewas.param.t ),       n, "TCA"))
        asso.n.df = rbind(asso.n.df, c((Unico.time[[config]]$ewas.param.t ),      n, "Unico-parametric"))
        asso.n.df = rbind(asso.n.df, c((Unico.time[[config]]$ewas.asymp.t ),      n, "Unico-asymptotic"))
    }
}

tensor.n.df = data.frame("time" = as.numeric(tensor.n.df[,"time"]),
                         "n" = as.numeric(tensor.n.df[,"n"]),
                         "method" = as.character(tensor.n.df[,"method"]))

asso.n.df = data.frame("time"    = as.numeric(asso.n.df[,"time"]),
                       "n"      = as.numeric(asso.n.df[,"n"]),
                       "method" = as.character(asso.n.df[,"method"]))



head(tensor.n.df)

head(asso.n.df)

runtime.list = list(tensor.k.df = tensor.k.df, 
                    asso.k.df = asso.k.df,
                
                    tensor.n.df = tensor.n.df,
                    asso.n.df = asso.n.df)

saveRDS(runtime.list, file.path(res.dir, "runtime.list.rds"))
