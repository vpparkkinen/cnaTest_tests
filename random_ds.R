if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

#library(cnaOpt)
library(cnasimtools)
library(doParallel)
library(ggplot2)
library(data.table)
library(frscore)
source("cnaTest_V2.R")

cores <- detectCores() - 2L
options(mc.cores = cores)


dsets_rand <- bs_dat_create(Nsets = 400, 
                       N = c(20, 30, 40, 50),
                       varnum = 5)

mods_rand <- mclapply(dsets_rand, \(x) frscored_cna(x, outcome = "A")[[1]])
#mods <- mclapply(dsets, \(x) frscored_cna(x)[[1]])
ffree <- lapply(mods_rand, \(x) is.null(x)) |> unlist()

pv <- mclapply(dsets_rand, \(data) cnaTest(d = data, 
                                    outcomes = "A",
                                    aggregFn = min)$p_value)

#pv <- mclapply(dsets, \(data) cnaTest(d = data, 
                                      #aggregFn = min)$p_value)


hist(unlist(pv))
res_rand <- data.table(ffree = ffree, pval = unlist(pv))

res_rand[pval <= 0.05 & ffree==F, .N] / res_rand[pval <= 0.05, .N]
hist(res_rand[ffree==T,pval], breaks = 100)

sig <- which(unlist(pv) <= 0.05)
sigres <- mclapply(dsets_rand[sig], \(x) frscored_cna(x, outcome = "A")[[1]])
sum(unlist(lapply(sigres, is.null))) / length(sigres)

unsig <- which(unlist(pv) >= 0.95)
unsigres <- mclapply(dsets_rand[unsig], \(x) frscored_cna(x, outcome = "A")[[1]])
sum(unlist(lapply(unsigres, is.null))) / length(unsigres)




