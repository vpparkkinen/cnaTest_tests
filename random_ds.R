if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

library(cnaOpt)
library(cnasimtools)
library(doParallel)
library(ggplot2)
library(data.table)
library(frscore)
source("cnaTest_V2.R")

cores <- detectCores() - 2L
options(mc.cores = cores)


dsets <- bs_dat_create(Nsets = 500, varnum = 5)

mods <- mclapply(dsets, \(x) frscored_cna(x, outcome = "A")[[1]])
mods <- mclapply(dsets, \(x) frscored_cna(x)[[1]])
nnull <- lapply(mods, \(x) !is.null(x)) |> unlist()

pv <- mclapply(dsets, \(data) cnaTest(d = data, 
                                    outcomes = "A",
                                    aggregFn = min)$p_value)

pv <- mclapply(dsets, \(data) cnaTest(d = data, 
                                      aggregFn = min)$p_value)


hist(unlist(pv))
res <- data.table(f_positive = nnull, pval = unlist(pv))

res[pval >= 0.95 & f_positive==T, .N] / 500
hist(res[f_positive==T,pval], breaks = 100)


