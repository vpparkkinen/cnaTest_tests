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
source("cnaTest_V2.R")

cores <- detectCores() - 2L
options(mc.cores = cores)


dsets <- replicate(100, 
                   noisyDat(6, set_N = sample(c(30, 40, 50, 60, 70), 1),
                            noisefraction = 0.1,
                            n.asf = 1L), 
                   simplify = F)
targets <- mclapply(dsets, \(x) attributes(x)$target)

grab_outcomes <- function(x){
  asfs <- unlist(strsplit(x, "\\)\\*\\("))
  asfs <- sapply(asfs, function(x) gsub("\\(|\\)", "", x), USE.NAMES = F)
  rhss <- sapply(asfs, rhs, USE.NAMES = FALSE)
  return(rhss)
}

outcomes <- lapply(targets, grab_outcomes)


mods <- mcmapply(\(x,y) cnaOpt(x=x, outcome=y),
                 x = dsets, y = outcomes, 
                 SIMPLIFY = F)

modsc <- mclapply(mods, `[[`, 2)

correct <- mcmapply(\(x, y) any(is.submodel(x = x, y = y)),
                    x = modsc, y = targets,
                    SIMPLIFY = FALSE)

pv <- mcmapply(\(data, out) cnaTest(d = data, 
                                  outcomes = out,
                                  aggregFn = min),
               data = dsets, out = outcomes,
               SIMPLIFY = FALSE)

pvals <- mclapply(pv, \(x) x$p_value)

res <- data.table(correct = unlist(correct), pval = unlist(pvals))

res[pval <= 0.05 & correct==FALSE, .N]



