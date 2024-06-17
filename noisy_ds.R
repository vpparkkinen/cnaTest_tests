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


dsets <- replicate(300, 
                   noisyDat(6, set_N = sample(c(20, 30, 40, 50), 1),
                            noisefraction = 0.2,
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


mods <- mcmapply(\(x,y) frscored_cna(x=x, outcome=y)[[1]],
                 x = dsets, y = outcomes, 
                 SIMPLIFY = F)

nnull <- lapply(mods, \(x) !is.null(x)) |> unlist()
mods <- lapply(mods, function(y) if (!is.null(y)){
            y[y$score >= quantile(y$score, 0.98, na.rm = T),]})
modsc <- lapply(mods[nnull], `[[`, 2)
correct <- vector(length = length(dsets))
correct[nnull] <- mcmapply(\(x, y) any(is.submodel(x = x, y = y)),
                    x = modsc, y = targets[nnull],
                    SIMPLIFY = FALSE)

correct[!nnull] <- TRUE

pv <- mcmapply(\(data, out) cnaTest(d = data, 
                                    outcomes = out,
                                    aggregFn = min)$p_value,
               data = dsets, out = outcomes,
               SIMPLIFY = FALSE)

pvals <- pv

res <- data.table(ffree = unlist(correct), pval = unlist(pvals), empty = !nnull)

res[ffree == FALSE & pval <= 1, .N] / res[ffree == FALSE | empty == TRUE, .N] 

res[pval >= .95 & ffree==T, .N]
hist(res[f_positive==T,pval], breaks = 100)


