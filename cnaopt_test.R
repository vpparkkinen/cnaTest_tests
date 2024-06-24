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


dsets_co <- replicate(100,
                   noisyDat(6, set_N = sample(c(10, 20, 30), 1),
                            noisefraction = 0.2,
                            n.asf = 1L),
                   simplify = F)
targets_co <- mclapply(dsets_co, \(x) attributes(x)$target)
dsets_co <- lapply(dsets_co, \(x) cbind(x, list(U= rbinom(nrow(x), 1, 0.5))))

grab_outcomes <- function(x){
  asfs <- unlist(strsplit(x, "\\)\\*\\("))
  asfs <- sapply(asfs, function(x) gsub("\\(|\\)", "", x), USE.NAMES = F)
  rhss <- sapply(asfs, rhs, USE.NAMES = FALSE)
  return(rhss)
}

outcomes_co <- lapply(targets_co, grab_outcomes)

mods_co <- mcmapply(\(x,y,...) cnaOpt(x=x, outcome=y,...),
                 x = dsets_co,
                 y = outcomes_co,
                 MoreArgs = list(reduce = "rreduce"),
                 SIMPLIFY = F)

mcmapply(\(x,y,...) cnaOpt(x=x, outcome=y,...),
                 x = dsets_co[90:100],
                 y = outcomes_co[90:100],
                 MoreArgs = list(reduce = "rreduce"),
                 SIMPLIFY = F)



prob_set <- list(dsets_co[[97]])
mcmapply(\(x,y,...) cnaOpt(x=x, outcome=y,...),
                 x = prob_set,
                 y = outcomes_co[[97]],
                 MoreArgs = list(reduce = "rreduce"),
                 SIMPLIFY = F)

cnaOpt(dsets_co[[97]], outcome = outcomes_co[[97]], reduce = "rreduce")
empty <- lapply(mods_co, \(x) is.null(x) | !is.data.frame(x))

ffree <- mapply(\(x, y) if(is.null(x) || !is.data.frame(x)) {TRUE} else any(is.submodel(x = x$condition, y = y)),
                     x = mods_co, y = targets_co,
                     SIMPLIFY = FALSE)
