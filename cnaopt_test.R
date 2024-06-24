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
#library(frscore)
source("cnaTest_V2.R")

cores <- detectCores() - 2L
options(mc.cores = cores)


dsets_co <- replicate(1000,
                   noisyDat(6, set_N = sample(c(30,40,70), 1),
                            noisefraction = 0.1,
                            n.asf = 1L),
                   simplify = F)
targets_co <- mclapply(dsets_co, \(x) attributes(x)$target)
# dsets_co <- lapply(dsets_co, \(x) cbind(x, list(U= rbinom(nrow(x), 1, 0.5))))

grab_outcomes <- function(x){
  asfs <- unlist(strsplit(x, "\\)\\*\\("))
  asfs <- sapply(asfs, function(x) gsub("\\(|\\)", "", x), USE.NAMES = F)
  rhss <- sapply(asfs, rhs, USE.NAMES = FALSE)
  return(rhss)
}

outcomes_co <- lapply(targets_co, grab_outcomes)

mods_co <- mcmapply(\(x,y,...) try(cnaOpt(x=x, outcome=y,...)),
                 x = dsets_co,
                 y = outcomes_co,
                 MoreArgs = list(reduce = "rreduce"),
                 SIMPLIFY = F)

# mcmapply(\(x,y,...) cnaOpt(x=x, outcome=y,...),
#                  x = dsets_co[90:100],
#                  y = outcomes_co[90:100],
#                  MoreArgs = list(reduce = "rreduce"),
#                  SIMPLIFY = F)

empty <- lapply(mods_co, \(x) is.null(x) | !is.data.frame(x))

ffree <- mapply(\(x, y) if(is.null(x) || !is.data.frame(x)) {TRUE} else {any(is.submodel(x = x$condition, y = y))},
                     x = mods_co, y = targets_co,
                     SIMPLIFY = FALSE)

pv_co <- mcmapply(\(data, out) cnaTest(d = data,
                                    outcomes = out,
                                    aggregFn = min)$p_value,
               data = dsets_co, out = outcomes_co,
               SIMPLIFY = FALSE)

res_co <- data.table(ffree = unlist(ffree), pval = unlist(pv_co), empy = unlist(empty))
res_co[ffree == T, .N]
res_co[ffree == T & pval > 0.05, .N]

hist(res_co[, pval])
