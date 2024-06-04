if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

library(cnasimtools)
library(moments)
library(doParallel)

source("cnaTest_V2.R")


n.cores <- detectCores() - 2
dcluster <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl = dcluster)


N <- NULL # set sample size
nfac <- 6 # how many factors in the DGS

idat <- replicate(200, randomDat(x = nfac, samplesize = N), simplify = FALSE)
targets <- lapply(idat, \(x) attributes(x)$target)


# pvals: all facs as outcomes, aggregate as min versus mean
pvals_agg_min <- unlist(mclapply(idat, \(x) cnaTest(x, aggregFn = min)$p_value))
pvals_agg_mean <- unlist(mclapply(idat, \(x) cnaTest(x, aggregFn = mean)$p_value))

skewness(pvals_agg_min)
skewness(pvals_agg_mean)

hist(pvals_agg_min, breaks = 100)
hist(pvals_agg_mean, breaks = 200)

quantile(pvals_agg_mean)
quantile(pvals_agg_min)


## sanity check -- pval distributions under null

fun <- function() {data.frame(setNames(replicate(6, 
                                      rbinom(30, 1, 0.5), 
                                      simplify = F), 
                                    LETTERS[1:6]))}
rdats <- replicate(100, expr = fun(), simplify = FALSE)

nullps_agg_mean <- unlist(mclapply(rdats, \(x) cnaTest(x, aggregFn = mean)$p_value))
nullps_agg_min <- unlist(mclapply(rdats, \(x) cnaTest(x, aggregFn = min)$p_value))

hist(nullps_agg_mean, breaks = 400)
hist(nullps_agg_min, breaks = 400)
skewness(nullps_agg_mean)
skewness(nullps_agg_min)
table(nullps_agg_min)
