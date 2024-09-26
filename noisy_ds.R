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


dsets <- replicate(1000,
                   noisyDat(6, set_N = sample(c(10, 20, 30), 1),
                            noisefraction = 0.3,
                            n.asf = 1L),
                   simplify = F)
targets <- mclapply(dsets, \(x) attributes(x)$target)
dsets <- lapply(dsets, \(x) cbind(x, list(U= rbinom(nrow(x), 1, 0.5))))

grab_outcomes <- function(x){
  asfs <- unlist(strsplit(x, "\\)\\*\\("))
  asfs <- sapply(asfs, function(x) gsub("\\(|\\)", "", x), USE.NAMES = F)
  rhss <- sapply(asfs, rhs, USE.NAMES = FALSE)
  return(rhss)
}

outcomes <- lapply(targets, grab_outcomes)


mods <- mcmapply(\(x,y,...) frscored_cna(x=x, outcome=y, ...),
                 x = dsets,
                 y = outcomes,
                 MoreArgs = list(fit.range = c(1, 0.7)),
                 SIMPLIFY = F)

nnull <- lapply(mods, \(x) !is.null(x)) |> unlist()
mods <- lapply(mods, function(y) if (!is.null(y)){
            y[[1]][y[[1]]$score >= quantile(y[[1]]$score, 0.98, na.rm = T),]})
# modsc <- lapply(mods[nnull], `[[`, 2)
# correct <- vector(length = length(dsets))
# correct[nnull] <- mcmapply(\(x, y) any(is.submodel(x = x, y = y)),
#                     x = modsc, y = targets[nnull],
#                     SIMPLIFY = FALSE)

# correct2 <- vector(length = length(dsets))
#
correct <- mapply(\(x, y) if(is.null(x)) TRUE else any(is.submodel(x = x$condition, y = y)),
                     x = mods, y = targets,
                     SIMPLIFY = FALSE)

# correct[!nnull] <- TRUE

pv_n <- mcmapply(\(data, out) cnaTest(d = data,
                                    outcomes = out,
                                    aggregFn = min)$p_value,
               data = dsets, out = outcomes,
               SIMPLIFY = FALSE)


res_noise <- data.table(ffree = unlist(correct), pval = unlist(pv_n), empty = !nnull)



saveRDS(res_noise, "20pnoise_1000sets.RDS")



#res_noise[ffree==T & empty == F, .N] / 400

res_noise[ffree == T & empty == F, .N]
#

# res_noise[ffree == T, pval] |> hist(breaks=100)
#
#
# res_noise[pval <= 1 & ffree==T, .N] / res_noise[pval <= 1, .N]
(res_noise[ffree == T & empty == F & pval <= 0.05, .N] + res_noise[pval > 0.05, .N]) / 1000

res_noise[ffree == T & empty == F & pval <= 0.05, .N] / res_noise[pval <= 0.05, .N]

res_noise[ffree == T & empty == F & pval >= 0.95, .N] / res_noise[pval >= 0.95, .N]

(res_noise[ffree == T & empty == F & pval >= 0.95, .N] + res_noise[pval < 0.95, .N]) / 1000

res_noise[pval == 1 & ffree == T, .N] / res_noise[pval == 1, .N]

res_noise[ffree == T & empty == F & pval <= 0.05, .N]
# res_noise[pval > 0.05 & ffree == T & empty == F, .N]
# allr <- rbind(res_noise, res_rand[, c("ffree", "pval", "empty")])
#
# allr[ffree == TRUE & pval >= 0.95, .N] / allr[pval >= 0.95, .N]
#
# allr[ffree == F, pval] |> hist()
