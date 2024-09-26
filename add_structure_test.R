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
library(tidyr)
library(ggplot2)
#library(data.table)
#library(frscore)
source("cnaTest_V2.R")

cores <- detectCores() - 2L
options(mc.cores = cores)

nfac <- 6
N = 32
outcome = "A"
base_noise_dats <- bs_dat_create(Nsets = 10, N = N, varnum = nfac)

base_noise_pvals <- mclapply(base_noise_dats, 
                           \(x) cnaTest(x, 
                                        outcomes = outcome,
                                        aggregFn = min)$p_value)

# targets <- replicate(length(base_noise_dats),
#                             randomAsf(nfac, outcome = outcome), simplify = FALSE)

cleandats <- replicate(length(base_noise_dats), randomDat(nfac, 
                                                condtype = "asf",
                                                outcome = outcome),
                       simplify = FALSE)

clean_increment <- 8
inc_nrow_factor <- N / clean_increment
res <- vector("list", inc_nrow_factor)

for(i in seq_along(res)){
  noisy_remove_idx <- sample(1:N, clean_increment * i)
  clean_to_add <- lapply(cleandats, \(x) x[sample(1:N, clean_increment * i), ])
  noise <- lapply(cleandats, \(x) x[-noisy_remove_idx, ]) 
  noisy_and_clean <- mapply(rbind, noise, clean_to_add, SIMPLIFY = FALSE) 
  res[[i]] <- mclapply(noisy_and_clean,
                       \(x) cnaTest(x,
                                    outcomes = outcome,
                                    aggregFn = min)$p_value)
}

res_all <- lapply(res, unlist)
base_all <- unlist(base_noise_pvals)

comb_all <- cbind(base_all, as.data.frame(res_all))
colnames(comb_all) <- seq(from = 0, to = N, by = clean_increment)   

pval_plot <- ggplot(gather(comb_all), aes(value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~key, scales = "free_x")
