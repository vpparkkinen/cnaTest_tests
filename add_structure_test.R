if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

if(!require(cnasimtools)){
  remotes::install_github("vpparkkinen/cnasimtools")
}

library(cnasimtools) 
library(doParallel)
library(tidyr)
library(ggplot2)
#library(data.table)
#library(frscore)
source("cnaTest_V2.R")

cores <- detectCores() - 2L
options(mc.cores = cores)

## setup
nfac <- 6 # how many factors in data sets
num_of_datasets <- 1000 # how many data sets
N = 32 # sample size
outcome = "A" #outcome

# create data sets of pure noise
base_noise_dats <- bs_dat_create(Nsets = num_of_datasets, N = N, varnum = nfac)

# calculate p-vals for pure noise data sets
base_noise_pvals <- mclapply(base_noise_dats, 
                           \(x) cnaTest(x, 
                                        outcomes = outcome,
                                        aggregFn = min)$p_value)

# create data sets that conform to a particular target,
# to be used for injecting structure into pure noise
cleandats <- replicate(length(base_noise_dats), randomDat(nfac, 
                                                condtype = "asf",
                                                outcome = outcome),
                       simplify = FALSE)

# targets for the clean data, not really needed except for checking
# that things worked as intended
tars <- lapply(cleandats, \(x) attributes(x)$target)

## additional check that all worked so far, i.e. that the clean data
## do in fact perfectly conform to some or other target structure
# co <- mapply(\(x, y) condition(x, y), tars, cleandats, SIMPLIFY = FALSE)
# co <- lapply(co, \(x) attributes(x)$info[,c("consistency", "coverage")])
# co <- do.call(rbind, co)
## now check colSums or something to make sure all con/cov are one


clean_increment <- 4 # incrementally add this many rows of structure
inc_nrow_factor <- N / clean_increment

res <- vector("list", inc_nrow_factor)
noisy_and_clean <- vector("list", length(res))

# using the noisy data and the clean data sets,
# create sets of data sets with increasing amount of structure
# calculate p-values
for(i in seq_along(res)){
  noisy_remove_idx <- sample(1:N, clean_increment * i)
  clean_to_add <- lapply(cleandats, \(x) x[sample(1:N, clean_increment * i), ])
  noise <- lapply(base_noise_dats, \(x) x[-noisy_remove_idx, ]) 
  noisy_and_clean[[i]] <- mapply(rbind, noise, clean_to_add, SIMPLIFY = FALSE) 
  res[[i]] <- mclapply(noisy_and_clean[[i]],
                       \(x) cnaTest(x,
                                    outcomes = outcome,
                                    aggregFn = min)$p_value)
}

saveRDS(res, "results/increment_structure_pvals.RDS")
saveRDS(noisy_and_clean, "results/increment_structure_datasets.RDS")

## Additional check that everything is working so far
# cccheck <- vector("list", length(res))
# for(i in seq_along(res)){
#   cccheck[[i]] <- mapply(\(x, y) attributes(condition(x, y))$info[,c("consistency", "coverage")],
#                                             tars, noisy_and_clean[[i]],
#                                             SIMPLIFY = FALSE)
# }

# store the p-values, including ones for pure noise, in a data frame, and plot
res_all <- lapply(res, unlist)
base_all <- unlist(base_noise_pvals)
comb_all <- cbind(base_all, as.data.frame(res_all)) 
levels <- seq(from = 0, to = N, by = clean_increment)
colnames(comb_all) <- levels
comb_all <- pivot_longer(comb_all, names(comb_all))
comb_all$name <- factor(comb_all$name, levels = levels) # just for ordering
                                                        # the plots

# plot p-value distributions per amount of structure in the data sets
pval_plot <- ggplot(comb_all, aes(value)) +
  geom_histogram(bins = 200) +
  facet_wrap(~name)

#pval_plot
ggsave("results/pval_plot.pdf")

#comb_all |> filter(name == 28 & value <= 0.05) |> nrow() 


##### frscored_cna() over the sets of data sets

alldats <- c(list(base_noise_dats), noisy_and_clean)

fr_tempres <- vector("list", length(alldats))

for(i in seq_along(alldats)){
  fr_tempres[[i]] <- lapply(alldats[[i]],
                            \(x) frscored_cna(x, outcome = outcome))
}

gettopq <- function(mods, tquant){
  lapply(mods, function(y) if (!is.null(y)){
    y[[1]][y[[1]]$score >= quantile(y[[1]]$score, tquant, na.rm = T),]})
}

fr_topq <- lapply(fr_tempres, \(x) gettopq(x, 0.98))


checkcorrect <- function(mods, tars){
  mapply(\(x, y) if(is.null(x)) NA else any(is.submodel(x = x$condition, y = y)),
         x = mods, y = tars,
         SIMPLIFY = FALSE)
}

fr_cor <- lapply(fr_topq, \(x) checkcorrect(x, tars))

fr_cor <- lapply(fr_cor, unlist)

fr_res_cor <- lapply(fr_cor, \(x) sum(x, na.rm = T) / num_of_datasets)

fr_res_false_pos <- lapply(fr_cor, \(x) sum(!x, na.rm = T) / num_of_datasets)

fr_cor_df <- data.frame(fr_res_cor)
colnames(fr_cor_df) <- levels

saveRDS(fr_cor_df, "results/fr_cor_df.RDS")

fr_res_FP_df <- data.frame(fr_res_false_pos)
colnames(fr_res_FP_df) <- levels

saveRDS(fr_res_FP_df, "results/fr_res_FP_df.RDS")

