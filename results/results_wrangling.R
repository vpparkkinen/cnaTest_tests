if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}
library(data.table)


fr_cor <- readRDS("fr_cor_means.RDS")
fr_cor_long <- melt(data.table(fr_cor), measure.vars = names(fr_cor))
fr_fp <- readRDS("fr_FP_means.RDS")
fr_fp_long <- melt(data.table(fr_fp), measure.vars = names(fr_fp))
pvals <- readRDS("pvals_df_wide.RDS")

all_fr <- fr_cor_long[fr_fp_long, on = "variable"]
colnames(all_fr) <- c("signalrows", "cor", "f_pos")

pvals_df <- as.data.table(pvals)
#colnames(pvals_df) <- as.character(seq(4, 32, 4))

all_fr[,`:=`(mean_p = colMeans(pvals_df), 
            median_p = apply(pvals_df, 2, median),
            p_signif = apply(pvals_df, 2,
                             \(x) sum(x <= 0.05) / nrow(pvals_df)))]

all_fr_long <- melt(all_fr, id.vars = "signalrows")

all_fr_long |> ggplot(aes(x = variable, y = value, fill = variable,
                          label = sprintf("%0.2f", 
                                          round(value, digits = 2)))) +
  geom_bar(stat = "identity") +
  geom_text(size = 2, nudge_y = 0.05) +
  facet_wrap(~signalrows) +
  theme(panel.spacing.y = unit(0.5, "lines"))

long <- melt(pvals_df, measure.vars = names(pvals_df))
long[variable == 24 & value < 0.05, .N]

long |> ggplot(aes(value)) +
  geom_histogram(bins = 200) +
  facet_wrap(~variable)
  