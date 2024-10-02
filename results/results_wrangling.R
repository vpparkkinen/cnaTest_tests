if(is.na(Sys.getenv("RSTUDIO", unset = NA))){
  setwd(system2("pwd", stdout = TRUE)) # if not in RStudio, assume R runs in
} else {                               # a shell. otherwise assume RStudio
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}
library(data.table)
library(ggplot2)

fr_cor_64 <- readRDS("fr_cor_means_Sep30_23h54m.RDS")
fr_fp_64 <- readRDS("fr_FP_means_Sep30_23h54m.RDS")
pvals_64 <- readRDS("pvals_df_Sep30_23h54m.RDS")

fr_cor_32 <- readRDS("fr_cor_means.RDS")
fr_fp_32 <- readRDS("fr_fp_means.RDS")
pvals_32 <- readRDS("pvals_df_wide.RDS")

plot_all <- function(cor_means_df, fp_means_df, pvals_df){
  fr_cor_long <- melt(data.table(cor_means_df), 
                      measure.vars = names(cor_means_df))
  fr_fp_long <- melt(data.table(fp_means_df), 
                     measure.vars = names(fp_means_df))
  pvals_df <- as.data.table(pvals_df)
  all_fr <- fr_cor_long[fr_fp_long, on = "variable"]
  colnames(all_fr) <- c("signalrows", "cor", "f_pos")
  all_fr[,`:=`(mean_p = colMeans(pvals_df), 
               median_p = apply(pvals_df, 2, median),
               p_signif = apply(pvals_df, 2,
                                \(x) sum(x <= 0.05) / nrow(pvals_df)))]
  # /this to label any non-empty output in 0-signal condition false positive
  all_fr[signalrows == 0, f_pos := f_pos + cor]
  all_fr[signalrows == 0, cor := 0L]
  # this to label any non-empty output in 0-signal condition false positive/
  all_fr_long <- melt(all_fr, id.vars = "signalrows")
  out <- all_fr_long |> ggplot(aes(x = variable, y = value, fill = variable,
                                   label = sprintf("%0.2f", 
                                                   round(value, digits = 2)))) +
    geom_bar(stat = "identity") +
    geom_text(size = 2, nudge_y = 0.05) +
    facet_wrap(~signalrows) +
    theme(panel.spacing.y = unit(0.5, "lines"),
          axis.text.x = element_blank())
  return(out)
}

plot_pval_dist <- function(pvals_df){
  long <- melt(data.table(pvals_df), measure.vars = names(pvals_df))
  out <- long |> ggplot(aes(value)) + 
    geom_histogram(bins = 200) +
    facet_wrap(~variable)
  return(out)  
}


allplot_64 <- plot_all(fr_cor_64, fr_fp_64, pvals_64)

allplot_32 <- plot_all(cor_means_df = fr_cor_32, fp_means_df = fr_fp_32, pvals_df = pvals_32)

plot_pval_dist(pvals_32)

plot_pval_dist(pvals_64)


  