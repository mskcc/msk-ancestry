library(ggplot2)
library(reshape2)
library(RColorBrewer)

resultsdir = snakemake@params[["outdir"]]
samples = snakemake@params[["samples"]]

myColors <- brewer.pal(5,"Set1")
names(myColors) <- c("EAS","SAS","EUR","NAM","AFR")

print(samples)
