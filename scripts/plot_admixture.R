library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)

######## Functions ##############

setcolors <- function(popnames) {
  # Set colors for the ancestral populations
  popnames <- unique(popnames)
  if(length(popnames)<10){
    myColors <- brewer.pal(length(popnames),"Set1")
  } else {
    myColors <- brewer.pal(length(popnames),"Set3")
  }
  if( all(popnames %in% c("EAS","SAS","EUR","NAM","AFR") ) ) {
    names(myColors) <- c("EAS","SAS","EUR","NAM","AFR")
  } else if ( all(popnames %in% c("EAS","SAS","EUR","AMR","AFR") ) ) {
    names(myColors) <- c("EAS","SAS","EUR","AMR","AFR")
  } else {
    names(myColors) <- popnames
  }
  return(myColors)
}

sortsamples <- function(admixture) {
  dt <- admixture %>% 
    gather('pop', 'prob', names(admixture)[-1] ) %>%
    group_by(Sample) %>%
    mutate(likely_assignment = pop[which.max(prob)], assignment_prob = max(prob)) %>%
    select(Sample,likely_assignment, assignment_prob) %>% unique()
  dt <- merge(admixture,dt)
  ordering <- dt[order(dt$likely_assignment, -dt$assignment_prob),]$Sample
  return(ordering)
}

admixturebarplot <- function(admixture) {
  admix <- reshape2::melt(admixture, id.vars=c("Sample"))
  admix$Sample <- factor(admix$Sample, levels=sortsamples(admixture))
  p <- ggplot(admix, aes(x=Sample, y=value, fill=variable)) +
    geom_bar(stat="identity") + 
    theme(legend.position="bottom", legend.title = element_blank()) + 
    ylab("Ancestry fraction") + 
    scale_fill_manual(values=setcolors(names(admixture)[-1]))
  if (nrow(admixture) > 25) {
    p <- p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) + xlab("Samples")
  } else {
    p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x = element_blank())
  }
  return(p)
}

######### Main ###################
# get list of individual admixture results from the snakemake rule input
inputlist <- unlist(snakemake@input, recursive=FALSE)

# Merge ancestry results for all sample into one table
admixture <- bind_rows(lapply(inputlist, fread))

# Write merged ancestry result to the output file listed in the snakemake rule
write.table(admixture,file=snakemake@output[["txt"]], sep="\t", quote=F, row.names=F)

# plot admixture results #
if(nrow(admixture)<10){
  pdfwidth=5
} else {
  pdfwidth=10
}
pdf(snakemake@output[["pdf"]], width=pdfwidth, height=5)
theme_set(theme_bw())
admixturebarplot(admixture)
dev.off()

