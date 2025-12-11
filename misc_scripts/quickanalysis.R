summary_statistics <- read.delim("~/PhD/datafiles/summary_statistics.tsv", header=FALSE,
                                 col.names = c("sampleId","assemblyId","coverageP","ani",
                                               "coverageAvg","coveragePExp","coveragePRatio","Heterozygosity"))
summary_statistics2<-summary_statistics%>%select(sampleId,coveragePExp,coverageAvg,ani,Heterozygosity)%>%arrange(desc(coverageAvg))
print(xtable(summary_statistics2, type = "latex",digits = 6), file = "filename2.tex")

library(ape)
library(funfuns)
library(ggtree)
library(tidyverse)
mash_dat <-read_mash_tri("/home/vincentven/distance_triangle.tsv")
tr<-nj(mash_dat)
ape::write.tree(tr, file='mash_tree.newick')
tree_mutans <- drop.tip(tr, "GCF_002355215.1")
p1<-ggtree(tree_mutans) + geom_nodelab(geom='label',size = 2)+
  ggtitle("Neighbour joining tree\ncomputed from MASH distance matrix")

library(patchwork)


library(viridis)

fastani <- read.delim("/home/vincentven/PhD/Strep_Phylogenetics_09-21/data/mutans_refs/fani_GCFtax.out", header=FALSE, col.names = c("File","Alignment","ANI"))

#File preprocessing of taxonomic info linked to fastani output for species
fastani <- transform(fastani,ANI=1-(ANI/100)) #transform to 1-ANI
histdat<-fastani%>%filter(ANI!=0.00000)
p2<-ggplot(histdat, aes(x=ANI)) + geom_histogram(bins=500)+ggtitle("Histogram of ANI between samples")
p1+p2+plot_annotation(title = 'Quality Control of Modern Streptococcus Mutans samples',tag_levels = 1)



ani_rep <- read.table("porphy_ref_rep_ani.matrix",quote="\"", comment.char="",check.names = FALSE)
