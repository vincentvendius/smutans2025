library(tidyverse)
library(ggtree)
library("ggtree")
library("treeio")
library(BactDating)
library(ape)
library(aplot)









###if you want to look at statistics for the samples
sampleInfo <-read.delim("PhD/community_analysis_09-22/tax_pipeline/20230101.diseases.sampleInfo.tsv") %>% 
  separate(groupLabel,c("Country","age_period","misc"),sep="_",remove = FALSE,extra = "merge")%>%
  select("sampleId","clusterLabel","country","latitude","longitude","ageAverage","region","age_period")
sampleInfo<-rbind(sampleInfo,c("CGG100272",NA,"Iran",32.000,53.000,as.double(2480.0),"CentralAsia","IronAge"))
sampleInfo<-rbind(sampleInfo,c("CGG101233",NA,"Afghanistan",47.20,13.20,as.double(820),"CentralAsia","IronAge"))
sampleInfo<-rbind(sampleInfo,c("CGG100534",NA,"Denmark",57.05,9.92,as.double(279),"NorthernEurope","Historical"))
summary_statistics <- read.delim("~/PhD/datafiles/summary_statistics.tsv", header=FALSE,
                                 col.names = c("sampleId","assemblyId","coverageP","ani",
                                               "coverageAvg","coveragePExp","coveragePRatio","Heterozygosity"))%>%
  #filter(coverageAvg>5)%>%
  mutate(clonality=Heterozygosity*1000)%>%
  mutate(Clonal=ifelse(clonality>0.1,"Mixed","Clonal"))%>%mutate(Coverage=ifelse(coverageAvg>5,"High","Low"))%>%
  select("sampleId","clonality","Clonal","Coverage")
  #select(c("sampleId","coverageP","coverageAvg","coveragePRatio","clonality","ageAverage","region","age_period","clonal","Quality"))
strep_table <- read.delim("~/PhD/datafiles/strep_mutans_table3")%>%filter(contigId=="NZ_CP044221.1")%>%
  left_join(summary_statistics,by="sampleId")%>% replace_na(list(Clonal = "Unknown", Coverage = "Low"))%>%
  left_join(sampleInfo)%>%mutate(across(where(is.numeric), ~ round(., 2)))%>%
  select(c("sampleId","coverageP","coverageAvg","ani","clonality",
           "BayesianZ","ageAverage","region","age_period","Clonal"))




gunc_result <- read.delim("~/PhD/datafiles/gunc_result.tsv")%>%select("genome",
                                                                      "n_genes_called","n_contigs","contamination_portion",
                                                                      "reference_representation_score","pass.GUNC")


####GET METADATA####
classification <- read.delim("~/PhD/PhylTreeScripts/data/accession_tree2") #read sample classification
country <- read.delim("~/PhD/PhylTreeScripts/data/country2_dat.tsv",col.names = c("NAME","COUNTRY")) #get country names
serotypes<-read.delim("/home/vincentven/PhD/PhylTreeScripts/data/serotypes",col.names=c("NAME","SEROTYPE"))
source<-read.delim("/home/vincentven/PhD/PhylTreeScripts/data/isolation_source",col.names=c("NAME","SOURCE"))
dating <- read.delim("~/PhD/PhylTreeScripts/data/time_dat.tsv",col.names = c("NAME","DATE")) #get ages
dating<-rbind(dating,c("GCF_002355215N1",as.integer(2015)))
dating$DATE<-as.integer(dating$DATE)
combined<-full_join(dating,country,by="NAME") %>% full_join(serotypes,by="NAME")%>%full_join(source,by="NAME") #combine the metadata
combined_meta<-full_join(combined,select(classification,c("NAME","TYPE","COV")),by="NAME") #combine this with the samples classification
ancientmeta<-read.delim("~/PhD/community_analysis_09-22/tax_pipeline/20230101.diseases.sampleInfo.tsv", header=TRUE)%>%select("sampleId","country","ageAverage")
colnames(ancientmeta)<-c("NAME","COUNTRY","DATE") #get other sample metadata
ancientmeta<-rbind(ancientmeta,c("CGG100272","Iran",as.double(2480.0)))
ancientmeta<-rbind(ancientmeta,c("CGG101233","Afghanistan",as.double(2675)))
ancientmeta<-rbind(ancientmeta,c("CGG100534","Denmark",as.double(279)))
ancientmeta$DATE<-as.double(ancientmeta$DATE)
ancientmeta$SOURCE<-"Ancient teeth"
ancientmeta$SEROTYPE<-NA
combined_meta$DATE<-2022-as.double(combined_meta$DATE) #correct the dating with the other metadata
combined_final<-left_join(combined_meta,ancientmeta,by="NAME")%>%mutate(COUNTRY = coalesce(COUNTRY.x, COUNTRY.y)) %>% 
  mutate(DATE = coalesce(DATE.x, DATE.y)) %>%mutate(SOURCE = coalesce(SOURCE.x, SOURCE.y))%>%mutate(SEROTYPE = coalesce(SEROTYPE.x, SEROTYPE.y))%>%
  select(NAME, COUNTRY, DATE,TYPE,COV,SEROTYPE,SOURCE) %>% mutate_if(is.character, ~na_if(., '')) #combine all the sample metadata into one file


library("ggtree")
library("treeio")
library(ape)
library(patchwork)
#read the raxml tree into a ggtree file
### Get unadultered raxml tree for mutans with bootstrap values
#tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/cluster3/cluster_3_mutans.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpregub.support")
#tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/cluster3/cluster_3_mutans.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpostgub.support")
#tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/cluster3/cluster_3_mutans.final_bootstrapped_tree.tre")
#tree_mutans<-drop.tip(tree,"NEO105")
tree <- read.tree(paste("~/PhD/datafiles/placement_new/","merge.newick",sep=""))
tree_mutans <- drop.tip(tree, "GCF_002355215N1")
#tree_mutans$node.label<-as.numeric(tree_mutans$node.label)
missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
metadata_cluster <- read.delim("~/PhD/datafiles/PCA_mutans/PCA_postgub/postgub_clusters.txt")%>%select(sample_id,group)%>%left_join(metadata,by=c("sample_id"="NAME"))
tree_mutans_full<-tree_mutans%>%left_join(metadata_cluster,by=c("label"="sample_id"))
p3<-ggtree(tree_mutans_full, aes(color=as.numeric(label))) + 
  scale_color_gradientn(name="bootstrap support",colours=c("red", 'orange', 'green', 'cyan', 'blue'))+
  #geom_rootedge(rootedge = 0.0005) +geom_tippoint(
  #mapping = aes(
  #  shape = COUNTRY
  #))+geom_tiplab(aes(subset=TYPE=="HIGHCOV"),size = 3)+
  #ggtitle("RaXML-NG: SNPs-post-masking")
  #ggtitle("RaXML-NG: masked core genome") + theme(legend.position = "none")
  #ggtitle("Gubbins: core genome") + theme(legend.position = "none")
  ggtitle("Principal clusters for RaxML-NG: masked core genome")+ theme(legend.position = "none")
  #ggtitle("RaXML-NG: core genome")+ theme(legend.position = "none")
p3+geom_rootedge(rootedge = 0.0005) +geom_tippoint(
  mapping = aes(
    shape = COUNTRY
  ))+geom_tiplab(aes(subset=TYPE!="REF"),size = 3)
p+p1+p2+p3+plot_annotation(title="Bootstrap values for internal nodes for trees mutans")  
#p+p1+p2+plot_annotation(title="Bootstrap values for internal nodes for ancient mutans clade")  

cluster_meta<-metadata_cluster%>%select("group")%>%as.matrix
row.names(cluster_meta)<-metadata_cluster$sample_id

gheatmap(p3,cluster_meta, offset=0.00000001, width=0.1, 
         colnames=FALSE, legend_title="group") 
p3<-ggplot(tree_mutans, aes(x, y)) + geom_tree() + theme_tree() + geom_tiplab(size = 2)+geom_rootedge(rootedge = 0.01)
#looks like bad bootstrap values for a lot of the deeper clades, try again, but now with the core genomes included
tree <- read.tree("/home/vincentven/PhD/PhylTreeScripts/data/strep_mutans_dating.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxml.support")
tree_mutans <- drop.tip(tree, "GCF_002355215N1")
ggtree(tree_mutans) + geom_nodelab(geom='label',size = 2)+ggtitle("Boostrap values for nodes for core genome using GTR+G+I")

#now lets try and coorborate that with a bootstrap density tree
btrees <- read.tree(system.file("extdata/RAxML", 
                                "RAxML_bootstrap.H3", 
                                package="treeio")
)
btrees <- read.tree("/home/vincentven/PhD/PhylTreeScripts/data/raxml_lewis/streptococcus_mutans.GCF_009738105.1_ASM973810v1.raxml.bootstraps")
ggtree(btrees) + facet_wrap(~.id, ncol=10)
tree_mutans <- drop.tip(btrees, "GCF_002355215N1")

ggdensitree(btrees, alpha=.3, colour='steelblue') +
    geom_tiplab(size=3) + hexpand(.35)+ggtitle("Merged bootstrap trees with branches ordered by most common topology")

missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))


ggtree(tree_mutans_full, layout='circular') + 
  geom_rootedge(rootedge = 0.005)+geom_tippoint(
  mapping = aes(
    #shape = TYPE, 
    color = DATE))+
    #shape = TYPE
  #))+
  geom_tiplab(size = 2,offset=.01,)
###


#Visualize likelihood weights for placements
lowcov<-c("CGG100272","CGG100534","CGG_2_18521","NEO938","RISE373","NEO31","DA112","DA246","DA223","VK146","DA195","NEO65","RISE569",
          "RISE540","VK35","RISE349","DA104","RISE413","NEO170","NEO160","VK441","S_troglodytae")
#NEO065 og VK441 excluded, bringing list to 19.
lowcov<-c("RISE540","NEO65","DA246")
lowcov<-c("NEO170","NEO160","VK441","RISE349")
lowcov<-c("DA246")

lowcov<-c("S_troglodytae")
#check different placements on the tree
pdf(paste0("~/PhD/datafiles/placement_4/placementunrooted.pdf"))#, width=8, height=25)
pdf(paste0("~/PhD/datafiles/placement_4/troglodytae.pdf"))#, width=8, height=25)
for(i in 1:length(lowcov)){
  sample<-lowcov[i]
  tree <- read.tree(paste("~/PhD/datafiles/placement_5/",sample,
                          ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.newick",sep=""))
  d <- read.jplace(paste("~/PhD/datafiles/placement_5/",sample,
                         ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.jplace",sep=""))
  #tree_mutans <- drop.tip(tree, "GCF_002355215N1")
  tree_mutans<-tree
  #d<-drop.tip(d, "GCF_002355215N1") #remove the troglodytae root
  missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
  metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
  tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
  p<-ggtree(tree_mutans_full) +#, layout='circular' 
    geom_rootedge(rootedge = 0.005)+#geom_tippoint(
    #  mapping = aes(
        #shape = TYPE, 
        #color = DATE, 
        #shape = TYPE
     # ))+
     geom_tiplab(aes(subset=COV=="LOW"),size = 2)+ggtitle(paste(sample," placement",sep=""))#+ # (major allele)
    #coord_cartesian(xlim=c(0.43,0.45))
  #print(p)
  #jplace loglikelihood tree
  d@data$llr <- 0 #get another int parameter in the data for the tree
  d@data$llr[match(d@placements$node, d@data$node)] <- d@placements$like_weight_ratio #append the loglikelihood
  #so each node now has a corresponding value for its likelihood
  d@data$hasPlacement <- ifelse(d@data$llr > 0, "yes", "no")
  #also apppend node information that is just binary
  #d@phylo <- root(d@phylo, "NEO137") #root by neo137 which has shown to be basal to the other strains of mutans 
  idxL <- c("solid", "dashed", NA) #associate a line style to each possible binary value for whether
  #a node is associated with the placement of a sample 
  names(idxL) <- c("yes", "no", NA)
  h <- length(d@phylo$tip.label) %/% 10 #h is just a more compressed version of the tip lengths 
  missingtaxa<-setdiff(combined_final$NAME,c(d@placements$name,d@phylo$tip.label))
  metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
  d<-d%>%left_join(metadata,by=c("label"="NAME"))
  p <- ggtree(d, layout="equal_angle", aes(x = x,
                     y = y,
                     color = llr,
                     alpha = llr,
                     size = llr, layout="equal_angle",
                     linetype = hasPlacement))
  print(p+#%>%collapse(node=432)%>%collapse(node=352) +
          #geom_label(aes(x=branch,label=node),size=2,col='black') +
          geom_tiplab(aes(subset=TYPE!="REF"),size = 1.5,
                      hjust = 0,
                      vjust = 0.5,
                      alpha = 1,
                      color = "black",
                      show.legend = FALSE) +
          scale_colour_gradient2(name = "Likelihood weight ratio",
                                 low = "grey80",
                                 mid = "steelblue4",
                                 high = "firebrick1",
                                 midpoint = 0.5,
                                 limits = c(0, 1)) +
          scale_size_continuous(name = "Likelihood weight ratio",
                                range = c(0.15, 1),
                                limits = c(0, 1)) +
          scale_alpha_continuous(range = c(1, 1),
                                 guide = "none") +
          scale_linetype_manual(values = c(2, 1)) +
          theme_tree() +
          theme(legend.position = "right") +
          guides(fill = guide_legend(override.aes = list(size = 2.5))) +
          ggtitle("Placement likelihood for ",sample))#+coord_cartesian(xlim=c(0.43,0.45))

}
dev.off()


###Get the tree in a pretty format with placements 
tree <- read.tree(paste("~/PhD/datafiles/placement_5/","merge.newick",sep=""))
tree$node.label<-1:tree$Nnode
tree$node<-1:tree$Nnode
tree_mutans <- drop.tip(tree, "GCF_002355215N1")
#tree_mutans<-rooted
missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
p<-ggtree(tree_mutans_full) + geom_nodelab(geom='label',size = 2)
library(castor)
count_tips_per_node(tree_mutans_full@phylo)[153]
count_tips_per_node(tree_mutans_full@phylo)[94]
count_tips_per_node(tree_mutans_full@phylo)[78]
View(p$data)
p
#ggtree(tree_mutans_full)+geom_label(aes(x=branch, label=node), fill='lightgreen') 
p<-ggtree(tree_mutans_full) + 
  #geom_label(aes(x=branch),fill='lightgreen') +
  geom_rootedge(rootedge = 0.005)+geom_tippoint(
    mapping = aes(subset=TYPE!="REF",
      #shape = COUNTRY, 
      color = DATE,size = 2,show_guide = FALSE,guide="none",show.legend=FALSE
    ))+geom_tiplab(aes(subset=TYPE=="HIGHCOV"),size = 3,offset=.0008)+guides(size=FALSE)+
  ggtitle("RaxML-ng phylogenetic tree including placement samples",subtitle = "Collapsed nodes with only modern samples: 103 (Green circle), 49 (Red rhombus), 10 (Blue Square)")
p2<-p%>%collapse(node=437)+geom_point2(aes(subset=(node==437)), shape=21, size=8, fill='green')
p2 <- collapse(p2, node=378) + 
  geom_point2(aes(subset=(node==378)), shape=23,size=5, fill='red')
p2<- collapse(p2, node=362) + 
  geom_point2(aes(subset=(node==362)), shape=22,size=5, fill='blue')
#geom_label(aes(x=branch,label=node),fill='lightgreen')
p2

#We see two clades with ancient high coverage samples.
#let us try and extract the clade with four samples in it

#now i am trying to find where to root my tree
tree <- read.tree(paste("~/PhD/datafiles/placement_5/","S_troglodytae",
                        ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.newick",sep=""))
tree<-root(tree,"Troglodytae",resolve.root = TRUE)%>%drop.tip("Troglodytae")
tree$node.label<-1:tree$Nnode
subtree<-get_subtree_at_node(tree, )
#tree <- read.tree(paste("~/PhD/datafiles/placement_3/","merge.newick",sep=""))
#tree$node.label<-1:tree$Nnode
#tree<- drop.tip(tree, "GCF_002355215N1")
#ggtree(tree_mutans_full) + geom_nodelab(geom='label',size = 2)+geom_tiplab()
subtree<-get_subtree_at_node(tree, 130)
tree_mutans<-subtree$subtree
missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
ggtree(tree) + geom_nodelab(geom='label',size = 2)+geom_tiplab()
p<-ggtree(tree_mutans_full) + 
  #geom_label(aes(x=branch),fill='lightgreen') +
  geom_rootedge(rootedge = 0.005)+geom_tippoint(
    mapping = aes(
      shape = TYPE, 
      color = DATE, 
      #shape = TYPE
    ))+geom_tiplab()#geom_tiplab(aes(subset=TYPE=="HIGHCOV"),size = 3,offset=.001)+ggtitle(paste("High quality ancient genomes"," major allele",sep=""))
p


#Csv file describing branch statistics with regards to recombination
BranchStats <- read.delim("~/PhD/PhylTreeScripts/data/gubbins/strep_mutans_dating.GCF_009738105.1_ASM973810v1.genotyped.SNP.per_branch_statistics.csv")
### Do molecular dating using high coverage samples with gubbins
library(BactDating)


tree<-loadGubbins("~/PhD/PhylTreeScripts/data/gubbins/strep_mutans_dating.GCF_009738105.1_ASM973810v1.genotyped.SNP")
is.rooted(tree)
tree <- drop.tip(tree, "GCF_002355215N1")
is.rooted(tree)
ggtree(phyl_postgub) + geom_nodelab(geom='label',size = 2)

tree <- read.tree(paste("~/PhD/datafiles/placement_5/","merge.newick",sep=""))
tree_trog<-root(tree,"GCF_002355215N1")#,node = "Node_263",resolve.root = TRUE)
#tree<-loadGubbins("~/PhD/PhylTreeScripts/data/gubbins_withOutgroup/strep_mutans_dating.GCF_009738105.1_ASM973810v1.genotyped.SNP")
tree_2<-read.tree("~/PhD/PhylTreeScripts/data/gubbins/strep_mutans_dating.GCF_009738105.1_ASM973810v1.genotyped.SNP.node_labelled.final_tree.tre")
#tree <- read.tree(paste("~/PhD/datafiles/placement_major/","merge.newick",sep=""))
tree_2 <- drop.tip(tree_2, "GCF_002355215N1")
library(castor)
subtree<-get_subtree_with_tips(tree,only_tips = c("RISE373","NEO938","CGG100272","NEO137","VK63","NEO105","CGG_2_18521","CGG101233"))


tree_mutans<-tree_trog
tree_mutans<-tree
missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
p<-ggtree(tree_mutans_full) + 
  #geom_label(aes(x=branch),fill='lightgreen') +
  geom_rootedge(rootedge = 0.005)+geom_tippoint(
    mapping = aes(
      #shape = TYPE, 
      color = DATE, 
      shape = TYPE
    ))+geom_tiplab(aes(subset=TYPE!="REF"),size = 3)
p+ggtitle("Gubbins output with ancient samples highlighted")
ggplot(tree_mutans, aes(x, y)) + geom_tree() + theme_tree() + geom_tiplab(size = 2)+geom_rootedge(rootedge = 0.01)

#fix dating metadata
combined_final2<-combined_final
combined_final2$YEAR<-2023-combined_final2$DATE
combined_final2<-filter(combined_final2,NAME%in%intersect(combined_final2$NAME,tree$tip.label))
tiplabels<-cbind(tree$tip.label,1:length(tree$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)
#compute molecular clock
rooted=initRoot(tree,tiplabels$YEAR)
#rooted<-drop.tip(rooted,"GCF_002355215N1")

combined_final2$YEAR<-2023-combined_final2$DATE
combined_final2<-filter(combined_final2,NAME%in%intersect(combined_final2$NAME,rooted$tip.label))
tiplabels<-cbind(rooted$tip.label,1:length(rooted$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)

r=roottotip(rooted,tiplabels$YEAR)
res=bactdate(rooted,tiplabels$YEAR,showProgress = TRUE)#,nbIts=1e5
plot(res,'trace')
plot(res,'treeCI')
print(res$rootdate)
res2=bactdate(rooted,showProgress = TRUE,rep(2015,length(tree$tip.label)))
modelcompare(res,res2)



#lets try and look at the mutans references only
tree<-loadGubbins("~/PhD/PhylTreeScripts/data/gubbins_Clade/mutans")
combined_final2<-combined_final
tree_mutans<-tree
tree_mutans$tip.label<-sub("\\.","N",tree_mutans$tip.label)
missingtaxa<-setdiff(combined_final2$NAME,tree_mutans$tip.label)
metadata<-combined_final2 %>% filter(!NAME%in%missingtaxa)
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
p<-ggtree(tree_mutans_full) + 
  #geom_label(aes(x=branch),fill='lightgreen') +
  geom_rootedge(rootedge = 0.005)+geom_tippoint(
    mapping = aes(
      #shape = TYPE, 
      color = DATE, 
      shape = TYPE
    ))+geom_tiplab(aes(subset=TYPE!="REF"),size = 5,offset=.001)
p
combined_final2$YEAR<-2023-combined_final2$DATE
combined_final2<-filter(combined_final2,NAME%in%intersect(combined_final2$NAME,tree_mutans$tip.label))
tiplabels<-cbind(tree_mutans$tip.label,1:length(tree_mutans$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)
rooted=initRoot(tree_mutans,tiplabels$YEAR)
r=roottotip(rooted,tiplabels$YEAR)
res=bactdate(rooted,tiplabels$YEAR,showProgress = TRUE)
res2=bactdate(tree,showProgress = TRUE, rep(2015,length(rooted$tip.label)))
modelcompare(res,res2)

#looking at the gubbins inferred clade with the ancient samples
#masked 
tree <- read.tree("~/PhD/PhylTreeScripts/data/gubbins/ancient_clade.final_tree.tre")

tree <- read.tree(paste("~/PhD/datafiles/placement_new/","merge.newick",sep=""))
tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_v2_all_snps.GCF_009738105.1_ASM973810v1.raxmlmasked.bestTree.raxml.support")
tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_v2_all_snps.final_bootstrapped_tree.tre")

tree <- drop.tip(tree, "GCF_002355215N1")
is.rooted(tree)
tree$tip.label<-gsub("'","",tree$tip.label)
combined_final2<-combined_final
tree_mutans<-tree
tree_mutans$tip.label<-sub("\\.","N",tree_mutans$tip.label)
missingtaxa<-setdiff(combined_final2$NAME,tree_mutans$tip.label)
metadata<-combined_final2 %>% filter(!NAME%in%missingtaxa)
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
p<-ggtree(tree_mutans_full) + 
  #geom_label(aes(x=branch),fill='lightgreen') +
  geom_rootedge(rootedge = 0.00005)+geom_tippoint(
    mapping = aes(
      #shape = TYPE, 
      color = DATE, 
      shape = COUNTRY
    ))+geom_tiplab(aes(subset=TYPE!="REF"),size = 5,offset=.001)+
  ggtitle("RAXML-NG tree with placements and recombination corrected core")
  #ggtitle("Clade with four ancient samples extracted from Gubbins Phylogenetic Tree")
p
combined_final2$YEAR<-2023-combined_final2$DATE
combined_final2<-filter(combined_final2,NAME%in%intersect(combined_final2$NAME,tree_mutans$tip.label))
tiplabels<-cbind(tree_mutans$tip.label,1:length(tree_mutans$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)
rooted=initRoot(tree_mutans,tiplabels$YEAR)
rooted=tree
r=roottotip(rooted,tiplabels$YEAR)
res=bactdate(rooted,tiplabels$YEAR,showProgress = TRUE)
res2=bactdate(rooted,showProgress = TRUE, rep(2015,length(rooted$tip.label)))
modelcompare(res,res2)
#just more bullshit


#Try and generate metadata for gubbinsplot
tree<-loadGubbins("~/PhD/PhylTreeScripts/data/gubbins/strep_mutans_dating.GCF_009738105.1_ASM973810v1.genotyped.SNP")
#only subset of metadata which fits with the tips included in tree
missingtaxa<-setdiff(combined_final$NAME,tree$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)%>%select("NAME","COUNTRY","DATE","TYPE","SEROTYPE","SOURCE")
#get only most frequent countries
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
source_counts<-metadata%>%group_by(SOURCE)%>%count(SOURCE)%>%arrange(n)%>%tail(n=6)%>%select(SOURCE)
metadata$SOURCE <-ifelse(metadata$SOURCE %in% source_counts$SOURCE,metadata$SOURCE,"Other")
metadata$TYPE <-ifelse(metadata$TYPE=="REF","Modern","Ancient")
#arrange correctly for r script format
metadata$id<-metadata$NAME
metadatacsv<-metadata%>%select("id","COUNTRY","TYPE","SEROTYPE","SOURCE")
write.csv(metadatacsv,"~/PhD/PhylTreeScripts/data/gubbins/metadata2.csv",row.names=FALSE,quote=TRUE)





############NEW ITERATION WITH NO ADDED SAMPLES TO POSTMASKED RAXML############



library("ggtree")
library("treeio")
library(ape)
library(patchwork)
#read the raxml tree into a ggtree file
### Get unadultered raxml tree for mutans with bootstrap values

tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_asm5_4.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpostgub.support")
tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_asm5_4.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpregub.support")
tree<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_asm5_4.final_bootstrapped_tree.tre")

tree_mutans<-tree
missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
metadata_cluster <- read.delim("/home/vincentven/PhD/datafiles/PCA_mutans/PCA_asm5_3/chromopainter.tsv")%>%left_join(metadata,by=c("sample_id"="NAME"))
tree_mutans_full<-tree_mutans%>%left_join(metadata_cluster,by=c("label"="sample_id"))

p2<-ggtree(tree_mutans_full, aes(color=as.numeric(label))) + 
  scale_color_viridis(name="bootstrap support")+
  #scale_color_gradientn(name="bootstrap support",colours=c("red", 'orange', 'green', 'cyan', 'blue'))+
  geom_rootedge(rootedge = 0.0005) +#geom_tippoint(
  #geom_tiplab(aes(subset=TYPE=="HIGHCOV"),size = 3,offset=0.001)+
  #ggtitle("RaXML-NG: SNPs-post-masking")
  ggtitle("RaXML-NG: masked \ncore genome")
  #ggtitle("Gubbins: \ncore genome") + theme(legend.position = "none")
  #ggtitle("Principal clusters for RaxML-NG: masked core genome")+ theme(legend.position = "none")
  #ggtitle("RaXML-NG: \ncore genome")#+ theme(legend.position = "none")
  #ggtitle("RaXML-NG: \ncore genome")+ theme(legend.position = "none")




p1+geom_tiplab(aes(subset=TYPE=="MAG"))
p+p1+p2+plot_annotation(title="Bootstrap values for internal nodes for trees mutans")  
#p+p1+p2+plot_annotation(title="Bootstrap values for internal nodes for ancient mutans clade") 



chromo_clus <- read.delim("/home/vincentven/PhD/datafiles/PCA_mutans/PCA_asm5_3/chromopainter.tsv")

p6<-ggplot(chromo_clus,aes(y=sample_id, x="Cluster",fill=group))+geom_tile()+
  scale_fill_manual(values=c("Cluster_A_ancient"="black","Cluster_A"="#ffce54",
                             "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568",
                             "Cluster_D_ancient"="black","Cluster_E"="#ed5564"),
                             name="Chromopainter \nCluster")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)
p6
#now i am trying to find where to root my tree
tree <- read.tree(paste("~/PhD/datafiles/placement_5/","S_troglodytae",
                        ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.newick",sep=""))
tree$node.label<-1:tree$Nnode
subtree<-get_subtree_at_node(tree, 194)
#find which tips to root on
cladelabs<-subtree$subtree$tip.label[2:length(subtree$subtree$tip.label)]
tree_mutans_rooted<-root(tree_mutans_full,cladelabs,resolve.root = TRUE)
library("viridis")

#get unrooted tree where you can see where my root comes from

p2<-ggtree(tree_mutans_full,  branch.length="none",
           ladderize = FALSE, aes(color=as.numeric(label))) + 
  scale_color_viridis(name="bootstrap support")+
  #geom_tiplab(aes(subset=TYPE=="MAG"),size = 3,offset=0.001)
  #geom_tiplab(aes(color=group),size = 3,offset=0.001)+
  ggtitle("RaXML-NG: masked \ncore genome")
p2  #mapping = aes(
#  shape = COUNTRY

p6%>%insert_left(p,width=7)#%>%insert_right(p5)



#run bactdating
combined_final2$YEAR<-2023-combined_final2$DATE
combined_final2<-filter(combined_final2,NAME%in%intersect(combined_final2$NAME,tree_mutans$tip.label))
tiplabels<-cbind(tree_mutans$tip.label,1:length(tree_mutans$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)
rooted=initRoot(tree_mutans,tiplabels$YEAR)
r=roottotip(rooted,tiplabels$YEAR)
res=bactdate(rooted,tiplabels$YEAR,showProgress = TRUE)
res2=bactdate(rooted,showProgress = TRUE, rep(2015,length(rooted$tip.label)))
modelcompare(res,res2)




# Load the 'ape' package
library(ape)
library(dendextend)
library(phytools)
# Read the two phylogenetic trees
tree2<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_asm5_4.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpregub.support")
# Create a plot with the first tree
#fix dating metadata
combined_final2<-combined_final
combined_final2$YEAR<-2023-combined_final2$DATE
combined_final2<-filter(combined_final2,NAME%in%intersect(combined_final2$NAME,tree$tip.label))
tiplabels<-cbind(tree2$tip.label,1:length(tree2$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)
#compute molecular clock
rooted1=initRoot(tree2,tiplabels$YEAR)
is.rooted(rooted1)
is.ultrametric(rooted1)
ultra1<-force.ultrametric(rooted1)
is.ultrametric(ultra1)
dend1 <- as.dendrogram(ultra1)

tree1<-read.tree("/home/vincentven/PhD/datafiles/treecomparison/strep_mutans_asm5_2.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpostgub.support")
tiplabels<-cbind(tree1$tip.label,1:length(tree1$tip.label))
colnames(tiplabels)<-c("NAME","SEQUENCE")
tiplabels<-as.data.frame(as.matrix(tiplabels))%>%left_join(combined_final2)
rooted2=initRoot(tree1,tiplabels$YEAR)
is.rooted(rooted2)
is.ultrametric(rooted2)
ultra2<-force.ultrametric(rooted2)
is.ultrametric(ultra2)
dend2 <- as.dendrogram(ultra2)

# Create the tanglegram
bojangle<-tanglegram(dend1, dend2, main = "Tanglegram between RaXML-NG with\n masked (right) and unmasked (left) MSA")
plot(bojangle)
#how to plot this bullshit
#tree_ultrametric <- ladderize(chronopl(tree_rooted)$tree)









scov<-read.delim("PhD/datafiles/GCF_009738105.1_ASM973810v1.strep_mutans_asm5_4.concat_gene_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("id","gene","start","end","avgcov","sampleId"))
scds<-read.delim("PhD/datafiles/GCF_009738105.1_ASM973810v1.strep_mutans_asm5_4.concat_cds_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("name","desc","start","end","avgcov","sampleId"))
combined_cov<-full_join(scov,scds,by=c("start","end","avgcov","sampleId"))
scov<-combined_cov
scov$sampleId<-gsub("\\.","N",scov$sampleId)
#complete_sample<-left_join(sample_info,sample_meta,by="sampleId") #get the metadata on to my sample classification table
genes<-c("dltC","ffh","cnm","cnaB","cbpA","cbm")  #define the genes i want to look for
#,"gtfA","gtfB","gtfC"
#"spaP","spaA","spaB","spaC"
genes_of_interest<-scov[which(scov$gene%in%genes),]
combined_final$sampleId<-combined_final$NAME
pathogen<-left_join(scov,combined_final,by="sampleId")

pathogen$normavgcov<-as.factor(ifelse(pathogen$avgcov>0.1,1,0))#ascertain based on average gene coverage if the gene is therr or not
pathogen2<-pathogen%>%left_join(strep_table,by="sampleId")
pathogen2<-bind_rows(pathogen,mutans_samples) #bind age column to table for metadata in heatmap
pathogen3<-pathogen2[grep("mub",pathogen2$gene),]
missing_from_ancient<-pathogen%>%filter(TYPE!="REF")%>%filter(DATE>600)%>%group_by(gene)%>%
  summarise(freq=mean(as.double(normavgcov)-1))%>%filter(freq<0.1)%>%select(gene)%>%as.list
not_missing_from_modern<-pathogen%>%filter(TYPE=="REF")%>%group_by(start)%>%
  summarise(freq=mean(as.double(normavgcov)-1))%>%filter(freq>0.3)%>%select(start)%>%as.list
rare_cds<-intersect(not_missing_from_modern$start,missing_from_ancient$start)
pathogen4<-pathogen[pathogen$start%in%rare_cds,]
genemat<-pathogen4%>%filter(name!="NA")%>%#filter(id!="F5989_RS00160")%>%
  select("sampleId","start","normavgcov")%>%#unique%>%
  spread(key=start,value=normavgcov)#%>%as.matrix
genenames<-pathogen4%>%filter(name!="NA")%>%filter(sampleId=="NEO105")%>%select("sampleId","start","gene","desc")%>%arrange(start)##%>%as.list
rownames<-genemat$sampleId
genemat<-select(genemat,-sampleId)%>%as.matrix
row.names(genemat)<-rownames

metmat<-metadata%>%select("NAME","SOURCE")
rownames<-metmat$NAME
metmat<-select(metmat,-NAME)%>%as.matrix
row.names(metmat)<-rownames

hey<- gheatmap(p8,taxabunmat, offset=0, width=5, 
               colnames=TRUE, legend_title="Gene coverage",
               colnames_angle=90,colnames_offset_y=-1,
               custom_column_labels=genenames$desc,hjust=1) + 
  scale_x_ggtree() + #ggtitle("hello")+
  scale_y_continuous(expand=c(0, 300))
hey
#save with width of 2000 and height 1500
library(ggnewscale)
hey2<-hey + new_scale_fill()
gheatmap(hey2,metmat,offset=53.75,width=.1)




#
library(aplot)


p1 <- ggplot(metadata, aes(NAME, DATE)) + geom_col()+
  #geom_text(aes(label=NAME, y= DATE+.1)) +
  coord_flip() + theme_tree2() + theme(legend.position='none')+ggtitle("Sample Age")
p1
g <- ggtree(tree_mutans_full) +geom_rootedge(rootedge=0.0005)+
  ggtitle("RaxML-ng tree")#+ geom_tiplab(align=TRUE)
p2 <- ggplot(pathogen3, aes(y=sampleId, x=reorder(id,start))) + ggtitle("mub operon coverage") + 
  geom_tile(aes(fill=normavgcov))  + theme_minimal()+ #+ scale_fill_viridis_c()
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  xlab(NULL) + ylab(NULL) + scale_x_discrete(label=genenames$gene) 
p2
p3<-ggplot(metadata,aes(y=NAME, x="Country",fill=factor(COUNTRY,levels = c("Brazil","USA","United Kingdom",
                                                                           "Denmark","China","Other"))))+geom_tile()+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  scale_fill_manual(values=c("Denmark"="#fc7d0b","United Kingdom"="#c85200","USA"="#5fa2ce",
                             "Brazil"="#1170aa","China"="#50A14F","Other"="#57606c"),na.value="#FFFFFF",
                    name="Origin")+
  xlab(NULL) + ylab(NULL)
p3
p4<-ggplot(metadata,aes(y=NAME, x="Source",fill=factor(SOURCE,levels = c("oral cavity","plaque","Ancient teeth",
                                                                          "blood","stool","Other"))))+geom_tile()+
  scale_fill_manual(values=c("Ancient teeth"="#f0bd27","plaque"="#ffda66","oral cavity"="#e39802",
                             "stool"="#51b364","blood"="#b60a1c","Other"="#57606c"),na.value="#FFFFFF",
                    name="Source")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)
p4
p5<-ggplot(metadata,aes(y=NAME, x="Serotype",fill=SEROTYPE))+geom_tile()+
  scale_fill_manual(values=c("C"="#1170aa","E"="#c85200","F"="#c8d0d9",
                             "K"="#ffbc79"),na.value="#FFFFFF",
                    name="Serotype")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)
p5


#do the tree with the clusters 
p1%>%insert_left(p3,width=0.1)%>%insert_left(p4,width=0.1)%>%insert_left(p5,width=0.1)%>%insert_left(p6,width=0.1)%>%insert_left(p8,width=0.8)

p2%>%insert_right(p1,width=0.6)%>%insert_left(p3,width=0.1)%>%insert_left(p4,width=0.1)%>%insert_left(g,width=1)#%>%insert_left(p5,width=0.1)
test<-pathogen3%>%filter(sampleId=="NEO105")%>%select(gene)%>%as.list


melt(cbind(sample=rownames(qualityScores), qualityScores))




######---------------Sample analysis and chromopainter comparisons-----------####
###Get the tree in a pretty format with placements 
library("treeio")

tree <- read.tree(paste("~/PhD/datafiles/placement_5/","merge.newick",sep=""))
tree<-root(tree,"GCF_002355215N1",resolve.root = TRUE)
missingtaxa<-setdiff(combined_final$NAME,tree$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
country_counts<-metadata%>%group_by(SOURCE)%>%count(SOURCE)%>%arrange(n)%>%tail(n=6)%>%select(SOURCE)
metadata$SOURCE <-ifelse(metadata$SOURCE %in% country_counts$SOURCE,metadata$SOURCE,"Other")

p8<-ggtree(tree,branch.length = "none")+geom_rootedge(rootedge=0.0005)+
  ggtitle("Placement merged\nPhylogenetic Tree")
chromo_clus$groupsimp<-gsub("_ancient","",chromo_clus$group)
p9<-ggplot(chromo_clus,aes(y=sample_id, x="Cluster",fill=groupsimp))+geom_tile()+
  scale_fill_manual(values=c("ancient_low_cov"="black","Cluster_A_ancient"="black","Cluster_A"="#ffce54",
                             "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568",
                             "Cluster_D_ancient"="black","Cluster_E"="#ed5564"),
                    name="Chromopainter \nCluster")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)
p9<-p9+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
             legend.key.height = unit(0.3, 'cm'), #change legend key height
             legend.key.width = unit(0.3, 'cm'), #change legend key width
             legend.title = element_text(size=8), #change legend title font size
             legend.text = element_text(size=8)) #change legend text font size
p4<-p4+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
             legend.key.height = unit(0.3, 'cm'), #change legend key height
             legend.key.width = unit(0.3, 'cm'), #change legend key width
             legend.title = element_text(size=8), #change legend title font size
             legend.text = element_text(size=8)) #change legend text font size
p3<-p3+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
             legend.key.height = unit(0.3, 'cm'), #change legend key height
             legend.key.width = unit(0.3, 'cm'), #change legend key width
             legend.title = element_text(size=8), #change legend title font size
             legend.text = element_text(size=8)) #change legend text font size
p5<-p5+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
             legend.key.height = unit(0.3, 'cm'), #change legend key height
             legend.key.width = unit(0.3, 'cm'), #change legend key width
             legend.title = element_text(size=8), #change legend title font size
             legend.text = element_text(size=8)) #change legend text font size



pathogen2<-pathogen%>%left_join(strep_table,by="sampleId")
pathogen5<-pathogen2%>%mutate(normavgcov=ifelse(is.na(coveragePExp),avgcov,ifelse(avgcov/coveragePExp>1,1,avgcov/coveragePExp)))%>%
  mutate(genelength=end-start)%>%mutate(bincov=ifelse(normavgcov>0.4,1,0))%>%
  mutate(anc_cov=ifelse(TYPE!="REF"&bincov==1,3,bincov))%>%
  mutate(gene_coverage=ifelse(TYPE!="REF"&bincov==0,2,anc_cov))
pathogen3<-pathogen5[grep("mub",pathogen2$gene),]
missing_from_ancient<-pathogen5%>%filter(TYPE!="REF")%>%group_by(start)%>%
  summarise(freq=mean(as.double(bincov)))%>%filter(freq<0.2)%>%as.list
not_missing_from_modern<-pathogen5%>%filter(TYPE=="REF")%>%group_by(start)%>%
  summarise(freq=mean(as.double(bincov)))%>%filter(freq>0.3)%>%as.list
rare_cds<-intersect(not_missing_from_modern$start,missing_from_ancient$start)
pathogen4<-pathogen5[pathogen5$start%in%rare_cds,]%>%filter(gene!="protein_coding")




genemat<-pathogen4%>%filter(name!="NA")#%>%select("sampleId","start","desc","normavgcov")
genenames<-pathogen4%>%filter(name!="NA")%>%filter(sampleId=="NEO105")%>%select("sampleId","start","gene","desc")%>%arrange(start)##%>%as.list
unique(genenames$gene)
p10 <- ggplot(genemat, aes(y=sampleId, x=reorder(id,start))) + ggtitle("peripheal coding sequences") + 
  geom_tile(aes(fill=as.character(gene_coverage)))+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_manual(name="Gene coverage",
                    labels=c("No", "Yes","No_ancient","Yes_ancient"),
                    values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=genenames$gene) 
p10<-p10+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
             legend.key.height = unit(0.3, 'cm'), #change legend key height
             legend.key.width = unit(0.3, 'cm'), #change legend key width
             legend.title = element_text(size=8), #change legend title font size
             legend.text = element_text(size=8)) #change legend text font size

genes<-c("dltC","ffh","gtfA","gtfB","gtfC","gtfD","cnm","spaP","spaA","spaB","spaC","cnaB","cbpA","cbm","mutA","mutA'","mutB","mutC","mutD","fruA","fruB","manL",
         "levD","levE","levF","levG","eepA","gbpA","gbpB","gbpC","gbpD","levQ","srtA","brpA","fabM","ftf","spx","comX","comS","comR","ftsN","aguA","sfcA",
         "atpE","atpB","atpF","atpH","atpA","atpG","atpC","atpD","comEC","comEA","comFC","comFA","comA","comGA","comGB","comGC","comGD","comGE","comGF","comGG",
         "comP","comEB","comB","comF","scrA","scrB","scrAB","lacD","lacR","lacA","lacB","lacC","lacF","lacE","lacG","lacX","lacI",
         "levD","msmK","msmG","msmF","msmE","msmR","sppA","frwB","manX","agaB","manY","manZ","wapA", "pepO", "covR","covS","covRS","VicK","VicR","VicX","Smu63c","atlA","ccpA")
genes<-c("gtfA","gtfB","gtfC","spaP","lacD","ccdA","lacD","lacR","lacA","lacB","lacC","lacF","lacE","lacG","lacX","lacI")
#genes<-c("covR","covS","covRS","VicK","VicR","VicX","gbpA","gbpB","gbpC","gbpD","atlA")
genes<-c("msmK","msmG","msmF","msmE","msmR","mutA","mutA'","mutB","mutC","yhgE","esaA","essA","essB","essC","sunT","bacA",
         "mubY","mubX","mubT","mubE","mubG","mubH","mubA","mubB","mubC","mubD","mubZ","mubM","mubI","mubJ","mubP","arcC","aguA","ptcA",
         "macB","lanM","gtfA","gtfB","gtfC","gtfD")
genes<-c("mubY","mubX","gtfA","gtfB","gtfC","msmE")
genes<-c("mubY","mubX","mubT","mubE","mubG","mubH","mubA","mubB","mubC","mubD","mubZ","mubM","mubI","mubJ","mubP","arcC","aguA","ptcA")
pathogen3<-filter(pathogen5,gene%in%genes)
genemat<-pathogen3%>%filter(name!="NA")#%>%select("sampleId","start","desc","normavgcov")
genenames<-pathogen3%>%filter(name!="NA")%>%filter(sampleId=="NEO105")%>%select("sampleId","start","gene","desc")%>%arrange(start)##%>%as.list
unique(genenames$gene)
p11 <- ggplot(genemat, aes(y=sampleId, x=reorder(id,start))) + ggtitle("selected virulence factors\n(Reference mapping)") + 
  #geom_tile(aes(fill=as.character(gene_coverage)))+
  geom_tile(aes(fill=normavgcov))+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Gene coverage",trans = 'reverse',option = "magma")+
  #scale_fill_manual(name="Gene coverage")+
                    #labels=c("No", "Yes","No_ancient","Yes_ancient"),
                    #values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=genenames$gene) 
p11<-p11+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
               legend.key.height = unit(0.3, 'cm'), #change legend key height
               legend.key.width = unit(0.3, 'cm'), #change legend key width
               legend.title = element_text(size=8), #change legend title font size
               legend.text = element_text(size=8)) #change legend text font size
p11
p9%>%insert_left(p8,width=4)%>%insert_right(p3)%>%insert_right(p4)%>%insert_right(p10,widt=6)
pdf("~/PhD/datafiles/mappedgenes.pdf")
p9%>%insert_left(p8,width=8)%>%insert_right(p11,widt=15)
dev.off()
pdf("~/PhD/datafiles/finaltree.pdf")
p9%>%insert_left(p8,width=7)%>%insert_right(p3)%>%insert_right(p5)%>%insert_right(p4)%>%insert_right(p1,width=10)
dev.off()
p9%>%insert_left(p8,width=2)%>%insert_right(p10,width=10)



#looking for atlA
pathogen6<-filter(pathogen5,desc%in%c("KxYKxGKxW signal peptide domain-containing protein",
                                      #"YSIRK-type signal peptide-containing protein",
                                      #"GbpC/Spa domain-containing protein",
                                      #"CHAP domain-containing protein",
                                      #"cell wall metabolism sensor histidine kinase VicK",
                                      #"GBS Bsp-like repeat-containing protein",
                                      #"response regulator YycF","peptidase T",
                                      #"MBL fold metallo-hydrolase",
                                      "LPXTG cell wall anchor domain-containing protein",
                                      "dextranase DexA",
                                      "F0F1 ATP synthase subunit beta"))
pathogen6<-filter(pathogen5,desc=="zinc ABC transporter substrate-binding protein AdcA")
pathogen6<-filter(pathogen5,gene%in%c("atpA","atpF","atpB","atpD","ftsE","ftsX","ftsY","ftsZ","ftsA","ftsL","ftsH","dexA","dexB"))
#looking for lac operon
pathogen6<-pathogen5[grep("lac",pathogen5$gene),]
#looking for com genes
pathogen6<-filter(pathogen5,start>1500413)%>%filter(start<1508500)
#looking for bacteriocin blp family 2
pathogen6<-filter(pathogen5,start>1561928)%>%filter(start<1571840)
pathogen6<-pathogen5[grep("bacteriocin",pathogen5$desc),]%>%filter(gene=="protein_coding")
#look for manL
pathogen6<-filter(pathogen5,gene=="ptsP")
genemat<-pathogen6%>%filter(name!="NA")#%>%select("sampleId","start","desc","normavgcov")
genenames<-pathogen6%>%filter(name!="NA")%>%filter(sampleId=="NEO105")%>%select("sampleId","start","gene","desc")%>%arrange(start)##%>%as.list
unique(genenames$gene)
p15 <- ggplot(genemat, aes(y=sampleId, x=reorder(id,start))) + ggtitle("ftsN,gbpA,dexA,atpD coding sequences\n(Reference mapping)") + 
  #geom_tile(aes(fill=as.character(gene_coverage)))+
  geom_tile(aes(fill=normavgcov))+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Gene coverage",trans = 'reverse',option = "magma")+
  #scale_fill_manual(name="Gene coverage")+
  #labels=c("No", "Yes","No_ancient","Yes_ancient"),
  #values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=genenames$desc) 
p15<-p15+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
               legend.key.height = unit(0.3, 'cm'), #change legend key height
               legend.key.width = unit(0.3, 'cm'), #change legend key width
               legend.title = element_text(size=8), #change legend title font size
               legend.text = element_text(size=8)) #change legend text font size
p9%>%insert_left(p8,width=8)%>%insert_right(p15,widt=15)
#atlA is everywhere

#vicK and yycF
pathogen6<-filter(pathogen5,desc%in%c("LacI family DNA-binding transcriptional regulator",
                                      "glucan 1%2C6-alpha-glucosidase DexB",
                                      "sn-glycerol-3-phosphate ABC transporter ATP-binding protein UgpC",
                                      "sucrose phosphorylase",
                                      "carbohydrate ABC transporter permease",
                                      "sugar ABC transporter permease",
                                      "sugar-binding protein MsmE",
                                      "alpha-galactosidase",
                                      "AraC family transcriptional regulator"))
pathogen6<-filter(pathogen5,start>463172)%>%filter(start<474111)
pathogen6<-pathogen5%>%mutate(alphaparam=ifelse(TYPE!="REF",1,0.6))%>%
  mutate(alphaparam2=ifelse(sampleId=="GCF_009738105N1",1,alphaparam))
genemat<-pathogen6%>%filter(name!="NA")#%>%select("sampleId","start","desc","normavgcov")
genenames<-pathogen6%>%filter(name!="NA")%>%filter(sampleId=="NEO105")%>%select("sampleId","start","gene","desc")%>%arrange(start)##%>%as.list
unique(genenames$gene)
p15 <- ggplot(genemat, aes(y=sampleId, x=reorder(id,start),alpha=alphaparam2)) + 
  ggtitle("Whole genome reference mapping")+#ggtitle("msm operon\n(Reference mapping)") + 
  #geom_tile(aes(fill=as.character(gene_coverage)))+
  geom_tile(aes(fill=normavgcov))+
  scale_alpha_identity(guide = "none")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 10, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Gene coverage",trans = 'reverse',option = "magma")+
  #scale_fill_manual(name="Gene coverage")+
  #labels=c("No", "Yes","No_ancient","Yes_ancient"),
  #values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=genenames$desc) 
p15<-p15+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
               legend.key.height = unit(0.3, 'cm'), #change legend key height
               legend.key.width = unit(0.3, 'cm'), #change legend key width
               legend.title = element_text(size=8), #change legend title font size
               legend.text = element_text(size=8)) #change legend text font size
pdf("~/PhD/datafiles/msm_mapping.pdf")
p9%>%insert_left(p8,width=8)%>%insert_right(p15,widt=350)

dev.off()



#########################HEATMAP ECO#####################
#eco pan ancient sample heatmap
combined_taxabun <- read.table("~/PhD/datafiles/combined_taxabun.tsv", quote="\"", comment.char="",col.names = c("reference","cov","expcov","sample"))%>%
  filter(!sample%in%c("NEO137","NEO105","CGG101233","VK63"))%>%filter(!reference%in%c("CGG100272","CGG100534","NEO938"))
combined_taxabun$NAME<-gsub("\\.","N",combined_taxabun$reference)
#combined_taxabun$count<-1
#sample_counts<-group_by(combined_taxabun,NAME)%>%summarise(count=sum(count))


#trying to look at ECO with heatmap3, utter failue, too many zero values
taxabunmat<-combined_taxabun%>%select("sample","NAME","expcov")%>%
  spread(key=NAME,value=expcov)
rownames<-taxabunmat$sample
taxabunmat<-select(taxabunmat,-sample)%>%as.matrix
row.names(taxabunmat)<-rownames
#colnames(taxabunmat)<-c(colnames(taxabunmat),setdiff(rooted2$tip.label,))
appendmat<-matrix(0,length(row.names(taxabunmat)),length(setdiff(rooted2$tip.label,colnames(taxabunmat))))
colnames(appendmat)<-setdiff(rooted2$tip.label,colnames(taxabunmat))
row.names(appendmat)<-row.names(taxabunmat)
taxabunmat<-cbind(taxabunmat,appendmat)
taxabunmat[is.na(taxabunmat)] <- 0
rooted2=tree_mutans_rooted@phylo
is.rooted(rooted2)
is.ultrametric(rooted2)
ultra2<-force.ultrametric(rooted2)
is.ultrametric(ultra2)
dend2 <- as.dendrogram(ultra2)
any(is.na(taxabunmat))
any(is.null(taxabunmat))
heatmap3(t(taxabunmat),distfun=dist,Rowv = dend2)



##ggplot method


p5<-ggplot(combined_taxabun,aes(y=NAME, x=sample,fill=expcov))+geom_tile(color="black")+
  scale_color_viridis(name="Mapping expected coverage")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)+ggtitle("")
p5


p5<-ggplot(combined_taxabun,aes(y=NAME, x=sample,fill=expcov)) +
  geom_tile(color = "black") +
  scale_fill_viridis(name="Mapping expected coverage")+
  #scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)+ggtitle("Expected coverage against\nancient sample")

p11<-ggtree(tree_mutans_rooted)+#, aes(color=as.numeric(label))) + 
  
  #scale_color_gradientn(name="bootstrap support",colours=c("red", 'orange', 'green', 'cyan', 'blue'))+
  geom_rootedge(rootedge = 0.0005) +#geom_tippoint(
  #geom_tiplab(aes(subset=TYPE=="HIGHCOV"),size = 3,offset=0.001)+
  #ggtitle("RaXML-NG: SNPs-post-masking")
  ggtitle("RaXML-NG:\n masked \ncore genome")
#ggtitle("Gubbins: \ncore genome") + theme(legend.position = "none")
#ggtitle("Principal clusters for RaxML-NG: masked core genome")+ theme(legend.position = "none")
#ggtitle("RaXML-NG: \ncore genome")+ theme(legend.position = "none")
p11

p9%>%insert_left(p11,width=4)%>%insert_right(p5,width=10)
