#### mutans paper 
######load packages
library("treeio")
library(tidyverse)
library(ape)
library(aplot)
library(gplots)
library(ggforce)
library(ggnewscale)
library(ggrepel)
library(ggtree)
library(readxl)
library(patchwork)
library(viridis)
library(phytools)
library(TraMineR)
world <- map_data("world")
set.seed(16)


#Load functions
edit_phylo_branch<- function(tree,highlighted){
  d <- tree$data
  # 1. Determine the desired aligned distance (max x among all tips)
  target_x <- max(d$x[d$isTip])
  
  # 2. For each highlighted tip, compute how much extra length is needed
  extend_df <- d %>%
    filter(isTip, label %in% highlighted) %>%
    mutate(
      xend = target_x,               # aligned target
      extend = xend - x              # needed extension
    ) %>%
    filter(extend > 0)               # only for shorter tips
  
  # 3. Plot: tree + dotted extensions + adjusted labels
  corrected_length <- tree +
    geom_tree() +
    
    # dotted extension lines
    geom_segment(
      data = extend_df,
      aes(x = x, y = y, xend = xend, yend = y),
      linetype = "dotted",
      size = 0.5
    )# +
  
  # labels for aligned tips (placed at xend)
  #geom_text(
  #  data = extend_df,
  #  aes(x = xend, y = y, label = label),
  #  hjust = 0,
  #  vjust = 0.5
  #)
  return(corrected_length)
}

phylo_bootstrap <- function(phylo,highlighted){
  tree<-ggtree(phylo, aes(color=as.numeric(label))) + 
    scale_color_viridis(name="bootstrap support")+
    geom_tree()+geom_point2(aes(subset = (label %in% highlighted & isTip)),
                            color = "red", size = 3)#+
    #geom_tiplab(aes(label = ifelse(label %in% highlighted, label, "")),
    #            hjust = -0.1)
  return(tree)
}

####Phylogenetic trees####
phyl_pregub<-read.tree("~/PhD/misc_scripts/data/treecomparison/strep_mutans_asm5_4.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpregub.support")
phyl_gub<-read.tree("~/PhD/misc_scripts/data/treecomparison/strep_mutans_asm5_4.final_bootstrapped_tree.tre")
phyl_rec<-loadGubbins(prefix='~/PhD/datafiles/gubbins/gubbins_asm5_3/strep_mutans_asm5_4')
phyl_postgub<-read.tree("~/PhD/misc_scripts/data/treecomparison/strep_mutans_asm5_4.GCF_009738105.1_ASM973810v1.raxml.bestTree.raxmlpostgub.support")
phyl_postgub_cf<-read.tree("~/PhD/misc_scripts/data/treecomparison/strep_mutans_asm5_4.cf.tree")
phyl_postgub_trog <- read.tree(paste("~/PhD/misc_scripts/data/placement_5/","S_troglodytae",
                   ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.newick",sep=""))
phyl_postgub_magplacements <- read.tree("~/PhD/misc_scripts/data/placement_5/placementmag/merge.newick")
phyl_postgub_withplacements <- read.tree("~/PhD/misc_scripts/data/placement_5/merge.newick")
phyl_withmags<-root(phyl_postgub_magplacements,"GCF_002355215N1",resolve.root = TRUE)%>%drop.tip("GCF_002355215N1")
phyl_final<-root(phyl_postgub_withplacements,"GCF_002355215N1",resolve.root = TRUE)%>%drop.tip("GCF_002355215N1")

anc_phylo<-c("NEO105","NEO137","VK63","CGG101233")
#read the raxml tree into a ggtree file
### Get unadultered raxml tree for mutans with bootstrap values
tree_pregub<-phylo_bootstrap(phyl_pregub,anc_phylo)+
  ggtitle("RaXML-NG: \ncore genome") +   theme(legend.position = "none")


tree_gub<-phylo_bootstrap(phyl_gub,anc_phylo)+
  ggtitle("Gubbins: \ncore genome") +    theme(legend.position = "none")

tree_postgub<-phylo_bootstrap(phyl_postgub,anc_phylo)+
  ggtitle("RaXML-NG: masked \ncore genome")


# Identify rows in edge matrix that correspond to those tips
highlighted<-c("CGG100272","CGG100534","CGG_2_18521","NEO938","RISE373","NEO31","DA112","DA246","DA223","VK146","DA195","NEO65","RISE569",
               "RISE540","VK35","RISE349","DA104","RISE413","NEO170","NEO160","VK441","S_troglodytae","NEO137","NEO105","VK63","CGG101233")
tip_rows <- which(phyl_final$edge[,2] %in% match(highlighted, phyl_final$tip.label))
phyl_final$edge.length[tip_rows] <- 0
phyl_final$edge.length<-ifelse(phyl_final$edge.length>0.05,new_length,phyl_final$edge.length) #identify branches made very long by cluster of ancient lowcov samples
tree_phyl_final <- ggtree(phyl_final)+geom_rootedge(rootedge=0.0005)+#,branch.length = "none"
  ggtitle("Placement merged\nPhylogenetic Tree")#+xlim(0,0.1)
tree_final<-edit_phylo_branch(tree_phyl_final,highlighted)

tip_rows <- which(phyl_withmags$edge[,2] %in% match(highlighted, phyl_withmags$tip.label))
phyl_withmags$edge.length[tip_rows] <- 0
phyl_withmags$edge.length<-ifelse(phyl_withmags$edge.length>0.05,new_length,phyl_withmags$edge.length) #identify branches made very long by cluster of ancient lowcov samples
tree_withmags<- ggtree(phyl_withmags)+ geom_rootedge(rootedge=0.0005)+#,branch.length= "none"
  ggtitle("Phylogenetic tree with \nplacement of MAGs")

tree_withmags<-edit_phylo_branch(tree_withmags,highlighted)



####GET METADATA####

#ancient sampleinfo
sampleInfo <-read.delim("~/PhD/misc_scripts/data/20230101.diseases.sampleInfo.tsv") %>% 
  separate(groupLabel,c("Country","age_period","misc"),sep="_",remove = FALSE,extra = "merge")%>%
  select("sampleId","clusterLabel","country","latitude","longitude","ageAverage","region","age_period")
sampleInfo<-rbind(sampleInfo,c("CGG100272",NA,"Iran",38.2498,48.2933,as.double(2480.0),"CentralAsia","IronAge"))
sampleInfo<-rbind(sampleInfo,c("CGG101233",NA,"Afghanistan",34.818813, 67.838997,as.double(820),"CentralAsia","Historical"))
sampleInfo<-rbind(sampleInfo,c("CGG100534",NA,"Denmark",57.05,9.92,as.double(279),"NorthernEurope","Historical"))
sampleInfo$latitude<-as.double(sampleInfo$latitude)
sampleInfo$longitude<-as.double(sampleInfo$longitude)
ancientmeta<-sampleInfo%>%select("sampleId","country","ageAverage")
colnames(ancientmeta)<-c("NAME","COUNTRY","DATE") #get other sample metadata
ancientmeta$DATE<-as.double(ancientmeta$DATE)
ancientmeta$SOURCE<-"Ancient teeth"
ancientmeta$SEROTYPE<-NA



####Table 2 and 3####
summary_statistics <- read.delim("~/PhD/misc_scripts/data/summary_statistics.tsv", header=FALSE,
                                 col.names = c("sampleId","assemblyId","coverageP","ani",
                                               "coverageAvg","coveragePExp","coveragePRatio","Heterozygosity"))%>%
  mutate(clonality=ifelse(coverageAvg>5,Heterozygosity*1000,NA))%>%mutate(Clonal=ifelse(clonality>0.1,"Mixed","Clonal"))%>%mutate(Coverage=ifelse(coverageAvg>5,"High","Low"))%>%
  select("sampleId","clonality","Clonal","Coverage")
#summary_statistics<-rbind(summary_statistics,c("VK35",as.double(0.000377591*1000),"Mixed","Low"))
summary_statistics<-rbind(summary_statistics,c("VK35",NA,NA,"Low"))

summary_statistics$clonality<-as.double(summary_statistics$clonality)

strep_table <- read.delim("~/PhD/misc_scripts/data/strep_mutans_table3")%>%filter(contigId=="NZ_CP044221.1")%>%
  left_join(summary_statistics,by="sampleId")%>% replace_na(list(Clonal = "Unknown", Coverage = "Low"))%>%
  left_join(sampleInfo)%>%mutate(across(where(is.numeric), ~ round(., 2)))
strep_table$age_period<-gsub("VikingAge","Historical",strep_table$age_period)
strep_table$age_period<-gsub("IronAge","Iron Age",strep_table$age_period)
strep_table$age_period<-gsub("BronzeAge","Bronze Age",strep_table$age_period)
strep_table$ageAverage<-as.double(strep_table$ageAverage)
strep_table$region<-gsub("CentralAsia","Central Asia",strep_table$region)
strep_table$region<-gsub("NorthernEurope","Northern Europe",strep_table$region)
strep_table$region<-gsub("SoutheastAsia","Southeast Asia",strep_table$region)
strep_table$region<-gsub("CentralEasternEurope","Central Eastern Europe",strep_table$region)
strep_table$region<-gsub("NorthAsia", "North Asia",strep_table$region)
strep_table$region<-gsub("SouthernEurope","Southern Europe",strep_table$region)
strep_table$region<-gsub("WesternAsia","Western Asia",strep_table$region)
strep_table$region<-gsub("WesternEurope","Western Europe",strep_table$region)





#Phylogenetic tree sample metadata
classification <- read.delim("~/PhD/misc_scripts/data/accession_tree2") #read sample classification
country <- read.delim("~/PhD/misc_scripts/data/country2_dat.tsv",col.names = c("NAME","COUNTRY")) #get country names
serotypes<-read.delim("~/PhD/misc_scripts/data/serotypes",col.names=c("NAME","SEROTYPE"))
source<-read.delim("~/PhD/misc_scripts/data/isolation_source",col.names=c("NAME","SOURCE"))
dating <- read.delim("~/PhD/misc_scripts/data/time_dat.tsv",col.names = c("NAME","DATE")) #get ages
dating$DATE<-as.integer(dating$DATE)
chromo_clus <- read.delim("~/PhD/misc_scripts/data/PCA_asm5_3/chromopainter.tsv",col.names = c("NAME","group"))

combined<-full_join(dating,country,by="NAME") %>% full_join(serotypes,by="NAME")%>%full_join(source,by="NAME") #combine the metadata
combined_meta<-full_join(combined,select(classification,c("NAME","TYPE","COV")),by="NAME")%>%full_join(chromo_clus,by="NAME") #combine this with the samples classification
combined_meta$DATE<-2022-as.double(combined_meta$DATE) #correct the dating with the other metadata
strep_metadata<-left_join(combined_meta,ancientmeta,by="NAME")%>%mutate(COUNTRY = coalesce(COUNTRY.x, COUNTRY.y)) %>% 
  mutate(DATE = coalesce(DATE.x, DATE.y)) %>%mutate(SOURCE = coalesce(SOURCE.x, SOURCE.y))%>%mutate(SEROTYPE = coalesce(SEROTYPE.x, SEROTYPE.y))%>%
  select(NAME, COUNTRY, DATE,TYPE,COV,SEROTYPE,SOURCE,group) %>% mutate_if(is.character, ~na_if(., ''))%>% #combine all the sample metadata into one file
  filter(NAME%in%phyl_final$tip.label) #take out modern assemblies or ancient assemblies which were not included in the final version
strep_metadata$YEAR<-2022-strep_metadata$DATE
#column simplifications
country_counts<-strep_metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
strep_metadata$COUNTRY <-ifelse(strep_metadata$COUNTRY %in% country_counts$COUNTRY,strep_metadata$COUNTRY,"Other")
source_counts<-strep_metadata%>%group_by(SOURCE)%>%count(SOURCE)%>%arrange(n)%>%tail(n=6)%>%select(SOURCE)
strep_metadata$SOURCE <-ifelse(strep_metadata$SOURCE %in% source_counts$SOURCE,strep_metadata$SOURCE,"Other")
strep_metadata$groupsimp<-gsub("_ancient","",strep_metadata$group) #do not distinguish between ancient and modern in chromopainter clusters


####Pangenomics####
scov<-read.delim("~/PhD/misc_scripts/data/mutans_gene_cov/GCF_009738105.1_ASM973810v1.strep_mutans_asm5_4.concat_gene_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("id","gene","start","end","avgcov","sampleId"))
scds<-read.delim("~/PhD/misc_scripts/data/mutans_gene_cov/GCF_009738105.1_ASM973810v1.strep_mutans_asm5_4.concat_cds_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("name","desc","start","end","avgcov","sampleId"))
combined_cov<-full_join(scov,scds,by=c("start","end","avgcov","sampleId"))
combined_cov$NAME<-gsub("\\.","N",combined_cov$sampleId) 
gene_cov<-left_join(combined_cov,strep_table,by="sampleId")%>% #combined gene and cds, and correct for average read depth when evaluating gene coverage
  mutate(normavgcov=ifelse(is.na(coveragePExp),avgcov,ifelse(avgcov/coveragePExp>1,1,avgcov/coveragePExp)))%>%
  mutate(alphaparam=ifelse(is.na(coveragePExp),0.8,1))#parameter that highlights ancient samples by increasing transparency of modern assemblies

#BAKTA annotations
CGG100272 <- read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/CGG100272-contigs.tsv", header=TRUE,skip=5)
CGG100272$NAME<-"CGG100272"
VK63 <- read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/VK63-contigs.tsv", header=TRUE,skip=5)
VK63$NAME<-"VK63"
NEO137 <- read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/NEO137-contigs.tsv", header=TRUE,skip=5)
NEO137$NAME<-"NEO137"
NEO105 <- read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/NEO105-contigs-filtered.tsv", header=TRUE,skip=5)
NEO105$NAME<-"NEO105"
NEO938 <- read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/NEO938-contigs.tsv", header=TRUE,skip=5)
NEO938$NAME<-"NEO938"
CGG100534<-read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/CGG100534-contigs.tsv", header=TRUE,skip=5)
CGG100534$NAME<-"CGG100534"
CGG101233<-read.delim2("~/PhD/misc_scripts/data/MAGs-ref_annotations/CGG101233-contigs.tsv", header=TRUE,skip=5)
CGG101233$NAME<-"CGG101233"
BAKTA_modern <- read.delim("~/PhD/misc_scripts/data/mutans_gene_cov/cds_modern.tsv", 
                         col.names = c("sampleId","X.Sequence.Id","Type","Start","Stop","Strand","Locus.Tag","Gene","Product","DbXrefs"))
BAKTA_modern$NAME<-sub("\\.","N",cds_modern$sampleId)
BAKTA_comb<-BAKTA_modern%>%select(!sampleId)%>%rbind(NEO938)%>%rbind(NEO137)%>%rbind(NEO105)%>%rbind(VK63)%>%rbind(CGG100534)%>%
  rbind(CGG101233)%>%rbind(CGG100272) #get bakta annotations for both modern assemblies and ancient mags
#define if a gene is annotated ,it has gene coverage, if the same gene is not found in another sample, it does not have gene coverage
#bakta only includes genes it finds for each sample
BAKTA_comb$gene_coverage<-as.integer(1)
BAKTA_comb<-complete(BAKTA_comb, NAME,Gene, fill = list(gene_coverage = as.integer(0)))%>%subset(Gene != "")%>%
  mutate_if(is.integer, ~replace(., is.na(.), 0))

BAKTA_all<-left_join(BAKTA_comb,strep_metadata,by="NAME") #add metadata
in_ancient<-BAKTA_all%>%filter(TYPE=="MAG")%>%group_by(Gene)%>% #Genes present in one or more ancient MAGs
  summarise(freq=mean(as.double(gene_coverage)))%>%filter(freq>0.1)%>%select(Gene)%>%as.list
non_ref_genes<-BAKTA_all%>%filter(NAME=="GCF_009738105N1")%>%filter(gene_coverage==0)%>%select(Gene) #Genes not present in reference genome
gene_freq<-BAKTA_all%>%group_by(Gene)%>% #Get frequency of each gene across samples
  summarise(freq=mean(as.double(gene_coverage)))
ancient_genes<-intersect(non_ref_genes$Gene,in_ancient$Gene) #genes found in an ancient mag, but not seen in the reference genome
BAKTA_all_fxd<-left_join(BAKTA_all,gene_freq,by="Gene")%>% #add distinction between gene presence in ancient and modern samples to make the ancient samples stand out
  mutate(anc_cov=ifelse(TYPE=="MAG"&gene_coverage==1,3,gene_coverage))%>%
  mutate(anc_gencov=ifelse(TYPE=="MAG"&gene_coverage==0,2,anc_cov))


##########Data Tables
#Table 1
#handcrafted on the keyboard

#Table 2
table1<-strep_table%>%select(c("sampleId","ageAverage","region","age_period","Clonal","coverageP","coverageAvg","ani","clonality",
                               "BayesianZ"))%>%arrange(desc(as.numeric(coverageAvg)))
colnames(table1)<-c("Sample Id","Mean Sample Age Estimate","Geographical Origin","Archaeological Context","Strain multiplicity",
                    "Percentage Reference Coverage","Average Read Depth","Average Nucleotide Identity","Multi Allele Rate","metaDMG Bayesian Z score")
write.csv(table1,"~/PhD/Submission folder/Mutans paper folder/table2.csv",row.names=FALSE,quote=TRUE)

#Table 3
gunc_result <- read.delim("~/PhD/misc_scripts/data/gunc_result.tsv")%>%select("genome",
                                                                              "n_genes_called","n_contigs","contamination_portion",
                                                                              "reference_representation_score","pass.GUNC")%>%arrange(desc(as.numeric(n_genes_called)))
colnames(gunc_result)<-c("MAG","Genes Called","Contigs","Contamination Portion","Reference Representation Score","GUNC pass")
write.csv(gunc_result,"~/PhD/Submission folder/Mutans paper folder/table3.csv",row.names=FALSE,quote=TRUE)


#Table 4
mag_stats <- read.delim("~/PhD/misc_scripts/data/mag_stats.tsv")%>%arrange(desc(as.numeric(total_length)))
write.csv(mag_stats,"~/PhD/Submission folder/Mutans paper folder/table4.csv",row.names=FALSE,quote=TRUE)

#Supplemental Table 1
gene_cov_anc<-gene_cov%>%
  filter(!is.na(coveragePExp))%>%
  select(!c("alphaparam","age_period","Clonal","Coverage","clusterLabel","country",
            "latitude","longitude","ageAverage","region","NAME","pContigsCovered",
            "covPosRelEntropy","editDistMode","nSoftClipAvg","contigId",
            "contigL","nReads","coverageSd","coverageBp","coverageP",
            "coveragePRatio","coverageCv","coverageEvennessScore","editDistAvg",
            "editDistAvgDecay","editDistDecayEnd","readLAvg","mqAvg","ani",
            "BayesianDMax","BayesianDMaxSd","BayesianZ","clonality"))%>%
  mutate(genelength=end-start,.after=end)%>%rename(gene_covP=avgcov)%>%
  mutate(gene_covP_ratio=normavgcov,.after=gene_covP)%>%select(!c("normavgcov"))
write.csv(gene_cov_anc,"~/PhD/Submission folder/Mutans paper folder/supp_table1.csv",row.names=FALSE,quote=TRUE)


#Supplemental Table 2
#this is the file loaded as jackie_dat in this script



####metadata columns####
#settheme
four_legend_theme<-  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
                         legend.key.height = unit(0.4, 'cm'), #change legend key height
                         legend.key.width = unit(0.4, 'cm'), #change legend key width
                         legend.title = element_text(size=10), #change legend title font size
                         legend.text = element_text(size=10), #change legend text font size
                         axis.text.y = element_blank(),axis.text.x = element_text(size = 9, angle = 45, hjust = 1))
three_legend_theme<-theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
                         legend.key.height = unit(0.4, 'cm'), #change legend key height
                         legend.key.width = unit(0.4, 'cm'), #change legend key width
                         legend.title = element_text(size=10), #change legend title font size
                         legend.text = element_text(size=10), #change legend text font size
                         axis.text.y = element_blank(),axis.text.x = element_text(size = 9, angle = 45, hjust = 1))

chrclus_col<-ggplot(strep_metadata,aes(y=NAME, x="Cluster",fill=group))+geom_tile()+
  scale_fill_manual(values=c("ancient_low_cov"="black","Cluster_A_ancient"="black","Cluster_A"="#ffce54",
                             "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568",
                             "Cluster_D_ancient"="black","Cluster_E"="#ed5564"),na.value="#FFFFFF",
                    name="Chromopainter \nCluster",na.translate = F)+
  four_legend_theme+
  xlab(NULL) + ylab(NULL)

chrclus_col_mags<-ggplot(strep_metadata%>%filter(NAME%in%phyl_withmags$tip.label),aes(y=NAME, x="Cluster",fill=group))+geom_tile()+
  scale_fill_manual(values=c("Cluster_A_ancient"="black","Cluster_A"="#ffce54","ancient_low_cov"="black",
                             "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568",
                             "Cluster_D_ancient"="black","Cluster_E"="#ed5564"),
                    name="Chromopainter \nCluster")+
  three_legend_theme+
  xlab(NULL) + ylab(NULL)

chrclus_col_pregub<-ggplot(strep_metadata%>%filter(NAME%in%phyl_pregub$tip.label),aes(y=NAME, x="Cluster",fill=group))+geom_tile()+
  scale_fill_manual(values=c("Cluster_A_ancient"="black","Cluster_A"="#ffce54",
                             "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568",
                             "Cluster_D_ancient"="black","Cluster_E"="#ed5564"),
                    name="Chromopainter \nCluster")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)


country_col<-ggplot(strep_metadata,aes(y=NAME, x="Origin",fill=factor(COUNTRY,levels = c("Brazil","USA","United Kingdom",
                                                                           "Denmark","China","Other"))))+geom_tile()+
  scale_fill_manual(values=c("Denmark"="#fc7d0b","United Kingdom"="#c85200","USA"="#5fa2ce",
                             "Brazil"="#1170aa","China"="#50A14F","Other"="#57606c"),na.value="#FFFFFF",
                    name="Origin",na.translate = F)+
  four_legend_theme+
  xlab(NULL) + ylab(NULL)


source_col<-ggplot(strep_metadata,aes(y=NAME, x="Source",fill=factor(SOURCE,levels = c("oral cavity","plaque","Ancient teeth",
                                                                         "blood","stool","Other"))))+geom_tile()+
  scale_fill_manual(values=c("Ancient teeth"="#f0bd27","plaque"="#ffda66","oral cavity"="#e39802",
                             "stool"="#51b364","blood"="#b60a1c","Other"="#57606c"),na.value="#FFFFFF",
                    name="Source",na.translate = F)+
  four_legend_theme+
  xlab(NULL) + ylab(NULL)

serotype_col<-ggplot(strep_metadata,aes(y=NAME, x="Serotype",fill=SEROTYPE))+geom_tile()+
  scale_fill_manual(values=c("C"="#1170aa","E"="#c85200","F"="#c8d0d9",
                             "K"="#ffbc79"),na.value="#FFFFFF",
                    name="Serotype",na.translate = F)+
  four_legend_theme+
  xlab(NULL) + ylab(NULL)

age_bar <- ggplot(strep_metadata, aes(y=NAME,DATE)) + geom_col()+
   theme_tree2() + theme(legend.position='none')+ggtitle("Sample Age") #coord_flip() +

## Genes
genes<-c("mubY","mubX","mubT","mubE","mubG","mubH","mubA","mubB","mubC","mubD","mubZ","mubM","mubI","mubJ","mubP","arcC","aguA","ptcA")
virfactors<-filter(gene_cov,gene%in%genes)

vir_bar <- ggplot(virfactors, aes(y=NAME, x=reorder(id,start))) + 
  ggtitle("virulence factors rare in ancient\nS. mutans (Reference mapping)") + 
  geom_tile(aes(fill=normavgcov))+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 8, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Observed/Expected\nbreadth of coverage",trans = 'reverse',option = "magma")+
  scale_x_discrete(label=virfactors$gene) 

#com operon coordinates
com_operon_mat<-filter(gene_cov,start>1500413)%>%filter(start<1508500)

com_operon_bar<-ggplot(com_operon_mat, aes(y=NAME, x=reorder(id,start),alpha=alphaparam)) + ggtitle("comG/comY operon\n(Reference mapping)") + 
  geom_tile(aes(fill=normavgcov))+
  scale_alpha_identity(guide = "none")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Observed/Expected\nbreadth of coverage",trans = 'reverse',option = "magma")+
  #labels=c("No", "Yes","No_ancient","Yes_ancient"),
  #values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=com_operon_mat$gene)


#virulence genes that have different coverage pattern in BAKTA and through reference mapping
genes<-c("mubY","mubX","gtfA","gtfB","gtfC","msmE")
selected_virfactors<-filter(gene_cov,gene%in%genes)%>%filter(NAME%in%phyl_withmags$tip.label)

sel_vir_bar <- ggplot(selected_virfactors, aes(y=NAME, x=reorder(id,start))) + ggtitle("selected virulence factors\n(Reference mapping)") + 
  #geom_tile(aes(fill=as.character(gene_coverage)))+
  geom_tile(aes(fill=normavgcov))+
  three_legend_theme +
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Observed/Expected\nbreadth of coverage",trans = 'reverse',option = "magma")+
  #scale_fill_manual(name="Gene coverage")+
  #labels=c("No", "Yes","No_ancient","Yes_ancient"),
  #values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=selected_virfactors$gene) 


sel_vir_BAKTA<-filter(BAKTA_all_fxd,Gene%in%genes)
sel_vir_bakta_bar <- ggplot(sel_vir_BAKTA, aes(y=NAME, x=factor(Gene,level=unique(selected_virfactors$gene)),fill=as.character(anc_gencov))) +
  ggtitle("selected virulence\nfactors (BAKTA)")+ geom_tile() + 
  three_legend_theme + 
  scale_fill_manual(name="Gene presence/absence",
                    labels=c("No", "Yes","No_ancient","Yes_ancient"),values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  xlab(NULL) + ylab(NULL)  


ancient_BAKTA<-filter(BAKTA_all_fxd,Gene%in%ancient_genes) #these are the genes which cannot be found through reference mapping, which have nevertheless been seen in ancient Streptococcus mutans
non_ref_BAKTA<-filter(BAKTA_all_fxd,Gene%in%non_ref_genes$Gene)#these are the genes which are not part of the reference genome in general
non_ref_coverage_mat <- ggplot(ancient_BAKTA, aes(y=NAME, x=reorder(Gene,freq),fill=as.character(anc_gencov))) + 
  ggtitle("Genes absent from reference genome GCF_009738105.1\nSorted by frequency of appearance")+
  geom_tile()  + theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
                       legend.key.size = unit(0.4, 'cm'), #change legend key size
                       legend.key.height = unit(0.4, 'cm'), #change legend key height
                       legend.key.width = unit(0.4, 'cm'), #change legend key width
                       legend.title = element_text(size=11), #change legend title font size
                       legend.text = element_text(size=10)) +#change legend text font size                                                       
  scale_fill_manual(name="Gene coverage",
                    labels=c("No", "Yes","No ancient","Yes ancient"),values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  xlab(NULL) + ylab(NULL) 


jackie_dat<-read_excel("~/PhD/misc_scripts/data/Supplementary table.xlsx")
cat_counts<-jackie_dat%>%select("Gene","Category")%>%distinct()%>%
  filter(!Category=="Translation")%>%filter(!Category=="Function Unknown")%>%filter(!Category=="nuclease")%>%#filter(!Gene%in%coregenes)%>%
  group_by(Category)%>%count(Category)%>%arrange(n)%>%tail(n=11)#%>%select(Category)
jackie_dat$simpcat <-ifelse(jackie_dat$Category %in% cat_counts$Category,jackie_dat$Category,"Other")
jackie_bacteriocin<-jackie_dat%>%filter(Category=="Bacteriocin protein/transcriptional regulator")
bacteriocin_BAKTA<-filter(BAKTA_all_fxd,Gene%in%jackie_bacteriocin$Gene)

jackie_bacteriocin_mat <- ggplot(bacteriocin_BAKTA, aes(y=NAME, x=reorder(Gene,freq),fill=as.character(anc_gencov))) + 
  ggtitle("Bacteriocins absent from reference genome GCF_009738105.1\nSorted by frequency of appearance")+
  geom_tile()  + theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
                       legend.key.size = unit(0.4, 'cm'), #change legend key size
                       legend.key.height = unit(0.4, 'cm'), #change legend key height
                       legend.key.width = unit(0.4, 'cm'), #change legend key width
                       legend.title = element_text(size=11), #change legend title font size
                       legend.text = element_text(size=10)) +#change legend text font size                                                                        
  scale_fill_manual(name="Gene coverage",
                    labels=c("No", "Yes","No_ancient","Yes_ancient"),values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  xlab(NULL) + ylab(NULL)  





#### Figures #### 

#Figure 1
plot_world<-ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    colour = "black", fill = "lightgray", size = 0.1,alpha=0.99,#width=3,height=1
  ) +  theme_void()+geom_point(
    data = strep_table,#position=position_jitter(width=1,height=1,seed=1),
    aes(longitude, latitude,color=ageAverage,label=sampleId,
        shape=factor(age_period, 
                     levels = c("Historical","Iron Age", "Bronze Age","Neolithic","Mesolithic"))),
    alpha = 0.99,
  ) +ylim(15,61) + xlim(-15,110)+#xlim(-50,165) #ylim(41, 73)  #+ ylim(37, 73) + xlim(-50,165) 
  scale_shape_manual(values = c(15,16,17,18,19))+
  scale_color_viridis(option = "viridis")+
  labs(title="Mutans sample overview",y="latitude",x="longitude")+
  geom_text_repel(data = strep_table,aes(longitude, latitude,label=sampleId),size=3,
                  min.segment.length = unit(0.2, 'lines'))+#,position=position_jitter(width=1,height=1,seed=1))+
  guides(col=guide_colourbar("Sample age"),shape=guide_legend("Archaeological\ncontext"))

plot_mutans <-ggplot(strep_table, aes(y=coveragePExp, x=coverageAvg, 
                                      shape=factor(age_period, 
                                                   levels = c("Historical","Iron Age", "Bronze Age","Neolithic","Mesolithic")))) +
  geom_point(aes(color=ageAverage),size=2) + 
  geom_vline(xintercept=7, color = "springgreen4",linetype="dashed",)+  scale_x_log10() + #scale_y_log10() + 
  annotate("text", x=25, y=0.5, label="assembly viable", angle=0)+
  scale_color_viridis(option = "viridis")+
  scale_shape_manual(values = c(15,16,17,18,19))+
  #geom_abline(color = "navyblue",linetype="dashed",yintercept=0.1) +
  #scale_color_manual(values=c('#999999','#E69F00'))+
  #geom_text(label=File,size=3, nudge_x = 0.09,nudge_y = -0.04) +
  #geom_label(label=File,position=position_jitter(height=0.2),size=3) +
  labs(title="Reference coverage",x="depth",y="breadth")+
  annotation_logticks(sides="b")+
  geom_text_repel(data = strep_table,aes(label=sampleId),size=3,min.segment.length = unit(0.2, 'lines')) +
  theme(legend.position='none')

plot_world / plot_mutans+ plot_layout(heights = c(4,2),guides = 'collect')+plot_annotation(tag_levels="A")


#Figure 2
#cmd command
#cd /home/vincentven/PhD/misc_scripts/data/PCA_asm5_3/
#Rscript plot_pca.R full.strep_mutans_asm5_4.eigenvec full metadata_chromopainter3.tsv ancientlist full.strep_mutans_asm5_4.pca.plot.pdf

#Figure 3
age_bar%>%insert_left(source_col,width=0.1)%>%insert_left(serotype_col,width=0.1)%>%insert_left(country_col,width=0.1)%>%
  insert_left(chrclus_col,width=0.1)%>%insert_left(tree_final,width=0.45)

#Figure 4
#beg martin

#Figure 5
vir_bar%>%insert_left(chrclus_col,width=0.1)%>%insert_left(tree_final,width=0.4)


#Figure S1
#bespoke picture

#Figure S2
#beg martin sikora

#Figure S3
tree_pregub_legend<-ggtree(midpoint.root(phyl_pregub), aes(color=as.numeric(label))) + 
  scale_color_viridis(name="bootstrap support")+
  ggtitle("RaXML-NG: core genome",subtitle="With midpoint rooting")
chrclus_col_pregub%>%insert_left(edit_phylo_branch(tree_pregub_legend,anc_phylo),width=7)#%>%insert_right(p5)

#Figure S4
gubbins_metadata<-strep_metadata%>%filter(NAME%in%phyl_gub$tip.label)%>%select("NAME","group","COUNTRY","SEROTYPE","SOURCE")
#arrange correctly for r script format
colnames(gubbins_metadata)<-c("id","Cluster","Origin","Serotype","Source")
write.csv(gubbins_metadata,"~/PhD/misc_scripts/data/gubbins_plot/metadata3.csv",row.names=FALSE,quote=TRUE)
#cmd command to create plot
#"Rscript /home/vincentven/PhD/misc_scripts/gubbins_plots.R --tree /home/vincentven/PhD/misc_scripts/data/gubbins_plot/strep_mutans_asm5_4.final_tree.tre \
#--rec /home/vincentven/PhD/misc_scripts/data/gubbins_plot/strep_mutans_asm5_4.recombination_predictions.gff --output plot2.pdf \
#--meta /home/vincentven/PhD/misc_scripts/data/gubbins_plot/metadata3.csv \
#--meta-width 0.1 --output-height 10 --legend-height 0.35 --meta-label-size 4 --heatmap-height 0.085"

#Figure S5
####Bootstrap comparisons####

tree_pregub+tree_gub+tree_postgub+plot_annotation(title="Bootstrap values for internal nodes for trees mutans",subtitle="Trees not rooted, ancient samples highlighted")  

#Figure S6
db_meta <- read.delim("~/PhD/misc_scripts/data/kraken_db_metadata",header=FALSE,col.name=c("File","Species"))
fastani <- read.delim("~/PhD/misc_scripts/data/fani_GCFtax.out", header=FALSE, col.names = c("File","Alignment","taxFile","taxAlignment","ANI"))
#File preprocessing of taxonomic info linked to fastani output for species
fastani <- transform(fastani,ANI=1-(ANI/100)) #transform to 1-ANI
faninames <-fastani[ , which(names(fastani) %in% c("File","taxFile"))] #take out the taxonomy info
colnames(faninames) <- c("File","Tax")
taxinfo <- faninames[!duplicated(faninames$File),] #remove duplicate names
taxinfo <- merge(taxinfo,db_meta) #get species info
row.names(taxinfo)<-taxinfo$File  
taxinfo <- taxinfo[,c("Species"), drop = F]
taxinfo[,c("Species")] = as.character(taxinfo$Species)

c21 <- c(
  "black","dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
   "khaki2",
  "maroon", "orchid1", "deeppink1",  "steelblue4",
  "darkturquoise", "green1",  "yellow3",
   "brown"
)

#convert fastani output to distance matrix and do hierachical clustering
ani_long <- fastani[,which(names(fastani)%in% c("File","Alignment","ANI"))]
ani_wide <- spread(ani_long,Alignment,ANI,fill=0.25) #fill missing values as having ANI=0.75
row.names(ani_wide)<-ani_wide$File #convert to matrix
ani_wide<-ani_wide[,-1]
ani_wide<-as.matrix(ani_wide)

repseqs<-dissrep(ani_wide,criterion="density",coverage=0.99)#extract representatives using neighbourhood distance with 0.99 desired coverage level
repindex<-repseqs[1:length(repseqs)]

#plot the representative subset
#use version of taxinfo with no infrequent taxa defined as other
taxinfo.freq = as.data.frame(table(taxinfo)) #find the frequency of different tax names
taxinfo.freq[,c("taxinfo")] = as.character(taxinfo.freq$taxinfo)
other = taxinfo.freq$Species[taxinfo.freq$Freq < 5] #define infrequent taxa as other for better overview
taxinfo$Species = ifelse(taxinfo$Species %in% other, "other", taxinfo$Species)
ani_rep<-ani_wide[repindex,repindex]
clustered = hclust(as.dist(ani_rep), method = "average") #compute hierarchical clustering

rpdb <- ggtree(clustered, layout = "circular") + 
  ggtitle("Hierarchical clustering of 115 representative subset of \ncustom kraken database for Streptococcus",subtitle="S. mutans marked with red diamond")+
  geom_point2(aes(subset=(node==8)), shape=23,size=5, fill='red')
rpdbring <- gheatmap(rpdb,taxinfo, width = 0.08,offset=0.01,colnames=F) + 
  theme(legend.key.size = unit(0.3, 'cm'))+scale_fill_manual(name="Species Name",values=c21)     
rpdbring





#Figure S7
mag_stats_1kb <- read_tsv("~/PhD/misc_scripts/data/mutans_bins_1kb.txt")  %>% 
  #mutate('100-red' = 100 - percent_redundancy) %>% 
  pivot_longer(cols = c("total_length","num_contigs", "N50", "GC_content", "percent_completion", "percent_redundancy")) %>% 
  mutate(name = factor(name, levels = c("total_length","num_contigs", "N50", "percent_completion", "percent_redundancy", "GC_content"))) 
mag_stats_plot_1 <- ggplot(mag_stats_1kb,aes(x = Sample,
                                             y = value)) +
  facet_grid(name ~ ., scale="free",space="free_x") +
  geom_point(aes(fill=Sample, color = `Assembly-strategy`, shape = `Binning-strategy`), 
             position = (position = position_jitterdodge(jitter.width = 0.1, 
                                                         dodge.width =  0.8))) +
  ylab("assembly & binning strategy") + ggtitle("Mutans bin statistics using 1kb contig threshold") +
  theme_bw()

#Figure S8
mag_stats_2_5kb <- read_tsv("~/PhD/misc_scripts/data/mutans_bins_2_5kb.txt")  %>% 
  #mutate('100-red' = 100 - percent_redundancy) %>% 
  pivot_longer(cols = c("total_length","num_contigs", "N50", "GC_content", "percent_completion", "percent_redundancy")) %>% 
  mutate(name = factor(name, levels = c("total_length","num_contigs", "N50", "percent_completion", "percent_redundancy", "GC_content"))) 
mag_stats_plot_2 <- ggplot(mag_stats_2_5kb,aes(x = Sample,
                                               y = value)) +
  facet_grid(name ~ ., scale="free",space="free_x") +
  geom_point(aes(fill=Sample, color = `Assembly-strategy`, shape = `Binning-strategy`), 
             position = (position = position_jitterdodge(jitter.width = 0.1, 
                                                         dodge.width =  0.8))) +
  ylab("assembly & binning strategy") + ggtitle("Mutans bin statistics using 1kb contig threshold") +
  theme_bw()




#Figure S9
modules <- read.delim("~/PhD/misc_scripts/data/pangenome/Mutans_enriched_modules.txt")%>%filter(p_MODERN>0.8)%>%filter(p_MAG<0.8)%>%rename(module=accession) #get metabolic modules 
mod_info<-read.delim("~/PhD/misc_scripts/data/pangenome/modules_info.txt")%>%filter(module%in%modules$module)
mod_cov <- read.delim("~/PhD/misc_scripts/data/pangenome/Mutans_metabolism_mat-module_stepwise_completeness-MATRIX.txt")%>%
  filter(module%in%modules$module)#how complete are these modules /how many of the constituent genes are present
colnames(mod_cov)<-gsub("MUTREF","Modern Mutans Assembly ",colnames(mod_cov))
colnames(mod_cov)<-gsub("_"," ",colnames(mod_cov))
#convert dataframe to matrix
module_matrix<- mod_cov[,-1]
row.names(module_matrix)<-mod_cov[,1]
module_matrix<-as.matrix(module_matrix)
module_matrix <- replace(module_matrix,is.na(module_matrix),0)

heatmap.2(t(module_matrix), 
          col = heat.colors(20),    # Choose the color palette (you can change this)
          scale = "column",         # Scale the values by column
          hclustfun = function(x) hclust(x, method="average"),
          Rowv = TRUE, Colv = TRUE,     # Turn off row and column clustering
          margins = c(5, 10),       # Set margins for row and column labels
          main = "Coverage for metabolic modules (%)", # Set the title of the plot
          xlab = "Selected Metabolic modules", ylab = NA,
          colRow = NA, colCol = "black", # Turn off color for row dendrogram and use black for column dendrogram
          trace = "none",           # Turn off trace
          key = TRUE, keysize = 1.0, cexRow = 0.8, cexCol = 0.8, # Adjust key parameters
          density.info = "none",    # Turn off density plot
          dendrogram = "both")    # Show only column dendrogram
#ColSideColors = column_metadata)  # Add metadata sidebar for columns

#Figure S10
com_operon_bar%>%insert_left(chrclus_col,width=0.1)%>%insert_left(tree_final,width=0.4)

#Figure S11
sel_vir_bar%>%insert_left(sel_vir_bakta_bar,width=0.8)%>%insert_left(chrclus_col_mags,width=0.2)%>%insert_left(tree_withmags)

#Figure S12
non_ref_coverage_mat%>%insert_left(chrclus_col_mags,width=0.05)%>%insert_left(tree_withmags,width=0.4)

#Figure S13
jackie_bacteriocin_mat%>%insert_left(chrclus_col_mags,width=0.05)%>%insert_left(tree_withmags,width=0.4)




####Supplemental Section Figures
#Figure SA1
#gunc plot function from: https://grp-bork.embl-community.io/gunc/usage.html 

#Figure SB1
#molecular clock from bactdating
set.seed(16)
tiplabels<-as.data.frame(phyl_rec$tip.label)%>%rename("NAME"="phyl_rec$tip.label")%>%left_join(strep_metadata)
rtt<-initRoot(phyl_rec,tiplabels$YEAR,useRec=T)
rtt_pp<-roottotip(rtt,tiplabels$YEAR)
#Figure SB2
bd_gub<-bactdate(phyl_rec,tiplabels$YEAR,useRec=T,showProgress = TRUE)
bd_gubcon<-bactdate(phyl_rec,showProgress = TRUE,useRec=T, rep(2022,length(phyl_rec$tip.label)))
modelcompare(bd_gub,bd_gubcon)
plot(bd_gub,'trace')
#plot(bd_gub,'treeCI')

#Table SB1 and SB2
#On the subtree of the recombination tree
tree_cluster<-extract.clade(phyl_rec,node="Node_256")
tiplabels<-as.data.frame(tree_cluster$tip.label)%>%rename("NAME"="tree_cluster$tip.label")%>%left_join(strep_metadata)
bd_cluster<-bactdate(tree_cluster,tiplabels$YEAR,showProgress = TRUE,useRec=T)
bd_control<-bactdate(tree_cluster,showProgress = TRUE,useRec=T, rep(2022,length(tree_cluster$tip.label)))
modelcompare(bd_cluster,bd_control)
rtt<-initRoot(tree_cluster,tiplabels$YEAR,useRec=T)
rtt_pp<-roottotip(rtt,tiplabels$YEAR)




#Figure SC1
#Print out the sCF tree from IQ tree
tree_postgub_cf<-ggtree(midpoint.root(phyl_postgub_cf), aes(color=as.numeric(label))) + 
  scale_color_viridis(name="sCF")+
  #scale_color_gradientn( #old color gradient which colors anything below 33 the same
  #  name = "sCF",
  #  colours=as.numeric(label),
  #colours = c(rep("black",33),tail(viridis(100),n=66)),
  #  limits = c(0, 100)
  #)+
  ggtitle("IQ-Tree: masked \ncore genome")
tree_postgub_root<-ggtree(midpoint.root(phyl_postgub), aes(color=as.numeric(label))) + 
  scale_color_viridis(name="bootstrap support")+
  ggtitle("RaXML-NG: masked core genome")

tree_postgub_cf+tree_postgub_root

#Figure SD1
#picture credit to Nicolas and Signe


#### Additional File 2
#Description: Phylogenetic placement of a low coverage ancient genomes using EPA-NG on the phylogenetic tree inferred with the high-coverage modern and ancient samples.
#Visualize likelihood weights for placements
lowcov<-c("CGG100272","CGG100534","CGG_2_18521","NEO938","RISE373","NEO31","DA112","DA246","DA223","VK146","DA195","NEO65","RISE569",
          "RISE540","VK35","RISE349","DA104","RISE413","NEO170","NEO160","VK441","S_troglodytae")
#NEO065 og VK441 excluded, bringing list to 19.
#lowcov<-c("RISE540","NEO65","DA246")
#lowcov<-c("S_troglodytae")
#check different placements on the tree
#pdf(paste0("~/PhD/datafiles/placement_4/placementunrooted.pdf"))#, width=8, height=25)
#pdf(paste0("~/PhD/datafiles/placement_4/troglodytae.pdf"))#, width=8, height=25)
for(i in 1:length(lowcov)){
  sample<-lowcov[i]
  tree <- read.tree(paste("~/PhD/misc_scripts/data/placement_5/",sample,
                          ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.newick",sep=""))
  d <- read.jplace(paste("~/PhD/misc_scripts/data/placement_5/",sample,
                         ".strep_mutans_asm5_4",".GCF_009738105.1_ASM973810v1.epa_ng.jplace",sep=""))
  #tree_mutans <- drop.tip(tree, "GCF_002355215N1")
  tree_mutans<-tree
  #d<-drop.tip(d, "GCF_002355215N1") #remove the troglodytae root
  missingtaxa<-setdiff(strep_metadata$NAME,tree_mutans$tip.label)
  metadata<-strep_metadata %>% filter(!NAME%in%missingtaxa)
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
  missingtaxa<-setdiff(strep_metadata$NAME,c(d@placements$name,d@phylo$tip.label))
  metadata<-strep_metadata %>% filter(!NAME%in%missingtaxa)
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
#dev.off()

#Additional File 3
#outputted by R script in reference mapping pipeline
#Description: plot showing percentage of reference covered by mapping an ancient sample against each species of genus Streptococcus included in this study, coloured by ANI. The references are sorted by phylogeny, with species name on the left and the NCBI RefSeq assembly name on the right. 

#Additional File 4
#Made during anvio mag construction 
#Description: Anvi’o visualization of our ancient MAGs.
#Each layer in these circular plots represents a statistics for each contig. The contigs are clustered together in a tree by 4-mer frequencies.
#The green layer represents the GC content.
#The black layers describe the number of reads from each library in the sample which can be mapped to the contig.
#The outermost black layer (damage) describes the average local damage for each read used to assemble this contig as calculated through metaDMG
#The familiarity layer is a binary layer which describes whether a gene cluster annotated within the contig can be found in any of the other 6 MAGs
#through a pangenome analysis before manual refinement is started.
#Any contig marked as contaminant is mapped to the NCBI database of genes, to confirm that these cannot be mapped to S. mutans.



###trash
#On the final tree with root-to-tip correlation rooting
tiplabels<-as.data.frame(phyl_gub$tip.label)%>%rename("NAME"="phyl_gub$tip.label")%>%left_join(strep_metadata)
rtt<-initRoot(phyl_gub,tiplabels$YEAR)
rtt_pp<-roottotip(rtt,tiplabels$YEAR)
#compared with bactdatings rooting
bd_pp <- bactdate(phyl_gub,tiplabels$YEAR,showProgress = TRUE,updateRoot= T)
plot(bd_pp,'trace')
plot(bd_pp,'treeCI')
#On the final tree rooted with S troglodytae outgroup rooting
trogroot<-root(phyl_postgub_trog,"Troglodytae",resolve.root = TRUE)%>%drop.tip("Troglodytae")
tiplabels<-as.data.frame(trogroot$tip.label)%>%rename("NAME"="trogroot$tip.label")%>%left_join(strep_metadata)
bd_trogroot<-bactdate(trogroot,tiplabels$YEAR,showProgress = TRUE,updateRoot= F)
#compare the two models
modelcompare(bd_pp,bd_trogroot)

#On the recombination tree from gubbins gubbins tree
bd_gub<-bactdate(phyl_rec,tiplabels$YEAR,useRec=T,showProgress = TRUE,nbIts = 500000)
#bd_control<-bactdate(phyl_rec,rep(2022,length(phyl_rec$tip.label)),useRec=T,showProgress = TRUE,nbIts = 500000)
#modelcompare(bd_gub,bd_control)