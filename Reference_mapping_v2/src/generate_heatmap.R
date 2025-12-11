library(gplots) 
library(ggplot2)
library(viridis)
library(tidyverse)




################## ################## ##################
plot_style=list(
  xlab("Gene"),
  ylab("Sample"),
  theme(strip.background = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=0),
        strip.text.y.left = element_text(angle = 90,hjust=1),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")))
################## ################## ##################
norm_cov<-function(x,avgcovall){
  genepresence<-ifelse((x/avgcovall)>0.1,1,0)
  return(genepresence)
}
findavg<-function(x,mutans_samples){
  return(mutans_samples[which(mutans_samples$sample==x),]$coverageAvg)
}


#sample_info <- read.delim("~/PhD/Strep_mapping/data/sumdat/topmap_allinfo.tsv", header=FALSE,col.names =c("sampleId","assemblyId","coverageP",
#                                                                                                          "ani","coverageAvg","coverageEvenessScore",
#                                                                                                          "heterogenity","classification")
#)
sample_info <- read.delim("~/PhD/datafiles/summary_statistics.tsv", header=FALSE,col.names =c("sampleId","assemblyId","coverageP",
                                                                                                          "ani","coverageAvg","coverageEvenessScore",
                                                                                                          "unknown","heterogenity")#,"classification")
)
sample_meta <-read.delim("PhD/community_analysis_09-22/tax_pipeline/20230101.diseases.sampleInfo.tsv") %>% select("sampleId","clusterLabel","country","latitude","longitude","ageAverage","region")
sample_meta<-rbind(sample_meta,c("CGG100272",NA,"Iran",NA,NA,as.double(2480.0),"CentralAsia"))
sample_meta<-rbind(sample_meta,c("CGG101233",NA,"Austria",NA,NA,as.double(2675),"CentralEasternEurope"))
sample_meta$ageAverage<-as.double(sample_meta$ageAverage)
testcov3<-read.delim("PhD/datafiles/cov_all_highcov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                     col.names=c("id","gene","start","end","avgcov","min","max","sampleId","covP"))

complete_sample<-left_join(sample_info,sample_meta,by="sampleId") #get the metadata on to my sample classification table
mutans_dominant<-which(complete_sample$assemblyId=="GCF_000007465.2") #find the samples that were mapped to mutans
bad_samples<-which(complete_sample$coverageAvg<1) #filter out those samples where i did not classify it as worth looking into
mutans<-setdiff(mutans_dominant,bad_samples) #find out which mutans samples are bad samples
mutans_samples<-complete_sample[mutans,] #and save those that pass the filter
genes<-c("dltC","ffh","gtfA","gtfB","gtfC","cnm","spaP","spaA","spaB","spaC","cnaB","cbpA","cbm")  #define the genes i want to look for






mutans_genes<-testcov3[which(testcov3$sampleId%in%mutans_samples$sampleId),] #find out which samples are good mutans samples
genes_of_interest<-mutans_genes[which(mutans_genes$gene%in%genes),]
pathogen<-left_join(mutans_genes,sample_meta,by="sampleId")
#pathogen<-left_join(genes_of_interest,sample_meta,by="sampleId") #filter out the genes of interest and add sample metadata to the table
mutans_samples$gene<-"Age" #define a gene as age, so i can add this metadata to the heatmap as a column
avgcovall<-sapply(pathogen$sampleId,findavg,mutans_samples)  #find mean coverage across reference genome for normalization
pathogen$normavgcov<-as.factor(norm_cov(pathogen$avgcov,avgcovall))#ascertain based on average gene coverage if the gene is therr or not
pathogen2<-bind_rows(pathogen,mutans_samples) #bind age column to table for metadata in heatmap
pathogen4<-pathogen2[grep("mub",pathogen2$gene),]
pathogen5<-pathogen2[pathogen2$gene%in%genes,]
pathogen3<-bind_rows(pathogen5,pathogen4)
pdf("virulence_heatmap3.pdf", width = 8, height = 10)
ggplot(filter(pathogen,gene!='Age'),aes(x=sampleId,y=gene))+ #mathematical racism
  geom_tile(color = "white",aes(fill=normavgcov))+ #define heatmap
  facet_grid(~region,scale='free',space='free',switch='y')+ #discriminate people by birthplace
  #scale_fill_manual()+
  plot_style + #do it with style!
  geom_point(data = filter(pathogen2,gene=='Age'),  #add the age of the people
           aes(color = ageAverage), size = 8)+
  #scale_fill_viridis(option="magma")+
  scale_fill_discrete(name="Gene coverage",
                      labels=c("No", "Yes"))
dev.off()


#adding metadata to sample files
units <- read.delim("~/PhD/community_analysis_09-22/tax_pipeline/units.tsv")
sample_all<-left_join(complete_sample,units,by="sampleId")
write.table(sample_all,file="~/PhD/Strep_mapping/data/sumdat/topmap_meta.tsv",sep = "\t") #write it out to a file





############################Doing it on all mutans samples available#######################################
library(tidyverse)
library(ggtree)
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
ancientmeta<-rbind(ancientmeta,c("CGG101233","Austria",2675))
ancientmeta$DATE<-as.double(ancientmeta$DATE)
ancientmeta$SOURCE<-"Archaeology"
combined_meta$DATE<-2022-as.double(combined_meta$DATE) #correct the dating with the other metadata
combined_final<-left_join(combined_meta,ancientmeta,by="NAME")%>%mutate(COUNTRY = coalesce(COUNTRY.x, COUNTRY.y)) %>% 
  mutate(DATE = coalesce(DATE.x, DATE.y)) %>%mutate(SOURCE = coalesce(SOURCE.x, SOURCE.y))%>%
  select(NAME, COUNTRY, DATE,TYPE,COV,SEROTYPE,SOURCE) %>% mutate_if(is.character, ~na_if(., '')) #combine all the sample metadata into one file
#country_counts<-combined_final%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
#combined_final$COUNTRY <-ifelse(combined_final$COUNTRY %in% country_counts$COUNTRY,combined_final$COUNTRY,"Other")
scov<-read.delim("PhD/datafiles/GCF_009738105.1_ASM973810v1.strep_mutans_v3_all_snps.concat_gene_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                     col.names=c("id","gene","start","end","avgcov","sampleId"))
scds<-read.delim("PhD/datafiles/GCF_009738105.1_ASM973810v1.strep_mutans_v3_all_snps.concat_cds_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("name","desc","start","end","avgcov","sampleId"))
combined_cov<-full_join(scov,scds,by=c("start","end","avgcov","sampleId"))
scov<-combined_cov
#complete_sample<-left_join(sample_info,sample_meta,by="sampleId") #get the metadata on to my sample classification table
genes<-c("dltC","ffh","cnm","cnaB","cbpA","cbm")  #define the genes i want to look for
#,"gtfA","gtfB","gtfC"
#"spaP","spaA","spaB","spaC"
genes_of_interest<-scov[which(scov$gene%in%genes),]
combined_final$sampleId<-combined_final$NAME
pathogen<-left_join(scov,combined_final,by="sampleId")
#pathogen<-left_join(genes_of_interest,sample_meta,by="sampleId") #filter out the genes of interest and add sample metadata to the table
mutans_samples$gene<-"Age" #define a gene as age, so i can add this metadata to the heatmap as a column
pathogen$normavgcov<-as.factor(ifelse(pathogen$avgcov>0.1,1,0))#ascertain based on average gene coverage if the gene is therr or not
pathogen2<-pathogen
pathogen2<-bind_rows(pathogen,mutans_samples) #bind age column to table for metadata in heatmap
pathogen3<-pathogen2[grep("mub",pathogen2$gene),]
#pathogen6<-pathogen2[grep("spa",pathogen2$gene),]
#pathogen7<-pathogen2[grep("gtf",pathogen2$gene),]
#pathogen5<-pathogen2[pathogen2$gene%in%genes,]
#pathogen3<-bind_rows(pathogen5,pathogen4,pathogen6,pathogen7)

ggplot(filter(pathogen3,gene!='Age'),aes(x=sampleId,y=reorder(id,start),fill=normavgcov,label=gene))+ #mathematical racism
  geom_tile(color = "white")+ #define heatmap
  
  #facet_grid(~TYPE,scale='free',space='free',switch='y')+ #discriminate people by birthplace
  #scale_fill_manual()+
  plot_style+ #do it with style!
  geom_point(data = filter(pathogen2,gene=='Age'),  #add the age of the people
             aes(color = ageAverage), size = 8)+
  #scale_fill_viridis(option="magma")+
  scale_fill_discrete(name="Gene coverage",
                      labels=c("No", "Yes"))+
  guides(fill = guide_colourbar(label = FALSE,
                                ticks = FALSE))



###Get the tree in a pretty format with placements 
tree <- read.tree(paste("~/PhD/datafiles/placement_3/placementall/","merge.newick",sep=""))
library("treeio")
tree_mutans <- drop.tip(tree, "GCF_002355215N1")
#tree_mutans<-rooted
missingtaxa<-setdiff(combined_final$NAME,tree_mutans$tip.label)
metadata<-combined_final %>% filter(!NAME%in%missingtaxa)
country_counts<-metadata%>%group_by(COUNTRY)%>%count(COUNTRY)%>%arrange(n)%>%tail(n=6)%>%select(COUNTRY)
metadata$COUNTRY <-ifelse(metadata$COUNTRY %in% country_counts$COUNTRY,metadata$COUNTRY,"Other")
country_counts<-metadata%>%group_by(SOURCE)%>%count(SOURCE)%>%arrange(n)%>%tail(n=6)%>%select(SOURCE)
metadata$SOURCE <-ifelse(metadata$SOURCE %in% country_counts$SOURCE,metadata$SOURCE,"Other")
tree_mutans_full<-tree_mutans%>%left_join(metadata,by=c("label"="NAME"))
p<-ggtree(tree_mutans_full) + 
  #geom_label(aes(x=branch),fill='lightgreen') +
  geom_rootedge(rootedge = 0.0005)+geom_tippoint(
    mapping = aes(
      subset=TYPE!="REF",
      #shape = COUNTRY, 
      color = DATE,
      size=4,offset=0.001)
    )
   #+geom_tiplab(aes(subset=TYPE!="REF"),size = 4,offset=.001)
  #ggtitle("RaxML-ng phylogenetic tree including placement samples",subtitle = "Nodes with only modern samples: 441 (Green circle), 361 (Red square)")
p
missing_from_ancient<-pathogen%>%filter(TYPE!="REF")%>%group_by(start)%>%
  summarise(freq=mean(as.double(normavgcov)-1))%>%filter(freq<0.5)%>%select(start)%>%as.list
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

metmat<-metadata%>%select("NAME","COUNTRY")
rownames<-metmat$NAME
metmat<-select(metmat,-NAME)%>%as.matrix
row.names(metmat)<-rownames

hey<- gheatmap(p,genemat, offset=0, width=5, 
         colnames=TRUE, legend_title="Gene coverage",
         colnames_angle=90,colnames_offset_y=-1,
         custom_column_labels=genenames$desc,hjust=1) + 
  scale_x_ggtree() + ggtitle("hello")+
  scale_y_continuous(expand=c(0, 300))
hey
#save with width of 2000 and height 1500
library(ggnewscale)
hey2<-hey + new_scale_fill()
gheatmap(hey2,metmat,offset=0.0175,width=.1)

library(aplot)

gheatmap(hey,genemat)
p+geom_tiplab(aes(subset=TYPE!="REF"),size = 4,offset=.001)

p%>%insert_left(p)%>%insert_right(p2)

ggplot(p)

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
p3<-ggplot(metadata,aes(y=NAME, x="Country",fill=COUNTRY))+geom_tile()+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)
p3
p4<-ggplot(metadata,aes(y=NAME, x="Source",fill=SOURCE))+geom_tile()+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 45, hjust = 1))+
  xlab(NULL) + ylab(NULL)
p4
p2%>%insert_right(p1,width=0.6)%>%insert_left(p3,width=0.1)%>%insert_left(p4,width=0.1)%>%insert_left(g,width=1)
test<-pathogen3%>%filter(sampleId=="NEO105")%>%select(gene)%>%as.list


melt(cbind(sample=rownames(qualityScores), qualityScores))






##############################Porphyromonas Gingivalis######################################################

sample_info <- read.delim("~/PhD/datafiles/porphyromonas/summary.tsv", header=FALSE,col.names =c("sampleId","assemblyId","coverageP",
                                                                                              "ani","coverageAvg","coverageEvenessScore",
                                                                                              "unknown","heterogenity")#,"classification")
)
sample_meta <-read.delim("PhD/community_analysis_09-22/tax_pipeline/20230101.diseases.sampleInfo.tsv") %>% select("sampleId","clusterLabel","country","latitude","longitude","ageAverage","region")
sample_meta<-rbind(sample_meta,c("CGG100272",NA,"Iran",NA,NA,as.double(2480.0),"CentralAsia"))
sample_meta<-rbind(sample_meta,c("CGG101233",NA,"Austria",NA,NA,as.double(2675),"CentralEasternEurope"))
sample_meta$ageAverage<-as.double(sample_meta$ageAverage)
testcov3<-read.delim("~/PhD/datafiles/porphyromonas/cov_porphy.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                     col.names=c("id","gene","start","end","avgcov","min","max","sampleId","covP"))

complete_sample<-left_join(sample_info,sample_meta,by="sampleId") #get the metadata on to my sample classification table
porphy_dominant<-which(complete_sample$assemblyId=="GCF_002754155.1") #find the samples that were mapped to mutans
bad_samples<-which(complete_sample$coverageAvg<1) #filter out those samples where i did not classify it as worth looking into
porphy<-setdiff(porphy_dominant,bad_samples) #find out which mutans samples are bad samples
porphy_samples<-complete_sample[porphy,] #and save those that pass the filter
genes<-c("dltC","ffh","gtfA","gtfB","gtfC","cnm","spaP","spaA","spaB","spaC","cnaB","cbpA","cbm")  #define the genes i want to look for






porphy_genes<-testcov3[which(testcov3$sampleId%in%porphy_samples$sampleId),] #find out which samples are good mutans samples
genes_of_interest<-porphy_genes[which(mutans_genes$gene%in%genes),]
pathogen<-left_join(porphy_genes,sample_meta,by="sampleId")
#pathogen<-left_join(genes_of_interest,sample_meta,by="sampleId") #filter out the genes of interest and add sample metadata to the table
porphy_samples$gene<-"Age" #define a gene as age, so i can add this metadata to the heatmap as a column
avgcovall<-sapply(pathogen$sampleId,findavg,porphy_samples)  #find mean coverage across reference genome for normalization
pathogen$normavgcov<-as.factor(norm_cov(pathogen$avgcov,avgcovall))#ascertain based on average gene coverage if the gene is therr or not
pathogen2<-bind_rows(pathogen,porphy_samples) #bind age column to table for metadata in heatmap
rare_genes<-pathogen%>%group_by(start)%>%
  summarise(freq=mean(as.double(normavgcov)-1))%>%filter(freq<0.9)%>%select(start)%>%as.list
pathogen4<-pathogen[pathogen$start%in%rare_genes$start,]
pathogen4<-pathogen2[grep("mub",pathogen2$gene),]
pathogen5<-pathogen2[pathogen2$gene%in%genes,]
pathogen3<-bind_rows(pathogen5,pathogen4)
#pdf("virulence_heatmap3.pdf", width = 8, height = 10)
ggplot(filter(pathogen4,gene!='Age'),aes(x=sampleId,y=id))+ #mathematical racism
  geom_tile(color = "white",aes(fill=normavgcov))+ #define heatmap
  facet_grid(~region,scale='free',space='free',switch='y')+ #discriminate people by birthplace
  #scale_fill_manual()+
  plot_style + #do it with style!
  geom_point(data = filter(pathogen2,gene=='Age'),  #add the age of the people
             aes(color = ageAverage), size = 8)+
  scale_y_discrete(labels=pathogen4$gene)+
  #scale_fill_viridis(option="magma")+
  scale_fill_discrete(name="Gene coverage",
                      labels=c("No", "Yes"))
#dev.off()



#################VIRULENCE FACTOR HUNT#####################################
CGG100272 <- read.delim2("~/Downloads/MAGs-ref_annotations/CGG100272-contigs.tsv", header=TRUE,skip=5)
CGG100272$NAME<-"CGG100272"
VK63 <- read.delim2("~/Downloads/MAGs-ref_annotations/VK63-contigs.tsv", header=TRUE,skip=5)
VK63$NAME<-"VK63"
NEO137 <- read.delim2("~/Downloads/MAGs-ref_annotations/NEO137-contigs.tsv", header=TRUE,skip=5)
NEO137$NAME<-"NEO137"
NEO105 <- read.delim2("~/Downloads/MAGs-ref_annotations/NEO105-contigs-filtered.tsv", header=TRUE,skip=5)
NEO105$NAME<-"NEO105"
NEO938 <- read.delim2("~/Downloads/MAGs-ref_annotations/NEO938-contigs.tsv", header=TRUE,skip=5)
NEO938$NAME<-"NEO938"
CGG100534<-read.delim2("~/Downloads/MAGs-ref_annotations/CGG100534-contigs.tsv", header=TRUE,skip=5)
CGG100534$NAME<-"CGG100534"
CGG101233<-read.delim2("~/Downloads/MAGs-ref_annotations/CGG101233-contigs.tsv", header=TRUE,skip=5)
CGG101233$NAME<-"CGG101233"
genes<-c("dltC","ffh","gtfA","gtfB","gtfC","gtfD","cnm","spaP","spaA","spaB","spaC","cnaB","cbpA","cbm","mutA","mutA'","mutB","mutC","mutD","fruA","fruB","manL",
         "levD","levE","levF","levG","eepA","gbpA","gbpB","gbpC","gbpD","levQ","srtA","brpA","fabM","ftf","spx","comX","comS","comR","ftsN","aguA","sfcA",
         "atpE","atpB","atpF","atpH","atpA","atpG","atpC","atpD","comEC","comEA","comFC","comFA","comA","comGA","comGB","comGC","comGD","comGE","comGF","comGG",
         "comP","comEB","comB","comF","scrA","scrB","scrAB","lacD","lacR","lacA","lacB","lacC","lacF","lacE","lacG","lacX","lacI",
         "levD","msmK","msmG","msmF","msmE","msmR","sppA","frwB","manX","agaB","manY","manZ","wapA", "pepO", "covR","covS","covRS","VicK","VicR","VicX","Smu63c","atlA","ccpA","ptsI")
genes<-c("mutA","mutA'","mutB","mutC","esaA","essA","essB","essC","bacA","macB","lanM")
         #"msmK","msmG","msmF","msmE","msmR","mubY","mubX","mubT","mubE","mubG","mubH","mubA","mubB","mubC","mubD","mubZ","mubM","mubI","mubJ","mubP","arcC","aguA","ptcA",
         #,"gtfA","gtfB","gtfC","gtfD")
genes<-c("mubY","mubX","gtfA","gtfB","gtfC","msmE","spaP")

#looking for genes that are missing
#smu630 or atlA: These should be there as well
genes<-c("covR","covS","covRS","VicK","VicR","VicX","gbpA","gbpB","gbpC","gbpD","atlA")
genes<-c("vicK","yycF","ompR")
genes<-c("mutA","mutA'","mutB","mutC","mutD","mutE","mutF","mutG")
#modern<-read.delim2("~/Downloads/MAGs-ref_annotations/GCF_009738105.1_ASM973810v1_genomic.tsv", header=TRUE,skip=5)
#setdiff(modern$Gene,CGG100272$Gene)
#combgenes<-rbind(CGG100272,modern)
#View(combgenes)





cds_modern <- read.delim("~/PhD/datafiles/cds_modern.tsv", 
                         col.names = c("sampleId","X.Sequence.Id","Type","Start","Stop","Strand","Locus.Tag","Gene","Product","DbXrefs"))
cds_modern$NAME<-sub("\\.","N",cds_modern$sampleId)
comb_genes<-cds_modern%>%select(!sampleId)%>%rbind(NEO938)%>%rbind(NEO137)%>%rbind(NEO105)%>%rbind(VK63)%>%rbind(CGG100534)%>%
  rbind(CGG101233)%>%rbind(CGG100272)#%>%rbind(c("GCF_002355215N1",NA,NA,NA,NA,NA,NA,NA,NA,NA))
comb_genes$gene_coverage<-as.integer(1)
comb_genes<-complete(comb_genes, NAME,Gene, fill = list(gene_coverage = as.integer(0)))%>%subset(Gene != "")%>%
  mutate_if(is.integer, ~replace(., is.na(.), 0))
virulence<-comb_genes%>%filter(Gene%in%genes)
unique(virulence$Gene)
setdiff(genes,unique(virulence$Gene))
#virulencetest<-comb_genes%>%filter(Gene%in%genes)



#Not many genes unique among the ancient
setdiff(CGG100272$Gene,cds_modern$Gene)
setdiff(NEO137$Gene,cds_modern$Gene)
setdiff(CGG101233$Gene,cds_modern$Gene)
setdiff(NEO105$Gene,cds_modern$Gene)
setdiff(CGG100534$Gene,cds_modern$Gene)
setdiff(NEO938$Gene,cds_modern$Gene)
setdiff(VK63$Gene,cds_modern$Gene)
#setdiff(cds_modern$Gene,CGG100272$Gene)

#good command for finding virulence genes
unique(comb_genes$NAME[grep("Listeria sRNA rli28",comb_genes$Product,ignore.case=TRUE)])#Glucosyltransferase-S
unique(comb_genes$Gene[grep("muk",comb_genes$Gene,ignore.case=TRUE)])#Glucosyltransferase-SI
unique(comb_genes$Gene[grep("muk",comb_genes$Product,ignore.case=TRUE)])
unique(comb_genes$Product[grep("muk",comb_genes$Gene,ignore.case=TRUE)])

genes_met<-left_join(comb_genes,combined_final,by="NAME")
in_ancient<-genes_met%>%filter(TYPE=="MAG")%>%group_by(Gene)%>%
  summarise(freq=mean(as.double(gene_coverage)))%>%filter(freq>0.1)%>%select(Gene)%>%as.list
non_ref_genes<-genes_met%>%filter(NAME=="GCF_009738105N1")%>%filter(gene_coverage==0)%>%select(Gene)
#non_core_genes<-genes_met%>%filter(TYPE=="REF")%>%group_by(Gene)%>%
#  summarise(freq=mean(as.double(gene_coverage)))%>%filter(freq<0.75)%>%select(Gene)%>%as.list
#true_start<-genes_met%>%group_by(Gene)%>%
#  summarise(cons_start=max(as.double(Start)))
gene_freq<-genes_met%>%group_by(Gene)%>%
  summarise(freq=mean(as.double(gene_coverage)))
ancient_genes<-intersect(non_ref_genes$Gene,in_ancient$Gene)
#remove core genes

genes_met_fxd<-left_join(genes_met,gene_freq,by="Gene")%>%
  mutate(anc_cov=ifelse(TYPE=="MAG"&gene_coverage==1,3,gene_coverage))%>%
  mutate(anc_gencov=ifelse(TYPE=="MAG"&gene_coverage==0,2,anc_cov))
ancient_met<-filter(genes_met_fxd,Gene%in%ancient_genes)
vir_met<-filter(genes_met_fxd,Gene%in%genes)
#genelabel<-ancient_met%>%select("cons_start","Gene")%>%arrange(cons_start)%>%unique
ancient_met<-filter(genes_met_fxd,Gene%in%jackiedat$Gene)




### plot with phylo tree with mag placements
library(aplot)
tree <- read.tree(paste("~/PhD/datafiles/placement_5/placementmag/","merge.newick",sep=""))
tree<-root(tree,"GCF_002355215N1",resolve.root = TRUE)
#ultra<-force.ultrametric(tree)
#is.ultrametric(ultra)
#dend_mut <- as.dendrogram(ultra)
g2<-ggtree(tree,branch.length = "none")+geom_rootedge(rootedge=0.0005)+
    ggtitle("Placement merged\nPhylogenetic Tree")
p9<-p9+theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=6), #change legend title font size
        legend.text = element_text(size=6)) #change legend text font size
p12 <- ggplot(ancient_met, aes(y=NAME, x=reorder(Gene,freq),fill=as.character(anc_gencov))) + ggtitle("UA159 absent Bacteriocins")+
  geom_tile()  + theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                       legend.key.size = unit(0.2, 'cm'), #change legend key size
                       legend.key.height = unit(0.2, 'cm'), #change legend key height
                       legend.key.width = unit(0.2, 'cm'), #change legend key width
                       legend.title = element_text(size=10), #change legend title font size
                       legend.text = element_text(size=10)) +#change legend text font size                                                       
  scale_fill_manual(name="Gene coverage",
  labels=c("No", "Yes","No_ancient","Yes_ancient"),values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  xlab(NULL) + ylab(NULL)  
p13 <- ggplot(vir_met, aes(y=NAME, x=reorder(Gene,freq),fill=as.character(anc_gencov))) + ggtitle("selected virulence factors\n(BAKTA)")+
  geom_tile() + 
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1),
                       legend.key.size = unit(0.3, 'cm'), #change legend key size
                       legend.key.height = unit(0.3, 'cm'), #change legend key height
                       legend.key.width = unit(0.3, 'cm'), #change legend key width
                       legend.title = element_text(size=8), #change legend title font size
                       legend.text = element_text(size=8)) +#change legend text font size                                                       
  scale_fill_manual(name="Gene presence/absence",
                    labels=c("No", "Yes","No_ancient","Yes_ancient"),values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  #scale_x_discrete(label=genenames$gene[3:8])+
  xlab(NULL) + ylab(NULL)  
p13
#pdf("~/PhD/datafiles/raregenes.pdf")#, width=8, height=25)
pdf("~/PhD/datafiles/nonrefbacteriocin.pdf",width=10)
p12%>%insert_left(p9,width=0.04)%>%insert_left(g2,width=0.4)
dev.off()
#p1%>%insert_left(p3,width=0.1)







p15 <- ggplot(vir_met, aes(y=NAME, x=factor(Gene,level=genenames$gene[3:8]),fill=as.character(anc_gencov))) + ggtitle("selected virulence\nfactors (BAKTA)")+
  geom_tile() + 
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)) +#change legend text font size                                                       
  scale_fill_manual(name="Gene presence/absence",
                    labels=c("No", "Yes","No_ancient","Yes_ancient"),values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  #scale_x_discrete(label=genenames$gene[3:8])+
  xlab(NULL) + ylab(NULL)  
p15
p14 <- ggplot(genemat, aes(y=sampleId, x=id)) + ggtitle("selected virulence\nfactors (Reference mapping)") + 
  #geom_tile(aes(fill=as.character(gene_coverage)))+
  geom_tile(aes(fill=normavgcov))+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 7, angle = 55, hjust = 1))+
  xlab(NULL) + ylab(NULL)+
  scale_fill_viridis(name="Gene coverage",trans = 'reverse',option = "magma")+
  #scale_fill_manual(name="Gene coverage")+
  #labels=c("No", "Yes","No_ancient","Yes_ancient"),
  #values=c("3"="black","2"="lightgrey","1"="blue","0"="white"))+
  scale_x_discrete(label=genenames$gene) 
p14<-p14+theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
               legend.key.height = unit(0.3, 'cm'), #change legend key height
               legend.key.width = unit(0.3, 'cm'), #change legend key width
               legend.title = element_text(size=8), #change legend title font size
               legend.text = element_text(size=8)) #change legend text font size
p14
pdf("~/PhD/datafiles/bakta_mapping.pdf")
p15%>%insert_left(p9,width=0.1)%>%insert_left(g2,width=0.9)%>%insert_right(p14)
dev.off()

#create table for jackie
virtabledat<-vir_met%>%select(c("Gene","Product"))%>%subset(!is.na(Product))%>%unique
count_anc<-ancient_met%>%filter(anc_gencov==3)%>%group_by(Gene)%>%
  summarise(Ancient_count=sum(gene_coverage))
count_mod<-ancient_met%>%filter(anc_gencov!=3)%>%group_by(Gene)%>%
  summarise(Modern_freq=mean(as.double(gene_coverage)))
crispr<-unique(ancient_met$Product[grep("CRISPR",ancient_met$Product,ignore.case=TRUE)])
transpose<-unique(ancient_met$Product[grep("transposase",ancient_met$Product,ignore.case=TRUE)])
genes<-c("mutA","mutA'","mutB","mutC","mutD","esaA","essA","essB","essC","bacA","macB")
COG_categories <- read.delim("~/PhD/datafiles/COG_categories", header=FALSE,col.names = c("catref","COG_Category"))
pertabledat<-ancient_met%>%subset(!is.na(Product))%>%
  filter(!Product%in%transpose)%>%filter(!Product%in%crispr)%>%
  mutate(COGcheck=ifelse(grepl("COG",DbXrefs),DbXrefs,NA))%>%
  separate(COGcheck,c("COG1","COG2"),sep=",",remove = FALSE,extra = "drop")%>%
  mutate(COG_Function=gsub("COG:","",ifelse(grepl("COG:COG",COG1),COG1,COG2)))%>%
  mutate(catref=gsub(" ","",gsub("COG:","",ifelse(grepl("COG:COG",COG1),COG2,COG1))))%>%
  left_join(COG_categories,by="catref")%>%
  mutate(Abranches_Category=ifelse(Gene%in%genes,"Bacteriocin",NA))%>%
  mutate(Category=ifelse(is.na(COG_Category),Abranches_Category,COG_Category))%>%
  select(c("Gene","Product","Category"))%>%
  left_join(count_anc,by="Gene")%>%left_join(count_mod,by="Gene")%>%unique%>%group_by(Gene)



#after getting results back from jackie, i check whether bacteriocins are enriched in my old samples
jackiedat <- read_excel("Downloads/nonref_gene_tableJA.xlsx")%>%filter(Category=="Bacteriocin/antimicrobial")%>%
  select("Gene","Ancient_count","Modern_freq")%>%unique
coregenes<-c("gbpA","atpD","ftsN","dexA") #these genes are actually core genes which we only find a few times here
#gene number cwlC, ftsN and rodA, and  are both cell envelope homeostasis genes
#ftsN should not have a product annotation
#gbpA and other "should be more common genes" for which i can find homologous functions in the reference genome need to be added the "-like" suffix
# anything related to recombinatio, nuclease and repair need to be put together in the same category "Nucleic acid replcation, recombination and repair
#Mg and Co transporter also falls within the transporter/exporter category
#remove any question marks! We know
#Listeria sRNA rli28 is just called smallRNA
#Purine metabolsin goes to central metabolism (purQ)
  
library(xlsx)
read


write.xlsx(as.data.frame(virtabledat),"~/PhD/datafiles/virulence_gene_table.xlsx",  
           col.names = TRUE, row.names = FALSE, append = FALSE)
write.xlsx(as.data.frame(pertabledat),"~/PhD/datafiles/nonref_gene_table.xlsx",  
           col.names = TRUE, row.names = FALSE, append = FALSE)

grepl("transposase",ancient_met$Product,ignore.case=TRUE)



#get genes not in ancient, which cannot be found through reference mapping
genes_met<-left_join(comb_genes,combined_final,by="NAME")
in_ancient<-genes_met%>%filter(TYPE=="MAG")%>%group_by(Gene)%>%
  summarise(freq=mean(as.double(gene_coverage)))%>%filter(freq<0.01)%>%select(Gene)%>%as.list
non_core_genes<-genes_met%>%filter(TYPE=="REF")%>%group_by(Gene)%>%
  summarise(freq=mean(as.double(gene_coverage)))%>%filter(freq>0.2)%>%select(Gene)%>%as.list
gene_freq<-genes_met%>%group_by(Gene)%>%
  summarise(freq=mean(as.double(gene_coverage)))
ancient_genes<-intersect(non_core_genes$Gene,in_ancient$Gene)
#remove core genes

genes_met_fxd<-left_join(genes_met,gene_freq,by="Gene")%>%
  mutate(anc_cov=ifelse(TYPE=="MAG"&gene_coverage==1,3,gene_coverage))%>%
  mutate(anc_gencov=ifelse(TYPE=="MAG"&gene_coverage==0,2,anc_cov))
mubgenes<-unique(genes_met_fxd$Gene[grep("mub",genes_met_fxd$Gene,ignore.case=TRUE)])
casgenes<-unique(genes_met_fxd$Gene[grep("cas",genes_met_fxd$Gene,ignore.case=TRUE)])
trna<-unique(genes_met_fxd$Gene[grep("trn",genes_met_fxd$Gene,ignore.case=TRUE)])
#crispr<-unique(genes_met_fxd$Gene[grep("trn",genes_met_fxd$Gene,ignore.case=TRUE)])
ribosomal<-unique(genes_met_fxd$Gene[grep("ribosomal RNA",genes_met_fxd$Product,ignore.case=TRUE)])
ancient_met<-filter(genes_met_fxd,Gene%in%ancient_genes)%>%filter(!Gene%in%c(mubgenes,casgenes,crispr,ribosomal))



####Ref mapping CGG101233 against all other ancient strep
scov<-read.delim("PhD/datafiles/CGG101233.concat_gene_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("contig","start","end","id","gene","avgcov","sampleId"))
scds<-read.delim("PhD/datafiles/CGG101233.concat_cds_cov.tsv", header=FALSE, #i will look for the genes in this gene coverage analysis file
                 col.names=c("contig","start","end","id","desc","avgcov","sampleId"))
combined_cov<-full_join(scov,scds,by=c("contig","start","end","id","avgcov","sampleId"))
scov<-combined_cov
scov$sampleId<-gsub("\\.","N",scov$sampleId)
#complete_sample<-left_join(sample_info,sample_meta,by="sampleId") #get the metadata on to my sample classification table
genes<-c("dltC","ffh","cnm","cnaB","cbpA","cbm")  #define the genes i want to look for
#,"gtfA","gtfB","gtfC"
#"spaP","spaA","spaB","spaC"
genes_of_interest<-scov[which(scov$gene%in%genes),]
combined_final$sampleId<-combined_final$NAME
pathogen7<-left_join(scov,combined_final,by="sampleId")%>%left_join(strep_table,by="sampleId")
casgenes<-unique(pathogen7$gene[grep("cas",pathogen7$gene,ignore.case=TRUE)])
pathogen5<-pathogen7%>%mutate(normavgcov=ifelse(is.na(coveragePExp),avgcov,ifelse(avgcov/coveragePExp>1,1,avgcov/coveragePExp)))%>%
  mutate(genelength=end-start)%>%mutate(bincov=ifelse(normavgcov>0.4,1,0))%>%
  mutate(anc_cov=ifelse(TYPE!="REF"&bincov==1,3,bincov))%>%
  mutate(gene_coverage=ifelse(TYPE!="REF"&bincov==0,2,anc_cov))%>%filter(!gene%in%c(mubgenes,casgenes,crispr,ribosomal))
ref_genes<-genes_met%>%filter(NAME=="GCF_009738105N1")%>%filter(gene_coverage!=0)%>%select(Gene)
pathogen4<-pathogen5[!pathogen5$gene%in%ref_genes$Gene,]%>%filter(gene!="protein_coding")


pathogen4$ageAverage2<-as.double(pathogen4$ageAverage)
ggplot(pathogen4,aes(x=reorder(sampleId,ageAverage2),y=gene))+ #mathematical racism
  geom_tile(color = "white",aes(fill=normavgcov))+ #define heatmap
  #facet_grid(~age_period,scale='free',space='free',switch='y')+# #discriminate people by birthplace
  #scale_fill_manual()+
 # plot_style + #do it with style!
  theme(axis.text.x = element_text(size = 7, angle = 55, hjust = 1))
  #geom_point(data = strep_table,  #add the age of the people
  #           aes(color = ageAverage), size = 8)#+
  #scale_fill_viridis(option="magma")+
  #scale_fill_discrete(name="Gene coverage",
  #                    labels=c("No", "Yes"))





####matrix method
library(Matrix)

virulencetestmat<-comb_genes%>%select("NAME","Gene","gene_coverage")%>%unique%>%subset(Gene != "")
non_core_genes<-virulencetestmat%>%group_by(Gene)%>%
  summarise(freq=mean(as.double(gene_coverage)))%>%filter(freq<0.99)%>%select(Gene)%>%as.list

#remove core genes

virulencetestmat<-virulencetestmat[virulencetestmat$Gene%in%non_core_genes$Gene,]
genemat<-spread(virulencetestmat,key=Gene,value=gene_coverage)
rownames<-genemat$NAME
genemat<-select(genemat,-NAME)%>%as.matrix
row.names(genemat)<-rownames

heatmap3(genemat)


heatmap3(genemat,distfun=dist)





