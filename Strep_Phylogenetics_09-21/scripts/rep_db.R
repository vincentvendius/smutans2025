library(reshape2)
library(ggplot2)
library(tidyr)
library(ggtree)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("ggtree")



#Read in FASTANI alignment
#fastani <- read.delim("/home/vincentven/PhD/Strep_Phylogenetics_09-21/data/Porphy/fani_GCFtax.out", header=FALSE, col.names = c("taxFile","File","taxAlignment","Alignment","ANI"))
#db_meta <- read.delim("/home/vincentven/PhD/Strep_Phylogenetics_09-21/data/Porphy/kraken_db_metadata",header=FALSE,col.name=c("File","Species"))
#fastani <- read.delim("/home/vincentven/PhD/Strep_Phylogenetics_09-21/data/mutans_refs/fani_GCFtax.out", header=FALSE, col.names = c("taxFile","File","taxAlignment","Alignment","ANI"))
db_meta <- read.delim("/home/vincentven/PhD/Strep_Phylogenetics_09-21/data/kraken_db_metadata",header=FALSE,col.name=c("File","Species"))
fastani <- read.delim("/home/vincentven/PhD/Strep_Phylogenetics_09-21/data/fani_GCFtax.out", header=FALSE, col.names = c("File","Alignment","taxFile","taxAlignment","ANI"))
#fastani <- read.delim(snakemake@input[["fastout"]], header=FALSE, col.names = c("File","Alignment","taxFile","taxAlignment","ANI"))
#db_meta <- read.delim(snakemake@input[["dbmeta"]],header=FALSE,col.name=c("File","Species"))
#File preprocessing of taxonomic info linked to fastani output for species
fastani <- transform(fastani,ANI=1-(ANI/100)) #transform to 1-ANI
Filenames <-fastani[ , which(names(fastani) %in% c("File","taxFile"))] #take out the taxonomy info
colnames(Filenames) <- c("File","Tax")
Alnnames <- fastani[ , which(names(fastani) %in% c("Alignment","taxAlignment"))] #take out the taxonomy info for aligned files
colnames(Alnnames) <- c("File","Tax")
faninames <- rbind(Filenames,Alnnames) #combine them into a dictionary of GCF and tax names
taxinfo <- faninames[!duplicated(faninames$File),] #remove duplicate names
taxinfo <- merge(taxinfo,db_meta) #get species info
row.names(taxinfo)<-taxinfo$File  
taxinfo <- taxinfo[,c("Species"), drop = F]
taxinfo[,c("Species")] = as.character(taxinfo$Species)
taxinfo.freq = as.data.frame(table(taxinfo)) #find the frequency of different tax names
taxinfo.freq[,c("taxinfo")] = as.character(taxinfo.freq$taxinfo)
other = taxinfo.freq$taxinfo[taxinfo.freq$Freq < 2] #define infrequent taxa as other for better overview
taxinfo_all <- taxinfo #keep seperate variable as original file
#define infrequent species as other
taxinfo$Species = ifelse(taxinfo$Species %in% other, "other", taxinfo$Species) 

color_pal<-c("#000000","#69b3a2","#3357FF","#33FF4E","#FF33BB","#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
             "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
             "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
             "#8A7C64", "#599861","#c00000","#AA4371")
color_pal2<-c(rep("#000000",18),"#689030",rep("#000000",14))

#convert fastani output to distance matrix and do hierachical clustering
ani_long <- fastani[,which(names(fastani)%in% c("File","Alignment","ANI"))]
ani_wide <- spread(ani_long,Alignment,ANI,fill=0.25) #fill missing values as having ANI=0.75
row.names(ani_wide)<-ani_wide$File #convert to matrix
ani_wide<-ani_wide[,-1]
ani_wide<-as.matrix(ani_wide)
clustered = hclust(as.dist(ani_wide), method = "average") #compute hierarchical clustering
#clustered <- nj(as.dist(ani_wide))




#get representative sequences
library(TraMineR)
#This package sorts the genomes in the dist matrix by the amount of sequences in their neighbourhood
#https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.385.9760&rep=rep1&type=pdf

# the neighbourhood defined by a radius of 0.01
# Iteratively, a sorted list of samples by neighbour hood density is generated, 
#the top samples is selected as representative, and neighbour density is 
#calculated a new, untill the top sample is a singleton by the radius criterion
repseqs<-dissrep(ani_wide,criterion="density",coverage=1,nrep=115)
summary(repseqs)
repindex<-repseqs[1:length(repseqs)]

#plot the representative subset
#use version of taxinfo with no infrequent taxa defined as other
taxinfo_all.freq = as.data.frame(table(taxinfo_all)) #find the frequency of different tax names
taxinfo_all.freq[,c("taxinfo")] = as.character(taxinfo_all.freq$taxinfo_all)
other = taxinfo_all.freq$Species[taxinfo_all.freq$Freq < 2] #define infrequent taxa as other for better overview
taxinfo_all$Species = ifelse(taxinfo_all$Species %in% other, "other", taxinfo_all$Species)
ani_rep<-ani_wide[repindex,repindex]
clustered = hclust(as.dist(ani_rep), method = "average") #compute hierarchical clustering
#remove this to get plot of representative (ani_rep) instead of everything
#clustered = hclust(as.dist(ani_wide), method = "average") #compute hierarchical clustering

#create phylogenetic trees for the database and its representative subset
#phylo <- ggtree(clustered, layout = "circular") + 
#  ggtitle("Phylogenetic tree of custom kraken database for Streptococcus\nHierarchical clustering")
#phyloring <- gheatmap(phylo,taxinfo, width = 0.06,colnames=F) + 
#  theme(legend.key.size = unit(0.3, 'cm'))+scale_fill_manual(name="Species Name",values=color_pal)
#phyloring
#pdf(snakemake@output[["dbphylo"]])
rpdb <- ggtree(clustered, layout = "circular") + 
  #geom_tiplab()+#aes(subset=TYPE!="REF"),size = 5,offset=.001)+
  #ggtitle("Representative subset of the Streptococcus Genus \nbased on 80 percent AverageNucleotideDistance")
  ggtitle("Phylogenetic tree of representative subset of \ncustom kraken database for Streptococcus",subtitle="S. mutans marked with red diamond")+
  geom_point2(aes(subset=(node==8)), shape=23,size=5, fill='red')
rpdbring <- gheatmap(rpdb,taxinfo_all, width = 0.08,offset=0.01,colnames=F) + 
  theme(legend.key.size = unit(0.3, 'cm'))+scale_fill_manual(name="Species Name",values=color_pal)     
rpdbring
#dev.off()

#change ani_wide to ani_rep for representative subset
#i changed this to wide for the sake of Porphy
ani_dat<-as.data.frame(ani_wide,row.names = row.names(ani_wide))
write(row.names(ani_wide),"rep_db_list")

#this bullshit is because porphy has the wrong samples names from reference mapping
link <- read.table("~/PhD/Strep_Phylogenetics_09-21/data/Porphy/linkingGCFstrain", quote="\"", comment.char="")
colnames(ani_dat) <- link$V2[match(colnames(ani_dat), link$V1)]
row.names(ani_dat) <- link$V2[match(row.names(ani_dat), link$V1)]
write.table(ani_dat,"porphy_ref_rep_ani.matrix",row.names=TRUE)
#write(row.names(ani_rep),snakemake@output[["rpdb"]]) #write the representative samples to a file
#write.table(ani_dat,snakemake@output[["dist_matrix"]],row.names = TRUE)
