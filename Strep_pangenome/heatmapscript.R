modules <- read.delim("~/PhD/datafiles/pangenome/newanalysis/Mutans_enriched_modules.txt")%>%filter(p_MODERN>0.8)%>%filter(p_MAG<0.8)%>%rename(module=accession) #get metabolic modules 
mod_info<-read.delim("~/PhD/datafiles/pangenome/newanalysis/modules_info.txt")%>%filter(module%in%modules$module)
mod_cov <- read.delim("~/PhD/datafiles/pangenome/newanalysis/Mutans_metabolism_mat-module_stepwise_completeness-MATRIX.txt")%>%
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



#time to look at general functional enrichment
func_mod <- read.delim("~/PhD/datafiles/pangenome/newanalysis/Functional_enriched_modules.txt")
mag_func<-func_mod[grep("MAG",func_mod$associated_groups),]%>%filter(p_ANCIENT>0.42)%>%filter(p_MODERN<0.1)

func_mod_count <- read.delim("~/PhD/datafiles/pangenome/newanalysis/Functional_enriched_modules.tsv")%>%filter(X%in%mag_func$COG20_FUNCTION)
data_wide2<- func_mod_count[,-1]
row.names(data_wide2)<-func_mod_count[,1]
data_wide<-as.matrix(data_wide2)
data_wide <- replace(data_wide,is.na(data_wide),0)
heatmap.2(t(data_wide), 
          col = heat.colors(20),    # Choose the color palette (you can change this)
          scale = "column",         # Scale the values by column
          hclustfun = function(x) hclust(x, method="average"),
          Rowv = TRUE, Colv = TRUE,     # Turn off row and column clustering
          margins = c(5, 10),       # Set margins for row and column labels
          main = "Coverage for selected modules", # Set the title of the plot
          xlab = "Metabolic modules", ylab = NA,
          colRow = NA, colCol = "black", # Turn off color for row dendrogram and use black for column dendrogram
          trace = "none",           # Turn off trace
          #key = TRUE, keysize = 1.0, cexRow = 0.8, cexCol = 0.8, # Adjust key parameters
          density.info = "none",    # Turn off density plot
          dendrogram = "both")    # Show only column dendrogram
#ColSideColors = column_metadata)  # Add metadata sidebar for columns


heatmap.2(t(data_wide), 
          #col = colorRampPalette(brewer.pal(9, "Blues"))(20),  # Change the color palette (you can choose from other palettes)
          scale = "column",         # Scale the values by column
          Rowv = NA, Colv = NA,     # Turn off row and column clustering
          margins = c(5, 10),       # Set margins for row and column labels
          main = "Annotated Heatmap", # Set the title of the plot
          xlab = "Columns", ylab = "Rows",
          colRow = NA, colCol = "black", # Turn off color for row dendrogram and use black for column dendrogram
          trace = "none",           # Turn off trace
          #key = TRUE, keysize = 1.0, cexRow = 0.8, cexCol = 0.8, # Adjust key parameters
          density.info = "none",    # Turn off density plot
          dendrogram = "column",    # Show only column dendrogram
          #ColSideColors = column_metadata,  # Add metadata sidebar for columns
          #lmat = cbind(c(0, 3), c(2, 1)),  # Add space for metadata legend
          #lhei = c(1, 4),           # Adjust the height of the legend
          #cex.axis = 0.8,           # Adjust axis label size
          #cex.lab = 0.8,            # Adjust axis title size
          #lty = c(1, 1),            # Set legend line types
          #lwid = c(2, 2),           # Set legend line widths
          #key.title = "Legend",     # Set legend title
          #key.xlab = "Values",      # Set legend x-axis label
          key.ylab = NULL)          # Remove legend y-axis label



ggplot(long_count, aes(y=Sample, x=X,fill=Count))+geom_tile()
long_count<-func_mod_count%>%gather("Sample","Count",2:33)

#Only on HQ mags
func_mod <- read.delim("~/PhD/datafiles/pangenome/newanalysis/HQ/Functional_enriched_modules.txt")
mag_func<-func_mod[grep("MAG",func_mod$associated_groups),]#%>%filter(p_MODERN<0.5)%>%filter(p_ANCIENT<0.8)

func_mod_count <- read.delim("~/PhD/datafiles/pangenome/newanalysis/Functional_enriched_modules.tsv")#%>%filter(X%in%mag_func$COG20_FUNCTION)

#HQ mags + all modern mutans from phylogenetic analysis
func_mod <- read.delim("~/PhD/datafiles/pangenome/hq_analysis/testFunctional_enriched_modules.txt")
mag_func<-func_mod%>%filter(p_MODERN>0.6)#%>%filter(p_ANCIENT<0.5)


#metabolic analysis
modules <- read.delim("~/PhD/datafiles/pangenome/hq_analysis/Mutans_enriched_modules.txt")

mod_info<-read.delim("~/PhD/datafiles/pangenome/hq_analysis/modules_info.txt")
mod_tab<-modules%>%filter(p_MODERN>0.8)%>%filter(p_ANCIENT<0.8)%>%rename(module=accession)%>%
  left_join(modern_enriched)%>%select("module","name","category","subcategory","sample_ids","p_GENUS","p_MODERN","p_ANCIENT")




#With all phylo samples, mags and 20 other species
#Mutans Metabolic modules
modules_all <- read.delim("~/PhD/datafiles/pangenome/magphylogenus/Mutans_enriched_modules.txt")%>%filter(p_MODERN>0.9)%>%filter(p_GENUS<0.5)

func_mod <- read.delim("~/PhD/datafiles/pangenome/magphylogenus/Functional_enriched_modules_FUNCTION.txt")

#Mutans Core functions
mag_func<-filter(func_mod,p_GENUS<0.1)%>%filter(p_MODERN>0.5)
#Modern unique
mag_func1<-filter(func_mod,p_MODERN>0.5)%>%filter(p_ANCIENT<0.3)%>%
  select(COG20_FUNCTION,gene_clusters_ids,p_GENUS,p_MODERN,p_ANCIENT)
write.xlsx(as.data.frame(mag_func),"~/PhD/datafiles/Modern_enriched_functions.xlsx",  
           col.names = TRUE, row.names = TRUE, append = FALSE)
#Ancient unique
mag_func_ancient<-func_mod[grep("ANCIENT",func_mod$associated_groups),]%>%filter(p_MODERN<0.13)%>%filter(p_ANCIENT>0.1)%>%
  select(COG20_FUNCTION,gene_clusters_ids,p_GENUS,p_MODERN,p_ANCIENT)
write.xlsx(as.data.frame(mag_func),"~/PhD/datafiles/Ancient_enriched_functions.xlsx",  
           col.names = TRUE, row.names = TRUE, append = FALSE)

func_mod_all <- read.delim("~/PhD/datafiles/pangenome/magphylogenus3/Functional_enriched_modules_id.txt")
#Mutans Core functions
mag_func_all<-filter(func_mod_all,p_GENUS<0.01)%>%filter(p_MODERN>0.999)
#Ancient unique
mag_func_all<-func_mod_all[grep("ANCIENT",func_mod_all$associated_groups),]%>%filter(p_MODERN<0.3)%>%filter(p_ANCIENT>0.2)



#Get gene cluster names 
genenames<-filter(func_mod_all,p_MODERN>0.95)%>%filter(p_GENUS<0.05)%>%pull(IDENTITY)%>%as.data.frame
#which of these core mutans gene clusters are uncommon in the ancient samples
filter(func_mod_all,p_MODERN>0.99999)%>%filter(p_GENUS<0.000001)%>%filter(p_ANCIENT<0.8)
#not many
#how many of these gene clusters are not only mutans but also modern mutans unique
modmutgene<-strsplit(paste0(mag_func$gene_clusters_ids,collapse=", "),", ")[[1]]
intersect(modmutgene,genenames)
#none

#look at the gene clusters, and what they represent in the reference mutans genome UA159
core_genes <- read.delim("~/PhD/datafiles/pangenome/magphylogenus3/core_genes.tsv", col.names=c("name","gene","str","end","cov"))%>%
  filter(cov==1)

core_cds <- read.delim("~/PhD/datafiles/pangenome/magphylogenus3/core_cds.tsv", col.names =c("id","desc","str","end","cov"))%>%
  filter(cov==1)
combined_core<-full_join(core_cds,core_genes,by=c("str","end","cov"))
library(xlsx)
write.xlsx(as.data.frame(combined_core),"~/PhD/datafiles/Mutans_unique_genes_v2.xlsx",  
           col.names = TRUE, row.names = TRUE, append = FALSE)



#----------------New iteration---------#
#Sending updated files raw to collaborators
func_cat <- read.delim("~/PhD/datafiles/pangenome/magphylogenus4/Functional_enriched_function.txt")

write.xlsx(as.data.frame(func_cat),"~/PhD/datafiles/pangenome/magphylogenus4/pangenome_kegg.xlsx",  
           col.names = TRUE, row.names = TRUE, append = FALSE)

enriched<-filter(func_cat,p_MODERN<0.1)%>%filter(p_ANCIENT>0.2)
#save this one

#Modern unique
mag_func1<-filter(func_cat,p_MODERN>0.8)%>%filter(p_GENUS<0.1)%>%
  select(COG20_FUNCTION,gene_clusters_ids,p_GENUS,p_MODERN,p_ANCIENT)

mag_func_ancient<-filter(func_cat,p_MODERN<0.01)%>%filter(p_ANCIENT>0.1)%>%
  select(COG20_FUNCTION,gene_clusters_ids,p_GENUS,p_MODERN,p_ANCIENT)


write.xlsx(as.data.frame(mag_func1),"~/PhD/datafiles/pangenome/magphylogenus3/pangenome_mutansv2.xlsx",  
           col.names = TRUE, row.names = TRUE, append = FALSE)
#Ancient unique
mag_func_ancient<-filter(func_cat,p_MODERN<0.01)%>%filter(p_ANCIENT>0.1)%>%
  select(COG20_FUNCTION,gene_clusters_ids,p_GENUS,p_MODERN,p_ANCIENT)
write.xlsx(as.data.frame(mag_func_ancient),"~/PhD/datafiles/pangenome/magphylogenus3/pangenome_ancient.xlsx",  
           col.names = TRUE, row.names = TRUE, append = FALSE)


#---------Comparing different mcl parameters
mcl<-"10"
ancients<-c("CGG100272","NEO137","NEO105","NEO938","VK63","CGG101233","CGG100534")
pansum <- read.delim(paste0("~/Documents/pangenome summaries/mcl",mcl,"_gene_clusters_summary.txt"))

length(unique(pansum$gene_cluster_id))
length(unique(pansum$COG20_FUNCTION))


newhomogen<-pansum%>%group_by(gene_cluster_id)%>%summarize(het=mean(combined_homogeneity_index))
newhomogen$mcl<-mcl
homogen<-rbind(homogen,newhomogen)


p<-ggplot(homogen, aes(x=het,color=mcl)) +
  geom_density()+
  ggtitle("Gene cluster combined homogenity density")+xlab("combined homogenity score")
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
  #           linetype="dashed")
p

ancientpan<-filter(pansum,genome_name%in%ancients)
length(unique(ancientpan$gene_cluster_id))
length(unique(ancientpan$COG20_FUNCTION))

anewhomogen<-ancientpan%>%group_by(gene_cluster_id)%>%summarize(het=mean(combined_homogeneity_index))
anewhomogen$mcl<-mcl
ahomogen<-rbind(ahomogen,anewhomogen)
p<-ggplot(ahomogen, aes(x=het,color=mcl)) +
  geom_density()+
  ggtitle("Gene cluster combined homogenity density for ancient samples")+xlab("combined homogenity score")
#geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
#           linetype="dashed")
p
