library(tidyverse)
library(ggforce)
library(patchwork)
library(ggrepel)
library(viridis)

sampleInfo <- read.delim("~/PhD/datafiles/20220603.diseases.sampleInfo.tsv", header=TRUE)
samples10e5kmers <- read.table("~/PhD/Strep_mapping/data/sample_info/samples10e5kmers.txt", quote="\"", comment.char="",col.names = c("sampleId"))
anvio_candidates <- read.delim("~/PhD/Strep_mapping/data/sumdat/anvio_candidates.tsv", col.names=c("contig","coverage","ani","sampleId"))
connected<-intersect(anvio_candidates$sampleId,neo.impute.1000g.sampleInfo$sampleId)
unconnected<-setdiff(anvio_candidates$sampleId,neo.impute.1000g.sampleInfo$sampleId)
tax_info<-neo.impute.1000g.sampleInfo[which(neo.impute.1000g.sampleInfo$sampleId %in% connected),]

sampleInfo <-read.delim("PhD/community_analysis_09-22/tax_pipeline/20230101.diseases.sampleInfo.tsv") %>% 
  separate(groupLabel,c("Country","age_period","misc"),sep="_",remove = FALSE,extra = "merge")%>%
  select("sampleId","clusterLabel","country","latitude","longitude","ageAverage","region","age_period")

sampleInfo<-rbind(sampleInfo,c("CGG100534",NA,"Denmark",57.047,9.98,as.double(279),"NorthernEurope","Historical"))
sampleInfo<-rbind(sampleInfo,c("CGG100272",NA,"Iran",32.000,53.000,as.double(2480.0),"CentralAsia","IronAge"))
sampleInfo<-rbind(sampleInfo,c("CGG101233",NA,"Afghanistan",34.818813, 67.838997,as.double(2675),"CentralAsia","IronAge"))
sampleInfo$ageAverage<-as.double(sampleInfo$ageAverage)
sampleInfo$latitude<-as.double(sampleInfo$latitude)
sampleInfo$longitude<-as.double(sampleInfo$longitude)
classification <- read.delim("~/PhD/PhylTreeScripts/data/accession_tree2")%>%
  filter(TYPE!="REF")%>%select("NAME","TYPE") #get names of ancient samples
colnames(classification)<-c("sampleId","quality")
tax_info2<-left_join(classification,sampleInfo,by="sampleId")#%>%mutate(COUNTRY = coalesce(COUNTRY.x, COUNTRY.y)) %>% 
  #mutate(DATE = coalesce(DATE.x, DATE.y)) %>%mutate(SOURCE = coalesce(SOURCE.x, SOURCE.y))%>%
  #select(NAME, COUNTRY, DATE,TYPE,COV,SEROTYPE,SOURCE) %>% mutate_if(is.character, ~na_if(., '')) #combine all the sample metadata into one file

tax_info2$age_period<-gsub("VikingAge","Historical",tax_info2$age_period)
strep_table$age_period<-gsub("VikingAge","Historical",strep_table$age_period)
strep_table$ageAverage<-as.double(strep_table$ageAverage)



world <- map_data("world")
#pdf("/home/vincentven/PhD/anvio_sample_origin", width=20, height=7)
#this plot needs to be plotted with a 46 to 125 height to width ratio = 2.717391
#geom_jitter(position = position_jitter(seed = 1)) +
#  geom_text_repel(aes(label = ifelse(wt > 3, wt, "")), position = position_jitter(seed = 1))
plot_world<-ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    colour = "black", fill = "white", size = 0.1,alpha=0.99,#width=3,height=1
  ) +  theme_void()+geom_point(
    data = tax_info2,#position=position_jitter(width=1,height=1,seed=1),
    aes(longitude, latitude,color=ageAverage,label=sampleId,
        shape=factor(age_period, 
        levels = c("Historical","IronAge", "BronzeAge","Neolithic","Mesolithic"))),
        alpha = 0.99,
  ) +ylim(15,61) + xlim(-15,110)+#xlim(-50,165) #ylim(41, 73)  #+ ylim(37, 73) + xlim(-50,165) 
  scale_shape_manual(values = c(15,16,17,18,19))+
  scale_color_viridis(option = "viridis")+
  labs(title="Mutans sample overview",y="latitude",x="longitude")+
  geom_text_repel(data = tax_info2,aes(longitude, latitude,label=sampleId),size=3,
                  min.segment.length = unit(0.2, 'lines'))+#,position=position_jitter(width=1,height=1,seed=1))+
  guides(col=guide_colourbar("Sample age"),shape=guide_legend("Archaeological\ncontext"))
#dev.off()
plot_world
plot_mutans <-ggplot(strep_table, aes(y=coveragePExp, x=coverageAvg, 
                                      shape=factor(age_period, 
                                      levels = c("Historical","IronAge", "BronzeAge","Neolithic","Mesolithic")))) +
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
  #
plot_mutans






NEOproject<- tax_info[grep("NEO",tax_info$sampleId),]
RISEproject<- tax_info[as.double(grep("RISE",tax_info$sampleId)),]
DAproject<- tax_info[as.double(grep("DA",tax_info$sampleId)),]



sample_table<-sample_metadata%>%filter(age_period!="Modern"&material!="dental calculus")%>%filter(label%in%row.names(amawdmg@sam_data))
sample_table$latitude<-as.double(sample_table$latitude)
sample_table$longitude<-as.double(sample_table$longitude)
sampleInfo<-sample_table%>%
  mutate(country=replace(country,country=="Faroes","Faroe Islands"))%>%
  mutate(country=replace(country,country=="Great Britain","UK"))%>%
  mutate(country=replace(country,country=="CzechRepublic","Czech Republic"))
country_df<-sampleInfo%>%count(country)
country_df$Samples<-country_df$n
world_map<-map_data("world")
sampleInfo$Estimated_age<-sampleInfo$ageAverage
valid_countries<-which(world_map$region%in%unique(country_df$country))
new_map<-world_map[valid_countries,]
#pdf("/home/vincentven/PhD/world_sample_origin", width=10, height=4)
p1<-ggplot(country_df) +
  geom_map(aes(map_id = country, fill = Samples,), map = new_map) +#ylim(25,85) + xlim(-73,200)+
  geom_polygon(data = new_map, aes(x = long, y = lat, group = group), colour = 'black', fill = NA) + 
  theme_void() + 
  coord_fixed()+
  ggtitle("Sample origin by country") + 
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + 
  scale_fill_viridis(option="plasma")

p2<-ggplot(sampleInfo, aes(x="",y=ageAverage), height = 50) + 
  scale_y_log10()+
  xlab("Samples")+ylab("Estimated Age")+
  ggtitle("Sample dating")+
  geom_sina()+ coord_flip()

p1 / p2 +plot_layout(heights = c(5, 1))+plot_annotation(title = 'Sample Overview',tag_levels = "A")
#dev.off()


pdf("/home/vincentven/PhD/notes/plots/mutans_sample_overview",width=10, height=6)
plot_world / plot_mutans+ plot_layout(heights = c(4,2),guides = 'collect')+plot_annotation(tag_levels="A")
dev.off()

