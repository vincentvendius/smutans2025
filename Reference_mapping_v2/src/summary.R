
################################
## summary report for species ##
################################

## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(viridis))
suppressMessages(library(purrr))
suppressMessages(library(Rsamtools))
suppressMessages(library(doParallel))
suppressMessages(library(tidyr))
suppressMessages(library(dendextend))
suppressMessages(library(cowplot))
suppressMessages(library(ggdendro))
suppressMessages(library(gridExtra))
suppressMessages(library(reshape2))

## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
dbPath <- args[2]
sampleId <- args[3]
threads <- as.integer(args[4])
bamPath <- args[5]
outPlot <- args[6]
ani_rep <- read.table(args[7],quote="\"", comment.char="",check.names = FALSE)
outBestmap <- args[8]
out_tsv <- args[9]
## --------------------------------------------------------------------------------
## helpers

registerDoParallel(threads)

## helper, extract first field of prefix for matching in file name assemblyId parsing
prefix1 <- strsplit(prefix, "\\.") %>%
    map_chr(1)
seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "cccccccc")


## --------------------------------------------------------------------------------
## genome coverage summary
#list.files(path=paste(bampath,"/bam",sampleId,sep=""),pattern="genomecov",full.names=TRUE)
f1 <- list.files(path = paste(bamPath,"bam", prefix,  sampleId, sep = "/"),
                 pattern = "genomecov",
                 full.names = TRUE)
dGc <- foreach(ff = f1) %dopar% {
    r1 <- read_tsv(ff,
                   col_type = "ciddd",
                   col_names = c("contigId", "dp", "count", "l", "p"))
    if(nrow(r1) == 0){
        NULL
    } else {

        ## find contigs with coverage
        r11 <- r1 %>%
            filter(dp > 0,
                   contigId != "genome") %>%
            distinct(contigId)

        nContigs <- r1 %>%
            pull(contigId) %>%
            unique() %>%
            length() -1

        ## coverage summary stats
        r2 <- r1 %>%
            filter(contigId %in% r11$contigId) %>%
            group_by(contigId) %>%
            summarise(contigL = l[1],
                      coverageAvg = sum(dp * count) / contigL,
                      coverageSd = sqrt(sum((dp - coverageAvg)^2*count)/sum(count)),
                      coverageBp = sum(count[dp>0]), coverageP = 1 - p[1],
                      coveragePExp = 1 - exp(-coverageAvg),
                      coveragePRatio = coverageP / coveragePExp,
                      coverageCv = coverageSd / coverageAvg,
                      coverageEvennessScore = 1 - sum((ceiling(coverageAvg) - dp[dp <= ceiling(coverageAvg)])*count[dp <= ceiling(coverageAvg)] / (ceiling(coverageAvg) * contigL)))

        r21 <- r1 %>%
            filter(contigId == "genome") %>%
            summarise(contigId = "genome",
                      contigL = l[1],
                      coverageAvg = sum(dp * count) / contigL,
                      coverageSd = sqrt(sum((dp - coverageAvg)^2*count)/sum(count)),
                      coverageBp = sum(count[dp>0]), coverageP = 1 - p[1],
                      coveragePExp = 1 - exp(-coverageAvg),
                      coveragePRatio = coverageP / coveragePExp,
                      coverageCv = coverageSd / coverageAvg,
                      coverageEvennessScore = 1 - sum((ceiling(coverageAvg) - dp[dp <= ceiling(coverageAvg)])*count[dp <= ceiling(coverageAvg)] / (ceiling(coverageAvg) * contigL))) %>%
            mutate(nContigs = nContigs,
                   pContigsCovered = nrow(r11) / nContigs)

        r2 <- bind_rows(r2, r21)## split fields separated by '.' in input file names
        r3 <- ff %>%
            gsub(".*\\/", "", .) %>%
            strsplit("\\.") %>%
            unlist()
        #idx <- match(prefix1, r3) ## index of field where prefix starts
        idx <- match("srt",r3)
        r2 <- r2 %>%
            mutate(genusId = prefix1, #r3[1],
                   assemblyId = paste(r3[1:(idx - 1)], collapse = "."))

        select(r2, genusId, assemblyId, contigId, contigL, coverageAvg:pContigsCovered)
    }
}
dGc <- bind_rows(dGc)


## --------------------------------------------------------------------------------
## bam stats

f1 <- list.files(path = paste(bamPath,"bam",prefix, sampleId, sep = "/"), pattern = "filtered.bam$", full.names = TRUE)
dBam <- foreach(ff = f1) %dopar% {

    p1 <- ScanBamParam(what = c("rname", "mapq", "qwidth", "cigar"),
                       tag=c("NM"),
                       flag = scanBamFlag(isDuplicate = FALSE))
    r1 <- scanBam(ff, param=p1)

    if(length(r1[[1]]$rname)== 0){
        NULL
    } else {
        r2 <- as_tibble(r1[[1]][c("rname", "mapq", "qwidth", "cigar")]) %>%
            mutate(nm = r1[[1]]$tag$NM,
                   rname = as.character(rname))

        ## parse cigar string to count number of soft-clipped bases; map function converts all non-S containing cigar fields into NA
        r21 <- r2 %>%
            pull(cigar) %>%
            gsub("([A-Z])", "\\1.", .) %>%
            strsplit("\\.") %>%
            map(function(x) as.integer(gsub("S", "", x)) %>%
                            sum(na.rm = TRUE)) %>%
            as.integer()
        r2 <- r2 %>%
            mutate(nSoftClip = r21)

        ## generate final summary tibble
        r3 <- gsub(".*\\/", "", f1)

        r41 <- r2 %>%
            group_by(rname) %>%
            count(nm) %>%
            mutate(p = n / sum(n)) %>%
            select(rname, nm, n, p) %>%
            summarise(nReads = sum(n),
                      editDistMode = nm[which.max(n)],
                      editDistAvg = sum(n * nm) / sum(n),
                      editDistAvgDecay = mean(diff(n)) / sum(n),
                      editDistDecayEnd = ifelse(is.na(which(diff(n) > 0)[1]), max(nm), which(diff(n) > 0)[1])) %>%
            ungroup()
        colnames(r41)[1] <- "contigId"

        r411 <- r2 %>%
            count(nm) %>%
            mutate(p = n / sum(n)) %>%
            select(nm, n, p) %>%
            summarise(nReads = sum(n),
                      editDistMode = nm[which.max(n)],
                      editDistAvg = sum(n * nm) / sum(n),
                      editDistAvgDecay = mean(diff(n)) / sum(n),
                      editDistDecayEnd = ifelse(is.na(which(diff(n) > 0)[1]), max(nm), which(diff(n) > 0)[1])) %>%
            mutate(contigId = "genome")
        r41 <- bind_rows(r41, r411)

        r42 <- r2 %>%
            group_by(rname) %>%
            summarise(readLAvg = mean(qwidth),
                      mqAvg = mean(mapq),
                      nSoftClipAvg = mean(nSoftClip),
                      ani = 1 - sum(nm)/sum(qwidth)) %>%
            ungroup() %>%
            select(-rname)
        r421 <- r2 %>%
            summarise(readLAvg = mean(qwidth),
                      mqAvg = mean(mapq),
                      nSoftClipAvg = mean(nSoftClip),
                      ani = 1 - sum(nm)/sum(qwidth))
        r42 <- bind_rows(r42, r421)

        ## split fields separated by '.' in input file names
        r5 <- ff %>%
            gsub(".*\\/", "", .) %>%
            strsplit("\\.") %>%
            unlist()

        #idx <- match(prefix1, r5) ## index of field where prefix starts
        idx <- match("dedup",r5)
        r6 <- bind_cols(r41, r42) %>%
            mutate(genusId = prefix1, #r5[1],
                   assemblyId = paste(r5[1:(idx - 1)], collapse = ".")) %>%
            select(genusId, assemblyId, contigId, nReads, editDistMode:ani)
        r6
    }
}
dBam <- bind_rows(dBam)

dB1 <- dBam %>%
    filter(contigId == "genome") %>%
    group_by(genusId) %>%
    mutate(aniRank = rank(dplyr::desc(ani))) %>%
    ungroup() %>%
    select(assemblyId, contigId, aniRank)

dB2 <- dBam %>%
    filter(contigId == "genome",
           nReads >= 100) %>%
    group_by(genusId) %>%
    mutate(aniRank100 = rank(dplyr::desc(ani))) %>%
    ungroup() %>%
    select(assemblyId, contigId, aniRank100)

dBam <- dBam %>%
    left_join(dB1) %>%
    left_join(dB2)




## --------------------------------------------------------------------------------
## final result table

s1 <- filter(seqInfo, assemblyId %in% dGc$assemblyId) %>%
    select(assemblyId, taxId, taxIdSpecies, taxNameSpecies)

dGc$assemblyId<-gsub(".srt","",dGc$assemblyId) #tilgiv mig sikora

res <- left_join(dGc, dBam, by = c("genusId", "assemblyId", "contigId")) %>%
    #left_join(dDamage, by = c("genusId", "assemblyId", "contigId")) %>%
    mutate(sampleId = sampleId, flag = "") %>%
    left_join(s1, by = "assemblyId") %>%
    #left_join(dKraken) %>%
    #select(sampleId, genusId, taxIdSpecies, taxNameSpecies, assemblyId:dam3pAvgDecay, krakenNClade:krakenKmerRank, flag) %>%
    filter(!is.na(nReads))



#idx <- res$dam5p >= 0.1 & res$dam5pAvgDecay < 0
idx1 <- res$aniRank100 < 2
idx2 <- res$coveragePRatio >= 0.8
#idx3 <- res$krakenKmerRank < 2 & !is.na(res$krakenKmerRank)

#res$flag[idx] <- "damage_0.1"
res$flag[idx1] <- "aniRank100_1"
res$flag[idx2] <- "coveragePRatio_0.8"

#res$flag[idx & idx1] <- "damage_0.1;aniRank100_1"
#res$flag[idx & idx2] <- "damage_0.1;coveragePRatio_0.8"
#res$flag[idx & idx3] <- "damage_0.1;krakenKmerRank_1"
res$flag[idx1 & idx2] <- "aniRank100_1;coveragePRatio_0.8"
#res$flag[idx1 & idx3] <- "aniRank100_1;krakenKmerRank_1"
#res$flag[idx2 & idx3] <- "coveragePRatio_0.8;krakenKmerRank_1"
#res$flag[idx & idx1 & idx2] <- "damage_0.1;aniRank100_1;coveragePRatio_0.8"
#res$flag[idx & idx1 & idx3] <- "damage_0.1;aniRank100_1;krakenKmerRank_1"
#res$flag[idx & idx2 & idx3] <- "damage_0.1;coveragePRatio_0.8;krakenKmerRank_1"
#res$flag[idx1 & idx2 & idx3] <- "aniRank100_1;coveragePRatio_0.8;krakenKmerRank_1"

#res$flag[idx & idx1 & idx2 & idx3] <- "damage_0.1;aniRank100_1;coveragePRatio_0.8;krakenKmerRank_1"
#res$flag[idx & idx1 & idx2 & idx3] <- "damage_0.1;aniRank100_1;coveragePRatio_0.8;krakenKmerRank_1"
write_tsv(res, file = out_tsv, na = "NaN")

#only include whole genome analysis
genomes <- res[which(res$contigId=="genome"),]
#get the best mapping based on coverage
bestmap <- cbind(c(sampleId),genomes[which(genomes$coverageP==max(genomes$coverageP)),c("assemblyId","coverageP","ani","coverageAvg","coveragePExp","coveragePRatio")])
write_tsv(bestmap, file = outBestmap, na = "NaN",col_names=FALSE)
#version of dist matrix which only includes references for which we have results
#increases robustness of script
#print(ani_rep)
#print(genomes$assemblyId)
ani_rep2<-ani_rep[which(row.names(ani_rep) %in% unique(genomes$assemblyId)),which(colnames(ani_rep) %in% unique(genomes$assemblyId))]
#print(ani_rep2)
#create dendogram from distance matrix
clustered_rep = hclust(as.dist(ani_rep2), method = "average") #compute hierarchical clustering
dend<-as.dendrogram(clustered_rep)
dend_data<-dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))

#include taxname into the dendogram for tick labels
ani_dat<-dend_data$labels
colnames(ani_dat)<-c("x","y","assemblyId")
ani_dat_species<-ani_dat %>% left_join(genomes[c("assemblyId","taxNameSpecies")])

#Position the dendogram labels
ref_pos_table <- with(
    ani_dat_species, 
    data.frame(y_center = x, assemblyId = as.character(assemblyId),taxNameSpecies=as.character(taxNameSpecies), height = 1))

# Limits for the vertical axes
ref_axis_limits <- with(
    ref_pos_table, 
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
    0.1 * c(-1, 1) # extra spacing: 0.1

plt_dendr <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_x_reverse(expand = c(0, 0)) + 
    scale_y_continuous(breaks = ref_pos_table$y_center, 
                       labels = ref_pos_table$taxNameSpecies, 
                       limits = ref_axis_limits,
                       expand = c(0, 0)) + 
    labs(x = "Distance", y = NULL, colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid.minor = element_blank())+ ggtitle(paste0("Reference phylogeny and mapping coverage for: \n",sampleId))
plt_dendr


genomes$assemblyId<-factor(genomes$assemblyId,levels=labels(dend)) #align genomes_plot with dendogram assembly id
refcov<-ggplot(data=genomes, aes(x=assemblyId, y=coverageP,fill=ani)) +
    #geom_point()+
    geom_bar(stat="identity",position="dodge")+
    labs(x="",y="Coverage") +
    scale_fill_continuous(low="blue", high="red") +
    scale_y_continuous(limits=c(0,1),expand=c(0,0)) + #limits=c(0,1),max(genomes$mqAvg)+1
    scale_x_discrete(position="top") +
    guides(fill=guide_legend(title="ANI")) + coord_flip() +
    theme_bw() + theme(panel.grid.minor = element_blank()) 
refcov

pdf(outPlot,width=17, height=20)
plot_grid(plt_dendr, refcov, align = "h", rel_widths = c(0.6, 1),nrow=1)
dev.off()

