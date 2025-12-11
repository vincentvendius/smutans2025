####################
## plot plink pca ##
####################


## --------------------------------------------------------
## libraries

library(tidyverse)
library(ggrepel)

## --------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
eig_file <- args[1]
panel <- args[2]
sample_file <- args[3]
hl_file <- args[4]
out_file <- args[5]


## --------------------------------------------------------
## helpers

sample_info <- read_tsv(sample_file)
samples_hl <- scan(hl_file, what = character())

th <- theme_bw() +
    theme(
        panel.grid.major = element_line(
            linetype = "dotted",
            linewidth = 0.25
        ),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")
    )

s <- sample_info %>%
    distinct(group, color, fill, shape)

pal_c <- s$color
names(pal_c) <- s$group
pal_f <- s$fill
names(pal_f) <- s$group
pal_s <- s$shape
names(pal_s) <- s$group


## --------------------------------------------------------
## process data

d <- read_table(eig_file, col_names = FALSE)

n_pcs <- ncol(d) - 2
colnames(d)[1:2] <- c("fam_id", "sample_id")
colnames(d)[3:(n_pcs + 2)] <- paste("PC", 1:n_pcs, sep = "")

## var explained
v_file <- gsub("eigenvec", "eigenval", eig_file)
v1 <- scan(v_file)
v2 <- v1 / sum(v1)


## --------------------------------------------------------
## plot PCA

labs <- matrix(paste("PC", 1:n_pcs, sep = ""), nrow = 2)
l1 <- formatC(v2 * 100, format = "f", digits = 2)
mm <- matrix(paste(labs, " (", l1[1:n_pcs], "%)", sep = ""), nrow = 2)


#d1 <- d %>%
#    left_join(sample_info) %>%
#    mutate(group = factor(group,
#        levels = unique(sample_info$group)
#    ))
#
#d1_bg <- d1 %>%
#    filter(!sample_id %in% samples_hl)
#
#d1_hl <- d1 %>%
#    filter(sample_id %in% samples_hl)
#
#pdf(out_file,
#    width = 7,
#    height = 4.5
#)
#
#walk(seq_len(ncol(l2)), ~ {
#    p <- ggplot(d1, aes(
#        x = !!sym(labs[1, .x]),
#        y = !!sym(labs[2, .x])
#    ))
#
#    print(p +
#        geom_point(size = NA) +
#        geom_point(
#            aes(
#                color = group,
#                fill = group,
#                shape = group
#            ),
#            size = 1.5,
#            alpha = 0.9,
#            data = d1_bg
#        ) +
#        geom_point(
#            aes(
#                fill = group,
#		color = group,
#                shape = group
#            ),
#            size = 2,
#            alpha = 1,
#            data = d1_hl
#        ) +
#        geom_text_repel(
#            aes(
#                label = strain
#            ),
#            size = 1.5,
#            segment.size = 0.25,
#            segment.color = "grey",
#            data = d1_hl,
#            max.overlaps = Inf
#        ) +
#        xlab(l2[1, .x]) +
#        ylab(l2[2, .x]) +
#        scale_color_manual(values = pal_c) +
#        scale_fill_manual(values = pal_f) +
#        scale_shape_manual(values = pal_s) +
#        th)
#})
#dev.off()


d5 <- d %>%
    left_join(sample_info) %>%
    mutate(group = factor(group,
        levels = unique(sample_info$group)
    ))

d5_bg <- d5 %>%
    filter(!sample_id %in% samples_hl)

d5_hl <- d5 %>%
    filter(sample_id %in% samples_hl)

pdf(out_file,
    width = 7,
    height = 4.5
)

walk(seq_len(ncol(mm)), ~ {
#    p <- ggplot(d1, aes(
#    ))
    p <- ggplot(d5, aes(
        x = !!sym(labs[1, .x]),
        y = !!sym(labs[2, .x])
#        x = !!sym(mm[1, .x]),
#        y = !!sym(mm[2, .x])
    ))
    print(p +
        geom_point(
            aes(
                fill = group,
                shape = group,
                color = group
            ),
            size = 1.5,
            alpha = 0.6,
            data = d5_bg
        ) +
        geom_point(
            aes(
                fill = group,
		color = group,
                shape = group
            ),
            size = 2.5,
            alpha = 1,
            data = d5_hl
        ) +
        geom_text_repel(
            aes(
                label = strain
            ),
            size = 1.5,
            segment.size = 0.25,
            segment.color = "grey",
            data = d5_hl,
            max.overlaps = Inf
        ) + ggtitle("PCA of S. mutans",subtitle="High coverage ancient and modern assemblies colored by Chromopainter cluster")+
        scale_fill_manual(values=c("ancient_low_cov"="black","Cluster_A_ancient"="#ffce54","Cluster_A"="#ffce54",
                                   "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568","Cluster_D_ancient"="#a0d568",
                                   "Cluster_E"="#ed5564"),guide=NULL)+

        scale_color_manual(values=c("ancient_low_cov"="black","Cluster_A_ancient"="black","Cluster_A"="#ffce54",
                                    "Cluster_B"="#ac92eb","Cluster_C"="#4fc1e8","Cluster_D"="#a0d568","Cluster_D_ancient"="black",
                                    "Cluster_E"="#ed5564"),name="Chromopainter \nCluster") +
		    #scale_shape_manual(values=c(
		    #  "Cluster_A_ancient"=15,"Cluster_D_ancient"="#a0d568","Ancient Whole_tooth"=16,"Ancient Mixed"=17,"Ancient Bone"=18),na.value="#FFFFFF",
		    #name="Pre-Lundbeck Samples",na.translate = F)+
        scale_shape_manual(values = c(21,21,21,21,21,24,24,24),guide=NULL) +
        th)
})
dev.off()
