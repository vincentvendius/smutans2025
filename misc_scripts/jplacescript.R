####################
## plot EPA trees ##
####################


## --------------------------------------------------------------------------------
## libraries
#Rscript src/plotEpa.R {input.jplace} {wildcards.sample}.epa {output.pdf}
suppressMessages(library(tidyverse))
suppressMessages(library(viridis))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(ape))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
jplace <- args[1]
sampleId <- args[2]
outPdf <- args[3]


## --------------------------------------------------------------------------------
## read and process placement data

d <- read.jplace(jplace)
d@data$llr <- 0
d@data$llr[match(d@placements$node, d@data$node)] <- d@placements$like_weight_ratio
d@data$hasPlacement <- ifelse(d@data$llr > 0, "yes", "no")

d@phylo <- root(d@phylo, "Akhmeta_7_GEO_MH607143_Vani_2010")


## --------------------------------------------------------------------------------
## plot trees

## placements
idxL <- c("solid", "dashed", NA)
names(idxL) <- c("yes", "no", NA)

h <- length(d@phylo$tip.label) %/% 10

pdf(outPdf, width = 8, height = h)
p <- ggtree(d, aes(x = x,
                   y = y,
                   color = llr,
                   alpha = llr,
                   size = llr,
                   linetype = hasPlacement))
print(p +
        geom_tiplab(size = 1.5,
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
        ggtitle(sampleId))

## graft
tr <- read.tree(gsub("jplace$", "newick", jplace))
tr <- root(tr, "Akhmeta_7_GEO_MH607143_Vani_2010")

idxX <- match(sampleId, tr$tip.label)
idx1 <- match(idxX, tr$edge[,2])
tr$edge.length[idx1] <- 0.002

p <- ggtree(tr, aes(x = x,
                    y = y),
            size = 0.25,
            colour = "grey")
print(p +
        geom_tiplab(aes(color = label == sampleId),
                    size = 1.5,
                    hjust = 0,
                    vjust = 0.5,
                    alpha = 1,
                    show.legend = FALSE) +
        scale_color_manual(values = c("grey", "black")) +
        theme_tree() +
        theme(legend.position = "right") +
        guides(fill = guide_legend(override.aes = list(size = 2.5))) +
        ggtitle(sampleId))
dev.off()