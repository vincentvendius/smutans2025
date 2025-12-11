"~/PhD/datafiles/mutans_bins_1kb.txt" %>%  
  read_tsv()  %>% 
  mutate('100-red' = 100 - percent_redundancy) %>% 
  pivot_longer(cols = c("total_length","num_contigs", "N50", "GC_content", "percent_completion", "100-red", "percent_redundancy")) %>% 
  mutate(name = factor(name, levels = c("total_length","num_contigs", "N50", "percent_completion", "100-red", "percent_redundancy", "GC_content"))) -> data

data %>% 
  ggplot(aes(x = `Assembly-strategy`,
             y = value)) +
  facet_grid(name ~ ., scale="free",space="free_x") +
  geom_boxplot(outlier.colour = NA, alpha = 0.2,
               position = position_dodge(width = 0.8),
               aes(group = interaction(`Assembly-strategy`, `Binning-strategy`))) +
  geom_point(aes(fill=Sample, group = interaction(Sample, `Assembly-strategy`, `Binning-strategy`), 
                 color = Sample, shape = `Binning-strategy`), 
             position = (position = position_jitterdodge(jitter.width = 0.1, dodge.width =  0.8))) +
  ylab("assembly - binning strategy") +
  theme_bw()  -> fig

fig

data %>% 
  ggplot(aes(x = `Binning-strategy`,
             y = value)) +
  facet_grid(name ~ ., scale="free",space="free_x") +
  geom_boxplot(outlier.colour = NA, alpha = 0.2,
               position = position_dodge(width = 0.8),
               aes(group = interaction(`Assembly-strategy`, `Binning-strategy`))) +
  geom_point(aes(fill=Sample, group = interaction(Sample, `Assembly-strategy`, `Binning-strategy`), 
                 color = Sample, shape = `Binning-strategy`), 
             position = (position = position_jitterdodge(jitter.width = 0.1, dodge.width =  0.8))) +
  ylab("assembly - binning strategy") +
  theme_bw()  -> fig

fig

data %>% 
  ggplot(aes(x = Sample,
             y = value)) +
  facet_grid(name ~ ., scale="free",space="free_x") +
  # geom_boxplot(outlier.colour = NA, alpha = 0.2,
  #              position = position_dodge(width = 0.8),
  #              aes(group = interaction(Sample))) +
  geom_point(aes(fill=Sample,
                 # group = interaction(Sample, `Assembly-strategy`, `Binning-strategy`), 
                 color = `Assembly-strategy`, shape = `Binning-strategy`), 
             position = (position = position_jitterdodge(jitter.width = 0.1, 
                                                         dodge.width =  0.8))) +
  ylab("assembly & binning strategy") + ggtitle("Mutans bin statistics using 1kb contig threshold") +
  theme_bw()  -> fig

fig

mag_stats <- read.delim("~/PhD/datafiles/mag_stats.tsv")