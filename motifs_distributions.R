#analysis of motifs distributions 
library(tidyverse)
library(scales)
data <- read_csv("human_proteome_motifs_across_domains.csv")

#focus only on the motif "YXX[LIMFV] and find the number of motifs in each domain type: ordered and disordered for each protein 

YXX <- data %>% 
  filter(motif_type == "YXX[LIMFV]") %>%
  group_by(sequence_id,domain_type) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  mutate(sum = D + O)

#heatmap of atleast one motif in the ordered or disordered region 
YXX_general_plot <- YXX %>%
  filter(sum > 0) %>%
  ggplot(aes(x = O, y = D)) + 
  scale_x_continuous("Motif Count in Ordered Regions", breaks = pretty_breaks()) +
  scale_y_continuous("Motif Count in Disordered Regions") +
  labs(title = "Motif Counts") +
  geom_hex(bins = 10) +
  stat_binhex(aes(label=..count..), geom="text", bins=10, colour="white")

#heatmap of atleast one motif in each domain 
YXX_domain_plot <- YXX %>%
  filter(D > 0) %>%
  filter(O > 0) %>%
  ggplot(aes(x = O, y = D)) + 
  scale_x_continuous("Motif Count in Ordered Regions", breaks = pretty_breaks()) +
  scale_y_continuous("Motif Count in Disordered Regions") +
  labs(title = "Motif Counts where each Domain has at least one") +
  geom_hex(bins = 10) +
  stat_binhex(aes(label=..count..), geom="text", bins=10, colour="white") 

