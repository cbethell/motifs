#analysis of motifs distributions 
library(tidyverse)
data <- read_csv("human_proteome_motifs_across_domains.csv")

#focus only on the motif "YXX[LIMFV] and find the number of motifs in each domain type: ordered and disordered for each protein 

YXX <- data %>% 
  filter(motif_type == "YXX[LIMFV]") %>%
  group_by(sequence_id,domain_type) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0))

YXX_plot <- YXX %>%
  as_tibble() %>%
  ggplot(aes(x = O, y = D)) + 
  geom_hex()

#cluster to include sequence ids (since amount of observations seems to be too large)?

#find proteins that exist in only one of the datasets (does not overlap)
