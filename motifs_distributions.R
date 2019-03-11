#analysis of motifs distributions 
library(tidyverse)
data <- read_csv("human_proteome_motifs_across_domains.csv")

#focus only on the motif "YXX[LIMFV] and find the number of motifs in each domain type: ordered and disordered for each protein 

YXX <- data %>% 
  filter(motif_type == "YXX[LIMFV]") %>%
  group_by(sequence_id,domain_type) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  mutate(o_count = ifelse(domain_type == "O", count, 0)) %>%
  mutate(d_count = ifelse(domain_type == "D", count, 0))

YXX %>%
  as_tibble() %>% 
  ggplot(aes(x = domain_type == "O", y = domain_type == "D")) +
  geom_hex(aes(colour = o_count, fill = d_count))

YXX %>%
  as_tibble() %>%
  ggplot(aes(x = o_count, y = d_count)) + 
  geom_hex(aes(colour = domain_type))

#cluster to include sequence ids (since amount of observations seems to be too large)?

#find proteins that exist in only one of the datasets (does not overlap)
