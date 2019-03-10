#analysis of motifs distributions 
library(tidyverse)

data <- read_csv("human_proteome_motifs_across_domains.csv")

#focus only on the motif "YXX[LIMFV] and find the number of motifs in each domain type: ordered and disordered for each protein 

YXX <- data %>% 
  filter(motif_type == "YXX[LIMFV]") %>%
  group_by(sequence_id,domain_type) %>%
  summarize(count=n()) 

### One solution is geom_hex: https://ggplot2.tidyverse.org/reference/geom_hex.html . This makes a heatmap which can be colored by number of proteins at each point. Try that out for Monday.

YXX %>%
  as_tibble() %>% 
  ggplot(aes(x = sequence_id, y = domain_type)) +
  geom_hex(aes(color = count)) 



#cluster to include sequence ids (since amount of observations seems to be too large)?

#find proteins that exist in only one of the datasets (does not overlap)
