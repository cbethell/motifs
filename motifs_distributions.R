#analysis of motifs distributions 
data <- read.csv("human_proteome_motifs_across_domains.csv",sep=",",stringsAsFactors=FALSE,strip.white=TRUE,fill=TRUE,header=TRUE)

#focus only on the motif "YXX[LIMFV]
YXX = subset(data, motif_type == "YXX[LIMFV]")
total_ordered_count <- sum(YXX$domain_type == "O")
total_disordered_count <- sum(YXX$domain_type == "D")

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)

YXX %>%
  drop_na() %>%
  summary()

#### TODO: Try to consolidate the code that produces `total_count_YXX` into a single stream of piping using the dplyr function group_by

#find the motifs that are in the ordered regions
ordered_YXX <- YXX %>%
  filter(domain_type == "O")

#find the motifs that are in the disordered regions
disordered_YXX <- YXX %>%
  filter(domain_type == "D")

#number of motifs found in ordered regions for each protein
ordered_id_count <- ordered_YXX %>%
  group_by(sequence_id) %>%
  summarize(count=n())

#number of motifs found in disordered regions for each protein 
disordered_id_count <- disordered_YXX %>%
  group_by(sequence_id) %>%
  summarize(count=n())

#shows ordered vs disordered for proteins that had motifs in both regions but missing those that are in disordered and not in ordered (need to fix)
total_count_YXX <- ordered_id_count %>% 
  left_join(disordered_id_count, by = "sequence_id") %>%
  as_tibble() %>%
  rename(o_count = count.x) %>%
  rename(d_count = count.y) %>%
  mutate(d_count = replace_na(d_count, 0)) %>%
  mutate(sig = o_count - d_count >= abs(-5)) #look for most significant 


### TODO: abs(7) is 7... what are you trying to take the absolute value of? 

total_count_YXX %>%
  filter(sig == TRUE) %>%
  ggplot(aes(x = o_count, y = d_count)) +
  geom_point()

### TODO: use filter rather than logical indexing. The first goal here should be to make a scatterplot, not a line plot. X axis = number of ordered motifs per protein, Y axis = number of disordered motifs per protein. At this stage we don't care about the specific proteins. Note how nothing is really meaningfully colored in this plot
##ggplot(data = sig[sig$sig==TRUE,], mapping = aes(x = o_count, y = d_count, color = sequence_id, group = 1)) + 
##  geom_line()
#find proteins that exist in only one of the datasets (does not overlap)
