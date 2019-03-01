#analysis of motifs distributions 
data <- read.csv("human_proteome_motifs_across_domains.csv",sep=",",stringsAsFactors=FALSE,strip.white=TRUE,fill=TRUE,header=TRUE)
YXX = subset(data, motif_type == "YXX[LIMFV]")
total_ordered_count <- sum(YXX$domain_type == "O")
total_disordered_count <- sum(YXX$domain_type == "D")

library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)

YXX %>%
  drop_na() %>%
  summary()

ordered_YXX <- YXX %>%
  filter(domain_type == "O")

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

#shows ordered vs disordered for proteins that had motifs in both regions
total_count_YXX <- ordered_id_count %>% 
  join(disordered_id_count, type = "inner", by = "sequence_id")

names(total_count_YXX)[2] <- "o_count" 
names(total_count_YXX)[3] <- "d_count" 
total_count_YXX <- as.tibble(total_count_YXX)

#too large
ggplot(data = total_count_YXX, mapping = aes(x = o_count, y = d_count, color = sequence_id)) +
  geom_line() 
