#analysis of motifs distributions 
library(tidyverse)
library(scales)
data <- read_csv("human_proteome_motifs_across_domains.csv")

### Also note.. some of the motifs here are subsets of the others, i.e. "DLL" is in "XL[LI]". Let's make a faceted plot without some of the duplicates, i.e. you'd exclude DLL but keep the motif that contains it. You should look into this manually, pick the specific ones to keep.



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

ggsave(file = "YXX_general_plot.pdf", YXX_general_plot)


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

ggsave(file = "YXX_domain_plot.pdf", YXX_domain_plot)


motif_domain_plot <- data %>%
  group_by(motif_type,sequence_id,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0) %>%
  filter(O > 0) %>%
  group_by(motif_type) %>%
  filter(motif_type != "DLL") %>%
  filter(motif_type != "LLX") %>%
  ggplot(aes(x = O, y = D)) + 
  scale_x_continuous("Motif Count in Ordered Regions", breaks = pretty_breaks()) +
  scale_y_continuous("Motif Count in Disordered Regions", breaks = pretty_breaks()) +
  labs(title = "Motif Counts where each Domain has at least one") +
  geom_hex(bins = 5) +
  stat_binhex(aes(label=..count..), geom="text", bins=5, colour="white") +
  facet_wrap(~motif_type, ncol = 3)
  
ggsave(file = "motif_domain_plot.pdf", motif_domain_plot)

endo_data <- read_csv("endocytosis_involved_proteins.csv")

data <- rename(data, UNIPROT_ID = sequence_id)

endo_proteins <- endo_data %>%
  select(UNIPROT_ID) %>%
  mutate(endo = UNIPROT_ID %in% data$UNIPROT_ID) %>%
  filter(endo == TRUE) %>%
  inner_join(data, by = c('UNIPROT_ID')) %>%
  group_by(UNIPROT_ID,motif_type,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0) %>%
  filter(O > 0) %>%
  group_by(motif_type)
  

non_endo_proteins <- data %>%
  select(UNIPROT_ID) %>%
  mutate(non_endo = UNIPROT_ID %in% endo_data$UNIPROT_ID) %>%
  filter(non_endo == FALSE) %>%
  inner_join(data, by = c('UNIPROT_ID')) %>%
  group_by(UNIPROT_ID,motif_type,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0) %>%
  filter(O > 0) %>%
  group_by(motif_type)


#compare proteins in endocytosis file to proteins not in endocytosis file(but in human proteome orginal file)
# Focus further on these proteins which are specifically involved in forming clathrin-coated pits for endocytosis: AP2A1, AP2A2, AP2B1, AP2S1, AP2M1, AAK1.You can use summarize or tally functions to determine the counts of motifs per structural domain in each of these proteins.