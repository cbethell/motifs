#analysis of motifs distributions 
library(tidyverse)
library(scales)
library(cowplot)
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

#plot showing count of motifs in disordered regions for each motif type
motif_disordered_plot <- data %>%
  group_by(motif_type,sequence_id,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0, O > 0) %>%
  group_by(motif_type) %>%
  filter(motif_type != "DLL") %>%
  filter(motif_type != "LLX") %>%
  ggplot(aes(x = motif_type, y = D)) + 
  theme(axis.text.x = element_text(size = 5)) +
  geom_hex(bins = 10) +
  stat_binhex(aes(label=..count..), geom="text", bins=10, colour="white",size = 3)

ggsave(file = "motif_disordered_plot.pdf", motif_disordered_plot)

#plot showing count of motifs in ordered regions for each motif type
motif_ordered_plot <- data %>%
  group_by(motif_type,sequence_id,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0, O > 0) %>%
  group_by(motif_type) %>%
  filter(motif_type != "DLL") %>%
  filter(motif_type != "LLX") %>%
  ggplot(aes(x = motif_type, y = O)) + 
  theme(axis.text.x = element_text(size = 5)) +
  geom_hex(bins = 10) +
  stat_binhex(aes(label=..count..), geom="text", bins=10, colour="white",size = 2.5)

ggsave(file = "motif_ordered_plot.pdf", motif_ordered_plot)

#cowplot to compare the disordered and ordered plots
combined_domain_plot <- plot_grid(motif_disordered_plot, motif_ordered_plot, labels = "auto")
ggsave(file = "combined_domain_plot.pdf", combined_domain_plot)

#read in endocytosis protein data
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
  filter(D > 0, O > 0) %>%
  group_by(motif_type) %>%
  ungroup()
  

non_endo_proteins <- data %>%
  select(UNIPROT_ID) %>%
  mutate(non_endo = UNIPROT_ID %in% endo_data$UNIPROT_ID) %>%
  filter(non_endo == FALSE) %>%
  inner_join(data, by = c('UNIPROT_ID')) %>%
  group_by(UNIPROT_ID,motif_type,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0, O > 0) %>%
  group_by(motif_type) %>%
  ungroup()

definite_active_motifs <- data %>%
  filter(UNIPROT_ID == "O95782" | UNIPROT_ID == "O94973" | UNIPROT_ID == "P63010" | UNIPROT_ID == "P53680" | UNIPROT_ID == "Q96CW1" | UNIPROT_ID == "Q2M2I8
") %>%
  group_by(UNIPROT_ID,domain_type) %>%
  summarize(count=n()) %>%
  spread(domain_type, count)%>%
  ggplot(aes(x = O, y = D)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label=UNIPROT_ID))

ggsave(file = "active_endo_motifs.pdf", definite_active_motifs)
