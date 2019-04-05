#analysis of motifs distributions 
library(tidyverse)
library(scales)
library(cowplot)
library(ggrepel)
data <- read_csv("human_proteome_motifs_across_domains.csv")

#focus only on the motif "YXX[LIMFV] and find the number of motifs in each domain type: ordered and disordered for each protein 

#function to filter and count the number of motifs in the regions of each motif expression
filter_count <- function(data, motif) {
  data %>%
  group_by(sequence_id,domain_type,motif_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>% 
  arrange(desc(D)) %>%
  filter(motif_type == !!motif)
}

#function to plot the counts in the regions of each motif expression
filter_plot <- function(data,motif) {
  data %>%
  filter(D > 0) %>%
  ggplot(aes(x = O, y = D)) + 
  scale_x_continuous("Motif Count in Ordered Regions", breaks = pretty_breaks(), limits=c(0,10)) +
  scale_y_continuous("Motif Count in Disordered Regions", limits=c(0,18)) +
  labs(title = paste0("Motif Counts in ",motif," Expression")) +
  geom_hex(bins = 15) +
  geom_abline(color="red") + 
  stat_binhex(aes(label=..count..), geom="text", bins=15, colour="white", size=3.5) 
}
#filter for YXX[LIMFV] motif expression
YXX <- filter_count(data, quo("YXX[LIMFV]"))

#heatmap of atleast one motif in the disordered region of YXX[LIMFV]
YXX_domain_plot <- filter_plot(YXX, "YXX[LIMFV]")
ggsave(file = "YXX_plot.pdf", YXX_domain_plot, width=7, height=6)


#filter for NPF motif expression
NPF <- filter_count(data, quo("NPF"))

#heatmap of atleast one motif in the disordered region of NPF
NPF_domain_plot <- filter_plot(NPF, "NPF")
ggsave(file = "NPF_plot.pdf", NPF_domain_plot, width=7, height=6)

#filter for X[DE]XXXL[LI] motif expression
XD <- filter_count(data, quo("X[DE]XXXL[LI]")) 

#heatmap of atleast one motif in the disordered region of X[DE]XXXL[LI]
XD_domain_plot <- filter_plot(XD, "X[DE]XXXL[LI]")
ggsave(file = "XD_plot.pdf", XD_domain_plot, width=7, height=6)

#filter for NPXY motif expression
NPXY <- filter_count(data, quo("NPXY"))

#heatmap of atleast one motif in the disordered region of NPXY
NPXY_domain_plot <- filter_plot(NPXY, "NPXY")
ggsave(file = "NPXY_plot.pdf", NPXY_domain_plot, width=7, height=6)

#plot showing count of motifs in disordered regions for each motif type
motif_disordered_plot <- data %>%
  group_by(motif_type,sequence_id,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0) %>%
  filter(motif_type %in% c("NPF","NPXY","X[DE]XXXL[LI]","YXX[LIMFV]")) %>%
  ggplot(aes(x = motif_type, y = D)) + 
  ylim(0,18) +
  labs(title = ("Motif Counts in the Disordered Regions of each Motif Expression")) +
  theme(axis.text.x = element_text(size = 5),plot.title = element_text(size = 8)) +
  geom_hex(bins = 10,aes(colour = factor(motif_type))) +
  stat_binhex(aes(label=..count..), geom="text", bins=10, colour="white",size = 3)

ggsave(file = "motif_disordered_plot.pdf", motif_disordered_plot)

#plot showing count of motifs in ordered regions for each motif type
motif_ordered_plot <- data %>%
  group_by(motif_type,sequence_id,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  filter(D > 0) %>%
  group_by(motif_type) %>%
  filter(motif_type %in% c("NPF","NPXY","X[DE]XXXL[LI]","YXX[LIMFV]")) %>%
  ggplot(aes(x = motif_type, y = O)) + 
  labs(title = ("Motif Counts in the Ordered Regions of each Motif Expression")) +
  theme(axis.text.x = element_text(size = 5),plot.title = element_text(size = 8)) +
  geom_hex(bins = 10,aes(colour = factor(motif_type))) +
  stat_binhex(aes(label=..count..), geom="text", bins=10, colour="white",size = 3)

ggsave(file = "motif_ordered_plot.pdf", motif_ordered_plot)

#cowplot to compare the disordered and ordered plots
combined_domain_plot <- plot_grid(motif_disordered_plot, motif_ordered_plot, labels = "auto")
ggsave(file = "combined_domain_plot.pdf", combined_domain_plot, width = 10)

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
  filter(D > 0) %>%
  group_by(motif_type) %>%
  filter(motif_type %in% c("NPF","NPXY","X[DE]XXXL[LI]","YXX[LIMFV]")) %>%
  arrange(desc(D)) %>%
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
  filter(motif_type %in% c("NPF","NPXY","X[DE]XXXL[LI]","YXX[LIMFV]")) %>%
  arrange(desc(D)) %>%
  ungroup()


definite_active_motifs <- data %>%
  filter(UNIPROT_ID == "O95782" | UNIPROT_ID == "O94973" | UNIPROT_ID == "P63010" | UNIPROT_ID == "P53680" | UNIPROT_ID == "Q96CW1" | UNIPROT_ID == "Q2M2I8") %>%
  inner_join(endo_data, by = c('UNIPROT_ID')) %>%
  group_by(UNIPROT_ID,UNIPROT_name,motif_type,domain_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0))%>%
  group_by(motif_type) %>%
  filter(motif_type %in% c("NPF","X[DE]XXXL[LI]","YXX[LIMFV]", "NPXY")) 

ggsave(file = "active_endo_motifs.pdf", definite_active_motifs)


