#analysis of motifs distributions 
library(tidyverse)
library(scales)
library(cowplot)
library(ggrepel)
library(ggpubr)
data <- read_csv("human_proteome_motifs_across_domains.csv")


#function to filter and count the number of motifs in the regions of each motif expression
data %<>%
  group_by(sequence_id,domain_type,motif_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>% 
  filter (D > 0)

#function to plot the counts in the regions of each motif expression
filter_plot <- function(data,motif,numx,numy) {
  data %>%
    group_by(D,O) %>%
    tally() %>%
    mutate(xmin = O-0.5, ymin = D-0.5, xmax = xmin + 1, ymax = ymin+1) %>%
    ggplot() + 
    geom_rect(color = "white", aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = n)) + 
    geom_text(color = "white", aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label = n)) + 
    scale_x_continuous("Ordered Motif Count", breaks = 0:numx) +
    scale_y_continuous("Disordered Motif Count", breaks = 0:numy) +
    grids(linetype = "dashed") +
    scale_fill_continuous(name = "Protein Count")
}

filter_STplot <- function(data,motif,number,numx,numy) {
  data %>%
    ggplot(aes(x = O, y = D)) + 
    scale_x_continuous("Ordered Motif Count", breaks = pretty_breaks(n = numx), limits=c(0,numx)) +
    scale_y_continuous("Disordered Motif Count", breaks = pretty_breaks(n = numy), limits=c(0,numy)) +
    geom_hex(bins = number) +
    geom_abline(color="green") + 
    stat_binhex(aes(label=..count..), geom="text", bins=number, colour="white", size=2.5) +
    scale_fill_gradient(high = "black", low = "red",name = "Frequency",na.value=NA) +
    grids(linetype = "dashed")
}
#filter for YXX[LIMFV] motif expression
YXX <- filter(data, motif_type == "YXX[LIMFV]")


#heatmap of atleast one motif in the disordered region of YXX[LIMFV]
YXX_domain_plot <- filter_plot(YXX, "YXX[LIMFV]", 10, 18) + geom_abline(color = "green")
ggsave(file = "YXX_plot.pdf", YXX_domain_plot, width=7, height=6)


#filter for NPF motif expression
NPF <- filter(data, motif_type == "NPF")

#heatmap of atleast one motif in the disordered region of NPF
NPF_domain_plot <- filter_plot(NPF, "NPF", 3, 5) + geom_abline(color = "green")
ggsave(file = "NPF_plot.pdf", NPF_domain_plot, width=7, height=6)


#filter for X[DE]XXXL[LI] motif expression
XD <- filter(data, motif_type == "X[DE]XXXL[LI]")

#heatmap of atleast one motif in the disordered region of X[DE]XXXL[LI]
XD_domain_plot <- filter_plot(XD, "X[DE]XXXL[LI]", 10, 18) + geom_abline(color = "green")
ggsave(file = "XD_plot.pdf", XD_domain_plot, width=7, height=6)

#filter for NPXY motif expression
NPXY <- filter(data, motif_type == "NPXY")

#heatmap of atleast one motif in the disordered region of NPXY
NPXY_domain_plot <- filter_plot(NPXY, "NPXY", 10, 18)
ggsave(file = "NPXY_plot.pdf", NPXY_domain_plot, width=7, height=6)


FXDX <- filter(data, motif_type == "FXDX[FILMV]")
FXDX_domain_plot <- filter_plot(FXDX, "FXDX[FILMV]", 5, 5) + geom_abline(color = "green")
ggsave(file = "FXDX_plot.pdf", FXDX_domain_plot, width=7, height=6)

DX <- filter(data, motif_type == "DX[FW]")
DX_domain_plot <- filter_plot(DX, "DX[FW]", 8, 20) + geom_abline(color = "green")
ggsave(file = "DX_plot.pdf", DX_domain_plot, width=7, height=6)

ST <- filter(data, motif_type == "[ST]XXXX[LI]")
ST_domain_plot <- filter_STplot(ST, "[ST]XXXX[LI]", 15, 10, 30) 
ggsave(file = "ST_plot.pdf", ST_domain_plot, width=7, height=6)


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

ggsave(file = "motif_disordered_plot.pdf", motif_disordered_plot,width=7, height=6)

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

ggsave(file = "motif_ordered_plot.pdf", motif_ordered_plot,width=7, height=6)

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
  filter(motif_type %in% c("DX[FW]","X[DE]XXXL[LI]","YXX[LIMFV]", "NPXY","[ST]XXXX[LI]","FXDX[FILMV]", "NPF")) %>%
  ggplot(aes(x = O, y = D)) + 
  scale_x_continuous("Motif Count in Ordered Regions", breaks = pretty_breaks()) +
  scale_y_continuous("Motif Count in Disordered Regions", breaks = pretty_breaks()) +
  geom_point(stat = "identity") +
  geom_abline(color="red") + 
  geom_text_repel(aes(label=UNIPROT_name)) +
  facet_wrap(~motif_type, ncol = 2, scales = "free")

ggsave(file = "active_endo_motifs.pdf", definite_active_motifs,width=7, height=6)


NPF_count <- NPF %>%
  count(D > O)
NPF_sum <- sum(NPF_count$`D > O` == TRUE)

YXX_count <- YXX %>%
  count(D > O)
YXX_sum <- sum(YXX_count$`D > O` == TRUE)

NPXY_count <- NPXY %>%
  count (D > O)
NPXY_sum <- sum(NPXY_count$`D > O` == TRUE)

XD_count <- XD %>%
  count (D > O)
XD_sum <- sum(XD_count$`D > O` == TRUE)

FXDX_count <- FXDX %>%
  count (D > O)
FXDX_sum <- sum(FXDX_count$`D > O` == TRUE)

DX_count <- DX %>%
  count (D > O)
DX_sum <- sum(DX_count$`D > O` == TRUE)

ST_count <- ST %>%
  count (D > O)
ST_sum <- sum(ST_count$`D > O` == TRUE)
