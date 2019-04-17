library(tidyverse)



data <- read_csv("human_proteome_motifs_across_domains.csv")

filter_count <- function(data, motif) {
  data %>%
  group_by(sequence_id,domain_type,motif_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>%
  arrange(desc(D)) %>%
  filter(motif_type == !!motif)
}

data %>% 
    filter_count("NPF") %>%
    filter(D > 0) -> d2

#### maaaaybe not darkorchid1? 
### Try to figure out other color schemes too, because it's more fun that way.
d2 %>% 
    group_by(D, O) %>% 
    tally() %>%
    mutate(xmin = O-0.5, ymin = D-0.5, xmax = xmin + 1, ymax = ymin+1) %>%
    ggplot() + geom_rect(color = "white", aes(xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, fill = n)) + 
    geom_text(color = "darkorchid1", aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label = n)) + 
    scale_x_continuous(name = "Ordered Motif Count", breaks=0:3) + 
    scale_y_continuous(name = "Disordered Motif Count",breaks=1:6) + 
    ggtitle("Number of proteins") 