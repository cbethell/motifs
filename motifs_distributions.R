#analysis of motifs distributions 
data <- read.csv("human_proteome_motifs_across_domains.csv",sep=",",stringsAsFactors=FALSE,strip.white=TRUE,fill=TRUE,header=TRUE)

#focus only on the motif "YXX[LIMFV]
YXX = subset(data, motif_type == "YXX[LIMFV]")
total_ordered_count <- sum(YXX$domain_type == "O")
total_disordered_count <- sum(YXX$domain_type == "D")

## SJS Libraries must always be loaded at the top of a file. Move to top, AND redo all above code using the tidyverse.
library(tidyverse)
### SJS don't load these libraries - they are loaded in tidyverse already. this is what tidyverse is
#library(tidyr)
#library(dplyr)
#library(ggplot2)


### SJS nothing gets saved here... What is this line doing?
YXX %>%
  drop_na() %>%
  summary()

#### TODO: Try to consolidate the code that produces `total_count_YXX` into a single stream of piping using the dplyr function group_by


#####################################################################
### SJS lines of code in this comment chunk can ALL be done in a SINGLE piped command using group_by. This removes the need for joining below, see my comments there

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
##########################################################################


#shows ordered vs disordered for proteins that had motifs in both regions but missing those that are in disordered and not in ordered (need to fix)
total_count_YXX <- ordered_id_count %>% 
  left_join(disordered_id_count, by = "sequence_id") %>%
  as_tibble() %>%
  rename(o_count = count.x) %>%  ### SJS the fact that you have two counts here tells me there may be an issue with your joining. You can address this entirely by only having ONE line of piping code that collects motifs per domain. See comment above in the comment "chunk". In other words, you don't actually need to join anything.
  rename(d_count = count.y) %>%
  mutate(d_count = replace_na(d_count, 0)) %>%  
  mutate(sig = o_count - d_count >= abs(-5)) #look for most significant 
### SJS What is negative -5? This is hardcoded, and I don't understand where the number -5 comes from. 
## This is also NOT how statistical significance is assessed. This is highly arbitrary, definitely remove it. Also, it means you aren't actually visualizing all motifs - just those at an arbitrary cutoff.

### TODO: abs(7) is 7... what are you trying to take the absolute value of? 


### This plot is a great start, but look at it a bit - notice how it's almost impossible to see any information? That's because points are all on top of each other
### One solution is geom_hex: https://ggplot2.tidyverse.org/reference/geom_hex.html . This makes a heatmap which can be colored by number of proteins at each point. Try that out for Monday.
total_count_YXX %>%
  filter(sig == TRUE) %>%
  ggplot(aes(x = o_count, y = d_count)) +
  geom_point()

### TODO: use filter rather than logical indexing. The first goal here should be to make a scatterplot, not a line plot. X axis = number of ordered motifs per protein, Y axis = number of disordered motifs per protein. At this stage we don't care about the specific proteins. Note how nothing is really meaningfully colored in this plot
##ggplot(data = sig[sig$sig==TRUE,], mapping = aes(x = o_count, y = d_count, color = sequence_id, group = 1)) + 
##  geom_line()
#find proteins that exist in only one of the datasets (does not overlap)
