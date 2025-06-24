library(readr)
library(janitor)
library(dplyr)
library(tidyr)
#Read the divsum file
kimura_div <- read_table("/repeats/HiTE/derLae1_hic.FINAL.fasta.divsum", 
                         col_types = cols(Div = col_integer(), 
                                          DNA = col_double(), X22 = col_skip()), 
#Use grep "Coverage for each repeat class and divergence (Kimura)" *.divsum -n to determine this number
                         skip = 1434)
#nombres is for plotting
nombres <- colnames(kimura_div)

#Read table with scaffold sizes chrom sizes
assembly <- read_delim("derLa1_hic.Final.chrom.sizes", delim = "\t", 
                       escape_double = FALSE, col_names = FALSE, 
                       col_types = cols(X2 = col_integer()), 
                       trim_ws = TRUE)
#Just to get the genome size
genome_size <- sum(assembly$X2)

#Make percentage of genome size
kimura_div <- kimura_div %>% mutate(across(c(-Div), \(x) 100*x/genome_size))

#Transform to long dataframe for graph
kimura_long <- kimura_div %>% pivot_longer(!Div, names_to = "Family", values_to = "percent")

library(stringr)
kimura_long <- kimura_long %>% mutate(
  Family = str_replace(Family, "I-Jockey", "Jockey"),
  Element = lapply(str_split(Family, "-"), first) %>% unlist)

library(ggplot2)
library(RColorBrewer)
kimura_long_short_names <- kimura_long %>% group_by(Div, Element) %>%
  summarise(Percent = sum(percent))

#Dont graph these elements that are not detected efficiently
families_2_plot <- kimura_long_short_names %>%
  filter(!(Element %in% c("tRNA", "snRNA", "LINE/Tad1", "rRNA", "SINE?", "Simple_repeat"))) %>% pull(Element) %>% unique()

#use scale_fill_hue(c = 50, direction =-1) to get colors around the color wheel, otherwise similar colors are plotted together and it is impossible to differentiate the elements
colors <- scale_fill_hue(c = 50, direction =-1)$palette(length(families_2_plot)) %>% as.list()
names(colors) <- families_2_plot
intersplot <- kimura_long_short_names %>%
  filter(Element %in% families_2_plot) %>% 
  ggplot() + geom_bar(aes(Div, Percent, fill = Element), stat = "identity") + 
  scale_fill_manual(values = colors) +
  theme_bw() + xlim(c(-1,50))

#DNA transposons
dna_elem <- str_detect(families_2_plot, "DNA*")

dnaplot <- kimura_long_short_names %>% 
  filter(Element %in% families_2_plot[dna_elem]) %>% 
  ggplot() + geom_bar(aes(Div, Percent, fill = Element), stat = "identity") + 
  #scale_fill_manual(values = colors[dna_elem]) +
  scale_fill_hue(c = 50, direction = -1) + 
  theme_bw() + xlim(c(-1,50))

#LINEs
LINEs <- str_detect(families_2_plot, "LINE")
lineplot <- kimura_long_short_names %>%
  filter(Element %in% families_2_plot[LINEs]) %>% 
  ggplot() + geom_bar(aes(Div, Percent, fill = Element), stat = "identity") + 
  scale_fill_hue(c = 50, direction = -1) + 
  theme_bw() + xlim(c(-1,50))

resto <- families_2_plot %in% c("LTR/ERVK", "LTR/ERV", "LTR/Gypsy", "LTR/Pao", "RC/Helitron", "SINE", "SINE/tRNA", "PLE/Naiad")

ltrsineplot <- kimura_long_short_names %>%
  filter(Element %in% families_2_plot[resto]) %>% 
  ggplot() + geom_bar(aes(Div, Percent, fill = Element), stat = "identity") + 
  scale_fill_hue(c = 50, direction = -1) + 
  theme_bw() + xlim(c(-1,50))

#Figure 3A-D
ggsave(plot = intersplot, filename = "hite_interspersed_repeat_landscape.pdf", height = 5, width = 10)
ggsave(plot = dnaplot, filename = "hite_dnaTE_repeat_landscape.pdf", height = 5, width = 10)
ggsave(plot = lineplot, filename = "hite_line_repeat_landscape.pdf", height = 5, width = 10)
ggsave(plot = ltrsineplot, filename = "hite_ltrsine_repeat_landscape.pdf", height = 5, width = 10)
