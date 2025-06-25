library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(plyranges)
library(ggbio)
library(janitor)
library(stringr)
library(purrr)
library(readr)

#Handy table of scaffold lengths
crom_sizes <- read.table("../csvs/cromosome_sizes.csv",
           col.names = c("scaffold", "length"))

#lengths to assembly
crom_sizes <- mutate(crom_sizes, start_in_ass = cumsum(lag(length, default = 1)),
           end_in_ass = length + start_in_ass - 1)

prefijo <- "../Alignments/toSelfgap/self_"

#Only the scaffolds that were not 100% repeat were aligned
nonrepeats <- read.table("../Alignments/Scaffolds_NOT100p_repeatslowercase.txt", header = F)
colnames(nonrepeats) <- c("scaffold")
#nonrepeats <- crom_sizes  %>% mutate(scaffold = toupper(scaffold)) %>% right_join(nonrepeats)
nonrepeats <- crom_sizes  %>% right_join(nonrepeats)

#laevew and laevep just short for whole (assembly) and parts (individual scaffolds)
alignments <- nonrepeats %>% tibble() %>% mutate(dots = map(scaffold, \(x) read_tsv(paste0(prefijo,x,".csv"), col_names = c("laevew", "laevep"), skip = 1)))

#Filter no rows
alignments <- alignments %>% mutate(filas = unlist(map(dots, nrow))) %>% filter(filas > 0)

alignments <- alignments %>% mutate(dots = map(dots, \(x) remove_empty(x, "rows")))

alignments <- alignments %>% mutate(dots = map(dots, \(x) mutate(x, ort_block_id = (row_number() - 1) %/% 2, axis = (row_number()-1) %% 2)))


alignments <- alignments %>% mutate(dots = map(dots, \(x) pivot_wider(x, names_from = axis, values_from = c(laevew, laevep))))

alignments <- alignments %>% unnest(dots)

alignments <- alignments %>% mutate(laeve_ass_0 = laevep_0 + start_in_ass -1, 
	laeve_ass_1 = laevep_1 + start_in_ass -1, width = sqrt((laevep_0 - laevep_1)^2 + (laevew_0 - laevew_1)^2))

minalsize <- 200
p <- filter(alignments, width >= minalsize) %>% ggplot() + geom_point(aes(x= laeve_ass_0, y = laevew_0, alpha = width, col = log(width)), size = 0.1) + scale_color_gradient2(low = "lightgray", high = "black", mid = "gray", space = "Lab", midpoint = log(10), limit = c(0, log(500))) + theme_bw()


endoflaeveass <- crom_sizes$end_in_ass[31]

p <- p + geom_hline(yintercept = endoflaeveass, col = "blue") + geom_vline(xintercept = endoflaeveass, col = "blue")


crombound <- crom_sizes$start_in_ass[2:31]

p <- p + geom_hline(yintercept = crombound, col = "red", linetype = "dashed") + geom_vline(xintercept = crombound, col = "red", linetype = "dashed")

p <- filter(alignments, width >= minalsize) %>% ggplot() + geom_bin2d(aes(x= laeve_ass_0, y = laevew_0), bins = 3000) + theme_bw() + scale_fill_gradient2(low = "lightgray", high = "red", mid = "blue", space = "Lab", midpoint = 2000)
