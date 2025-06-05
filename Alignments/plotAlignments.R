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
library(hexbin)

#faSize -detailed /genome/derLae1_hic.FINAL.fasta > /csvs/cromosome_sizes.csv
crom_sizes <- read.table("../csvs/cromosome_sizes.csv",
           col.names = c("scaffold", "length"))

#Read the length of the individual scaffolds, transform them into an assembly-wide coordinate system.
#More robust would be to first arrange crom_sizes by length in descending order, but this is not necessary as the numbering was assigned by size by 3D-DNA
crom_sizes <- mutate(crom_sizes, start_in_ass = cumsum(lag(length, default = 1)),
           end_in_ass = length + start_in_ass - 1)

#Read the names of the scaffolds that are not completely covered by repeats and were aligned
nonrepeats <- read.table("../Alignments/Scaffolds_NOT100p_repeats.txt", header = F)
colnames(nonrepeats) <- c("scaffold")

nonrepeats <- crom_sizes  %>% mutate(scaffold = toupper(scaffold)) %>% right_join(nonrepeats)

#Use map to read each csv file with the alignment segment information
alignments <- nonrepeats %>% tibble() %>% mutate(dots = map(scaffold, \(x) read_tsv(paste0("../Alignments/toFulica/",x,".csv"), col_names = c("fulica", "laeve"), skip = 1)))

#If there was no significant alignment for a given scaffold, it is empty, remove
alignments <- alignments %>% mutate(filas = unlist(map(dots, nrow))) %>% filter(filas > 0)

#Remove the NA rows
alignments <- alignments %>% mutate(dots = map(dots, \(x) remove_empty(x, "rows")))

#Just transformations to keep track of each block and position between each row.
alignments <- alignments %>% mutate(dots = map(dots, \(x) mutate(x, ort_block_id = (row_number() - 1) %/% 2, axis = (row_number()-1) %% 2)))

#Pivot to have four columns basically start and end for A fulica and start and end for D laeve
#Basically, we want each alignment segment in its own row
alignments <- alignments %>% mutate(dots = map(dots, \(x) pivot_wider(x, names_from = axis, values_from = c(fulica, laeve))))
#Take out the nested datasets
alignments <- alignments %>% unnest(dots)

#The main thing is to transform coordinates from scaffold specific to the whole assembly. Get info on the length of the alignment as well
alignments <- alignments %>% mutate(laeve_ass_0 = laeve_0 + start_in_ass -1, 
	laeve_ass_1 = laeve_1 + start_in_ass -1, width = sqrt((laeve_0 - laeve_1)^2 + (fulica_0 - fulica_1)^2))
#alignments %>% ggplot() + geom_histogram(aes(width)) + scale_x_log10()
#Based on the histogram, we keep only alignments larger than 150
minalsize <- 150
p <- filter(alignments, width >= minalsize) %>% ggplot() + 
#In the final figure, the color is a function of the alignment segment size
geom_point(aes(x= fulica_0, y = laeve_ass_0, alpha = log(width), col = log10(width)), size = 0.2) + scale_color_gradient2(low = "lightgray", high = "black", mid = "black", space = "Lab", midpoint = log10(2000), limit = c(0, log10(2781))) + theme_bw()

#Get the point in D laeve assembly where the 31 chromosomes end
endoflaeveass <- crom_sizes$end_in_ass[31]

#Draw the blue line separating chromosomes from unplaced scaffolds
p <- p + geom_hline(yintercept = endoflaeveass, col = "blue")

#Easy to get with faSize or similar
afulica_crom_sizes <- read.table("/afulica/Achatina_crom_sizes.txt", col.names = c("scaffold", "length"))

#Same operations for a fulica chromosome sizes
afulica_crom_sizes <- mutate(afulica_crom_sizes, start_in_ass = cumsum(lag(length, default = 1)),
           end_in_ass = length + start_in_ass - 1)
           endofulicass <- afulica_crom_sizes$end_in_ass[31]

p <- p + geom_vline(xintercept = endofulicass, col = "blue")

#Vectors of positions of chromosome boundaries
crombound <- crom_sizes$start_in_ass[2:31]
afulicacrombound <- afulica_crom_sizes$start_in_ass[2:31]

#Paint chromosome boundaries in a red dashed line, vertical if D laeve, horizontal if a fulica
p <- p + geom_hline(yintercept = crombound, col = "red", linetype = "dashed") + geom_vline(xintercept = afulicacrombound, col = "red", linetype = "dashed")



