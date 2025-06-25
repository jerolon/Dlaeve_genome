library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(plyranges)
library(ggbio)
library(stringr)

#Table of chromosome sizes
crom_sizes <- read.table("../csvs/cromosome_sizes.csv",
           col.names = c("scaffold", "length"))

crom_sizes <- mutate(crom_sizes, start_in_ass = cumsum(lag(length, default = 1)),
           end_in_ass = length + start_in_ass - 1)


#Chromosomes Assembly. 32 Scaffolds are in mcols. Ranges are assembly based: 0 to length of Genome.
cromsinass = GRanges(seqnames = c("assembly"),
ranges=IRanges(start = crom_sizes$start_in_ass[1:32], end = crom_sizes$end_in_ass[c(1:31,nrow(crom_sizes))]))
 mcols(cromsinass) <- DataFrame(scaffold = c(crom_sizes$scaffold[1:31], "Unplaced"))

#31 scaffolds as seqnames. Ranges are Chromosome based
gr = GRanges(seqnames = c(crom_sizes$scaffold[1:31]) ,
ranges=IRanges(start = c(1), end = crom_sizes$length[1:31]))
 seqlengths(gr) <- width(gr)

#31 scaffolds as seqnames. Ranges are Chromosome based
grfull = GRanges(seqnames = cromsinass$scaffold,
ranges=IRanges(start = c(1), end = width(cromsinass)))
seqlengths(grfull) <- width(grfull)


#High quality complete repeats from HiTE
HiTE_intact <- read_gff("../repeats/HiTE/HiTE_intact.sorted.gff3")

#Must transform all HiC scaffolds higher than 31 to Unplaced. 
meta <- mcols(HiTE_intact)
new_range <- data.frame(old_scaffold = as.vector(seqnames(HiTE_intact)), start = start(HiTE_intact), end = end(HiTE_intact), strand=strand(HiTE_intact))

#Get the shift value which must be 0 for chromosomes but the cumlength for the unplaced.
shiftby <- crom_sizes$start_in_ass[32]
#Must also be 0 for scaffold32 since it is the first of unplaced
crom_sizes <- mutate(crom_sizes, shift_start = if_else(row_number() <= 32, 0, start_in_ass - shiftby))
#Just name all scaffolds after 31 to "Unplaced"
crom_sizes <- mutate(crom_sizes, new_chrom = if_else(row_number() <= 31, scaffold, "Unplaced"))
#Add this info to the ranges we got from HiTE repeats
new_range <- left_join(new_range, crom_sizes, by = join_by(old_scaffold == scaffold))
#and transform so that scaffold 33 starts not at zero, but at the end of scaffold 32, and so on.
new_range <- transmute(new_range, new_chrom, start = start + shift_start, end = end + shift_start, strand)

#make a new ranges object for the HiTe repeats with just 31 chromosomes and an "unplaced" supra-scaffold
trHiTE_intact <- GRanges(seqname=new_range$new_chrom, 
		     	 ranges = IRanges(start=new_range$start,
                                           end = new_range$end),
 			 strand= new_range$strand)
#Add the metadata with family and other info
mcols(trHiTE_intact) <- meta

#Homogenize seq levels
seqlevels(trHiTE_intact) <- seqlevels(grfull)
seqinfo(trHiTE_intact) <- seqinfo(grfull)

#How many of each class
table(trHiTE_intact$classification)

#Make 1Mb windows along the chromosome
windows <- tileGenome(seqlengths = seqlengths(grfull), tilewidth = 1e6, cut.last.tile.in.chrom=T)
#Keep only windows that actually overlap some transposon
non_zero <- subsetByOverlaps(windows, trHiTE_intact)
#find the overlap between windows and transposons
overlaps <- findOverlaps(non_zero, trHiTE_intact)
overlaps_by_window <- split(subjectHits(overlaps), queryHits(overlaps))

#Columns to stores the coverage of different elements
non_zero$helitron <- 0
non_zero$dna <- 0
non_zero$line <- 0
non_zero$LTR <- 0
non_zero$unknown <- 0
non_zero$sine <- 0


#This takes a loooooong time
#For each window
for(i in seq_along(non_zero)){
  #get the hits to transposons
hits <- overlaps_by_window[[i]]
  #Get all the hits to a particular class. "Reduce" in order to not count elements that overlap twice.
  #We could get away with dividing by 1,000,0000 but some windows at the end of chromosomes are not necessarily 1Mb
  #Get the density of a given element by dividing its extent by the window width
non_zero$dna[i] <- sum(width(GenomicRanges::reduce(filter(trHiTE_intact[hits], str_detect(classification, "DNA"))))) / width(non_zero[i])
non_zero$line[i] <- sum(width(GenomicRanges::reduce(filter(trHiTE_intact[hits], str_detect(classification, "LINE"))))) / width(non_zero[i])
non_zero$LTR[i] <- sum(width(GenomicRanges::reduce(filter(trHiTE_intact[hits], str_detect(classification, "LTR"))))) / width(non_zero[i])
non_zero$unknown[i] <- sum(width(GenomicRanges::reduce(filter(trHiTE_intact[hits], classification == "Unknown")))) / width(non_zero[i])
non_zero$sine[i] <- sum(width(GenomicRanges::reduce(filter(trHiTE_intact[hits], str_detect(classification, "SINE"))))) / width(non_zero[i])
non_zero$helitron[i] <- sum(width(GenomicRanges::reduce(filter(trHiTE_intact[hits], classification == "RC/Helitron")))) / width(non_zero[i])
}


#Write this files automatically or individually to use later.
library(readr)
repeat_classes <- colnames(mcols(non_zero))
lapply(repeat_classes, function(clase){
write_tsv(x = (as.data.frame(non_zero) %>% select(seqnames, start, end, clase)),
file = paste0(clase, "_HiTE_mb_density.txt"),
col_names = FALSE)
})

#Individual 
tosave <- as.data.frame(non_zero) %>% select(seqnames, start, end, line)
write_tsv(tosave, col_names = FALSE, file = "line_HiTE_mb_density.txt")





#########################
# Process tandem repeats
#########################

#Read the sorted bed file 
trf_raw <- read_table("derLae.TRFsorted.bed", col_names = 
c("scaffold", "start_bed", "end_bed", "program", "period_size", "n_copies", "consensus_size",
"p_matches_adjacent", "p_indels_adjacent", "aln_score", "p_adenine", "p_cytosine", "p_guanine", "p_timine", "entropy", "consensus"))

#Just a trick to remove overlaps, group by start and take the one with the largest end.
#Because the file is huge
trf_filter <- group_by(trf_raw, scaffold, start_bed) %>% arrange(desc(end_bed)) %>% slice_head(n=1)
#Transform to gRnages in order to reduce
#+1 because bed is zero based
trf <- GRanges(seqnames = trf_filter$scaffold, ranges=IRanges(start = trf_filter$start_bed+1, end = trf_filter$end_bed))

#read the sequences
library(Biostrings)
genoma <- readDNAStringSet("../derLae1_hic.FINAL.fasta")

#Only non-redundant intervals
trf_reduced <- reduce(trf)

#Get the sequences of the tandem repeat intervals
trf_seq <- getSeq(genoma, trf_reduced)

#Telomeres
findTelomereRep <- function(secuencia){
telo <- Biostrings::matchPattern("TTAGGG", secuencia)
return(length(telo))
}

#Get the number of times the motif was found
ntelo <- sapply(trf_seq, findTelomereRep)
#Get the GC content of the repeat
gc_trf <- letterFrequency(trf_seq, "GC", as.prob = TRUE)

#So that 0 is 50%, negative AT-rich
trf_reduced$gc <- gc_trf[,1] - 0.5 
trf_reduced$n_telomere_repeats <- ntelo

trf_df <- as.data.frame(trf_reduced) %>% select(-strand)

#reuse_crom_sizes to transform to Unplaced
trf_df <- left_join(trf_df, crom_sizes, by = join_by(seqnames == scaffold))

#Transform the Unplaced coordinates
trf_df <- trf_df %>% mutate(start = start + shift_start, end = end + shift_start)

#Keep only repeats greater than 100kb. Unless the original Unplaced scaffold length was less than that

#The files with four columns can be used directly by circos.
tosave <- mutate(trf_df, tokeep = if_else(length > 100000 & width < 100000, "NO", "YES")) %>% filter(tokeep == "YES") %>% select(new_chrom, start, end, gc)
write_tsv(tosave, col_names = FALSE, file = "tandemsGC.txt")

#Save telomeres
tosave <- filter(trf_df, n_telomere_repeats > 0) %>% select(new_chrom, start, end, n_telomere_repeats)
write_tsv(tosave, col_names = FALSE, file = "telomeres.txt")


telomeres <- GRanges(seqnames = trf_df$new_chrom, ranges = IRanges(start = trf_df$start, end = trf_df$end))
seqlevels(telomeres) <- seqlevels(grfull)
seqlengths(telomeres) <- seqlengths(grfull)
mcols(telomeres) <- trf_df %>% select(n_telomere_repeats, width)

#pirnas
#Assume pirnas has table S6 where coordinates are in string form like "chr12:54031931-54185160" so we must separate them
pirnas <- separate_wider_delim(pirnas, coordinates, delim = ":", names = c("seqnames", "range")) %>% separate_wider_delim(range, "-", names = c("start", "end"))

pirnas <- na.omit(pirnas)
pirnas$exists <- 1

tosave <- pirnas %>% select(seqnames, start, end, exists)
write_tsv(tosave, col_names = FALSE, file = "pirnas_clusters.txt")

############################
#Plot linear with R
############################
library(readr)
library(patchwork)

#If re-reading the files
dens_df <- tibble(family=colnames(mcols(non_zero)))
#Karyotype
library(purrr)
#Read the density files again
dens_df <- mutate(dens_df, data = map(family, function(x){read_tsv(paste0("//repeats/", x, "_HiTE_mb_density.txt"), col_names = c("seqnames", "start", "end", "density"))}))

dens_df <- unnest(dens_df, data)

chroms <- seqlevels(grfull)

tandemsGC <- read_tsv("/repeats/tandemsGC.txt", col_names = c("scaffold", "start", "end", "gc"))

ttaggg <- read_tsv("/telomeres.txt", col_names = c("scaffold", "start", "end", "n_repeats"))

pirnas <- read_tsv("/repeats/pirnas_clusters.txt", col_names = c("scaffold", "start", "end", "n_clusters"))

#Transform these densities to negative so that they are plotted as a streamchart
dens2 <- dens_df %>% mutate(
    density = ifelse(family %in% c("dna","helitron","sine"),
                    density, -density)
  )


chroms <- seqlevels(grfull)

library(scales) #To get the x scale in Mb units

#Function to plot a single scaffold
make_scaffold_panel <- function(scaf) {
  #p1 is the telomere track
  chrom_size <- seqlengths(grfull)[scaf]
  p1 <- ttaggg %>%
    filter(scaffold == scaf) %>%
    mutate(ldensidad = log2(n_repeats + 1)) %>%
    ggplot(aes(xmin = start, xmax = end, ymin = -1, ymax = pmin(ldensidad,log2(50)))) +
    geom_rect(fill = "red") +
    labs(title=scaf, y=NULL, x=NULL) +
    theme_minimal() + xlim(0,chrom_size) +
    scale_x_continuous(labels = comma_format(), expand   = c(0, 0)) +
    theme(
    legend.position = "none", axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    plot.margin       = unit(c(0, 0, -50, 0), "pt"), plot.title      = element_blank())
  
  #Tandem GC track
  p2 <-  ggplot(tandemsGC %>% filter(scaffold == scaf)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=0.1, fill = gc>-0.05)) +
    scale_fill_manual(values=c(`TRUE`="magenta", `FALSE`="blue")) +
    geom_point(
      data = pirnas %>% filter(scaffold == scaf),
      aes(x = (start + end) / 2, y = 0),
      shape = 17,      # filled triangle
      color = "red",
      size  = 1
    ) +
    labs(y=NULL, x=NULL) +
    theme_minimal() + xlim(0,chrom_size) +
    scale_x_continuous(labels = comma_format(), expand   = c(0, 0)) +
    theme(
    legend.position = "none",axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    plot.margin = unit(c(-50, 0, -50, 0), "pt"), plot.title = element_blank())
  
  #Repeat density track
  p3 <- dens2 %>%
    filter(seqnames == scaf) %>%
    ggplot(aes(x=start, y=density, fill=family)) +
    geom_area(position = "stack") + 
    scale_x_continuous(labels = scales::label_number(scale=1e-6, suffix=" Mb", accuracy=1),
    expand   = c(0, 0)) +
    xlim(0,chrom_size) +
    theme_minimal() +
    theme(legend.position = "none",
     plot.margin     = unit(c(-100, 0, 0, 0), "pt"), plot.title      = element_blank())+
    labs(y=NULL, x="")
    
  
  # Stack them, giving relative heights
  p1 / p2 / p3 +
  plot_layout(heights = c(0.5, 0.5, 2), ncol = 1) &
  theme(plot.background = element_blank())
}

#Manually and tortuosly give the scaffold to chromosome translation
chroms <- c("HiC_scaffold_1","HiC_scaffold_5","HiC_scaffold_4","HiC_scaffold_7","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_8","HiC_scaffold_10","HiC_scaffold_6","HiC_scaffold_13","HiC_scaffold_12","HiC_scaffold_14","HiC_scaffold_16","HiC_scaffold_9","HiC_scaffold_11","HiC_scaffold_17", 
            rev(c("Unplaced", "HiC_scaffold_15","HiC_scaffold_18","HiC_scaffold_21","HiC_scaffold_20","HiC_scaffold_19","HiC_scaffold_23","HiC_scaffold_22","HiC_scaffold_26","HiC_scaffold_25","HiC_scaffold_24","HiC_scaffold_29","HiC_scaffold_27","HiC_scaffold_28","HiC_scaffold_31","HiC_scaffold_30")))

#Apply the function to get the plots for all chromosomes
panels     <- setNames(lapply(chroms, make_scaffold_panel), chroms)

#This is arbitrary, you can plot any chromosome you like so that they fit in a page
main_figure <- wrap_plots(panels[paste0("HiC_scaffold_",c(1,7,13,19,22,31))], ncol = 2, byrow = F)

pdf("repeats_scaffolds1_7_13_19_22_31.pdf", width = 8, height = 5.33)
main_figure
dev.off()


sup_fig1 <- wrap_plots(panels[paste0("HiC_scaffold_",c(5,4,2,3,8,10,6,12,14,16))], ncol = 2, byrow = F)

pdf("repeats_two_col_10scaffolds.pdf", width = 8, height = 5.33)
sup_fig1
dev.off()


sup_fig2 <- wrap_plots(panels[paste0("HiC_scaffold_",c(9,11,17,15,18,21,20,23,26,25,24,29,27,28,30))], ncol = 3, byrow = F)

pdf("repeats_three_col_15scaffolds.pdf", width = 8, height = 5.33)
sup_fig2
dev.off()

#Unplaced, same, but with legends
scaf <- "Unplaced"
  chrom_size <- seqlengths(grfull)[scaf]
  p1 <- ttaggg %>%
    filter(scaffold == scaf) %>%
    mutate(ldensidad = log2(n_repeats + 1)) %>%
    ggplot(aes(xmin = start, xmax = end, ymin = -1, ymax = pmin(ldensidad,log2(50)))) +
    geom_rect(fill = "red") +
    labs(title=scaf, y="Telomeric\nlog10(n+1)", x=NULL) +
    theme_minimal() + xlim(0,chrom_size) +  
    scale_x_continuous(labels = comma_format(), expand = c(0,0)) +
    theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    plot.margin       = unit(c(0, 0, -50, 0), "pt"), plot.title      = element_blank())
  
  # 2. Tandem GC track
  p2 <-  ggplot(tandemsGC %>% filter(scaffold == scaf)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=0.1, fill = gc>-0.05)) +
    scale_fill_manual(values=c(`TRUE`="magenta", `FALSE`="blue"),
    labels = c(`FALSE` = "AT-rich", `TRUE` = "GC-rich"),
    name   = "") +
    geom_point(
      data = pirnas %>% filter(scaffold == scaf),
      aes(x = (start + end) / 2, y = 0),
      shape = 17,      # filled triangle
      color = "red",
      size  = 1
    ) +
    labs(y=NULL, x=NULL) +
    theme_minimal() + xlim(0,chrom_size) +  
    scale_x_continuous(labels = comma_format(), expand = c(0,0)) +
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    plot.margin = unit(c(-50, 0, -50, 0), "pt"), plot.title = element_blank())
  
  # 3. Repeat density track (segments)
  p3 <- dens2 %>%
    filter(seqnames == scaf) %>%
    ggplot(aes(x=start, y=density, fill=family)) +
    geom_area(position = "stack") + 
    scale_x_continuous(
    labels   = scales::label_number(scale=1e-6, suffix=" Mb", accuracy=1),
    expand   = c(0, 0)) +
    xlim(0,chrom_size) +
    #Limit of the big bacterial chromosomes
    geom_vline(xintercept = 46461513, linetype = 2) + 
    theme_minimal() +
    theme(plot.margin     = unit(c(-100, 0, 0, 0), "pt"), 
    plot.title      = element_blank())+
    labs(y="Repeat\nsigned density", x="Genomic Coordinates (bp)")
    
  
  # Stack them, giving relative heights (e.g. 1:1:2)
  (unplaced_plot <- p1 / p2 / p3 +
  plot_layout(heights = c(0.5, 0.5, 2), ncol = 1) &
  theme(plot.background = element_blank()))

pdf("repeats_unplaced.pdf", width = 8, height = 2)
unplaced_plot
dev.off()
