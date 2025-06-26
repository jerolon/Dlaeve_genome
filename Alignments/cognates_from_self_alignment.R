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

chrom.sizes <- read_tsv("../statistics/derLae_Full_completeStats.txt", col_names = F)
colnames(chrom.sizes) <- c("chromosome", "size", "repeats", "real", "N", "real2")

#Only want chromosomes
chrom.sizes <- chrom.sizes[1:31,]
dir_path <- "toself_inner/"
tsv_files <- paste0(dir_path, chrom.sizes$chromosome, ".tsv")

lastz_segments <- map_dfr(tsv_files, read_tsv)
lastz_segments <- clean_names(lastz_segments)

colnames(lastz_segments)[1] <- "chr_query"
colnames(lastz_segments)[4] <- "chr2_subject"

lastz_segments <- filter(lastz_segments, chr_query != chr2_subject)
lastz_segments <- filter(lastz_segments, chr_query %in% chrom.sizes$chromosome, chr2_subject %in% chrom.sizes$chromosome)


#Nested segments aligned to each query chromosome
nested_gr <- lastz_segments %>%
  group_by(chr_query) %>%
  summarise(
    gr = list(
      GRanges(
        seqnames = chr2_subject,
        ranges = IRanges(start = start2, end = end2),
        strand = strand2
      )
    ),
    .groups = "drop"
  )

#Reduce ranges, we dont want to count twice if the same sequence is in two alignments of the same chromosome
nested_gr <- nested_gr %>% mutate(
    gr_reduced = map(gr, ~ GenomicRanges::reduce(.x, ignore.strand= TRUE))
  )

#seqname_widths just has the length of the alignments to the other chromosomes
nested_gr <- nested_gr %>%
  mutate(
    seqname_widths = map(gr_reduced, ~ {
      as.data.frame(.x) %>%
        group_by(seqnames) %>%
        summarise(total_width = sum(width), .groups = "drop")
    })
  )

#Join each nested dataframe with the info on size and number of repeats  
nested_gr <- nested_gr %>%
  mutate(
    seqname_widths = map(seqname_widths, ~ {
      left_join(.x, chrom.sizes, by = c("seqnames" = "chromosome"))
    })
  )

#Just the columns we need. Unnest into a dataframe
cognates_df <- nested_gr %>%
  select(chr_query, seqname_widths) %>%
  unnest(seqname_widths)

#Normalize by the target chromosome non-repeat or "real" sequence length
cognates_df <- cognates_df %>%
   mutate(alignment_density = total_width / real)

#Get the maximum, to get the cognate of each chromosome
unique_pairs_df <- cognates_df %>% group_by(chr_query) %>% slice_max(order_by=alignment_density) %>% select(chr_query, seqnames) %>% ungroup()

#Remove redundant pairs
unique_pairs_df <- unique_pairs_df %>%
  transmute(
    chr1 = pmin(as.character(chr_query), as.character(seqnames)),
    chr2 = pmax(as.character(chr_query), as.character(seqnames))
  ) %>% 
  distinct(chr1, chr2)

write_tsv(unique_pairs_df, "chromosome_cognates.tsv")

############# Get alignment segments as genomic ranges
homologous_segments <- merge(
  lastz_segments,
  unique_pairs_df,
  by.x = c("chr_query", "chr2_subject"),
  by.y = c("chr1", "chr2")
)


nested_homolog_gr <- homologous_segments %>%
  group_by(chr_query, chr2_subject) %>%
  summarise(
    gr_query = list(
      GRanges(
        seqnames = chr_query,
        ranges = IRanges(start = start1, end = end1),
        strand = strand1
      )
    ),
    gr_subj = list(
      GRanges(
        seqnames = chr2_subject,
        ranges = IRanges(start = start2, end = end2),
        strand = strand2
      )
    ),
    .groups = "drop"
  )
  
nested_homolog_gr <- nested_homolog_gr %>% mutate(
    gr_q_reduced = map(gr_query, ~ GenomicRanges::reduce(.x, ignore.strand= TRUE)),
    gr_s_reduced = map(gr_subj, ~ GenomicRanges::reduce(.x, ignore.strand= TRUE))
  )
  

library(GenomicFeatures)

#From the annotation section
txdb <- loadDb("Annotation/dlaeve_txDb_R_GenomicFeatures.db")

compute_overlap_bp <- function(homolog_range, gene_ranges) {
  # Find overlaps
  olaps <- findOverlaps(homolog_range, gene_ranges, ignore.strand = TRUE)
  
  # Get intersected ranges
  intersected <- pintersect(homolog_range[queryHits(olaps)],
                            gene_ranges[subjectHits(olaps)],
                            ignore.strand = TRUE)
  
  sum(width(GenomicRanges::reduce(intersected)))  # Total genic bp
}


genes_gr <- genes(txdb)

nested_homolog_gr <- nested_homolog_gr %>%
  mutate(
    bp_q_total   = map_int(gr_q_reduced, ~ sum(width(.x))),
    bp_q_genic   = map_int(gr_q_reduced, ~ compute_overlap_bp(.x, genes_gr)),
    bp_q_intergenic = bp_q_total - bp_q_genic,

    bp_s_total   = map_int(gr_s_reduced, ~ sum(width(.x))),
    bp_s_genic   = map_int(gr_s_reduced, ~ compute_overlap_bp(.x, genes_gr)),
    bp_s_intergenic = bp_s_total - bp_s_genic
  )


library(ggforce)

plot_df <- nested_homolog_gr %>%
  mutate(
    total_bp = bp_q_total + bp_s_total,
    genic_bp = bp_q_genic + bp_s_genic,
    intergenic_bp = total_bp - genic_bp
  ) %>%
  dplyr::select(chr_query, chr2_subject, total_bp, genic_bp, intergenic_bp) %>%
  pivot_longer(
    cols = c(genic_bp, intergenic_bp),
    names_to = "type",
    values_to = "bp"
  ) %>%
  group_by(chr_query, chr2_subject) %>%
  mutate(
    fraction = bp / sum(bp),
    cumulative = cumsum(fraction),
    start = lag(cumulative, default = 0),
    end = cumulative,
    mid = (start + end) / 2
  ) %>%
  ungroup()

library(forcats)
plot_df <- plot_df %>%
  mutate(
    chr_query = fct_reorder(chr_query, as.integer(str_extract(chr_query, "\\d+"))),
    chr2_subject = fct_reorder(chr2_subject, as.integer(str_extract(chr2_subject, "\\d+")))
  )

p <- ggplot(plot_df) +
  geom_arc_bar(aes(
    x0 = as.numeric(chr_query),
    y0 = as.numeric(chr2_subject),
    r0 = 0,
    r = sqrt(total_bp) / 2200,  # scale radius for better visibility
    start = 2 * pi * start,
    end = 2 * pi * end,
    fill = type
  ), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("genic_bp" = "forestgreen", "intergenic_bp" = "goldenrod")) +
  coord_fixed() +
  scale_x_continuous(
    breaks = 1:length(levels(plot_df$chr_query)),
    labels = levels(plot_df$chr_query)
  ) +
  scale_y_continuous(
    breaks = 1:length(levels(plot_df$chr2_subject)),
    labels = levels(plot_df$chr2_subject)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Genic vs Intergenic Content per Homologous Chromosomes",
    x = "Query Chromosome",
    y = "Subject Chromosome",
    fill = "Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
