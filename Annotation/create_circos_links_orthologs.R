library(GenomicFeatures)
library(dplyr)
library(readr)
library(tidyr)

croms <- names(seqlengths(dlaevm))
crom_sizes <- read.table("../csvs/cromosome_sizes", col.names = c("chrom", "length"))
croms <- filter(crom_sizes, chrom %in% croms)
dlaevm <- GenomicFeatures::makeTxDbFromGFF("../Dlaeve_annotated.gff3", organism = "Deroceras laeve", chrominfo= croms)

dl_genes <- genes(dlaevm)
dl_exons_by_tx <- exonsBy(dlaevm, by = "tx")

#Mono exon multi exon
monoexon <- width(dl_exons_by_tx) %>% 
	lapply(length) %>% unlist() %>%
	 `==`(1) %>% sum

multiexon <- width(dl_exons_by_tx) %>% 
	lapply(length) %>% unlist() %>%
	 `>=`(2) %>% sum

#Median number of exons
width(dl_exons_by_tx) %>% lapply(length) %>% unlist %>% median
width(dl_exons_by_tx) %>% lapply(length) %>% unlist %>% mean
width(dl_exons_by_tx) %>% lapply(length) %>% unlist %>% max

#Length of exons
width(dl_exons_by_tx) %>% unlist %>% sum
width(dl_exons_by_tx) %>% unlist %>% mean

#Introns
dl_intrns_by_tx <- intronsByTranscript(dlaevm)

nintron <- sapply(dl_intrns_by_tx, length)
nintron[which(nintron > 0)] %>% median
nintron[which(nintron > 0)] %>% median

width(dl_intrns_by_tx) %>% unlist %>% mean

saveDb(dlaevm, "dlaeve_txDb_R_GenomicFeatures.db")

dltxdb <- loadDb("dlaeve_txDb_R_GenomicFeatures.db")
jaccard <- read.table("full_jaccard_orthologs.csv", header = T)

transcritos <- transcripts(dltxdb)
crom_sizes <- read.table("../csvs/cromosome_sizes", col.names = c("scaffold", "length"))

crom_sizes <- mutate(crom_sizes, start_in_ass = cumsum(lag(length, default = 1)),
	end_in_ass = length + start_in_ass)

shiftby <- crom_sizes$start_in_ass[32]
crom_sizes <- mutate(crom_sizes, shift_start = if_else(row_number() <= 32, 0, start_in_ass - shiftby))
crom_sizes <- mutate(crom_sizes, new_chrom = if_else(row_number() <= 31, scaffold, "Unplaced"))

trans_df <- data.frame(genes =  transcritos$tx_name, scaffold = seqnames(transcritos),
  start = start(transcritos), end = end(transcritos))

trans_df <- left_join(trans_df, crom_sizes)

trans_df <- transmute(trans_df, genes, new_chrom, start = start + shift_start, end = end + shift_start)

df <- jaccard %>% select(-jaccard) %>% filter(gene1 != "gene1")


df <- left_join(df, trans_df, by = join_by(gene1 == genes))
df <- left_join(df, trans_df, by = join_by(gene2 == genes))
df <- df %>% select(-c(gene1, gene2))

write_tsv(df, "links.txt", col_names = FALSE)
