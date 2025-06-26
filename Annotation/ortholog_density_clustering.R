library(GenomicFeatures)
library(dplyr)
library(readr)
library(igraph)
library(plyranges)
library(stringr)
#Load the database of genomic features to get localization of genes
dltxdb <- loadDb("../dlaeve_txDb_R_GenomicFeatures.db")

#Load jaccard to get orthology relationships
jaccard <- read.table("/mnt/Timina/alfredvar/jmiranda/50-Genoma_Referencia_Dlaeve/59-scripts/full_jaccard_orthologs.csv", header = F)

transcritos <- transcripts(dltxdb)
colnames(jaccard) <- c("gene1", "gene2", "jaccard_similarity")

trans_df <- data.frame(genes =  transcritos$tx_name, scaffold = seqnames(transcritos),
  start = start(transcritos), end = end(transcritos))
  
df <- left_join(jaccard, trans_df, by = join_by(gene1 == genes))
df <- left_join(df, trans_df, by = join_by(gene2 == genes))

#Filter cis-chromosomal orthologs
df <- df %>% filter(scaffold.x == scaffold.y)

df <- rowwise(df) %>% mutate(span = max(start.x, start.y, end.x, end.y) - min(start.x, start.y, end.x, end.y), 
	widths = end.x - start.x + end.y - start.y, genomic_distance = span - widths) 

df <- mutate(df, genomic_distance = if_else(genomic_distance <= 1, 2, genomic_distance))
df$jaccard_similarity <- as.numeric(df$jaccard_similarity)

gene_connections <- select(df, gene1, gene2)
weights <- 1/log2(df$genomic_distance) * (df$jaccard_similarity)
graph <- graph_from_data_frame(gene_connections, directed = FALSE)
E(graph)$weight <- weights
clusters <- cluster_louvain(graph, weights = E(graph)$weight)
membership <- membership(clusters)
print(membership)
groups <- split(clusters$names, clusters$membership)

#Is this necessary or just split by membership?
group_ranges <- lapply(groups, \(x){filter(transcritos, tx_name %in% x)})

#Strand: direct or inverted repeats?
strands <- lapply(group_ranges, \(x){as.vector(strand(x))})

onlypairs <- strands[sapply(strands, function(x) length(x) == 2)]
logicalstrand <- lapply(onlypairs, \(x){x == "+"})
#True is when the strands are different
sapply(logicalstrand, \(x){xor(x[1],x[2])})

plot(graph, layout = layout_with_mds, vertex.size = 2, vertex.label = NA)
plot(graph, layout = layout_with_mds, vertex.size = 2, vertex.label.cex = 0.2)
#Index of communities that are segmental duplications
seg_dups <- which(sizes(clusters) == 2)
#Actual clusters (bigger than 2)
gene_clusters <- which(sizes(clusters) > 2)
gene_clusters_list <- communities(clusters)[gene_clusters]
#gene ids that belong to segmental duplications
segmental_duplications_list <- communities(clusters)[seg_dups]
genes_seg_dups <- unlist(segmental_duplications_list)
#Graph with communities bigger than 2 only
non_segmental_dup <- delete_vertices(graph, genes_seg_dups)

#function to get a range object from a list of ids
get_range_from_txid <- function(listoftx){
transcritos %>% filter(tx_name %in% listoftx)
}


#segmental_duplications_ranges <- lapply(segmental_duplications_list, get_range_from_txid)
#Without ignore.strand you can derive head-tail tail-tail or head head because head to tail have a single strand,
#tail tail are + then -  and head-head are - then +
#tosave <- lapply(segmental_duplications_ranges, range, ignore.strand = T) %>% GRangesList %>% unlist()
#Segmental duplications are already marked in the circos plot


gene_clusters_ranges <- lapply(gene_clusters_list, get_range_from_txid)
tosave <- lapply(gene_clusters_ranges, range, ignore.strand = T) %>% GRangesList %>% unlist()

#Get gene functions
gene_functions <- read_tsv("../51-Functions_Domains_Orthologs/functions.tsv", skip = 6)
gene_functions <- na.omit(gene_functions)
library(janitor)
gene_functions <- clean_names(gene_functions)

get_gene_function <- function(genelist){
filter(gene_functions, number_query %in% genelist) %>% select(cog_category, description, preferred_name, pfa_ms) %>% unique()
} 

(funciones_genes <- lapply(gene_clusters_list, get_gene_function))

tibble_2_1row <- function(input_tibble){
input_tibble %>% 
summarise(across(everything(), ~ str_c(unique(.), collapse = ", ")))
}

functions_df <- lapply(funciones_genes, tibble_2_1row) %>% bind_rows()

tosave <- cbind(as.data.frame(tosave), functions_df)
tosave <- tibble(tosave)
tosave$nmembers <- lapply(gene_clusters_list, length)

tosavegood <- tosave %>% filter(width < 5e+06, nmembers > 4) %>% arrange(desc(width), seqnames) %>% select(seqnames, start, enscription, width, cog_category, preferred_name)

write_tsv(tosavegood, "gene_clusters.txt", col_names = F) 

#Homeobox genes

hoxgenes <- filter(gene_functions, str_starts(preferred_name, "HOX")) %>% select(number_query, preferred_name)
hoxranges <- get_range_from_txid(hoxgenes$number_query) %>% as.data.frame()

hox_df inner_join(hoxgenes, hoxranges, by = join_by(number_query == tx_name))

write_tsv(select(hox_df, seqnames, start, end, preferred_name), "hox_clusters.txt", col_names = F) 
