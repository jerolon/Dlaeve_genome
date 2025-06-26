library(dplyr)
library(readr)
library(purrr)
library(Matrix)
library("doFuture")
library(doParallel)
library(foreach)
library(stringr)
message("start")

del_ortho <- readr::read_tsv("orthologs_non_redundant.tsv", skip = 10, col_names = c("query", "orth_type", "species", "orthologs"))
del_ortho2 <- del_ortho |> select(query, orthologs)
del_ortho2 <- del_ortho2 |> mutate(ortos= map(orthologs, \(x){str_split(x, ",") %>% unlist}))
del_ortho2 <- del_ortho2 |> group_by(query) |>  summarise(uniqueOrthos = list(unique(unlist(ortos))))

#All the unique hits from EggNOGG, to give them a unique ID
all_orthologs <- unique(unlist(del_ortho2$uniqueOrthos))
#Here, we match a given ortholog to its unique number from the all_orthologs vector
ortholog_index <- match(unlist(del_ortho2$uniqueOrthos), all_orthologs)
#Number of dlaeve genes with at least one hit to EGGnog
ngenes <- length(del_ortho2$query)
# Construct sparse matrix
#Here, each d laeve gene is given an ID from 1 to ngenes, we make as many rows as unique orthologs that gene has
rows <- rep(1:ngenes, sapply(del_ortho2$uniqueOrthos, length)) # Gene indices
#The columns are the ortholog indives
cols <- ortholog_index                                   # Ortholog indices
values <- rep(1, length(rows))                          # Presence indicator
sparse_matrix <- sparseMatrix(i = rows, j = cols, x = values,
                               dims = c(ngenes, length(all_orthologs)),
                               dimnames = list(del_ortho2$query, all_orthologs))

head(sparse_matrix)

# Step 2: Jaccard Similarity Function (Efficient for Sparse Matrices)
jaccard_sparse <- function(i, j, matrix) {
  intersect_count <- sum(matrix[i, ] & matrix[j, ])
  union_count <- sum(matrix[i, ] | matrix[j, ])
  return(intersect_count / union_count)
}


message("Calc combinations")
gene_combinations <- combn(nrow(sparse_matrix), 2)
num_combinations <- ncol(gene_combinations)
message(num_combinations)

# Step 3: Set up Parallel Processing
#registerDoFuture()
options(future.globals.maxSize = 943718400)
num_cores <- 43
#cl <- makeCluster(num_cores)
plan(multicore, workers = num_cores )
# Step 4: Parallel Computation
message("export to each worker")
#plan(future.batchtools::batchtools_sge)
#clusterExport(cl, varlist = c("jaccard_sparse", "sparse_matrix", "gene_combinations"), envir = environment())
message("start loop")
for(x in 1:99){
count = x * 1000000
  start <- (x-1) * 1000000 + 1
results <- foreach(k = start:count, .combine = 'rbind', .options.future = list(chunk.size = 200L)) %dofuture% {
  i <- gene_combinations[1, k]
  j <- gene_combinations[2, k]
  similarity <- jaccard_sparse(i, j, sparse_matrix)
  return(c(Gene1 = i,
             Gene2 = j,
             Jaccard = similarity))
  }

message("done")
filtered_results <- results[results[, "Jaccard"] > 0.5, ]
# Convert to Data Frame
nombres <- rownames(sparse_matrix)
results_df <- data.frame(gene1 = nombres[filtered_results[,1]],
gene2 = nombres[filtered_results[,2]],
jaccard = filtered_results[,3])

# View Results
n <- str_pad(x, 2, pad="0")
write_tsv(results_df, paste0(n, "_jaccard_orthologs.csv"))
}

#The last cycle
results <- foreach(k = 100000000:num_combinations, .combine = 'rbind', .options.future = list(chunk.size = 200L)) %dofuture% {
  i <- gene_combinations[1, k]
  j <- gene_combinations[2, k]
  similarity <- jaccard_sparse(i, j, sparse_matrix)
  return(c(Gene1 = i,
             Gene2 = j,
             Jaccard = similarity))
  }
write_tsv(results_df, "last_jaccard_orthologs.csv")
message("done")
