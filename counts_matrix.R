library(tidyverse)

counts <- read.delim("I:/psivakumar/MSA_RNASeq/analysis/MSA_RNASeq_counts.txt", header = TRUE, skip = 1)

counts_matrix <- counts %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-c(1:5))
  
sample_names <- names(counts_matrix) %>%
  as.list() %>%
  sapply(strsplit, split = "_") %>%
  sapply(`[[`, 2) %>%
  sapply(strsplit, split = "bams.") %>%
  sapply(`[[`, 2)

names(counts_matrix) <- paste("ID_", gsub(pattern = ".", replacement = "_", sample_names, fixed = T), sep = "")

counts.matrix <- matrix(as.numeric(as.character(unlist(counts_matrix))), nrow = nrow(counts_matrix))
colnames(counts.matrix) <- names(counts_matrix)
row.names(counts.matrix) <- row.names.data.frame(counts_matrix)

#save(counts.matrix, file = "I:/psivakumar/MSA_RNASeq/analysis/MSA_RNASeq_counts_matrix.RData")

#################################################
#checking if marking and ignoring duplicates == to removing them from bams entirely

counts <- read.delim("I:/psivakumar/MSA_RNASeq/analysis/analysis_rmDups/MSA_RNASeq_rmDups_counts.txt", header = TRUE, skip = 1)

#markDups matrix
markDups_counts_matrix <- counts %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(ends_with("Processed.out.bam"))

markDups_sample_names <- names(markDups_counts_matrix) %>%
  as.list() %>%
  sapply(strsplit, split = "_") %>%
  sapply(`[[`, 3) %>%
  sapply(strsplit, split = "Dups.") %>%
  sapply(`[[`, 2)

names(markDups_counts_matrix) <- paste("ID_", gsub(pattern = ".", replacement = "_", markDups_sample_names, fixed = T), sep = "") 

markDups_counts.matrix <- matrix(as.numeric(as.character(unlist(markDups_counts_matrix))), nrow = nrow(markDups_counts_matrix))
colnames(markDups_counts.matrix) <- names(markDups_counts_matrix)
row.names(markDups_counts.matrix) <- row.names.data.frame(markDups_counts_matrix)

#rmDups matrix
rmDups_counts_matrix <- counts %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(ends_with("removedDups.bam"))

rmDups_sample_names <- names(rmDups_counts_matrix) %>%
  as.list() %>%
  sapply(strsplit, split = "_") %>%
  sapply(`[[`, 3) %>%
  sapply(strsplit, split = "Dups.") %>%
  sapply(`[[`, 2)

names(rmDups_counts_matrix) <- paste("ID_", gsub(pattern = ".", replacement = "_", rmDups_sample_names, fixed = T), sep = "") 

rmDups_counts.matrix <- matrix(as.numeric(as.character(unlist(rmDups_counts_matrix))), nrow = nrow(rmDups_counts_matrix))
colnames(rmDups_counts.matrix) <- names(rmDups_counts_matrix)
row.names(rmDups_counts.matrix) <- row.names.data.frame(rmDups_counts_matrix)

identical(rmDups_counts.matrix, markDups_counts.matrix)
#counts.matrix <- rmDups_counts.matrix

save(counts.matrix, file = "I:/psivakumar/MSA_RNASeq/analysis/analysis_rmDups/MSA_RNASeq_rmDups_counts_matrix.RData")
