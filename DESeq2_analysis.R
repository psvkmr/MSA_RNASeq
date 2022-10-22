library(tidyverse)
library(biomaRt)
library(DESeq2)
library(gProfileR)
library(gplots)

load("I:/psivakumar/MSA_RNASeq/analysis/MSA_RNASeq_counts_matrix.RData")
#load("I:/psivakumar/MSA_RNASeq/analysis/analysis_rmDups/MSA_RNASeq_rmDups_counts_matrix.RData")
sample_info <- read.delim("I:/psivakumar/MSA_RNASeq/analysis/MSA_rna_seq_master.txt")

#load gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biomart_ids <- getBM(mart = mart, attributes = c("ensembl_gene_id", "external_gene_name"))

#modify sample info df
sample_info[c(7, 21), 2] <- "done"
sample_info <- filter(sample_info, RNA.seq == "done")
sample_info <- arrange(sample_info, Patient.ID.number)

#convert sample_info ids to df ids
partial_ids <-    
  sapply(colnames(counts.matrix), strsplit, split = "_") %>% 
  sapply('[[', 2)

sample_info$df_id <- c()

for (i in 1:nrow(sample_info)) {
  sample_info$df_id[[i]] <- agrep(partial_ids[[i]], colnames(counts.matrix), value = TRUE, max = list(sub = 0, ins = 0))
}

mutant_ids <- filter(sample_info, Group == "MSA")$df_id

#Removed ID-01129-53593635 due to low counts
MSA_counts <- counts.matrix[, colnames(counts.matrix) != "ID_01129_53593635"]

#remove 0 counts
MSA_counts <- MSA_counts[apply(MSA_counts, 1, function(x) all(x != 0)), ]

#pseudocounts for analysis
pseudocounts <- log2(MSA_counts + 1)
pseudocounts_df <- pseudocounts %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_gene_id")

#hist of counts distribution for one sample
msa_counts_hist <- 
  ggplot(pseudocounts_df, aes(x = `ID_01020-53530554`)) + 
  geom_histogram(binwidth = 0.6)

msa_melt_counts <- pseudocounts_df %>%
  gather(key = sample, value = counts, -ensembl_gene_id) %>%
  left_join(dplyr::select(sample_info, c(genotype = "Group", "df_id")), by = c("sample" = "df_id"))

msa_counts_boxplot <-
  ggplot(msa_melt_counts, aes(sample, y = counts, fill = genotype)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#data frame of variants
msa_col_df <- dplyr::select(sample_info, sample = df_id, genotype = Group) %>%
  filter(sample %in% colnames(MSA_counts))

#design setting
design_setting <- formula(~ genotype)

#model matrix
modelled_matrix <- model.matrix(design_setting, msa_col_df)

#create dds
msa_dds <- DESeqDataSetFromMatrix(countData = MSA_counts, 
                                  colData = msa_col_df, 
                                  design = design_setting)
#remove sex genes
sex_genes <- getBM(mart = mart, attributes = c("ensembl_gene_id", "chromosome_name")) %>%
  filter(chromosome_name == "X" | chromosome_name == "Y")
msa_dds <- msa_dds[!rownames(msa_dds) %in% sex_genes$ensembl_gene_id, ]

#normalise data with variance stabilising transformation for PCA
msa_dds_vst <- varianceStabilizingTransformation(msa_dds)

msa_dists_vst <- msa_dds_vst %>%
  assay() %>%
  t() %>%
  dist()

#Euclidean distances heatmap
create_heatmap <- function(x) {
  plot_heatmap <- function() heatmap(as.matrix(x))
}

msa_dds_heatmap <- create_heatmap(msa_dists_vst)

msa_pca <-
  BiocGenerics::plotPCA(msa_dds_vst, intgroup = "genotype") +
  theme_classic()

#normalised_counts
size_factors <- estimateSizeFactors(msa_dds)
size_factors_df <- sizeFactors(size_factors)
normalised_counts <- DESeq2::counts(size_factors, normalized = TRUE)

#plot normalised counts (assume biomart ids loaded)
normalised_counts_scatter_plotter <- function(gene_name) {
  ifelse(grepl("ENSG", toupper(gene_name)), ensemblID <- gene_name, ensemblID <- filter(biomart_ids, external_gene_name == toupper(gene_name))$ensembl_gene_id)
  df <- rownames_to_column(as.data.frame(normalised_counts))
  df <- df[grep(pattern = ensemblID, x = df$rowname), ] %>%
    gather(key = sample, value = count, -rowname) %>%
    mutate(genotype = ifelse(sample %in% mutant_ids, "MSA", "CTRL"))
  ggplot(df, aes(sample, count, colour = genotype)) +
    geom_point() +
    scale_color_manual(values = c("dark grey", "red")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle(paste0(ensemblID, " _ ", filter(biomart_ids, ensembl_gene_id == ensemblID)$external_gene_name))
}

#deseq results
msa_deseq_est <- DESeq(msa_dds)
msa_deseq_res <- results(msa_deseq_est)
msa_deseq_summary <- summary(msa_deseq_res)

msa_deseq <- as.data.frame(msa_deseq_res) %>%
  rownames_to_column(var = "EnsemblID") %>%
  left_join(biomart_ids, by = c("EnsemblID" = "ensembl_gene_id")) %>%
  arrange(padj) 

msa_volcano <- msa_deseq %>%
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(colour = ifelse(padj < 0.1, 
                                 ifelse(log2FoldChange > 0, "Up", "Down"), 
                                 "None"))) +
  scale_colour_manual(name = "Differential expression",
                      limits = c("Up", "Down", "None"),
                      values = c("#d40000", "#2c5aa0", "#cccccc")) +
  xlim(-2, 2) +
  ylim(0, 5) +
  theme_classic()
  

msa_deseq_with_counts <-
  left_join(msa_deseq, rownames_to_column(as.data.frame(normalised_counts), var = "EnsemblID"), 
            by = "EnsemblID")

msa_pvalue_hist <-
  ggplot(na.omit(msa_deseq), aes(x = pvalue)) +
  geom_histogram(bins = 100)

#GO terms
go_terms <-
  gprofiler(query = filter(na.omit(msa_deseq), padj < 0.1)$EnsemblID, 
            organism = "hsapiens", 
            significant = TRUE,
            min_set_size = 5,
            min_isect_size = 2,
            correction_method = "fdr", 
            custom_bg = na.omit(msa_deseq)$EnsemblID) %>%
  arrange(p.value)

msa_go <- 
  ggplot(go_terms, aes(term.name, y = -log10(p.value))) +
  geom_col() +
  coord_flip()

#DGE heatmap + clustering
sigf_genes <- filter(msa_deseq, padj < 0.1)$EnsemblID

#Incomplete
create_deseq_heatmap <- function(x) {
  plot_deseq_heatmap <- function() {
    y <- x[sigf_genes, ]
    y <- y - rowMeans(y)
    heatmap.2(y, 
              scale = "none", 
              na.rm = TRUE, 
              trace = "none", 
              labRow = NA, 
              srtCol = 45, 
              offsetCol = 0)
  }
}

#write.table(msa_deseq, "I:/psivakumar/MSA_RNASeq/analysis/deseq_MSA_blood_RNASeq.tab", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(msa_deseq, "I:/psivakumar/MSA_RNASeq/analysis/analysis_rmDups/deseq_MSA_blood_RNASeq_rmDups.tab", row.names = FALSE, quote = FALSE, sep = '\t')
