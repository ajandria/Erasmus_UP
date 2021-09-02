
setwd('/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/nextflowScripts')

library(magrittr)
library(dplyr)
library(tibble)
library(stringr)
library(DESeq2)

meta <- read.csv('/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/sampleTable2.csv')

meta$condition <- factor(meta$condition, levels = c('NonAsthma', 'Asthma'))

featurecounts_matrix <- read.csv('/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/nextflowScripts/featureCounts/featureCounts_matrix.txt',
         sep = "\t", skip = 1L) %>% 
  dplyr::select(-c('Chr','Start','End','Strand','Length')) %>% 
  column_to_rownames(var = 'Geneid')

colnames(featurecounts_matrix) <- meta$samples

dds <- DESeqDataSetFromMatrix(countData = featurecounts_matrix,
                              colData = meta,
                              design = ~condition)

keep <- rowSums(counts(dds)) >= 20

dds_f <- dds[keep,]

dds_deseq <- DESeq(dds_f)

res <- results(dds_deseq)

summary(res)

dds_norm <- estimateSizeFactors(dds)
counts_norm <- counts(dds_norm, normalized=TRUE)

res_sroted <- res %>% 
  data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  arrange(id)

counts_norm_sorted <- counts_norm %>% 
  data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  arrange(id)

write.csv(res_sroted, file = '/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/deseq_featureCounts_dge.csv')

write.csv(counts_norm_sorted, file = '/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/deseq_featureCounts_normalised_counts.csv')








