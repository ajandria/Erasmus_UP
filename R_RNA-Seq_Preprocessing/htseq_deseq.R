directory <- "/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/htseq/htseq"
sampleFiles <- list.files(directory)
class(sampleFiles)
sampleFiles <- sort(sampleFiles)
sampleFiles
meta <- read.csv('/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/sampleTable2.csv')
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=meta$condition)
sampleTable
sampleTable$condition <- factor(sampleTable$condition, levels=c('NonAsthma','Asthma'))
sampleTable$condition
library(DESeq2)
dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory,design=~condition)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
res

dds_norm <- estimateSizeFactors(dds)
counts_norm <- counts(dds_norm, normalized=TRUE)

write.csv(res, file='/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/htseq/deseq_htseq.csv')
write.csv(counts_norm, file='/Users/andrzejeljaszewicz/Erasmus_UP/liam_data_dge_analysis/htseq/deseq_htseq_counts_norm.csv')
