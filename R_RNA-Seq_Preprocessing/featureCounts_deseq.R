
# Dependencies ------------------------------------------------------------
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(readr))

# Setup -------------------------------------------------------------------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# Read merged featureCounts matrix in
merged_featureCounts_matrix <-
  '/Users/andrzejeljaszewicz/Erasmus_UP/scripts_main/R_miRNA_main/miRNA_Witek_raw_counts.csv'

# Read merged featureCounts matrix in
featureCounts_matrix <-
  read_delim(
    merged_featureCounts_matrix
  ) %>%
  column_to_rownames(var = 'Geneid')

# Trim the colnames if neeeded
colnames(featureCounts_matrix) <- stringr::str_subset(colnames(featureCounts_matrix), '.fastq.gz_fastp.fastq.gzAligned.out.bam')

# Read meta or create the condition
# Meta file has to have 'samples' and 'condition' columns
meta_file <-
  ''

# Check if meta file already exist and set factors, but if not do it manually
if (file.exists(meta_file)) {
  meta <- read.csv(meta_file)
  # Set factor leves
  meta$condition <-
    factor(meta$condition, levels = c('NonAsthma', 'Asthma'))
} else {
  meta <- data.frame(samples = colnames(featureCounts_matrix),
                     condition = factor(rep(1:2, each = 100), levels = c(1,2)))
}

# Specify outpit directories and file names
de_results_out <- '/Users/andrzejeljaszewicz/Erasmus_UP/scripts_main/R_miRNA_main/miRNA_Witek_de_results.csv'

normalised_counts_out <- '/Users/andrzejeljaszewicz/Erasmus_UP/scripts_main/R_miRNA_main/miRNA_Witek_normalised_counts.csv'

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# Function for already merged featureCounts matrix ------------------------

#featurecounts_matrix <- read.csv('/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/nextflowScripts/featureCounts/featureCounts_matrix.txt',
#         sep = "\t", skip = 1L) %>%
#  dplyr::select(-c('Chr','Start','End','Strand','Length')) %>%
#  column_to_rownames(var = 'Geneid')

# Load --------------------------------------------------------------------

# Assign colnames as samples
colnames(featureCounts_matrix) <- meta$samples

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = featureCounts_matrix,
                              colData = meta,
                              design = ~ condition)

# Prefilter for low count genes
keep <- rowSums(counts(dds)) >= 20

# Perform filtering
dds_f <- dds[keep, ]

# Perform DE
dds_deseq <- DESeq(dds_f)

# Obtain DE results from DE object
res <- results(dds_deseq)

# Summarize
summary(res)

# Normalized counts using DESeq2's median of ratios method allowing for sample-sample comparisions
# (only other tool capable of this kind of normalization is EdgeR)
dds_norm <- estimateSizeFactors(dds)
counts_norm <- counts(dds_norm, normalized = TRUE)

# Prepare DE results for export
res_sroted <- res %>%
  data.frame() %>%
  rownames_to_column(var = 'id') %>%
  arrange(id)

# Prepare normalized matrix for export
counts_norm_sorted <- counts_norm %>%
  data.frame() %>%
  rownames_to_column(var = 'id') %>%
  arrange(id)

# Export data as .csv
write_csv(res_sroted, file = de_results_out)

write_csv(counts_norm_sorted, file = normalised_counts_out)
