
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon_with_inferential_replicates

# Dependencies ------------------------------------------------------------
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tximportData))
suppressPackageStartupMessages(library(dplyr))
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# Setup -------------------------------------------------------------------
meta <- readxl::read_excel('/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/Liam SIBRO data.xls', sheet = 2) %>% 
  mutate(Identifier = stringr::str_sub(Identifier, 0, 6))

metaCondition <- c('Non-asthma'="NonAsthma", 'Asthma'="Asthma")
meta$`Asthma Status` <- as.character(metaCondition[meta$`Asthma Status`])

dir <- '/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/nextflowScripts'

sample_names <- list.files('/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/nextflowScripts/salmon_fastp')

files <- file.path(dir, "salmon_fastp", sample_names, "quant.sf")

names(files) <- sample_names

all(file.exists(files))

txdb_gff <- makeTxDbFromGFF(file = '/Users/andrzejeljaszewicz/Erasmus_UP/nextflowTranscriptomics/geneModel/genes.gtf')

k <- keys(txdb_gff, keytype = 'TXNAME')

tx2gene <- select(txdb_gff, k, 'GENEID', 'TXNAME')

txi.salmon <- tximport(files, type = 'salmon', tx2gene = tx2gene)

sampleTable <- data.frame(samples = stringr::str_sub(sample_names, 35),
                          Identifier = stringr::str_sub(sample_names, 35, 40))

sampleTable2 <- sampleTable %>% 
  full_join(meta) %>% 
  filter(!is.na(samples)) %>% 
  dplyr::rename('condition' = "Asthma Status")

sampleTable2$condition <- factor(sampleTable2$condition, levels = c('NonAsthma', 'Asthma'))

# Load --------------------------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable2, ~condition)

keep <- rowSums(counts(dds)) >= 20

dds_f <- dds[keep,]

dds_deseq <- DESeq(dds_f)

res <- results(dds_deseq)

summary(res)

plotMA(res, ylim=c(-10,10))


plotMA(lfcShrink(res))



####
tpm_mat1 <- txi.salmon$counts/txi.salmon$length

tpm_mat2 <- t( t(tpm_mat1) * 1e6 / colSums(tpm_mat1) )
####

# Export ------------------------------------------------------------------
data_results_export <- data.frame(res) %>% 
  tibble::rownames_to_column(var = 'id') %>% 
  dplyr::select(id, everything()) %>% 
  mutate('temp' = NA) %>% 
  dplyr::rename('Compared samples: condition Asthma vs Non.asthma' = 'temp')

#write.csv(data_results_export, file="dge_liam.csv")
#write.csv(data.frame(tpm_mat2) %>% tibble::rownames_to_column(var = 'id'), file="TPM_liam_f.csv")
#write.csv(data.frame(counts(dds, normalized = F))%>% tibble::rownames_to_column(var = 'id'), file="raw_counts_liam_f.csv")

write.table(data_results_export,'dge_liam', row.names = FALSE)
write.table(data.frame(counts(dds, normalized = F))%>% tibble::rownames_to_column(var = 'id'),'raw_counts_liam_f.txt', row.names = FALSE)
write.table(data.frame(tpm_mat2) %>% tibble::rownames_to_column(var = 'id'),'TPM_liam_f.txt', row.names = FALSE)

# Raw counts QC -----------------------------------------------------------

dds_norm <- estimateSizeFactors(dds)
counts_norm <- counts(dds_norm, normalized=TRUE)

colnames(counts_norm) <- stringr::str_sub(sample_names, 35)

#write_csv(data.frame(counts_norm) %>% tibble::rownames_to_column(var='id'), file = 'deseq2_counts_norm.csv')

pdf('QC_DESeq2_normalized_counts.pdf')

# PCA -----------------------------

PCA_values <- prcomp(t(counts_norm))

library(factoextra)

options(ggrepel.max.overlaps = Inf)

fviz_eig(PCA_values)

groups_pca <- data.frame(condition = as.factor(sampleTable2$condition),
                         sample = sampleTable2$samples)

fviz_pca_ind(PCA_values,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T    # Avoid text overlapping,
)

#' @reorder PCA without 6 outliers

fviz_pca_ind(PCA_values,
             col.ind = groups_pca$condition, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             #addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

PCA_values_cut <- data.frame(counts_norm) %>% 
  dplyr::select(-c('SIB002','SIB018','SIB011','SIB016')) %>% 
  t() %>% 
  prcomp()

groups_pca_cut <- groups_pca %>% 
  filter(sample %in% rownames(PCA_values_cut$x))

fviz_pca_ind(PCA_values_cut,
             col.ind = groups_pca_cut$condition, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             #addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             title = "without 'SIB002','SIB018','SIB011','SIB016'"
)

PCA_values_cut2 <- data.frame(counts_norm) %>% 
  dplyr::select(-c('SIB002','SIB018','SIB011','SIB016','SIB012','SIB014_replicate_1')) %>% 
  t() %>% 
  prcomp()

groups_pca_cut2 <- groups_pca %>% 
  filter(sample %in% rownames(PCA_values_cut2$x))

fviz_pca_ind(PCA_values_cut2,
             col.ind = groups_pca_cut2$condition, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             #addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             title = "without 'SIB002','SIB018','SIB011','SIB016','SIB012','SIB014_replicate_1'"
             
)

# Distance Matrix --------------------

#' @reorder samples; generate corr plot can be with dist matrix

mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]

library("RColorBrewer")
sampleDists <- dist(t(counts_norm))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
rownames(sampleDistMatrix) <- paste(sampleTable2$condition,sampleTable2$samples, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
sampleDistMatrix<-sampleDistMatrix[order(rownames(sampleDistMatrix)),order(colnames(sampleDistMatrix))]
#sampleDistMatrix <- mat_ord(sampleDistMatrix)
pheatmap(sampleDistMatrix,
         #clustering_distance_rows=sampleDists,
         #clustering_distance_cols=sampleDists,
         col=colors,
         show_rownames = T, show_colnames = T,
         cluster_rows = F,
         cluster_cols = F)


# Boxplots ------------------------------------------------------

#' @reorder samples

library(reshape2)

df = melt(counts_norm, variable.name = "Samples", 
          value.name = "count") %>% # reshape the matrix
  dplyr::rename('samples' = 'Var2') %>% 
  full_join(sampleTable2) %>% 
  arrange(samples) %>% 
  arrange(condition)

df$samples <- factor(df$samples, levels = unique(df$samples))

ggplot(df, aes(x=samples, y=log2(count+1), fill=condition))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Density -----------------------------------------------------------------

ggplot(df, aes(x = log(count+1), colour = samples, fill = samples)) + 
  #ylim(c(0, 0.17)) + 
  geom_density(alpha = 0.2, size = 1.25) + 
  facet_wrap(~ condition) + 
  theme(legend.position = "top") + 
  xlab(expression(log[2](count + 1)))+ 
  theme(legend.title = element_text(size = 7, face = 'bold'), 
        legend.text = element_text(size = 7, face = 'bold')) +
  theme(legend.key.size = unit(1,"line"))



# correlation -------------------------------------------------------------

dist_mat <- dist(t(counts_norm),
          diag = TRUE,
          upper = TRUE) %>% 
  as.matrix()

cor_mat <- cor(dist_mat)


rownames(cor_mat) <- paste(sampleTable2$condition,sampleTable2$samples, sep="-")
colnames(cor_mat) <- rownames(cor_mat)
cor_mat2<-cor_mat[order(rownames(cor_mat)),order(colnames(cor_mat))]

library(corrplot)

corrplot(cor_mat2, method = 'color', tl.col = 'black', tl.cex = 0.8,
                   col=colorRampPalette(c("red","white","blue"))(300))

#dev.off()


# cluster -----------------------------------------------------------------
d_iris <- dist(t(counts_norm)) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
iris_species <- rev(levels(sampleTable2$condition))

library(dendextend)
dend <- as.dendrogram(hc_iris)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:150)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=2)#, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
labels_colors(dend) <-
  colorspace::rainbow_hcl(2)[sort_levels_values(
    as.numeric(sampleTable2$condition)[order.dendrogram(dend)]
  )]

# We shall add the flower type to the labels:
labels(dend) <- paste(as.character(sampleTable2$condition)[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
par(mar = c(2,2,2,14)+0.1, xpd = NA)
plot(dend, 
     main = "Clustered data set", 
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = iris_species, fill = colorspace::rainbow_hcl(2))
#the_bars <- ifelse(as.character(sampleTable2$condition), "grey", "gold")
sampleTable2$colorCondition <- ifelse(as.character(sampleTable2$condition)=='Asthma', '#E495A5','#39BEB1')
colored_bars(colors = sampleTable2$colorCondition, dend = dend, rowLabels = "",
             add = TRUE, horiz = TRUE)

dev.off()
############

# Create the dendrogram, use default options
dend_mtcars <- t(counts_norm) %>% 
  dist %>% 
  hclust() %>% 
  as.dendrogram

# Set the plot margin: bottom, left, top & right
par(mar = c(10, 3, 3, 4) + 0.1,
    xpd = NA) # allow content to go into outer margin 

# Plot
plot(dend_mtcars)

# Setup the color bar based on $am & $vs
the_bars <- ifelse(as.character(sampleTable2$condition), "grey", "gold")
colored_bars(colors = the_bars, dend = dend_mtcars, rowLabels = c("Corectness"))

# Add the legend manually
legend("topright", legend = c('0', '1'), pch = 15, pt.cex = 3, cex = 1.5, bty = 'n',
       inset = c(-0.1, 0), # place outside
       title = "Status", 
       col = c('beige', 'firebrick3'))




# 3d pca ------------------------------------------------------------------

library(rgl)
library(pca3d)

pca3d(PCA_values, group=groups_pca$condition,show.group.labels=T,
      fancy = T)
#rglwidget(width = 1920, height = 1080)

# PCA ---------------------------------------------------------------------
library(ggplot2)

vsd <- vst(dds)

plotPCA(vsd, intgroup=c("condition"), ntop= 500) + 
  geom_text(aes(label=vsd$Identifier),vjust=2,check_overlap = T, size = 3)+
  theme_bw()
plotPCA(vsd, intgroup=c("condition"), ntop= 1500) + 
  geom_text(aes(label=vsd$Identifier),vjust=2,check_overlap = T, size = 3)+
  theme_bw()
plotPCA(vsd, intgroup=c("condition"), ntop= nrow(vsd)) + 
  geom_text(aes(label=vsd$Identifier),vjust=2,check_overlap = T, size = 3)+
  theme_bw()







