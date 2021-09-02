
library(magrittr)
library(tidyr)
library(stringr)
library(dplyr)

refFlat <- read.delim('/Users/andrzejeljaszewicz/noChr_refFlat.txt', header = F)

gtf <- read.delim('/Users/andrzejeljaszewicz/mirna_lincrna_Homo_sapiens.GRCh37.87.chr.gtf', header = F)

gtf_name1 <- sub(".*gene_name ", "", gtf$V9)

gtf_name2 <- sub(";.*", "", gtf_name1)


refFlat_mirna_lincrna <- refFlat %>% 
  filter(V1 %in% gtf_name2)

write.table(refFlat_mirna_lincrna, file = 'refFlat_mirna_lincrna.txt', sep = "\t",
            quote = F, col.names = F, row.names = F)
