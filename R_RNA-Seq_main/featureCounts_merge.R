library(purrr)
library(tidyverse)
library(readr)

setwd("/Users/andrzejeljaszewicz/Erasmus_UP/mirnas/input_mirnas/mirna_input_featureCounts")

f_files<- list.files("/Users/andrzejeljaszewicz/Erasmus_UP/mirnas/input_mirnas/mirna_input_featureCounts")

read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names =T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
  return(cnt)
}

raw_counts<- map(f_files, read_in_feature_counts)

raw_counts_df<- purrr::reduce(raw_counts, inner_join) 

write_csv(raw_counts_df, '/Users/andrzejeljaszewicz/Erasmus_UP/mirnas/input_mirnas/input_mirnas_liam_featureCounts_matrix.txt')
