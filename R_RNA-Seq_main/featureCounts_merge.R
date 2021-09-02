
# Setup -------------------------------------------------------------------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

featureCounts_directory <- 
  '/Users/andrzejeljaszewicz/Erasmus_UP/scripts_main/R_miRNA_main/featureCounts_Witek_mirnas'

output_file_dir <- 
  '/Users/andrzejeljaszewicz/Erasmus_UP/scripts_main/R_miRNA_main/miRNA_Witek_raw_counts.csv'

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# Dependencies ------------------------------------------------------------
library(purrr)
library(tidyverse)
library(readr)

# Load --------------------------------------------------------------------

# Change working dir so the listed files can be found
setwd(featureCounts_directory)

# List files in directory
f_files <-
  list.files(
    featureCounts_directory
  )

# Create function to merge featureCounts matrices for individual samples
read_in_feature_counts <- function(file) {
  cnt <- read_tsv(file, col_names = T, comment = "#")
  cnt <- cnt %>% dplyr::select(-Chr,-Start,-End,-Strand,-Length)
  return(cnt)
}

# Merge and create matrix
raw_counts <- map(f_files, read_in_feature_counts)

# Create object matrix
raw_counts_df <- purrr::reduce(raw_counts, inner_join)

# Save file as .csv
write_csv(
  raw_counts_df,
  output_file_dir
)
