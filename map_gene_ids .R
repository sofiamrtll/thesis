
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# 
#BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)
#install.packages("tidyverse")
#install.packages("readxl")
#install.packages("writexl")
library(readxl)
library(writexl)

# taken from: https://bioinformatics.stackexchange.com/questions/15573/what-is-the-best-way-to-programmatically-convert-ensembl-ids-from-one-release-to

mart37 = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 host    = "https://grch37.ensembl.org",
                 path    = "/biomart/martservice",
                 dataset = "hsapiens_gene_ensembl")

# default host is grch38://
mart38 = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 path    = "/biomart/martservice",
                 dataset = "hsapiens_gene_ensembl")


#read.xlsx("HEK293.xlsx", rowIndex = 2:57908, colIndex = 1)
# example input list
#ensembl_gene_ids <- c("ENSG00000198062", "ENSG00000130538", "ENSG00000198445", "ENSG00000172967", "ENSG00000215568")
#example = readxl_example("C:/Users/sofia/Desktop/THESIS/PROJECT/HEK293.xlsx")
ensembl_gene_ids2 <- read_excel("/cloud/project/expressed_gene_ids_HEK.xlsx")
ensembl_gene_ids3 = data.frame( ensembl_gene_ids2)
ensembl_gene_ids3
ensembl_gene_ids = ensembl_gene_ids3[,1]
ensembl_gene_ids
columns = c('ensembl_gene_id', 'ccds')

out37 <- getBM(attributes = columns,
               filters = c('ensembl_gene_id'),
               values = ensembl_gene_ids,
               mart = mart37) %>% #empty string CCDS ID's
  filter(ccds != "")

#NOTE:: to check which filters can be used to query mart37, run this command:
#View(listFilters(mart37))

# query on CCDS as the key to map gene ids between releases
# https://www.ensembl.org/info/genome/genebuild/ccds.html
out38 <- getBM(attributes = columns,
               filters = c('ccds'),
               values = list(out37$ccds),
               mart = mart38)

merge2 = merge(x = out37, y = out38, by = 'ccds', all = TRUE)
merge = data.frame(merge2)

write_xlsx(merge ,"/cloud/project/output.xlsx")

