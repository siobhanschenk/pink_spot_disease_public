library(phyloseq)
library(tidyverse)
library(zoo)

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/imput_data")

tax = read.table("cascadia_all_taxonomy_noNAs_SILVAv138_16s.txt", header=T, row.names = "row_names")
seq = readRDS("seqtab_nochim_cascadia_all_studies_full_length.RDS")

seq = seq |> 
  as.data.frame() |> 
  rownames_to_column(var="rn")

diseaseonly = subset(seq, seq$rn %in% c(samplist))

diseaseonly = diseaseonly |> 
  rownames_to_column(var="otu_sum") |> 
  column_to_rownames(var="rn") |> 
  as.matrix() |> 
    t()

class(diseaseonly) <-"numeric"

sum(diseaseonly)

diseaseonly.sum= as.data.frame(rowSums(diseaseonly))
diseaseonly.sum


both = cbind(diseaseonly.sum, diseaseonly)

both = subset(both, both$`rowSums(diseaseonly)`>0)
dim(both)
