##### set up #####
library(tidyverse)
library(phyloseq)
library(ggh4x)
library(plyr)
library(qualpalr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
                                   axis.text.x = element_text(colour = "black", face = "bold", size = 12),
                                   legend.text = element_text(size = 8, face ="bold", colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   axis.ticks = element_blank(),
                                   panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank()))

setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/imput_data")

## read in data
all = readRDS("all_studies_filtered_phyloseq.RDS")
all

all@sam_data$rd_no_off = sample_sums(all)

chl = read.csv("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/imput_data/filtering_chloroplast_loss.csv")
chl = chl[,c("X", "no_chloroplast_reads", "rd_post_merge")]

##### subset samples ######

cs = subset_samples(all, study=="cascadia")

meta = as.data.frame(as.matrix(cs@sam_data)) |> rownames_to_column(var="X")
meta = meta[,c(-3)]

reads = right_join(chl, meta)


#### get number of chloroplast per library #####
reads$number_chloroplast_reads = as.numeric(reads$rd_post_merge)-as.numeric(reads$no_chloroplast_reads)
reads$ra_number_chloroplast_reads = reads$number_chloroplast_reads/as.numeric(reads$rd_post_merge)

reads = separate(reads, col="X", into=c("illumina_name", "library_id"), sep="-_-")

reads = subset(reads, reads$library_id !="NA")

reads$illumina_name = paste0(reads$illumina_name, "-_-", reads$library_id)

write.csv(reads, "cascadia_read_number_chloroplast_blocker_table.csv")


ggplot(reads, aes(x=sample_target, y=ra_number_chloroplast_reads, color=library))+
  geom_boxplot()

ggplot(reads, aes(x=sample_target, y=as.numeric(rd_no_off), color=library))+
  geom_boxplot()
