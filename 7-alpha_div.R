##### set up #####
library(tidyverse)
library(phyloseq)
library(plyr)
library(qualpalr)
library(ggh4x)
library(vegan)
library(microbiome)
library(pairwiseAdonis)
library(car)
library(ggpubr)
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
                                   panel.grid.minor = element_blank()
                             ))


## load functions
dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  ## if your metadta is empty after running this, you need to use 
  otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="ASVid")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "ASVid", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}


setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/imput_data")
path="C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/output/"


## replace metadata
all = readRDS("all_studies_rarefied_phyloseq.RDS")

## metadata
meta = read.csv("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/pink-spot_metadata.csv")

meta = meta |> column_to_rownames(var="rn")

all = phyloseq(sample_data(meta),
               otu_table(all@otu_table),
               tax_table(all@tax_table))
all

##### SUBSET SAMPLES #####

nur = subset_samples(all,
                     pink_health %in% c( "pink_healthy", "pink_diseased") &
                       kelp_species=="Saccharina latissima" & treatment =="control")


nur@sam_data$rd_filt = sample_sums(nur)

#### CALUCLATE DIVERSITY INDICIES ######

## extract data
meta = as.data.frame(as.matrix(nur@sam_data))
otu = as.matrix(as.data.frame(as.matrix(nur@otu_table)))

## get species shan for each sample
meta$rich = vegan::specnumber(otu)
meta$shan = vegan::diversity(otu)

###### GET MEAN VALUES ######
meta.sum = ddply(meta, c("pink_health", "year_sampled"),
                 summarise,
                 meanrich = mean(rich),
                 meanshan = mean(shan))


meta = full_join(meta, meta.sum)

##### PLOT #####
spool = subset(meta, meta$sample_target=="spool")

spool$pink_health = factor(spool$pink_health, levels=c("pink_healthy","pink_diseased"))

ggplot(spool, aes(x=pink_health, y=as.numeric(shan), shape=location, color=pink_health))+
  geom_point(size = 4, stroke=1.1, position=position_jitter(width=0.2, height=0))+
  scale_color_manual(values=c("#be9125", "magenta"))+
  labs(color="substrate", shape="year")+
  ylim(0,4.7)+
  facet_nested(.~year_sampled, scales="free", space="free")+
  scale_shape_manual(values=c(19, 15, 17, 8))

ggsave(paste0(path,'alpha_div.pdf'), height=4.3, width=5.7, units="in")


##### ANOVA #####
spool2021 = subset(spool, spool$year_sampled=="2021")
spool2022 = subset(spool, spool$year_sampled=="2022")

a1=aov(shan~pink_health, data=spool2021)
summary(a1)
#plot(a1)
leveneTest(a1)
TukeyHSD(a1)


a1=aov(shan~pink_health, data=spool2022)
summary(a1)
#plot(a1)
leveneTest(a1)
TukeyHSD(a1)
