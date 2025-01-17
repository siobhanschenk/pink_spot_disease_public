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

## read in data
all = readRDS("all_studies_rarefied_phyloseq.RDS")
all

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
                      kelp_species=="Saccharina latissima" & treatment =="control" &
                      location =="prince_rupert")


nur@sam_data$rd_filt = sample_sums(nur)

##### CREATE A COUMN THAT IS THE RELATIVE ABUNDANCE OF ALGICOLA ######
ag = subset_taxa(nur, genus=="Algicola")
ag@sam_data$rd_algi = sample_sums(ag)

metafull = as.data.frame(as.matrix(nur@sam_data))
agmeta = as.data.frame(as.matrix(ag@sam_data))
agmeta = select(agmeta, rd_algi,rd_filt)
agmeta$ra_ag = as.numeric(agmeta$rd_algi)/as.numeric(agmeta$rd_filt)

metaboth = merge(metafull, agmeta, by=0)

metaboth = metaboth |> column_to_rownames(var="Row.names")

nuralgi = phyloseq(sample_data(metaboth),
                   otu_table(nur@otu_table, taxa_are_rows = F),
                   tax_table(nur@tax_table))

#### NMDS NURSERY - WITH ALGICOLA ####
#nuralgi2021 = subset_samples(nuralgi, year_sampled=="2021")
nuralgi2022 = subset_samples(nuralgi, year_sampled=="2022")


nur.ord2022 = ordinate(nuralgi2022, "NMDS","bray")
nur.ord2022
jwa2022=plot_ordination(nuralgi2022, nur.ord2022, color="pink_health", shape="location")+
  geom_point(size=6, stroke=1.5)+
  scale_color_manual(values=c("magenta", "#be9125"))+
  labs(color="substrate", shape="year")+
  scale_shape_manual(values=c(15, 17, 8))
jwa2022

#### PERMANOVS WITH ALGICOLA SET UP ####
# get data frames
otu.all = as.data.frame(as.matrix(nur@otu_table))
meta.all = as.data.frame(as.matrix(nur@sam_data))
## merge the data
all.mo = merge(meta.all, otu.all, by=0)
## break up by year
#all.mo2021 = subset(all.mo, all.mo$year_sampled=="2021")
all.mo2022 = subset(all.mo, all.mo$year_sampled=="2022")
## count colum numbers
metacols = ncol(meta.all)+1

#### PERMANOVA WITH ALGICOLA 2022 #####
## run permanova 
perm.all= adonis2(all.mo2022[,-c(1:metacols)] ~ as.character(pink_health), 
                  data=all.mo2022, method = "bray")
perm.all

# ## post hoc test
all.pairwise = pairwise.adonis(all.mo2022[,-c(1:metacols)], as.factor(all.mo2022$pink_health))
all.pairwise


nur2022 = subset_samples(nur, year_sampled=="2022")
bdt <- phyloseq::distance(nur2022, method = "bray")
sample_df <- data.frame(sample_data(nur2022))
st <- betadisper(bdt, sample_df$pink_health) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bst=permutest(st) 
bst



#### NMDS NURSERY - NO ALGICOLA ####
nurNOalgi = subset_taxa(nuralgi,!(genus %in% c("Algicola")))


#na2021 = subset_samples(nurNOalgi, year_sampled=="2021")
na2022 = subset_samples(nurNOalgi, year_sampled=="2022")

nur.ord2022 = ordinate(na2022, "NMDS","bray")
nur.ord2022
jna2022=plot_ordination(na2022, nur.ord2022, color="pink_health", shape="location")+
  geom_point(size=6, stroke=1.5)+
  scale_color_manual(values=c("magenta", "#be9125"))+
  labs(color="substrate", shape="year")+
  scale_shape_manual(values=c(15, 17, 8))
jna2022


#### PERMANOVA - NO ALGICOLA 2022 ######
# get data frames
otu2022 = as.data.frame(as.matrix(na2022@otu_table))
meta2022 = as.data.frame(as.matrix(na2022@sam_data))
## merge the data
all.mo = merge(meta2022, otu2022, by=0)
## count colum numbers
metacols = ncol(meta2022)+1

## run permanova 
perm.all= adonis2(all.mo[,-c(1:metacols)] ~ as.character(pink_health), 
                  data=all.mo, method = "bray")
perm.all

# betadispersion test
bdt <- phyloseq::distance(na2022, method = "bray")
sample_df <- data.frame(sample_data(na2022))

st <- betadisper(bdt, sample_df$pink_health) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
bst=permutest(st) 
bst

# ## post hoc test
all.pairwise = pairwise.adonis(otu2022, as.factor(meta2022$pink_health))
all.pairwise


#### POST HOC NO ALGI 2021 #####
nur2021 = subset_samples(nurNOalgi, year_sampled=="2021")
otu.all = as.data.frame(as.matrix(nur2021@otu_table))
meta.all = as.data.frame(as.matrix(nur2021@sam_data))

# ## post hoc test
all.pairwise = pairwise.adonis(otu.all, as.factor(meta.all$pink_health))
all.pairwise


##### ARRANGE PLOTS ######
ggarrange(jwa2021, jwa2022,jna2021, jna2022)

ggsave(filename = paste0(path, "spool_nmds_bray.pdf"), height =9.4, width=12.9, units = "in")
