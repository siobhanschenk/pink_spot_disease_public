##### set up #####
library(tidyverse)
library(phyloseq)
library(plyr)
library(qualpalr)
library(ggh4x)
library(vegan)
library(ggpattern)
library(ggpubr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12),
                                   axis.text.x = element_text(colour = "black", size = 8),
                                   legend.text = element_text(size = 8, colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   # axis.ticks = element_blank(),
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
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
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

## read in data
setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/imput_data")

## read in data
all = readRDS("all_studies_filtered_phyloseq.RDS")

## metadata
meta = read.csv("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/pink-spot_metadata.csv")

meta = meta |> column_to_rownames(var="rn")

all = phyloseq(sample_data(meta),
               otu_table(all@otu_table),
               tax_table(all@tax_table))
all


##### SUBSET SAMPLES #####

spool = subset_samples(all,
                    pink_health %in% c( "pink_healthy", "pink_diseased") &
                      kelp_species=="Saccharina latissima" & treatment =="control")

other = subset_samples(all, treatment=="control" & sample_target %in% c("seawater", "airline") & diseased_tank =="pink_symptomatic")
#View(other@sam_data)

cs = merge_phyloseq(spool, other)
cs

cs = prune_taxa(taxa_sums(cs)>0, cs)
cs

cs@sam_data$rd_no_off = sample_sums(cs)
mean(cs@sam_data$rd_no_off)
sum(cs@sam_data$rd_no_off)

samplist = rownames(cs@sam_data)

#### SET UP GROUPS FOR NURSERY KELP ####

## group at genus level
cs.gen = tax_glom(cs, taxrank = "genus")
## calculate read depth
cs.gen@sam_data$rd_no_off_with_algicola = sample_sums(cs.gen)

## only algicola
only.algicola = subset_taxa(cs.gen, genus=="Algicola")

## calculate read depth
cs.gen@sam_data$rd_no_off = sample_sums(cs.gen)

## remove algicola
cs.gen = subset_taxa(cs.gen, genus!="Algicola")

## get data out of phyloseq
kelpsubdf = dephyloseq(cs.gen)


## make groups
grouplist = c(unique(kelpsubdf$pink_health))


## calculate relative abundance
kelpsubdf$ra = as.numeric(kelpsubdf$asv_abundance)/as.numeric(kelpsubdf$rd_no_off)

## make plotnames
kelpsubdf$plotnames = paste0(kelpsubdf$order,"; ", kelpsubdf$genus)

## summarize data by taxaplot group type. 
kelpsub.sum = ddply(kelpsubdf, c("pink_health", "plotnames", "sample_target"),
                summarise,
                sum = sum(ra))

## sort data by relative abundance. This is how the loop will pick the mos tabundant taxa
sorted = kelpsub.sum[order(-kelpsub.sum$sum),]


#### PICK TOP TAXA FOR EACH GROUP ####
## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in grouplist) {
  for(j in i) {
    
    ## subset dataframe by samples
    #!# Remeber to change te substrate to your group! 
    sample = subset(sorted, sorted$pink_health %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:10),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}


##### SUBSET AND JOIN DATAFRAMES #####
toplist = c(unique(top.df$plotnames))

## calculate read depth
cs.gen@sam_data$rd_no_off = sample_sums(cs.gen)
## get out of phyloseq
kelpdf = dephyloseq(cs.gen)
## make plotnames column
kelpdf$plotnames = paste0(kelpdf$order,"; ", kelpdf$genus)
## only keep target taxa
kelpdf.sub = subset(kelpdf, kelpdf$plotnames %in% c(toplist))


## add info about groups from taxaplot
kelpsubtop = full_join(kelpsubdf, top.df)
kelpsubtop = kelpsubtop[,c("Row.names", "plotnames", "pink_health", "sample_target")]

kelpdf.taxaplot = right_join(kelpsubtop, kelpdf.sub)

## calculate relative abundance
kelpdf.taxaplot$ra = as.numeric(kelpdf.taxaplot$asv_abundance)/as.numeric(kelpdf.taxaplot$rd_no_off)

## group the dataframe to get mean RA for each taxa
tp = ddply(kelpdf.taxaplot, c("plotnames", "pink_health", "Row.names", "tank_id",
                              "dual_swab", "sample_target","year_sampled", "location", "area_swabbed", "date_sampled"),
           summarise,
           meanra = mean(ra),
           sdra = sd(ra))

## calculate others
tpothers = ddply(tp,  c("pink_health", "Row.names","tank_id", "dual_swab",
                        "sample_target","year_sampled", "location", "area_swabbed","date_sampled"),
                 summarise,
                 sumra = sum(meanra))
tpothers$meanra = 1-tpothers$sumra
tpothers$plotnames = "Others"

alldata = full_join(tp, tpothers)

####GET COLORS FOR TAXAPLOT #####

# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))

# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(3)

# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")

# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")

# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")

# 6. Join color list and taxa names
topcolors = cbind(tops, hex)

# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"

# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

## ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$meanra),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## order disease
alldata$pink_health = factor(alldata$pink_health, levels=c("pink_healthy", "pink_diseased", "Symp.Env"))

## set plot_names levels
plot.order = unlist(plot.order)


##### MAKE PLOT #####
#alldata = subset(alldata, alldata$pink_health!="NA")


## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))

alldata$pink_health = ifelse(alldata$sample_target %in% c("seawater", "airline"), "Symp.Env.", alldata$pink_health)

bp=ggplot(alldata, aes(x=paste(location, date_sampled, as.character(Row.names)),
               y=meanra,
               fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~year_sampled+pink_health+sample_target, scales="free", space="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(ncol=2))
bp

only.algicola.df = dephyloseq(only.algicola)
only.algicola.df$ra = as.numeric(only.algicola.df$asv_abundance)/as.numeric(only.algicola.df$rd_no_off_with_algicola)
only.algicola.df$presabs = ifelse(only.algicola.df$ra==0, "abs", "pres")

only.algicola.df$sample_target = factor(only.algicola.df$sample_target, levels=c("spool", "airline", "seawater"))
only.algicola.df$pink_health = factor(only.algicola.df$pink_health, levels=c("pink_healthy", "pink_diseased","airline", "seawater"))

ggplot(only.algicola.df, aes(x=paste(location, date_sampled, as.character(Row.names)),
                    y=paste(genus, "   spacing"),
                    size=sqrt(ra), alpha=presabs, color=genus))+
  geom_point()+
  facet_nested(.~year_sampled+pink_health, scales="free", space="free")+
  scale_color_manual(values=c("magenta"))+
  scale_size(range=c(2,8),
             breaks=c(0,0.2, 0.4, 0.6, 0.8),labels=c("0","0.2","0.4","0.6","0.8"),guide="legend")

only.algicola.df$ra = as.numeric(ifelse(only.algicola.df$ra>0, only.algicola.df$ra, "NA"))

dp=ggplot(only.algicola.df, aes(x=paste(location, date_sampled, as.character(Row.names)),
                             y=paste(genus, "   spacing"),
                             fill=sqrt(ra)))+
  geom_tile()+
  facet_nested(.~year_sampled+pink_health, scales="free", space="free")+
  scale_fill_gradient(low="#b4b4b4", high="black", na.value = "white")+
  theme(axis.text.x = element_blank())

ggarrange(dp,bp, ncol=1, heights=c(0.4, 1))

ggsave("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/output/taxaplot_NoAlgicola.pdf",
       width=11.5, height=8.5, units="in")






