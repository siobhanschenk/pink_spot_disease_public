##### set up #####
library(tidyverse)
library(phyloseq)
library(plyr)
library(qualpalr)
library(ggh4x)
library(vegan)
library(ggpattern)
library(DESeq2)
library(ggpubr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
                                   axis.text.x = element_text(colour = "black", face = "bold", size = 8),
                                   legend.text = element_text(size = 8, face ="bold", colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 12, colour = "black", face = "bold"),
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


setwd("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/imput_data")
path="C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/output/"

## read in data
all = readRDS("all_studies_filtered_phyloseq.RDS")
all
all@sam_data$rd_no_off = sample_sums(all)

## metadata
meta = read.csv("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/pink-spot_metadata.csv")

meta = meta |> column_to_rownames(var="rn")

all = phyloseq(sample_data(meta),
               otu_table(all@otu_table),
               tax_table(all@tax_table))
all

##### SUBSET SAMPLES #####

cs = subset_samples(all,
                      pink_health %in% c( "pink_healthy", "pink_diseased") &
                      kelp_species=="Saccharina latissima" & treatment =="control")
#View(cs@sam_data)

project_data = prune_taxa(taxa_sums(cs)>0, cs)

tax = as.data.frame(as.matrix(project_data@tax_table))


alpha <- 0.01 #your significance threshold for MULTIPLE TEST CORRECTED pvals

####DESeq2 analysis preparation####
## add 1 to every cell
project_data@otu_table = project_data@otu_table +1

#OPTIONAL/IF NEEDED:assign variables as either factors or numeric
sample_data(project_data)$pink_health <- factor(sample_data(project_data)$pink_health, 
                                                      levels=c("pink_diseased", "pink_healthy"))

#IMPORTANT: factors need to be watched carefully, recoding variables as factors can re-arrange data w/in the data frame if you have previously ordered it in a specific way

####DESeq: taxa correlated with a binary factor ("variable3"), controlling for two other factors ("variable1" and "variable2") ####
#parameters up front

## OPTION 1: DESEQ RUN USING LRT MODEL ##
dds.var3 <- phyloseq_to_deseq2(project_data, design = ~ pink_health)
dds.var3 <- DESeq(dds.var3, test = "Wald", fitType = "parametric") 
resultsNames(dds.var3) #check results names for the name of the contrast you'd like to examine


#### save output "pink_health_pink_healthy_vs_pink_diseased" #####

## sig output
res.PDvsPH <- results(dds.var3, cooksCutoff = FALSE, alpha = alpha, pAdjustMethod = "BH", 
                    name="pink_health_pink_healthy_vs_pink_diseased", 
                    altHypothesis = "greaterAbs") 
summary(res.PDvsPH, alpha = alpha)

#reorder the results table by adjusted p-value and remove any "NA" entries
res_p_ordered <- res.PDvsPH[order(res.PDvsPH$padj, na.last = NA), ]   
#filter out any results which have a adjust p-value less than alpha 
res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj <= alpha), ]  
res.PDvsPH <- cbind(as(res_p_ordered_filt, "data.frame"), 
                                  as(tax_table(project_data)[rownames(res_p_ordered_filt), ], "matrix"))
res.PDvsPH$comparison = "res.PDvsPH"

res.PDvsPH = res.PDvsPH %>%
  rownames_to_column(var="asv")

write.csv(res.PDvsPH,  paste0(path,"sig_deseq_output.csv"))

###### OCCURENCE OF SEIG ASV IN HC AND PD ######
cs = subset_taxa(cs, asv_id %in% c(res.PDvsPH$asv_id))
spool = dephyloseq(cs)

## add presence anscence column
spool$ra = as.numeric(spool$asv_abundance)/as.numeric(spool$rd_no_off)
spool$presabs = as.numeric(ifelse(spool$ra>0, "1", "0"))

pd = subset(spool, spool$pink_health=="pink_diseased")
pd = ddply(pd, c("asv_id"),
           summarise,
           nsamples = sum(presabs))
pd = subset(pd, pd$nsamples>3)

## keep taxa that made the cut
sigtax = subset(spool, spool$asv_id %in% c(pd$asv_id))

## include deseq output
sigtax = left_join(sigtax, res.PDvsPH)

##### heatmap #####
sigtax$ra[sigtax$ra==0]<-NA
sigtax$plotnames = paste0(sigtax$genus, " (", sigtax$asv_id,")")

## order pink health
sigtax$pink_health = factor(sigtax$pink_health, levels=c("pink_healthy", "pink_diseased"))

## order by log2 fold cange
sigtax = sigtax[order(sigtax$log2FoldChange),]
## get list of factors in order
taxa.order = as.list(c(unique(sigtax$plotnames)))
## set plot_names levels
plot.order = unlist(taxa.order)
sigtax$plotnames = factor(sigtax$plotnames, levels=c(plot.order))

## make healtmap
heat=ggplot(sigtax, aes(x=paste(as.character(Row.names), location, date_sampled), 
                        y=fct_rev(plotnames), fill=sqrt(ra)))+
  geom_tile()+
  facet_nested(.~year_sampled+pink_health, space="free", scales="free")+
  scale_fill_gradient2(low="#5ef4ff", mid="#c85bd6", midpoint=0.5, high="#ff0000", na.value = "white")+
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y="", x="Samples")
heat

point=ggplot(sigtax, aes(x=log2FoldChange, y=fct_rev(plotnames)))+
  geom_point(cex=4)

ggarrange(point, heat)

ggsave(paste0(path, "deseq_heatmap.pdf"), height=7, width=11.2, units="in")
