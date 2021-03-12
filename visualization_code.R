############################################################################
# March 12th, 2021
# Heather Deel
# Compilation of code used to generate figures in bone_2020_deel manuscript
############################################################################

### libraries ###
library(tidyverse)
library("lme4")
library("ggplot2")
library("HLMdiag")
library("DHARMa")
library("car")
library("Matrix")
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(randomcoloR)
library(tidyr)
library(ggpubr)
library(randomcoloR)
library(tibble)
library(ggpubr)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(plyr)
library(rstatix)
library(rlang)

### Figure 1 ###

# 16S linear mixed effects

setwd("/Users/alex/Dropbox/PMI_3_analyses/bone/01_16S/07_longitudinal/72f27eb0-57fb-42ee-a20f-68912e6e2064/data")

# heather pathway
setwd("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/07_longitudinal/72f27eb0-57fb-42ee-a20f-68912e6e2064/data")

faithdata<-read.delim("raw-data.tsv",header=TRUE,sep="\t")

#remove two potential outliers
#11553.STAFS.2016.065.R11 [44]
#11553.SHSU.STAFS2016.011.L08 [10]

faithdata2<-faithdata[-c(44,10),]

#look at linearity
ggplot(faithdata2,aes(faith_pd,residual,color=season))+
  geom_point()

#homogeneity of variance
faithdata2$residualsquare <- faithdata2$residual^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.Model <- lm(residualsquare ~ host_subject_id , data=faithdata2) #ANOVA of the squared residuals
anova(Levene.Model) #p-value is 0.04; residuals do not have equal variance across subjects

#model with random effects and random intercepts
Model.F<-lmer(faith_pd ~ ADD_0*season
              + (1+ADD_0|host_subject_id),
              data=faithdata2, REML=TRUE)
#model with random intercepts
Model.F2<-lmer(faith_pd ~ ADD_0*season
               + (1|host_subject_id), 
               data=faithdata2, REML=TRUE)
library(lmtest)
lrtest(Model.F,Model.F2)

Plot.Model.F <- plot(Model.F)

#normality
require("lattice")
qqmath(Model.F, id=0.05)

#plot faith's pd
library(ggpubr)
A<-ggscatter(faithdata2, x = "ADD_0", y = "faith_pd",
             add = "reg.line",                         # Add regression line
             conf.int = TRUE,                          # Add confidence interval
             color = "season", palette = "jco",           # Color by groups "cyl"
             shape = "season"                             # Change point shape by groups "cyl"
)+
  xlab("Accumulated Degree Day")+
  ylab("Faith's pd")+
  theme(axis.text.x=element_text(angle = 90, vjust=0.5, hjust = 1))+
  scale_y_continuous(breaks=c(10,20,30,40,50),limits=c(0,50))

A

# add correlation coefficient and p-value to plot
A + stat_cor(method = "pearson")

ggsave("Faithpd16S.png", height=5, width=5, device="png", dpi=500)

# 18S linear  mixed effects

faithdata3<-read.delim("/Users/alex/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/05_core_metrics/core-metrics-results_hostID/longitudinal/39bc961c-e5f6-443d-a217-b656136174c7/data/raw-data.tsv",header=TRUE,sep="\t")

# pathway for Heather
faithdata3<- read.delim("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/05_core_metrics/core-metrics-results_hostID/longitudinal/39bc961c-e5f6-443d-a217-b656136174c7/data/raw-data.tsv",header=TRUE,sep="\t")

#remove one potential outliers
#STAFS.2016.064.R12 [22]
faithdata3<-faithdata3[-c(22),]

#look at linearity
ggplot(faithdata3,aes(residual,faith_pd,color=season))+
  geom_point()

#homogeneity of variance
faithdata3$residualsquare <- faithdata3$residual^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.Model <- lm(residualsquare ~ host_subject_id, data=faithdata3) #ANOVA of the squared residuals
anova(Levene.Model) #p-value is 0.26

Model.F<-lmer(faith_pd ~ ADD*season
              +(1|host_subject_id), #random effects for subject ID and (+) FTND by (|) subject, which are allowed to vary independently from eachother (the 1 and 0 notation)
              data=faithdata3, REML=FALSE)
Plot.Model.F <- plot(Model.F)

#normality
require("lattice")
qqmath(Model.F, id=0.05)

#plot faith's pd
library(ggpubr)
g1 <- subset(faithdata3, host_subject_id == "STAFS2016.064")
B<-ggscatter(faithdata3, x = "ADD", y = "faith_pd",
             add = "reg.line",                         # Add regression line
             conf.int = TRUE,                          # Add confidence interval
             color = "season", palette = "jco",           # Color by groups "cyl"
             shape = "season"                             # Change point shape by groups "cyl"
)+
  xlab("Accumulated Degree Day")+
  ylab("Faith's pd")+
  theme(axis.text.x=element_text(angle = 90, vjust=0.5, hjust = 1))+
  geom_point(data=g1, colour="red", shape="triangle") +
  scale_y_continuous(breaks=c(10,20,30,40,50),limits=c(0,50))

B

ggsave("Faithpd18S.png", height=5, width=5, device="png", dpi=500)

plotmerge<-ggarrange(A,B, nrow=1, common.legend=TRUE, legend="bottom",labels=c("A","B"))

plotmerge

ggsave("FaithpdCombined.png", height=4, width=8, device="png", dpi=500)

### Figure 2A ###

#make a phyloseq object
physeq_16S<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/feature_tables/frag_ins_filtered_noChloMito_17098_table.qza", 
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/02_metadata/maps/map3_R.txt"
)

#taxa table as formatted wouldn't work with my qiime2 library, so reformatted and merged
taxa<-read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/taxonomy/taxonomy_no_prefixes.qza")$data
tax<-taxa %>% separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE,
                       convert = FALSE, extra = "warn", fill = "warn")

#getting rid of blanks and NAs and calling them "Unclassified"
tax[tax==" "]<-NA
tax[is.na(tax)]<-"Unclassified"

#formatting taxonomy table for merging
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)
tax<-tax_table(tax)

#remake physeq to include taxonomy
physeq_16S<-merge_phyloseq(physeq_16S,tax)


#combine at genus level and transform abundance to relative abundance
physeq_16S_genus <-physeq_16S %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept. Works across entire data set
ftax<-physeq_16S_genus %>% filter_taxa(function(x) mean(x) < 0.005, TRUE) %>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_16S_genus<-physeq_16S_genus %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
genus_16S_Avg <- physeq_16S_genus %>%
  group_by(season,ADD_0,sample_simple,Phylum, Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(season,sample_simple,ADD_0) %>%
  unite(Taxa, c(Phylum, Class, Order, Family, Genus), sep = ",", remove = FALSE)

#assign to rare taxa
genus_16S_Avg$Phylum[genus_16S_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
genus_16S_Avg$Taxa[genus_16S_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal phylum levels 
levels(as.factor(genus_16S_Avg$Phylum)) #seven levels 

# use below palette for sequential colors
palette<-c("red","green","yellow","pink","blue","orange", "monochrome")

#to load in a saved palette
palette <- readRDS(file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette.rds")

#use these same colors for the same taxa; build a palette
myColors1 <- palette
names(myColors1)<-levels(as.factor(genus_16S_Avg$Phylum))

#create column Phy_Colors; just adding the palette to my data frame
genus_16S_Avg$Phy_Colors<-genus_16S_Avg$Phylum
genus_16S_Avg$Phy_Colors<-as.factor(genus_16S_Avg$Phy_Colors)
levels(genus_16S_Avg$Phy_Colors)<-c(palette)

###get a small summarised dataframe for color for loop
p<-genus_16S_Avg %>% select_all() %>% group_by(Phylum,Phy_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each phylum and save to object called phylum
phylum<-c()

for(i in 1:nrow(p)) {
  print(p$Phylum[i])
  a<-randomColor(count = p$n_colors[i], hue = p$Phy_Colors[i], luminosity = "bright")
  phylum<-c(phylum,a)
}

#get unique taxa to pair with colors
gen_colors<-genus_16S_Avg %>% select_all() %>% 
  group_by(Phylum,Phy_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named phylum
names(phylum)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
genus_16S_Avg$Taxa<- factor(genus_16S_Avg$Taxa, levels = gen_colors$u_genera)

#plot
ggplot(data=genus_16S_Avg, aes(x = reorder(sample_simple, ADD_0), y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  facet_grid(~season, scales = 'free', space = "free") +
  scale_fill_manual(values= phylum) +
  theme_classic()+
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Accumulated Degree days (ADD)")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=21))

ggsave("RelAbundGenus_luminositybright_final.png", units="in", width = 7.5, height = 5.5,  dpi=500, device="png")  

saveRDS(palette, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette.rds")

### Figure 2B ###

#make initial phyloseq object
physeq_18S_initial<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/03_feature_tables/18S_plate52_ribonly_table_full_taxfiltered_214940.qza",
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/02_metadata/map_18S_R.txt"
)

#get taxonomy files
taxonomy1 = read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/02_taxonomy/18S_plate52_ribonly_full-taxonomy_no_prefixes.qza")$data #this  is the updated taxonomy file using PR2

#separate the taxa by semi-colon according  to PR2 groups; warning is because of extra semi-colon at the end of some rows
#change what is in c() to reflect SILVA rather than PR2
tax<-taxonomy1 %>% separate("Taxon", c("L0", "L1","L2", "L3", "L4", "L5", "L6", "L7","L8","L9","L10","L11","L12","L13","L14"), sep = ";", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn")

#remove confidence and transform to matrix; the 10 will need to be changed
tax<-tax[,-17]
tax[is.na(tax)]<-"Unclassified"

#remove feature id column
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)

#make into tax table
tax2<-tax_table(tax)

#merge to make final phyloseq object
physeq_18S<-merge_phyloseq(physeq_18S_initial,tax2)

#combine at genus level and transform abundance to relative abundance
physeq_18S_L8 <-physeq_18S %>%
  tax_glom(taxrank = "L8") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept and prevalence of 2%
ftax<-physeq_18S_L8 %>% filter_taxa(function(x) mean(x) < 0.005, TRUE) %>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_18S_L8<-physeq_18S_L8 %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
physeq_18S_L8_Avg <- physeq_18S_L8 %>%
  group_by(season,ADD_0,sample_simple,L2,L3,L4,L5,L6,L7,OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(season,sample_simple,ADD_0) %>%
  unite(Taxa, c(L2, L3, L4, L5, L6, L7), sep = ",", remove = FALSE)

#assign to rare taxa
physeq_18S_L8_Avg$L2[physeq_18S_L8_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
physeq_18S_L8_Avg$Taxa[physeq_18S_L8_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal L2 levels 
levels(as.factor(physeq_18S_L8_Avg$L2)) #7 levels 
physeq_18S_L8_Avg$L2<-fct_relevel(physeq_18S_L8_Avg$L2, "RareTaxa", after=6)

# use below palette for sequential colors
palette<-c("red","green","yellow","blue","purple", "monochrome", "monochrome")

#use these same colors for the same taxa; build a palette
myColors1 <- palette
names(myColors1)<-levels(as.factor(physeq_18S_L8_Avg$L2))

#create column Phy_Colors; just adding the palette to my data frame
physeq_18S_L8_Avg$L2_Colors<-physeq_18S_L8_Avg$L2
physeq_18S_L8_Avg$L2_Colors<-as.factor(physeq_18S_L8_Avg$L2_Colors)
levels(physeq_18S_L8_Avg$L2_Colors)<-c(palette)

###get a small summarised dataframe for color for loop
p<-physeq_18S_L8_Avg %>% select_all() %>% group_by(L2,L2_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each L2 and save to object called L2
level2 <- c()

#if more random colors rather than sequential
for(i in 1:nrow(p)) {
  print(p$L2[i])
  a<-randomColor(count = p$n_colors[i], hue = p$L2_Colors[i], luminosity = "bright")
  level2<-c(level2,a)
}

#get unique taxa to pair with colors
gen_colors<-physeq_18S_L8_Avg %>% select_all() %>% 
  group_by(L2,L2_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named level2
names(level2)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
physeq_18S_L8_Avg$Taxa<- factor(physeq_18S_L8_Avg$Taxa, levels = gen_colors$u_genera)

#plot
ggplot(data=physeq_18S_L8_Avg, aes(x = reorder(sample_simple, ADD_0), y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  facet_grid(~season, scales = 'free', space = "free") +
  scale_fill_manual(values= level2) +
  theme_classic()+
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Accumulated Degree days (ADD)")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=21))

ggsave("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/04_results_other/taxonomy/RelAbund_luminositybright_final.png", units="in", width = 7.5, height = 5.5,  dpi=500, device="png")  

saveRDS(palette, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/04_results_other/taxonomy/palette.rds")

### Figure 3A and 3B ###
# spring alone
# import metadata and weighted unifrac qza files
metadata <- read_csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_skin_soil_spring_summer_allQIITAIDs_R.csv")
metadata

pcoa_weighted_spring <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-STsamplesonly-spring-4548/weighted_unifrac_pcoa_results.qza")

# just looking at stuff 
pcoa_weighted_spring$uuid
head(pcoa_weighted_spring$data$ProportionExplained)
pcoa_weighted_spring$data$Vectors[1:5, 1:3]

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

### plot weighted unifrac
pcoa_weighted_spring$data$Vectors %>%
  #rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color = sample_type_timepoint)) +
  geom_point(size = 2.5) +
  xlab(paste("PC1: ", round(100*pcoa_weighted_spring$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_weighted_spring$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"),
        #axis.title=element_text(size=20, face="bold", color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        plot.title=element_text(size=20, color = "black"),
        title = black.bold.text, axis.title = black.bold.text,
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.text=element_text(size=20, color = "black"),
        legend.title =element_text(size=20, color = "black"),
        legend.background = element_rect(fill = "white")) +
  scale_color_manual(name = "Sample Type",
                     labels = c("Rib", "Fresh Skin", "Advanced Skin", "Fresh Soil", "Advanced Soil"),
                     values=c("red","deepskyblue1", "deepskyblue4", "gold", "lightsalmon4")) +
  labs(title = "Spring")

# summer alone
# import metadata and weighted unifrac qza files
metadata <- read_csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_skin_soil_spring_summer_allQIITAIDs_R.csv")
metadata

pcoa_weighted_summer <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-STsamplesonly-summer-4548/weighted_unifrac_pcoa_results.qza")

# just looking at stuff 
pcoa_weighted_summer$uuid
head(pcoa_weighted_summer$data$ProportionExplained)
pcoa_weighted_summer$data$Vectors[1:5, 1:3]

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

### plot weighted unifrac
pcoa_weighted_summer$data$Vectors %>%
  #rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color = sample_type_timepoint)) +
  geom_point(size = 2.5) +
  xlab(paste("PC1: ", round(100*pcoa_weighted_summer$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_weighted_summer$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"),
        #axis.title=element_text(size=20, face="bold", color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        plot.title=element_text(size=20, color = "black"),
        title = black.bold.text, axis.title = black.bold.text,
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.text=element_text(size=20, color = "black"),
        legend.title =element_text(size=20, color = "black"),
        legend.background = element_rect(fill = "white")) +
  scale_color_manual(name = "Sample Type",
                     labels = c("Rib", "Fresh Skin", "Advanced Skin", "Fresh Soil", "Advanced Soil"),
                     values=c("red","deepskyblue1", "deepskyblue4", "gold", "lightsalmon4")) +
  labs(title = "Summer")

### Figure 4 ###
# both seasons
# import metadata and unweighted unifrac qza files
metadata <- read_tsv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/02_metadata/maps/map3.txt")
metadata

pcoa <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/core-metrics/core-metrics-results-frag-ins/unweighted_unifrac_pcoa_results.qza")

# just looking at stuff 
pcoa$uuid
head(pcoa$data$ProportionExplained)
pcoa$data$Vectors[1:5, 1:3]

### plot beta diversity
# unweighted unifrac
pcoa$data$Vectors %>%
  rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=season, size=ADD)) +
  geom_point() +
  xlab(paste("PC1: ", round(100*pcoa$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        legend.title =element_text(size=20),
        plot.title=element_text(size=20)) +
  scale_color_manual(values=c("deepskyblue","tomato")) +
  ggtitle("Spring and Summer")

### spring

# import qza
pcoa_spring <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/core-metrics/core-metrics-frag-ins-spring/unweighted_unifrac_pcoa_results.qza")

# plot unweighted unifrac pcoa
pcoa_spring$data$Vectors %>%
  rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, size=ADD)) +
  geom_point(color = "deepskyblue") +
  xlab(paste("PC1: ", round(100*pcoa_spring$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_spring$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        legend.title =element_text(size=20),
        plot.title=element_text(size=20)) +
  ggtitle("Spring")

### summer

# import qza
pcoa_summer <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/core-metrics/core-metrics-frag-ins-summer/unweighted_unifrac_pcoa_results.qza")

# plot unweighted unifrac pcoa
pcoa_summer$data$Vectors %>%
  rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, size=ADD)) +
  geom_point(color = "tomato") +
  xlab(paste("PC1: ", round(100*pcoa_summer$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_summer$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        legend.title =element_text(size=20),
        plot.title=element_text(size=20)) +
  ggtitle("Summer")

### Figure S2A ###

#make a phyloseq object
physeq_16S_neg<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/feature_tables/frag_ins_filtered_noChloMito_table_negextractions_only.qza", 
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/02_metadata/maps/map3_R_negs.txt"
)

#taxa table as formatted wouldn't work with my qiime2 library, so reformatted and merged
taxa<-read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/taxonomy/taxonomy_no_prefixes.qza")$data
tax<-taxa %>% separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE,
                       convert = FALSE, extra = "warn", fill = "warn")

#getting rid of blanks and NAs and calling them "Unclassified"
tax[tax==" "]<-NA
tax[is.na(tax)]<-"Unclassified"

#formatting taxonomy table for merging
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)
tax<-tax_table(tax)

#remake physeq to include taxonomy
physeq_16S_neg<-merge_phyloseq(physeq_16S_neg,tax)

#combine at genus level and transform abundance to relative abundance
physeq_16S_neg_genus <-physeq_16S_neg %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept. Works across entire data set
ftax<-filter_taxa(physeq_16S_neg_genus, function(x) mean(x) < 0.005, TRUE) %>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_16S_neg_genus<-physeq_16S_neg_genus %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
genus_16S_neg_Avg <- physeq_16S_neg_genus %>%
  group_by(Sample,Phylum, Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(Sample) %>%
  unite(Taxa, c(Phylum, Class, Order, Family, Genus), sep = ",", remove = FALSE)

#assign to rare taxa 
genus_16S_neg_Avg$Phylum[genus_16S_neg_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
genus_16S_neg_Avg$Taxa[genus_16S_neg_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal phylum levels 
levels(as.factor(genus_16S_neg_Avg$Phylum)) #five levels 

# use below palette for sequential colors
palette_neg<-c("red","green","yellow","blue","monochrome")

#to load in a saved palette
palette_neg <- readRDS(file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette_neg.rds")

#use these same colors for the same taxa; build a palette
myColors_neg <- palette_neg
names(myColors_neg)<-levels(as.factor(genus_16S_neg_Avg$Phylum))

#create column Phy_Colors; just adding the palette to my data frame
genus_16S_neg_Avg$Neg_Colors<-genus_16S_neg_Avg$Phylum
genus_16S_neg_Avg$Neg_Colors<-as.factor(genus_16S_neg_Avg$Neg_Colors)
levels(genus_16S_neg_Avg$Neg_Colors)<-c(palette_neg)

###get a small summarised dataframe for color for loop
p<-genus_16S_neg_Avg %>% select_all() %>% group_by(Phylum,Neg_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each phylum and save to object called phylum
phylum<-c()

#if more random colors rather than sequential
library(randomcoloR) 
for(i in 1:nrow(p)) {
  print(p$Phylum[i])
  a<-randomColor(count = p$n_colors[i], hue = p$Neg_Colors[i], luminosity = "bright")
  phylum<-c(phylum,a)
}

#get unique taxa to pair with colors
gen_colors<-genus_16S_neg_Avg %>% select_all() %>% 
  group_by(Phylum,Neg_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named phylum
names(phylum)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
genus_16S_neg_Avg$Taxa<- factor(genus_16S_neg_Avg$Taxa, levels = gen_colors$u_genera)

#plot
ggplot(data=genus_16S_neg_Avg, aes(x = Sample, y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  scale_fill_manual(values=phylum) +
  theme_classic()+
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Negative Extraction Samples")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=8))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 1))

ggsave("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/RelAbundGenus_negs_luminositybright_final.png", units="in", width = 10, height = 5.5,  dpi=500, device="png")  

saveRDS(palette_neg, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette_neg.rds")

### Figure S2B ###

#make initial phyloseq object
physeq_18S_neg_initial<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/03_feature_tables/18S_plate52_ribonly_table_full_taxfiltered_negextractions_only.qza",
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/02_metadata/map_18S_negs_R.txt"
)

#get taxonomy files
taxonomy1 = read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/02_taxonomy/18S_plate52_ribonly_full-taxonomy_no_prefixes.qza")$data #this  is the updated taxonomy file using PR2

#separate the taxa by semi-colon according  to PR2 groups; warning is because of extra semi-colon at the end of some rows
#change what is in c() to reflect SILVA rather than PR2
tax<-taxonomy1 %>% separate("Taxon", c("L0", "L1","L2", "L3", "L4", "L5", "L6", "L7","L8","L9","L10","L11","L12","L13","L14"), sep = ";", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn")

#remove confidence and transform to matrix; the 10 will need to be changed
tax<-tax[,-17]
tax[is.na(tax)]<-"Unclassified"

#remove feature id column
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)

#make into tax table
tax2<-tax_table(tax)

#merge to make final phyloseq object
physeq_18S_neg<-merge_phyloseq(physeq_18S_neg_initial,tax2)

#combine at genus level and transform abundance to relative abundance
physeq_18S_neg_L8 <-physeq_18S_neg %>%
  tax_glom(taxrank = "L8") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept and prevalence of 2%
ftax<-physeq_18S_neg_L8 %>% filter_taxa(function(x) mean(x) < 0.005, TRUE) %>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_18S_neg_L8<-physeq_18S_neg_L8 %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
physeq_18S_neg_L8_Avg <- physeq_18S_neg_L8 %>%
  group_by(Sample,L2,L3,L4,L5,L6,L7,OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(Sample) %>%
  unite(Taxa, c(L2, L3, L4, L5, L6, L7), sep = ",", remove = FALSE)

#assign to rare taxa alternative
physeq_18S_neg_L8_Avg$L2[physeq_18S_neg_L8_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
physeq_18S_neg_L8_Avg$Taxa[physeq_18S_neg_L8_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal L2 levels 
levels(as.factor(physeq_18S_neg_L8_Avg$L2)) #5 levels 
physeq_18S_neg_L8_Avg$L2<-fct_relevel(physeq_18S_neg_L8_Avg$L2, "RareTaxa", after=4)

# use below palette for sequential colors
palette_neg<-c("red","yellow","blue","monochrome", "monochrome")

#to load in a saved palette
palette_neg <- readRDS(file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/04_results_other/taxonomy/palette_neg.rds")

#use these same colors for the same taxa; build a palette
myColors_neg <- palette_neg
names(myColors_neg)<-levels(as.factor(physeq_18S_neg_L8_Avg$L2))

#create column L2_Colors; just adding the palette to my data frame
physeq_18S_neg_L8_Avg$L2_Colors<-physeq_18S_neg_L8_Avg$L2
physeq_18S_neg_L8_Avg$L2_Colors<-as.factor(physeq_18S_neg_L8_Avg$L2_Colors)
levels(physeq_18S_neg_L8_Avg$L2_Colors)<-c(palette_neg)

###get a small summarised dataframe for color for loop
p<-physeq_18S_neg_L8_Avg %>% select_all() %>% group_by(L2,L2_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each L2 and save to object called L2
level2 <- c()

#if more random colors rather than sequential
for(i in 1:nrow(p)) {
  print(p$L2[i])
  a<-randomColor(count = p$n_colors[i], hue = p$L2_Colors[i], luminosity = "bright")
  level2<-c(level2,a)
}

#get unique taxa to pair with colors
gen_colors<-physeq_18S_neg_L8_Avg %>% select_all() %>% 
  group_by(L2,L2_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named level2
names(level2)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
physeq_18S_neg_L8_Avg$Taxa<- factor(physeq_18S_neg_L8_Avg$Taxa, levels = gen_colors$u_genera)

#plot
ggplot(data=physeq_18S_neg_L8_Avg, aes(x = Sample, y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  scale_fill_manual(values= level2) +
  theme_classic()+
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Negative Extraction Samples")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=8))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 1))

ggsave("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/04_results_other/taxonomy/RelAbund_negs_luminositybright_final.png", units="in", width = 10, height = 5.5,  dpi=500, device="png")  

saveRDS(palette_neg, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/04_results_other/taxonomy/palette_neg.rds")

### Figure S3 ###

# spring

ST_shannon_spring <- read.csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/shannon_R_spring.csv")

# do a boxplot with pairwise comparison
# pairwise test is Wilcoxon (which is also nonparametric)
# the wilcoxon test p-values are corrected for multiple comparisons
# adjusted p-value is shown

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

# statistical test
stat.test_spring <- ST_shannon_spring %>%
  wilcox_test(shannon_entropy ~ sample_type_timepoint) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_xy_position(x = "supp") %>% 
  add_significance("p.adj")

stat.test_spring

ggboxplot(ST_shannon_spring, x = "sample_type_timepoint", y = "shannon_entropy", color = "sample_type_timepoint") +
  stat_pvalue_manual(stat.test_spring, label = "{p.adj}{p.adj.signif}", tip.length = 0, hide.ns = TRUE) +
  scale_x_discrete(labels = c("rib" = "Rib", 
                              "skin_early" = "Fresh Skin",
                              "skin_late" = "Advanced Skin",
                              "soil_early" = "Fresh Soil",
                              "soil_late" = "Advanced Soil")) +
  scale_color_manual(values=c("red","deepskyblue1", "deepskyblue4", "gold", "lightsalmon4")) +
  labs(x="Sample Type", y="Shannon", title = "Spring") +
  theme(legend.position = "none",
        title = black.bold.text, axis.title = black.bold.text,
        axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"))

# summer

ST_shannon_summer <- read.csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/shannon_R_summer.csv")

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

# statistical test
stat.test_summer <- ST_shannon_summer %>%
  wilcox_test(shannon_entropy ~ sample_type_timepoint) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_xy_position(x = "supp") %>% 
  add_significance("p.adj")

stat.test_summer

ggboxplot(ST_shannon_summer, x = "sample_type_timepoint", y = "shannon_entropy", color = "sample_type_timepoint") +
  stat_pvalue_manual(stat.test_summer, label = "{p.adj}{p.adj.signif}", tip.length = 0, hide.ns = TRUE) +
  scale_x_discrete(labels = c("rib" = "Rib", 
                              "skin_early" = "Fresh Skin",
                              "skin_late" = "Advanced Skin",
                              "soil_early" = "Fresh Soil",
                              "soil_late" = "Advanced Soil")) +
  scale_color_manual(values=c("red","deepskyblue1", "deepskyblue4", "gold", "lightsalmon4")) +
  labs(x="Sample Type", y="Shannon", title = "Summer") +
  theme(legend.position = "none",
        title = black.bold.text, axis.title = black.bold.text,
        axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"))

### Figure S4A ###

# spring, skin and soil
#make a phyloseq object
physeq_ST_spring<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-STsamplesonly-spring-4548/rarefied_table.qza", 
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_skin_soil_spring_nosoilctrl_noribs_allQIITAIDs_R.txt"
)

#taxa table as formatted wouldn't work with my qiime2 library, so reformatted and merged
taxa<-read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/taxonomy/bone_SFWall_spring_taxonomy_no_prefixes.qza")$data
tax<-taxa %>% separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE,
                       convert = FALSE, extra = "warn", fill = "warn")

#getting rid of blanks and NAs and calling them "Unclassified"
tax[tax==" "]<-NA
tax[is.na(tax)]<-"Unclassified"

#formatting taxonomy table for merging
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)
tax<-tax_table(tax)

#remake physeq to include taxonomy
physeq_ST_spring<-merge_phyloseq(physeq_ST_spring,tax)

#combine at genus level and transform abundance to relative abundance
physeq_ST_spring_genus <-physeq_ST_spring %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept. Works across entire data set
ftax<-physeq_ST_spring_genus %>% filter_taxa(function(x) mean(x) < 0.002, TRUE) #%>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_ST_spring_genus<-physeq_ST_spring_genus %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
genus_ST_spring_Avg <- physeq_ST_spring_genus %>%
  group_by(day_of_collection_numeric, sample_type_timepoint, sample_simple, Phylum, Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(sample_type_timepoint, sample_simple, day_of_collection_numeric) %>%
  unite(Taxa, c(Phylum, Class, Order, Family, Genus), sep = ",", remove = FALSE)

#assign to rare taxa
genus_ST_spring_Avg$Phylum[genus_ST_spring_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
genus_ST_spring_Avg$Taxa[genus_ST_spring_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal phylum levels 
levels(as.factor(genus_ST_spring_Avg$Phylum)) #seven levels 

# use below palette for sequential colors
palette<-c("purple","red", "green","yellow","pink","orange", "blue","monochrome")

#to load in a saved palette
palette <- readRDS(file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/palette.rds")

#use these same colors for the same taxa; build a palette
myColors1 <- palette
names(myColors1)<-levels(as.factor(genus_ST_spring_Avg$Phylum))

#create column Phy_Colors; just adding the palette to my data frame
genus_ST_spring_Avg$Phy_Colors<-genus_ST_spring_Avg$Phylum
genus_ST_spring_Avg$Phy_Colors<-as.factor(genus_ST_spring_Avg$Phy_Colors)
levels(genus_ST_spring_Avg$Phy_Colors)<-c(palette)

###get a small summarised dataframe for color for loop
p<-genus_ST_spring_Avg %>% select_all() %>% group_by(Phylum,Phy_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each phylum and save to object called phylum
phylum<-c()

#if more random colors rather than sequential
library(randomcoloR) 
for(i in 1:nrow(p)) {
  print(p$Phylum[i])
  a<-randomColor(count = p$n_colors[i], hue = p$Phy_Colors[i], luminosity = "bright")
  phylum<-c(phylum,a)
}

#get unique taxa to pair with colors
gen_colors<-genus_ST_spring_Avg %>% select_all() %>% 
  group_by(Phylum,Phy_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named phylum
names(phylum)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
genus_ST_spring_Avg$Taxa<- factor(genus_ST_spring_Avg$Taxa, levels = gen_colors$u_genera)

# change facet grid labels
facet_labels <- list('skin_early'="Fresh Skin",
                     'skin_late'="Advanced Skin",
                     'soil_early'="Fresh Soil",
                     'soil_late'="Advanced Soil")

facet_labeller <- function(variable,value){
  return(facet_labels[value])}

#plot
ggplot(data=genus_ST_spring_Avg, aes(x = reorder(sample_simple, day_of_collection_numeric), y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  facet_grid(~sample_type_timepoint, scales="free", space="free", labeller = facet_labeller) +
  scale_fill_manual(values= phylum) +
  theme_classic()+
  #scale_y_continuous(limits=c(0,1.25))+
  # Remove x axis title
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Source Samples")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=10))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=21))

ggsave("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/RelAbundGenus_sources_spring_luminositybright_final.png", units="in", width = 25, height = 10,  dpi=500, device="png")  

saveRDS(palette, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette.rds")

### Figure S4B ###

# summer, skin and soil
#make a phyloseq object
physeq_ST_summer<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-STsamplesonly-summer-4548/rarefied_table.qza", 
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_skin_soil_summer_nosoilctrl_noribs_allQIITAIDs_R.txt"
)

#taxa table as formatted wouldn't work with my qiime2 library, so reformatted and merged
taxa<-read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/taxonomy/bone_SFWall_spring_taxonomy_no_prefixes.qza")$data
tax<-taxa %>% separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE,
                       convert = FALSE, extra = "warn", fill = "warn")

#getting rid of blanks and NAs and calling them "Unclassified"
tax[tax==" "]<-NA
tax[is.na(tax)]<-"Unclassified"

#formatting taxonomy table for merging
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)
tax<-tax_table(tax)

#remake physeq to include taxonomy
physeq_ST_summer<-merge_phyloseq(physeq_ST_summer,tax)

#combine at genus level and transform abundance to relative abundance
physeq_ST_summer_genus <-physeq_ST_summer %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept. Works across entire data set
ftax<-physeq_ST_summer_genus %>% filter_taxa(function(x) mean(x) < 0.003, TRUE) #%>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_ST_summer_genus<-physeq_ST_summer_genus %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
genus_ST_summer_Avg <- physeq_ST_summer_genus %>%
  group_by(day_of_collection_numeric, sample_type_timepoint, sample_simple, Phylum, Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(sample_type_timepoint, sample_simple, day_of_collection_numeric) %>%
  unite(Taxa, c(Phylum, Class, Order, Family, Genus), sep = ",", remove = FALSE)

#assign to rare taxa alternative
genus_ST_summer_Avg$Phylum[genus_ST_summer_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
genus_ST_summer_Avg$Taxa[genus_ST_summer_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal phylum levels 
levels(as.factor(genus_ST_summer_Avg$Phylum)) #seven levels 

# use below palette for sequential colors
palette<-c("purple","red", "green","yellow","orange", "blue","monochrome")

#to load in a saved palette
palette <- readRDS(file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/palette.rds")

#use these same colors for the same taxa; build a palette
myColors1 <- palette
names(myColors1)<-levels(as.factor(genus_ST_summer_Avg$Phylum))

#create column Phy_Colors; just adding the palette to my data frame
genus_ST_summer_Avg$Phy_Colors<-genus_ST_summer_Avg$Phylum
genus_ST_summer_Avg$Phy_Colors<-as.factor(genus_ST_summer_Avg$Phy_Colors)
levels(genus_ST_summer_Avg$Phy_Colors)<-c(palette)

###get a small summarised dataframe for color for loop
p<-genus_ST_summer_Avg %>% select_all() %>% group_by(Phylum,Phy_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each phylum and save to object called phylum
phylum<-c()

#if more random colors rather than sequential
library(randomcoloR) 
for(i in 1:nrow(p)) {
  print(p$Phylum[i])
  a<-randomColor(count = p$n_colors[i], hue = p$Phy_Colors[i], luminosity = "bright")
  phylum<-c(phylum,a)
}

#get unique taxa to pair with colors
gen_colors<-genus_ST_summer_Avg %>% select_all() %>% 
  group_by(Phylum,Phy_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named phylum
names(phylum)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
genus_ST_summer_Avg$Taxa<- factor(genus_ST_summer_Avg$Taxa, levels = gen_colors$u_genera)

# change facet grid labels
facet_labels <- list('skin_early'="Fresh Skin",
                     'skin_late'="Advanced Skin",
                     'soil_early'="Fresh Soil",
                     'soil_late'="Advanced Soil")

facet_labeller <- function(variable,value){
  return(facet_labels[value])}

#plot
ggplot(data=genus_ST_summer_Avg, aes(x = reorder(sample_simple, day_of_collection_numeric), y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  facet_grid(~sample_type_timepoint, scales="free", space="free", labeller = facet_labeller) +
  scale_fill_manual(values= phylum) +
  theme_classic()+
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Source Samples")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=10))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=21))

ggsave("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/RelAbundGenus_sources_summer_luminositybright_final.png", units="in", width = 25, height = 10,  dpi=500, device="png")  

saveRDS(palette, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette.rds")

### Figure S5 ###
# spring alone
# import metadata and weighted unifrac qza files
metadata <- read_csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_skin_soil_spring_summer_allQIITAIDs_R.csv")
metadata

pcoa_unweighted_spring <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-STsamplesonly-spring-4548/unweighted_unifrac_pcoa_results.qza")

# just looking at stuff 
pcoa_unweighted_spring$uuid
head(pcoa_unweighted_spring$data$ProportionExplained)
pcoa_unweighted_spring$data$Vectors[1:5, 1:3]

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

### plot unweighted unifrac
pcoa_unweighted_spring$data$Vectors %>%
  #rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color = sample_type_timepoint)) +
  geom_point(size = 2.5) +
  xlab(paste("PC1: ", round(100*pcoa_unweighted_spring$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_unweighted_spring$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"),
        #axis.title=element_text(size=20, face="bold", color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        plot.title=element_text(size=20, color = "black"),
        title = black.bold.text, axis.title = black.bold.text,
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.text=element_text(size=20, color = "black"),
        legend.title =element_text(size=20, color = "black"),
        legend.background = element_rect(fill = "white")) +
  scale_color_manual(name = "Sample Type",
                     labels = c("Rib", "Fresh Skin", "Advanced Skin", "Fresh Soil", "Advanced Soil"),
                     values=c("red","deepskyblue1", "deepskyblue4", "gold", "lightsalmon4")) +
  labs(title = "Spring")

# summer alone
# import metadata and weighted unifrac qza files
metadata <- read_csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_skin_soil_spring_summer_allQIITAIDs_R.csv")
metadata

pcoa_unweighted_summer <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-STsamplesonly-summer-4548/unweighted_unifrac_pcoa_results.qza")

# just looking at stuff 
pcoa_unweighted_summer$uuid
head(pcoa_unweighted_summer$data$ProportionExplained)
pcoa_unweighted_summer$data$Vectors[1:5, 1:3]

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

### plot weighted unifrac
pcoa_unweighted_summer$data$Vectors %>%
  #rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color = sample_type_timepoint)) +
  geom_point(size = 2.5) +
  xlab(paste("PC1: ", round(100*pcoa_unweighted_summer$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_unweighted_summer$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"),
        #axis.title=element_text(size=20, face="bold", color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        plot.title=element_text(size=20, color = "black"),
        title = black.bold.text, axis.title = black.bold.text,
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.text=element_text(size=20, color = "black"),
        legend.title =element_text(size=20, color = "black"),
        legend.background = element_rect(fill = "white")) +
  scale_color_manual(name = "Sample Type",
                     labels = c("Rib", "Early Skin", "Decomposer Skin", "Early Soil", "Decomposer Soil"),
                     values=c("red","deepskyblue1", "deepskyblue4", "gold", "lightsalmon4")) +
  labs(title = "Summer")

### Figure S6 ###

# import metadata and unweighted unifrac qza files
metadata <- read_tsv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/02_metadata/map_18S.txt")
metadata

pcoa <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/02_18S/Hiseq_original_run/forward_reads_only/01_qiime2_analysis/05_core_metrics/core-metrics-results/unweighted_unifrac_pcoa_results.qza")

# just looking at stuff 
pcoa$uuid
head(pcoa$data$ProportionExplained)
pcoa$data$Vectors[1:5, 1:3]

### plot beta diversity
# unweighted unifrac
pcoa$data$Vectors %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=season, size=ADD)) +
  geom_point() +
  xlab(paste("PC1: ", round(100*pcoa$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"),
        axis.title=element_text(size=20, face="bold", color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        plot.title=element_text(size=20, color = "black"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.text=element_text(size=20, color = "black"),
        legend.title =element_text(size=20, color = "black"),
        legend.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("deepskyblue","tomato")) +
  ggtitle("18S rRNA Unweighted UniFrac PCoA")

### Figure S8B (Figure S8A made in QIIME2 view and powerpoint) ###
#make a phyloseq object
physeq_ST_ctrl<-qza_to_phyloseq(
  features="/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/core_metrics/core-metrics-soilctrlonly-4548/rarefied_table.qza", 
  metadata = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/mapping_files/map_bone_SHSU_soil_ctrls_only_allQIITAIDs_R.txt"
)

#taxa table as formatted wouldn't work with my qiime2 library, so reformatted and merged
taxa<-read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/qiime2_analysis/taxonomy/bone_SFWall_spring_taxonomy_no_prefixes.qza")$data
tax<-taxa %>% separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE,
                       convert = FALSE, extra = "warn", fill = "warn")

#getting rid of blanks and NAs and calling them "Unclassified"
tax[tax==" "]<-NA
tax[is.na(tax)]<-"Unclassified"

#formatting taxonomy table for merging
rownames(tax)<-tax$Feature.ID
tax<-tax[,-1]
tax<-as.matrix(tax)
tax<-tax_table(tax)

#remake physeq to include taxonomy
physeq_ST_ctrl<-merge_phyloseq(physeq_ST_ctrl,tax)

#combine at genus level and transform abundance to relative abundance
physeq_ST_ctrl_genus <-physeq_ST_ctrl %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at level of interest
  phyloseq::transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance

#can use this for filtering rare taxa;  genera with a mean greater than 0.005 are kept. Works across entire data set
ftax<-physeq_ST_ctrl_genus %>% filter_taxa(function(x) mean(x) < 0.008, TRUE) #%>% filter_taxa(function(x) sum(x > 0) > (0.02*length(x)), TRUE)
#get names of taxa with mean less than 0.005
ftax.names<-taxa_names(ftax)

#melt the physeq table
physeq_ST_ctrl_genus<-physeq_ST_ctrl_genus %>% psmelt()

#Average relative abundance by category; obviously you don't have to get a mean relative abundance 
genus_ST_ctrl_Avg <- physeq_ST_ctrl_genus %>%
  group_by(season,day_of_collection_numeric,sample_simple, Phylum, Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::arrange(season,sample_simple,day_of_collection_numeric) %>%
  unite(Taxa, c(Phylum, Class, Order, Family, Genus), sep = ",", remove = FALSE)

#assign to rare taxa 
genus_ST_ctrl_Avg$Phylum[genus_ST_ctrl_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")
genus_ST_ctrl_Avg$Taxa[genus_ST_ctrl_Avg$OTU %in% ftax.names] <- as.character("RareTaxa")

#make palette; need colors to equal phylum levels 
levels(as.factor(genus_ST_ctrl_Avg$Phylum)) #seven levels 

# use below palette for sequential colors
palette<-c("purple","red", "green", "orange", "blue", "pink", "yellow", "monochrome")

#to load in a saved palette
palette <- readRDS(file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/palette.rds")

#use these same colors for the same taxa; build a palette
myColors1 <- palette
names(myColors1)<-levels(as.factor(genus_ST_ctrl_Avg$Phylum))

#create column Phy_Colors; just adding the palette to my data frame
genus_ST_ctrl_Avg$Phy_Colors<-genus_ST_ctrl_Avg$Phylum
genus_ST_ctrl_Avg$Phy_Colors<-as.factor(genus_ST_ctrl_Avg$Phy_Colors)
levels(genus_ST_ctrl_Avg$Phy_Colors)<-c(palette)

###get a small summarised dataframe for color for loop
p<-genus_ST_ctrl_Avg %>% select_all() %>% group_by(Phylum,Phy_Colors)  %>%
  summarize(n_colors=n_distinct(Taxa)) %>%
  mutate_if(is.factor, as.character)

#gets unique colors for each phylum and save to object called phylum
phylum<-c()

#if more random colors rather than sequential
library(randomcoloR) 
for(i in 1:nrow(p)) {
  print(p$Phylum[i])
  a<-randomColor(count = p$n_colors[i], hue = p$Phy_Colors[i], luminosity = "bright")
  phylum<-c(phylum,a)
}

#get unique taxa to pair with colors
gen_colors<-genus_ST_ctrl_Avg %>% select_all() %>% 
  group_by(Phylum,Phy_Colors) %>% arrange(avg.rel.abund) %>%
  summarize(u_genera=unique(Taxa))

#add unique taxa to our newly made palette named phylum
names(phylum)<-gen_colors$u_genera

#need to reorder factor levels for the Taxa column so we get a nice gradient
genus_ST_ctrl_Avg$Taxa<- factor(genus_ST_ctrl_Avg$Taxa, levels = gen_colors$u_genera)

#plot
ggplot(data=genus_ST_ctrl_Avg, aes(x = reorder(sample_simple, day_of_collection_numeric), y = avg.rel.abund)) + 
  geom_bar(aes(fill=Taxa),stat = "identity",color="black") +
  facet_grid(~season, scales = 'free', space = "free") +
  scale_fill_manual(values= phylum) +
  theme_classic()+
  theme(axis.title.x = element_text(size=8)) + 
  xlab("Soil Control Samples")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, nrow=13)) +
  ylab("Average Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        strip.text.x = element_text(face="bold", vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_blank(), legend.text =element_text(size=10))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=21))

ggsave("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/06_sourcetracker/ST_graphing_files/RelAbundGenus_soilctrl_luminositybright_final.png", units="in", width = 15, height = 10,  dpi=500, device="png")  

saveRDS(palette, file = "/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/04_results_other/taxa_plots/palette.rds")

### Figure S9 ###

#import metadata
metadata <- read.csv("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/02_metadata/maps/map3_R.csv", header = TRUE)

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

### precipitation graph with both seasons
precip_plot <- ggplot(data = metadata, aes(x = ADD_0, y = accumulated_precipitation_inches_base0, group = season)) +
  geom_line(aes(linetype=season)) +
  geom_point() +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5, size = 12), 
        axis.text.y = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        title = black.bold.text, axis.title = black.bold.text,
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title.align = 0, 
        plot.background = element_rect(fill='white'),
        legend.background = element_rect(fill='white'),
        legend.text = element_text(color = "black", size = 17),
        legend.key.size = unit(2,"line"),
        legend.title = element_text(color = "black", size = 18)) +
  labs(x = "Accumulated Degree Day", y = "Accumulated Precipitation (inches)")

### humidity graph with both seasons
hum_plot <- ggplot(data = metadata, aes(x = ADD_0, y = accumulated_percent_humidity_days_base0, group = season)) +
  geom_line(aes(linetype=season)) +
  geom_point() +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5, size = 12), 
        axis.text.y = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        title = black.bold.text, axis.title = black.bold.text,
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title.align = 0, 
        plot.background = element_rect(fill='white'),
        legend.background = element_rect(fill='white'),
        legend.text = element_text(color = "black", size = 17),
        legend.key.size = unit(2,"line"),
        legend.title = element_text(color = "black", size = 18)) +
  labs(x = "Accumulated degree day", y = "Accumulated Percent Humidity")

### combine the two graphs

precip_hum_spring_summer <- ggarrange(precip_plot, hum_plot,
                                      labels = c("A", "B"),
                                      ncol = 2)
precip_hum_spring_summer

### Figure S10 ###
# import qza
pcoa_host <- read_qza("/Users/heatherdeel/Dropbox/PMI_3_analyses/bone/01_16S/01_qiime2_analysis/core-metrics/core-metrics-results-frag-ins/weighted_unifrac_pcoa_results.qza")

black.bold.text = element_text(face = "bold", color = "black")
black.text = element_text(color = "black")

### plot weighted unifrac
pcoa_host$data$Vectors %>%
  rename("#SampleID"=SampleID) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=host_subject_id)) +
  geom_point(size = 3.5) +
  xlab(paste("PC1: ", round(100*pcoa_host$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*pcoa_host$data$ProportionExplained[2]), "%")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        legend.title =element_text(size=20),
        plot.title=element_text(size=20)) +
  scale_color_manual(name = "Host Number",
                     labels = c("007", "011", "024", "064", "065", "067"),
                     values = c("red","orange1","gold","green","blue","purple")) +
  ggtitle("16S rRNA Rib Communities by Host")
