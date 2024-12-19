setwd("M:/ORDIAmur_Phase_II/WP3-Experiment/WP3_ASV-Level-analysis/")

#load packages
library(ggplot2)
library(forcats)
library(colorspace)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

relabu_ASVs <- read.table("WP3_selected_genera_ASVs.txt", header=T, sep="\t", dec=".")
View(relabu_ASVs)

#subset Nectria
WP3_Nectria <- subset(relabu_ASVs, Genus =="Nectria")
WP3_Nectria

Nectria <-ggplot(WP3_Nectria, aes(fill=ASV, y=RA, x=genotype)) + 
  geom_bar(position="stack", stat="identity", color="grey")+
  facet_grid(~fct_relevel(soil, 'HG', 'EH', 'RU')) +
  #scale_fill_manual()+
  labs(x= "", y= "RA [%]", title= "Nectria") +
  theme_bw(base_size = 16) +
  scale_x_discrete(guide = guide_axis(angle = 45))
Nectria
#ggsave("WP3_Nectria_ASVs_stackedbarchart.png", width = 9, height = 6, units = "in")


#subset Fusarium
WP3_Fusarium <- subset(relabu_ASVs, Genus =="Fusarium")
WP3_Fusarium

Fusarium <-ggplot(WP3_Fusarium, aes(fill=ASV, y=RA, x=genotype)) + 
  geom_bar(position="stack", stat="identity", color="grey")+
  facet_grid(~fct_relevel(soil, 'HG', 'EH', 'RU')) +
  #scale_fill_manual()+
  theme_bw(base_size = 16) +
  labs(x= "", y= "RA [%]", title= "Fusarium") +
  scale_x_discrete(guide = guide_axis(angle = 45))
Fusarium
#ggsave("WP3_Fusarium_ASVs_stackedbarchart.png", width = 12, height = 6, units = "in")

ggarrange(Nectria,Fusarium, ncol=1, nrow=2, common.legend=FALSE)     # Ohne Legende unter allen Plots
ggsave("WP3_stacked_barchartsNectria_Fusarium_ASV.png", width = 14, height = 10, units = "in")
