# differentially abundant taxa were identified using ANCOM-BC2. Post-analysis, taxa were filtered and only taxa that passed sensitivity analysis, appeared in n>4 replicates and that had an average relative abundance >0.25%  were filtered.
# for display, only four bacterial and fungal genera with differential abundance between genotypes were selected for display. A summary of all differentially abundant taxa can be found in supplemental tables S9 (bacteria) and S10) fungi of the publication.  
#load packages
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(gridExtra)
library(ggpubr)

differntial_taxa <- read.table("WP3_relabu_difftaxa_gt_all.txt", row.names=1, header=TRUE, sep="\t") #this dataframe includes the relative abundances of each differentially abundant taxon for each of the 8 replicates. Samples are in columns and taxa in rows.
View(differntial_taxa)

#bacteria
p3= ggplot(data=differntial_taxa, aes(x=rootstock, y=Allo_Neo_Para_Rhizobium, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +     
  coord_cartesian(ylim=c(0.0,4.0))+
  labs(x= "", y= "RA [%]", title= "Allo_Neo_Para_Rhizobium") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45)) 

p17= ggplot(data=differntial_taxa, aes(x=rootstock, y=Novosphingobium, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,4.0))+
  labs(x= "", y= "", title= "Novosphingobium") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

p52= ggplot(data=differntial_taxa, aes(x=rootstock, y=Sphingobium, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,4.0))+
  labs(x= "", y= "", title= "Sphingobium") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

p25= ggplot(data=differntial_taxa, aes(x=rootstock, y=Streptomyces, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,4.0))+
  labs(x= "", y= "", title= "Streptomyces") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

#fungi
p38= ggplot(data=differntial_taxa, aes(x=rootstock, y=Nectria, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,60.0))+
  labs(x= "", y= "RA [%]", title= "Nectria") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

p45= ggplot(data=differntial_taxa, aes(x=rootstock, y=Sarocladium, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,18.0))+
  labs(x= "", y= "", title= "Sarocladium") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

p46= ggplot(data=differntial_taxa, aes(x=rootstock, y=Solicoccozyma, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,18.0))+
  labs(x= "", y= "", title= "Solicoccozyma") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

p49= ggplot(data=differntial_taxa, aes(x=rootstock, y=Trichoderma, fill=site)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("#41917C", "#DFA398","#505778")) +        
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) + 
  theme_bw(base_size = 8) +
  #theme(axis.text.x=element_blank()) +     # delete x-axis labels 
  coord_cartesian(ylim=c(0.0,4.0))+
  labs(x= "", y= "", title= "Trichoderma") +                                       
  scale_x_discrete(guide = guide_axis(angle = 45))

ggarrange(p3,p17,p52,p25,p38,p45,p46,p49, ncol=4, nrow=2, common.legend=TRUE, legend="bottom")     
ggsave("WP3_boxplots_gt_responder_revised.png", width = 16.0, height = 10, units = "cm", dpi = 600)