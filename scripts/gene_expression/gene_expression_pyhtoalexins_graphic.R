#-------------------------------------------------------------------------------
#----------------------- gene expression ---------------------------------------
#-------------------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(emmeans)
library(lmerTest)
library(ggpubr)
library(ggh4x)

setwd("D:/Versuche/WP3 Versuch 2021/Manuskript")

dat_gene <- read.csv("gene_expression_WP3_BIS3_tidy.csv",
                     header=TRUE,
                     na.strings = "NA",
                     strip.white = TRUE,
                     dec = ".",
                     sep= ";") %>% 
        mutate(Boden=factor(Boden, levels=c("HG", 
                                            "HO",
                                            "EH",
                                            "PI",
                                            "RU",
                                            "ME"))) %>% 
        na.omit()

# Fit model
fit <- lm(log(BIS3) ~ Behandlung*Genotyp*Boden,
          dat_gene)


three_way_means <- emmeans(fit, 
                           specs= ~Behandlung*Genotyp*Boden,
                           type="response") %>% 
        data.frame()

three_way_means

# arith_means <- dat_gene %>% group_by(Behandlung, Genotyp, Boden) %>% 
#         transmute(arith_mean=mean(BIS3))%>% 
#         data.frame()
# 
# means_re <- full_join(three_way_means, arith_means) %>% 
#         mutate(ls_mean=response,
#                mean_diff=ls_mean-arith_mean) %>% 
#         select(Behandlung,
#                Genotyp,
#                Boden,
#                ls_mean,
#                arith_mean,
#                mean_diff)
# 
# write.csv2(means_re, "relative_expression_means.csv")

# ggplot(data=dat_gene, aes(x=Behandlung, y=BIS3))+
#         theme_bw()+
#         geom_pointrange(data=three_way_means,
#                         aes(y=response,
#                             ymin=lower.CL,
#                             ymax=upper.CL))+
#         geom_jitter(aes(color=Boden),
#                     width=0.2)+
#         facet_nested(~Boden+Genotyp)+
#         scale_color_manual(values=c("HG"="#DFA398", 
#                                    "HO"="#883b41",
#                                    "EH"="#41917C",
#                                    "PI"="#A9D284",
#                                    "RU"="#505778",
#                                    "ME"="#A0A5C0"))


# Arithmethic mean for gamma treatment is 0.3301233
filter(dat_gene, Behandlung=="G") %>% 
        transmute(total_mean <- mean(BIS3))

# ls mean for gamma treatmen is 0.205
emmeans(fit, 
        specs= ~Behandlung,
        type="response") 


dat_rel_expr <- full_join(dat_gene, three_way_means) %>% 
        mutate(rel_expr= BIS3 - 0.205,
               treat="ARD (RE)") %>% 
        filter(Behandlung=="ARD")

dat_rel_expr_means <- dat_rel_expr %>% group_by(Boden, Genotyp) %>% 
        transmute(mean_rel_expr=mean(rel_expr)) %>% 
        mutate(treat="ARD (RE)") %>% 
        data.frame()




# Relative expression "within each soil genotype combination"
# dat_rel_expr <- full_join(dat_gene, three_way_means) %>% 
#         mutate(rel_expr= BIS3 - response,
#                treat="ARD (RE)") %>% 
#         filter(Behandlung=="ARD")
# 
# dat_rel_expr_means <- dat_rel_expr %>% group_by(Boden, Genotyp) %>% 
#         transmute(mean_rel_expr=mean(rel_expr)) %>% 
#         mutate(treat="ARD (RE)") %>% 
#         data.frame()

fig_a <- ggplot(data=dat_rel_expr, 
                aes(x=Genotyp, y=rel_expr))+
        theme_bw()+
        geom_jitter(aes(fill=Boden),
                    width=0.15,
                    shape=21)+
        geom_point(data=dat_rel_expr_means,
                   aes(y=mean_rel_expr,
                       x=Genotyp),
                   shape=21,
                   size=2,
                   fill="white",
                   alpha=0.4)+
        facet_nested(~Boden)+
        scale_fill_manual(values=c("HG"="#DFA398", 
                                    "HO"="#883b41",
                                    "EH"="#41917C",
                                    "PI"="#A9D284",
                                    "RU"="#505778",
                                    "ME"="#A0A5C0"))+
        ylab(expression("Relative "*italic("BIS3")*" expression"))+
        theme(axis.title.y = element_text(size = 9))+
        xlab("")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        labs(fill="Soil origin")

fig_a

#-------------------------------------------------------------------------------
#------------- Phytoalexins ----------------------------------------------------
#-------------------------------------------------------------------------------

dat <- read.csv2("WP3_Phytoalexine_t8.csv",
                 header=TRUE,
                 strip.white = TRUE,
                 sep=";",
                 dec=".") %>% 
        mutate(Soil=factor(Soil, levels=c("HG", 
                                          "HO",
                                          "EH",
                                          "PI",
                                          "RU",
                                          "ME")))
str(dat)

dat$Table <- "1"
dat$Table[which(dat$Block=="C")] <- "2"
dat$Table[which(dat$Block=="D")] <- "2"

#-------------------------------------------------------------------------------

dat$Phyto1 <- dat$Phytoalexin+1

fit <- lmer(log(Phyto1) ~ Treatment*Genotype* Soil+
                    (1|Table) + (1|Table:Block)+
                    (1|Table:Treatment) + (1|Table:Block:Treatment)+
                    (1|Table:Genotype) + (1|Table:Block:Genotype)+
                    (1|Table:Soil) + (1|Table:Block:Soil),
            data=dat)

#-------------------------------------------------------------------------------

# genotyp:soil interaktion

means_geno_soil <- emmeans(fit,
                           specs= ~Genotype:Soil,
                           type="response") %>% 
        data.frame()

means_treat_soil <- emmeans(fit,
                            specs= ~Treatment:Soil,
                            type="response") %>% 
        data.frame()

means_treat_soil_geno <- emmeans(fit,
                                 specs= ~Genotype:Treatment:Soil,
                                 type="response") %>% 
        data.frame()

fig_b <- ggplot(dat, aes(x=Genotype,
                         y=Phytoalexin))+
        theme_bw()+
        # geom_rect(data=merge(dat, means_treat_soil),
        #           aes(ymin = lower.CL  , ymax = upper.CL,
        #               xmin = -Inf, xmax = Inf), fill = "grey", alpha = 0.006) +
        geom_jitter(aes(fill=Soil,
                        shape=Treatment),
                    width=0.15,
                    alpha=0.8)+
        # geom_hline(data=means_treat_soil,
        #            aes(yintercept = response),
        #            color="gray")+
  #geom_pointrange(data=means_geno_soil,
   #               aes(y=response,
    #                  ymin=lower.CL,
     #                 ymax=upper.CL),
      #            shape=24,
       #           fill="white",
        #          size=0.2)+
        geom_pointrange(data=means_treat_soil_geno,
                        aes(y=response,
                            ymin=lower.CL,
                            ymax=upper.CL,
                            shape=Treatment),
                        fill="white",
                        position = position_dodge2(w = 0.75),
                        size=0.2,
                        fatten=10,
                        stroke=0.5)+
        facet_nested(~Soil)+
        scale_shape_manual(values=c("ARD"=21,
                                    "G"=22),
                           labels=c("ARD"="ARD",
                                    "G"=expression(paste(gamma, "ARD"))))+
        scale_fill_manual(values=c("HG"="#DFA398", 
                                   "HO"="#883b41",
                                   "EH"="#41917C",
                                   "PI"="#A9D284",
                                   "RU"="#505778",
                                   "ME"="#A0A5C0"))+
        scale_x_discrete(labels=c("M26"="M.26",
                                  "MAL0130"="MAL0130",
                                  "MAL0739"="MAL0739"))+
        guides(fill = "none")+
        ylab(expression("Root phytoalexins ["*mu*"g g"^-1*" DW]"))+
        theme(axis.title.y = element_text(size = 9))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig_b

#-------------------------------------------------------------------------------


ggarrange(fig_a,
          fig_b,
          nrow=2,
          labels = c("A", "B"))

ggsave("fig_3_final.png", width=16, height=14, units="cm")
