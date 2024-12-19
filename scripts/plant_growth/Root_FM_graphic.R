#-------------------------------------------------------------------------------
#------------- ARD Vergleich Frischmasse ---------------------------------------
#-------------------------------------------------------------------------------
### Meta data

# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#         [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Berlin
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] emmeans_1.8.8   lmerTest_3.1-3  ggpubr_0.6.0    lme4_1.1-35.3   Matrix_1.6-5    ggh4x_0.2.8    
# [7] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.4    
# [13] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#         [1] gtable_0.3.5        rstatix_0.7.2       lattice_0.21-8      tzdb_0.4.0          numDeriv_2016.8-1.1
# [6] vctrs_0.6.5         tools_4.3.1         generics_0.1.3      parallel_4.3.1      pbkrtest_0.5.2     
# [11] sandwich_3.1-0      fansi_1.0.6         pkgconfig_2.0.3     lifecycle_1.0.4     farver_2.1.1       
# [16] compiler_4.3.1      textshaping_0.3.6   munsell_0.5.1       codetools_0.2-19    carData_3.0-5      
# [21] pillar_1.9.0        car_3.1-2           nloptr_2.0.3        MASS_7.3-60         boot_1.3-28.1      
# [26] abind_1.4-5         multcomp_1.4-25     nlme_3.1-162        tidyselect_1.2.1    mvtnorm_1.2-4      
# [31] stringi_1.8.3       labeling_0.4.3      splines_4.3.1       cowplot_1.1.3       grid_4.3.1         
# [36] colorspace_2.1-0    cli_3.6.2           magrittr_2.0.3      survival_3.5-5      utf8_1.2.4         
# [41] broom_1.0.5         TH.data_1.1-2       withr_3.0.0         scales_1.3.0        backports_1.4.1    
# [46] timechange_0.2.0    estimability_1.4.1  ggsignif_0.6.4      ragg_1.2.5          zoo_1.8-12         
# [51] hms_1.1.3           coda_0.19-4         rlang_1.1.3         Rcpp_1.0.12         xtable_1.8-4       
# [56] glue_1.7.0          rstudioapi_0.14     minqa_1.2.6         R6_2.5.1            systemfonts_1.0.4  

#-------------------------------------------------------------------------------

# load packages

library(tidyverse)
library(lme4)
library(ggpubr)
library(lmerTest)
library(emmeans)
library(forcats)

#-------------------------------------------------------------------------------
# Data import

setwd("")

dat <- read.csv2(".csv",
                 strip.white = TRUE) %>% 
        mutate(Genotype=factor(Genotype),
               Soil=factor(Soil),
               Treat=factor(Treat, levels=c("G", "ARD")),
               Plant=factor(Plant),
               Block=factor(Block),
               Table=factor(Table)) %>% 
        na.omit()


#-------------------------------------------------------------------------------

# Linear mixed model with tree-way interaction of Soil*Genotype*Treat as fixed (for stat. testing)
# and Table and Blocks as random effects (dependence structure of experimental design)

fit_RFM <- lmer(RFM ~ Soil*Genotype*Treat + (1|Table) + (1|Block), data=dat)
summary(fit_RFM)

# Random effects:
# Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.01137  0.1066  
# Table    (Intercept) 0.00000  0.0000  
# Residual             1.67502  1.2942  
# Number of obs: 765, groups:  Block, 4; Table, 2

# - The tables do not explain any variance
# - The Blocks explain a bit of the total variance
# - Most of the variance is residual

#-------------------------------------------------------------------------------

# Graphical display of the results


#Adjustments to ordering and labeling

inter3_means <- data.frame(emmeans(fit_RFM, specs=~Soil:Genotype:Treat))

dat$Genotype <- fct_relevel(dat$Genotype, 'M.26', 'MAL0130','MAL0739', 'EMR.2', 'G.202', 'G.935')
dat$Soil <- fct_relevel(dat$Soil, "HG","HO","EH", "PI", "RU", "ME")

inter3_means$Genotype <- fct_relevel(inter3_means$Genotype, 'M.26', 'MAL0130','MAL0739', 'EMR.2', 'G.202', 'G.935')
inter3_means$Soil <- fct_relevel(inter3_means$Soil, "HG","HO","EH", "PI", "RU", "ME")

Treatlabels <- c("γARD", "ARD", "γARD", "ARD","γARD", "ARD","γARD", "ARD","γARD", "ARD","γARD", "ARD")


# Observations, means, 95% pointwise CI

plot_RFM_3 <- ggplot(dat, aes(x=Treat, y=RFM, fontsize=14))+
  theme_bw()+
  geom_violin(
    aes(fill = Soil),  # Specify fill aesthetic here
    position = "dodge",
    trim = TRUE,
    scale = "area"
  ) +
  scale_fill_manual(values = c(
    HG = "#DFA398CC",
    HO = "#A20000CC",
    EH = "#41917CCC",
    PI = "#A9D284CC",
    RU = "#505778CC",
    ME = "#A0A5C0CC"
  )) +
  geom_pointrange(data=inter3_means,
                  aes(y=emmean,
                      ymax=upper.CL,
                      ymin=lower.CL),
                  size = 0.1,
                  stroke=0.2,
                  color="black",
                  fill="white",
                  shape=21)+
  #scale_color_manual(values = c(HG="#DFA398", HO="#A20000", EH="#41917C", PI="#A9D284", RU="#505778", ME="#A0A5C0"), aesthetics = c("colour", "fill"))+ 
  # facet_nested(Table + Block ~ Soil)
  facet_grid(Genotype ~ Soil)+
  scale_x_discrete(labels= Treatlabels)+
  ylab("root FM [g]")+
  xlab("Treatment")+
  ggtitle("B")+
  ylim(c(0, 15))+
  theme(
    strip.text.x = element_text(size=8),  # Adjust the size for Soil labels
    strip.text.y = element_text(size=8)   # Adjust the size for Genotype labels
  )

plot_RFM_3

ggsave("RFM_violin_final.jpg", 
       width=16,
       height=11,
       units="cm",
       dpi=500)






