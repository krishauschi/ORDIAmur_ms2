Data import
dat <- read.csv("gene_expression_WP3_BIS3_tidy.csv",
                header=TRUE,
                na.strings = "NA",
                strip.white = TRUE,
                dec = ".",
                sep= ";") %>%
  na.omit()
str(dat)
## ’data.frame’: 105 obs. of 5 variables:
## $ Genotyp : chr "M.26" "M.26" "M.26" "M.26" ...
## $ Boden : chr "HO" "HO" "HO" "HO" ...
## $ Behandlung: chr "ARD" "ARD" "ARD" "G" ...
## $ Pool : chr "M26_HO_ARD_01" "M26_HO_ARD_02" "M26_HO_ARD_03" "M26_HO_G_01" ...
## $ BIS3 : num 0.2621 0.6725 1.75 0.0896 0.0732 ...
## - attr(*, "na.action")= ’omit’ Named int [1:3] 13 31 38
## ..- attr(*, "names")= chr [1:3] "13" "31" "38"
Grafical overview
# Overview on response scale
ggplot(data=dat, aes(x=Behandlung, y=BIS3))+
  theme_bw()+
  geom_point(aes(color=Boden))+
  facet_grid(~Genotyp)

With increasing gene expression, also the variance increases. Log-transformation will stabilize the variance.
# Overview on log scale
ggplot(data=dat, aes(x=Behandlung, y=log(BIS3)))+
  theme_bw()+
  geom_point(aes(color=Boden))+
  facet_grid(~Genotyp)

Modelling
Modelling was done on log scale. Since, due to the pooling the randomization structure was not available
anymore and all reactions were run on the same plate, a simple linear model was fit using the treatments,
the genotypes and the different soils as well as their two- and threeway interactions as explanatory variables.
fit <- lm(log(BIS3) ~ Behandlung*Genotyp*Boden,
          dat)
Both analytical plots look ok.
# Residuals vs fitted
plot(fit,
     which=1)

# qq-Plot
plot(fit,
     which=2)

ANOVA
anova(fit)
## Analysis of Variance Table
##
## Response: log(BIS3)
## Df Sum Sq Mean Sq F value Pr(>F)
## Behandlung 1 74.529 74.529 193.7977 < 2.2e-16 ***
## Genotyp 2 29.836 14.918 38.7917 5.133e-12 ***
## Boden 5 5.264 1.053 2.7375 0.02582 *
## Behandlung:Genotyp 2 0.214 0.107 0.2788 0.75753
## Behandlung:Boden 5 1.241 0.248 0.6452 0.66602
## Genotyp:Boden 10 28.088 2.809 7.3036 9.252e-08 ***
## Behandlung:Genotyp:Boden 10 0.654 0.065 0.1701 0.99778
## Residuals 69 26.535 0.385
## ---
## Signif. codes: 0 ’***’ 0.001 ’**’ 0.01 ’*’ 0.05 ’.’ 0.1 ’ ’ 1
All three main effects as well as the genotype:soil interaction are significant.
Multiple mean comparisons
Treatment effect
The treatment effect (Behandlung) is the only main effect that is significant in the ANOVA but not involved
in interactions. Therefore, it can be deduced that the treatments differ significantly from each other and
this difference is constant over all combinations of soils and genotypes.
The mean comparisons were calculated as differences on log scale ˆμARD − ˆμG. Back-transformation to
response scale yields the proportion ˆμARD/ˆμG. Hence also the null hypothesis changes to μARD/μG = 1.
comp_treat <- emmeans(fit,
                      specs="Behandlung",
                      contr="pairwise",
                      adjust="mvt",
                      type="response")
comp_treat
## $emmeans
## Behandlung response SE df lower.CL upper.CL
## ARD 1.089 0.0957 69 0.914 1.30
## G 0.195 0.0164 69 0.164 0.23
##
## Results are averaged over the levels of: Genotyp, Boden
## Confidence level used: 0.95
## Intervals are back-transformed from the log scale
##
## $contrasts
## contrast ratio SE df null t.ratio p.value
## ARD / G 5.6 0.682 69 1 14.141 <.0001
##
## Results are averaged over the levels of: Genotyp, Boden
## Tests are performed on the log scale
The relative expression of BIS3 is estimated to be 5.6 times higher in ARD than in gamma irradiated soil.
Since the corresponding confidence interval does not cover the null hypothesis of one, it can be deduces that
the expression differs significantly between the two treatments.
plot(confint(comp_treat$contrasts))+
  theme_bw()+
  geom_vline(xintercept=1)+
  scale_x_continuous(breaks=c(1, 4, 5, 6, 7,8),
                     labels=c(1, 4, 5, 6, 7, 8),
                     limits=c(0, 8))

Genotype-soil interaction
# simple contrasts
comp_simple_gs <- emmeans(fit,
                          specs="Boden",
                          by="Genotyp",
                          contr="pairwise",
                          adjust="mvt",
                          type="response")
comp_simple_gs

## $emmeans
## Genotyp = M.26:
## Boden response SE df lower.CL upper.CL
## EH 1.010 0.2556 69 0.6092 1.673
## HG 0.978 0.2768 69 0.5561 1.720
## HO 0.258 0.0654 69 0.1558 0.428
## ME 0.972 0.2460 69 0.5864 1.610
## PI 1.337 0.3385 69 0.8068 2.215
## RU 0.613 0.1734 69 0.3483 1.078
##
## Genotyp = MAL0130:
## Boden response SE df lower.CL upper.CL
## EH 0.473 0.1198 69 0.2855 0.784
## HG 0.598 0.1514 69 0.3608 0.991
## HO 0.872 0.2468 69 0.4958 1.534
## ME 0.777 0.1967 69 0.4687 1.287
## PI 0.533 0.1349 69 0.3216 0.883
## RU 0.368 0.0932 69 0.2222 0.610
##
## Genotyp = MAL0739:
## Boden response SE df lower.CL upper.CL
## EH 0.121 0.0305 69 0.0728 0.200
## HG 0.176 0.0445 69 0.1061 0.291
## HO 1.108 0.2806 69 0.6688 1.836
## ME 0.156 0.0396 69 0.0943 0.259
## PI 0.220 0.0556 69 0.1325 0.364
## RU 0.140 0.0355 69 0.0845 0.232
##
## Results are averaged over the levels of: Behandlung
## Confidence level used: 0.95
## Intervals are back-transformed from the log scale
##
## $contrasts
## Genotyp = M.26:
## contrast ratio SE df null t.ratio p.value
## EH / HG 1.032 0.3920 69 1 0.084 1.0000
## EH / HO 3.910 1.4000 69 1 3.809 0.0040
## EH / ME 1.039 0.3720 69 1 0.107 1.0000
## EH / PI 0.755 0.2704 69 1 -0.784 0.9692
## EH / RU 1.648 0.6258 69 1 1.315 0.7753
## HG / HO 3.788 1.4386 69 1 3.507 0.0099
## HG / ME 1.006 0.3822 69 1 0.017 1.0000
## HG / PI 0.732 0.2778 69 1 -0.823 0.9622
## HG / RU 1.596 0.6391 69 1 1.169 0.8498
## HO / ME 0.266 0.0951 69 1 -3.702 0.0055
## HO / PI 0.193 0.0691 69 1 -4.593 0.0003
## HO / RU 0.421 0.1600 69 1 -2.275 0.2178
## ME / PI 0.727 0.2603 69 1 -0.891 0.9473
## ME / RU 1.586 0.6024 69 1 1.215 0.8278
## PI / RU 2.182 0.8287 69 1 2.055 0.3225
##
## Genotyp = MAL0130:
## contrast ratio SE df null t.ratio p.value
## EH / HG 0.791 0.2833 69 1 -0.654 0.9862
## EH / HO 0.542 0.2060 69 1 -1.611 0.5940
## EH / ME 0.609 0.2180 69 1 -1.385 0.7355
## EH / PI 0.888 0.3178 69 1 -0.333 0.9994
## EH / RU 1.285 0.4599 69 1 0.699 0.9814
## HG / HO 0.686 0.2603 69 1 -0.994 0.9182
## HG / ME 0.770 0.2756 69 1 -0.731 0.9774
## HG / PI 1.122 0.4017 69 1 0.321 0.9995
## HG / RU 1.624 0.5813 69 1 1.353 0.7539
## HO / ME 1.123 0.4264 69 1 0.305 0.9996
## HO / PI 1.636 0.6215 69 1 1.297 0.7854
## HO / RU 2.368 0.8993 69 1 2.270 0.2202
## ME / PI 1.458 0.5219 69 1 1.052 0.8981
## ME / RU 2.109 0.7552 69 1 2.085 0.3071
## PI / RU 1.447 0.5181 69 1 1.032 0.9053
##
## Genotyp = MAL0739:
## contrast ratio SE df null t.ratio p.value
## EH / HG 0.686 0.2456 69 1 -1.052 0.8982
## EH / HO 0.109 0.0390 69 1 -6.195 <.0001
## EH / ME 0.772 0.2763 69 1 -0.724 0.9783
## EH / PI 0.549 0.1966 69 1 -1.674 0.5532
## EH / RU 0.861 0.3083 69 1 -0.418 0.9983
## HG / HO 0.159 0.0568 69 1 -5.142 <.0001
## HG / ME 1.125 0.4027 69 1 0.328 0.9995
## HG / PI 0.800 0.2866 69 1 -0.622 0.9891
## HG / RU 1.255 0.4494 69 1 0.635 0.9880
## HO / ME 7.091 2.5387 69 1 5.471 <.0001
## HO / PI 5.046 1.8066 69 1 4.521 0.0003
## HO / RU 7.913 2.8331 69 1 5.777 <.0001
## ME / PI 0.712 0.2548 69 1 -0.950 0.9317
## ME / RU 1.116 0.3996 69 1 0.306 0.9996
## PI / RU 1.568 0.5615 69 1 1.257 0.8071
##
## Results are averaged over the levels of: Behandlung
## P value adjustment: mvt method for 15 tests
## Tests are performed on the log scale
plot(confint(comp_simple_gs$contrasts))+
  theme_bw()+
  geom_vline(xintercept=1)+
  scale_x_continuous(trans = "log10",
                     breaks=c(0.01, 0.1, 1, 10, 100, 200),
                     labels=c(0.01, 0.1, 1, 10, 100, 200))

Means for interaction contrasts
means_inter_gs <- emmeans(fit,
                          specs= ~ Boden:Genotyp,
                          type="response")
means_inter_gs
## Boden Genotyp response SE df lower.CL upper.CL
## EH M.26 1.010 0.2560 69 0.6092 1.673
## HG M.26 0.978 0.2770 69 0.5561 1.720
## HO M.26 0.258 0.0654 69 0.1558 0.428
## ME M.26 0.972 0.2460 69 0.5864 1.610
## PI M.26 1.337 0.3380 69 0.8068 2.215
## RU M.26 0.613 0.1730 69 0.3483 1.078
## EH MAL0130 0.473 0.1200 69 0.2855 0.784
## HG MAL0130 0.598 0.1510 69 0.3608 0.991
## HO MAL0130 0.872 0.2470 69 0.4958 1.534
## ME MAL0130 0.777 0.1970 69 0.4687 1.287
## PI MAL0130 0.533 0.1350 69 0.3216 0.883
## RU MAL0130 0.368 0.0932 69 0.2222 0.610
## EH MAL0739 0.121 0.0305 69 0.0728 0.200
## HG MAL0739 0.176 0.0445 69 0.1061 0.291
## HO MAL0739 1.108 0.2810 69 0.6688 1.836
## ME MAL0739 0.156 0.0396 69 0.0943 0.259
## PI MAL0739 0.220 0.0556 69 0.1325 0.364
## RU MAL0739 0.140 0.0355 69 0.0845 0.232
##
## Results are averaged over the levels of: Behandlung
## Confidence level used: 0.95
## Intervals are back-transformed from the log scale
# Interaction contrasts
comp_inter_gs <- contrast(means_inter_gs,
                          interaction=c(Genotyp="pairwise",
                                        Boden="pairwise"),
                          adj="mvt",
                          type="response")
comp_inter_gs
## Genotyp_pairwise Boden_pairwise ratio SE df null t.ratio p.value
## M.26 / MAL0130 EH / HG 1.3046 0.6809 69 1 0.509 1.0000
## M.26 / MAL0739 EH / HG 1.5046 0.7853 69 1 0.783 0.9993
## MAL0130 / MAL0739 EH / HG 1.1533 0.5840 69 1 0.282 1.0000
## M.26 / MAL0130 EH / HO 7.2089 3.7625 69 1 3.785 0.0118
## M.26 / MAL0739 EH / HO 35.9309 18.1933 69 1 7.073 <.0001
## MAL0130 / MAL0739 EH / HO 4.9843 2.6014 69 1 3.078 0.0887
## M.26 / MAL0130 EH / ME 1.7059 0.8638 69 1 1.055 0.9922
## M.26 / MAL0739 EH / ME 1.3463 0.6817 69 1 0.587 1.0000
## MAL0130 / MAL0739 EH / ME 0.7892 0.3996 69 1 -0.467 1.0000
## M.26 / MAL0130 EH / PI 0.8507 0.4308 69 1 -0.319 1.0000
## M.26 / MAL0739 EH / PI 1.3752 0.6963 69 1 0.629 0.9999
## MAL0130 / MAL0739 EH / PI 1.6165 0.8185 69 1 0.948 0.9966
## M.26 / MAL0130 EH / RU 1.2829 0.6696 69 1 0.477 1.0000
## M.26 / MAL0739 EH / RU 1.9136 0.9988 69 1 1.243 0.9734
## MAL0130 / MAL0739 EH / RU 1.4917 0.7553 69 1 0.790 0.9993
## M.26 / MAL0130 HG / HO 5.5257 2.9676 69 1 3.183 0.0669
## M.26 / MAL0739 HG / HO 23.8809 12.4640 69 1 6.080 <.0001
## MAL0130 / MAL0739 HG / HO 4.3218 2.2556 69 1 2.804 0.1674
## M.26 / MAL0130 HG / ME 1.3076 0.6825 69 1 0.514 1.0000
## M.26 / MAL0739 HG / ME 0.8948 0.4670 69 1 -0.213 1.0000
## MAL0130 / MAL0739 HG / ME 0.6843 0.3465 69 1 -0.749 0.9995
## M.26 / MAL0130 HG / PI 0.6521 0.3403 69 1 -0.819 0.9990
## M.26 / MAL0739 HG / PI 0.9140 0.4770 69 1 -0.172 1.0000
## MAL0130 / MAL0739 HG / PI 1.4016 0.7097 69 1 0.667 0.9998
## M.26 / MAL0130 HG / RU 0.9833 0.5281 69 1 -0.031 1.0000
## M.26 / MAL0739 HG / RU 1.2719 0.6831 69 1 0.448 1.0000
## MAL0130 / MAL0739 HG / RU 1.2934 0.6549 69 1 0.508 1.0000
## M.26 / MAL0130 HO / ME 0.2366 0.1235 69 1 -2.761 0.1836
## M.26 / MAL0739 HO / ME 0.0375 0.0190 69 1 -6.486 <.0001
## MAL0130 / MAL0739 HO / ME 0.1583 0.0826 69 1 -3.531 0.0254
## M.26 / MAL0130 HO / PI 0.1180 0.0616 69 1 -4.094 0.0046
## M.26 / MAL0739 HO / PI 0.0383 0.0194 69 1 -6.444 <.0001
## MAL0130 / MAL0739 HO / PI 0.3243 0.1693 69 1 -2.157 0.5296
## M.26 / MAL0130 HO / RU 0.1779 0.0956 69 1 -3.214 0.0618
## M.26 / MAL0739 HO / RU 0.0533 0.0278 69 1 -5.619 <.0001
## MAL0130 / MAL0739 HO / RU 0.2993 0.1562 69 1 -2.311 0.4248
## M.26 / MAL0130 ME / PI 0.4987 0.2525 69 1 -1.374 0.9479
## M.26 / MAL0739 ME / PI 1.0214 0.5172 69 1 0.042 1.0000
## MAL0130 / MAL0739 ME / PI 2.0482 1.0371 69 1 1.416 0.9375
## M.26 / MAL0130 ME / RU 0.7520 0.3925 69 1 -0.546 1.0000
## M.26 / MAL0739 ME / RU 1.4214 0.7418 69 1 0.674 0.9998
## MAL0130 / MAL0739 ME / RU 1.8901 0.9570 69 1 1.257 0.9713
## M.26 / MAL0130 PI / RU 1.5079 0.7870 69 1 0.787 0.9993
## M.26 / MAL0739 PI / RU 1.3915 0.7263 69 1 0.633 0.9999
## MAL0130 / MAL0739 PI / RU 0.9228 0.4673 69 1 -0.159 1.0000
##
## Results are averaged over the levels of: Behandlung
## P value adjustment: mvt method for 45 tests
## Tests are performed on the log scale
plot(comp_inter_gs)+
  geom_vline(xintercept=1)+
  theme_bw()+
  scale_x_continuous(trans = "log10",
                     breaks=c(0.01, 0.1, 1, 10, 100, 200),
                     labels=c(0.01, 0.1, 1, 10, 100, 200))

Three way interaction
means_inter_3way <- emmeans(fit,
                            specs="Behandlung" ,
                            by=c("Boden","Genotyp"),
                            contr="pairwise",
                            adjust="mvt",
                            type="response")
data.frame(confint(means_inter_3way$contrasts))
## contrast Boden Genotyp ratio SE df lower.CL upper.CL
## 1 ARD / G EH M.26 4.424774 2.240442 69 1.611390 12.150142
## 2 ARD / G HG M.26 8.496081 4.809682 69 2.746300 26.283871
## 3 ARD / G HO M.26 6.848490 3.467669 69 2.494046 18.805512
## 4 ARD / G ME M.26 10.374852 5.253209 69 3.778258 28.488675
## 5 ARD / G PI M.26 4.490065 2.273502 69 1.635168 12.329429
## 6 ARD / G RU M.26 6.085515 3.445046 69 1.967101 18.826431
## 7 ARD / G EH MAL0130 5.137379 2.601263 69 1.870903 14.106912
## 8 ARD / G HG MAL0130 9.594212 4.857939 69 3.493968 26.345087
## 9 ARD / G HO MAL0130 5.057922 2.863320 69 1.634938 15.647422
## 10 ARD / G ME MAL0130 4.681197 2.370280 69 1.704773 12.854266
## 11 ARD / G PI MAL0130 3.367890 1.705299 69 1.226500 9.248008
## 12 ARD / G RU MAL0130 4.339214 2.197120 69 1.580232 11.915203
## 13 ARD / G EH MAL0739 5.183138 2.624433 69 1.887567 14.232561
## 14 ARD / G HG MAL0739 6.549187 3.316119 69 2.385048 17.983645
## 15 ARD / G HO MAL0739 5.062168 2.563181 69 1.843513 13.900385
## 16 ARD / G ME MAL0739 5.459357 2.764294 69 1.988159 14.991043
## 17 ARD / G PI MAL0739 4.014033 2.032468 69 1.461809 11.022276
## 18 ARD / G RU MAL0739 6.228506 3.153745 69 2.268264 17.103076


