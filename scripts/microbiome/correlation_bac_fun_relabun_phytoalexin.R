setwd("M:/ORDIAmur_Phase_II/WP3-Experiment/WP3_ITS_analysis/Wp3_ITS_spearman")

library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(colorspace)

WP3_ITS <- read.table("WP3_spearman_corr_diff-taxa_revised.txt", header=T, sep="\t", dec=".")
View(WP3_ITS)
WP3_cor <- cor(WP3_ITS)

# Perform correlation significance test
testRes = cor.mtest(WP3_ITS, conf.level = 0.95)

# Open a PNG device
png(filename = "WP3_Fig7_revised.png", width = 6.5, height = 6, units = 'in', res = 600)

# Set up layout for adding a legend below the main plot
layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1))  # Main plot above, legend below

# Draw the correlation plot
corrplot(WP3_cor,
         p.mat = testRes$p, 
         insig = 'blank', 
         method = "circle", 
         type = 'upper', 
         diag = FALSE, 
         tl.cex = 0.7, 
         tl.col = "black", 
         tl.srt = 90, 
         cl.cex = 0.7, 
         cl.pos = 'r', 
         col = COL2(diverging = c("BrBG"), n = 200))

# Calculate exact sizes for the legend dots
max_dot_size <- 1  # Maximum size in corrplot corresponds to correlation = 1
legend_sizes <- max_dot_size * abs(c(-1, -0.5, 0, 0.5, 1))

# Add a custom legend below with calculated sizes
par(mar = c(0, 0, 0, 0))  # No margins
plot.new()
legend("center", 
       legend = c("-1", "-0.5", "0", "0.5", "1"), 
       title = "Correlation Size",
       pt.cex = legend_sizes,  # Use calculated sizes
       col = "black", 
       pch = 16,  # Solid circles
       horiz = TRUE,
       bty = "n",  # No box around legend
       x.intersp = 1.2)  # Space between legend points

# Close the PNG device
dev.off()