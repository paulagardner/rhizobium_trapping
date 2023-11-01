#!/usr/bin/env Rscript
# plot PCA
# cjfiscus
# 2023-03-15

# libs
library(pacman)
p_load(ggplot2, cowplot, ggthemes, scales, readxl, ggmap, maps, 
       mapdata, maptools, raster, rgeos, magrittr)

setwd("~/Desktop")

# import pca data
df<-read.delim("plink2.eigenvec")
eg<-read.table("plink2.eigenval")

# calc prop var explained
eg$pc<-1:nrow(eg)
eg$prop<-eg$V1/sum(eg$V1)

## plot prop var exp per pc
p1<-ggplot(eg, aes(x=pc, y=prop*100)) + geom_bar(stat="identity") +
  xlab("PC") +
  ylab("% variance explained") + theme_cowplot() +
  scale_x_continuous(breaks=seq(1,10,1))
ggsave("propvar.pdf", p1, height=4, width=4)

# plot pca with admix groups
#assign<-read.delim("cluster_assignments.txt")
#names(df)[2]<-"sample"
#df<-merge(df, assign, by="sample")

# pc 1 vs pc 2
p1<-ggplot(df, aes(x=PC1, y=PC2)) + geom_point() + theme_cowplot() +
  xlab(paste0("PC 1 (", round(eg[1,3]*100), "%)")) +
  ylab(paste0("PC 2 (", round(eg[2,3]*100), "%)"))
ggsave("pc1pc2.pdf", p1, height=4, width=4)

# plot by BIO8
wc<-read.delim("Varizonica_worldclim2_1.txt")
names(df)[1]<-"ID"
df<-merge(df, wc, by="ID")

p1<-ggplot(df, aes(x=PC1, y=PC2, color=BIO8)) + geom_point(alpha=0.7) + theme_cowplot() +
  xlab(paste0("PC 1 (", round(eg[1,3]*100), "%)")) +
  ylab(paste0("PC 2 (", round(eg[2,3]*100), "%)")) + 
  scale_colour_gradient_tableau(palette="Green-Gold")
ggsave("pc1pc2_bio8.pdf", p1, height=4, width=5)