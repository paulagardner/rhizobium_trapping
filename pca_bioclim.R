#!/usr/bin/env Rscript
# bioclim pca
# cjfiscus
# 2023-08-28

library(pacman)
p_load(ggplot2, cowplot, ggrepel, khroma, reshape2, ggbeeswarm)

setwd("~/Desktop")

# read in and prep data
df<-read.delim("vitis_worldclim2_1.txt")
row.names(df)<-df$ID
df$ID<-NULL
df<-na.omit(df)

# do pca
pca<-prcomp(df, center=T, scale=T)

##########
# extract variables
## per var exp
per.var<-as.data.frame(cbind(1:19, (pca$sdev^2/sum(pca$sdev^2))*100))
names(per.var)<-c("pc", "prop_var_exp")

scores<-as.data.frame(pca$x)

ld<-as.data.frame(pca$rotation)
##########
## plot per var exp
p1<-ggplot(per.var, aes(x=pc, y=prop_var_exp)) + 
  geom_bar(stat="identity") +
  theme_cowplot() +
  scale_x_continuous(breaks=seq(1,19,1)) + ylab("% variance explained") +
  xlab("principal component")
ggsave("per_var_exp.pdf", p1, height=4, width=4)

## plot loadings
temp<-c(paste0("BIO", seq(1:11)))
ld$lab<-sapply(strsplit(row.names(ld), split="_"), "[", 1)
ld$type<-ifelse(ld$lab %in% temp, "Temperature", "Precipitation")

p1<-ggplot(ld, aes(x=PC1, y=PC2, label=lab, color=type)) + 
  geom_point() +
  geom_text_repel() + 
  theme_cowplot() +
  theme(legend.position="top") + scale_color_bright()
ggsave("load_pc1_pc2.pdf", p1, height=4, width=4)

p1<-ggplot(ld, aes(x=PC1, y=PC3, label=lab, color=type)) + 
  geom_point() +
  geom_text_repel() + 
  theme_cowplot() +
  theme(legend.position="top") + scale_color_bright()
ggsave("load_pc1_pc3.pdf", p1, height=4, width=4)

## plot scores
scores$sample<-row.names(scores)

## apply kmeans clustering to determine cluster number
temp<-scores
temp$sample<-NULL

ss<-data.frame()
for (i in 1:10){
  sub<-kmeans(temp, centers=i, nstart=25)
  add<-as.data.frame(cbind(i, sub$tot.withinss))
  ss<-as.data.frame(rbind(ss, add))
}

## plot ss
p1<-ggplot(ss, aes(x=i, V2)) + 
  geom_line() + 
  geom_point() +
  theme_cowplot() +
  xlab("K clusters") + ylab("sum of squares within") +
  scale_x_continuous(breaks=seq(1,10))
ggsave("ss_within.pdf", p1, height=4, width=4)

## choose 2
sub<-kmeans(temp, centers=2, nstart=25)
assign<-as.data.frame(sub$cluster)
assign$ID<-row.names(assign)
names(assign)<-c("cluster", "sample")

## cluster in pc space
scores<-merge(scores, assign, by="sample")
scores$cluster<-as.factor(scores$cluster)

p1<-ggplot(scores, aes(x=PC1, y=PC2, color=cluster)) + geom_point(alpha=0.5) +
  theme_cowplot() +
  theme(legend.position="right") +
  xlab(paste0("PC 1 (", round(per.var[1,2], digits=2), "%)")) +
  ylab(paste0("PC 2 (", round(per.var[2,2], digits=2), "%)")) +
  scale_color_bright()
ggsave("pc1_pc2.pdf", p1, height=4, width=5)

p1<-ggplot(scores, aes(x=PC1, y=PC3, color=cluster)) + geom_point(alpha=0.5) +
  theme_cowplot() +
  theme(legend.position="right") +
  xlab(paste0("PC 1 (", round(per.var[1,2], digits=2), "%)")) +
  ylab(paste0("PC 3 (", round(per.var[3,2], digits=2), "%)")) +
  scale_color_bright()
ggsave("pc1_pc3.pdf", p1, height=4, width=5)

## plot summary by cluster
df$sample<-row.names(df)

m<-merge(df, scores, by="sample")

sub<-m[,c("sample", "cluster", "BIO1", "BIO12")]

sub<-melt(sub, id.vars=c("sample", "cluster"))

p1<-ggplot(sub, aes(x=cluster, y=value, color=cluster)) + 
  geom_quasirandom() + 
  scale_color_bright() + theme_cowplot() +
  facet_wrap(~variable, scales="free") +
  theme(strip.background = element_blank()) +
  stat_summary(fun = "median", colour = "black", size = 3, geom = "point")
ggsave("diff_climate.pdf", p1, height=4, width=6)

## write out data
write.table(scores, "scores_cluster.txt", sep="\t", quote=F, row.names=F)
