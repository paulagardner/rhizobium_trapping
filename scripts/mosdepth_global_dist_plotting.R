#!/usr/bin/env Rscript

library(ggplot2)

input_file <- "~/rhizobium_trapping/data/test_output.mosdepth.global.dist.txt"

#plotting formatting taken from github issue that suggested the dist output in the first place!: https://github.com/MultiQC/MultiQC/issues/924
df <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("chr", "coverage", "percent"))
df$diff <- c(df$percent[1], diff(df$percent))
df_total <- subset(df, chr == "total")
df_total$diff <- c(df_total$percent[1], diff(df_total$percent))

## Cumulative sum plot
ggplot(data = subset(df_total, percent > 0.02)) + geom_line(aes(x = coverage, y = percent))

## Histogram-like plot
ggplot(data = subset(df_total, percent > 0.02)) + geom_line(aes(x = coverage, y = diff))

### plots by chromosome: 
ggplot(data = subset(df, percent > 0.02)) + geom_line(aes(x = coverage, y = percent)) + facet_wrap(~chr, scales="free")