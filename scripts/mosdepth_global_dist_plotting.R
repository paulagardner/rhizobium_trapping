#!/usr/bin/env Rscript

library(ggplot2)

input_file <- "~/rhizobium_trapping/data/test_output.mosdepth.global.dist.txt"

df <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("chr", "coverage", "percent"))
df <- subset(df, chr == "total")
df$diff <- c(df$percent[1], diff(df$percent))

## Cumulative sum plot
ggplot(data = subset(df, percent > 0.02)) + geom_line(aes(x = coverage, y = percent))

## Histogram-like plot
ggplot(data = subset(df, percent > 0.02)) + geom_line(aes(x = coverage, y = diff))
