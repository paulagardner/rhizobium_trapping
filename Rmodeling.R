#!/usr/bin/env Rscript
# Load the model from file
library("tidyverse")
library("poolfstat")


matrix_rds <- read_rds('pairwisefsts.rds')
print("matrix loaded, checking")
ls()
# Use the model however you want

