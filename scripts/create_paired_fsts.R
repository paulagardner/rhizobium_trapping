#options(error = function() { traceback(2L)
#  dump.frames(dumpto = "last.dump", to.file = TRUE)
#})

#Allele Data ####
#what we'll use to calculate pooled data
library(poolfstat)
#tidyverse to save as rds at the end
library(tidyverse)

infile = 'data/sync_all_bams.sync' #name of sync file to analyze
gff_filename = 'data/genome.gff3' #Name of annotation file
gene_output_filename = 'gene_tests.tsv' #Name of file for gene-level tests
snp_output_filename = 'snp_tests.tsv' #name of file for gene-level tests
poolsizes = rep(870021, 24) #pool size setting

#Vector of treatments for testing allele frequency differences. MUST
#be in the same order as the columns of the sync file. Currently, this
#script is set up to test just one tratment factor. 


pooldata = popsync2pooldata(sync.file="~/Fudge016/TrappingCh2/rhizobium_trapping/data/sync_all_bams.sync", poolsizes=rep(872001,24))
print("pooldata object made")

test_pca <- randomallele.pca(pooldata, scale = TRUE, return.snploadings = FALSE, plot.pcs = c(1, 2),)
write_rds(p.fst, '~/Fudge016/TrappingCh2/rhizobium_trapping/data/test_pca.rds')
print("PCA calculated")



#use this to compute global and per SNP FSTs:
snp.fsts <- computeFST(pooldata, method="Anova")

print("test stop")

#use this to compute pairwise FSTs:
pair.fst <- compute.pairwiseFST(pooldata,method="Anova")

print("calculating pairwise fsts finished")

#r is very much not happy trying to represent pair.fst through print.default below
#print("pariwise fst matrix snippet", pair.fst)
#save pairwise matrix:

#new poolfstat data structure makes it so that pair.fst is a s4 object, not s3, so different symbol to access the matrix variable in the dataset 
p.fst <- pair.fst@PairwiseFSTmatrix
ls()
#print("Creating matrices finished")

write_rds(p.fst, '~/Fudge016/TrappingCh2/rhizobium_trapping/data/pairwisefsts.rds')
#p.fst<- as.matrix(p.fst)
#heatmap(p.fst, Rowv=NA, Colv=NA, symm=TRUE, scale="none")

#print("headmap made")

