Comparing our methods to Jorrin and Imperial 2015 paper:

Methods we have done so far:
Filter reads with Trimmomatic (Kbase)
Align reads to Rlv reference genome
Transform output with Sam Tools
Detect SNPs (They used varscan, we used PoPoolation)

Methods they did, which we should do next:
Filter SNPs bassed on requirements (FST > 0.1, p < 0.05)
Calculate euclidean distances and multidimensional scaling based on output
