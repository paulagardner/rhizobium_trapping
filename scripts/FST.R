library(poolfstat)

pooldata <- vcf2pooldata(
  vcf.file = "C:/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/2025revisionsrepo/scripts/missingness_filtered_output.vcf",
  poolsizes = rep(50, 24), #100 is a conservative estimate of how many bacterial genomes may have been pooled in each sample
  poolnames = c(
    "AUMerit_Roseau_Cold.sorted", "AUMerit_Roseau_Warm.sorted",
    "AUMerit_Rosemount_Cold.sorted", "AUMerit_Rosemount_Warm.sorted",
    "AUMerit_StPaul_Cold.sorted", "AUMerit_StPaul_Warm.sorted",
    "Hungvillosa_Roseau_Cold.sorted", "Hungvillosa_Roseau_Warm.sorted",
    "Hungvillosa_Rosemount_Cold.sorted", "Hungvillosa_Rosemount_Warm.sorted",
    "Hungvillosa_StPaul_Cold.sorted", "Hungvillosa_StPaul_Warm.sorted",
    "MSP4045_Roseau_Cold.sorted", "MSP4045_Roseau_Warm.sorted",
    "MSP4045_Rosemount_Cold.sorted", "MSP4045_Rosemount_Warm.sorted",
    "MSP4045_StPaul_Cold.sorted", "MSP4045_StPaul_Warm.sorted",
    "PupleBounty_StPaul_Warm.sorted", "PurpleBounty_Roseau_Cold.sorted",
    "PurpleBounty_Roseau_Warm.sorted", "PurpleBounty_Rosemount_Cold.sorted",
    "PurpleBounty_Rosemount_Warm.sorted", "PurpleBounty_StPaul_Cold.sorted"
  ),
  min.rc = 1,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.01,
  nlines.per.readblock = 1e+06,
)

PairwiseFST <- compute.pairwiseFST(pooldata)

print(PairwiseFST@PairwiseFSTmatrix)

library(ggplot2)
library(reshape2)

# Convert the matrix to a long format for ggplot
fst_melted <- melt(PairwiseFST@PairwiseFSTmatrix)

# Create the heatmap
ggplot(fst_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Pairwise FST Heatmap", x = "Population", y = "Population", fill = "FST") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

write.csv(PairwiseFST@PairwiseFSTmatrix, file = "pairwise_FST.csv", row.names = TRUE)
