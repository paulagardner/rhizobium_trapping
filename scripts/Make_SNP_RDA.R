# library(vcfR)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
require(tidyr)
library(tidyverse)
library(wesanderson) # color schemes
library(patchwork) # display two figures at once


#R version being used in documentation-env: 4.3.3
print("packages loaded")
setwd("~/Fudge016/TrappingCh2/rhizobium_trapping/data")
getwd()
print("working directory set to data")

# Step 1: Import snp_freq_diffs_test_mincov2_rc
# How to import data in Rstudio:
# Click file -> "import dataset" -> select file -> comment = none
# can also do this through wd

# load in as a df. not sure why read.table doesn't maintain headers as I'd like
data <- read.delim("snp_freq_diffs_test_mincov2_rc", sep = "\t", header = T)
print(head(data))




#####################PUT ALLELE FREQUENCIES IN DECIMAL FORMAT ###########################################
convert_and_replace <- function(data, start_column_index) {
  # Get the columns from the specified index to the end
  fraction_columns <- names(data)[start_column_index:ncol(data)]
  
  for (col in fraction_columns) {
    # Extract numerator and denominator
    fractions <- strsplit(as.character(data[[col]]), "/")
    
    # Convert to decimal and replace the original entry
    data[[col]] <- sapply(fractions, function(fraction) as.numeric(fraction[1]) / as.numeric(fraction[2]))
  }
  
  return(data)
}

# Specify the start column index
start_column_index <- 10  # Change this to the appropriate index
snp_data <- convert_and_replace(data, start_column_index) #convert fraction rows to decimal rows

# Display the result
print("result_df made")
print(head(snp_data))


############################ENSURE HEADER NAMES ARE POPULATED VIA BAM FILE .TXT THAT POPOOLATION, SAMTOOLS WORK FROM#########################
lines <- readLines("sorted_bam_files.txt")

# Process each line
sample_names <- sub("^.*bamfiles/", "", lines)  # Remove everything before and including "bamfiles/"
sample_names <- sub("\\.sorted$", "", sample_names)  # Remove ".sorted"

# Print the processed lines to check the result
print("df of sample names from .bams")
print(sample_names)



# Example data frames
# sample_names 
# result_df
samples_df <- data.frame(sample_names)

# Prefix to add
prefix <- "major_alleles_"

# Use a for loop to rename the major allele columns in the snp data
for (i in 1:nrow(samples_df)) {
  snp_data_col_index <- i + 9  # Adjust index to start renaming from column 10
  colnames(snp_data)[snp_data_col_index] <- paste0(prefix, samples_df$sample_names[i])
}

print(head(snp_data))



# #### ALL Replicons ####

# # Import file
# load("result_df.RData")
# # For some reason, this is called result_df

# # drop rows with NA
# snp_data <- data.frame(result_df) 
# head(snp_data)
# snp_data <- snp_data %>% drop_na() 

# # Add new names
# colnames(snp_data)[1] = "Chr"
# colnames(snp_data)[2] = "Pos"
# colnames(snp_data)[3] = "rc"
# colnames(snp_data)[4] = "allele_count"
# colnames(snp_data)[5] = "allele_states"
# colnames(snp_data)[6] = "deletion_sum"
# colnames(snp_data)[7] = "snp_type"
# colnames(snp_data)[8] = "major_alleles(maa)"
# colnames(snp_data)[9] = "minor_alleles(mia)"
# colnames(snp_data)[10] = "AUMerit_Roseau_Cold"
# colnames(snp_data)[11] = "AUMerit_Roseau_Warm"
# colnames(snp_data)[12] = "AUMerit_Rosemount_Cold"
# colnames(snp_data)[13] = "AUMerit_Rosemount_Warm"
# colnames(snp_data)[14] = "AUMerit_StPaul_Cold"
# colnames(snp_data)[15] = "AUMerit_StPaul_Warm"
# colnames(snp_data)[16] = "Hungvillosa_Roseau_Cold"
# colnames(snp_data)[17] = "Hungvillosa_Roseau_Warm"
# colnames(snp_data)[18] = "Hungvillosa_Rosemount_Cold"
# colnames(snp_data)[19] = "Hungvillosa_Rosemount_Warm"
# colnames(snp_data)[20] = "Hungvillosa_StPaul_Cold"
# colnames(snp_data)[21] = "Hungvillosa_StPaul_Warm"
# colnames(snp_data)[22] = "MSP4045_Roseau_Cold"
# colnames(snp_data)[23] = "MSP4045_Roseau_Warm"
# colnames(snp_data)[24] = "MSP4045_Rosemount_Cold"
# colnames(snp_data)[25] = "MSP4045_Rosemount_Warm"
# colnames(snp_data)[26] = "MSP4045_StPaul_Cold"
# colnames(snp_data)[27] = "MSP4045_StPaul_Warm"
# colnames(snp_data)[28] = "PupleBounty_StPaul_Warm"
# colnames(snp_data)[29] = "PupleBounty_Roseau_Cold"
# colnames(snp_data)[30] = "PupleBounty_Roseau_Warm"
# colnames(snp_data)[31] = "PupleBounty_Rosemount_Cold"
# colnames(snp_data)[32] = "PupleBounty_Rosemount_Warm"
# colnames(snp_data)[33] = "PupleBounty_StPaul_Cold"
# colnames(snp_data)[34] = "mia_AUMerit_Roseau_Cold"
# colnames(snp_data)[35] = "mia_AUMerit_Roseau_Warm"
# colnames(snp_data)[36] = "mia_AUMerit_Rosemount_Cold"
# colnames(snp_data)[37] = "mia_AUMerit_Rosemount_Warm"
# colnames(snp_data)[38] = "mia_AUMerit_StPaul_Cold"
# colnames(snp_data)[39] = "mia_AUMerit_StPaul_Warm"
# colnames(snp_data)[40] = "mia_Hungvillosa_Roseau_Cold"
# colnames(snp_data)[41] = "mia_Hungvillosa_Roseau_Warm"
# colnames(snp_data)[42] = "mia_Hungvillosa_Rosemount_Cold"
# colnames(snp_data)[43] = "mia_Hungvillosa_Rosemount_Warm"
# colnames(snp_data)[44] = "mia_Hungvillosa_StPaul_Cold"
# colnames(snp_data)[45] = "mia_Hungvillosa_StPaul_Warm"
# colnames(snp_data)[46] = "mia_MSP4045_Roseau_Cold"
# colnames(snp_data)[47] = "mia_MSP4045_Roseau_Warm"
# colnames(snp_data)[48] = "mia_MSP4045_Rosemount_Cold"
# colnames(snp_data)[49] = "mia_MSP4045_Rosemount_Warm"
# colnames(snp_data)[50] = "mia_MSP4045_StPaul_Cold"
# colnames(snp_data)[51] = "mia_MSP4045_StPaul_Warm"
# colnames(snp_data)[52] = "mia_PupleBounty_StPaul_Warm"
# colnames(snp_data)[53] = "mia_PupleBounty_Roseau_Cold"
# colnames(snp_data)[54] = "mia_PupleBounty_Roseau_Warm"
# colnames(snp_data)[55] = "mia_PupleBounty_Rosemount_Cold"
# colnames(snp_data)[56] = "mia_PupleBounty_Rosemount_Warm"
# colnames(snp_data)[57] = "mia_PupleBounty_StPaul_Cold"

# # I actually think we don't need minor alleles!

# small_snp_data <- snp_data[,c(1:2,10:33)]
# small_snp_data$SNP <- paste(small_snp_data$Chr,small_snp_data$Pos,sep="_")
# head(small_snp_data)

# # Make last row the first row
# small_snp_data <- small_snp_data %>%
#   relocate(SNP, before="Chr")


# # Get rid of columns I don't need
# small_snp_data <- small_snp_data[,-c(2:3)]
# head(small_snp_data)

# # Make column 1 the row names

# RDA_input <- small_snp_data[,-1]
# rownames(RDA_input) <- small_snp_data[,1]
# head(RDA_input)

# RDA_input <- t(RDA_input)
# RDA_input <- as.data.frame(RDA_input)
# rownames(RDA_input)

# # extract and define metadata from snp matrix row names
# meta <- data.frame("ID" = rownames(RDA_input))
# meta <- separate(meta,col=ID, into = c("Ecotype","Site","Temp"), sep = "_")
# meta$ID <- rownames(RDA_input)
# meta <- as.data.frame(meta)

# SNP_rda <- rda(RDA_input ~ meta$Ecotype + meta$Site + meta$Temp)

# smry <- summary(SNP_rda)
# all_rda <- SNP_rda
# variance <- SNP_rda$CA$eig/SNP_rda$tot.chi*100
# # all_varp <- varpart(RDA_input, ~ meta$Ecotype, ~ meta$Site, ~ meta$Temp) # Var part not relevant

# # Summary:
# # Partitioning of variance:
# #                 Inertia   Proportion
# # Total           40411     1.0000
# # Constrained     13105     0.3243
# # Unconstrained   27307     0.6757

# # RDA1 explains 0.5347, RDA2 explains 0.1781, RDA3 explains 0.1369


# # What if we take out ecotype, as Liana suggested?
# SNP_rda_no_ecotype <- rda(RDA_input ~ meta$Site + meta$Temp)
# smry_no_ecotype <- summary(SNP_rda_no_ecotype)
# rda.model.no.ecotype<-anova(SNP_rda_no_ecotype, step=1000, perm.max=1000, by= "terms")

# # Summary:
# # Partitioning of variance:
# #                 Inertia Proportion
# # Total           40411     1.0000
# # Constrained      9137     0.2261
# # Unconstrained   31275     0.7739

# # Not sure if this is better or worse


# # What if we just look at site?
# SNP_rda_only_site <- rda(RDA_input ~ meta$Site)
# rda.model.only.site <-anova(SNP_rda_only_site, step=1000, perm.max=1000, by= "terms")

# summary_rda_only_site <- data.frame(Term=row.names(rda.model.only.site),
#                                      Df=rda.model.only.site$Df,
#                                      Prop.Var = round(rda.model.only.site$Variance/sum(rda.model.only.site$Variance),3),
#                                      Fstat=round(rda.model.only.site$F,2),
#                                      Pvalue=round(rda.model.only.site$`Pr(>F)`,3),
#                                      Radj=as.numeric(RsquareAdj(SNP_rda_only_site))[2])

# df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
# df2  <- data.frame(smry$biplot[,1:2])   # loadings for PC1 and PC2
# df1$Ecotype <- meta$Ecotype
# df1$Site <- meta$Site
# df1$Temp <- meta$Temp

# rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Site) +
#   geom_point(aes(color = Site, shape=Temp),size=2.5) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   scale_color_manual(values = wes_palette("Darjeeling1"), name = "Sites") +
#   labs(title = "RDA of Allele Frequencies", subtitle = "All sites") +
#   theme_bw()

# rownames(df2) <- gsub("sam\\$", "", rownames(df2))

# rda.biplot.ecotype <- rda.plot +
# #  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2),
# #               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
# #  geom_text(data=df2,
# #            aes(x=RDA1,y=RDA2,label=rownames(df2),
# #                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))),
# #            color="red", size=3) +
#   labs(title = "RDA of Allele Frequencies",
#        x = paste("RDA1 (", round(variance[1], 2),"%)"),
#        y = paste("RDA1 (", round(variance[2], 2),"%)"))

# rda.biplot.ecotype


# rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Site) +
#   geom_point(aes(color = Site)) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   scale_color_manual(values = wes_palette("Darjeeling1"), name = "Site")

# rownames(df2) <- gsub("sam\\$", "", rownames(df2))

# rda.biplot.site <- rda.plot +
# #  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2),
# #               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
# #  geom_text(data=df2,
# #            aes(x=RDA1,y=RDA2,label=rownames(df2),
# #                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))),
# #            color="red", size=3) +
#   labs(title = "RDA of Allele Frequencies", subtitle = "All Sites",
#        x = paste("RDA1 (", round(variance[1], 2),"%)"),
#        y = paste("RDA1 (", round(variance[2], 2),"%)"))

# rda.biplot.site

# rda.temp <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Temp) +
#   geom_point(aes(color = Temp)) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   scale_color_manual(values = wes_palette("Darjeeling1"), name = "Temp")

# rownames(df2) <- gsub("sam\\$", "", rownames(df2))

# rda.biplot.temp <- rda.temp +
#   geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2),
#                color="red", arrow=arrow(length=unit(0.01,"npc"))) +
#   geom_text(data=df2,
#             aes(x=RDA1,y=RDA2,label=rownames(df2),
#                 hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))),
#             color="red", size=3) +
#   labs(title = "RDA of Allele Frequencies", subtitle = "All Replicons (MAF = 0.05)",
#        x = paste("RDA1 (", round(variance[1], 2),"%)"),
#        y = paste("RDA1 (", round(variance[2], 2),"%)"))

# rda.biplot.temp


# # Run with anova ####
# all_results <- anova(all_rda)

# # Code from Liana script
# rda.model.all<-anova(all_rda, step=1000, perm.max=1000, by= "terms")

# summary_rda_all <- data.frame(Term=row.names(rda.model.all),
#                                  Df=rda.model.all$Df,
#                                  Prop.Var = round(rda.model.all$Variance/sum(rda.model.all$Variance),3),
#                                  Fstat=round(rda.model.all$F,2),
#                                  Pvalue=round(rda.model.all$`Pr(>F)`,3),
#                                  Radj=as.numeric(RsquareAdj(SNP_rda))[2])

# summary_rda_no_ecotype <- data.frame(Term=row.names(rda.model.no.ecotype),
#                               Df=rda.model.no.ecotype$Df,
#                               Prop.Var = round(rda.model.no.ecotype$Variance/sum(rda.model.no.ecotype$Variance),3),
#                               Fstat=round(rda.model.no.ecotype$F,2),
#                               Pvalue=round(rda.model.no.ecotype$`Pr(>F)`,3),
#                               Radj=as.numeric(RsquareAdj(SNP_rda_no_ecotype))[2])

# # Permutation test for rda under reduced model
# # Terms added sequentially (first to last)
# # Permutation: free
# # Number of permutations: 999

# # Model: rda(formula = RDA_input ~ meta$Ecotype + meta$Site + meta$Temp)
# # Df Variance      F Pr(>F)  
# # meta$Ecotype  3   3968.1 0.8235  0.628  
# # meta$Site     2   7811.8 2.4317  0.021 *
# # meta$Temp     1   1324.7 0.8247  0.480  
# # Residual     17   27306.8                

# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# # Roseau ####

# RDA_input %>%
#   filter(grepl("Roseau",rownames(RDA_input))) -> RDA_input_Roseau

# rownames(RDA_input_Roseau)
# # Success!!
# # Now for meta:

# meta_Roseau <- meta %>%
#   filter(grepl("Roseau",ID))

# # Run RDA
# SNP_rda_Roseau <- rda(RDA_input_Roseau ~ meta_Roseau$Ecotype + meta_Roseau$Temp)
# # SNP_rda_Roseau <- rda(RDA_input_Roseau ~ meta_Roseau$Temp)

# # Anova
# rda.model.Roseau<-anova(SNP_rda_Roseau, step=1000, perm.max=1000, by= "terms")


# # Model: rda(formula = RDA_input_Roseau ~ meta_Roseau$Ecotype + meta_Roseau$Temp)
# #                      Df  Variance F       Pr(>F)  
# # meta_Roseau$Ecotype  3   6265.4   1.3029  0.116  
# # meta_Roseau$Temp     1   2698.6   1.6835  0.071 .
# # Residual             3   4809.0 

# # Temp approaching significance: p = 0.071

# # LESS significant when I take out ecotype as a term...
# # Model: rda(formula = RDA_input_Roseau ~ meta_Roseau$Temp)
# #                   Df  Variance  F       Pr(>F)
# # meta_Roseau$Temp  1   2698.6    1.4621  0.129
# # Residual          6   11074.4  


# smry_Roseau <- summary(SNP_rda_Roseau)
# #variance_Roseau <- SNP_rda_Roseau$CA$eig/SNP_rda_Roseau$tot.chi*100
# #varp_Roseau <- varpart(RDA_input_Roseau, ~ meta_Roseau$Ecotype, ~ meta_Roseau$Temp)

# summary_rda_Roseau <- data.frame(Term=row.names(rda.model.Roseau),
#                               Df=rda.model.Roseau$Df,
#                               Prop.Var = round(rda.model.Roseau$Variance/sum(rda.model.Roseau$Variance),3),
#                               Fstat=round(rda.model.Roseau$F,2),
#                               Pvalue=round(rda.model.Roseau$`Pr(>F)`,3),
#                               Radj=as.numeric(RsquareAdj(SNP_rda_Roseau))[2])

# # Get ready to plot
# df1_Roseau  <- data.frame(smry_Roseau$sites[,1:2])       # PC1 and PC2
# df2_Roseau  <- data.frame(smry_Roseau$biplot[,1:2])   # loadings for PC1 and PC2
# df1_Roseau$Ecotype <- meta_Roseau$Ecotype
# df1_Roseau$Temp <- meta_Roseau$Temp

# rda.plot.Roseau <- ggplot(df1_Roseau, aes(x=RDA1, y=RDA2), group = Temp) +
#   geom_point(aes(color = Temp),size=3) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   scale_color_manual(values=c("lightskyblue2", "salmon1")) +
#   labs(title = "Site: Roseau") +
#   theme_bw()
  
# rda.biplot.Roseau <- rda.plot.Roseau +
#   #  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2),
#   #               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
#   #  geom_text(data=df2,
#   #            aes(x=RDA1,y=RDA2,label=rownames(df2),
#   #                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))),
#   #            color="red", size=3) +
#   labs(title = "RDA of Allele Frequencies", subtitle = "Site: Roseau",
#        x = paste("RDA1 (", round(variance_Roseau[1], 2),"%)"),
#        y = paste("RDA1 (", round(variance_Roseau[2], 2),"%)"))

# rda.biplot.Roseau


# # St Paul ####
# RDA_input_StPaul <- RDA_input %>%
#   filter(grepl("StPaul",rownames(RDA_input)))

# rownames(RDA_input_StPaul)
# # Success!!
# # Now for meta:

# meta_StPaul <- meta %>%
#   filter(grepl("StPaul",ID))

# # Run RDA
# SNP_rda_StPaul <- rda(RDA_input_StPaul ~ meta_StPaul$Ecotype + meta_StPaul$Temp)
# # SNP_rda_StPaul <- rda(RDA_input_StPaul ~ meta_StPaul$Temp)

# smry_StPaul <- summary(SNP_rda_StPaul)
# #variance_StPaul <- SNP_rda_StPaul$CA$eig/SNP_rda_StPaul$tot.chi*100
# #varp_StPaul <- varpart(RDA_input_StPaul, ~ meta_StPaul$Ecotype, ~ meta_StPaul$Temp)

# # Get ready to plot
# df1_StPaul  <- data.frame(smry_StPaul$sites[,1:2])       # PC1 and PC2
# df2_StPaul  <- data.frame(smry_StPaul$biplot[,1:2])   # loadings for PC1 and PC2
# df1_StPaul$Ecotype <- meta_StPaul$Ecotype
# df1_StPaul$Temp <- meta_StPaul$Temp

# rda.plot.StPaul <- ggplot(df1_StPaul, aes(x=RDA1, y=RDA2), group = Temp) +
#   geom_point(aes(color = Temp),size=3) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   scale_color_manual(values=c("lightskyblue2", "salmon1")) +
#   labs(title = "Site: St. Paul") +
#   theme_bw()

# rda.biplot.StPaul <- rda.plot.StPaul +
#   #  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2),
#   #               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
#   #  geom_text(data=df2,
#   #            aes(x=RDA1,y=RDA2,label=rownames(df2),
#   #                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))),
#   #            color="red", size=3) +
#   labs(title = "RDA of Allele Frequencies", subtitle = "Site: St. Paul",
#        x = paste("RDA1 (", round(variance_StPaul[1], 2),"%)"),
#        y = paste("RDA1 (", round(variance_StPaul[2], 2),"%)"))

# rda.biplot.StPaul


# # Permanova

# rda.model.StPaul<-anova(SNP_rda_StPaul, step=1000, perm.max=1000, by= "terms")

# # Model: rda(formula = RDA_input_StPaul ~ meta_StPaul$Ecotype + meta_StPaul$Temp)
# #                      Df Variance F       Pr(>F)
# # meta_StPaul$Ecotype  3  14349.0  0.8217  0.579
# # meta_StPaul$Temp     1   2505.9  0.4305  0.694
# # Residual             3  17462.0  

# # p even HIGHER when we take out ecotype term
# # Model: rda(formula = RDA_input_StPaul ~ meta_StPaul$Temp)
# #                   Df    Variance  F       Pr(>F)
# # meta_StPaul$Temp  1     2506      0.4726  0.702
# # Residual          6     31811  

# summary_rda_StPaul <- data.frame(Term=row.names(rda.model.StPaul),
#                                  Df=rda.model.StPaul$Df,
#                                  Prop.Var = round(rda.model.StPaul$Variance/sum(rda.model.StPaul$Variance),3),
#                                  Fstat=round(rda.model.StPaul$F,2),
#                                  Pvalue=round(rda.model.StPaul$`Pr(>F)`,3),
#                                  Radj=as.numeric(RsquareAdj(SNP_rda_StPaul))[2])



# # Rosemount ####
# RDA_input_Rosemount <- RDA_input %>%
#   filter(grepl("Rosemount",rownames(RDA_input)))

# rownames(RDA_input_Rosemount)
# # Success!!
# # Now for meta:

# meta_Rosemount <- meta %>%
#   filter(grepl("Rosemount",ID))

# # Run RDA
# SNP_rda_Rosemount <- rda(RDA_input_Rosemount ~ meta_Rosemount$Ecotype + meta_Rosemount$Temp)
# # SNP_rda_Rosemount <- rda(RDA_input_Rosemount ~ meta_Rosemount$Temp)

# smry_Rosemount <- summary(SNP_rda_Rosemount)
# variance_Rosemount <- SNP_rda_Rosemount$CA$eig/SNP_rda_Rosemount$tot.chi*100
# varp_Rosemount <- varpart(RDA_input_Rosemount, ~ meta_Rosemount$Ecotype, ~ meta_Rosemount$Temp)

# # Get ready to plot
# df1_Rosemount  <- data.frame(smry_Rosemount$sites[,1:2])       # PC1 and PC2
# df2_Rosemount  <- data.frame(smry_Rosemount$biplot[,1:2])   # loadings for PC1 and PC2
# df1_Rosemount$Ecotype <- meta_Rosemount$Ecotype
# df1_Rosemount$Temp <- meta_Rosemount$Temp

# rda.plot.Rosemount <- ggplot(df1_Rosemount, aes(x=RDA1, y=RDA2), group = Temp) +
#   geom_point(aes(color = Temp),size=3) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   geom_vline(xintercept=0, linetype="dotted") +
#   scale_color_manual(values=c("lightskyblue2", "salmon1")) +
#   labs(title = "Site: Rosemount") +
#   theme_bw()

# rda.biplot.Rosemount <- rda.plot.Rosemount +
#   #  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2),
#   #               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
#   #  geom_text(data=df2,
#   #            aes(x=RDA1,y=RDA2,label=rownames(df2),
#   #                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))),
#   #            color="red", size=3) +
#   labs(title = "RDA of Allele Frequencies", subtitle = "Site: Rosemount",
#        x = paste("RDA1 (", round(variance_Rosemount[1], 2),"%)"),
#        y = paste("RDA1 (", round(variance_Rosemount[2], 2),"%)"))

# rda.biplot.Rosemount

# # Permanova

# rda.model.Rosemount<-anova(SNP_rda_Rosemount, step=1000, perm.max=1000, by= "terms")

# # Model: rda(formula = RDA_input_Rosemount ~ meta_Rosemount$Ecotype + meta_Rosemount$Temp)
# #                         Df   Variance F       Pr(>F)  
# # meta_Rosemount$Ecotype  3    20142    0.9990  0.402  
# # meta_Rosemount$Temp     1    18720    2.7855  0.025 *
# # Residual                3    20162 

# summary_rda_Rosemount <- data.frame(Term=row.names(rda.model.Rosemount),
#                                  Df=rda.model.Rosemount$Df,
#                                  Prop.Var = round(rda.model.Rosemount$Variance/sum(rda.model.Rosemount$Variance),3),
#                                  Fstat=round(rda.model.Rosemount$F,2),
#                                  Pvalue=round(rda.model.Rosemount$`Pr(>F)`,3),
#                                  Radj=as.numeric(RsquareAdj(SNP_rda_Rosemount))[2])

# # Big figure
# rda.plot.Roseau + rda.plot.Rosemount + rda.plot.StPaul



