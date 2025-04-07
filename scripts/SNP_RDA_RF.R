# library(vcfR)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
require(tidyr)
library(tidyverse)
library(wesanderson) # color schemes
library(patchwork) # display two figures at once

data <- read.delim("snp_freq_diffs_test_mincov2_rc", sep = "\t", header = T)

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
result_df <- convert_and_replace(data, start_column_index)

# Display the result
print(result_df)

#### ALL Replicons ####
# drop rows with NA
snp_data <- data.frame(result_df) 
head(snp_data)
snp_data <- snp_data %>% drop_na() 

# Add new names
colnames(snp_data)[1] = "Chr"
colnames(snp_data)[2] = "Pos"
colnames(snp_data)[3] = "rc"
colnames(snp_data)[4] = "allele_count"
colnames(snp_data)[5] = "allele_states"
colnames(snp_data)[6] = "deletion_sum"
colnames(snp_data)[7] = "snp_type"
colnames(snp_data)[8] = "major_alleles(maa)"
colnames(snp_data)[9] = "minor_alleles(mia)"
colnames(snp_data)[10] = "AUMerit_Roseau_Cold"
colnames(snp_data)[11] = "AUMerit_Roseau_Warm"
colnames(snp_data)[12] = "AUMerit_Rosemount_Cold"
colnames(snp_data)[13] = "AUMerit_Rosemount_Warm"
colnames(snp_data)[14] = "AUMerit_StPaul_Cold"
colnames(snp_data)[15] = "AUMerit_StPaul_Warm"
colnames(snp_data)[16] = "Hungvillosa_Roseau_Cold"
colnames(snp_data)[17] = "Hungvillosa_Roseau_Warm"
colnames(snp_data)[18] = "Hungvillosa_Rosemount_Cold"
colnames(snp_data)[19] = "Hungvillosa_Rosemount_Warm"
colnames(snp_data)[20] = "Hungvillosa_StPaul_Cold"
colnames(snp_data)[21] = "Hungvillosa_StPaul_Warm"
colnames(snp_data)[22] = "MSP4045_Roseau_Cold"
colnames(snp_data)[23] = "MSP4045_Roseau_Warm"
colnames(snp_data)[24] = "MSP4045_Rosemount_Cold"
colnames(snp_data)[25] = "MSP4045_Rosemount_Warm"
colnames(snp_data)[26] = "MSP4045_StPaul_Cold"
colnames(snp_data)[27] = "MSP4045_StPaul_Warm"
colnames(snp_data)[28] = "PupleBounty_StPaul_Warm"
colnames(snp_data)[29] = "PupleBounty_Roseau_Cold"
colnames(snp_data)[30] = "PupleBounty_Roseau_Warm"
colnames(snp_data)[31] = "PupleBounty_Rosemount_Cold"
colnames(snp_data)[32] = "PupleBounty_Rosemount_Warm"
colnames(snp_data)[33] = "PupleBounty_StPaul_Cold"
colnames(snp_data)[34] = "mia_AUMerit_Roseau_Cold"
colnames(snp_data)[35] = "mia_AUMerit_Roseau_Warm"
colnames(snp_data)[36] = "mia_AUMerit_Rosemount_Cold"
colnames(snp_data)[37] = "mia_AUMerit_Rosemount_Warm"
colnames(snp_data)[38] = "mia_AUMerit_StPaul_Cold"
colnames(snp_data)[39] = "mia_AUMerit_StPaul_Warm"
colnames(snp_data)[40] = "mia_Hungvillosa_Roseau_Cold"
colnames(snp_data)[41] = "mia_Hungvillosa_Roseau_Warm"
colnames(snp_data)[42] = "mia_Hungvillosa_Rosemount_Cold"
colnames(snp_data)[43] = "mia_Hungvillosa_Rosemount_Warm"
colnames(snp_data)[44] = "mia_Hungvillosa_StPaul_Cold"
colnames(snp_data)[45] = "mia_Hungvillosa_StPaul_Warm"
colnames(snp_data)[46] = "mia_MSP4045_Roseau_Cold"
colnames(snp_data)[47] = "mia_MSP4045_Roseau_Warm"
colnames(snp_data)[48] = "mia_MSP4045_Rosemount_Cold"
colnames(snp_data)[49] = "mia_MSP4045_Rosemount_Warm"
colnames(snp_data)[50] = "mia_MSP4045_StPaul_Cold"
colnames(snp_data)[51] = "mia_MSP4045_StPaul_Warm"
colnames(snp_data)[52] = "mia_PupleBounty_StPaul_Warm"
colnames(snp_data)[53] = "mia_PupleBounty_Roseau_Cold"
colnames(snp_data)[54] = "mia_PupleBounty_Roseau_Warm"
colnames(snp_data)[55] = "mia_PupleBounty_Rosemount_Cold"
colnames(snp_data)[56] = "mia_PupleBounty_Rosemount_Warm"
colnames(snp_data)[57] = "mia_PupleBounty_StPaul_Cold"



small_snp_data <- snp_data[,c(1:2,10:33)]
small_snp_data$SNP <- paste(small_snp_data$Chr,small_snp_data$Pos,sep="_")
head(small_snp_data)

# Make last row the first row
small_snp_data <- small_snp_data %>%
  relocate(SNP, before="Chr")


# Get rid of columns I don't need
small_snp_data <- small_snp_data[,-c(2:3)]
head(small_snp_data)

# Make column 1 the row names

RDA_input <- small_snp_data[,-1]
rownames(RDA_input) <- small_snp_data[,1]

RDA_input <- t(RDA_input)
RDA_input <- as.data.frame(RDA_input)

# extract and define metadata from snp matrix row names
meta <- data.frame("ID" = rownames(RDA_input))
meta <- separate(meta,col=ID, into = c("Ecotype","Site","Temp"), sep = "_")
meta$ID <- rownames(RDA_input)
meta <- as.data.frame(meta)

SNP_rda <- rda(RDA_input ~ meta$Ecotype + meta$Site + meta$Temp)
smry <- summary(SNP_rda)
all_rda <- SNP_rda
variance <- SNP_rda$CA$eig/SNP_rda$tot.chi*100

#permanova
rda.model.all<-anova(all_rda, step=1000, perm.max=1000, by= "terms")
print(rda.model.all)

summary_rda_all <- data.frame(Term=row.names(rda.model.all),
                                 Df=rda.model.all$Df,
                                 Prop.Var = round(rda.model.all$Variance/sum(rda.model.all$Variance),3),
                                 Fstat=round(rda.model.all$F,2),
                                 Pvalue=round(rda.model.all$`Pr(>F)`,3),
                                 Radj=as.numeric(RsquareAdj(SNP_rda))[2])

print(summary_rda_all)


##################### Roseau ####################################

RDA_input %>%
  filter(grepl("Roseau",rownames(RDA_input))) -> RDA_input_Roseau

rownames(RDA_input_Roseau)

meta_Roseau <- meta %>%
  filter(grepl("Roseau",ID))

# Run RDA
SNP_rda_Roseau <- rda(RDA_input_Roseau ~ meta_Roseau$Ecotype + meta_Roseau$Temp)


smry_Roseau <- summary(SNP_rda_Roseau)
variance_Roseau <- SNP_rda_Roseau$CA$eig/SNP_rda_Roseau$tot.chi*100

####permanova
rda.model.Roseau<-anova(SNP_rda_Roseau, step=1000, perm.max=1000, by= "terms")
print(rda.model.Roseau)



summary_rda_Roseau <- data.frame(Term=row.names(rda.model.Roseau),
                              Df=rda.model.Roseau$Df,
                              Prop.Var = round(rda.model.Roseau$Variance/sum(rda.model.Roseau$Variance),3),
                              Fstat=round(rda.model.Roseau$F,2),
                              Pvalue=round(rda.model.Roseau$`Pr(>F)`,3),
                              Radj=as.numeric(RsquareAdj(SNP_rda_Roseau))[2])
print(summary_rda_Roseau)


# Get ready to plot
df1_Roseau  <- data.frame(smry_Roseau$sites[,1:2])       # PC1 and PC2
df2_Roseau  <- data.frame(smry_Roseau$biplot[,1:2])   # loadings for PC1 and PC2
df1_Roseau$Ecotype <- meta_Roseau$Ecotype
df1_Roseau$Temp <- meta_Roseau$Temp

rda.plot.Roseau <- ggplot(df1_Roseau, aes(x=RDA1, y=RDA2), group = Temp) +
  geom_point(aes(color = Temp),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_color_manual(values=c("lightskyblue2", "salmon1")) +
  labs(title = "Site: Roseau", x = paste("RDA1 (", round(variance_Roseau[1], 2),"%)"),
       y = paste("RDA2 (", round(variance_Roseau[2], 2),"%)")) +
  theme_bw()


# St Paul ####
RDA_input_StPaul <- RDA_input %>%
  filter(grepl("StPaul",rownames(RDA_input)))

rownames(RDA_input_StPaul)


meta_StPaul <- meta %>%
  filter(grepl("StPaul",ID))

# Run RDA
SNP_rda_StPaul <- rda(RDA_input_StPaul ~ meta_StPaul$Ecotype + meta_StPaul$Temp)


smry_StPaul <- summary(SNP_rda_StPaul)
variance_StPaul <- SNP_rda_StPaul$CA$eig/SNP_rda_StPaul$tot.chi*100




# Permanova

rda.model.StPaul<-anova(SNP_rda_StPaul, step=1000, perm.max=1000, by= "terms")
print(rda.model.StPaul)

###Summary
summary_rda_StPaul <- data.frame(Term=row.names(rda.model.StPaul),
                                 Df=rda.model.StPaul$Df,
                                 Prop.Var = round(rda.model.StPaul$Variance/sum(rda.model.StPaul$Variance),3),
                                 Fstat=round(rda.model.StPaul$F,2),
                                 Pvalue=round(rda.model.StPaul$`Pr(>F)`,3),
                                 Radj=as.numeric(RsquareAdj(SNP_rda_StPaul))[2])
print(summary_rda_StPaul)


# Get ready to plot
df1_StPaul  <- data.frame(smry_StPaul$sites[,1:2])       # PC1 and PC2
df2_StPaul  <- data.frame(smry_StPaul$biplot[,1:2])   # loadings for PC1 and PC2
df1_StPaul$Ecotype <- meta_StPaul$Ecotype
df1_StPaul$Temp <- meta_StPaul$Temp

rda.plot.StPaul <- ggplot(df1_StPaul, aes(x=RDA1, y=RDA2), group = Temp) +
  geom_point(aes(color = Temp),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_color_manual(values=c("lightskyblue2", "salmon1")) +
  labs(title = "Site: St. Paul", x = paste("RDA1 (", round(variance_StPaul[1], 2),"%)"),
       y = paste("RDA2 (", round(variance_StPaul[2], 2),"%)")) +
  theme_bw()



# Rosemount ####
RDA_input_Rosemount <- RDA_input %>%
  filter(grepl("Rosemount",rownames(RDA_input)))

rownames(RDA_input_Rosemount)
# Success!!
# Now for meta:

meta_Rosemount <- meta %>%
  filter(grepl("Rosemount",ID))

# Run RDA
SNP_rda_Rosemount <- rda(RDA_input_Rosemount ~ meta_Rosemount$Ecotype + meta_Rosemount$Temp)

smry_Rosemount <- summary(SNP_rda_Rosemount)
variance_Rosemount <- SNP_rda_Rosemount$CA$eig/SNP_rda_Rosemount$tot.chi*100

# Permanova

rda.model.Rosemount<-anova(SNP_rda_Rosemount, step=1000, perm.max=1000, by= "terms")

print(rda.model.Rosemount)
##Summary
summary_rda_Rosemount <- data.frame(Term=row.names(rda.model.Rosemount),
                                 Df=rda.model.Rosemount$Df,
                                 Prop.Var = round(rda.model.Rosemount$Variance/sum(rda.model.Rosemount$Variance),3),
                                 Fstat=round(rda.model.Rosemount$F,2),
                                 Pvalue=round(rda.model.Rosemount$`Pr(>F)`,3),
                                 Radj=as.numeric(RsquareAdj(SNP_rda_Rosemount))[2])
print(summary_rda_Rosemount)

# Get ready to plot
df1_Rosemount  <- data.frame(smry_Rosemount$sites[,1:2])       # PC1 and PC2
df2_Rosemount  <- data.frame(smry_Rosemount$biplot[,1:2])   # loadings for PC1 and PC2
df1_Rosemount$Ecotype <- meta_Rosemount$Ecotype
df1_Rosemount$Temp <- meta_Rosemount$Temp

rda.plot.Rosemount <- ggplot(df1_Rosemount, aes(x=RDA1, y=RDA2), group = Temp) +
  geom_point(aes(color = Temp),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_color_manual(values=c("lightskyblue2", "salmon1")) +
  labs(title = "Site: Rosemount", x = paste("RDA1 (", round(variance_Rosemount[1], 2),"%)"),
       y = paste("RDA2 (", round(variance_Rosemount[2], 2),"%)")) +
  theme_bw()



# Big figure
rda.plot.Roseau
rda.plot.StPaul
rda.plot.Rosemount
rda.plot.Roseau + rda.plot.Rosemount + rda.plot.StPaul
print(summary_rda_all)