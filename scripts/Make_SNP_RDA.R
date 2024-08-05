# library(vcfR)
#confirming which packages are actually used 
library(vegan)
library(ggplot2)
library(tidyverse)
require(tidyr) #for separate()
require(dplyr) #for relocate()


library(ggpubr)
library(tidyverse)
library(wesanderson) # color schemes
library(patchwork) # display two figures at once


#R version being used in documentation-env: 4.3.3
setwd("C:/Users/paula/Desktop/TrappingCh2/rhizobium_trapping/data")
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


############################ENSURE HEADER NAMES ARE POPULATED VIA A BAM FILE .TXT THAT POPOOLATION, SAMTOOLS WORK FROM#########################
lines <- readLines("sorted_bam_files.txt")

# Process each line
sample_names <- sub("^.*bamfiles/", "", lines)  # Remove everything before and including "bamfiles/"
sample_names <- sub("\\.sorted$", "", sample_names)  # Remove ".sorted"

# Print the processed lines to check the result
print("df of sample names from .bams made")
sample_names_df <- data.frame(sample_names)

# suffix to add
suffix <- "_majoralleles"

# Use a for loop to rename the major allele columns in the snp data
for (i in 1:nrow(sample_names_df)) {
  snp_data_col_index <- i + 9  # Adjust index to start renaming from column 10
  colnames(snp_data)[snp_data_col_index] <- paste0( sample_names_df$sample_names[i], suffix)
}




###########################################       RDA ANALYSIS ############################################

# #### ALL Replicons ####
# drop rows with NA
snp_data <- snp_data %>% drop_na() 

# minor alleles not used in analysis


#drop the columns between 'pos' and the major allele data, so we only have chromosome, position, and
#snp frequency
small_snp_data <- snp_data[,c(1:2,10:33)]

####make a new column combining the chromosome and position information
small_snp_data$SNP <- paste(small_snp_data$X..chr, small_snp_data$pos, sep="_")
# now move it to beginning
small_snp_data <- small_snp_data %>% relocate(SNP)
# Get rid of the columns we referenced to make that row
small_snp_data <- small_snp_data[,-c(2:3)]

#now, make that column into rownames (we could do this with tidyverse, too, but here's the native way)
RDA_input <- small_snp_data[,-1]
rownames(RDA_input) <- small_snp_data[,1]
#result: a df that has all major allele proportion values for all SNPS for all replicates 
RDA_input <- t(RDA_input)
RDA_input <- as.data.frame(RDA_input)



# extract and define metadata from snp matrix row names
meta <- data.frame("ID" = rownames(RDA_input))
meta <- separate(meta,col=ID, into = c("Ecotype","Site","Temp"), sep = "_", extra="drop") 
#extra: removes the suffix added in via the for() loop

meta$ID <- rownames(RDA_input)
meta <- as.data.frame(meta)
print(head(meta))
 



# print("data to perform RDA made")
SNP_rda <- rda(RDA_input ~ meta$Ecotype + meta$Site + meta$Temp)
smry <- summary(SNP_rda)
all_rda <- SNP_rda
all_variance <- SNP_rda$CA$eig/SNP_rda$tot.chi*100

# all_varp <- varpart(RDA_input, ~ meta$Ecotype, ~ meta$Site, ~ meta$Temp) # Var part not relevant

#SAMTOOLS 1.10 (earlier)
# Summary:
# Partitioning of variance:
#                 Inertia   Proportion
# Total           40411     1.0000
# Constrained     13105     0.3243
# Unconstrained   27307     0.6757

# # RDA1 explains 0.5347, RDA2 explains 0.1781, RDA3 explains 0.1369

#SAMTOOLS 1.16(july)
#Partitioning of variance:
#               Inertia Proportion
#Total           39.05     1.0000
#Constrained     11.74     0.3006
#Unconstrained   27.31     0.6994

#RDA 1,2,3
#Proportion Explained  0.1506 0.06925 0.03705


# # What if we take out ecotype, as Liana suggested?
SNP_rda_no_ecotype <- rda(RDA_input ~ meta$Site + meta$Temp)
smry_no_ecotype <- summary(SNP_rda_no_ecotype)
rda.model.no.ecotype<-anova(SNP_rda_no_ecotype, step=1000, perm.max=1000, by= "terms")

#SAMTOOLS 1.10
# # Summary:
# # Partitioning of variance:
# #                 Inertia Proportion
# # Total           40411     1.0000
# # Constrained      9137     0.2261
# # Unconstrained   31275     0.7739

# # Not sure if this is better or worse

##SAMTOOLS 1.16
#Partitioning of variance:
#  Inertia Proportion
#Total          39.048     1.0000
#Constrained     8.441     0.2162
#Unconstrained  30.607     0.7838


# # What if we just look at site?
SNP_rda_only_site <- rda(RDA_input ~ meta$Site)
rda.model.only.site <-anova(SNP_rda_only_site, step=1000, perm.max=1000, by= "terms")

summary_rda_only_site <- data.frame(Term=row.names(rda.model.only.site),
                                      Df=rda.model.only.site$Df,
                                      Prop.Var = round(rda.model.only.site$Variance/sum(rda.model.only.site$Variance),3),
                                      Fstat=round(rda.model.only.site$F,2),
                                      Pvalue=round(rda.model.only.site$`Pr(>F)`,3),
                                      Radj=as.numeric(RsquareAdj(SNP_rda_only_site))[2])

df1  <- data.frame(smry$sites[,1:2])    # PC1 and PC2
df2  <- data.frame(smry$biplot[,1:2])   # loadings for PC1 and PC2
df1$Ecotype <- meta$Ecotype
df1$Site <- meta$Site
df1$Temp <- meta$Temp

rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Site) +
   geom_point(aes(color = Site, shape=Temp),size=2.5) +
   geom_hline(yintercept=0, linetype="dotted") +
   geom_vline(xintercept=0, linetype="dotted") +
   scale_color_manual(values = wes_palette("Darjeeling1"), name = "Sites") +
   labs(title = "RDA of Allele Frequencies", subtitle = "All sites", x = paste("RDA1 (", round(all_variance[1], 2),"%)"),
        y = paste("RDA1 (", round(all_variance[2], 2),"%)")) 
   
rda.plot

rownames(df2) <- gsub("sam\\$", "", rownames(df2))

#############do we want/need this? 
rda.temp <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Temp) +
   geom_point(aes(color = Temp)) +
   geom_hline(yintercept=0, linetype="dotted") +
   geom_vline(xintercept=0, linetype="dotted") +
   scale_color_manual(values = wes_palette("Darjeeling1"), name = "Temp")

rownames(df2) <- gsub("sam\\$", "", rownames(df2))


# Run with anova ####
all_results <- anova(all_rda)

# Code from Liana script
rda.model.all<-anova(all_rda, step=1000, perm.max=1000, by= "terms")

summary_rda_all <- data.frame(Term=row.names(rda.model.all),
                                 Df=rda.model.all$Df,
                                 Prop.Var = round(rda.model.all$Variance/sum(rda.model.all$Variance),3),
                                 Fstat=round(rda.model.all$F,2),
                                 Pvalue=round(rda.model.all$`Pr(>F)`,3),
                                 Radj=as.numeric(RsquareAdj(SNP_rda))[2])

summary_rda_no_ecotype <- data.frame(Term=row.names(rda.model.no.ecotype),
                              Df=rda.model.no.ecotype$Df,
                              Prop.Var = round(rda.model.no.ecotype$Variance/sum(rda.model.no.ecotype$Variance),3),
                              Fstat=round(rda.model.no.ecotype$F,2),
                              Pvalue=round(rda.model.no.ecotype$`Pr(>F)`,3),
                              Radj=as.numeric(RsquareAdj(SNP_rda_no_ecotype))[2])

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



############# Roseau ####

RDA_input %>%
  filter(grepl("Roseau",rownames(RDA_input))) -> RDA_input_Roseau

rownames(RDA_input_Roseau)
# Success!!
# Now for meta:

meta_Roseau <- meta %>%
  filter(grepl("Roseau",ID))

# Run RDA
SNP_rda_Roseau <- rda(RDA_input_Roseau ~ meta_Roseau$Ecotype + meta_Roseau$Temp)
# SNP_rda_Roseau <- rda(RDA_input_Roseau ~ meta_Roseau$Temp)

# Anova
rda.model.Roseau<-anova(SNP_rda_Roseau, step=1000, perm.max=1000, by= "terms")


# Model: rda(formula = RDA_input_Roseau ~ meta_Roseau$Ecotype + meta_Roseau$Temp)
#                      Df  Variance F       Pr(>F)
# meta_Roseau$Ecotype  3   6265.4   1.3029  0.116
# meta_Roseau$Temp     1   2698.6   1.6835  0.071 .
# Residual             3   4809.0

# Temp approaching significance: p = 0.071

# LESS significant when I take out ecotype as a term...
# Model: rda(formula = RDA_input_Roseau ~ meta_Roseau$Temp)
#                   Df  Variance  F       Pr(>F)
# meta_Roseau$Temp  1   2698.6    1.4621  0.129
# Residual          6   11074.4


smry_Roseau <- summary(SNP_rda_Roseau)
variance_Roseau <- SNP_rda_Roseau$CA$eig/SNP_rda_Roseau$tot.chi*100
#varp_Roseau <- varpart(RDA_input_Roseau, ~ meta_Roseau$Ecotype, ~ meta_Roseau$Temp)
summary_rda_Roseau <- data.frame(Term=row.names(rda.model.Roseau),
                              Df=rda.model.Roseau$Df,
                              Prop.Var = round(rda.model.Roseau$Variance/sum(rda.model.Roseau$Variance),3),
                              Fstat=round(rda.model.Roseau$F,2),
                              Pvalue=round(rda.model.Roseau$`Pr(>F)`,3),
                              Radj=as.numeric(RsquareAdj(SNP_rda_Roseau))[2])

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
       y = paste("RDA1 (", round(variance_Roseau[2], 2),"%)")) +
  theme_bw()

rda.plot.Roseau


# # St Paul ####
RDA_input_StPaul <- RDA_input %>%
  filter(grepl("StPaul",rownames(RDA_input)))

rownames(RDA_input_StPaul)
# Success!!
# Now for meta:

meta_StPaul <- meta %>%
  filter(grepl("StPaul",ID))

#Run RDA
SNP_rda_StPaul <- rda(RDA_input_StPaul ~ meta_StPaul$Ecotype + meta_StPaul$Temp)
# SNP_rda_StPaul <- rda(RDA_input_StPaul ~ meta_StPaul$Temp)

smry_StPaul <- summary(SNP_rda_StPaul)
variance_StPaul <- SNP_rda_StPaul$CA$eig/SNP_rda_StPaul$tot.chi*100
#varp_StPaul <- varpart(RDA_input_StPaul, ~ meta_StPaul$Ecotype, ~ meta_StPaul$Temp)

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
       y = paste("RDA1 (", round(variance_StPaul[2], 2),"%)")) +
  theme_bw()


# # Permanova

rda.model.StPaul<-anova(SNP_rda_StPaul, step=1000, perm.max=1000, by= "terms")

# Model: rda(formula = RDA_input_StPaul ~ meta_StPaul$Ecotype + meta_StPaul$Temp)
#                      Df Variance F       Pr(>F)
# meta_StPaul$Ecotype  3  14349.0  0.8217  0.579
# meta_StPaul$Temp     1   2505.9  0.4305  0.694
# Residual             3  17462.0

# p even HIGHER when we take out ecotype term
# Model: rda(formula = RDA_input_StPaul ~ meta_StPaul$Temp)
#                   Df    Variance  F       Pr(>F)
# meta_StPaul$Temp  1     2506      0.4726  0.702
# Residual          6     31811

summary_rda_StPaul <- data.frame(Term=row.names(rda.model.StPaul),
                                 Df=rda.model.StPaul$Df,
                                 Prop.Var = round(rda.model.StPaul$Variance/sum(rda.model.StPaul$Variance),3),
                                 Fstat=round(rda.model.StPaul$F,2),
                                 Pvalue=round(rda.model.StPaul$`Pr(>F)`,3),
                                 Radj=as.numeric(RsquareAdj(SNP_rda_StPaul))[2])



# # Rosemount ####
RDA_input_Rosemount <- RDA_input %>%
  filter(grepl("Rosemount",rownames(RDA_input)))

rownames(RDA_input_Rosemount)
# Success!!
# Now for meta:

meta_Rosemount <- meta %>%
  filter(grepl("Rosemount",ID))

# Run RDA
SNP_rda_Rosemount <- rda(RDA_input_Rosemount ~ meta_Rosemount$Ecotype + meta_Rosemount$Temp)
# SNP_rda_Rosemount <- rda(RDA_input_Rosemount ~ meta_Rosemount$Temp)

smry_Rosemount <- summary(SNP_rda_Rosemount)
variance_Rosemount <- SNP_rda_Rosemount$CA$eig/SNP_rda_Rosemount$tot.chi*100
varp_Rosemount <- varpart(RDA_input_Rosemount, ~ meta_Rosemount$Ecotype, ~ meta_Rosemount$Temp)

# # Get ready to plot
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
       y = paste("RDA1 (", round(variance_Rosemount[2], 2),"%)")) +
  theme_bw()


# # Permanova
rda.model.Rosemount<-anova(SNP_rda_Rosemount, step=1000, perm.max=1000, by= "terms")

# # Model: rda(formula = RDA_input_Rosemount ~ meta_Rosemount$Ecotype + meta_Rosemount$Temp)
# #                         Df   Variance F       Pr(>F)
# # meta_Rosemount$Ecotype  3    20142    0.9990  0.402
# # meta_Rosemount$Temp     1    18720    2.7855  0.025 *
# # Residual                3    20162

summary_rda_Rosemount <- data.frame(Term=row.names(rda.model.Rosemount),
                                 Df=rda.model.Rosemount$Df,
                                 Prop.Var = round(rda.model.Rosemount$Variance/sum(rda.model.Rosemount$Variance),3),
                                 Fstat=round(rda.model.Rosemount$F,2),
                                 Pvalue=round(rda.model.Rosemount$`Pr(>F)`,3),
                                 Radj=as.numeric(RsquareAdj(SNP_rda_Rosemount))[2])

# # Big figure
rda.plot.Roseau + rda.plot.Rosemount + rda.plot.StPaul








#########################write functions to do all this plotting ################

# Define your functions
metadata_generator <- function(site_name) {
  metadata = meta %>% filter(grepl(site_name,ID))
}

RDA_input_generator <- function(site_name) {
  RDA_input %>% filter(grepl(site_name,rownames(RDA_input))) 
}

tasks <- list(
  list(func = metadata_generator, args = list(site_name = "Roseau")), #SNP RDA for everything 
  list(func = metadata_generator, args = list(site_name = "Rosemount")),
  list(func = metadata_generator, args = list(site_name = "StPaul")),
  list(func = RDA_input_generator, args = list(site_name = "Roseau")), #SNP RDA for everything
  list(func = RDA_input_generator, args = list(site_name = "Rosemount")),
  list(func = RDA_input_generator, args = list(site_name = "StPaul"))
  ) 
# Function to call a given function with arguments from a list
call_with_args <- function(tasks) {
  do.call(tasks$func, tasks$args)
}

# Iterate over functions and argument sets and call the functions
results <- lapply(tasks, call_with_args)



RDAs <- function(RDA_input, metadata) {
  #what the function does
  rda(RDA_input ~ sample_metadata$Ecotype + sample_metdata$Temp)
  return(result)
}









########################break down analyses by site #############

#make blank lists to store our variables to reference when doing the RDAs
RDA_inputs_list <- list()
metadata_list <- list()
#create inputs and metadata to do RDA analyses with
site_names <- c("", "Rosemount", "Roseau", "StPaul") #having this be blank is how I'm using the unfiltered data

for (site_name in site_names) {
  if (site_name == "") {site_name <- "AllData"} #this gives a name to the unfiltered data
  
  if (site_name == "AllData") { 
  RDA_inputs_list[[site_name]] <- RDA_input
  metadata_list[[site_name]] <- meta} #this is handling the unfiltered case
  
  else {
  RDA_inputs_list[[site_name]] <- RDA_input %>% filter(grepl(site_name,rownames(RDA_input)))
  metadata_list[[site_name]] <- meta %>% filter(grepl(site_name, ID))
  }
}

#the actual RDA analysis
RDA_analyses_list <- list()

for (i in seq_along(metadata_list)) {
  metadata_df <- metadata_list[[i]] #seq_along and df renaming *may* not be necessary]
  data_input <- RDA_inputs_list[[i]]
  metadata_df_name <- names(metadata_list)[i]
  
  if (length(unique(metadata_df$Site)) > 1) {
    # print(head(RDA_inputs_list[[i]]))
    RDA_analyses_list[[metadata_df_name]] <- rda(data_input ~ metadata_df$Ecotype 
                                                 + metadata_df$Site 
                                                 + metadata_df$Temp)

    print("analysis")
  } 
  else {
    RDA_analyses_list[[metadata_df_name]] <- rda(data_input ~ metadata_df$Ecotype
                                                 + metadata_df$Temp) #no site, as there is only
    #one site in these samples (see: if statement)
    print("outside loop")
    
  
  }
  
}
  

###PLOTS
summaries_list <- list()
variances_list <- list()
RDA_plots <- list()
# df_list <- list()


for (i in seq_along(RDA_analyses_list)) {
  item <- RDA_analyses_list[[i]]
  item_name <- names(RDA_analyses_list[i])
  summaries_list[[item_name]] <- summary(item)
  variances_list[[item_name]] <- item$CA$eig/item$tot.chi*100
  df1 <- data.frame(summaries_list[[item_name]]$sites[,1:2])
  df1$Ecotype <- metadata_list[[item_name]]$Ecotype
  df1$Site <- metadata_list[[item_name]]$Site
  df1$Temp <- metadata_list[[item_name]]$Temp
  

  if (length(unique(metadata_list[[item_name]]$Site)) >1) {
    RDA_plots[[item_name]] <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Site) +
      geom_point(aes(color = Site, shape=Temp),size=2.5) +
      geom_hline(yintercept=0, linetype = 'dotted') + 
      geom_vline(xintercept=0, linetype = 'dotted') +
      scale_color_manual(values = wes_palette("Darjeeling1"), name = "Sites") +
      labs(title = "RDA of Allele Frequencies", subtitle = "All sites", 
           x= paste("RDA1", round(variances_list[[item_name]][1], 2),"%"),
           y= paste("RDA2", round(variances_list[[item_name]][2], 2), "%"))
  }
  else {
    RDA_plots[[item_name]] <- ggplot(df1, aes(x=RDA1, y=RDA2), group = Temp) +
      geom_point(aes(color = Temp),size=3) +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      scale_color_manual(values=c("lightskyblue2", "salmon1")) +
      labs(title = paste("Site:",item_name), x = paste("RDA1 (", round(variances_list[[item_name]][1], 2),"%)"),
           y = paste("RDA1 (", round(variances_list[[item_name]][2], 2),"%)")) +
      theme_bw()
  }
  
}



