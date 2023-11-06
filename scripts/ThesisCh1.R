library(tidyverse)
library(agricolae)
library(wesanderson)
library(poolfstat)
library(emmeans)
library(FactoMineR)
library(factoextra)
require("psych")
require("data.table")

trapping <- read.csv("TrappingData.csv")

fstmatrix <- read.csv("FSTmatrix.csv")

phen.data <- read.csv("PhenotypeMatrix.csv")

# PLANT DATA ####

harvest <- trapping %>%
  filter(TimePoint =="Harvest")

harvest$efficiency = harvest$BiomassDryWt/harvest$NodDryWt

warm.harvest <- harvest %>%
  filter(Temp =="Warm")
  
cold.harvest <- harvest %>%
  filter(Temp =="Cold")

# Growth data, across all soils
warmGrowth <- trapping %>%
  filter(Temp=="Warm") %>%
  filter(TimePoint %in% c("Week2","Week3","Week4","Week5","Week6","Week7","Week8")) %>%
  group_by(Ecotype,TimePoint) %>%
  summarise(avg = mean(LeafArea,na.rm=TRUE),sd=sd(LeafArea)) %>%
  arrange(avg)
warmGrowth %>%
  ggplot(aes(x=TimePoint,y=avg,group=Ecotype,color=Ecotype)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = "Time Point",y="Average leaf area (cm^2)")

coldGrowth <- trapping %>%
  filter(Temp=="Cold") %>%
  filter(TimePoint %in% c("Week2","Week3","Week4","Week5","Week6","Week7","Week8")) %>%
  group_by(Ecotype,TimePoint) %>%
  summarise(avg = mean(LeafArea),sd=sd(LeafArea)) %>%
  arrange(avg)
coldGrowth %>%
  ggplot(aes(x=TimePoint,y=avg,group=Ecotype,color=Ecotype)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = "Time Point",y="Average leaf area (cm^2)")


### Water uptake vs. leaf area ###
Week3 <- trapping %>%
  filter(TimePoint=="Week3")
plot(WaterUptake ~ LeafArea, data=Week3)
cor(Week3$WaterUptake,Week3$LeafArea)
# R2 = 0.66

Week4 <- trapping %>%
  filter(TimePoint=="Week4")
plot(WaterUptake ~ LeafArea, data=Week4)
cor(Week4$WaterUptake,Week4$LeafArea)
# R2 = 0.88

Week5 <- trapping %>%
  filter(TimePoint=="Week5")
plot(WaterUptake ~ LeafArea, data=Week5)
cor(Week5$WaterUptake,Week5$LeafArea,use="complete.obs")
# R2 = 0.93

Week6 <- trapping %>%
  filter(TimePoint=="Week6")
plot(WaterUptake ~ LeafArea, data=Week6)
cor(Week6$WaterUptake,Week6$LeafArea)
# R2 = 0.95

Week7 <- trapping %>%
  filter(TimePoint=="Week7")
plot(WaterUptake ~ LeafArea, data=Week7)
cor(Week6$WaterUptake,Week6$LeafArea)
# R2 = 0.95


# Same graphs, water uptake
waterUptakeCold <- trapping %>%
  filter(Temp=="Cold") %>%
  filter(TimePoint %in% c("Week4","Week5","Week6","Week7","Week8")) %>%
  group_by(Ecotype,TimePoint) %>%
  summarise(avg = mean(WaterUptake,na.rm=TRUE),sd=sd(WaterUptake,na.rm=TRUE)) %>%
  arrange(avg)

waterUptakeCold %>%
  ggplot(aes(x=TimePoint,y=avg,group=Ecotype,color=Ecotype)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = "Time Point",y="Water use (g)")

waterUptakeWarm <- trapping %>%
  filter(Temp=="Warm") %>%
  filter(TimePoint %in% c("Week3","Week4","Week5","Week6","Week7","Week8")) %>%
  group_by(Ecotype,TimePoint) %>%
  summarise(avg = mean(WaterUptake,na.rm=TRUE),sd=sd(WaterUptake,na.rm=TRUE)) %>%
  arrange(avg)

waterUptakeWarm %>%
  ggplot(aes(x=TimePoint,y=avg,group=Ecotype,color=Ecotype)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = "Time Point",y="Water use (g)")

# HARVEST ####

# Biomass:
harvest %>%
  ggplot(aes(x=Ecotype,y=BiomassDryWt,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  labs(x = "Ecotype",y="Biomass Dry Weight (g)")

aov1 <- aov(BiomassDryWt ~ Ecotype + Soil + Temp,data=harvest)
summary(aov1)
# Ecotype effect: p = 0.001
# Soil effect: p = 0.6 (no effect)
# Temp effect: p < 2e -16

aov1a <- aov(BiomassDryWt ~ Ecotype*Temp + Soil,data=harvest)
summary(aov1a)
# Very weak interaction effect. p = 0.055

lm.cold <- aov(BiomassDryWt ~ Ecotype + Soil,data=cold.harvest)
summary(lm.cold)
# Ecotype effect: p = 3.43e-05
# Soil effect: p = 0.84

cold.means <- HSD.test(lm.cold, "Ecotype")
cold.means

#                 BiomassDryWt groups
# AU Merit         0.5750000      a
# MSP 4045         0.3350000      b
# Hungvillosa      0.3225000      b
# Purple Bounty    0.3083333      b

lm.warm <- aov(BiomassDryWt ~ Ecotype + Soil,data=warm.harvest)
summary(lm.warm)
# Ecotype effect: p = 0.0112 (less of an effect than at low temps)
# Soil effect: p = 0.84

warm.means <- HSD.test(lm.warm, "Ecotype")
warm.means

#              BiomassDryWt    groups
# AU Merit          2.824167      a
# Hungvillosa       2.541667     ab
# MSP 4045          2.051667     ab
# Purple Bounty     1.833333      b

####

# Root dry weight:
harvest %>%
  ggplot(aes(x=Ecotype,y=RootDryWt,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  labs(x = "Ecotype",y="Root Dry Weight (g)")

aov2 <- aov(RootDryWt ~ Ecotype + Soil + Temp,data=harvest)
summary(aov2)
# Ecotype effect: p = 0.0297
# Soil effect: p = 0.7 (no effect)
# Temp effect: p < 1.44 e-11

lm.cold.root <- aov(RootDryWt ~ Ecotype + Soil,data=cold.harvest)
summary(lm.cold.root)
# Ecotype effect: p = 0.05 (not really an effect)
# Soil effect: p = 0.79 (no effect)

cold.means <- HSD.test(lm.cold.root, "Ecotype")
cold.means

lm.warm.root <- aov(RootDryWt ~ Ecotype + Soil,data=warm.harvest)
summary(lm.warm.root)
# Ecotype effect: p = 0.1 (less of an effect than at low temps)
# Soil effect: p = 0.65 (no effect)

warm.means <- HSD.test(lm.warm, "Ecotype")
warm.means

####
# Nod dry weight:
harvest %>%
  ggplot(aes(x=Ecotype,y=NodDryWt,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  labs(x = "Ecotype",y="Nodule Dry Weight (g)")

aov3 <- aov(NodDryWt ~ Ecotype + Soil + Temp,data=harvest)
summary(aov3)
# Ecotype effect: p < 2.05 e-08
# Soil effect: p = 0.153 (no effect)
# Temp effect: p < 2e-16

aov3a <- aov(NodDryWt ~ Ecotype*Temp + Soil,data=harvest)
summary(aov3a)
# Ecotype effect: p < 0.0001
# Soil effect: p = 0.13
# Temp effect: p < 0.0001
# Interaction: p = 0.02

lm.cold.nod <- aov(NodDryWt ~ Ecotype + Soil,data=cold.harvest)
summary(lm.cold.nod)
# Ecotype effect: p = 1.48 e-05
# Soil effect: p = 0.476

cold.means.nod <- HSD.test(lm.cold.nod, "Ecotype")
cold.means.nod

#                NodDryWt     groups
# AU Merit      0.024333333      a
# Hungvillosa   0.014416667      b
# MSP 4045      0.012916667      b
# Purple Bounty 0.009416667      b

# Note: Same groups as biomass for cold

####

lm.warm.nod <- aov(NodDryWt ~ Ecotype + Soil,data=warm.harvest)
summary(lm.warm.nod)
# Ecotype effect: p = 3e-05
# Soil effect: p = 0.135

warm.means.nod <- HSD.test(lm.warm.nod, "Ecotype")
warm.means.nod

#                NodDryWt   groups
# AU Merit      0.06325000      a
# Hungvillosa   0.04366667      b
# MSP 4045      0.04191667     bc
# Purple Bounty 0.02466667      c

# Interesting: more means separation than at cold temperatures.


# Efficiency (return on nodule construction cost)
harvest$efficiency = harvest$BiomassDryWt/harvest$NodDryWt
harvest %>%
  ggplot(aes(x=Ecotype,y=efficiency,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  labs(x = "Ecotype",y="Efficiency (Biomass g/Nodule g)")

harvest %>%
  ggplot(aes(x=Soil,y=efficiency,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  labs(x = "Soil",y="Efficiency (Biomass g/Nodule g)")

aov4 <- aov(efficiency ~ Ecotype + Soil + Temp,data=harvest)
summary(aov4)
# Ecotype effect: p = 0.005
# Soil effect: p = 0.017 (first time we see an effect!)
# Temp effect: p < 7.06 e-11

aov4a <- aov(efficiency ~ Ecotype*Temp + Soil,data=harvest)
summary(aov4a)
# No interaction between ecotype and temp

aov4b <- aov(efficiency ~ Soil*Temp + Ecotype,data=harvest)
summary(aov4b)
# Soil: p = 0.01
# Temp: p < 0.001
# Ecotype: p = 0.004
# Soil*Temp: p = 0.03

warm.harvest$efficiency = warm.harvest$BiomassDryWt/warm.harvest$NodDryWt
lm.warm.eff <- aov(efficiency ~ Ecotype + Soil,data=warm.harvest)
summary(lm.warm.eff)
# Ecotype effect: p = 0.07 (NS)
# Soil effect: p = 0.01 *

warm.means.eff <- HSD.test(lm.warm.eff, "Ecotype")
warm.means.eff

warm.means.eff <- HSD.test(lm.warm.eff, "Soil")
warm.means.eff

cold.harvest$efficiency = cold.harvest$BiomassDryWt/cold.harvest$NodDryWt
lm.cold.eff <- aov(efficiency ~ Ecotype + Soil,data=cold.harvest)
summary(lm.cold.eff)
# Ecotype effect: p = 0.002 (*)
# Soil effect: p = 0.78 (NS)

cold.means.eff <- HSD.test(lm.cold.eff, "Ecotype")
cold.means.eff

# Biomass vs. Nod Mass ####

# First look at cold vs. warm
harvest %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "Nodule Dry Weight (g)",y="Shoot dry weight (g)")

# Analysis by ecotype ####
m.interaction <- lm(NodDryWt ~ BiomassDryWt*Ecotype + Temp, data = harvest)
anova(m.interaction)

# We used linear regression to compare the relationship of nodule dry mass to shoot dry mass for each temperature
# Significant interaction of biomass dry weight and temperature to influence nod dry weight (p = 0.01)

# Obtain slopes:
m.interaction$coefficients
m.lst <- lstrends(m.interaction, "Temp", var="BiomassDryWt")
pairs(m.lst)

# There was a significant difference in the slopes between cold and warm temperatures (p = 0.01).
# But how do I list the slopes? Do I need a dummy variable to figure out the slopes? Check stats...

harvest <- as.data.table(harvest)
# Calculate Pearson's R
m.correlations <- harvest[, cor(NodDryWt,BiomassDryWt), by = Temp]
m.correlations
# Correlation in cold: 0.85
# Correlation in warm: 0.62

# Compare R values with Fisher's R to Z
paired.r(m.correlations[Temp=="Cold", V1], m.correlations[Temp=="Warm", V1], 
         n = harvest[Temp %in% c("Warm", "Cold"), .N])

# Correlation is significantly different in cold vs. warm (p < 0.001)

# Theoretically I could go through and do this by ecotype


# ==== #
# Break down by ecotype:
MSP4045 <- harvest %>%
  filter(Ecotype =="MSP 4045")

MSP4045 %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "Nodule Dry Weight (g)",y="Shoot dry weight (g)")

Hungvillosa <- harvest %>%
  filter(Ecotype =="Hungvillosa")

Hungvillosa %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "Nodule Dry Weight (g)",y="Shoot dry weight (g)")

PurpleBounty <- harvest %>%
  filter(Ecotype =="Purple Bounty")

PurpleBounty %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "Nodule Dry Weight (g)",y="Shoot dry weight (g)")


AUMerit <- harvest %>%
  filter(Ecotype =="AU Merit")

AUMerit %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Temp)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "Nodule Dry Weight (g)",y="Shoot dry weight (g)")

# Then look at four ecotypes within each temp:
warm.harvest %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Ecotype,color=Ecotype)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

warm.harvest %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Ecotype,color=Ecotype)) +
  geom_point() +
  theme_classic()

cold.harvest %>%
  ggplot(aes(x=NodDryWt,y=BiomassDryWt,shape=Ecotype,color=Ecotype)) +
  geom_point() +
  theme_classic()

####
# Root shoot ratio
harvest$rootshoot = harvest$BiomassDryWt/harvest$RootDryWt
harvest %>%
  ggplot(aes(x=Ecotype,y=rootshoot,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  ylab("Shoot Biomass / Root Biomass")

aov5 <- aov(rootshoot ~ Ecotype + Soil + Temp,data=harvest)
summary(aov5)
# Ecotype effect: p = 0.139 (NS)
# Soil effect: p = 0.819 (NS)
# Temp effect: p < 4.27 e-15

warm.harvest$rootshoot = warm.harvest$BiomassDryWt/warm.harvest$RootDryWt
lm.warm.ratio <- aov(rootshoot ~ Ecotype + Soil,data=warm.harvest)
summary(lm.warm.ratio)
# Ecotype effect: p = 0.581 (NS)
# Soil effect: p = 0.91 (NS)

cold.harvest$rootshoot = cold.harvest$BiomassDryWt/cold.harvest$RootDryWt
lm.cold.ratio <- aov(rootshoot ~ Ecotype + Soil,data=cold.harvest)
summary(lm.cold.ratio)
# Ecotype effect: p = 0.581 (NS)
# Soil effect: p = 0.91 (NS)

harvest %>%
  ggplot(aes(x=Soil,y=rootshoot,fill=Temp)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Temp) +
  scale_fill_manual(values=c("lightskyblue2", "salmon1")) +
  ylab("Shoot Biomass / Root Biomass")

# Compare variances of cold vs. warm ###
# Null hypothesis assumes ratio of variances = 1
# Low p value = reject the hypothesis that variances are equal

# Nodule dry weight:
harvest %>%
  ggplot(aes(x=Temp, y=BiomassDryWt)) + 
  geom_boxplot()
var.test(warm.harvest$BiomassDryWt,cold.harvest$BiomassDryWt)
# p <<<<< 0.01

# Nodule dry weight:
harvest %>%
  ggplot(aes(x=Temp, y=NodDryWt)) + 
  geom_boxplot()
var.test(warm.harvest$NodDryWt,cold.harvest$NodDryWt)
# p <<<<< 0.01

# Efficiency:
harvest %>%
  ggplot(aes(x=Temp, y=efficiency)) + 
  geom_boxplot()
var.test(warm.harvest$efficiency,cold.harvest$efficiency)
# p <<<<< 0.01


# SEQUENCING DATA ####

kaiju <- read.csv("KBaseNotes.csv")
names(kaiju)[6] <- "PctRhizobium" 

aov6 <- aov(PctRhizobium ~ Ecotype + Soil + Temp,data=kaiju)
summary(aov6)
# None have a significant effect, but ecotype is the closest (p = 0.09)


#####################
# Allele Data ####

library(poolfstat)

infile = 'sync_all_bams.sync'            # Name of sync file to analyze
gff_filename = 'genome.gff3'      # Name of annotation file
gene_output_filename = 'gene_tests.tsv' # Name of file for gene-level tests
snp_output_filename = 'snp_tests.tsv' # Name of file for gene-level tests
poolsizes = rep(870021, 24)     # Pool size setting

# Vector of treatments for testing allele frequency differences. MUST
# be in the same order as the columns of the sync file. Currently, this
# script is set up to test just one treatment factor.
#trts = c('cold', 'cold', 'cold', 'cold',
#         'warm', 'warm', 'warm', 'warm') 

pooldata=popsync2pooldata(sync.file="sync_all_bams.sync",poolsizes=rep(870021,24))

# Use this to compute global and per SNP FSTs:
snp.fsts <- computeFST(pooldata, method="Anova")

# Use this to compute pairwise FSTs:
pair.fst <- computePairwiseFSTmatrix(pooldata,method = "Anova")

# Save pairwise matrix:
p.fst <- pair.fst$PairwiseFSTmatrix
# write.csv(p.fst, file = "PairwiseFST.csv")

p.fst <- as.matrix(p.fst)
heatmap(p.fst, Rowv=NA,Colv=NA,symm=TRUE,scale="none")

p.fst <- as.data.frame(p.fst)
p.fst <- data.frame(t(apply(p.fst,1,sort)))
p.fst <- unique(p.fst)

ggpairs(p.fst, title="correlogram with ggpairs()") 
ggcorr(p.fst) 

corrgram(p.fst,order=FALSE,
         main = "Pairwise FST Test",
         lower.panel=panel.shade,upper.panel=NULL,
         diag.panel=panel.minmax)

corrgram(p.fst,order=FALSE,
         main = "Pairwise FST Test",
         lower.panel=panel.shade,upper.panel=panel.cor,
         diag.panel=panel.minmax,
         col.regions=colorRampPalette(c("lightgoldenrod1", "gold1","coral3","firebrick4")))

# PCA ####

fstmatrix <- as_tibble(fstmatrix)
data.pca <- princomp(fstmatrix)
?princomp

# Mantel Test ####
library(cluster)
# Make sure data is ready for Mantel test:
fstmatrix <- as.matrix(fstmatrix)

# Get phenotype data ready for dissimilarity matrix
phen.data <- as.data.frame(phen.data)
phen.data.daisy <- subset(phen.data, select = -Treatment)

# Make dissimilarity matrix for phenotype data
install.packages("devtools")
phen.dist <- dist(phen.data, method = "euclidean")
# phen.matrix <- daisy(phenmatrixdaisy, metric = "euclidean",stand=FALSE)

# Make distance objects
fst.dist <- as.dist(fstmatrix)

library(ade4)
mantel.rtest(fst.dist, phen.dist, nrepet = 9999)
dim(as.matrix(fstmatrix))
dim(as.matrix(phen.dist))

phen.matrix <- as.matrix(phen.dist)

# WHERE TO PICK UP:
# both vegan and ade4 package need matrices to be dist objects
# investigate dist objects - how to make them?
## need to use dist function, I think, to make dist object (instead of daisy)
# then run mantel test
