
# Before running code: ----------------------------------------------------
## Weight sheet and EA data should be copied into the IRMS batch file as additional tabs
## - Check weight sheet - ensure all rows have weight entered, even for blanks (enter 1 for blanks); 
##                      - ensure the correct description shows up in the "description" column for all rows
## - Check EA data - often, a line of zeros will be printed at bottom of EA data. This should usually be deleted.


#### LOAD PACKAGES #### (M. Gieske)

library(readxl) # Reads in .xls and .xlsx files
library(openxlsx) # Writes .xlsx files
library(dplyr) # Lets you rename columns (among other things)
library(ggplot2)
library(openxlsx)
library(tidyverse)

# Import data --------------------------------------------------------

# Read in the weigh sheet
weigh_sheet <- readxl::read_excel("WeightSheet.xlsx", 
                                  sheet = "Weight Sheet",
                                  skip = 5, 
                                  col_names = c("Analysis order", "Tray position",  
                                                "Sample ID", "description", "Target wt (mg)", "Actual wt (mg)")) %>%
  .[!is.na(.$`Actual wt (mg)`), ] # Remove any rows where there is no weight data

#Read in the raw EA data from the VarioPYRO Cube and select the desired columns
raw_data <- readxl::read_excel("CN_RawData.xlsx", sheet = "EA data")

# Look in the Environment window (upper right).  You should see the two files that were read in.
# Check to make sure they have the same number of "observations" (rows)

# Extract the data you need -----------------------------------------------
# Combine the weigh sheet and EA data while dropping columns you don't need
raw1 <- cbind(weigh_sheet[ , c(1, 3, 4, 6)], 
              raw_data[ , c(1, 2, 4, 5, 7, 8, 10)])  

# Check whether the sample weights and descriptions line up. If all the values in the selected columns match up, it will say "TRUE" in the window below.
# "FALSE" means that one or more values do not match up.  You should look and decide whether the mismatch is a real problem.  (Hint: click on the file name in the environment
# window to open it so you can scroll through it). Reminder: sometimes when loading the data into varioPYRO, the sample ID is copy/pasted
# into varioPYRO and other times the description is copy/pasted.
identical(raw1$`Actual wt (mg)`, raw1$`Weight (mg)`) 
identical(raw1$Name, raw1$`Sample ID`)

# Get the data you want for the "Processed data" sheet in your Excel workbook
#   and put the columns in the correct order
raw2 <- raw1[ , c(1:3, 5, 7:11)]

#### EXTRACT THE STANDARDS ####

acet <- raw2[grepl("ACETANILIDE", raw2$`Sample ID`, fixed=TRUE)==TRUE, ]
peach <- raw2[grepl("Peach leaves", raw2$`Sample ID`, fixed=TRUE)==TRUE, ]

# Put the three sets of standards together
stds <- rbind(acet, peach)

# Pull out the replicates
reps <- raw2[raw2$description=="unknown_replicate", ]
rep_ids <- unlist(strsplit(reps$`Sample ID`, "_"))
reps2 <- raw2[raw2$`Sample ID` %in% rep_ids, ]
reps3 <- rbind(reps, reps2)

#### GET R2, STD DEV, SLOPE, MEAN, EXPECTED VALUE, CALIBRATION VALUE #### (S. Pey, 3.24.2020) 

# Create a matrix to store the expected standard values, mean standard values and calculated
# calibration values for each of the four standards (Peach leaves, rosemount soil, acentanilide,
# and USGS 40)

qc <- matrix(NA, 12, 2)
row.names(qc) <- c("Acetanilide (r squared)", "Acetanilide (slope)", "Acetanilide (std dev)", 
                   "Acetanilide (mean)", "Aetenalide (expected)","Acetanilide (calibration)",
                   "Peach leaves (r squared)", "Peach leaves (slope)", "Peach leaves (std dev)", 
                   "Peach leaves (mean)", "Peach leaves (expected)","Peach leaves (calibration)")
colnames(qc) <- c("N [%]","C [%]") 

## [ACETANILIDE] Calculate the r-squared, standard deviation, slope and mean standard values 
## for %N anc %C. Add the expected %N and %C values to the matrix.
## Calculate the calibration value and add it to the matrix.

qc[1, 1] <- summary(lm(`N [%]` ~ `Analysis order`, data = acet))$r.squared
qc[1, 2] <- summary(lm(`C [%]` ~ `Analysis order`, data = acet))$r.squared

qc[2, 1] <- summary(lm(`N [%]` ~ `Analysis order`, data = acet))$coefficients[2, 1]
qc[2, 2] <- summary(lm(`C [%]` ~ `Analysis order`, data = acet))$coefficients[2, 1]

qc[3, 1] <- sd(acet$`N [%]`)
qc[3, 2] <- sd(acet$`C [%]`)

qc[4, 1] <- mean(acet$`N [%]`)
qc[4, 2] <- mean(acet$`C [%]`)

qc[5,1] <- 10.36 # Expected %N, acetanilide
qc[5,2] <- 71.09 # Expected %C acetanilide

qc[6, 1] <- qc[4,1] / qc[5,1] # %N calibration, acetanilide
qc[6, 2] <- qc[4,2] / qc[5,2] # %C calibration, acetanilide

## [PEACH LEAVES] Calculate the r-squared, standard deviation, slope and mean standard values 
## for %N and %C. Add the expected %N and %C values to the matrix.
## Calculate the calibration value and add it to the matrix.
qc[7, 1] <- summary(lm(`N [%]` ~ `Analysis order`, data = peach))$r.squared
qc[7, 2] <- summary(lm(`C [%]` ~ `Analysis order`, data = peach))$r.squared

qc[8, 1] <- summary(lm(`N [%]` ~ `Analysis order`, data = peach))$coefficients[2, 1]
qc[8, 2] <- summary(lm(`C [%]` ~ `Analysis order`, data = peach))$coefficients[2, 1]

qc[9, 1] <- sd(peach$`N [%]`)
qc[9, 2] <- sd(peach$`C [%]`)

qc[10, 1] <- mean(peach$`N [%]`)
qc[10, 2] <- mean(peach$`C [%]`)

qc[11,1] <- 2.28   # Expected %N, peach leaves
qc[11,2] <- 50.40  # Expected %C, peach leaves

qc[12, 1] <- qc[10, 1] / qc[11,1]  # %N calibration, peach leaves
qc[12, 2] <- qc[10, 2] / qc[11,2]  # %C calibration, peach leaves

# Round all values to two decimal places
#qc <- round(qc, 2)
# Change matrix to data frame
qc_df <- data.frame(rownames(qc), qc, row.names = NULL)  # Add rownames form the matrix as first column (MF 6/19/2020)
# data.frame() changes colname format. Rename columns (MF 6/2020)
names(qc_df) <- c("QC Value", "N [%]", "C [%]")  

# %N calibration
N.cal <- c(qc[6,1]) # [ACETANILIDE]
N.cal <- c(qc[12,1]) # [PEACH LEAVES]


# %C calibration
C.cal <- c(qc[6,2]) # [ACETANILIDE]
C.cal <- c(qc[12,2]) # [PEACH LEAVES]

## Add the calibration values to the raw2 dataframe. Double check that the correct
## calibration values appear in the table. If the calibration value is incorrect,
## make sure that you have selected the correct standard to use for calibration in 
## the above section of code. 

raw2$N.cal <- N.cal
raw2$C.cal <- C.cal

## Apply calibration standards to the data and add the corrected data values into 
## the raw2 dataframe.
raw2$corrected.N <- (raw2$`N [%]` / raw2$N.cal)
raw2$corrected.C <- (raw2$`C [%]` / raw2$C.cal)

## Calculate and add C/N ratio values into the raw2 dataframe. ### (C. Loopstra, 6/26/2020)
corrected.CNratio <- (raw2$corrected.C/raw2$corrected.N)
raw2$corrected.CNratio <- corrected.CNratio

## Change raw2 column names to their final column names
colnames(raw2) = c("Analysis order", "Sample ID", 
                   "Description", "Weight (mg)",
                   "N Area", "C Area", "N [%]", "C [%]", "C/N ratio",
                   "N calibration value", "C calibration value", 
                   "Final Corrected %N", "Final Corrected %C","Final Corrected CN ratio")
#### SAVE THE DATA ####

sheet_list <- list("Raw data" = raw_data,
                   "Processed data" = raw2,
                   "Standards" = stds,
                   "QC" = qc_df,
                   "Replicates" = reps3)

write.xlsx(sheet_list, file = "ProcessedCNdata.xlsx")  

