# Easy way to import data: click "import dataset" -> select file -> comment = none

#will have to load in the dataframe
data <- as.data.frame(snp_freq_diffs_test_mincov2_rc)
smalldata <- data[9:10,]
#for (column in small.data[,10:34])
#  for (entry in column)
#    print(entry)
#    fraction_vector <- unlist(strsplit(entry, "/"))
#    numerator <- as.numeric(fraction_vector[1])
#    denominator <- as.numeric(fraction_vector[2])
#    if (denominator != 0) {
#      result <- numerator / denominator
#    } 
#    print(result)

# Function to perform division on a fraction string
divideFraction <- function(fraction_string) {
  fraction_vector <- unlist(strsplit(as.character(fraction_string), "/"))
  numerator <- as.numeric(fraction_vector[1])
  denominator <- as.numeric(fraction_vector[2])
  if (denominator != 0) {
    return(numerator / denominator)
  } else {
    return(NA)  # Handle division by zero as needed
  }
}


# Specify the starting column index to update onward
start_column_index <- 10  # Starting from "fraction" column onward

# Try just the first line
divideFraction_test1 <- function(fraction_string) {
  fraction_vector <- unlist(strsplit(as.character(fraction_string), "/"))
}

test1 <- lapply(smalldata[, start_column_index:ncol(smalldata)], divideFraction_test1)

# Here is the weird thing: it is not differentiating among rows.
# For example, V10 is 22 22 22 22
# The right values are being extracted, but they are all clumped together

# Try it without if_else just to see
divideFraction_test2 <- function(fraction_string) {
  fraction_vector <- unlist(strsplit(as.character(fraction_string), "/"))
  numerator <- as.numeric(fraction_vector[1])
  denominator <- as.numeric(fraction_vector[2]) 
  }

test2 <- lapply(smalldata[, start_column_index:ncol(smalldata)], divideFraction_test2)

# This seemed to work


# Apply the function only to the specified columns
test3 <- lapply(smalldata[, start_column_index:ncol(smalldata)], divideFraction)
# in some sense, this worked. 
# For the first line, you expect 1, 1, then a number smaller than 1
# We get that: 1, 1, 0.98
# But where is the other row?
# Is the issue that we aren't applying to all rows?


smalldata2 <- lapply(smalldata[, start_column_index:ncol(smalldata)], divideFraction)

smalldata[, start_column_index:ncol(smalldata)] <- lapply(smalldata[, start_column_index:ncol(smalldata)], divideFraction)



# Display the updated data frame
print(smalldata)



# Example usage: Assuming your data frame is named your_data and you want to start from column 3
smalldata <- divide_and_update(smalldata, start_column = 3)

# Fuck it we ball: try it on the whole data frame
test4 <- lapply(data[, start_column_index:ncol(smalldata)], divideFraction)
# All NAs, just one column. Need to learn to apply what we did on test3 to a whole data frame.


# Try in tidyverse ####
data.narrow <- data %>%
  select(V1,V10:V50)

data.narrow %>%
  select(V10:V50) %>%
  divideFraction_test()
#Error in if (denominator != 0) { : missing value where TRUE/FALSE needed
#  In addition: Warning messages:
#    1: In divideFraction(.) : NAs introduced by coercion
#  2: In divideFraction(.) : NAs introduced by coercion


data.output <- data.narrow %>%
  select(V10:V50) %>%
  divideFraction_test2()
# Warning messages:
#  1: In divideFraction_test2(.) : NAs introduced by coercion
# 2: In divideFraction_test2(.) : NAs introduced by coercion

# Still didn't work! Need Chat GPT to help!
