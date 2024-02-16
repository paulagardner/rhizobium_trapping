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

# Apply the function only to the specified columns
smalldata[, start_column_index:ncol(smalldata)] <- lapply(smalldata[, start_column_index:ncol(smalldata)], divideFraction)

# Display the updated data frame
print(smalldata)

}

# Example usage: Assuming your data frame is named your_data and you want to start from column 3
smalldata <- divide_and_update(smalldata, start_column = 3)
