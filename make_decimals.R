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

# Function to divide "16/16" and update the data frame
divide_and_update <- function(data, start_column) {
  # Get the number of rows in the data frame
  num_rows <- nrow(data)
  
  # Loop through each row and perform the division
  for (i in 1:num_rows) {
    # Extract the "16/16" string from the specified columns
    fraction_vector <- unlist(strsplit(as.character(data[i, start_column]), "/"))
    numerator <- as.numeric(fraction_vector[1])
    denominator <- as.numeric(fraction_vector[2])
    
    # Check if the denominator is not zero to avoid division by zero
    if (denominator != 0) {
      # Perform the division and update the original data frame with the result
      data[i, start_column] <- numerator / denominator
    } else {
      # Handle the case where the denominator is zero (you may want to do something specific)
      data[i, start_column] <- NA
    }
  }
  # Return the updated data frame
  return(data)
}

# Example usage: Assuming your data frame is named your_data and you want to start from column 3
smalldata <- divide_and_update(smalldata, start_column = 3)
