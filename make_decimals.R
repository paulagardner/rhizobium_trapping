# Easy way to import data: click "import dataset" -> select file -> comment = none

#will have to load in the dataframe
data <- as.data.frame(snp_freq_diffs_test_mincov2_rc)

##### maybe worked? 

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
