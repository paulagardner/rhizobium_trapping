# Load necessary libraries
library(poolfstat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(wesanderson)  # for color palettes
library(cowplot)      # for arranging plots if needed

# Step 1: Load the pool data
pooldata <- vcf2pooldata(
  vcf.file = "C:/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/2025revisionsrepo/scripts/missingness_filtered_output.vcf",
  poolsizes = rep(50, 24),
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
  nlines.per.readblock = 1e+06
)

# Step 2: Compute pairwise FST
PairwiseFST <- compute.pairwiseFST(pooldata)
fst_matrix <- PairwiseFST@PairwiseFSTmatrix

# Step 3: Clean sample names and extract metadata
clean_names <- gsub("\\.sorted$", "", rownames(fst_matrix))
rownames(fst_matrix) <- clean_names
colnames(fst_matrix) <- clean_names

split_parts <- strsplit(clean_names, "_")
ecotypes <- sapply(split_parts, `[`, 1)
locations <- sapply(split_parts, `[`, 2)
temperatures <- sapply(split_parts, `[`, 3)
names(ecotypes) <- names(locations) <- names(temperatures) <- clean_names

# Step 4: Reorder and melt function
reorder_and_melt <- function(matrix, group_vector) {
  ordered_names <- names(sort(group_vector))
  reordered_matrix <- matrix[ordered_names, ordered_names]
  melted <- melt(reordered_matrix)
  melted$Var1 <- factor(melted$Var1, levels = ordered_names)
  melted$Var2 <- factor(melted$Var2, levels = ordered_names)
  return(melted)
}

fst_by_ecotype <- reorder_and_melt(fst_matrix, ecotypes)
fst_by_location <- reorder_and_melt(fst_matrix, locations)
fst_by_temperature <- reorder_and_melt(fst_matrix, temperatures)

# Step 5: Axis background color function
axis_coloring <- function(names_vec, palette_name) {
  factor_vec <- factor(names_vec, levels = unique(names_vec))
  col_vec <- wes_palette(palette_name, length(levels(factor_vec)), type = "continuous")
  color_map <- setNames(col_vec, levels(factor_vec))
  color_assign <- color_map[as.character(factor_vec)]
  return(color_assign)
}

# Step 6: Plotting function
plot_heatmap <- function(melted, title, group_vector, palette_name) {
  # Reorder group vector based on factor levels in plot
  sample_levels <- levels(melted$Var1)
  group_vector <- group_vector[sample_levels]
  
  axis_colors <- axis_coloring(group_vector, palette_name)
  
  ggplot(melted, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "lightskyblue2", high = "salmon1", name = "FST") +
    labs(title = title, x = "Sample", y = "Sample") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = axis_colors),
      axis.text.y = element_text(color = axis_colors)
    )
}


# Step 7: Generate plots
p1 <- plot_heatmap(fst_by_ecotype, "Pairwise FST grouped by Ecotype", ecotypes, "Darjeeling1")
p2 <- plot_heatmap(fst_by_location, "Pairwise FST grouped by Location", locations, "Royal1")
p3 <- plot_heatmap(fst_by_temperature, "Pairwise FST grouped by Temperature", temperatures, "Rushmore1")

# Step 8: Display each plot
print(p1)
print(p2)
print(p3)