# ==============================================================================
# Preprocess Lake Virome Abundance Matrix and Rarefaction Analysis
# ==============================================================================

# Set working directory
wd <- "D:/onedrive/zzy/isme_revised/1.overview"
setwd(wd)

# Load required packages
library(dplyr)
library(tidyr)
library(rtk)
library(ggplot2)

cat("Reading and cleaning raw data...\n")

# Read merged abundance file
df <- read.delim("all_abundance.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Remove ".temp" suffix from sample names
df$Sample <- gsub("\\.temp$", "", df$Sample)

cat("Converting long format to abundance matrix...\n")

# Pivot table to wide format
abundance_matrix <- df %>%
  select(Sample, Contig, `Trimmed Mean`) %>%
  pivot_wider(names_from = Sample, 
              values_from = `Trimmed Mean`, 
              values_fill = 0)

# Convert to standard dataframe with Contig as rownames
abundance_df <- as.data.frame(abundance_matrix)
rownames(abundance_df) <- abundance_df$Contig
abundance_df$Contig <- NULL

# Export cleaned abundance matrix
write.csv(abundance_df, "abundance_matrix.csv", quote = FALSE)

cat("Generating Presence/Absence matrix...\n")

# Create 0/1 matrix
pa_df <- abundance_df
pa_df[pa_df > 0] <- 1

cat("Calculating accumulation curve (100 permutations)...\n")

# Run rtk calculation
curve_matrix <- collectors.curve(pa_df, times = 100, bin = 1)

# Format conversion to handle different rtk outputs safely
if (is.list(curve_matrix) && !is.data.frame(curve_matrix)) {
  curve_matrix <- do.call(rbind, curve_matrix)
}

if (is.null(dim(curve_matrix))) {
  curve_matrix <- as.matrix(curve_matrix)
}

curve_matrix <- as.data.frame(curve_matrix)

# Transpose matrix if necessary
if (nrow(curve_matrix) != 100 && ncol(curve_matrix) == 100) {
  curve_matrix <- as.data.frame(t(curve_matrix))
}

cat("Calculating Mean and 95% Confidence Intervals...\n")

# Calculate statistics
plot_data <- data.frame(
  Sample_Size = 1:ncol(curve_matrix),
  Mean_vOTUs  = apply(curve_matrix, 2, mean, na.rm = TRUE),
  SD          = apply(curve_matrix, 2, sd, na.rm = TRUE)
)

# Calculate 95% CI
plot_data$CI_Lower <- pmax(0, plot_data$Mean_vOTUs - 1.96 * plot_data$SD)
plot_data$CI_Upper <- plot_data$Mean_vOTUs + 1.96 * plot_data$SD

# Plotting using ggplot2
p <- ggplot(plot_data, aes(x = Sample_Size, y = Mean_vOTUs)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), 
              fill = "#4A90E2", alpha = 0.3, color = NA) +
  geom_line(color = "#003366", linewidth = 1.2) +
  theme_classic() +
  labs(x = "Number of samples", 
       y = "Number of detected viral OTUs") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.line = element_line(linewidth = 0.8)
  )

# Display plot
print(p)

# Export as high-resolution PDF
ggsave("Rarefaction_Curve.pdf", plot = p, width = 8, height = 6, dpi = 300)

cat("Analysis and plotting finished!\n")