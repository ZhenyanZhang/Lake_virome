# Set working directory
wd <- "D:/onedrive/zzy/isme_revised/1.overview/pam"
setwd(wd)

# Load required packages
library(ggplot2)
library(reshape2)

cat("1. Loading data...\n")
data_df <- read.csv("shannon_data.csv", row.names = 1, check.names = FALSE)
cluster_info <- read.csv("PAM_JSD_cluster_assignments.csv", row.names = 1)

cat("2. Formatting data for plotting...\n")
# Build plotting dataframe
plot_df <- data.frame(
  Sample = rownames(data_df),
  Value = data_df[, 1], 
  Virometype = paste0("Virometype ", cluster_info[rownames(data_df), "Best_Cluster_Label"])
)

# Extract metric name for labels
metric_name <- colnames(data_df)[1]

cat("3. Generating boxplot...\n")
virometype_colors <- c("#E64B35", "#4DBBD5", "#00A087")

p_box <- ggplot(plot_df, aes(x = Virometype, y = Value, fill = Virometype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 19, stroke = 0.3) +
  scale_fill_manual(values = virometype_colors) +
  labs(
    x = "Virometype",
    y = paste(metric_name, "Value"),
    title = paste(metric_name, "across Virometypes")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
  theme(
    text = element_text(family = "sans", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none" 
  )

print(p_box)

cat("4. Saving plot...\n")
output_name <- paste0("Virometype_", metric_name, "_Boxplot.pdf")
# ggsave(output_name, p_box, width = 6, height = 5)

cat(paste("Done! Output saved as:", output_name, "\n"))