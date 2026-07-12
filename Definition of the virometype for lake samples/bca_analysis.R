# Set working directory
wd <- "D:/onedrive/zzy/isme_revised/redo/overview/PAM"
setwd(wd)

# Load required packages
library(ade4)
library(ggplot2)

# Load distance matrix
dist_df <- read.csv("Distance_Matrix_JSD.csv", row.names = 1, check.names = FALSE)
dist_matrix <- as.dist(dist_df)

# Load cluster assignments
cluster_info <- read.csv("PAM_JSD_cluster_assignments.csv", row.names = 1)
virometype_labels <- as.factor(cluster_info$Best_Cluster_Label)
levels(virometype_labels) <- paste0("Virometype ", levels(virometype_labels))
best_k <- length(levels(virometype_labels))

# Perform BCA dimension reduction
pco_res <- dudi.pco(dist_matrix, scannf = FALSE, nf = best_k)
bca_res <- bca(pco_res, fac = virometype_labels, scannf = FALSE, nf = 2)

plot_data <- data.frame(bca_res$ls)
colnames(plot_data) <- c("BCA1", "BCA2")
plot_data$Virometype <- virometype_labels
eig_pct <- bca_res$eig / sum(bca_res$eig) * 100

# Generate BCA spider plot
centroids <- aggregate(cbind(BCA1, BCA2) ~ Virometype, data = plot_data, FUN = mean)
colnames(centroids)[2:3] <- c("Centroid_X", "Centroid_Y")
plot_data <- merge(plot_data, centroids, by = "Virometype")

p_bca <- ggplot(plot_data, aes(x = BCA1, y = BCA2, color = Virometype, fill = Virometype)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, color = NA, show.legend = FALSE) +
  geom_segment(aes(x = BCA1, y = BCA2, xend = Centroid_X, yend = Centroid_Y),
               alpha = 1, linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 1.5, alpha = 0.8, shape = 19, stroke = 0.2) +
  labs(
    x = sprintf("PCoA1"),
    y = sprintf("PCoA2"),
    title = "Lake Virometype Landscape (Spider Plot)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),                       
    text = element_text(family = "sans", size = 12),    
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )
print(p_bca)

# Monte-Carlo permutation test
rt <- randtest(bca_res, nrepet = 999)
print(rt)