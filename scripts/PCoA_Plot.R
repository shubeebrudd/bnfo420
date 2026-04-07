library(vegan)
library(ggplot2)

df <- read.csv("top_5_genera_per_sample.csv")

n_samples <- nrow(df)
abund_matrix <- matrix(0, nrow = n_samples, ncol = 3)
rownames(abund_matrix) <- df$Sample
colnames(abund_matrix) <- c("Desulfovibrio", "Klebsiella", "Proteus")

for(i in 1:n_samples) {
  # Genus 1
  if(df$Genus_label_1[i] == "Desulfovibrio") {
    abund_matrix[i, "Desulfovibrio"] <- df$Abundance_1[i]
  } else if(df$Genus_label_1[i] == "Klebsiella") {
    abund_matrix[i, "Klebsiella"] <- df$Abundance_1[i]
  } else if(df$Genus_label_1[i] == "Proteus") {
    abund_matrix[i, "Proteus"] <- df$Abundance_1[i]
  }
  
  # Genus 2
  if(df$Genus_label_2[i] == "Desulfovibrio") {
    abund_matrix[i, "Desulfovibrio"] <- df$Abundance_2[i]
  } else if(df$Genus_label_2[i] == "Klebsiella") {
    abund_matrix[i, "Klebsiella"] <- df$Abundance_2[i]
  } else if(df$Genus_label_2[i] == "Proteus") {
    abund_matrix[i, "Proteus"] <- df$Abundance_2[i]
  }
  
  # Genus 3
  if(df$Genus_label_3[i] == "Desulfovibrio") {
    abund_matrix[i, "Desulfovibrio"] <- df$Abundance_3[i]
  } else if(df$Genus_label_3[i] == "Klebsiella") {
    abund_matrix[i, "Klebsiella"] <- df$Abundance_3[i]
  } else if(df$Genus_label_3[i] == "Proteus") {
    abund_matrix[i, "Proteus"] <- df$Abundance_3[i]
  }
}

total_abund <- rowSums(abund_matrix)
valid <- total_abund > 0

abund_matrix <- abund_matrix[valid, ]
samples_valid <- df$Sample[valid]
diet_valid    <- df$Diet[valid]

cat("Samples kept:", sum(valid), "\n")
cat("Samples removed (all zeros):", sum(!valid), "\n")

#Bray-Curtis distance
dist_bray <- vegdist(abund_matrix, method = "bray")

pcoa_result <- cmdscale(dist_bray, k = 2, eig = TRUE)

pcoa_coords <- data.frame(
  Sample = samples_valid,
  PC1    = pcoa_result$points[,1],
  PC2    = pcoa_result$points[,2],
  Diet   = diet_valid
)

# Calculate variance 
eigenvals <- pcoa_result$eig
var_explained <- (eigenvals / sum(eigenvals)) * 100
pc1_label <- paste0("PCoA1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PCoA2 (", round(var_explained[2], 1), "%)")

ggplot(pcoa_coords, aes(x = PC1, y = PC2, color = Diet, shape = Diet)) +
  geom_point(size = 4.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  labs(
    title = "PCoA of Gut Microbiome Samples",
    subtitle = "Bray-Curtis distance, colored by Diet",
    x = pc1_label,
    y = pc2_label
  ) +
  scale_color_manual(values = c("Vegetarian" = "#1b9e77", "Omnivore" = "#d95f02")) +
  scale_shape_manual(values = c("Vegetarian" = 16, "Omnivore" = 17)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

