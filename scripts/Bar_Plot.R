
df <- read.csv("top_5_genera_per_sample.csv", stringsAsFactors = FALSE)

long_df <- data.frame(
  Sample     = rep(df$Sample, 3),
  Diet       = rep(df$Diet, 3),
  Genus      = c(df$Genus_label_1, df$Genus_label_2, df$Genus_label_3),
  Abundance  = c(df$Abundance_1, df$Abundance_2, df$Abundance_3)
)

diets  <- unique(long_df$Diet)
genera <- unique(long_df$Genus)

Diet_vec   <- character()
Genus_vec  <- character()
Mean_vec   <- numeric()
SE_vec     <- numeric()

for(d in diets){
  for(g in genera){
    vals <- long_df$Abundance[long_df$Diet == d & long_df$Genus == g]
    m  <- mean(vals, na.rm = TRUE)
    s  <- sd(vals, na.rm = TRUE)
    n  <- length(vals)
    se <- if(n > 1) s / sqrt(n) else 0
    
    Diet_vec  <- c(Diet_vec, d)
    Genus_vec <- c(Genus_vec, g)
    Mean_vec  <- c(Mean_vec, m)
    SE_vec    <- c(SE_vec, se)
  }
}

summary_df <- data.frame(
  Diet   = Diet_vec,
  Genus  = Genus_vec,
  Mean   = Mean_vec,
  SE     = SE_vec,
  stringsAsFactors = FALSE
)

means_matrix <- tapply(summary_df$Mean, list(summary_df$Genus, summary_df$Diet), function(x) x)
se_matrix    <- tapply(summary_df$SE,   list(summary_df$Genus, summary_df$Diet), function(x) x)

genus_colors <- c("Desulfovibrio" = "steelblue", 
                  "Klebsiella"    = "darkorange", 
                  "Proteus"       = "forestgreen")

bar_positions <- barplot(
  height = means_matrix,
  beside = TRUE,
  col    = genus_colors[rownames(means_matrix)],
  ylim   = c(0, max(summary_df$Mean + summary_df$SE, na.rm = TRUE) * 1.25),
  ylab   = "Mean Relative Abundance",
  xlab   = "Diet Group",
  main   = "Top Genera Abundance by Diet",
  las    = 1,                    
  cex.main = 1.4,
  cex.lab  = 1.1,
  legend.text = rownames(means_matrix),
  args.legend = list(x = "topright", bty = "n", cex = 0.95)
)

grid(nx = NA, ny = NULL, col = "gray75", lty = "dotted", lwd = 1)

arrows(
  x0 = bar_positions,
  y0 = means_matrix - se_matrix,
  x1 = bar_positions,
  y1 = means_matrix + se_matrix,
  angle = 90,
  code  = 3,          
  length = 0.06,
  lwd   = 1.6,
  col   = "black"
)

text(
  x      = bar_positions,
  y      = means_matrix + se_matrix + 30,
  labels = round(means_matrix, 0),
  cex    = 0.75,
  col    = "black"
)

