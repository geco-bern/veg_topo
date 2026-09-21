library(terra)
source(here::here("R/config.R"))

r_H_TWI <- rast(cor_twi_vegh_mosaic_file)
mi <- rast(mi_5km_file) * 0.0001
fused <- rast(fused_5km_file)

stacked <- c(r_H_TW, mi, fused)

# ============= Figure 1: All data =================================================

# Equal frequency binning
n_bins <- 100
mi_values_all <- values(mi)[,1]

# Calculate quantiles as break points
probs <- seq(0, 1, length.out = n_bins + 1)
breaks_all <- quantile(mi_values_all, probs = probs, na.rm = TRUE)

# Remove duplicates
if (length(unique(breaks_all)) < length(breaks_all)) {
  cat("Warning: Duplicate quantile break points detected\n")
  breaks_all <- unique(breaks_all)
  cat("Actual number of bins:", length(breaks_all) - 1, "\n")
}

# Binning
mi_bins_all <- cut(mi_values_all, breaks = breaks_all, include.lowest = TRUE)

# Extract r_H_TWI values
r_H_TW_values_all <- values(r_H_TW)[,1]

# Create data frame
df_all <- data.frame(
  r_H_TWI = r_H_TW_values_all,
  mi_bin = mi_bins_all
)
df_all <- na.omit(df_all)

# Plot Figure 1: All data
dev.new()
par(mar = c(8, 5, 4, 5))

boxplot(r_H_TWI ~ mi_bin, data = df_all,
        xlab = "",
        ylab = "",
        las = 2,
        cex.axis = 0.6,
        outline = FALSE,
        main = "All Pixels")

# Add x-axis title
mtext("MI Bins (equal frequency)", side = 1, line = 6, cex = 0.9)

# Add y-axis title
mtext("r_H_TWI Values", side = 2, line = 4, cex = 0.9)

# Add y = 0 reference line
abline(h = 0, col = "red", lwd = 2, lty = 2)

# Add sample size information
text(x = nlevels(df_all$mi_bin) * 0.9,
     y = max(df_all$r_H_TW, na.rm = TRUE) * 0.95,
     labels = paste("n =", format(nrow(df_all), scientific = FALSE)),
     cex = 0.8, adj = 1)

# ============= Figure 2: fused < 0.05 data ======================================

# Extract fused values
fused_values <- values(fused)[,1]

# Create filter condition
filter_fused <- fused_values < 0.05 & !is.na(fused_values)

# Filtered mi and r_H_TWI values
mi_values_filtered <- mi_values_all[filter_fused]
r_H_TW_values_filtered <- r_H_TW_values_all[filter_fused]

# Perform equal frequency binning on filtered data (based on filtered mi values)
breaks_filtered <- quantile(mi_values_filtered, probs = probs, na.rm = TRUE)

# Remove duplicates
if (length(unique(breaks_filtered)) < length(breaks_filtered)) {
  breaks_filtered <- unique(breaks_filtered)
}

# Binning
mi_bins_filtered <- cut(mi_values_filtered, breaks = breaks_filtered, include.lowest = TRUE)

# Create filtered data frame
df_filtered <- data.frame(
  r_H_TWI = r_H_TW_values_filtered,
  mi_bin = mi_bins_filtered
)
df_filtered <- na.omit(df_filtered)

# Output statistics
cat("\n============ Data Statistics ================================================\n")
cat("Original data pixel count:", nrow(df_all), "\n")
cat("fused < 0.05 pixel count:", nrow(df_filtered), "\n")
cat("Filtered percentage:", round(nrow(df_filtered)/nrow(df_all)*100, 2), "%\n")

# Plot Figure 2: fused < 0.05 data
dev.new()
par(mar = c(8, 5, 4, 5))

boxplot(r_H_TWI ~ mi_bin, data = df_filtered,
        xlab = "",
        ylab = "",
        las = 2,
        cex.axis = 0.6,
        outline = FALSE,
        main = expression("Pixels with fused < 0.05"))

# Add x-axis title
mtext("MI Bins (equal frequency)", side = 1, line = 6, cex = 0.9)

# Add y-axis title
mtext("r_H_TWI Values", side = 2, line = 4, cex = 0.9)

# Add y = 0 reference line
abline(h = 0, col = "red", lwd = 2, lty = 2)

# Add sample size information
text(x = nlevels(df_filtered$mi_bin) * 0.9,
     y = max(df_filtered$r_H_TW, na.rm = TRUE) * 0.95,
     labels = paste("n =", format(nrow(df_filtered), scientific = FALSE)),
     cex = 0.8, adj = 1)

# ============= Save figures ====================================================

# Save Figure 1
png(here::here("data/figures/4_04_boxplot_by_mi_bin_all_land.png"), width = 3000, height = 2000, res = 300)
par(mar = c(8, 5, 4, 5))
boxplot(r_H_TWI ~ mi_bin, data = df_all,
        xlab = "",
        ylab = "",
        las = 2,
        cex.axis = 0.6,
        outline = FALSE,
        main = "All Pixels")
mtext("MI Bins (equal frequency)", side = 1, line = 6, cex = 0.9)
mtext("r_H_TWI Values", side = 2, line = 4, cex = 0.9)
abline(h = 0, col = "red", lwd = 2, lty = 2)
text(x = nlevels(df_all$mi_bin) * 0.9,
     y = max(df_all$r_H_TW, na.rm = TRUE) * 0.95,
     labels = paste("n =", format(nrow(df_all), scientific = FALSE)),
     cex = 0.8, adj = 1)
dev.off()

# Save Figure 2
png(here::here("data/figures/4_04_boxplot_by_mi_bin_natural_land.png"), width = 3000, height = 2000, res = 300)
par(mar = c(8, 5, 4, 5))
boxplot(r_H_TWI ~ mi_bin, data = df_filtered,
        xlab = "",
        ylab = "",
        las = 2,
        cex.axis = 0.6,
        outline = FALSE,
        main = expression("Pixels with fused < 0.05"))
mtext("MI Bins (equal frequency)", side = 1, line = 6, cex = 0.9)
mtext("r_H_TWI Values", side = 2, line = 4, cex = 0.9)
abline(h = 0, col = "red", lwd = 2, lty = 2)
text(x = nlevels(df_filtered$mi_bin) * 0.9,
     y = max(df_filtered$r_H_TW, na.rm = TRUE) * 0.95,
     labels = paste("n =", format(nrow(df_filtered), scientific = FALSE)),
     cex = 0.8, adj = 1)
dev.off()
