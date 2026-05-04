# Generate bambamR logo using base R graphics only
# Run once: Rscript inst/extdata/create_logo.R

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

png("man/figures/logo.png", width = 518, height = 600, bg = "transparent")

par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

# Flat-top hexagon centered at (0.5, 0.5)
angles <- seq(0, 2 * pi, length.out = 7)
hx <- 0.5 + 0.49 * cos(angles)
hy <- 0.5 + 0.49 * sin(angles)

# Fill + border
polygon(hx, hy, col = "#2166AC", border = "#164578", lwd = 5)

# DNA helix motif (centered vertically at 0.48)
n <- 100
t <- seq(0, 3 * pi, length.out = n)
cx <- seq(0.15, 0.85, length.out = n)
cy1 <- 0.48 + sin(t) * 0.10
cy2 <- 0.48 - sin(t) * 0.10

lines(cx, cy1, col = "#DCEEFB", lwd = 3.5)
lines(cx, cy2, col = "#F8D7DA", lwd = 3.5)

# Base-pair rungs
bp <- seq(5, n, by = 9)
segments(cx[bp], cy2[bp], cx[bp], cy1[bp],
         col = "#FFFFFF55", lwd = 1.5)

# Bar chart motif at the bottom
bar_x <- seq(0.30, 0.70, length.out = 6)
bar_h <- c(0.04, 0.07, 0.11, 0.06, 0.13, 0.08)
bar_base <- 0.18
for (i in seq_along(bar_x)) {
  rect(bar_x[i] - 0.028, bar_base,
       bar_x[i] + 0.028, bar_base + bar_h[i],
       col = "#A8CCE8", border = NA)
}

# Package name at the top
text(0.5, 0.76, "bambamR", cex = 5, col = "white",
     font = 2, family = "sans")

# Tagline at bottom
text(0.5, 0.12, "RNA-seq toolkit", cex = 1.5, col = "#8BB8DE",
     font = 3, family = "sans")

dev.off()

cat("Logo saved to man/figures/logo.png\n")
cat("Size:", file.size("man/figures/logo.png"), "bytes\n")
