#* 13: Render Figures
#+ 13.0: Clean up any corrupted PDFs from previous runs
if (file.exists("Outputs/Figures/Figs1-3.pdf")) {
  file.remove("Outputs/Figures/Figs1-3.pdf")
}
#+ 13.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 1A
  draw_plot(p1A, x = 1.18, y = 6.774, 
                 width = 2.98, height = 3.055) +
  # 1B
  draw_plot(p1B, x = 4.23, y = 7.285,
                 width = 3.16, height = 2.6) +
  # Labels
  figure_labels(list(
    A = c(1.26, 10.00),
    B = c(4.35, 10.00),
    "Figure 1" = c(0.49, 10.43)
  ))
#+ 13.2: Figure 2
fig2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 2A
  draw_plot(p2A,
    x = 0.53, y = 3.699,
    width = 4.6, height = 6.28
  ) +
  # 2B
  draw_plot(p2B,
    x = 4.55, y = 7.016,
    width = 2.7945, height = 2.03
  ) +
  # 2C
  draw_plot(p2C,
    x = 4.55, y = 5.016,
    width = 2.8, height = 2.03
  ) +
  # Labels
  figure_labels(list(
    A = c(0.75, 10.00),
    B = c(4.6, 8.97),
    C = c(4.6, 6.97),
    "Figure 2" = c(0.49, 10.43)
  ))
#+ 13.3: Figure 3
fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(p3A,
    x = 0.75, y = 6.838,
    width = 7, height = 3.5
  ) +
  draw_plot(p3B,
    x = 0.7258, y = 2.889133,
    width = 3.346663, height = 3.74233
  ) +
  draw_plot(p3C,
    x = 3.96875, y = 2.935833,
    width = 4, height = 3.75
  ) +
  draw_plot(p3D,
    x = 0.7182667, y = 0.4446667,
    width = 2.75, height = 2.386667
  ) +
  draw_plot(p3E,
    x = 3.3, y = 0.388,
    width = 2.5, height = 2.5
  ) +
  draw_plot(p3F,
    x = 5.335+0.25, y = 0.388,
    width = 2.5, height = 2.5
  ) +
  # Labels
  figure_labels(list(
    A = c(0.785, 10.125),
    B = c(0.785, 6.563333),
    C = c(4.14, 6.563333),
    D = c(0.785, 2.75),
    E = c(3.424, 2.75),
    F = c(5.709, 2.75),
    "Figure 3" = c(0.49, 10.43)
  ))
#+ 13.4: Print All Main Figures
#- 13.4.1: As PNGs
print_to_png(fig1, "Fig1.png", output_dir = "Outputs/Figures/PNG")
print_to_png(fig2, "Fig2.png", output_dir = "Outputs/Figures/PNG")
print_to_png(fig3, "Fig3.png", output_dir = "Outputs/Figures/PNG")
#- 13.4.2: As TIFFs
print_to_tiff(fig1, "Fig1.tiff", output_dir = "Outputs/Figures/TIFF")
print_to_tiff(fig2, "Fig2.tiff", output_dir = "Outputs/Figures/TIFF")
print_to_tiff(fig3, "Fig3.tiff", output_dir = "Outputs/Figures/TIFF")
#+ 13.5: Save All Main Figures as PDF from PNGs
# Close any open graphics devices
while (!is.null(dev.list())) { 
  dev.off() 
}
pdf("Outputs/Figures/Figs1-3.pdf", width = 8.5, height = 11)
# Page 1: Fig1
img1 <- readPNG("Outputs/Figures/Fig1.png")
grid::grid.newpage()
grid::grid.raster(img1, width = grid::unit(8.5, "inches"), height = grid::unit(11, "inches"))
# Page 2: Fig2
img2 <- readPNG("Outputs/Figures/Fig2.png")
grid::grid.newpage()
grid::grid.raster(img2, width = grid::unit(8.5, "inches"), height = grid::unit(11, "inches"))
# Page 3: Fig3
img3 <- readPNG("Outputs/Figures/Fig3.png")
grid::grid.newpage()
grid::grid.raster(img3, width = grid::unit(8.5, "inches"), height = grid::unit(11, "inches"))
dev.off()
cat("PDF compiled: Outputs/Figures/Figs1-3.pdf\n")
