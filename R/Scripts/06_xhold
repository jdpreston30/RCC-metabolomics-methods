# #+ 13.4: Print All Main Figures
# #- 13.4.1: As PNGs
# print_to_png(fig1, "Fig1.png", output_dir = "Outputs/Figures/PNG")
# #+ 13.5: Save All Main Figures as PDF from PNGs
# # Close any open graphics devices
# while (!is.null(dev.list())) { 
#   dev.off() 
# }
# pdf("Outputs/Figures/Figs1-3.pdf", width = 8.5, height = 11)
# # Page 1: Fig1
# img1 <- readPNG("Outputs/Figures/Fig1.png")
# grid::grid.newpage()
# grid::grid.raster(img1, width = grid::unit(8.5, "inches"), height = grid::unit(11, "inches"))
# # Page 2: Fig2
# img2 <- readPNG("Outputs/Figures/Fig2.png")
# grid::grid.newpage()
# grid::grid.raster(img2, width = grid::unit(8.5, "inches"), height = grid::unit(11, "inches"))
# # Page 3: Fig3
# img3 <- readPNG("Outputs/Figures/Fig3.png")
# grid::grid.newpage()
# grid::grid.raster(img3, width = grid::unit(8.5, "inches"), height = grid::unit(11, "inches"))
# dev.off()
# cat("PDF compiled: Outputs/Figures/Figs1-3.pdf\n")
