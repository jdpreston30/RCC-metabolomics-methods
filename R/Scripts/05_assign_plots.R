#* 5: Assign Plots to Figure Panels
#+ 4.1: Figure 1 Plots
p1A <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig1A.png")), interpolate = TRUE)
p1B <- num_features_ratio
p1C <- mean_intensity_ratio
p1D <- ratio_table_grob
#+ 4.2: Figure 2 Plots