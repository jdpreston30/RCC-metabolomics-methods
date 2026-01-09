#* 5: Assign Plots to Figure Panels
#+ 4.1: Figure 1 Plots
p1A <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig1A.png")), interpolate = TRUE)
p1B <- num_features_ratio
p1C <- mean_intensity_ratio
p1D <- ratio_table_grob
#+ 4.2: Figure 2 Plots
p2A <- plot_features_time
p2B <- plot_intensity_time
#+ 4.6: Supplementary Figure 1 Plots
s1A <- plot_intensity_time_normal
s1B <- plot_features_time_normal