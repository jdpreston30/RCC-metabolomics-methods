#* 6: Render Figures
#+ 6.0: Clean up any corrupted PDFs from previous runs
if (file.exists("Outputs/Figures/Figs1-3.pdf")) {
  file.remove("Outputs/Figures/Figs1-3.pdf")
}
#+ 6.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 1A
  draw_plot(p1A, x = 0.645, y = 8.986667, 
                 width = 7.225, height = 1.02) +
  # 1B
  draw_plot(p1B, x = 0.7, y = 5.62,
                 width = 3.5, height = 3) +
  # 1C
  draw_plot(p1C, x = 4.283333333, y = 5.62,
                 width = 3.5, height = 3) +
  # 1D
  draw_plot(p1D, x = 1.75, y = 2.993333333,
                 width = 1.786666667, height = 3) +
  # Labels
  figure_labels(list(
    A = c(0.635, 10.00),
    B = c(0.635, 8.5),
    C = c(4.218333333, 8.5),
    D = c(0.635, 5.203333333),
    "Figure 1" = c(0.49, 10.43)
  ))
#+ 6.5: Print Final Figures
print_to_png(fig1, "Fig1.png", output_dir = "Outputs/Figures/PNG")
