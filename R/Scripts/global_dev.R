main_uft_20p_tumors <- main_uft_20p %>%
  filter(tissue_type == "Tumor") |>
  select(-tissue_type) |>
  filter(tx != "CCM")

#* 2.1: Global Metabolic Deviation Analysis
#+ 2.1.1: Calculate deviations from baseline using utility function
deviation_results <- calculate_global_deviations(
  data = main_uft_20p_tumors,
  baseline_tx = "GS",
  feature_pattern = "^(HILIC|C18)_"
)

# Extract results
global_deviations_long <- deviation_results$deviations_long %>%
  mutate(
    tx = factor(tx, levels = c("Ice", "RT")),
    time = factor(time, levels = c("15", "30", "60", "90"))
  ) %>%
  rename(delta_from_GS = delta_from_baseline)

global_mean_deviation <- deviation_results$mean_deviation %>%
  mutate(
    tx = factor(tx, levels = c("Ice", "RT")),
    time = factor(time, levels = c("15", "30", "60", "90"))
  )

feature_cols <- deviation_results$feature_cols

#+ 2.1.2: Statistical analysis using utility function
analysis_results <- analyze_global_deviations(
  mean_deviation_data = global_mean_deviation,
  response_var = "mean_delta",
  print_results = TRUE
)

# Extract results
global_model <- analysis_results$model
global_anova <- analysis_results$anova
global_summary <- analysis_results$summary_stats

#+ 2.1.3: Visualization - Global metabolic trajectory
{
  # Add baseline point at time = 0 (GS by definition is 0 deviation)
  baseline_global <- tibble(
    tx = factor(c("Ice", "RT"), levels = c("Ice", "RT")),
    time = factor("0", levels = c("0", "15", "30", "60", "90")),
    mean = 0,
    se = 0,
    sd = 0
  )
  
  global_summary_with_baseline <- bind_rows(baseline_global, global_summary) %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  # Create plot
  tx_colors <- c("Ice" = "#0958b3ff", "RT" = "#5d0f10")
  
  plot_global_trajectory <- ggplot(global_summary_with_baseline, 
                                   aes(x = time_numeric, y = mean, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_color_manual(values = tx_colors, labels = c("Ice", "RT")) +
    scale_x_continuous(
      breaks = c(0, 15, 30, 60, 90),
      labels = c("0", "15", "30", "60", "90")
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(size = 8, face = "bold", color = "black"),
      axis.title = element_text(size = 10, face = "bold", color = "black"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold")
    ) +
    labs(
      title = "Global Metabolic Deviation from Baseline",
      x = "Time (minutes)",
      y = "Mean Log₂ Fold Change from Baseline"
    )
  
  print(plot_global_trajectory)
  
  # Save figure
  print_to_png(plot_global_trajectory, "global_metabolic_trajectory.png", 
               width = 6, height = 4.5, output_dir = "Outputs/Figures", auto_open = TRUE)
}
