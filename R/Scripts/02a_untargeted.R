main_uft_full_0s_tumor <- main_uft_full_0s %>%
  filter(tissue_type == "Tumor") |>
  select(-tissue_type)

#* 2a.1: Time Course Analysis - Mean Intensity and Feature Detection
#+ 2a.1.1: Calculate metrics per sample and mode
{
  # Filter out CCM
  time_course_data <- main_uft_full_0s_tumor %>%
    filter(tx != "CCM")
  
  # Get HILIC and C18 feature columns
  hilic_cols <- colnames(time_course_data)[str_detect(colnames(time_course_data), "^HILIC_")]
  c18_cols <- colnames(time_course_data)[str_detect(colnames(time_course_data), "^C18_")]
  
  # Calculate metrics for each sample
  time_course_metrics <- time_course_data %>%
    rowwise() %>%
    mutate(
      # HILIC metrics
      hilic_mean_intensity = mean(c_across(all_of(hilic_cols))[c_across(all_of(hilic_cols)) > 0], na.rm = TRUE),
      hilic_num_features = sum(c_across(all_of(hilic_cols)) > 0, na.rm = TRUE),
      # C18 metrics
      c18_mean_intensity = mean(c_across(all_of(c18_cols))[c_across(all_of(c18_cols)) > 0], na.rm = TRUE),
      c18_num_features = sum(c_across(all_of(c18_cols)) > 0, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Sample_ID, pt_tissue_ID, tx, time, 
           hilic_mean_intensity, hilic_num_features, 
           c18_mean_intensity, c18_num_features)
  
  # Convert to long format for plotting
  time_course_long <- time_course_metrics %>%
    pivot_longer(
      cols = c(hilic_mean_intensity, hilic_num_features, c18_mean_intensity, c18_num_features),
      names_to = c("Mode", "metric"),
      names_pattern = "(.+)_(mean_intensity|num_features)",
      values_to = "value"
    ) %>%
    mutate(
      Mode = toupper(Mode),
      tx = factor(tx, levels = c("GS", "Ice", "RT")),
      time_numeric = ifelse(tx == "GS", 0, as.numeric(time))
    )
}

#+ 2a.1.2: Calculate summary statistics by treatment, time, and mode
{
  time_course_summary <- time_course_long %>%
    group_by(tx, time_numeric, Mode, metric) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

#+ 2a.1.3: Create plot for number of features over time
{
  # Filter for num_features metric
  features_data <- time_course_summary %>%
    filter(metric == "num_features")
  
  # Color mapping
  tx_colors <- c("GS" = "gray50", "Ice" = "#0958b3ff", "RT" = "#5d0f10")
  
  plot_features_time <- ggplot(features_data, 
                               aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    scale_color_manual(values = tx_colors, labels = c("GS (Baseline)", "Ice", "RT")) +
    scale_x_continuous(
      breaks = c(0, 15, 30, 60, 90),
      labels = c("0", "15", "30", "60", "90")
    ) +
    facet_wrap(~Mode, scales = "free_y", ncol = 2) +
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
      strip.text = element_text(size = 9, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold")
    ) +
    labs(
      title = "Feature Detection Over Time by Storage Condition",
      x = "Time (minutes)",
      y = "Number of Detected Features"
    )
  
  print(plot_features_time)
  
  # Save figure
  print_to_png(plot_features_time, "time_course_num_features.png", 
               width = 8, height = 4.5, output_dir = "Outputs/Figures", auto_open = TRUE)
}

#+ 2a.1.4: Create plot for mean intensity over time
{
  # Filter for mean_intensity metric
  intensity_data <- time_course_summary %>%
    filter(metric == "mean_intensity")
  
  plot_intensity_time <- ggplot(intensity_data, 
                                aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    scale_color_manual(values = tx_colors, labels = c("GS (Baseline)", "Ice", "RT")) +
    scale_x_continuous(
      breaks = c(0, 15, 30, 60, 90),
      labels = c("0", "15", "30", "60", "90")
    ) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~Mode, scales = "free_y", ncol = 2) +
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
      strip.text = element_text(size = 9, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold")
    ) +
    labs(
      title = "Mean Feature Intensity Over Time by Storage Condition",
      x = "Time (minutes)",
      y = "Mean Intensity (Raw, non-zero features)"
    )
  
  print(plot_intensity_time)
  
  # Save figure
  print_to_png(plot_intensity_time, "time_course_mean_intensity.png", 
               width = 8, height = 4.5, output_dir = "Outputs/Figures", auto_open = TRUE)
}

