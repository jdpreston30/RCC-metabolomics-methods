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
      n_samples = n(),
      .groups = "drop"
    )
  
  # DIAGNOSTIC: Check sample sizes per group
  cat("\n=== Sample Sizes by Treatment × Time ===\n")
  time_course_summary %>%
    filter(metric == "num_features", Mode == "HILIC") %>%
    select(tx, time_numeric, n_samples) %>%
    print(n = Inf)
  
  # DIAGNOSTIC: Check if TM11 is an outlier in feature detection
  cat("\n=== TM11 Feature Detection Compared to Other Patients ===\n")
  tm11_check <- time_course_metrics %>%
    filter(tx %in% c("Ice", "RT")) %>%
    group_by(tx, time) %>%
    summarise(
      mean_hilic_all = mean(hilic_num_features, na.rm = TRUE),
      mean_hilic_tm11 = mean(hilic_num_features[pt_tissue_ID == "TM11"], na.rm = TRUE),
      mean_hilic_others = mean(hilic_num_features[pt_tissue_ID != "TM11"], na.rm = TRUE),
      tm11_diff = mean_hilic_tm11 - mean_hilic_others,
      .groups = "drop"
    )
  print(tm11_check)
}

#+ 2a.1.3: Create plot for number of features over time
{
  # Filter for num_features metric
  features_data <- time_course_summary %>%
    filter(metric == "num_features")
  
  # Get GS baseline values at time=0 to connect Ice and RT lines
  gs_baseline <- features_data %>%
    filter(tx == "GS") %>%
    select(Mode, metric, baseline_value = mean_value)
  
  # Add baseline point to Ice and RT groups
  ice_rt_data <- features_data %>%
    filter(tx != "GS") %>%
    left_join(gs_baseline, by = c("Mode", "metric"))
  
  # Create baseline rows for Ice and RT
  baseline_points <- ice_rt_data %>%
    distinct(tx, Mode, metric, baseline_value) %>%
    mutate(time_numeric = 0, mean_value = baseline_value, se = 0) %>%
    select(-baseline_value)
  
  # Combine baseline with Ice/RT data
  plot_data_features <- bind_rows(baseline_points, ice_rt_data %>% select(-baseline_value))
  
  # Color mapping
  tx_colors <- c("Ice" = "#0958b3ff", "RT" = "#5d0f10")
  
  plot_features_time <- ggplot(plot_data_features, 
                               aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    scale_color_manual(values = tx_colors, labels = c("Ice", "RT")) +
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
  
  # Get GS baseline values at time=0 to connect Ice and RT lines
  gs_baseline_int <- intensity_data %>%
    filter(tx == "GS") %>%
    select(Mode, metric, baseline_value = mean_value)
  
  # Add baseline point to Ice and RT groups
  ice_rt_data_int <- intensity_data %>%
    filter(tx != "GS") %>%
    left_join(gs_baseline_int, by = c("Mode", "metric"))
  
  # Create baseline rows for Ice and RT
  baseline_points_int <- ice_rt_data_int %>%
    distinct(tx, Mode, metric, baseline_value) %>%
    mutate(time_numeric = 0, mean_value = baseline_value, se = 0) %>%
    select(-baseline_value)
  
  # Combine baseline with Ice/RT data
  plot_data_intensity <- bind_rows(baseline_points_int, ice_rt_data_int %>% select(-baseline_value))
  
  plot_intensity_time <- ggplot(plot_data_intensity, 
                                aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    scale_color_manual(values = tx_colors, labels = c("Ice", "RT")) +
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

#* 2a.2: Combined HILIC + C18 Analysis
#+ 2a.2.1: Calculate combined metrics (all features together)
{
  # Filter out CCM
  time_course_combined <- main_uft_full_0s_tumor %>%
    filter(tx != "CCM")
  
  # Get all feature columns (HILIC and C18 combined)
  all_feature_cols <- colnames(time_course_combined)[str_detect(colnames(time_course_combined), "^(HILIC|C18)_")]
  
  # Calculate metrics for each sample across ALL features
  time_course_combined_metrics <- time_course_combined %>%
    rowwise() %>%
    mutate(
      combined_mean_intensity = mean(c_across(all_of(all_feature_cols))[c_across(all_of(all_feature_cols)) > 0], na.rm = TRUE),
      combined_num_features = sum(c_across(all_of(all_feature_cols)) > 0, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Sample_ID, pt_tissue_ID, tx, time, combined_mean_intensity, combined_num_features) %>%
    mutate(
      tx = factor(tx, levels = c("GS", "Ice", "RT")),
      time_numeric = ifelse(tx == "GS", 0, as.numeric(time))
    )
}

#+ 2a.2.2: Calculate summary statistics for combined data
{
  time_course_combined_summary <- time_course_combined_metrics %>%
    pivot_longer(
      cols = c(combined_mean_intensity, combined_num_features),
      names_to = "metric",
      names_prefix = "combined_",
      values_to = "value"
    ) %>%
    group_by(tx, time_numeric, metric) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

#+ 2a.2.3: Create combined plot for number of features
{
  # Filter for num_features metric
  combined_features_data <- time_course_combined_summary %>%
    filter(metric == "num_features")
  
  # Get GS baseline value at time=0
  gs_baseline_combined <- combined_features_data %>%
    filter(tx == "GS") %>%
    pull(mean_value)
  
  # Create baseline rows for Ice and RT
  baseline_combined_features <- tibble(
    tx = factor(c("Ice", "RT"), levels = c("Ice", "RT")),
    time_numeric = 0,
    metric = "num_features",
    mean_value = gs_baseline_combined,
    se = 0
  )
  
  # Combine baseline with Ice/RT data
  plot_data_combined_features <- bind_rows(
    baseline_combined_features,
    combined_features_data %>% filter(tx != "GS")
  )
  
  # Color mapping
  tx_colors <- c("Ice" = "#0958b3ff", "RT" = "#5d0f10")
  
  plot_combined_features <- ggplot(plot_data_combined_features, 
                                   aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
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
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold")
    ) +
    labs(
      title = "Total Feature Detection Over Time (HILIC + C18)",
      x = "Time (minutes)",
      y = "Number of Detected Features"
    )
  
  print(plot_combined_features)
  
  # Save figure
  print_to_png(plot_combined_features, "time_course_combined_num_features.png", 
               width = 6, height = 4.5, output_dir = "Outputs/Figures", auto_open = TRUE)
}

#+ 2a.2.4: Create combined plot for mean intensity
{
  # Filter for mean_intensity metric
  combined_intensity_data <- time_course_combined_summary %>%
    filter(metric == "mean_intensity")
  
  # Get GS baseline value at time=0
  gs_baseline_combined_int <- combined_intensity_data %>%
    filter(tx == "GS") %>%
    pull(mean_value)
  
  # Create baseline rows for Ice and RT
  baseline_combined_intensity <- tibble(
    tx = factor(c("Ice", "RT"), levels = c("Ice", "RT")),
    time_numeric = 0,
    metric = "mean_intensity",
    mean_value = gs_baseline_combined_int,
    se = 0
  )
  
  # Combine baseline with Ice/RT data
  plot_data_combined_intensity <- bind_rows(
    baseline_combined_intensity,
    combined_intensity_data %>% filter(tx != "GS")
  )
  
  plot_combined_intensity <- ggplot(plot_data_combined_intensity, 
                                    aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    scale_color_manual(values = tx_colors, labels = c("Ice", "RT")) +
    scale_x_continuous(
      breaks = c(0, 15, 30, 60, 90),
      labels = c("0", "15", "30", "60", "90")
    ) +
    scale_y_continuous(labels = scales::comma) +
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
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold")
    ) +
    labs(
      title = "Mean Feature Intensity Over Time (HILIC + C18)",
      x = "Time (minutes)",
      y = "Mean Intensity (Raw, non-zero features)"
    )
  
  print(plot_combined_intensity)
  
  # Save figure
  print_to_png(plot_combined_intensity, "time_course_combined_mean_intensity.png", 
               width = 6, height = 4.5, output_dir = "Outputs/Figures", auto_open = TRUE)
}
