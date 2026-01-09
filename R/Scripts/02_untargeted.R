main_uft_20p_tumors <- main_uft_20p %>%
  filter(tissue_type == "Tumor") |>
  select(-tissue_type) |>
  filter(tx != "CCM")

#* 2.1: Global Metabolic Deviation Analysis
#+ 2.1.1: Extract GS baselines and calculate deviations from baseline
{
  # Get GS baseline samples (one per patient at time = 0)
  GS_baseline_global <- main_uft_20p_tumors %>%
    filter(tx == "GS") %>%
    select(pt_tissue_ID, starts_with("HILIC"), starts_with("C18"))
  
  # Get Ice/RT samples only
  global_ice_rt <- main_uft_20p_tumors %>%
    filter(tx != "GS")
  
  # Get feature columns
  feature_cols <- colnames(main_uft_20p_tumors)[str_detect(colnames(main_uft_20p_tumors), "^(HILIC|C18)_")]
  
  # Calculate deviations from GS for each patient - VECTORIZED APPROACH
  # Convert to matrices for fast computation
  gs_matrix <- GS_baseline_global %>%
    select(all_of(feature_cols)) %>%
    as.matrix()
  
  sample_matrix <- global_ice_rt %>%
    select(all_of(feature_cols)) %>%
    as.matrix()
  
  # Match patient IDs to get correct GS baseline for each sample
  patient_match <- match(global_ice_rt$pt_tissue_ID, GS_baseline_global$pt_tissue_ID)
  
  # Subtract GS from samples (vectorized across all features at once)
  delta_matrix <- sample_matrix - gs_matrix[patient_match, ]
  
  # Convert back to tibble and add delta column names
  delta_df <- as_tibble(delta_matrix)
  colnames(delta_df) <- paste0(feature_cols, "_delta")
  
  # Combine with metadata
  global_with_deviations <- bind_cols(
    global_ice_rt %>% select(Sample_ID, pt_tissue_ID, tx, time),
    delta_df
  )
  
  # Convert to long format for analysis
  delta_cols <- paste0(feature_cols, "_delta")
  
  global_deviations_long <- global_with_deviations %>%
    select(Sample_ID, pt_tissue_ID, tx, time, all_of(delta_cols)) %>%
    pivot_longer(
      cols = all_of(delta_cols),
      names_to = "feature_delta",
      values_to = "delta_from_GS"
    ) %>%
    mutate(
      feature = str_remove(feature_delta, "_delta"),
      tx = factor(tx, levels = c("Ice", "RT")),
      time = factor(time, levels = c("15", "30", "60", "90"))
    ) %>%
    select(-feature_delta)
  
  # Calculate mean deviation per sample (average across all features)
  global_mean_deviation <- global_deviations_long %>%
    group_by(Sample_ID, pt_tissue_ID, tx, time) %>%
    summarise(
      mean_delta = mean(delta_from_GS, na.rm = TRUE),
      median_delta = median(delta_from_GS, na.rm = TRUE),
      sd_delta = sd(delta_from_GS, na.rm = TRUE),
      n_features = sum(!is.na(delta_from_GS)),
      .groups = "drop"
    )
}

#+ 2.1.2: Statistical analysis of global deviation
{
  library(lmerTest)
  
  # 2-way repeated measures: treatment × time with patient random effect
  global_model <- lmer(mean_delta ~ tx * time + (1|pt_tissue_ID), data = global_mean_deviation)
  global_anova <- anova(global_model, type = "III")
  
  print("Global Metabolic Deviation Analysis:")
  print(global_anova)
  
  # Get mean deviation by treatment and time
  global_summary <- global_mean_deviation %>%
    group_by(tx, time) %>%
    summarise(
      mean = mean(mean_delta, na.rm = TRUE),
      se = sd(mean_delta, na.rm = TRUE) / sqrt(n()),
      sd = sd(mean_delta, na.rm = TRUE),
      .groups = "drop"
    )
}

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
