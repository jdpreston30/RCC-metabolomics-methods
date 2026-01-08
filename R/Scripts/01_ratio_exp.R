#* 1: Ratio Experiment Analysis
#+ 1.1: Total Feature Counts
#- 1.1.1: Calculate feature counts per sample
sample_feature_counts <- ratio_uft_full_0s %>%
  pivot_longer(
    cols = starts_with(c("HILIC_", "C18_")),
    names_to = "feature",
    values_to = "intensity"
  ) %>%
  mutate(column_esi = if_else(str_starts(feature, "HILIC_"), "HILIC", "C18")) %>%
  filter(intensity != 0) %>%
  group_by(Sample_ID, ratio, name, tissue_type, column_esi) %>%
  summarise(
    total_features = n(),
    .groups = "drop"
  )
#- 1.1.2: ANOVA for feature counts per sample
feature_anova_samples <- aov(total_features ~ ratio * tissue_type * column_esi,
  data = sample_feature_counts
)
#- 1.1.3: Plot feature counts by ratio
num_features_ratio <- plot_num_features(sample_feature_counts,
                               y_min = 6.1e3, y_max = 8.5e3,
                               anova_results = feature_anova_samples, anova_y = 6800)
#+ 1.2: Mean Intensity (Detected Features Only)
#- 1.2.1: Calculate mean intensity per sample (excluding zeros)
sample_mean_intensity <- ratio_uft_full_0s %>%
  pivot_longer(
    cols = starts_with(c("HILIC_", "C18_")),
    names_to = "feature",
    values_to = "intensity"
  ) %>%
  mutate(column_esi = if_else(str_starts(feature, "HILIC_"), "HILIC", "C18")) %>%
  filter(intensity != 0) %>%
  group_by(Sample_ID, ratio, name, tissue_type, column_esi) %>%
  summarise(
    mean_intensity = mean(intensity),
    .groups = "drop"
  )
#- 1.2.2: ANOVA for mean intensity per sample
intensity_anova_samples <- aov(mean_intensity ~ ratio * tissue_type * column_esi,
  data = sample_mean_intensity
)
#- 1.2.3: Plot mean intensity by ratio
mean_intensity_ratio <- plot_mean_intensity(sample_mean_intensity,
                               anova_results = intensity_anova_samples, y_min = 2e7, y_max = 1.2e8, anova_y = 4.9e7)
#+ 1.3: Ratio Performance Analysis
#- 1.3.1: Merge datasets
correlation_data <- sample_feature_counts %>%
  left_join(sample_mean_intensity,
    by = c("Sample_ID", "ratio", "name", "tissue_type", "column_esi")
  )
#- 1.3.2: Calculate performance metrics by ratio and method
ratio_method_stats <- correlation_data %>%
  group_by(ratio, column_esi) %>%
  summarise(
    n_samples = n(),
    mean_features = mean(total_features),
    sd_features = sd(total_features),
    mean_intensity_val = mean(mean_intensity),
    sd_intensity_val = sd(mean_intensity),
    .groups = "drop"
  ) %>%
  mutate(
    features_cv = sd_features / mean_features * 100,
    intensity_cv = sd_intensity_val / mean_intensity_val * 100
  ) %>%
  rename(
    mean_intensity = mean_intensity_val,
    sd_intensity = sd_intensity_val
  )
#- 1.3.3: Rank ratios across both methods
ratio_performance <- ratio_method_stats %>%
  group_by(ratio) %>%
  summarise(
    avg_features = mean(mean_features),
    avg_intensity = mean(mean_intensity),
    max_features_cv = max(features_cv),
    max_intensity_cv = max(intensity_cv),
    intensity_stable = all(intensity_cv < 20),
    .groups = "drop"
  ) %>%
  mutate(
    features_score = (avg_features - min(avg_features)) / (max(avg_features) - min(avg_features)) * 100,
    intensity_score = (avg_intensity - min(avg_intensity)) / (max(avg_intensity) - min(avg_intensity)) * 100,
    cv_penalty = max_intensity_cv,
    composite_score = (features_score + intensity_score) / 2 - cv_penalty
  ) %>%
  arrange(desc(composite_score))
#- 1.3.4: Format performance table for figure
ratio_performance_table <- ratio_performance %>%
  mutate(
    ratio = factor(ratio, levels = c("15:2", "15:1", "20:1", "30:1", "60:1"))
  ) %>%
  arrange(ratio) %>%
  transmute(
    Ratio = ratio,
    Features = round(avg_features),
    Intensity = sprintf("%.1e", avg_intensity) %>% 
      gsub("e\\+0*", "E", .) %>% 
      gsub("E0+$", "", .),
    `CV (Features)` = sprintf("%.1f%%", max_features_cv),
    `CV (Intensity)` = sprintf("%.1f%%", max_intensity_cv)
  )
#- 1.3.5: Create table grob for plotting
ratio_table_grob <- gridExtra::tableGrob(
  ratio_performance_table,
  rows = NULL,
  theme = gridExtra::ttheme_default(
    base_family = "Arial",
    core = list(
      fg_params = list(cex = 0.55, fontface = "plain"),
      bg_params = list(fill = "white", col = "black", lwd = 0.5)
    ),
    colhead = list(
      fg_params = list(cex = 0.55, fontface = "bold"),
      bg_params = list(fill = "gray90", col = "black", lwd = 0.5)
    )
  )
)
#_ Add thick outer border
ratio_table_grob <- gtable::gtable_add_grob(
  ratio_table_grob,
  grobs = grid::rectGrob(gp = grid::gpar(col = "black", lwd = 2, fill = NA)),
  t = 1, l = 1, b = nrow(ratio_table_grob), r = ncol(ratio_table_grob),
  z = Inf, name = "outer-border"
)
#! 20:1 was selected as optimal, balancing high feature detection (7496 features) and intensity (57M) with excellent reproducibility (CV <17%) while maximizing sample efficiency."