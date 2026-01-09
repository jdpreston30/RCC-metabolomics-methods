#* 2.1: Untargeted Features Over Time
#+ 2.0: Create data subsets
#- 2.0.1: Tumor samples
main_uft_full_0s_tumor <- main_uft_full_0s %>%
  filter(tissue_type == "Tumor") |>
  select(-tissue_type)
#- 2.0.2: Normal samples
main_uft_full_0s_normal <- main_uft_full_0s %>%
  filter(tissue_type == "Normal") |>
  select(-tissue_type)
#+ 2.1: Calculate TC metrics using utility function (Tumor)
#- 2.1.1: Run TC metrics with statistics
tc_tumor <- calculate_time_course_metrics(
  data = main_uft_full_0s_tumor,
  exclude_treatments = "CCM",
  tx_levels = c("GS", "Ice", "RT"),
  combine_modes = TRUE,
  run_diagnostics = TRUE,
  diagnostic_patient = "TM11",
  run_statistics = TRUE,
  baseline_tx = "GS",
  print_statistics = TRUE
)
#- 2.1.2: Extract results
tc_tumor_metrics <- tc_tumor$metrics_wide
tc_tumor_long <- tc_tumor$metrics_long
tc_tumor_summary <- tc_tumor$summary_stats
tc_tumor_stats <- tc_tumor$statistics
#- 2.1.3: Create plot of features over time
plot_features_time <- plot_time_course_with_baseline(
  summary_data = tc_tumor_summary,
  target_metric = "num_features",
  y_label = "Total Features Detected",
  use_facets = FALSE,
  loess_span = 0.8,
  anova_results = tc_tumor_stats$combined_num_features$anova,
  anova_x = 0,
  anova_y = 15000
)
#- 2.1.4: Create plot of mean intensity over time
plot_intensity_time <- plot_time_course_with_baseline(
  summary_data = tc_tumor_summary,
  target_metric = "mean_intensity",
  y_label = "Mean Intensity (Detected Features)",
  use_facets = FALSE,
  y_axis_format = "scientific",
  loess_span = 0.8,
  anova_results = tc_tumor_stats$combined_mean_intensity$anova,
  anova_x = 0,
  anova_y = Inf,
  y_min = 6e7,
  y_max = 7.5e7
)
#+ 2.2: Calculate TC metrics using utility function (Normal)
#- 2.2.1: Run TC metrics with statistics
tc_normal <- calculate_time_course_metrics(
  data = main_uft_full_0s_normal,
  exclude_treatments = "CCM",
  tx_levels = c("GS", "Ice", "RT"),
  combine_modes = TRUE,
  run_diagnostics = TRUE,
  diagnostic_patient = "TM11",
  run_statistics = TRUE,
  baseline_tx = "GS",
  print_statistics = TRUE
)
#- 2.2.2: Extract results
tc_normal_metrics <- tc_normal$metrics_wide
tc_normal_long <- tc_normal$metrics_long
tc_normal_summary <- tc_normal$summary_stats
tc_normal_stats <- tc_normal$statistics
#- 2.2.3: Create plot of features over time
plot_features_time_normal <- plot_time_course_with_baseline(
  summary_data = tc_normal_summary,
  target_metric = "num_features",
  y_label = "Total Features Detected",
  use_facets = FALSE,
  loess_span = 0.5,
  anova_results = tc_normal_stats$combined_num_features$anova,
  anova_x = 0,
  anova_y = Inf
)
#- 2.2.4: Create plot of mean intensity over time
plot_intensity_time_normal <- plot_time_course_with_baseline(
  summary_data = tc_normal_summary,
  target_metric = "mean_intensity",
  y_label = "Mean Intensity (Detected Features)",
  use_facets = FALSE,
  y_axis_format = "scientific",
  loess_span = 0.5,
  anova_results = tc_normal_stats$combined_mean_intensity$anova,
  anova_x = 0,
  anova_y = Inf,
  y_min = 5.8e7,
  y_max = 6.8e7
)