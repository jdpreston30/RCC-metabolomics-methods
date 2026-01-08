#- 1.1.2: Count non-zero features per ratio/tissue/platform
feature_counts <- ratio_uft_annotated %>%
  pivot_longer(
    cols = starts_with(c("HILIC_", "C18_")),
    names_to = "feature",
    values_to = "intensity"
  ) %>%
  mutate(platform = if_else(str_starts(feature, "HILIC_"), "HILIC", "C18")) %>%
  filter(intensity != 0) %>%
  group_by(ratio, tissue_type, platform) %>%
  summarise(
    total_features = n(),
    .groups = "drop"
  )
#- 1.1.3: Calculate mean intensity per ratio/tissue/platform
mean_intensities <- ratio_uft_annotated %>%
  pivot_longer(
    cols = starts_with(c("HILIC_", "C18_")),
    names_to = "feature",
    values_to = "intensity"
  ) %>%
  mutate(platform = if_else(str_starts(feature, "HILIC_"), "HILIC", "C18")) %>%
  filter(intensity != 0) %>%
  group_by(ratio, tissue_type, platform) %>%
  summarise(
    mean_intensity = mean(intensity),
    .groups = "drop"
  )
#- 1.1.4: Display results
feature_counts
mean_intensities
#+ 1.2: Statistical analysis of ratio performance
#- 1.2.1: ANOVA for feature counts
feature_anova <- aov(total_features ~ ratio + tissue_type + platform, data = feature_counts)
summary(feature_anova)
#- 1.2.2: Post-hoc pairwise comparisons for ratios
emmeans_features <- emmeans(feature_anova, ~ratio)
pairs(emmeans_features)
#- 1.2.3: Calculate coefficient of variation per ratio
cv_summary <- feature_counts %>%
  group_by(ratio) %>%
  summarise(
    mean_features = mean(total_features),
    sd_features = sd(total_features),
    cv_features = (sd_features / mean_features) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_features))
cv_summary
#- 1.2.4: Tissue type independence validation
tissue_feature_test <- t.test(total_features ~ tissue_type, data = feature_counts)
tissue_intensity_test <- t.test(mean_intensity ~ tissue_type, data = mean_intensities)
cat("Tissue type effect on features: p =", round(tissue_feature_test$p.value, 4), "\n")
cat("Tissue type effect on intensity: p =", round(tissue_intensity_test$p.value, 4), "\n")
#- 1.3.1: Tissue type comparison validation plots
p_tissue_features <- ggplot(feature_counts, aes(x = tissue_type, y = total_features, fill = tissue_type)) +
  geom_boxplot() +
  facet_wrap(~platform) +
  scale_fill_manual(values = c("Normal" = "#5d0f10", "Tumor" = "#dc9c52")) +
  labs(
    title = "Feature Detection: No Tissue Type Effect",
    subtitle = paste0("t-test p = ", round(tissue_feature_test$p.value, 3)),
    y = "Total Features Detected",
    x = "Tissue Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p_tissue_features
p_tissue_intensity <- ggplot(mean_intensities, aes(x = tissue_type, y = mean_intensity, fill = tissue_type)) +
  geom_boxplot() +
  facet_wrap(~platform, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#5d0f10", "Tumor" = "#dc9c52")) +
  labs(
    title = "Mean Intensity: No Tissue Type Effect",
    subtitle = paste0("t-test p = ", round(tissue_intensity_test$p.value, 3)),
    y = "Mean Intensity",
    x = "Tissue Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p_tissue_intensity
#- 1.3.2: Feature count comparison across ratios
p_features <- ggplot(feature_counts, aes(x = ratio, y = total_features, fill = tissue_type)) +
  geom_col(position = "dodge") +
  facet_wrap(~platform) +
  scale_fill_manual(values = c("Normal" = "#5d0f10", "Tumor" = "#dc9c52")) +
  labs(
    title = "Feature Detection by Ratio",
    y = "Total Features Detected",
    fill = "Tissue Type"
  ) +
  theme_minimal()
p_features
#- 1.3.3: Mean intensity comparison
p_intensity <- ggplot(mean_intensities, aes(x = ratio, y = mean_intensity, color = tissue_type)) +
  geom_point(size = 3) +
  geom_line(aes(group = tissue_type)) +
  facet_wrap(~platform, scales = "free_y") +
  scale_color_manual(values = c("Normal" = "#5d0f10", "Tumor" = "#dc9c52")) +
  labs(
    title = "Signal Intensity by Ratio",
    y = "Mean Intensity",
    color = "Tissue Type"
  ) +
  theme_minimal()
p_intensity
#- 1.3.4: Tissue-averaged correlation plot
combined_avg <- feature_counts %>%
  left_join(mean_intensities, by = c("ratio", "tissue_type", "platform")) %>%
  group_by(ratio, platform) %>%
  summarise(
    avg_features = mean(total_features),
    avg_intensity = mean(mean_intensity),
    .groups = "drop"
  )
p_correlation <- ggplot(combined_avg, aes(x = avg_features, y = avg_intensity, color = ratio)) +
  geom_point(size = 4) +
  geom_text(aes(label = ratio), vjust = -1, size = 3) +
  facet_wrap(~platform, scales = "free") +
  labs(
    title = "Feature Count vs Signal Intensity (Tissue-Averaged)",
    subtitle = "No trade-off: optimal ratios maximize both metrics",
    x = "Average Features Detected",
    y = "Average Mean Intensity"
  ) +
  theme_minimal()
p_correlation
#' Plot Number of Detected Features by Ratio
#'
#' Creates publication-style plots showing the number of detected features
#' across different ratios, with separate displays for normal vs tumor tissue.
#' Displays individual samples as points with bar overlays for group means.
#'
#' @param feature_counts Tibble with columns: ratio, tissue_type, platform, total_features
#' @param base_family Font family for plots (default: "Arial")
#' @param text_scale Scaling factor for all text elements (default: 1.0)
#' @param y_min Minimum value for y-axis (default: 0)
#' @param y_max Maximum value for y-axis (default: NULL, auto-calculated)
#'
#' @return ggplot object
#'
#' @export
plot_num_features <- function(feature_counts,
                              base_family = "Arial",
                              text_scale = 1.0,
                              y_min = 0,
                              y_max = NULL) {
  library(dplyr)
  library(ggplot2)
  library(forcats)

  # Define ratio order
  ratio_order <- c("15:2", "15:1", "20:1", "30:1", "60:1")

  # Prepare data with ordered factors
  plot_data <- feature_counts %>%
    mutate(
      ratio = factor(ratio, levels = ratio_order),
      tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")),
      platform = factor(platform, levels = c("HILIC", "C18"))
    ) %>%
    arrange(ratio, tissue_type)

  # Calculate summary statistics for bars
  summary_data <- plot_data %>%
    group_by(ratio, tissue_type) %>%
    summarise(
      mean_features = mean(total_features, na.rm = TRUE),
      .groups = "drop"
    )

  # Color mapping
  tissue_colors <- c("Normal" = "#dc9c52", "Tumor" = "#5d0f10")
  tissue_colors_light <- c("Normal" = "#f0d7a8", "Tumor" = "#8b5a5e")

  # Publication-style theme
  theme_pub <- theme_minimal(base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
      axis.line.y.left = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(size = 9 * text_scale, face = "bold", color = "black"),
      axis.text.x = element_text(size = 9 * text_scale, face = "bold", color = "black", angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_text(size = 9 * text_scale, face = "bold", color = "black"),
      axis.title = element_text(size = 10 * text_scale, face = "bold", color = "black"),
      axis.title.x = element_text(size = 10 * text_scale, face = "bold", color = "black", margin = margin(t = 5)),
      axis.title.y = element_text(size = 10 * text_scale, face = "bold", color = "black", margin = margin(r = 5)),
      plot.title = element_text(size = 11 * text_scale, face = "bold", hjust = 0.5, color = "black"),
      legend.position = "none",
      strip.text = element_text(size = 10 * text_scale, face = "bold", color = "black")
    )

  # Calculate y-axis limits with 4-tick structure
  if (is.null(y_max)) {
    y_max_data <- max(plot_data$total_features, na.rm = TRUE)
    top_tick_magnitude <- ceiling(log10(y_max_data))
    top_tick_base <- ceiling(y_max_data / (10^(top_tick_magnitude - 1))) * (10^(top_tick_magnitude - 1))

    if (top_tick_base < y_max_data * 1.1) {
      top_tick_base <- top_tick_base + (10^(top_tick_magnitude - 1))
    }

    tick_spacing <- (top_tick_base - y_min) / 3
    y_breaks <- c(y_min, y_min + tick_spacing, y_min + 2 * tick_spacing, top_tick_base)
    y_limit <- top_tick_base * 1.017
  } else {
    # Use provided y_max
    tick_spacing <- (y_max - y_min) / 3
    y_breaks <- c(y_min, y_min + tick_spacing, y_min + 2 * tick_spacing, y_max)
    y_limit <- y_max * 1.017
  }

  # Create the plot
  p <- ggplot(plot_data, aes(x = ratio, y = total_features, fill = tissue_type, color = tissue_type)) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.6,
      linewidth = 0.36,
      fatten = 1.2,
      alpha = 0.8,
      outlier.shape = NA
    ) +
    geom_point(
      aes(color = tissue_type),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
      size = 0.9,
      alpha = 1,
      shape = 16
    ) +
    scale_fill_manual(values = tissue_colors_light) +
    scale_color_manual(values = tissue_colors) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
      limits = c(y_min, y_limit),
      breaks = y_breaks,
      labels = function(x) {
        sapply(x, function(val) {
          if (is.na(val)) {
            return("")
          }
          if (val == 0) {
            return("0")
          }
          sci_notation <- sprintf("%.0e", val)
          cleaned <- gsub("e\\+0*", "E", sci_notation)
          cleaned <- gsub("E0+$", "", cleaned)
          return(cleaned)
        })
      }
    ) +
    facet_wrap(~platform, labeller = as_labeller(c("HILIC" = "HILIC+", "C18" = "C18-"))) +
    labs(
      x = "Ratio",
      y = "Total Features Detected",
      title = "Feature Detection by Extraction Ratio"
    ) +
    theme_pub

  return(p)
}
