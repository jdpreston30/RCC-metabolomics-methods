#' Plot Mean Intensity of Detected Features by Ratio
#'
#' Creates publication-style plots showing mean intensity of detected features
#' (excluding zeros/non-detects) across different ratios, with separate displays
#' for normal vs tumor tissue. Displays individual samples as points with bar
#' overlays for group means.
#'
#' @param intensity_data Tibble with columns: ratio, tissue_type, column_esi, mean_intensity
#' @param base_family Font family for plots (default: "Arial")
#' @param text_scale Scaling factor for all text elements (default: 1.0)
#' @param y_min Minimum value for y-axis (default: 0)
#' @param y_max Maximum value for y-axis (default: NULL, auto-calculated)
#' @param anova_results Optional ANOVA results object from aov() (default: NULL, no stats displayed)
#' @param anova_y Vertical position for ANOVA annotation (default: Inf, top of plot)
#'
#' @return ggplot object
#'
#' @export
plot_mean_intensity <- function(intensity_data,
                                base_family = "Arial",
                                text_scale = 1.0,
                                y_min = 0,
                                y_max = NULL,
                                anova_results = NULL,
                                anova_y = Inf) {
  
  library(dplyr)
  library(ggplot2)
  library(forcats)
  
  # Define ratio order
  ratio_order <- c("15:2", "15:1", "20:1", "30:1", "60:1")
  
  # Prepare data with ordered factors
  plot_data <- intensity_data %>%
    mutate(
      ratio = factor(ratio, levels = ratio_order),
      tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")),
      column_esi = factor(column_esi, levels = c("HILIC", "C18"))
    ) %>%
    arrange(ratio, tissue_type)
  
  # Calculate summary statistics for bars
  summary_data <- plot_data %>%
    group_by(ratio, tissue_type) %>%
    summarise(
      mean_of_means = mean(mean_intensity, na.rm = TRUE),
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
      axis.text = element_text(size = 8 * text_scale, face = "bold", color = "black"),
      axis.text.x = element_text(size = 8 * text_scale, face = "bold", color = "black", angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_text(size = 8 * text_scale, face = "bold", color = "black"),
      axis.title = element_text(size = 10 * text_scale, face = "bold", color = "black"),
      axis.title.x = element_text(size = 10 * text_scale, face = "bold", color = "black", margin = margin(t = 5)),
      axis.title.y = element_text(size = 10 * text_scale, face = "bold", color = "black", margin = margin(r = 5)),
      plot.title = element_text(size = 11 * text_scale, face = "bold", hjust = 0.5, color = "black"),
      legend.position = "none",
      strip.text = element_text(size = 10 * text_scale, face = "bold", color = "black")
    )
  
  # Calculate y-axis limits with 4-tick structure
  if (is.null(y_max)) {
    y_max_data <- max(plot_data$mean_intensity, na.rm = TRUE)
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
  p <- ggplot(plot_data, aes(x = ratio, y = mean_intensity, fill = tissue_type, color = tissue_type)) +
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
          sci_notation <- sprintf("%.1e", val)
          cleaned <- gsub("e\\+0*", "E", sci_notation)
          cleaned <- gsub("E0+$", "", cleaned)
          return(cleaned)
        })
      }
    ) +
    facet_wrap(~ column_esi, labeller = as_labeller(c("HILIC" = "HILIC+", "C18" = "C18-"))) +
    labs(
      x = "Ratio",
      y = "Mean Intensity (Detected Features)"
    ) +
    theme_pub
  
  # Add ANOVA annotation if provided (only on HILIC+ facet, bottom left)
  if (!is.null(anova_results)) {
    # Extract and format ANOVA results
    anova_summary <- summary(anova_results)[[1]]
    
    # Create lookup table for term display names and order
    term_lookup <- tibble::tribble(
      ~term, ~display, ~order,
      "tissue_type", "Tissue", 1,
      "ratio", "Ratio", 2,
      "column_esi", "Method", 3,
      "ratio:tissue_type", "Tissue \u00d7 Ratio", 4,
      "tissue_type:column_esi", "Tissue \u00d7 Method", 5,
      "ratio:column_esi", "Ratio \u00d7 Method", 6,
      "ratio:tissue_type:column_esi", "Tissue \u00d7 Ratio \u00d7 Method", 7
    )
    
    # Format annotation
    anova_annotation <- tibble::tibble(
      term = stringr::str_trim(rownames(anova_summary)),
      p_value = anova_summary[, "Pr(>F)"]
    ) %>%
      filter(term != "Residuals", !is.na(p_value)) %>%
      left_join(term_lookup, by = "term") %>%
      mutate(
        p_text = if_else(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value)),
        line = paste0(display, ": ", p_text)
      ) %>%
      arrange(order)
    
    anova_text <- paste(anova_annotation$line, collapse = "\n")
    
    # Add annotation to HILIC+ facet (left panel), bottom left
    annotation_data <- data.frame(
      column_esi = factor("HILIC", levels = c("HILIC", "C18")),
      x = 0.7,
      y = anova_y,
      label = anova_text
    )
    
    p <- p + 
      geom_text(data = annotation_data,
                aes(x = x, y = y, label = label),
                hjust = 0, vjust = 1,
                size = 1.6 * text_scale,
                family = base_family,
                color = "black",
                lineheight = 0.9,
                inherit.aes = FALSE)
  }
  
  return(p)
}
