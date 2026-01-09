#' Plot Metabolite Trajectory from Gold Standard
#'
#' Creates publication-style line plots showing fold change trajectories over time
#' for Ice vs RT storage conditions relative to patient-specific gold standard baseline.
#'
#' @param deviation_summary Tibble with columns: identified_name, feature, tx, time, mean_delta
#' @param identified_name Character string of metabolite name to plot
#' @param base_family Font family for plots (default: "Arial")
#' @param text_scale Scaling factor for all text elements (default: 1.0)
#' @param y_min Minimum value for y-axis (default: NULL, auto-calculated)
#' @param y_max Maximum value for y-axis (default: NULL, auto-calculated)
#'
#' @return ggplot object
#'
#' @export
plot_trajectory <- function(deviation_summary,
                            identified_name,
                            feature = NULL,
                            base_family = "Arial",
                            text_scale = 1.0,
                            y_min = NULL,
                            y_max = NULL) {
  
  library(dplyr)
  library(ggplot2)
  
  # Filter data for selected metabolite by identified_name and optionally by feature
  plot_data <- deviation_summary %>%
    filter(identified_name == !!identified_name)
  
  # If feature specified, filter further to that specific feature
  if (!is.null(feature)) {
    plot_data <- plot_data %>% filter(feature == !!feature)
  }
  
  plot_data <- plot_data %>%
    mutate(
      time_numeric = as.numeric(as.character(time)),
      tx = factor(tx, levels = c("Ice", "RT"))
    )
  
  # Check if data exists
  if (nrow(plot_data) == 0) {
    stop(paste("No data found for identified_name:", identified_name))
  }
  
  # Add baseline point at time = 0 with delta = 0 for both treatments (GS baseline)
  baseline_data <- tibble(
    identified_name = identified_name,
    feature = unique(plot_data$feature)[1],
    tx = factor(c("Ice", "RT"), levels = c("Ice", "RT")),
    time = factor("0", levels = c("0", "15", "30", "60", "90")),
    mean_delta = 0,
    time_numeric = 0
  )
  
  # Combine baseline with plot data
  plot_data <- bind_rows(baseline_data, plot_data) %>%
    arrange(tx, time_numeric)
  
  # Color mapping (matching plot_num_features style)
  tx_colors <- c("Ice" = "#0958b3ff", "RT" = "#5d0f10")  # Dark blue and dark red
  
  # Publication-style theme (matching plot_num_features)
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
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9 * text_scale, face = "bold", color = "black")
    )
  
  # Calculate y-axis limits with 4-tick structure (matching plot_num_features)
  y_data_min <- min(plot_data$mean_delta, na.rm = TRUE)
  y_data_max <- max(plot_data$mean_delta, na.rm = TRUE)
  
  if (is.null(y_min) || is.null(y_max)) {
    # Find range and extend symmetrically around zero
    y_range <- max(abs(y_data_min), abs(y_data_max))
    
    # Round to nice number
    magnitude <- floor(log10(y_range))
    base_value <- ceiling(y_range / (10^magnitude)) * (10^magnitude)
    
    if (is.null(y_min)) y_min <- -base_value
    if (is.null(y_max)) y_max <- base_value
  }
  
  # Create 4-tick breaks
  tick_spacing <- (y_max - y_min) / 3
  y_breaks <- c(y_min, y_min + tick_spacing, y_min + 2 * tick_spacing, y_max)
  y_limit_lower <- y_min * 1.017
  y_limit_upper <- y_max * 1.017
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = time_numeric, y = mean_delta, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.75) +
    geom_point(size = 3, shape = 19) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_color_manual(values = tx_colors, labels = c("Ice", "RT")) +
    scale_x_continuous(
      breaks = c(0, 15, 30, 60, 90),
      labels = c("0", "15", "30", "60", "90")
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = c(y_limit_lower, y_limit_upper)
    ) +
    labs(
      title = identified_name,
      x = "Time (minutes)",
      y = expression(bold("Log"[2]*" Fold Change from Baseline"))
    ) +
    theme_pub
  
  return(p)
}
