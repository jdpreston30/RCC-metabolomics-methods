#' Create Time Course Plot with Baseline Connection
#'
#' Creates a time course plot showing metabolomics metrics over time with
#' baseline connections from GS treatment to Ice/RT treatments. Supports
#' faceting by mode (HILIC/C18) and customizable styling.
#'
#' @param summary_data A tibble containing summary statistics with columns:
#'   - tx: Treatment factor
#'   - time_numeric: Numeric time values
#'   - Mode: Mode factor (for faceting, optional)
#'   - metric: Metric name
#'   - mean_value: Mean values to plot
#'   - se: Standard error (optional, for error bars)
#' @param target_metric Character string specifying which metric to plot
#' @param baseline_treatment Character string for baseline treatment (default: "GS")
#' @param treatment_colors Named vector of colors for treatments (default: c("Ice" = "#0958b3ff", "RT" = "#5d0f10"))
#' @param time_breaks Numeric vector of time breaks for x-axis (default: c(0, 15, 30, 60, 90))
#' @param plot_title Character string for plot title
#' @param y_label Character string for y-axis label
#' @param x_label Character string for x-axis label (default: "Time (minutes)")
#' @param use_facets Logical, whether to facet by Mode (default: TRUE)
#' @param facet_scales Character string for facet scales (default: "free_y")
#' @param add_error_bars Logical, whether to add error bars (default: FALSE)
#' @param error_bar_type Character string for error bar type: "se" for SEM or "sd" for SD (default: "se")
#' @param loess_span Numeric value for loess span (default: 0.75)
#' @param y_axis_format Character string for y-axis formatting: "comma", "scientific", or NULL (default: NULL)
#' @param y_min Minimum value for y-axis (default: NULL, auto-calculated)
#' @param y_max Maximum value for y-axis (default: NULL, auto-calculated)
#' @param anova_results Optional ANOVA results from lmerTest (default: NULL)
#' @param anova_x Horizontal position for ANOVA annotation (default: 0, left side)
#' @param anova_y Vertical position for ANOVA annotation (default: Inf, top of plot)
#'
#' @return A ggplot object
#'
#' @export
plot_time_course_with_baseline <- function(summary_data,
                                           target_metric,
                                           baseline_treatment = "GS",
                                           treatment_colors = c("Ice" = "#0958b3ff", "RT" = "#5d0f10"),
                                           time_breaks = c(0, 15, 30, 60, 90),
                                           plot_title = "Metabolomics Time Course",
                                           y_label = "Value",
                                           x_label = "Time (minutes)",
                                           use_facets = TRUE,
                                           facet_scales = "free_y",
                                           add_error_bars = FALSE,
                                           error_bar_type = "se",
                                           loess_span = 0.75,
                                           y_axis_format = NULL,
                                           y_min = NULL,
                                           y_max = NULL,
                                           anova_results = NULL,
                                           anova_x = 0,
                                           anova_y = Inf) {
  
  # Filter for target metric
  filtered_data <- summary_data %>%
    filter(metric == target_metric)
  
  if (use_facets) {
    # Get baseline values by Mode
    baseline_data <- filtered_data %>%
      filter(tx == baseline_treatment) %>%
      select(Mode, metric, baseline_value = mean_value)
    
    # Add baseline point to non-baseline treatments
    treatment_data <- filtered_data %>%
      filter(tx != baseline_treatment) %>%
      left_join(baseline_data, by = c("Mode", "metric"))
    
    # Create baseline rows for each treatment
    baseline_points <- treatment_data %>%
      distinct(tx, Mode, metric, baseline_value) %>%
      mutate(time_numeric = 0, mean_value = baseline_value, se = 0) %>%
      select(-baseline_value)
    
    # Combine baseline with treatment data
    plot_data <- bind_rows(baseline_points, treatment_data %>% select(-baseline_value))
    
  } else {
    # Get baseline value (single value)
    baseline_value <- filtered_data %>%
      filter(tx == baseline_treatment) %>%
      pull(mean_value)
    
    # Create baseline rows for each treatment
    treatment_levels <- unique(filtered_data$tx[filtered_data$tx != baseline_treatment])
    baseline_points <- tibble(
      tx = factor(treatment_levels, levels = levels(filtered_data$tx)),
      time_numeric = 0,
      metric = target_metric,
      mean_value = baseline_value,
      se = 0
    )
    
    # Combine baseline with treatment data
    plot_data <- bind_rows(
      baseline_points,
      filtered_data %>% filter(tx != baseline_treatment)
    )
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = time_numeric, y = mean_value, color = tx, group = tx)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = loess_span) +
    geom_point(size = 3, shape = 19) +
    scale_color_manual(values = treatment_colors, labels = names(treatment_colors)) +
    scale_x_continuous(
      breaks = time_breaks,
      labels = as.character(time_breaks)
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
      axis.line.y.left = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(size = 8, face = "bold", color = "black"),
      axis.text.x = element_text(size = 8, face = "bold", color = "black", angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_text(size = 8, face = "bold", color = "black"),
      axis.title = element_text(size = 10, face = "bold", color = "black"),
      axis.title.x = element_text(size = 10, face = "bold", color = "black", margin = margin(t = 5)),
      axis.title.y = element_text(size = 10, face = "bold", color = "black", margin = margin(r = 5)),
      plot.title = element_blank(),
      legend.position = c(0.95, 0.05),
      legend.justification = c("right", "bottom"),
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold", color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(2, 2, 2, 2),
      legend.spacing.y = unit(0.1, "lines"),
      legend.key.size = unit(0.8, "lines"),
      legend.key.height = unit(0.8, "lines")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, linewidth = 0.8))) +
    labs(
      x = x_label,
      y = y_label
    )
  
  # Add faceting if requested
  if (use_facets) {
    p <- p +
      facet_wrap(~Mode, scales = facet_scales, ncol = 2) +
      theme(
        strip.text = element_text(size = 10, face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5)
      )
  }
  
  # Add error bars if requested
  if (add_error_bars) {
    # Determine which error column to use
    error_col <- if(error_bar_type == "sd") "sd" else "se"
    
    # Check if the column exists in plot_data
    if (error_col %in% colnames(plot_data)) {
      p <- p +
        geom_errorbar(
          aes(ymin = mean_value - .data[[error_col]], ymax = mean_value + .data[[error_col]]),
          width = 2,
          linewidth = 0.5
        )
    } else {
      warning(paste0("Error bar type '", error_bar_type, "' requested but '", error_col, "' column not found in data."))
    }
  }
  
  # Set y-axis limits if specified
  if (!is.null(y_min) || !is.null(y_max)) {
    # Get current y limits from data if not fully specified
    if (is.null(y_min)) {
      y_min <- min(plot_data$mean_value, na.rm = TRUE) * 0.95
    }
    if (is.null(y_max)) {
      y_max <- max(plot_data$mean_value, na.rm = TRUE) * 1.05
    }
  }
  
  # Add y-axis formatting if requested
  if (!is.null(y_axis_format)) {
    if (y_axis_format == "comma") {
      if (!is.null(y_min) || !is.null(y_max)) {
        p <- p + scale_y_continuous(labels = scales::comma, limits = c(y_min, y_max))
      } else {
        p <- p + scale_y_continuous(labels = scales::comma)
      }
    } else if (y_axis_format == "scientific") {
      p <- p + scale_y_continuous(
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
        },
        limits = if (!is.null(y_min) || !is.null(y_max)) c(y_min, y_max) else NULL
      )
    }
  } else if (!is.null(y_min) || !is.null(y_max)) {
    # Apply limits even without formatting
    p <- p + scale_y_continuous(limits = c(y_min, y_max))
  }
  
  # Add ANOVA annotation if provided
  if (!is.null(anova_results)) {
    # Format ANOVA results
    term_lookup <- tibble::tribble(
      ~term, ~display, ~order,
      "tx", "Temperature", 1,
      "time", "Time", 2,
      "tx:time", "Temperature \u00d7 Time", 3
    )
    
    anova_annotation <- tibble::tibble(
      term = rownames(anova_results),
      p_value = anova_results[, "Pr(>F)"]
    ) %>%
      filter(!is.na(p_value)) %>%
      left_join(term_lookup, by = "term") %>%
      mutate(
        p_text = if_else(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value)),
        line = paste0(display, ": ", p_text)
      ) %>%
      arrange(order)
    
    anova_text <- paste(anova_annotation$line, collapse = "\n")
    
    # Add annotation as text
    p <- p + 
      annotate("text",
               x = anova_x, y = anova_y,
               label = anova_text,
               hjust = 0, vjust = 1,
               size = 1.6,
               family = "Arial",
               color = "black",
               lineheight = 0.9)
  }
  
  return(p)
}