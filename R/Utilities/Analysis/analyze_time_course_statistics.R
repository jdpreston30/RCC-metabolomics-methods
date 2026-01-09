#' Statistical Analysis of Time Course Metabolomics Data
#'
#' Performs linear mixed-effects models to test treatment and time effects
#' on deviations from baseline, accounting for repeated measures within patients.
#' Calculates deviation from baseline (GS) for each patient, then tests whether
#' these deviations differ between treatments over time.
#'
#' @param metrics_data A tibble containing metrics data with columns:
#'   - pt_tissue_ID: Patient tissue identifier
#'   - tx: Treatment factor (expects "GS", "Ice", "RT")
#'   - time: Time factor
#'   - Metric columns (e.g., combined_num_features, combined_mean_intensity)
#' @param metrics Character vector of metric column names to analyze
#'   (default: c("combined_num_features", "combined_mean_intensity"))
#' @param baseline_tx Character string for baseline treatment (default: "GS")
#' @param print_results Logical, whether to print statistical results (default: TRUE)
#'
#' @return A list containing for each metric:
#'   - deviation_data: Data with calculated deviations from baseline
#'   - model: The fitted lmer model object
#'   - anova: ANOVA results from the model
#'   - emmeans_tx: Estimated marginal means by treatment
#'   - emmeans_time: Estimated marginal means by time
#'   - emmeans_interaction: Estimated marginal means by treatment × time
#'
#' @export
analyze_time_course_statistics <- function(metrics_data,
                                           metrics = c("combined_num_features", "combined_mean_intensity"),
                                           baseline_tx = "GS",
                                           print_results = TRUE) {
  
  # Ensure required packages are loaded
  require(lme4)
  require(lmerTest)
  require(emmeans)
  
  # Initialize results list
  results <- list()
  
  # Analyze each metric
  for (metric in metrics) {
    
    # Check if metric exists in data
    if (!metric %in% colnames(metrics_data)) {
      warning(paste("Metric", metric, "not found in data. Skipping."))
      next
    }
    
    if (print_results) {
      cat("\n")
      cat("================================================================================\n")
      cat(paste("Statistical Analysis:", metric, "(Deviation from Baseline)\n"))
      cat("================================================================================\n")
    }
    
    # Get baseline values for each patient
    baseline_data <- metrics_data %>%
      filter(tx == baseline_tx) %>%
      select(pt_tissue_ID, baseline_value = !!sym(metric))
    
    # Calculate deviations from baseline
    deviation_data <- metrics_data %>%
      filter(tx != baseline_tx) %>%
      left_join(baseline_data, by = "pt_tissue_ID") %>%
      mutate(
        deviation = !!sym(metric) - baseline_value,
        tx = factor(tx),
        time = factor(time)
      )
    
    # Create model formula for deviation
    model_formula <- as.formula("deviation ~ tx * time + (1|pt_tissue_ID)")
    
    # Fit mixed-effects model
    model <- lmer(model_formula, data = deviation_data)
    
    # Get ANOVA results
    anova_results <- anova(model, type = "III")
    
    if (print_results) {
      cat("\nANOVA Results (Type III) - Testing Deviations from Baseline:\n")
      print(anova_results)
    }
    
    # Estimated marginal means
    emmeans_tx <- emmeans(model, ~ tx)
    emmeans_time <- emmeans(model, ~ time)
    emmeans_interaction <- emmeans(model, ~ tx * time)
    
    if (print_results) {
      cat("\n--- Main Effect: Treatment (Mean Deviation from Baseline) ---\n")
      print(emmeans_tx)
      cat("\nPairwise Comparisons (Treatment):\n")
      print(pairs(emmeans_tx))
      
      cat("\n--- Main Effect: Time (Mean Deviation from Baseline) ---\n")
      print(emmeans_time)
      
      cat("\n--- Interaction: Treatment × Time (Deviations) ---\n")
      print(emmeans_interaction)
      cat("\nSimple Effects (Treatment Difference at each Time):\n")
      print(pairs(emmeans_interaction, by = "time"))
    }
    
    # Store results
    results[[metric]] <- list(
      deviation_data = deviation_data,
      model = model,
      anova = anova_results,
      emmeans_tx = emmeans_tx,
      emmeans_time = emmeans_time,
      emmeans_interaction = emmeans_interaction,
      pairwise_tx = pairs(emmeans_tx),
      simple_effects = pairs(emmeans_interaction, by = "time")
    )
  }
  
  return(results)
}
