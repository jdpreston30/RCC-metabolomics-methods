#' Statistical Analysis of Global Metabolic Deviations
#'
#' Performs statistical analysis on global metabolic deviation data using
#' mixed-effects models to test treatment and time effects with patient
#' as a random effect.
#'
#' @param mean_deviation_data A tibble containing mean deviation data with columns:
#'   - pt_tissue_ID: Patient tissue identifier
#'   - tx: Treatment factor
#'   - time: Time factor  
#'   - mean_delta: Mean deviation values
#' @param response_var Character string specifying response variable column name
#'   (default: "mean_delta")
#' @param print_results Logical, whether to print ANOVA results (default: TRUE)
#'
#' @return A list containing:
#'   - model: The fitted lmer model object
#'   - anova: ANOVA results from the model
#'   - summary_stats: Summary statistics by treatment and time
#'
#' @export
analyze_global_deviations <- function(mean_deviation_data, 
                                      response_var = "mean_delta",
                                      print_results = TRUE) {
  
  # Ensure required packages are loaded
  require(lme4)
  require(lmerTest)
  
  # Create model formula dynamically
  model_formula <- as.formula(paste(response_var, "~ tx * time + (1|pt_tissue_ID)"))
  
  # Fit 2-way repeated measures: treatment × time with patient random effect
  model <- lmer(model_formula, data = mean_deviation_data)
  anova_results <- anova(model, type = "III")
  
  if (print_results) {
    cat("Global Metabolic Deviation Analysis:\n")
    print(anova_results)
  }
  
  # Get mean deviation by treatment and time
  summary_stats <- mean_deviation_data %>%
    group_by(tx, time) %>%
    summarise(
      mean = mean(.data[[response_var]], na.rm = TRUE),
      se = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()),
      sd = sd(.data[[response_var]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # Return results as list
  return(list(
    model = model,
    anova = anova_results,
    summary_stats = summary_stats
  ))
}