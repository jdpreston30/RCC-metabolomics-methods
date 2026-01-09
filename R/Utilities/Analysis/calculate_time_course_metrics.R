#' Calculate Time Course Metabolomics Metrics
#'
#' Calculates mean intensity and feature detection metrics for metabolomics data
#' across different modes (HILIC and C18) or combined, and provides summary 
#' statistics by treatment and time groups.
#'
#' @param data A tibble containing metabolomics data with columns:
#'   - Sample_ID: Sample identifier
#'   - pt_tissue_ID: Patient tissue identifier
#'   - tx: Treatment factor
#'   - time: Time point
#'   - Feature columns starting with "HILIC_" and "C18_"
#' @param exclude_treatments Vector of treatments to exclude (default: "CCM")
#' @param tx_levels Character vector specifying treatment factor levels (default: c("GS", "Ice", "RT"))
#' @param hilic_pattern Regex pattern for HILIC features (default: "^HILIC_")
#' @param c18_pattern Regex pattern for C18 features (default: "^C18_")
#' @param combine_modes Logical, whether to combine HILIC and C18 into single metrics (default: FALSE)
#' @param run_diagnostics Logical, whether to run diagnostic checks (default: TRUE)
#' @param diagnostic_patient Character, patient ID for outlier diagnostics (default: "TM11")
#' @param log_transform_intensity Logical, whether to log2 transform mean intensity values (default: TRUE)
#' @param run_statistics Logical, whether to run statistical analysis (default: TRUE)
#' @param baseline_tx Character, baseline treatment for statistics (default: "GS")
#' @param print_statistics Logical, whether to print statistical results (default: TRUE)
#'
#' @return A list containing:
#'   - metrics_wide: Wide format metrics per sample
#'   - metrics_long: Long format metrics for plotting
#'   - summary_stats: Summary statistics by treatment, time (and mode if not combined)
#'   - diagnostics: List of diagnostic results (if run_diagnostics = TRUE)
#'   - statistics: List of statistical results (if run_statistics = TRUE)
#'
#' @export
calculate_time_course_metrics <- function(data,
                                          exclude_treatments = "CCM",
                                          tx_levels = c("GS", "Ice", "RT"),
                                          hilic_pattern = "^HILIC_",
                                          c18_pattern = "^C18_",
                                          combine_modes = FALSE,
                                          run_diagnostics = TRUE,
                                          diagnostic_patient = "TM11",
                                          log_transform_intensity = TRUE,
                                          run_statistics = TRUE,
                                          baseline_tx = "GS",
                                          print_statistics = TRUE) {
  
  # Filter out excluded treatments
  if (!is.null(exclude_treatments)) {
    time_course_data <- data %>%
      filter(!tx %in% exclude_treatments)
  } else {
    time_course_data <- data
  }
  
  # Get HILIC and C18 feature columns
  hilic_cols <- colnames(time_course_data)[str_detect(colnames(time_course_data), hilic_pattern)]
  c18_cols <- colnames(time_course_data)[str_detect(colnames(time_course_data), c18_pattern)]
  all_feature_cols <- c(hilic_cols, c18_cols)
  
  if (combine_modes) {
    # Calculate combined metrics for each sample across ALL features
    time_course_metrics <- time_course_data %>%
      rowwise() %>%
      mutate(
        # Store raw intensity (for plotting)
        combined_mean_intensity = mean(c_across(all_of(all_feature_cols))[c_across(all_of(all_feature_cols)) > 0], na.rm = TRUE),
        # Store log2 intensity (for statistics)
        combined_mean_intensity_log2 = log2(combined_mean_intensity),
        # Count features (no transformation)
        combined_num_features = sum(c_across(all_of(all_feature_cols)) > 0, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      select(Sample_ID, pt_tissue_ID, tx, time, combined_mean_intensity, combined_mean_intensity_log2, combined_num_features) %>%
      mutate(
        tx = factor(tx, levels = tx_levels),
        time_numeric = ifelse(tx == "GS", 0, as.numeric(time))
      )
    
    # Convert to long format for plotting
    time_course_long <- time_course_metrics %>%
      pivot_longer(
        cols = c(combined_mean_intensity, combined_num_features),
        names_to = "metric",
        names_prefix = "combined_",
        values_to = "value"
      ) %>%
      mutate(
        tx = factor(tx, levels = tx_levels),
        time_numeric = ifelse(tx == "GS", 0, as.numeric(time))
      )
    
    # Calculate summary statistics by treatment and time (no Mode)
    time_course_summary <- time_course_long %>%
      group_by(tx, time_numeric, metric) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        se = sd(value, na.rm = TRUE) / sqrt(n()),
        n_samples = n(),
        .groups = "drop"
      )
    
  } else {
    # Calculate metrics for each sample by mode
    time_course_metrics <- time_course_data %>%
      rowwise() %>%
      mutate(
        # HILIC metrics
        hilic_mean_intensity_raw = mean(c_across(all_of(hilic_cols))[c_across(all_of(hilic_cols)) > 0], na.rm = TRUE),
        hilic_mean_intensity = if(log_transform_intensity) log2(hilic_mean_intensity_raw) else hilic_mean_intensity_raw,
        hilic_num_features = sum(c_across(all_of(hilic_cols)) > 0, na.rm = TRUE),
        # C18 metrics
        c18_mean_intensity_raw = mean(c_across(all_of(c18_cols))[c_across(all_of(c18_cols)) > 0], na.rm = TRUE),
        c18_mean_intensity = if(log_transform_intensity) log2(c18_mean_intensity_raw) else c18_mean_intensity_raw,
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
        tx = factor(tx, levels = tx_levels),
        time_numeric = ifelse(tx == "GS", 0, as.numeric(time))
      )
    
    # Calculate summary statistics by treatment, time, and mode
    time_course_summary <- time_course_long %>%
      group_by(tx, time_numeric, Mode, metric) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        se = sd(value, na.rm = TRUE) / sqrt(n()),
        n_samples = n(),
        .groups = "drop"
      )
  }
  
  # Initialize diagnostics list
  diagnostics <- list()
  
  # Run diagnostics if requested
  if (run_diagnostics) {
    # Check sample sizes per group
    if (combine_modes) {
      sample_sizes <- time_course_summary %>%
        filter(metric == "num_features") %>%
        select(tx, time_numeric, n_samples)
    } else {
      sample_sizes <- time_course_summary %>%
        filter(metric == "num_features", Mode == "HILIC") %>%
        select(tx, time_numeric, n_samples)
    }
    
    # Check if specified patient is an outlier in feature detection
    if (combine_modes) {
      patient_check <- time_course_metrics %>%
        filter(tx %in% c("Ice", "RT")) %>%
        group_by(tx, time) %>%
        summarise(
          mean_all = mean(combined_num_features, na.rm = TRUE),
          !!paste0("mean_", tolower(diagnostic_patient)) := mean(combined_num_features[pt_tissue_ID == diagnostic_patient], na.rm = TRUE),
          mean_others = mean(combined_num_features[pt_tissue_ID != diagnostic_patient], na.rm = TRUE),
          !!paste0(tolower(diagnostic_patient), "_diff") := get(paste0("mean_", tolower(diagnostic_patient))) - mean_others,
          .groups = "drop"
        )
    } else {
      patient_check <- time_course_metrics %>%
        filter(tx %in% c("Ice", "RT")) %>%
        group_by(tx, time) %>%
        summarise(
          mean_hilic_all = mean(hilic_num_features, na.rm = TRUE),
          !!paste0("mean_hilic_", tolower(diagnostic_patient)) := mean(hilic_num_features[pt_tissue_ID == diagnostic_patient], na.rm = TRUE),
          mean_hilic_others = mean(hilic_num_features[pt_tissue_ID != diagnostic_patient], na.rm = TRUE),
          !!paste0(tolower(diagnostic_patient), "_diff") := get(paste0("mean_hilic_", tolower(diagnostic_patient))) - mean_hilic_others,
          .groups = "drop"
        )
    }
    
    diagnostics <- list(
      sample_sizes = sample_sizes,
      patient_outlier_check = patient_check
    )
    
    # Print diagnostics
    cat("\n=== Sample Sizes by Treatment × Time ===\n")
    print(sample_sizes, n = Inf)
    
    cat(paste0("\n=== ", diagnostic_patient, " Feature Detection Compared to Other Patients ===\n"))
    print(patient_check)
  }
  
  # Run statistical analysis if requested
  statistics <- NULL
  if (run_statistics && combine_modes) {
    # Load required packages for statistical analysis
    require(lmerTest)
    require(emmeans)
    
    # Determine which metrics to analyze (use log2 for intensity)
    metrics_to_analyze <- list(
      combined_num_features = "combined_num_features",
      combined_mean_intensity = "combined_mean_intensity_log2"
    )
    
    # Initialize statistics list
    statistics <- list()
    
    for (metric_name in names(metrics_to_analyze)) {
      metric_col <- metrics_to_analyze[[metric_name]]
      
      # Calculate deviations from baseline for this metric
      baseline_data <- time_course_metrics %>%
        filter(tx == baseline_tx) %>%
        select(pt_tissue_ID, baseline_value = !!sym(metric_col))
      
      deviation_data <- time_course_metrics %>%
        filter(tx != baseline_tx) %>%
        left_join(baseline_data, by = "pt_tissue_ID") %>%
        mutate(deviation = !!sym(metric_col) - baseline_value) %>%
        select(pt_tissue_ID, tx, time, deviation) %>%
        mutate(
          tx = factor(tx, levels = setdiff(tx_levels, baseline_tx)),
          time = factor(time)
        )
      
      # Fit mixed-effects model with lmerTest
      model <- lmerTest::lmer(deviation ~ tx * time + (1|pt_tissue_ID), data = deviation_data)
      
      # Type III ANOVA with Satterthwaite
      anova_result <- anova(model)
      
      # Estimated marginal means
      emm_tx <- emmeans(model, ~ tx)
      emm_time <- emmeans(model, ~ time)
      emm_tx_time <- emmeans(model, ~ tx * time)
      
      # Pairwise comparisons
      pairwise_tx <- pairs(emm_tx)
      pairwise_time_by_tx <- contrast(emm_tx_time, "pairwise", by = "time")
      
      # Store results (use original metric name, not _log2 suffix)
      statistics[[metric_name]] <- list(
        model = model,
        anova = anova_result,
        emmeans_tx = emm_tx,
        emmeans_time = emm_time,
        emmeans_tx_time = emm_tx_time,
        pairwise_tx = pairwise_tx,
        pairwise_time_by_tx = pairwise_time_by_tx
      )
      
      # Print results if requested
      if (print_statistics) {
        cat("\n", strrep("=", 80), "\n", sep = "")
        cat("Statistical Analysis: ", metric_name, " (log2 transformed for intensity)\n", sep = "")
        cat(strrep("=", 80), "\n\n", sep = "")
        
        cat("ANOVA Results (Type III):\n")
        print(anova_result)
        
        cat("\n\nEstimated Marginal Means by Treatment:\n")
        print(summary(emm_tx))
        
        cat("\n\nPairwise Comparisons (Treatment):\n")
        print(summary(pairwise_tx))
        
        cat("\n\nPairwise Comparisons (Treatment at each Time):\n")
        print(summary(pairwise_time_by_tx))
      }
    }
  }
  
  # Return results as list
  return(list(
    metrics_wide = time_course_metrics,
    metrics_long = time_course_long,
    summary_stats = time_course_summary,
    diagnostics = if(run_diagnostics) diagnostics else NULL,
    statistics = statistics,
    feature_columns = list(hilic = hilic_cols, c18 = c18_cols)
  ))
}