#' Calculate Global Metabolic Deviations from Baseline
#'
#' This function calculates metabolic deviations from GS (baseline) samples
#' for each patient across all features, returning both individual feature
#' deviations and summary statistics per sample.
#'
#' @param data A tibble containing metabolomics data with columns:
#'   - pt_tissue_ID: Patient tissue identifier 
#'   - tx: Treatment (should include "GS" for baseline)
#'   - time: Time point
#'   - Feature columns starting with "HILIC_" or "C18_"
#' @param baseline_tx Character string specifying baseline treatment (default: "GS")
#' @param feature_pattern Regex pattern to identify feature columns (default: "^(HILIC|C18)_")
#'
#' @return A list containing:
#'   - deviations_long: Long format tibble with individual feature deviations
#'   - mean_deviation: Summary tibble with mean deviation per sample
#'   - feature_cols: Vector of feature column names used
#'
#' @export
calculate_global_deviations <- function(data, 
                                        baseline_tx = "GS", 
                                        feature_pattern = "^(HILIC|C18)_") {
  
  # Get GS baseline samples (one per patient at time = 0)
  baseline_samples <- data %>%
    filter(tx == baseline_tx) %>%
    select(pt_tissue_ID, starts_with("HILIC"), starts_with("C18"))
  
  # Get non-baseline samples
  treatment_samples <- data %>%
    filter(tx != baseline_tx)
  
  # Get feature columns
  feature_cols <- colnames(data)[str_detect(colnames(data), feature_pattern)]
  
  # Calculate deviations from baseline for each patient - VECTORIZED APPROACH
  # Convert to matrices for fast computation
  baseline_matrix <- baseline_samples %>%
    select(all_of(feature_cols)) %>%
    as.matrix()
  
  sample_matrix <- treatment_samples %>%
    select(all_of(feature_cols)) %>%
    as.matrix()
  
  # Match patient IDs to get correct baseline for each sample
  patient_match <- match(treatment_samples$pt_tissue_ID, baseline_samples$pt_tissue_ID)
  
  # Subtract baseline from samples (vectorized across all features at once)
  delta_matrix <- sample_matrix - baseline_matrix[patient_match, ]
  
  # Convert back to tibble and add delta column names
  delta_df <- as_tibble(delta_matrix)
  colnames(delta_df) <- paste0(feature_cols, "_delta")
  
  # Combine with metadata
  data_with_deviations <- bind_cols(
    treatment_samples %>% select(Sample_ID, pt_tissue_ID, tx, time),
    delta_df
  )
  
  # Convert to long format for analysis
  delta_cols <- paste0(feature_cols, "_delta")
  
  deviations_long <- data_with_deviations %>%
    select(Sample_ID, pt_tissue_ID, tx, time, all_of(delta_cols)) %>%
    pivot_longer(
      cols = all_of(delta_cols),
      names_to = "feature_delta",
      values_to = "delta_from_baseline"
    ) %>%
    mutate(
      feature = str_remove(feature_delta, "_delta"),
      tx = factor(tx),
      time = factor(time)
    ) %>%
    select(-feature_delta)
  
  # Calculate mean deviation per sample (average across all features)
  mean_deviation <- deviations_long %>%
    group_by(Sample_ID, pt_tissue_ID, tx, time) %>%
    summarise(
      mean_delta = mean(delta_from_baseline, na.rm = TRUE),
      median_delta = median(delta_from_baseline, na.rm = TRUE),
      sd_delta = sd(delta_from_baseline, na.rm = TRUE),
      n_features = sum(!is.na(delta_from_baseline)),
      .groups = "drop"
    )
  
  # Return results as list
  return(list(
    deviations_long = deviations_long,
    mean_deviation = mean_deviation,
    feature_cols = feature_cols
  ))
}