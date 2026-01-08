#' Filter Feature Table by Unique Value Percentage (Standalone)
#'
#' This function calculates the percentage of unique values for each numeric column
#' in a feature table and retains only columns that meet or exceed a specified
#' uniqueness threshold. Non-numeric columns are always retained.
#' This is a standalone version that can be used in pipes independently of preprocess_FT().
#'
#' @param feature_table A data frame or tibble containing feature data
#' @param unique_threshold Numeric value between 0 and 1 (default 0.7). 
#'   Columns with unique percentage >= this threshold are retained.
#'   For example, 0.7 means 70% or more unique values are required to keep the column.
#' @param verbose Logical (default FALSE). If TRUE, prints information about 
#'   which columns were removed and their unique percentages.
#'
#' @return A filtered data frame/tibble with only columns meeting the uniqueness threshold
#'   plus all non-numeric columns.
#'
#' @details
#' The function:
#' 1. Identifies numeric columns in the feature table
#' 2. Calculates unique percentage as: (number of unique non-NA values) / (total non-NA values)
#' 3. Retains columns where unique percentage >= unique_threshold
#' 4. Always retains non-numeric columns (e.g., Patient, grouping variables)
#'
#' @examples
#' \dontrun{
#' # Keep columns with 70% or more unique values (default)
#' filtered_data <- feature_table %>% filter_unique_features()
#' 
#' # Keep columns with 80% or more unique values
#' filtered_data <- feature_table %>% filter_unique_features(unique_threshold = 0.8)
#' 
#' # With verbose output to see what was removed
#' filtered_data <- feature_table %>% filter_unique_features(verbose = TRUE)
#' }
#'
#' @export
filter_unique_features <- function(feature_table, unique_threshold = 0.8, verbose = FALSE) {
  
  # Validate inputs
  if (!is.data.frame(feature_table)) {
    stop("feature_table must be a data frame or tibble")
  }
  
  if (!is.numeric(unique_threshold) || unique_threshold < 0 || unique_threshold > 1) {
    stop("unique_threshold must be a numeric value between 0 and 1")
  }
  
  # Identify numeric columns
  numeric_cols <- sapply(feature_table, is.numeric)
  numeric_col_names <- names(feature_table)[numeric_cols]
  non_numeric_col_names <- names(feature_table)[!numeric_cols]
  
  if (verbose) {
    cat("Analyzing", length(numeric_col_names), "numeric columns for uniqueness...\n")
  }
  
  # Calculate unique percentages for numeric columns
  unique_percentages <- sapply(numeric_col_names, function(col_name) {
    col_data <- feature_table[[col_name]]
    # Remove NA values for calculation
    non_na_data <- col_data[!is.na(col_data)]
    
    if (length(non_na_data) == 0) {
      return(0)  # If all NA, uniqueness is 0
    }
    
    unique_count <- length(unique(non_na_data))
    total_count <- length(non_na_data)
    
    return(unique_count / total_count)
  })
  
  # Identify columns that meet the threshold
  columns_to_keep <- numeric_col_names[unique_percentages >= unique_threshold]
  columns_to_remove <- numeric_col_names[unique_percentages < unique_threshold]
  
  # Always keep non-numeric columns
  final_columns_to_keep <- c(non_numeric_col_names, columns_to_keep)
  
  if (verbose) {
    cat("Threshold:", unique_threshold * 100, "%\n")
    cat("Columns meeting threshold:", length(columns_to_keep), "/", length(numeric_col_names), "\n")
    
    if (length(columns_to_remove) > 0) {
      cat("Removed columns (unique %):\n")
      removed_percentages <- unique_percentages[columns_to_remove]
      for (i in seq_along(columns_to_remove)) {
        cat("  ", columns_to_remove[i], ": ", round(removed_percentages[i] * 100, 1), "%\n", sep = "")
      }
    }
  }
  
  # Filter the feature table
  filtered_table <- feature_table[, final_columns_to_keep, drop = FALSE]
  
  # Print summary report
  n_above_threshold <- length(columns_to_keep)
  n_below_threshold <- length(columns_to_remove)
  total_numeric_features <- length(numeric_col_names)
  
  pct_above <- round((n_above_threshold / total_numeric_features) * 100, 1)
  pct_below <- round((n_below_threshold / total_numeric_features) * 100, 1)
  
  cat("\n=== Feature Filtering Report ===\n")
  cat("Unique threshold:", unique_threshold * 100, "%\n")
  cat("Features ≥ threshold:", n_above_threshold, "(", pct_above, "% of total features)\n")
  cat("Features < threshold:", n_below_threshold, "(", pct_below, "% of total features)\n")
  cat("Total numeric features analyzed:", total_numeric_features, "\n")
  cat("Non-numeric columns retained:", length(non_numeric_col_names), "\n")
  cat("Final dataset dimensions:", nrow(filtered_table), "rows ×", ncol(filtered_table), "columns\n")
  cat("==============================\n\n")
  
  return(filtered_table)
}
