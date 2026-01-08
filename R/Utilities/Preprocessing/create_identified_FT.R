#' Create Identified Feature Table
#'
#' This function combines feature name parsing, library matching, and feature table
#' subsetting into one comprehensive workflow. It parses feature names from an
#' untargeted feature table, matches them against a reference library, and creates
#' a subset table containing only library-confirmed features.
#'
#' @param feature_table Data frame with feature columns starting with C18 or HILIC
#' @param reference_library Tibble with columns: CNAME, Column, mz_lib, time_lib, Adduct, KEGGID, Library_Isomer
#' @param mz_thresh_ppm M/z threshold in ppm (default: 5)
#' @param time_thresh_sec Retention time threshold in seconds (default: 30)
#' @return List with two elements:
#'   - TFT_confirmed: Subset of original table with only matched features
#'   - matched_features: Annotation metadata with unique_match column
#'
#' @examples
#' \dontrun{
#'   results <- create_identified_FT(
#'     feature_table = UFT,
#'     reference_library = idx_lib_simple,
#'     mz_thresh_ppm = 5,
#'     time_thresh_sec = 30
#'   )
#'   
#'   TFT_confirmed <- results$TFT_confirmed
#'   annotations <- results$matched_features
#' }
#'
#' @export
create_identified_FT <- function(feature_table, 
                                reference_library, 
                                mz_thresh_ppm = 5, 
                                time_thresh_sec = 30) {
  
  # Load required libraries
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  
  cat("=== Creating Identified Feature Table ===\n")
  cat("Parameters:\n")
  cat("  M/Z threshold:", mz_thresh_ppm, "ppm\n")
  cat("  RT threshold:", time_thresh_sec, "seconds\n\n")
  
  # Step 1: Parse feature names
  cat("Step 1: Parsing feature names...\n")
  
  # Get column names that start with C18 or HILIC
  feature_columns <- colnames(feature_table)[str_detect(colnames(feature_table), "^(C18|HILIC)")]
  
  # Create tibble with parsed information
  feature_info <- tibble(
    feature = feature_columns
  ) %>%
    mutate(
      # Split the feature name by underscores
      parts = str_split(feature, "_"),
      # Extract mode (first part)
      mode = map_chr(parts, ~.x[1]),
      # Extract m/z (second part, convert to numeric)
      mz = map_dbl(parts, ~as.numeric(.x[2])),
      # Extract RT (third part, convert to numeric)
      rt = map_dbl(parts, ~as.numeric(.x[3]))
    ) %>%
    # Remove the temporary parts column
    select(feature, mode, mz, rt)
  
  cat("  Parsed", nrow(feature_info), "features\n")
  cat("  Modes found:", paste(unique(feature_info$mode), collapse = ", "), "\n\n")
  
  # Step 2: Match features against reference library
  cat("Step 2: Matching features against library...\n")
  
  # Initialize results list
  all_matches <- list()
  
  # Iterate through each feature
  for (i in 1:nrow(feature_info)) {
    
    current_feature <- feature_info[i, ]
    
    # Filter library by matching column/mode
    # Map feature modes to library column names
    mode_mapping <- c("C18" = "C18", "HILIC" = "HILIC")
    
    if (!current_feature$mode %in% names(mode_mapping)) {
      next  # Skip if mode not recognized
    }
    
    library_subset <- reference_library %>%
      filter(Column == mode_mapping[current_feature$mode])
    
    if (nrow(library_subset) == 0) {
      next  # Skip if no library entries for this mode
    }
    
    # Calculate m/z threshold in Daltons (ppm-based)
    mz_threshold_da <- (mz_thresh_ppm * current_feature$mz) / 1000000
    
    # Find m/z matches
    mz_matches <- which(abs(library_subset$mz_lib - current_feature$mz) <= mz_threshold_da)
    
    if (length(mz_matches) == 0) {
      next  # Skip if no m/z matches
    }
    
    # Check retention time matches for each m/z match
    for (match_idx in mz_matches) {
      
      library_entry <- library_subset[match_idx, ]
      
      # Calculate retention time difference
      rt_diff <- abs(current_feature$rt - library_entry$time_lib)
      
      # Check if within retention time threshold
      if (rt_diff <= time_thresh_sec) {
        
        # Create match result
        match_result <- tibble(
          # Feature information
          feature = current_feature$feature,
          feature_mode = current_feature$mode,
          feature_mz = current_feature$mz,
          feature_rt = current_feature$rt,
          
          # Library information
          compound_name = library_entry$CNAME,
          library_mode = library_entry$Column,
          library_mz = library_entry$mz_lib,
          library_rt = library_entry$time_lib,
          adduct = library_entry$Adduct,
          kegg_id = library_entry$KEGGID,
          library_isomer = library_entry$Library_Isomer,
          
          # Match quality metrics
          mz_error_ppm = ((current_feature$mz - library_entry$mz_lib) / library_entry$mz_lib) * 1000000,
          rt_error_sec = current_feature$rt - library_entry$time_lib,
          rt_error_abs_sec = rt_diff
        )
        
        all_matches <- append(all_matches, list(match_result))
      }
    }
  }
  
  # Combine all matches
  if (length(all_matches) > 0) {
    matched_features <- bind_rows(all_matches)
    
    # Sort by feature name and match quality
    matched_features <- matched_features %>%
      arrange(feature, rt_error_abs_sec, abs(mz_error_ppm))
    
    cat("  Total matches found:", nrow(matched_features), "\n")
    cat("  Features with matches:", length(unique(matched_features$feature)), "\n")
    cat("  Unique compounds matched:", length(unique(matched_features$compound_name)), "\n\n")
    
  } else {
    cat("  No matches found with current thresholds.\n")
    return(list(TFT_confirmed = tibble(), matched_features = tibble()))
  }
  
  # Step 3: Add unique match indicator
  cat("Step 3: Adding unique match indicators...\n")
  
  # Count matches per feature
  match_counts <- matched_features %>%
    group_by(feature) %>%
    summarise(match_count = n(), .groups = "drop")
  
  # Add unique_match column
  matched_features <- matched_features %>%
    left_join(match_counts, by = "feature") %>%
    mutate(unique_match = match_count == 1) %>%
    select(-match_count)
  
  unique_matches <- sum(matched_features$unique_match)
  multiple_matches <- nrow(matched_features) - unique_matches
  
  cat("  Unique matches:", unique_matches, "\n")
  cat("  Multiple matches:", multiple_matches, "\n\n")
  
  # Step 4: Create TFT_confirmed
  cat("Step 4: Creating TFT_confirmed...\n")
  
  # Get unique feature names that had matches
  matched_feature_names <- unique(matched_features$feature)
  
  # Get all column names from feature table
  ft_columns <- colnames(feature_table)
  
  # Find which columns are the matched features
  feature_columns_to_keep <- ft_columns[ft_columns %in% matched_feature_names]
  
  # Get non-feature columns (like sample info, metadata, etc.)
  non_feature_columns <- ft_columns[!grepl("^(C18|HILIC)", ft_columns)]
  
  # Create TFT_confirmed with metadata columns + matched feature columns
  columns_to_select <- c(non_feature_columns, feature_columns_to_keep)
  
  TFT_confirmed <- feature_table %>%
    select(all_of(columns_to_select))
  
  # Print summary
  cat("  Original features:", sum(grepl("^(C18|HILIC)", colnames(feature_table))), "\n")
  cat("  Matched features kept:", length(feature_columns_to_keep), "\n")
  cat("  Non-feature columns:", length(non_feature_columns), "\n")
  cat("  Total columns in TFT_confirmed:", ncol(TFT_confirmed), "\n")
  cat("  Rows in TFT_confirmed:", nrow(TFT_confirmed), "\n\n")
  
  cat("=== Process Complete ===\n")
  
  # Return both results
  return(list(
    TFT_confirmed = TFT_confirmed,
    matched_features = matched_features
  ))
}