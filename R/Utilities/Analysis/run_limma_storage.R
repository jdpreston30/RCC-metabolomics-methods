#' Run LIMMA for Storage Experiment (Repeated Measures)
#'
#' Analyzes factorial repeated-measures metabolomics with treatment × time interaction.
#' Uses duplicateCorrelation for within-patient blocking and handles missing data with weights.
#'
#' @param data data.frame/tibble with columns: pt_tissue_ID, tx, time, and numeric features
#' @param tx_var character. Treatment column (e.g., "tx")
#' @param time_var character. Time column (e.g., "time")
#' @param block_var character. Patient ID column (e.g., "pt_tissue_ID")
#' @param tx_levels character(2). Reference first, e.g. c("Ice", "RT")
#' @param time_levels character vector. Time points in order, e.g. c("0", "15", "30", "60", "90")
#' @param meta_cols character vector. Metadata columns to exclude from features
#' @param export_prefix optional. If provided, exports CSVs with this prefix
#' @param export_dir character. Output folder. Default "Outputs/LIMMA"
#' @param metaboanalyst logical. Create MetaboAnalyst/Mummichog tables. Default TRUE
#'
#' @return list with:
#'   $tx_main: Treatment main effect (RT vs Ice across all times)
#'   $time_main: Time main effect (change over time, averaged across treatments)
#'   $interaction: Treatment × Time interaction (different time courses?)
#'   $summary_final: Combined stats + means for all tx × time combinations
#'   $mummichog: MetaboAnalyst-ready tables (if metaboanalyst=TRUE)
#'   $consensus_cor: Within-patient correlation estimate
#'
#' @export
run_limma_storage <- function(
    data,
    tx_var = "tx",
    time_var = "time",
    block_var = "pt_tissue_ID",
    tx_levels = c("Ice", "RT"),
    time_levels = c("0", "15", "30", "60", "90"),
    meta_cols = c("Sample_ID", "pt_tissue_ID", "tx", "time"),
    export_prefix = NULL,
    export_dir = "Outputs/LIMMA",
    metaboanalyst = TRUE) {
  
  stopifnot(all(c(tx_var, time_var, block_var) %in% names(data)))
  
  # Prepare factors
  df <- data
  df[[tx_var]] <- factor(df[[tx_var]], levels = tx_levels)
  df[[time_var]] <- factor(df[[time_var]], levels = time_levels)
  df[[block_var]] <- factor(df[[block_var]])
  
  # Remove samples with missing factor values (before creating design matrix)
  complete_rows <- complete.cases(df[, c(tx_var, time_var, block_var)])
  if (!all(complete_rows)) {
    n_dropped <- sum(!complete_rows)
    cat(sprintf("⚠️  Dropping %d samples with missing factor values\n", n_dropped))
    df <- df[complete_rows, ]
  }
  
  # Check for complete tx × time design
  design_check <- df |>
    dplyr::group_by(!!rlang::sym(tx_var), !!rlang::sym(time_var)) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  missing_combos <- design_check |>
    dplyr::filter(n == 0)
  
  if (nrow(missing_combos) > 0) {
    cat("\n⚠️  Missing tx × time combinations:\n")
    print(missing_combos)
  }
  
  # Check for very sparse combinations (n=1 might cause issues)
  sparse_combos <- design_check |>
    dplyr::filter(n < 2)
  
  if (nrow(sparse_combos) > 0) {
    cat("\n⚠️  Sparse combinations (n<2) - may cause non-estimability:\n")
    print(sparse_combos)
    cat("\nConsider dropping these timepoints or treatments for more robust estimates.\n\n")
  }
  
  # Feature columns
  feat_cols <- setdiff(names(df), meta_cols)
  df[feat_cols] <- lapply(df[feat_cols], function(x) suppressWarnings(as.numeric(x)))
  
  # Expression matrix (features × samples) with NA weights
  E <- t(as.matrix(df[, feat_cols, drop = FALSE]))
  W <- ifelse(is.na(E), 0, 1)
  E0 <- E
  E0[is.na(E0)] <- 0
  
  # Set column names to sample IDs
  if ("Sample_ID" %in% names(df)) {
    colnames(E0) <- df$Sample_ID
    colnames(W) <- df$Sample_ID
  }
  
  # Design: ~ tx * time
  fml <- stats::as.formula(paste("~", tx_var, "*", time_var))
  design <- model.matrix(fml, data = df)
  block <- df[[block_var]]
  
  # Diagnostic info
  cat(sprintf("Expression matrix: %d features × %d samples\n", nrow(E0), ncol(E0)))
  cat(sprintf("Design matrix: %d samples × %d coefficients\n", nrow(design), ncol(design)))
  cat(sprintf("Block vector: %d samples\n", length(block)))
  
  # Check dimensions
  if (nrow(design) != ncol(E0)) {
    stop(sprintf(
      "Dimension mismatch: design has %d rows but expression matrix has %d columns",
      nrow(design), ncol(E0)
    ))
  }
  
  # Expected coefficient names
  tx2 <- levels(df[[tx_var]])[2]
  coef_tx <- paste0(tx_var, tx2)  # e.g., "txRT"
  
  # Time coefficients (all except baseline)
  time_coefs <- paste0(time_var, time_levels[-1])  # e.g., "time15", "time30", etc.
  
  # Interaction coefficients
  interact_coefs <- paste0(coef_tx, ":", time_coefs)  # e.g., "txRT:time15"
  
  # Verify coefficients exist
  miss <- setdiff(c(coef_tx, time_coefs, interact_coefs), colnames(design))
  if (length(miss)) {
    stop("Expected coefficients not found in design: ", paste(miss, collapse = ", "))
  }
  
  # LIMMA with duplicateCorrelation + weights
  corfit <- limma::duplicateCorrelation(E0, design, block = block, weights = W)
  fit <- limma::lmFit(E0, design, block = block, correlation = corfit$consensus, weights = W)
  fit <- limma::eBayes(fit)
  
  # Identify estimable coefficients
  estimable <- !is.na(fit$coefficients[1, ])
  estimable_coefs <- colnames(design)[estimable]
  non_estimable <- colnames(design)[!estimable]
  
  if (length(non_estimable) > 0) {
    cat("\n⚠️  Non-estimable coefficients (will be excluded):\n")
    cat(paste("  -", non_estimable, collapse = "\n"), "\n\n")
  }
  
  # Filter coefficient lists to only estimable ones
  time_coefs <- time_coefs[time_coefs %in% estimable_coefs]
  interact_coefs <- interact_coefs[interact_coefs %in% estimable_coefs]
  
  if (length(time_coefs) == 0) {
    stop("No estimable time coefficients - design may be too sparse.")
  }
  
  # Tidy topTable function
  tidy_tt <- function(fit, coef_name, contrast_name) {
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    tibble::tibble(
      feature = rownames(tt),
      t.score = unname(tt[, "t"]),
      p.value = unname(tt[, "P.Value"]),
      FDR = unname(tt[, "adj.P.Val"]),
      logFC = unname(tt[, "logFC"])
    ) |>
      dplyr::mutate(
        contrast = contrast_name,
        Rank = dplyr::row_number()
      ) |>
      dplyr::select(feature, contrast, t.score, p.value, FDR, logFC, Rank)
  }
  
  # Extract main effects and interaction
  tx_tbl <- tidy_tt(fit, coef_tx, paste0(tx2, "_vs_", tx_levels[1]))
  
  # Check if F-statistics are available (they won't be with non-estimable coefficients)
  has_F_stats <- !is.null(fit$F) && !all(is.na(fit$F))
  
  # Time main effect
  if (length(time_coefs) > 0 && has_F_stats) {
    time_tbl <- limma::topTable(fit, coef = time_coefs, number = Inf, sort.by = "none") |>
      tibble::as_tibble() |>
      tibble::add_column(feature = rownames(limma::topTable(fit, coef = time_coefs, 
                                                             number = Inf, sort.by = "none")), 
                         .before = 1) |>
      dplyr::transmute(
        feature,
        contrast = "Time_effect",
        F.score = F,
        p.value = P.Value,
        FDR = adj.P.Val,
        Rank = dplyr::row_number()
      )
  } else if (length(time_coefs) > 0) {
    # Fallback: use first estimable time coefficient as proxy
    cat("⚠️  F-statistics unavailable; using first time coefficient as proxy\n")
    time_tbl <- tidy_tt(fit, time_coefs[1], "Time_effect") |>
      dplyr::mutate(F.score = t.score^2) |>
      dplyr::select(-logFC)
  } else {
    time_tbl <- tibble::tibble(
      feature = rownames(E0),
      contrast = "Time_effect",
      t.score = NA_real_,
      p.value = NA_real_,
      FDR = NA_real_,
      Rank = seq_along(feature)
    )
    cat("⚠️  No time effects testable (insufficient estimable coefficients)\n")
  }
  
  # Interaction
  if (length(interact_coefs) > 0 && has_F_stats) {
    interact_tbl <- limma::topTable(fit, coef = interact_coefs, number = Inf, sort.by = "none") |>
      tibble::as_tibble() |>
      tibble::add_column(feature = rownames(limma::topTable(fit, coef = interact_coefs, 
                                                             number = Inf, sort.by = "none")), 
                         .before = 1) |>
      dplyr::transmute(
        feature,
        contrast = "Interaction",
        F.score = F,
        p.value = P.Value,
        FDR = adj.P.Val,
        Rank = dplyr::row_number()
      )
  } else if (length(interact_coefs) > 0) {
    # Fallback: use first estimable interaction coefficient as proxy
    cat("⚠️  F-statistics unavailable; using first interaction coefficient as proxy\n")
    interact_tbl <- tidy_tt(fit, interact_coefs[1], "Interaction") |>
      dplyr::mutate(F.score = t.score^2) |>
      dplyr::select(-logFC)
  } else {
    interact_tbl <- tibble::tibble(
      feature = rownames(E0),
      contrast = "Interaction",
      t.score = NA_real_,
      p.value = NA_real_,
      FDR = NA_real_,
      Rank = seq_along(feature)
    )
    cat("⚠️  No interaction effects testable (insufficient estimable coefficients)\n")
  }
  
  # Summary table: core stats merged
  core_stats <- tx_tbl |>
    dplyr::select(
      feature,
      p.value_tx = p.value,
      FDR_tx = FDR,
      logFC_tx = logFC
    ) |>
    dplyr::left_join(
      time_tbl |> dplyr::select(feature, p.value_time = p.value, FDR_time = FDR),
      by = "feature"
    ) |>
    dplyr::left_join(
      interact_tbl |> dplyr::select(feature, p.value_interaction = p.value, FDR_interaction = FDR),
      by = "feature"
    )
  
  # Means for each feature by tx × time (wide format)
  means_tbl <- df |>
    dplyr::select(dplyr::all_of(c(tx_var, time_var, feat_cols))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(feat_cols),
      names_to = "feature",
      values_to = "value"
    ) |>
    dplyr::group_by(
      feature,
      !!rlang::sym(tx_var),
      !!rlang::sym(time_var)
    ) |>
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
      mean_col = paste0(
        "Mean_",
        as.character(.data[[tx_var]]), "_",
        as.character(.data[[time_var]])
      )
    ) |>
    tidyr::pivot_wider(names_from = mean_col, values_from = mean_value)
  
  summary_final <- core_stats |>
    dplyr::left_join(means_tbl, by = "feature")
  
  # MetaboAnalyst/Mummichog tables
  mummi <- NULL
  if (isTRUE(metaboanalyst)) {
    to_mummichog <- function(stats_tbl) {
      stats_tbl |>
        dplyr::transmute(
          Feature = feature,
          p.value = p.value,
          mode = dplyr::case_when(
            stringr::str_starts(Feature, "HILIC") ~ "positive",
            stringr::str_starts(Feature, "C18") ~ "negative",
            TRUE ~ NA_character_
          ),
          m.z = suppressWarnings(as.numeric(stringr::str_extract(Feature, "(?<=_)[0-9.]+"))),
          r.t = suppressWarnings(as.numeric(
            stringr::str_extract(Feature, "_[0-9.]+_([0-9.]+)") |>
              stringr::str_extract("[0-9.]+$")
          ))
        ) |>
        dplyr::filter(!is.na(mode) & !is.na(m.z) & !is.na(r.t)) |>
        dplyr::select(m.z, p.value, mode, r.t)
    }
    
    mummi <- list(
      tx_main = to_mummichog(tx_tbl),
      time_main = to_mummichog(time_tbl),
      interaction = to_mummichog(interact_tbl)
    )
  }
  
  # Optional export
  if (!is.null(export_prefix)) {
    dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(tx_tbl, file.path(export_dir, paste0(export_prefix, "_tx_main.csv")))
    readr::write_csv(time_tbl, file.path(export_dir, paste0(export_prefix, "_time_main.csv")))
    readr::write_csv(interact_tbl, file.path(export_dir, paste0(export_prefix, "_interaction.csv")))
    readr::write_csv(summary_final, file.path(export_dir, paste0(export_prefix, "_summary.csv")))
    
    if (isTRUE(metaboanalyst)) {
      readr::write_csv(mummi$tx_main, file.path(export_dir, paste0(export_prefix, "_mummichog_tx.csv")))
      readr::write_csv(mummi$time_main, file.path(export_dir, paste0(export_prefix, "_mummichog_time.csv")))
      readr::write_csv(mummi$interaction, file.path(export_dir, paste0(export_prefix, "_mummichog_interaction.csv")))
    }
  }
  
  list(
    tx_main = tx_tbl,
    time_main = time_tbl,
    interaction = interact_tbl,
    summary_final = summary_final,
    mummichog = mummi,
    consensus_cor = corfit$consensus
  )
}
