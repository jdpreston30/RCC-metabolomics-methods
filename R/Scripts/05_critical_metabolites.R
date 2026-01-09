#* 5: Critical Metabolites Deviation Analysis
#+ 5.0: Load cancer metabolites list
cancer_metabolites_metadata <- cancer_metabolites |>
  select(feature, identified_name) |>
  unique()
#+ 5.1: Prepare data - remove CCM and calculate deviations from GS
#- 5.1.1: Filter out CCM samples
TFT_cancer_no_CCM <- TFT_combined_cancer_metabolites |>
  filter(tx != "CCM")
#- 5.1.2: Separate GS baseline values
GS_baseline <- TFT_cancer_no_CCM |>
  filter(tx == "GS") |>
  select(pt_tissue_ID, starts_with("HILIC"), starts_with("C18"))
#- 5.1.3: Get Ice and RT samples only
TFT_cancer_with_deviations <- TFT_cancer_no_CCM |>
  filter(tx != "GS") |>
  # Join with GS baseline
  left_join(
    GS_baseline |>
      rename_with(~paste0(.x, "_GS"), .cols = starts_with(c("HILIC", "C18"))),
    by = "pt_tissue_ID"
  )
#- 5.1.4: Calculate deviations (delta from GS for each feature)
feature_cols <- colnames(TFT_cancer_no_CCM)[str_detect(colnames(TFT_cancer_no_CCM), "^(HILIC|C18)_")]
for (feat in feature_cols) {
  gs_col <- paste0(feat, "_GS")
  delta_col <- paste0(feat, "_delta")
  
  TFT_cancer_with_deviations[[delta_col]] <- TFT_cancer_with_deviations[[feat]] - TFT_cancer_with_deviations[[gs_col]]
}
#+ 5.2: Statistical analysis to identify candidate features
#- 5.2.1: Prepare long format for analysis
delta_cols <- colnames(TFT_cancer_with_deviations)[str_detect(colnames(TFT_cancer_with_deviations), "_delta$")]
cancer_deviations_long <- TFT_cancer_with_deviations |>
  select(Sample_ID, pt_tissue_ID, tx, time, all_of(delta_cols)) |>
  pivot_longer(
    cols = all_of(delta_cols),
    names_to = "feature",
    values_to = "delta_from_GS"
  ) |>
  mutate(
    feature = str_remove(feature, "_delta$"),
    time = factor(time, levels = c("15", "30", "60", "90")),
    tx = factor(tx, levels = c("Ice", "RT"))
  )
#- 5.2.2: Run 2-way repeated measures analysis for each feature
feature_stats <- cancer_deviations_long |>
  group_by(feature) |>
  summarise(
    test_result = tryCatch({
      # Fit linear mixed model: treatment × time with patient random effect
      model <- lmer(delta_from_GS ~ tx * time + (1|pt_tissue_ID), data = pick(everything()))
      anova_result <- anova(model, type = "III")
      
      # Extract p-values for main effects and interaction
      list(tibble(
        tx_F = anova_result$`F value`[2],
        tx_p = anova_result$`Pr(>F)`[2],
        time_F = anova_result$`F value`[3],
        time_p = anova_result$`Pr(>F)`[3],
        interaction_F = anova_result$`F value`[4],
        interaction_p = anova_result$`Pr(>F)`[4],
        n_obs = nrow(pick(everything())),
        mean_Ice = mean(delta_from_GS[tx == "Ice"], na.rm = TRUE),
        mean_RT = mean(delta_from_GS[tx == "RT"], na.rm = TRUE),
        mean_abs_delta = mean(abs(delta_from_GS), na.rm = TRUE)
      ))
    }, error = function(e) {
      list(tibble(
        tx_F = NA, tx_p = NA, time_F = NA, time_p = NA, 
        interaction_F = NA, interaction_p = NA, n_obs = 0,
        mean_Ice = NA, mean_RT = NA, mean_abs_delta = NA
      ))
    }),
    .groups = "drop"
  ) |>
  unnest(test_result) |>
  mutate(
    tx_FDR = p.adjust(tx_p, method = "BH"),
    time_FDR = p.adjust(time_p, method = "BH"),
    interaction_FDR = p.adjust(interaction_p, method = "BH")
  )
#- 5.2.3: Identify candidate features for visualization
candidates_tx <- feature_stats |>
  filter(tx_p < 0.05) |>
  arrange(tx_p)
candidates_time <- feature_stats |>
  filter(time_p < 0.05) |>
  arrange(time_p)
candidates_interaction <- feature_stats |>
  filter(interaction_p < 0.05) |>
  arrange(interaction_p)
candidates_any <- feature_stats |>
  filter(tx_p < 0.05 | time_p < 0.05 | interaction_p < 0.05) |>
  arrange(interaction_p, tx_p, time_p)
# Store for later use
feature_candidates <- candidates_any |>
  left_join(cancer_metabolites_metadata, by = "feature") |>
  select(identified_name)
#+ 5.3: Visualize deviation patterns
#- 5.3.1: Heatmap of mean deviations by treatment (Ice vs RT only)
deviation_summary <- cancer_deviations_long |>
  group_by(feature, tx, time) |>
  summarise(mean_delta = mean(delta_from_GS, na.rm = TRUE), .groups = "drop") |>
  left_join(cancer_metabolites_metadata, by = "feature") |>
  mutate(
    identified_name = map_chr(identified_name, metabolite_sentence_case)
  )  |>
  select(identified_name, everything())
#- 5.3.2: Generate trajectory plots for all metabolites
{
  # Create output directory
  output_dir <- "Outputs/target_test"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Get unique combinations of identified_name and feature
  # This ensures we create separate plots when one metabolite has multiple features (MMD)
  unique_combinations <- deviation_summary %>%
    distinct(identified_name, feature)
  # Initialize storage list for plots
  trajectory_plots <- list()
  # Loop through each unique metabolite-feature combination
  cat("\nGenerating trajectory plots for", nrow(unique_combinations), "metabolite-feature combinations...\n")
  for (i in seq_len(nrow(unique_combinations))) {
    metabolite <- unique_combinations$identified_name[i]
    feature_id <- unique_combinations$feature[i]
    # Count how many features this metabolite has already processed
    existing_count <- sum(str_detect(names(trajectory_plots), paste0("^", metabolite)))
    # Create unique name with suffix if multiple features per metabolite
    if (existing_count == 0) {
      plot_name <- metabolite
    } else {
      plot_name <- paste0(metabolite, "_", existing_count + 1)
    }
    # Generate plot
    tryCatch({
      plot_obj <- plot_trajectory(deviation_summary, identified_name = metabolite, feature = feature_id)
      # Store plot object
      trajectory_plots[[plot_name]] <- plot_obj
      # Save to PNG - pass just filename and specify output_dir, disable auto_open
      png_filename <- paste0(plot_name, ".png")
      print_to_png(plot_obj, png_filename, width = 5, height = 4, output_dir = output_dir, auto_open = FALSE)
      cat("  ✓", plot_name, "\n")
    }, error = function(e) {
      cat("  ✗", metabolite, "- Error:", e$message, "\n")
    })
  }
}
