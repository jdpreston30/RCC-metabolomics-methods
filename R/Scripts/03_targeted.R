#* 3: Storage Experiment LIMMA Analysis
#+ 3.1: Subset MSMICA data to just tumors; exclude CCM and GS samples
main_msmica_20p_tumors <- main_msmica_20p |>
  filter(tissue_type == "Tumor") |>
  filter(!str_detect(Sample_ID, "CCM|GS")) |>  # Exclude CCM and gold standard
  select(-tissue_type) |>
  mutate(
    time = factor(time, levels = c("15", "30", "60", "90")),  # No time 0 since GS excluded
    tx = factor(tx, levels = c("Ice", "RT"))  # Ice as reference
  )
#+ 3.2: Run LIMMA analysis
u()
limma_storage <- run_limma_storage(
  data = main_msmica_20p_tumors,
  tx_var = "tx",
  time_var = "time",
  block_var = "pt_tissue_ID",
  tx_levels = c("Ice", "RT"),
  time_levels = c("15", "30", "60", "90"),  # Updated - no time 0
  meta_cols = c("Sample_ID", "pt_tissue_ID", "tx", "time"),
  export_prefix = "storage_experiment",
  export_dir = "Outputs/LIMMA",
  metaboanalyst = TRUE
)
#+ 3.3: Summarize results
#- 3.3.1: Extract summary statistics
summary_stoarage <- limma_storage$summary_final |>
  select(feature, p.value_tx, p.value_time, p.value_interaction) |>
  unique()
#- 3.3.2: Join with superclass key
superclass_key <- main_msmica_key |>
  select(feature, Superclass) |>
  unique() |>
  left_join(summary_stoarage, by = "feature") |>
  filter(!is.na(Superclass))
#- 3.3.3: Create stacked bar plots by superclass and significance
# Treatment effect plot
tx_superclass_summary <- superclass_key |>
  mutate(significance = ifelse(p.value_tx < 0.05, "Significant", "Not Significant")) |>
  count(Superclass, significance) |>
  group_by(Superclass) |>
  mutate(total = sum(n)) |>
  ungroup() |>
  arrange(desc(total))
plot_tx_superclass <- ggplot(tx_superclass_summary, 
                              aes(y = reorder(Superclass, total), x = n, fill = significance)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Significant" = "#D32F2F", "Not Significant" = "gray70")) +
  labs(
    title = "Ice vs. RT Effect",
    y = NULL,
    x = "Number of Features",
    fill = "P-value (< 0.05)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )
# Time effect plot
time_superclass_summary <- superclass_key |>
  mutate(significance = ifelse(p.value_time < 0.05, "Significant", "Not Significant")) |>
  count(Superclass, significance) |>
  group_by(Superclass) |>
  mutate(total = sum(n)) |>
  ungroup() |>
  arrange(desc(total))
plot_time_superclass <- ggplot(time_superclass_summary, 
                                aes(y = reorder(Superclass, total), x = n, fill = significance)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Significant" = "#D32F2F", "Not Significant" = "gray70")) +
  labs(
    title = "Time Effect",
    y = NULL,
    x = "Number of Features",
    fill = "P-value (< 0.05)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )
# Interaction effect plot
interaction_superclass_summary <- superclass_key |>
  mutate(significance = ifelse(p.value_interaction < 0.05, "Significant", "Not Significant")) |>
  count(Superclass, significance) |>
  group_by(Superclass) |>
  mutate(total = sum(n)) |>
  ungroup() |>
  arrange(desc(total))
plot_interaction_superclass <- ggplot(interaction_superclass_summary, 
                                       aes(y = reorder(Superclass, total), x = n, fill = significance)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Significant" = "#D32F2F", "Not Significant" = "gray70")) +
  labs(
    title = "Interaction",
    y = NULL,
    x = "Number of Features",
    fill = "P-value (< 0.05)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )
# Display plots
print_to_png(plot_tx_superclass, "Outputs/Figures/storage_tx_superclass.png", 6, 4)
print_to_png(plot_time_superclass, "Outputs/Figures/storage_time_superclass.png", 6, 4)
print_to_png(plot_interaction_superclass, "Outputs/Figures/storage_interaction_superclass.png", 7, 4)
