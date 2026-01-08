apply_main_labels <- function(tbl, labels = main_labels) {
  tbl %>%
    filter(!str_starts(Sample_ID, "nist")) %>%
    filter(!str_starts(Sample_ID, "q3")) %>%
    filter(!str_starts(Sample_ID, "q4")) %>%
    left_join(labels, by = "Sample_ID") %>%
    mutate(tissue_type = if_else(str_starts(pt_tissue_ID, "N"), "Normal", "Tumor")) |>
    relocate(Sample_ID, pt_tissue_ID, tissue_type, tx, time, .before = 1)
}
