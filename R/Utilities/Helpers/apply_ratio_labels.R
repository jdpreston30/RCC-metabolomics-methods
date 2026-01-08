apply_ratio_labels <- function(tbl, labels = ratio_labels) {
  tbl %>%
    filter(!str_starts(Sample_ID, "B")) %>%
    left_join(labels, by = "Sample_ID") %>%
    mutate(tissue_type = if_else(str_starts(name, "N"), "Normal", "Tumor")) |>
    relocate(Sample_ID, name, tissue_type, ratio, .before = 1)
}
