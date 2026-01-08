#* 0c: Feature Table Import and Preprocessing
#+ 0c.0: Import Metadata
#- 0c.0.1: Labels and weights (ratio)
ratio_labels <- read_excel(resolve_file(config, "labels_weights", "ratio"), sheet = "RCC_Labels") |>
  select(Sample_ID, name, ratio)
#- 0c.0.2: Feature key (ratio)
ratio_msmica_key <- read_csv(resolve_file(config, "MSMICA_key", "ratio"))
#- 0c.0.3: Labels and weights (main)
main_labels <- read_excel(resolve_file(config, "labels_weights", "main"), sheet = "RCC_Labels") |>
  mutate(
    time = str_remove(time, "'"),
    Sample_ID = str_remove(Sample_ID, "'")
  )
#- 0c.0.2: Feature key (main)
main_msmica_key <- read_csv(resolve_file(config, "MSMICA_key", "main"))
#+ 0c.2: Import feature tables for ratio experiment; join metadata
#- 0c.2.1: MSMICA feature tables (ratio)
ratio_msmica_20p <- read_csv(resolve_file(config, "MSMICA_20p", "ratio")) |>
  apply_ratio_labels()
#- 0c.2.2: Untargeted feature tables (ratio)
ratio_uft_full <- read_csv(resolve_file(config, "UFT_full", "ratio")) |>
  apply_ratio_labels()
ratio_uft_full_0s <- read_csv(resolve_file(config, "UFT_full_0s", "ratio")) |>
  apply_ratio_labels()
ratio_uft_20p <- read_csv(resolve_file(config, "UFT_20p", "ratio")) |>
  apply_ratio_labels()
#+ 0c.3: Import feature tables for main experiment; join metadata
#! Note: TM11_RT_60 missing for all, did not collect (experimental error)
#- 0c.3.1: MSMICA feature tables (main)
main_msmica_20p <- read_csv(resolve_file(config, "MSMICA_20p", "main")) |>
  apply_main_labels()
#- 0c.3.2: Untargeted feature tables (main)
main_uft_full <- read_csv(resolve_file(config, "UFT_full", "main")) |>
  apply_main_labels()
main_uft_full_0s <- read_csv(resolve_file(config, "UFT_full_0s", "main")) |>
  apply_main_labels()
main_uft_20p <- read_csv(resolve_file(config, "UFT_20p", "main")) |>
  apply_main_labels()