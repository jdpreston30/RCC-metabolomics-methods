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
main_msmica_key <- read_csv(resolve_file(config, "MSMICA_key", "main")) |>
  rename(
    feature = Feature,
    identified_name = `Identified Name`,
    rt = `Retention Time`,
    exact_mass = `Exact Mass`,
    ion_mode = `Ion Mode`,
    MMD = `Multi-Mode Detection`,
    mean_intensity = `Mean Intensity`,
    mean_intensity_w_zeros = `Mean Intensity w/ Zeros`,
    identification_method = `Identification Method`,
    annotation_probability = `Annotation Probability`,
    KEGG = `KEGG ID`,
    CID = `PubChem CID`,
    inchi_key = InChIKey,
    HMDB = `HMDB ID`,
    alternative_parent = `Alternative Parent`,
    isomer = Isomer
  )
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
  apply_main_labels() |>
  filter_unique_features()
#- 0c.3.2: Untargeted feature tables (main)
main_uft_full <- read_csv(resolve_file(config, "UFT_full", "main")) |>
  apply_main_labels()
main_uft_full_0s <- read_csv(resolve_file(config, "UFT_full_0s", "main")) |>
  apply_main_labels()
main_uft_20p <- read_csv(resolve_file(config, "UFT_20p", "main")) |>
  apply_main_labels() |>
  filter_unique_features()
#+ 0c.4: Make merged TFT of MSMICA and library confirmed features (main experiment)
#- 0c.4.1: Bring in IROA IDX feature library
hfqe_lib <- read_csv(config$paths$HFQE_lib)
#- 0c.4.2: Run function to match features to library and create identified TFT_confirmed
identified <- create_identified_FT(
  feature_table = main_uft_20p,
  reference_library = hfqe_lib,
  mz_thresh_ppm = 5,
  time_thresh_sec = 30
)
#- 0c.4.3: Assign and apply unique filter fxn just in case
TFT_confirmed <- identified$TFT_confirmed |>
  filter_unique_features()
#- 0c.4.4: Build key
TFT_confirmed_key <- identified$matched_features %>%
  select(identified_name = compound_name, isomer = library_isomer, everything()) %>%
  mutate(MMD = "")
#+ 0c.5: Create merged library/annotated TFT
#- 0c.5.1: Subset features from annotations with source tracking
TFT_annot_features <- main_msmica_key %>%
  mutate(lib_conf = "N", source = "annotation") %>%
  select(feature, lib_conf, source) %>%
  unique()
#- 0c.5.2: Subset features from confirmed with source tracking
TFT_confirmed_features <- TFT_confirmed_key %>%
  mutate(lib_conf = "Y", source = "library") %>%
  select(feature, lib_conf, source) %>%
  unique()
#- 0c.5.3: Identify features present in both datasets
overlapping_features_lookup <- intersect(TFT_annot_features$feature, TFT_confirmed_features$feature)
#- 0c.5.4: Create source labels for overlapping and unique features
TFT_merged_features <- bind_rows(TFT_confirmed_features, TFT_annot_features) %>%
  arrange(feature, desc(lib_conf)) %>%
  group_by(feature) %>%
  summarise(
    lib_conf = first(lib_conf),  # Prioritize "Y" over "N" due to arrange
    source = if_else(
      n() > 1, 
      "both", 
      dplyr::first(source)
    ),
    .groups = "drop"
  )
#- 0c.5.5: Identify overlapping feature columns to handle duplicates properly
annot_features <- names(main_msmica_20p)[!names(main_msmica_20p) %in% c("Sample_ID", "pt_tissue_ID", "tissue_type", "tx", "time")]
confirmed_features <- names(TFT_confirmed)[!names(TFT_confirmed) %in% c("Sample_ID", "pt_tissue_ID", "tissue_type", "tx", "time")]
overlapping_features <- intersect(annot_features, confirmed_features)
#- 0c.5.6: Create base from annotations, removing overlapping features
TFT_combined_base <- main_msmica_20p %>%
  select(-any_of(overlapping_features))
#- 0c.5.7: Prepare confirmed features only (Patient + feature columns)
TFT_confirmed_features_only <- TFT_confirmed %>%
  select(Sample_ID, all_of(confirmed_features))
#- 0c.5.8: Merge datasets and apply uniqueness filter
TFT_combined_main <- TFT_combined_base %>%
  left_join(TFT_confirmed_features_only, by = "Sample_ID") %>%
  filter_unique_features(unique_threshold = 0.8)