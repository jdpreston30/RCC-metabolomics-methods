#* 0c: Feature Table Import and Preprocessing
#+ 0c.1: Import feature tables and metadata for ratio experiment
#- 0c.1.1: MSMICA feature tables (ratio)
ratio_msmica_20p <- read_csv(resolve_file(config, "MSMICA_20p", "ratio"))
ratio_msmica_key <- read_csv(resolve_file(config, "MSMICA_key", "ratio"))
#- 0c.1.2: Untargeted feature tables (ratio)
ratio_uft_full <- read_csv(resolve_file(config, "UFT_full", "ratio"))
ratio_uft_20p <- read_csv(resolve_file(config, "UFT_20p", "ratio"))
#+ 0c.2: Import feature tables and metadata for main experiment
#- 0c.2.1: MSMICA feature tables (main)
main_msmica_20p <- read_csv(resolve_file(config, "MSMICA_20p", "main"))
main_msmica_key <- read_csv(resolve_file(config, "MSMICA_key", "main"))
#- 0c.2.2: Untargeted feature tables (main)
main_uft_full <- read_csv(resolve_file(config, "UFT_full", "main"))
main_uft_20p <- read_csv(resolve_file(config, "UFT_20p", "main"))