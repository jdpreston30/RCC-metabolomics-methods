{
source("R/Utilities/Helpers/restore_renv.R")
source("R/Utilities/Helpers/load_dynamic_config.R")
config <- load_dynamic_config(computer = "auto", config_path = "All_Run/config_dynamic.yaml")
source("R/Scripts/00a_environment_setup.R")
source("R/Scripts/00b_setup.R")
source("R/Scripts/00c_FTs.R")
source("R/Scripts/01_ratio_exp.R")
source("R/Scripts/02_untargeted_TC.R")
source("R/Scripts/09_assign_plots.R")
source("R/Scripts/10_render_figures.R")
}