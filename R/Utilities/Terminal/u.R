#' Reload All Utility Functions and Config
#'
#' Quick helper to re-source all utility functions from R/Utilities/ and reload config
#' Useful during development when making changes to utility functions or config
#'
#' @export
#' @examples
#' u()  # or U()
u <- function() {
  utils_path <- "R/Utilities/"
  if (dir.exists(utils_path)) {
    purrr::walk(
      list.files(utils_path, pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE),
      source
    )
    cat("✅ Utility functions updated\n")
  } else {
    cat("⚠️  Utilities directory not found: ", utils_path, "\n")
  }
  
  #- Reload config
  config_path <- "All_Run/config_dynamic.yaml"
  if (file.exists(config_path)) {
    config <<- load_dynamic_config(computer = "auto", config_path = config_path)
    .GlobalEnv$config <- config
    cat("✅ Config reloaded\n")
  } else {
    cat("⚠️  Config file not found: ", config_path, "\n")
  }
  
  invisible(NULL)
}

#' @rdname u
#' @export
U <- u
