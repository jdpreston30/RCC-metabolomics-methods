#' Reload All Utility Functions and Config
#'
#' Quick helper to re-source all utility functions from R/Utilities/ and reload config
#' Optionally run a numbered script after updating
#' Useful during development when making changes to utility functions or config
#'
#' @param n Optional script number to run after updating (e.g., 2 for 02_*.R)
#' @export
#' @examples
#' u() # Just update utilities and config
#' u(2) # Update utilities/config, then run script 02
u <- function(n = NULL) {
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

  # If script number provided, run that script
  if (!is.null(n)) {
    r(n)
  }

  invisible(NULL)
}

#' Update, Run Script, and Render Figures
#'
#' Convenience wrapper that runs u(n) then f()
#' Updates utilities/config, runs the specified script, then renders figures
#'
#' @param n Script number to run (e.g., 5 for 05_*.R)
#' @export
#' @examples
#' uf(5) # Equivalent to u(5) then f()
uf <- function(n) {
  u(n)
  f()
  invisible(NULL)
}

#' Run Render Figures Script
#'
#' Quick helper to run 05_assign_plots.R and 06_render_figures.R
#' Useful during development for quickly re-rendering figures
#'
#' @export
#' @examples
#' f()
f <- function() {
  # First run assign_figures
  assign_script <- "R/Scripts/05_assign_plots.R"
  if (file.exists(assign_script)) {
    cat("▶️  Running:", basename(assign_script), "\n")
    source(assign_script)
  } else {
    cat("⚠️  Script not found: ", assign_script, "\n")
    return(invisible(NULL))
  }

  # Then run render_figures
  render_script <- "R/Scripts/06_render_figures.R"
  if (file.exists(render_script)) {
    cat("▶️  Running:", basename(render_script), "\n")
    source(render_script)
    cat("✅ Figures rendered\n")
  } else {
    cat("⚠️  Script not found: ", render_script, "\n")
  }

  invisible(NULL)
}

#' Run Numbered Script
#'
#' Quick helper to run any numbered script by its number
#' Automatically pads single digit numbers with leading zero
#'
#' @param n Script number (e.g., 1 for 01_FTs.R, 5 for 05_targeted_volcano_diverge.R)
#' @export
#' @examples
#' r(1) # Runs 01_FTs.R
#' r(5) # Runs 05_targeted_volcano_diverge.R
#' r(9) # Runs 06_render_figures.R
r <- function(n) {
  # Pad single digit numbers with leading zero
  script_num <- sprintf("%02d", n)

  # Find script file that matches the pattern
  scripts_dir <- "R/Scripts/"
  pattern <- paste0("^", script_num, ".*\\.R$")

  matching_files <- list.files(scripts_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)

  if (length(matching_files) == 0) {
    cat("⚠️  No script found matching number:", n, "\n")
    return(invisible(NULL))
  }

  if (length(matching_files) > 1) {
    cat("⚠️  Multiple scripts found matching number:", n, "\n")
    cat("   ", paste(basename(matching_files), collapse = ", "), "\n")
    return(invisible(NULL))
  }

  script_path <- matching_files[1]
  cat("▶️  Running:", basename(script_path), "\n")
  source(script_path)
  cat("✅ Script completed\n")

  invisible(NULL)
}

#' @rdname u
#' @export
U <- u