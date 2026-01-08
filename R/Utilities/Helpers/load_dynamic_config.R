#' Load and Resolve Dynamic Configuration Based on Computer Environment
#'
#' Loads a YAML configuration file and dynamically resolves computer-specific paths
#' and settings based on the detected or specified computer environment. Supports
#' automatic detection of laptop vs desktop environments and variable substitution
#' in path strings using template syntax.
#'
#' @param computer Character string specifying the computer environment:
#'   - "auto": Automatically detect based on username and system info (default)
#'   - "laptop": Use laptop-specific configuration (user: jdp2019)
#'   - "desktop": Use desktop-specific configuration (user: JoshsMacbook2015)
#' @param config_path Character string specifying the path to the YAML configuration file
#'   (default: "config_dynamic.yaml")
#'
#' @return Named list containing the resolved configuration with:
#'   - All computer-specific variables substituted in path strings
#'   - Computer detection metadata in \code{computer_used} field
#'   - Template variables like \code{{base_data_path}} resolved to actual paths
#'
#' @details
#' The function performs the following operations:
#' 1. Loads the raw YAML configuration file
#' 2. Auto-detects computer environment if requested using multiple methods:
#'    - Username detection (jdp2019 = laptop, JoshsMacbook2015 = desktop)
#'    - Computer name pattern matching
#'    - Path signature detection as fallback
#' 3. Extracts computer-specific variables from the configuration
#' 4. Recursively resolves template variables in all path strings
#' 5. Returns the final configuration with resolved paths
#'
#' @examples
#' \dontrun{
#'   # Auto-detect computer and load config
#'   config <- load_dynamic_config()
#'   
#'   # Manually specify computer type
#'   config <- load_dynamic_config(computer = "laptop")
#'   
#'   # Use custom config file
#'   config <- load_dynamic_config(config_path = "custom_config.yaml")
#' }
#'
#' @importFrom yaml read_yaml
#' @importFrom here here
#' @importFrom glue glue
#' @export
load_dynamic_config <- function(computer = "auto", config_path = "config_dynamic.yaml") {
  # Load raw configuration
  raw_config <- yaml::read_yaml(here::here(config_path))
  
  # Auto-detect computer if requested
  if (computer == "auto") {
    current_user <- Sys.getenv("USER")
    computer_name <- Sys.info()["nodename"]
    
    # Multiple detection methods for robustness
    if (current_user == "jdp2019") {
      computer <- "laptop"
      cat("🔍 Detected laptop via username:", current_user, "\n")
    } else if (current_user == "JoshsMacbook2015" || grepl("JoshsMacbook", computer_name) || grepl("JDP", computer_name)) {
      computer <- "desktop"
      cat("🔍 Detected desktop via username/computer name:", current_user, "/", computer_name, "\n")
    } else {
      # Fallback: check for specific path signatures
      if (dir.exists("/Users/jdp2019")) {
        computer <- "laptop"
        cat("🔍 Detected laptop via path signature\n")
      } else if (dir.exists("/Users/JoshsMacbook2015")) {
        computer <- "desktop"
        cat("🔍 Detected desktop via path signature\n")
      } else {
        stop("Could not auto-detect computer. Available: 'laptop' (jdp2019) or 'desktop' (JoshsMacbook2015). ",
             "Current user: ", current_user, ", Computer: ", computer_name,
             "\nPlease specify computer = 'laptop' or computer = 'desktop'")
      }
    }
  }
  
  # Validate computer selection
  if (!computer %in% names(raw_config$computers)) {
    stop("Invalid computer selection. Available options: ", 
         paste(names(raw_config$computers), collapse = ", "))
  }
  
  cat("🖥️ Using configuration for:", computer, "\n")
  
  # Get computer-specific variables
  comp_vars <- raw_config$computers[[computer]]
  
  # Create substitution variables including derived ones
  substitution_vars <- comp_vars
  substitution_vars$base_data_path <- glue::glue(
    raw_config$paths$base_data_path, 
    .envir = list2env(comp_vars)
  )
  
  # Recursively substitute variables in all path strings
  resolve_paths <- function(obj, vars) {
    if (is.list(obj)) {
      return(lapply(obj, resolve_paths, vars))
    } else if (is.character(obj) && length(obj) == 1) {
      # Only substitute if string contains template variables
      if (grepl("\\{.*\\}", obj)) {
        return(as.character(glue::glue(obj, .envir = list2env(vars))))
      } else {
        return(obj)
      }
    } else {
      return(obj)
    }
  }
  
  # Resolve all paths
  resolved_config <- raw_config
  resolved_config$paths <- resolve_paths(raw_config$paths, substitution_vars)
  
  # Remove the computers section from final config
  resolved_config$computers <- NULL
  
  # Add metadata about which computer was used
  resolved_config$computer_used <- computer
  
  return(resolved_config)
}

#' Resolve File Paths from Config Templates
#'
#' Resolves parameterized file path templates from the configuration using
#' data version (ratio or main) and template name. Substitutes {data_path},
#' {metabojanitor_subpath}, {user_home}, and {onedrive_path} variables.
#'
#' @param config List: Configuration object returned by load_dynamic_config()
#' @param template_name Character: Name of the file template (e.g., "MSMICA_20p", "UFT_full")
#' @param version Character: Data version to use - "ratio" or "main" (default: "main")
#'
#' @return Character string containing the resolved file path
#'
#' @examples
#' \dontrun{
#'   config <- load_dynamic_config()
#'   msmica_path <- resolve_file(config, "MSMICA_20p", "main")
#'   uft_ratio <- resolve_file(config, "UFT_full", "ratio")
#' }
#'
#' @importFrom glue glue
#' @export
resolve_file <- function(config, template_name, version = "main") {
  # Get the template
  if (!template_name %in% names(config$file_templates)) {
    stop("Template '", template_name, "' not found. Available: ",
         paste(names(config$file_templates), collapse = ", "))
  }
  
  template <- config$file_templates[[template_name]]
  
  # Get data path key from data_versions
  if (!version %in% names(config$data_versions)) {
    stop("Version '", version, "' not found. Available: ",
         paste(names(config$data_versions), collapse = ", "))
  }
  
  data_path_key <- config$data_versions[[version]]$data_path
  
  # Get computer-specific variables from paths section
  # (these were stored during config resolution)
  comp_vars <- list(
    data_path = config$paths[[data_path_key]],
    metabojanitor_subpath = "MetaboJanitoR"  # Default, can be overridden if stored in config
  )
  
  # Resolve the template with substitution variables
  resolved_path <- glue::glue(template, .envir = list2env(comp_vars))
  
  return(resolved_path)
}