# GitHub Copilot Instructions

**R Code Style & Project Structure Preferences:**

## Comments & Spacing
- **Hierarchical comment structure**: `#*` (major section), `#+` (subsection), `#-` (detail/sub-subsection), `#_` (individual item), `#!` (important note)
- **Numbering format**: 
  - Major sections: `#* 9: Validation Plots Adjustment`
  - Subsections: `#+ 9.2: Plots that failed after review`
  - Details: `#- 9.2.1: Methylparaben (CP2252)`
  - Items can use colon format: `#_ Compound Name (ID)`
- **Inline comments**: Regular explanatory comments within code blocks just use `#` without special symbols
- **Comments are brief and descriptive**: Typically compound names, short phrases, or action descriptions
- **Compact code style**: NO extra blank lines between sections or subsections - code should be dense
- **No blank lines after section headers**: Code immediately follows comment headers
- **Example structure**:
  ```r
  #* 9: Major Section
  if (!isTRUE(config$skip_something)) {
  #+ 9.2: Subsection
  #- 9.2.1: Compound Name (ID)
  variable <- function_call(args)
  # Regular inline comment explaining something specific
  another_variable <- function_call(args)
  #- 9.2.2: Next Compound (ID)
  more_code <- function_call(args)
  }
  ```

## Project Architecture
- Use YAML configuration files for dynamic, computer-specific paths (e.g., `config_dynamic.yaml`)
- Store utility functions in organized subdirectories: `R/Utilities/Helpers/`, `R/Utilities/Validation/`, etc.
- Use renv for package management with automatic restore on first run
- Pipeline structure: Main run script (`All_Run/run.R`) sources numbered scripts sequentially
- Load all utility functions at setup with: `purrr::walk(list.files(...), source)`
- Store configuration in `.GlobalEnv$config` or similar for global access
- Use `load_dynamic_config()` pattern for automatic computer detection and path resolution

## Functions
- Use roxygen2-style documentation for utility functions
- Prefer tidyverse functions when appropriate
- Functions should respect YAML configuration flags (e.g., skip logic based on config parameters)
- Separate concerns: e.g., plot generation vs PDF compilation in separate functions

## File Organization
- Scripts: Numbered for pipeline order (00a, 00b, 01, 02, etc.)
- Outputs: Organized in `Outputs/Figures/`, `Outputs/Tables/`, `Outputs/Validation/`
- Metadata: Static reference files in `metadata_files/`
- Use OneDrive for final output storage, local temp dirs for run-specific I/O
