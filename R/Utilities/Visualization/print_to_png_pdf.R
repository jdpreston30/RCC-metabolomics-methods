#' Print Plot to PNG and PDF with Auto-Refresh
#'
#' Saves ggplot objects as both high-resolution PNG and PDF files with automatic 
#' Preview.app integration on macOS. This function combines PNG and PDF output
#' in a single call, using magick to convert from PNG to PDF at high quality.
#'
#' @param plot ggplot object or other plot object compatible with ggplot2::ggsave()
#' @param filename Character string for the base filename (without extension).
#'   The function will create both .png and .pdf versions
#' @param width Numeric value for plot width in inches (default: 8.5 for letter size)
#' @param height Numeric value for plot height in inches (default: 11 for letter size)  
#' @param dpi Numeric value for resolution in dots per inch (default: 600 for publication quality)
#' @param output_dir_png Character string specifying the PNG output directory. If NULL (default),
#'   uses "Supplementary/Components/Figures/PNG"
#' @param output_dir_pdf Character string specifying the PDF output directory. If NULL (default),
#'   uses "Supplementary/Components/Figures/PDF"
#' @param auto_open Logical indicating whether to automatically open the PNG in
#'   Preview.app on first save (default: TRUE)
#'
#' @return Invisibly returns a list with PNG and PDF file paths
#'
#' @details
#' This function streamlines the workflow for creating both PNG and PDF versions
#' of supplementary figures. The PNG is created first using print_to_png(), then
#' converted to PDF using magick at high density (600 DPI default).
#'
#' @examples
#' \dontrun{
#'   # Basic usage for supplementary figure
#'   print_to_png_pdf(sup_fig1, "S1")
#'   
#'   # Custom output directories
#'   print_to_png_pdf(sup_fig2, "S2", 
#'                    output_dir_png = "Outputs/Figures",
#'                    output_dir_pdf = "Outputs/Figures")
#' }
#'
#' @importFrom ggplot2 ggsave
#' @importFrom magick image_read image_write
#' @export
print_to_png_pdf <- function(plot, filename, width = 8.5, height = 11, dpi = 600,
                              output_dir_png = "Supplementary/Components/Figures/PNG",
                              output_dir_pdf = "Supplementary/Components/Figures/PDF",
                              auto_open = TRUE) {
  
  # Ensure filename doesn't have extension
  filename_base <- sub("\\.(png|pdf)$", "", filename, ignore.case = TRUE)
  
  # Create PNG filename
  png_filename <- paste0(filename_base, ".png")
  pdf_filename <- paste0(filename_base, ".pdf")
  
  # Print to PNG using existing function
  png_path <- print_to_png(plot, png_filename, 
                           width = width, height = height, dpi = dpi,
                           output_dir = output_dir_png, auto_open = auto_open)
  
  # Setup PDF output directory
  if (!dir.exists(output_dir_pdf)) {
    dir.create(output_dir_pdf, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create PDF path
  pdf_path <- file.path(output_dir_pdf, pdf_filename)
  
  # Convert PNG to PDF using magick
  image_write(image_read(png_path), pdf_path, format = "pdf", density = dpi)
  
  # Return both paths invisibly
  invisible(list(png = png_path, pdf = pdf_path))
}
