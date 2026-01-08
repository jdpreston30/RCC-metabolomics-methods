#' Figure Labels Generator for Cowplot Layouts
#'
#' Generates figure panel labels (A, B, C, etc.) at specified coordinates
#' for multi-panel figure layouts using cowplot's coordinate system.
#'
#' @param labels Named list where names are label text and values are coordinate vectors c(x, y)
#' @param size Font size for labels (default: 14)
#' @param fontface Font face for labels (default: "bold")
#' @param fontfamily Font family for labels (default: "Arial")
#' @param hjust Horizontal justification (default: 0)
#'
#' @return List of cowplot::draw_label objects
#'
#' @details
#' Creates publication-ready figure panel labels at precise coordinates.
#' Coordinates use cowplot's system where (0,0) is bottom-left corner.
#'
#' Common positioning:
#' - Top-left panels: around (0.8, 9.7)
#' - Top-right panels: around (4.3, 9.7)
#' - Bottom panels: adjust y-coordinate accordingly
#'
#' @examples
#' \dontrun{
#' # Define label positions
#' labels <- list(
#'   A = c(0.8, 9.7),
#'   B = c(4.3, 9.7),
#'   C = c(0.8, 5.2)
#' )
#'
#' # Add to cowplot layout
#' final_plot + figure_labels(labels)
#' }
#'
#' @export
#+ Figure Labels Generator
figure_labels <- function(labels, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0) {
  # Convert single label to list format if needed
  if (is.character(labels)) {
    stop("Please provide labels as a named list with x and y coordinates, e.g., list(A = c(0.8, 9.7), B = c(3.7, 9.7))")
  }

  # Generate draw_label calls for each label
  label_layers <- list()
  for (name in names(labels)) {
    coords <- labels[[name]]
    label_layers[[length(label_layers) + 1]] <-
      cowplot::draw_label(name,
        x = coords[1], y = coords[2],
        size = size, fontface = fontface,
        fontfamily = fontfamily, hjust = hjust
      )
  }

  return(label_layers)
}

#' Grid Guide for Positioning (Cowplot Coordinates)
#'
#' Creates a grid overlay with dashed lines and coordinate labels to assist with
#' precise positioning of plot elements in cowplot figure layouts.
#'
#' @param x_max Maximum x-coordinate for the grid (default: 8.5)
#' @param y_max Maximum y-coordinate for the grid (default: 11)
#' @param interval Spacing between grid lines (default: 0.25)
#' @param label_interval Spacing between coordinate labels (default: 0.5)
#' @param margins Distance from edges for margin lines in inches (default: 0.5)
#'
#' @return List of ggplot2 geom objects for grid overlay
#'
#' @details
#' The function creates a coordinate grid system for precise element positioning:
#' - Red dashed lines at specified intervals
#' - Coordinate labels at specified intervals
#' - Optional black margin lines to show printable area
#'
#' Coordinates follow cowplot's system where (0,0) is bottom-left.
#'
#' @examples
#' \dontrun{
#' # Add grid to existing plot
#' plot + grdgd()
#'
#' # Custom grid settings
#' plot + grdgd(x_max = 10, y_max = 8, interval = 0.5)
#' }
#'
#' @export
#+ Grid Guide Generator
grdgd <- function(x_max = 8.5, y_max = 11, interval = 0.25, label_interval = 0.5, margins = 0.5) {
  guide_elements <- list(
    ggplot2::geom_vline(xintercept = seq(0, x_max, interval), color = "red", alpha = 0.3, linetype = "dashed"),
    ggplot2::geom_hline(yintercept = seq(0, y_max, interval), color = "red", alpha = 0.3, linetype = "dashed"),
    ggplot2::annotate("text", x = seq(0, x_max, label_interval), y = 0.2, label = seq(0, x_max, label_interval), size = 3, color = "red"),
    ggplot2::annotate("text", x = 0.2, y = seq(0, y_max, label_interval), label = seq(0, y_max, label_interval), size = 3, color = "red")
  )
  # Add margin lines if margins is specified
  if (!is.null(margins)) {
    # Add vertical margin lines (left and right)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_vline(xintercept = margins, color = "black", linewidth = 1)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_vline(xintercept = x_max - margins, color = "black", linewidth = 1)

    # Add horizontal margin lines (bottom and top)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_hline(yintercept = margins, color = "black", linewidth = 1)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_hline(yintercept = y_max - margins, color = "black", linewidth = 1)
  }
  return(guide_elements)
}
