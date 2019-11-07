
#' energy_density_wang
#'
#' @param magnitude Richter magnitude 
#' @param distance epicentral distance
#'
#' @return value of the energy density
#' @export
#'
#' @examples
#' energy_density_wang(5, 9)
energy_density_wang <- function(magnitude, distance) {
  # Wang & Manga 2010
  10^((log10(distance) + 1.4 - (0.48 * magnitude)) / -0.33)
}




