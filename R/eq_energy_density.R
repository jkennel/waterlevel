
#' energy_density_tsuboi
#'
#' @param magnitude Richter magnitude 
#' @param distance epicentral distance
#'
#' @return
#' @export
#'
#' @examples
#' energy_density_tsuboi(5, 1000)
energy_density_tsuboi <- function(magnitude, distance) {
  # https://en.wikipedia.org/wiki/Richter_magnitude_scale
  # Tsuboi, University of Tokyo empirical
  10^(magnitude + 0.83 - 1.73 * log10(distance))
}


#' energy_density_lillie
#'
#' @param magnitude Richter magnitude 
#' @param distance epicentral distance
#'
#' @return
#' @export
#'
#' @examples
#' energy_density_lillie(5, 1000)
energy_density_lillie <- function(magnitude, distance) {
  # https://en.wikipedia.org/wiki/Richter_magnitude_scale
  # Lillie empirical
  10^(magnitude + 2.48 - 2.76 * log10(distance))
}




