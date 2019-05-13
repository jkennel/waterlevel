#' slug_hvorslev
#'
#' Response to a slug test
#'
#' @param elapsed_time elapsed time (numeric)
#' @param hydraulic_conductivity formation hydraulic conductivity (numeric)
#' @param screen_length length of the screen (numeric)
#' @param radius_well well radius (numeric)
#' @param radius_filter radius of the filter pack (numeric)
#'
#' @return hvorslev slugtest drawdown
#' 
#' @references Hvorslev, M. J., Time lag and soil permeability in ground water
#'  measurements, Bull. 36, 50 pp., U.S. Army Corps of Eng. Waterways Exp. Stn.,
#'  Vicksburg, Miss., 1951
#'  
#' @export
#'
#' 
slug_hvorslev <- function(elapsed_time, hydraulic_conductivity, screen_length,
                     radius_well, radius_filter) {
  
  exp((-2 * hydraulic_conductivity * screen_length * elapsed_time) /
        (radius_well^2 * log(screen_length / radius_filter^2)))
  
}