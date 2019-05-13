#' level_shift
#'
#' Create a factor column with the data subseted according to break points.
#' 
#' @param x vector of numeric
#' @param level_breaks break points in the data
#'
#' @return vector of type factor 
#' @export
#'
#' @examples
#' level_shift(1:100, c(2, 40, 86))
level_shift <- function(x,
                        level_breaks = NULL) {
  
  cut(x, 
      breaks = c(min(x, na.rm = TRUE),
                 level_breaks,
                 max(x, na.rm = TRUE)),
      include.lowest = TRUE)
  
}

