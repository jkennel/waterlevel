#' Lagged Earth tide response
#'
#' `step_earthtide` creates a *specification* of a recipe step
#'  that are the Earth tide harmonics for a particular 
#'  location.  Requires the *earthtide* package
#'
#' @inheritParams recipes::step_lag
#' @inheritParams earthtide::calc_earthtide
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Defaults to "earthtide"
#' @param prefix A prefix for generated column names, default to 
#'  "earthtide_".
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any). For the
#'  `tidy` method, a tibble with columns `terms` which is
#'  the columns that will be affected and `holiday`.
#' @keywords datagen
#' @concept generate Earth tide harmonics
#' @export
#' 
#' 
#' @details `step_earthtide` calculates the Earth tide harmonics
#' @examples
#' library(earthtide)
#' data(transducer)
#' data(eterna_wavegroups)
#' 
#' rec <- recipe(wl ~ .,
#'               data = transducer[1:1000, list(datetime, wl)])
#'
#' wg <- na.omit(eterna_wavegroups[eterna_wavegroups$time == '1 month',])
#' with_et <- rec %>%
#'   step_earthtide(datetime, 
#'                  latitude = 34, 
#'                  longitude = -118.5, 
#'                  wave_groups = wg) %>% 
#'   prep() %>%
#'   juice()
#'
#' @seealso [step_earthtide()] [recipe()]
#'   [prep.recipe()] [bake.recipe()]
#' @importFrom recipes add_step step terms_select ellipse_check rand_id
step_earthtide <-
  function(recipe,
           ...,
           role = "earthtide",
           trained = FALSE,
           method = "gravity",
           astro_update = 1L,
           latitude = 0,
           longitude = 0,
           elevation = 0,
           azimuth = 0,
           gravity = 0,
           earth_radius = 6378136.3,
           earth_eccen = 0.0066943979514,
           cutoff = 1e-06, 
           wave_groups = NULL,
           catalog = "ksm04", 
           eop = NULL,
           scale = TRUE,
           prefix = "earthtide_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("earthtide")) {
    add_step(
      recipe,
      step_earthtide_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        method = method,
        astro_update = astro_update,
        latitude = latitude,
        longitude = longitude,
        elevation = elevation,
        azimuth = azimuth,
        gravity = gravity,
        earth_radius = earth_radius,
        earth_eccen = earth_eccen,
        cutoff = cutoff, 
        wave_groups = wave_groups,
        catalog = catalog, 
        eop = eop,
        scale = scale,
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_earthtide_new <-
  function(terms, role, trained, method, astro_update, latitude, longitude,
           elevation, azimuth, gravity, earth_radius, earth_eccen, cutoff,
           wave_groups, catalog, eop, scale, default, prefix, columns, skip, id) {
    step(
      subclass = "earthtide",
      terms = terms,
      role = role,
      trained = trained,
      method = method,
      astro_update = astro_update,
      latitude = latitude,
      longitude = longitude,
      elevation = elevation,
      azimuth = azimuth,
      gravity = gravity,
      earth_radius = earth_radius,
      earth_eccen = earth_eccen,
      cutoff = cutoff, 
      wave_groups = wave_groups,
      catalog = catalog, 
      eop = eop,
      scale = scale,
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_earthtide <- function(x, training, info = NULL, ...) {
  
  step_earthtide_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    method = x$method,
    astro_update = x$astro_update,
    latitude = x$latitude,
    longitude = x$longitude,
    elevation = x$elevation,
    azimuth = x$azimuth,
    gravity = x$gravity,
    earth_radius = x$earth_radius,
    earth_eccen = x$earth_eccen,
    cutoff = x$cutoff, 
    wave_groups = x$wave_groups,
    catalog = x$catalog, 
    eop = x$eop,
    scale = x$scale,
    default = x$default,
    prefix = x$prefix,
    columns = terms_select(x$terms, info = info),
    skip = x$skip,
    id = x$id
  )
  
}

#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#' @importFrom earthtide calc_earthtide
#' @importFrom recipes bake prep
#' 
#' @export
bake.step_earthtide <- function(object, new_data, ...) {

  bind_cols(new_data, 
            as_tibble(calc_earthtide(utc = new_data[[object$columns]],
                                     do_predict = FALSE,
                                     method = object$method,
                                     astro_update = object$astro_update, 
                                     latitude = object$latitude,
                                     longitude = object$longitude,
                                     elevation = object$elevation, 
                                     azimuth = object$azimuth,
                                     gravity = object$gravity, 
                                     earth_radius = object$earth_radius,
                                     earth_eccen = object$earth_eccen,
                                     cutoff = object$cutoff,
                                     wave_groups = object$wave_groups,
                                     catalog = object$catalog,
                                     eop = object$eop,
                                     scale = object$scale,
                                     return_matrix = TRUE)))
  
}

#' @importFrom recipes printer
print.step_earthtide <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("earthtide ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }