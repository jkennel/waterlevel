#' Lagged Earth tide response
#'
#' `step_lag_earthtide` creates a *specification* of a recipe step
#'  that are lagged versions of the Earth tide response at a particular 
#'  location.  Requires the *earthtide* package
#'
#' @inheritParams recipes::step_lag
#' @inheritParams earthtide::calc_earthtide
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Defaults to "lag_earthtide"
#' @param prefix A prefix for generated column names, default to 
#'  "lag_earthtide_".
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any). For the
#'  `tidy` method, a tibble with columns `terms` which is
#'  the columns that will be affected and `holiday`.
#' @keywords datagen
#' @concept generate a lagged form of the Earth tide response
#' @export
#' @details `step_lag_earthtide` calculates the earthtide response and then
#'  lags (leads) the terms.
#' @examples
#' data(transducer)
#'
#' rec <- recipe(wl ~ .,
#'               data = transducer[1:1000, list(datetime, wl)])
#'
#' with_et <- rec %>%
#'   step_lag_earthtide(datetime, latitude = 34, longitude = -118.5, lag = -1:1) %>% 
#'   prep() %>%
#'   juice()
#'
#' @seealso [step_lag_earthtide()] [recipe()]
#'   [prep.recipe()] [bake.recipe()]
#' @importFrom recipes add_step step terms_select ellipse_check rand_id
step_lag_earthtide <-
  function(recipe,
           ...,
           role = "lag_earthtide",
           trained = FALSE,
           method = "gravity",
           lag = 0,
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
           prefix = "lag_earthtide_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("lag_earthtide")) {
    add_step(
      recipe,
      step_lag_earthtide_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        lag = lag,
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
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_lag_earthtide_new <-
  function(terms, role, trained,lag, method, astro_update, latitude, longitude,
           elevation, azimuth, gravity, earth_radius, earth_eccen, cutoff,
           wave_groups, catalog, eop, default, prefix, columns, skip, id) {
    step(
      subclass = "lag_earthtide",
      terms = terms,
      role = role,
      trained = trained,
      lag = lag,
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
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_lag_earthtide <- function(x, training, info = NULL, ...) {
  
  step_lag_earthtide_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    lag = x$lag,
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
    default = x$default,
    prefix = x$prefix,
    columns = terms_select(x$terms, info = info),
    skip = x$skip,
    id = x$id
  )
  
}

#' @importFrom dplyr mutate
#' @importFrom rlang call2 sym
#' @importFrom purrr map2
#' @importFrom earthtide calc_earthtide
#' @importFrom recipes bake prep check_name
#' 
#' @export
bake.step_lag_earthtide <- function(object, new_data, ...) {
  
  earthtide_to_lag = calc_earthtide(utc = new_data[[object$columns]],
                 do_predict = TRUE,
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
                 return_matrix = TRUE)
  
  make_call <- function(col, lag_val) {
    call2(
      "shift_subset",
      x = sym(col),
      lag = lag_val,
      n_subset = 1L,
      n_shift = 0L,
      .ns = "waterlevel"
    )
  }
  
  
  grid    <- expand.grid(col = 'earthtide_to_lag', lag_val = object$lag,
                         stringsAsFactors = FALSE)
  calls   <- as.vector(map2(grid$col, grid$lag_val, make_call))
  newname <- paste0(object$prefix, grid$lag_val)
  newname <- gsub('-', 'n', newname)
  calls   <- check_name(calls, new_data, object, newname, TRUE)
  
  as_tibble(mutate(new_data, !!!calls))
  
  
}

#' @importFrom recipes printer
print.step_lag_earthtide <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("lag_earthtide ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }