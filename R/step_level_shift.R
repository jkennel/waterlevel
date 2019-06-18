#' Gaps where there are shifts in the regression
#'
#' `step_level_shift` creates a *specification* of a recipe step
#'  for creating level shifts
#'
#' @inheritParams recipes::step_lag
#' @inheritParams find_level_shift
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Defaults to "level_shift"
#' @param prefix A prefix for generated column names, default to 
#'  "level_shift_".
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any). For the
#'  `tidy` method, a tibble with columns `terms` which is
#'  the columns that will be affected and `holiday`.
#' @keywords datagen
#' @concept generate a set of level shifts
#' @export
#' @details `step_level_shift` calculates the earthtide response and then
#'  lags (leads) the terms.
#' @examples
#' data(transducer)
#'
#' transducer[1000:1200, wl := NA_real_]
#' 
#' rec <- recipe(wl ~ .,
#'               data = transducer[, list(datetime, wl, baro)])
#'
#' with_levels <- rec %>%
#'   step_level_shift(wl, datetime, time_interval = 120L) %>%
#'   prep() %>%
#'   juice()
#'
#' @seealso [recipe()]
#'   [prep.recipe()] [bake.recipe()]
#' @importFrom recipes add_step step terms_select ellipse_check rand_id
step_level_shift <-
  function(recipe,
           ...,
           role = "level_shift",
           trained = FALSE,
           time_interval = 1L,
           prefix = "level_shift_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("level_shift")) {
    add_step(
      recipe,
      step_level_shift_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        time_interval = time_interval,
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_level_shift_new <-
  function(terms, role, trained, dep_var, time_var, time_interval, default, prefix, columns, skip, id) {
    step(
      subclass = "level_shift",
      terms = terms,
      role = role,
      trained = trained,
      time_interval = time_interval,
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_level_shift <- function(x, training, info = NULL, ...) {
  
  step_level_shift_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    time_interval = x$time_interval,
    default = x$default,
    prefix = x$prefix,
    columns = terms_select(x$terms, info = info),
    skip = x$skip,
    id = x$id
  )
  
}

#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#' @importFrom recipes bake prep
#' @export
bake.step_level_shift <- function(object, new_data, ...) {
  
  as_tibble(mutate(new_data, 
                   level_shift = add_level_shift(new_data,
                                                 as.character(object$columns[1]),
                                                 as.character(object$columns[2]),
                                                 object$time_interval)))
  
  
}


#' @importFrom recipes printer
print.step_level_shift <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("level_shift ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }