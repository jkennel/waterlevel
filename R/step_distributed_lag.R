#' Distributed lag response
#'
#' `step_distributed_lag` creates a *specification* of a recipe step
#'  that are distributed lag versions of a particular variable. Uses FFT for 
#'  fast calculation with a large maximum lag and many observations
#'
#' @inheritParams recipes::step_lag
#' @inheritParams distributed_lag
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Defaults to "distributed_lag"
#' @param prefix A prefix for generated column names, default to 
#'  "distributed_lag_".
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any). For the
#'  `tidy` method, a tibble with columns `terms` which is
#'  the columns that will be affected and `holiday`.
#' @keywords datagen
#' @concept generate a distributed lag response
#' @export
#' @details `step_distributed_lag` calculates the earthtide response and then
#'  lags (leads) the terms.
#' @examples
#' data(transducer)
#'
#' rec <- recipe(wl ~ .,
#'               data = transducer[1:1000, list(datetime, wl, baro)])
#'
#' with_et <- rec %>%
#'   step_distributed_lag(baro, knots = c(1, 10, 100)) %>%
#'   step_naomit(everything()) %>% 
#'   prep() %>%
#'   juice()
#'
#' @seealso [step_lag_matrix()] [recipe()]
#'   [prep.recipe()] [bake.recipe()]
#' @importFrom recipes add_step step terms_select ellipse_check rand_id
step_distributed_lag <-
  function(recipe,
           ...,
           role = "distributed_lag",
           trained = FALSE,
           knots = 1,
           spline_fun = splines::ns,
           prefix = "distributed_lag_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("distributed_lag")) {
    add_step(
      recipe,
      step_distributed_lag_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        knots = knots,
        spline_fun = spline_fun,
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_distributed_lag_new <-
  function(terms, role, trained, knots, spline_fun, default, prefix, columns, skip, id) {
    step(
      subclass = "distributed_lag",
      terms = terms,
      role = role,
      trained = trained,
      knots = knots,
      spline_fun = spline_fun,
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_distributed_lag <- function(x, training, info = NULL, ...) {
  
  step_distributed_lag_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    knots = x$knots,
    spline_fun = x$spline_fun,
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
bake.step_distributed_lag <- function(object, new_data, ...) {
  
  
  
  bind_cols(new_data, 
            as_tibble(distributed_lag(new_data[[object$columns]],
                                      object$knots,
                                      object$spline_fun, 
                                      object$columns)))
  
}


#' @importFrom recipes printer
print.step_harmonic <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("harmonic ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }