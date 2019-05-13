#' Generate harmonics
#'
#' `step_harmonic` creates a *specification* of a recipe step
#'  for harmonics (sin and cos terms).
#'
#' @inheritParams recipes::step_lag
#' @inheritParams harmonic
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Defaults to "harmonic"
#' @param prefix A prefix for generated column names, default to 
#'  "harmonic_".
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any). For the
#'  `tidy` method, a tibble with columns `terms` which is
#'  the columns that will be affected and `holiday`.
#' @keywords datagen
#' @concept generate harmonics of frequency in cycles per day
#' @export
#' @details `step_harmonic` calculates the Earth tide harmonics
#' @examples
#' data(transducer)
#' 
#' rec <- recipe(wl ~ .,
#'               data = transducer[1:1000, list(datetime, wl)])
#'
#' with_et <- rec %>%
#'   step_harmonic(datetime, freq = c(1, 1.93, 2)) %>% 
#'   prep() %>%
#'   juice()
#'
#' @seealso [step_earthtide()] [recipe()]
#'   [prep.recipe()] [bake.recipe()]
#' @importFrom recipes add_step step terms_select ellipse_check rand_id
step_harmonic <-
  function(recipe,
           ...,
           role = "harmonic",
           trained = FALSE,
           freq = 1,
           prefix = "harmonic_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("harmonic")) {
    add_step(
      recipe,
      step_harmonic_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        freq = freq,
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_harmonic_new <-
  function(terms, role, trained, freq, default, prefix, columns, skip, id) {
    step(
      subclass = "harmonic",
      terms = terms,
      role = role,
      trained = trained,
      freq = freq,
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_harmonic <- function(x, training, info = NULL, ...) {
  
  step_harmonic_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    freq = x$freq,
    default = x$default,
    prefix = x$prefix,
    columns = terms_select(x$terms, info = info),
    skip = x$skip,
    id = x$id
  )
  
}

#' @importFrom dplyr bind_cols
#' @importFrom recipes bake prep
#' 
#' @export
bake.step_harmonic <- function(object, new_data, ...) {
  

  # if (!all(object$lag == as.integer(object$lag)))
  #   stop("step_lag requires 'lag' argument to be integer valued.",
  #        call. = FALSE)
  
  # make_call_sin <- function(col, freq) {
  #   lang(
  #     "sin_harmonic",
  #     x = sym(col),
  #     freq = freq,
  #     .ns = "waterlevel"
  #   )
  # }
  # 
  # make_call_cos <- function(col, freq) {
  #   lang(
  #     "cos_harmonic",
  #     x = sym(col),
  #     freq = freq,
  #     .ns = "waterlevel"
  #   )
  # }
  # 
  # 
  # grid <- expand.grid(col = object$columns, freq = object$freq,
  #                     stringsAsFactors = FALSE)
  # calls_sin <- as.vector(map2(grid$col, grid$freq, make_call_sin))
  # new_name_sin <- paste0(object$prefix, 'sin_', grid$freq, "_", grid$col)
  # calls_sin <- check_name(calls_sin, new_data, object, new_name_sin, TRUE)
  # 
  # calls_cos <- as.vector(map2(grid$col, grid$freq, make_call_cos))
  # new_name_cos <- paste0(object$prefix, 'cos_', grid$freq, "_", grid$col)
  # calls_cos <- check_name(calls_cos, new_data, object, new_name_cos, TRUE)
  # 
  # 
  # as_tibble(mutate(new_data, !!!calls_sin) %>%
  #           mutate(!!!calls_cos))
  
  bind_cols(new_data, 
            as_tibble(harmonic(new_data[[object$columns]], object$freq)))
  
}


#' @importFrom recipes printer
print.step_harmonic <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("harmonic ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }