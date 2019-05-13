#' Create a lagged (or lead) predictor
#'
#' `step_lag_matrix` creates a *specification* of a recipe step that
#'   will add new columns of lagged data. Lagged data will
#'   by default include NA values where the lag was induced.
#'   These can be removed with [step_naomit()], or you may
#'   specify an alternative filler value with the `default`
#'   argument.  This method is faster than [step_lag()] and allows
#'   for negative values.
#'   
#' @inheritParams lag_matrix
#' @inheritParams recipes::step_lag
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See [selections()] for more details.
#' @param role Defaults to "predictor"
#' @param trained A logical to indicate if the quantities for preprocessing
#'   have been estimated.
#' @param lag A vector of integers. They can be positive, negative or zero.
#'  Each specified column will be lagged for each value in the vector.
#' @param n_subset subset every n_subset values 
#' @param n_shift shift the data n_shift values
#' @param prefix A prefix for generated column names, default to "lag_".
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#' @param id A character string that is unique to this step to identify it.
#' @param skip A logical. Should the step be skipped when the
#'  recipe is baked by [bake.recipe()]? While all operations are baked
#'  when [prep.recipe()] is run, some operations may not be able to be
#'  conducted on new data (e.g. processing the outcome variable(s)).
#'  Care should be taken when using `skip = TRUE` as it may affect
#'  the computations for subsequent operations
#' @return An updated version of `recipe` with the
#'   new step added to the sequence of existing steps (if any).
#' @details The step assumes that the data are already _in the proper sequential
#'  order_ for lagging.
#' @export
#' @rdname step_lag_matrix
#'
#' @examples
#' data(transducer)
#'
#' rec <- recipe(wl ~ .,
#'               data = transducer[1:1000, list(datetime, wl, baro)])
#'
#' with_et <- rec %>%
#'   step_lag_matrix(baro, lag = -1:1) %>%
#'   step_naomit(everything()) %>% 
#'   prep() %>%
#'   juice()
#'
#' @seealso [recipe()] [step_lag()] [prep.recipe()] [bake.recipe()]
#'          [step_naomit()]
#'          
#' @importFrom recipes add_step step terms_select ellipse_check rand_id
step_lag_matrix <-
  function(recipe,
           ...,
           role = "lag_matrix",
           trained = FALSE,
           lag = 1,
           n_subset = 1,
           n_shift = 0,
           prefix = "lag_matrix_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("lag_matrix")) {
    add_step(
      recipe,
      step_lag_matrix_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        lag = lag,
        n_subset = n_subset,
        n_shift = n_shift,
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_lag_matrix_new <-
  function(terms, role, trained, lag, n_subset, n_shift, 
           default, prefix, columns, skip, id) {
    step(
      subclass = "lag_matrix",
      terms = terms,
      role = role,
      trained = trained,
      lag = lag,
      n_subset = n_subset,
      n_shift = n_shift,
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_lag_matrix <- function(x, training, info = NULL, ...) {
  
  step_lag_matrix_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    lag = x$lag,
    n_subset = x$n_subset,
    n_shift = x$n_shift,
    default = x$default,
    prefix = x$prefix,
    columns = terms_select(x$terms, info = info),
    skip = x$skip,
    id = x$id
  )
  
}

#' @importFrom dplyr bind_cols mutate
#' @importFrom rlang call2 sym
#' @importFrom tibble as_tibble
#' @importFrom purrr map2 
#' @importFrom recipes bake prep check_name
#' 
#' @export
bake.step_lag_matrix <- function(object, new_data, ...) {

  if (!all(object$lag == as.integer(object$lag)))
    stop("step_lag_matrix requires 'lag' argument to be integer valued.",
         call. = FALSE)

  make_call <- function(col, lag_val, n_subset, n_shift) {
    call2(
      "shift_subset",
      x = sym(col),
      lag = lag_val,
      n_subset = n_subset,
      n_shift = n_shift,
      .ns = "waterlevel"
    )
  }
  
  if(object$n_subset > 1) {
    # subset new_data
    inds <- seq(object$n_shift+1, 
                nrow(new_data) - object$n_shift, 
                object$n_subset)
    lag_mat <- lag_matrix(new_data[[object$columns]],
                          object$lag, 
                          object$n_subset,
                          object$n_shift,
                          object$columns)
    return(bind_cols(new_data[inds,], as_tibble(lag_mat)))
  }

  grid <- expand.grid(col = object$columns, lag_val = object$lag,
                      stringsAsFactors = FALSE)
  calls <- as.vector(map2(grid$col, 
                                 grid$lag_val,
                                 make_call,
                                 n_subset = object$n_subset,
                                 n_shift = object$n_shift))
  newname <- paste0(object$prefix, grid$lag_val, "_", grid$col)
  newname <- gsub('-', 'n', newname)
  calls <- check_name(calls, new_data, object, newname, TRUE)
  
  as_tibble(mutate(new_data, !!!calls))
  
  
}


#' @importFrom recipes printer
print.step_lag_matrix <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Lagging ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }