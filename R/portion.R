#' Extract Finalized Training Set
#'
#' As steps are estimated by `prep`, these operations are
#'  applied to the training set. Rather than running `bake`
#'  to duplicate this processing, this function will return
#'  variables from the processed training set.
#' @inheritParams recipes::bake.recipe
#' @param object A `recipe` object that has been prepared
#'   with the option `retain = TRUE`.
#' @details When preparing a recipe, if the training data set is
#'  retained using `retain = TRUE`, there is no need to `bake` the
#'  recipe to get the preprocessed training set.
#'
#'  `portion` will return the results of a recipes where _all steps_
#'  have been applied to the data, irrespective of the value of
#'  the step's `skip` argument.
#'
#' @export
#' 
#' @importFrom purrr is_empty
#' @importFrom rlang quos
#' @importFrom recipes fully_trained
#' @importFrom tidyselect everything
#' @importFrom dplyr select
#' 
#' 
#' @seealso [recipe()] [prep.recipe()] [bake.recipe()] [juice.recipe()]
portion <- function(object, ...) {
  if (!fully_trained(object))
    stop("At least one step has not been trained. Please ",
         "run `prep`.",
         call. = FALSE)
  
  if(!isTRUE(object$retained))
    stop("Use `retain = TRUE` in `prep` to be able to extract the training set",
         call. = FALSE)
  
  terms <- quos(...)
  if (is_empty(terms))
    terms <- quos(everything())
  keepers <- terms_select(terms = terms, info = object$term_info)
  
  new_data <- object$template[, names(object$template) %in% keepers]
  
  ## Since most models require factors, do the conversion from character
  if (!is.null(object$levels)) {
    var_levels <- object$levels
    var_levels <- var_levels[keepers]
    check_values <-
      vapply(var_levels, function(x)
        (!all(is.na(x))), c(all = TRUE))
    var_levels <- var_levels[check_values]
    # need to replace strings2factors
    # if (length(var_levels) > 0)
    #   new_data <- strings2factors(new_data, var_levels)
  }

  
  convert_tibble_of_matrices(new_data, object$term_info)
  
  
}


convert_tibble_of_matrices <- function(x, term_info) {
  
  role_group <- split(term_info$variable, term_info$role)
  
  as_tibble(lapply(role_group, to_list_of_matrix, x = x))
  
  
}

to_list_of_matrix <- function(col_group, x) {

    return(as.matrix(select(x, col_group)))
  
}




