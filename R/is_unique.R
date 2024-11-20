#' Check if a data.frame or vector contains unique values only
#'
#' @param obj an object
#'
#' @return boolean of yes unique or no not unique
#' @export
is_unique <- function(obj) {
  UseMethod("is_unique")
}

#' @rdname is_unique
#' @export
is_unique.default <- function(obj) {
  if (length(obj) == length(unique(obj))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
