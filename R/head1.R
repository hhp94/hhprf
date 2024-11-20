#' Print 1:6 Rows and Cols
#'
#' If a matrix/df is too big, Rstudio will try to print out everything and
#' freeze. This is a wrapper around `head(x, c(r, m))`.
#'
#' @param obj data.frame or matrix
#' @param r rows
#' @param m cols
#' @param ... args passed to [utils::head()]
#'
#' @return data.frame or matrix
#' @export
#'
#' @examples
#' head1(mtcars)
head1 <- function(obj, ...) {
  UseMethod("head1")
}

#' @rdname head1
#' @export
head1.default <- function(obj, ...) {
  utils::head(obj, ...)
}

#' @rdname head1
tabular_return <- function(obj, r, m, ...) {
  min_n <- min(r, nrow(obj))
  min_m <- min(m, ncol(obj))
  utils::head(obj, c(min_n, min_m), ...)
}

#' @rdname head1
#' @export
head1.data.frame <- function(obj, r = 6, m = 6, ...) {
  tabular_return(obj = obj, r = r, m = m, ...)
}

#' @rdname head1
#' @export
head1.list <- function(obj, r = 6, m = 6, ...) {
  lapply(obj, head1, ...)
}

#' @rdname head1
#' @export
head1.matrix <- function(obj, r = 6, m = 6, ...) {
  tabular_return(obj = obj, r = r, m = m, ...)
}
