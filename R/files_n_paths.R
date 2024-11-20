#' Getting Files and Path in a Dir
#'
#' @param dir Directory
#' @param ... Passed to [list.files()]
#'
#' @export
files_n_paths <- function(dir, ...) {
  d <- data.frame(
    files = list.files(dir, ...),
    paths = list.files(dir, ..., full.names = TRUE)
  )

  return(d)
}
