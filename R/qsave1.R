#' Modify qs::qsave By Saving to Temp First
#'
#' Saving files to the hdd can create errors. Saving to temp first might create
#' less errors
#'
#' @param obj obj to serialize with qs
#' @param file destination
#' @param overwrite overwrite or not
#' @param ... args to qsave
#'
#' @export
qsave1 <- function(obj, file, overwrite = TRUE, ...) {
  stopifnot(overwrite %in% c(TRUE, FALSE))

  tmp_fn <- fs::file_temp()
  qs2::qs_save(obj, tmp_fn, ...)

  message("Generating hash")
  tmp_file_hash <- tools::md5sum(tmp_fn)

  fs::file_copy(tmp_fn, file, overwrite = overwrite)

  new_file_hash <- tools::md5sum(file)

  # While the file doesn't exists or the hash doesn't match or has not been overwritten
  if (new_file_hash != tmp_file_hash) {
    warning("Warning, hash of new file doesn't match old file")
  }
}
