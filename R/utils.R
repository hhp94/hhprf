retry_dsq <- function(job_status, out = "joblist_retry.txt") {
  if (!file.exists(job_status)) {
    stop("Error: The file '", job_status, "' does not exist.")
  }

  job_file <- readr::read_tsv(
    job_status,
    col_names = FALSE,
    col_select = c(2, 7)
  )

  failed <- job_file[job_file$X2 != 0, ]

  if (nrow(failed) == 0) {
    message("No jobs failed.")
    return(invisible(character(0)))
  }

  message(sprintf("Found %d failed jobs. Writing to '%s'...", nrow(failed), out))
  readr::write_lines(failed$X7, out)

  invisible(failed$X7)
}
