#' Extract Failed Jobs to Create a Retry List
#'
#' @description Parses a job status TSV file to identify tasks that failed (i.e., those
#'   with a non-zero exit status in the second column) and writes their corresponding
#'   commands (from the seventh column) to a new text file for easy resubmission.
#'
#' @param job_status Character string. The path to the job status TSV file. The function
#'   expects the job exit code to be in the second column and the job command to be in
#'   the seventh column.
#' @param out Character string. The destination file path where the failed job commands
#'   will be saved. Defaults to `"joblist_retry.txt"`.
#'
#' @returns A character string containing the path to the output file (`out`). If no
#'   jobs failed, it returns an empty character vector (`character(0)`).
#'
#' @export
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
    return(character(0))
  }

  message(sprintf("Found %d failed jobs. Writing to '%s'...", nrow(failed), out))
  readr::write_lines(failed$X7, out)

  return(out)
}
