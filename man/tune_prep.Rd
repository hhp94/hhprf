% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_imputation.R
\name{tune_prep}
\alias{tune_prep}
\title{Setup Grid Search for Tuning Parameters}
\usage{
tune_prep(
  X,
  group_sample = NULL,
  group_feature = NULL,
  rep = 1,
  hp = NULL,
  prop = 0.05,
  min_na = 50,
  max_na = 5000,
  c_thresh = 0.9,
  r_thresh = 0.9,
  max_iter = 1000,
  transpose = FALSE,
  parallel = FALSE,
  progress = FALSE,
  temp_path = NULL,
  ampute = TRUE
)
}
\arguments{
\item{X}{A numeric matrix. The input data to be imputed. Rows correspond to features, and columns correspond to samples.}

\item{group_sample}{A data frame with two columns: \code{sample_id} and \code{group}. Maps samples to their respective groups. Default is \code{NULL}.}

\item{group_feature}{A data frame with two columns: \code{feature_id} and \code{group}. Maps features to their respective groups. Default is \code{NULL}.}

\item{rep}{A vector of unique elements. The length of the vector is the times the tuning is repeated.}

\item{hp}{A data frame. Optional hyperparameter grid to be used for tuning. Default is \code{NULL}.}

\item{prop}{A numeric value between 0 and 1. Proportion of values to inject as missing in the subset of the matrix. Default is \code{0.05}.}

\item{min_na}{An integer. Minimum number of \code{NA} injected. Default is \code{50}. Used to calculate prediction performance.}

\item{max_na}{An integer. Maximum number of \code{NA} injected. Default is \code{5000}. Used to calculate prediction performance.}

\item{c_thresh}{A numeric value between 0 and 1. Maximum allowed proportion of missing values per column. Default is \code{0.9}.}

\item{r_thresh}{A numeric value between 0 and 1. Maximum allowed proportion of missing values per row. Default is \code{0.9}.}

\item{max_iter}{A positive integer. Maximum number of iterations for injecting missing values. Default is 1000.}

\item{transpose}{Transpose before NA injection or not. Default to FALSE}

\item{parallel}{Logical. If \code{TRUE}, enables parallel processing using the \code{furrr} package. Default is \code{FALSE}.}

\item{progress}{Logical. If \code{TRUE}, displays a progress bar for the tuning process. Default is \code{FALSE}.}

\item{temp_path}{Storage location for matrices. Options:
\itemize{
\item \code{NULL}: Save in temporary directory (cleared at session end)
\item Character string: Save in specified directory path
\item \code{NA}: Store in memory (may significantly increase memory usage)
}}

\item{ampute}{Logical. inject \code{NA} or not (for tuning). Default to \code{TRUE}.}
}
\value{
A data frame
}
\description{
Biological assays often require methods to impute missing data, which may involve hyperparameters.
Additionally, multiple methods may be applied. This function sets up a grid search to tune
imputation methods by injecting \code{NA} multiple times, imputing the data,
and comparing the imputed data to the original data before the \code{NA} injection. The function
also sets up tuning stratified by feature groups (e.g., chromosomes) and sample groups
(e.g., sample batches).
}
