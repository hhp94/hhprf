% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_imputation.R
\name{impute_PCA}
\alias{impute_PCA}
\title{Wrapper for \code{\link[missMDA:methyl_imputePCA]{missMDA::methyl_imputePCA()}}}
\usage{
impute_PCA(X, hash = NULL, ncp, seed, na_loc, ampute, ...)
}
\arguments{
\item{X}{A numeric matrix. The input data to be imputed. Rows correspond to features, and columns correspond to samples.}

\item{hash}{hash of the saved object used to check hash.}

\item{ncp}{number of PCs}

\item{seed}{An integer seed for random number generation.}

\item{na_loc}{A vector of location for amputation. Generated by \code{\link[=inject_na]{inject_na()}}.}

\item{ampute}{Logical. inject \code{NA} or not (for tuning). Default to \code{TRUE}.}

\item{...}{arguments passed to \code{\link[missMDA:methyl_imputePCA]{missMDA::methyl_imputePCA()}}}
}
\description{
Wrapper for \code{\link[missMDA:methyl_imputePCA]{missMDA::methyl_imputePCA()}}
}
