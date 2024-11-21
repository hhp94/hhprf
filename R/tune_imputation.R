#' Setup Grid Search for Tuning Parameters
#'
#' Biological assays often require methods to impute missing data, which may involve hyperparameters.
#' Additionally, multiple methods may be applied. This function sets up a grid search to tune
#' imputation methods by injecting `NA` multiple times, imputing the data,
#' and comparing the imputed data to the original data before the `NA` injection. The function
#' also sets up tuning stratified by feature groups (e.g., chromosomes) and sample groups
#' (e.g., sample batches).
#'
#' @param X A numeric matrix. The input data to be imputed. Rows correspond to features, and columns correspond to samples.
#' @param group_sample A data frame with two columns: `sample_id` and `group`. Maps samples to their respective groups. Default is `NULL`.
#' @param group_feature A data frame with two columns: `feature_id` and `group`. Maps features to their respective groups. Default is `NULL`.
#' @param rep A positive integer. The number of repeats for NA injection. Default is `1`.
#' @param hp A data frame. Optional hyperparameter grid to be used for tuning. Default is `NULL`.
#' @param prop A numeric value between 0 and 1. Proportion of values to inject as missing in the subset of the matrix. Default is `0.05`.
#' @param c_thresh A numeric value between 0 and 1. Maximum allowed proportion of missing values per column. Default is `0.9`.
#' @param r_thresh A numeric value between 0 and 1. Maximum allowed proportion of missing values per row. Default is `0.9`.
#' @param max_iter A positive integer. Maximum number of iterations for injecting missing values. Default is 1000.
#' @param transpose Transpose before NA injection or not. Default to FALSE
#'
#' @return A data frame
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' row.names(X) <- 1:10
#' colnames(X) <- 1:10
#' group_sample <- data.frame(sample_id = colnames(X), group = rep(1:2, each = 5))
#' group_feature <- data.frame(feature_id = rownames(X), group = rep(1:2, each = 5))
#' tune_prep(X, group_sample, group_feature, rep = 3, prop = 0.2, c_thresh = 1, r_thresh = 1)
#' @export
tune_prep <- function(
    X,
    group_sample = NULL,
    group_feature = NULL,
    rep = 1,
    hp = NULL,
    prop = 0.05,
    c_thresh = 0.9,
    r_thresh = 0.9,
    max_iter = 1000,
    transpose = FALSE) {
  # input validation
  stopifnot("X has to be a matrix" = is.matrix(X))

  if (ncol(X) > nrow(X)) {
    warning("X is expected to be long")
  }
  if (!(length(rep) == 1 && is.numeric(rep) && rep == as.integer(rep) && rep >= 1)) {
    stop("`rep` has to be a positive integer")
  }
  if (!is.null(hp)) {
    stopifnot(
      "`hp` has to be a data frame" = is.data.frame(hp),
      "`hp` needs at least 1 row" = nrow(hp) > 0
    )
  }

  if (is.null(group_sample)) {
    group_sample <- data.frame(sample_id = colnames(X), group = "all_samples")
  }

  if (is.null(group_feature)) {
    group_feature <- data.frame(feature_id = rownames(X), group = "all_features")
  }

  stopifnot(
    "`group_sample` has to be a data frame with 2 columns: `sample_id` and `group`" = is.data.frame(group_sample),
    "`group_sample` needs at least 1 row" = nrow(group_sample) > 0,
    "`group_sample` must contain unique `sample_id`" = !anyDuplicated(group_sample$sample_id)
  )

  stopifnot(
    "`group_feature` has to be a data frame with 2 columns: `feature_id` and `group`" = is.data.frame(group_feature),
    "`group_feature` needs at least 1 row" = nrow(group_feature) > 0,
    "`group_feature` must contain unique `feature_id`" = !anyDuplicated(group_feature$feature_id)
  )

  stopifnot(
    "`X` has to have row names (feature_id)" = !is.null(row.names(X)),
    "`X` has to have col names (sample_id)" = !is.null(colnames(X)),
    "`X`'s row names must match `group_feature$feature_id`" = all(row.names(X) %in% group_feature$feature_id),
    "`X`'s column names must match `group_sample$sample_id`" = all(colnames(X) %in% group_sample$sample_id)
  )

  sample_group_counts <- table(group_sample$group)
  if (any(sample_group_counts == 1)) {
    stop("Each group sample must have more than 1 sample.")
  }

  feature_group_counts <- table(group_feature$group)
  if (any(feature_group_counts == 1)) {
    stop("Each group feature must have more than 1 feature.")
  }

  df <- tidyr::expand_grid(
    sample_group = unique(group_sample$group),
    feature_group = unique(group_feature$group),
    rep = seq_len(rep)
  )

  if (!is.null(hp)) {
    df$params <- list(hp)
    df <- tidyr::unnest(df, cols = dplyr::all_of("params"))
  }

  group_sample_collapsed <- tidyr::nest(
    group_sample,
    sample_id = dplyr::all_of("sample_id"),
    .by = "group"
  )
  group_feature_collapsed <- tidyr::nest(
    group_feature,
    feature_id = dplyr::all_of("feature_id"),
    .by = "group"
  )

  df <- dplyr::inner_join(df, group_sample_collapsed, by = c("sample_group" = "group"))
  df <- dplyr::inner_join(df, group_feature_collapsed, by = c("feature_group" = "group"))
  df <- dplyr::mutate(
    df,
    sample_id = lapply(sample_id, \(x){
      x$sample_id
    }),
    feature_id = lapply(feature_id, \(x){
      x$feature_id
    }),
    seed = sample.int(10000, size = dplyr::n())
  )

  # ampute
  df$na_loc <- purrr::map2(
    df$feature_id,
    df$sample_id,
    \(fid, sid) {
      inject_na(
        X = X,
        feature_id = fid,
        sample_id = sid,
        prop = prop,
        c_thresh = c_thresh,
        r_thresh = r_thresh,
        max_iter = max_iter,
        transpose = transpose
      )
    }
  )
  df <- tidyr::unnest(df, cols = dplyr::all_of("na_loc"))
  return(df)
}

#' Inject Missing Values
#'
#' Injects missing values `NA` into a subset of a matrix while adhering to column-wise and row-wise
#' missingness thresholds.
#'
#' @inheritParams tune_prep
#' @param feature_id A vector of feature IDs to subset rows of `X`.
#' @param sample_id A vector of sample IDs to subset columns of `X`.
#' @param debug return extra object for debugging.
#'
#' @return A vector of indices indicating the locations of injected `NA` in the subset matrix.
#'
#' @examples
#' X <- matrix(1:1000, nrow = 100, ncol = 10)
#' colnames(X) <- 1:10
#' rownames(X) <- 1:100
#' feature_id <- rownames(X)[1:5]
#' sample_id <- colnames(X)[1:5]
#' inject_na(X, feature_id, sample_id, prop = 0.1, c_thresh = 1, r_thresh = 1, transpose = FALSE)
#' @export
inject_na <- function(
    X,
    feature_id,
    sample_id,
    prop,
    c_thresh = 0.9,
    r_thresh = 0.9,
    max_iter = 1000,
    transpose,
    debug = FALSE) {
  stopifnot(is.matrix(X))
  stopifnot(is.logical(transpose), length(transpose) == 1)
  for (i in list(c_thresh, r_thresh)) {
    if (!is.numeric(i)) stop("Thresholds must be numeric values")
    if (length(i) != 1) stop("Thresholds must be single values, not vectors")
    if (i < 0) stop("Thresholds must be non-negative")
    if (i > 1) stop("Thresholds must be less than or equal to 1")
  }
  stopifnot("max_iter must be positive" = (is.numeric(max_iter) && length(max_iter) == 1 && max_iter > 0 && as.integer(max_iter) == max_iter))
  stopifnot(
    "prop must be a single numeric value between 0 and 1 (inclusive)" =
      is.numeric(prop) && length(prop) == 1 && 0 < prop && prop <= 1
  )
  # subset the matrix to the specified features and samples
  M <- if (transpose) {
    t(X[feature_id, sample_id])
  } else {
    X[feature_id, sample_id]
  }
  na_mat <- !is.na(M)
  not_na <- which(na_mat)
  na_size <- ceiling(length(not_na) * prop)
  if (na_size == 0) {
    stop("Increase the proportion of missing values (prop).")
  }
  # max allowed missing counts per column and row
  max_col_miss <- floor(nrow(na_mat) * c_thresh)
  max_row_miss <- floor(ncol(na_mat) * r_thresh)
  # initialize variables for the while loop
  c_miss <- TRUE
  r_miss <- TRUE
  na_loc <- NULL
  iter <- 0
  # inject NAs while ensuring missingness thresholds and iter are not exceeded
  while (c_miss || r_miss) {
    iter <- iter + 1
    if (iter > max_iter) {
      stop("NA injection failed. Lower 'prop' or increase 'c_thresh' and 'r_thresh'.")
    }
    na_mat_test <- na_mat
    na_loc <- sample(not_na, size = na_size)
    na_mat_test[na_loc] <- FALSE
    # calculate the counts of missing values in columns and rows
    col_miss_count <- nrow(na_mat_test) - colSums(na_mat_test)
    row_miss_count <- ncol(na_mat_test) - rowSums(na_mat_test)
    # check if any column or row exceeds the missingness thresholds
    c_miss <- any(col_miss_count > max_col_miss)
    r_miss <- any(row_miss_count > max_row_miss)
  }
  if (debug) {
    return(list(na_loc = na_loc, M = M))
  }
  return(
    dplyr::tibble(
      na_loc = list(na_loc),
      r_names = list(row.names(M)),
      c_names = list(colnames(M))
    )
  )
}

# impute::impute.knn ----
#' Tune [impute::impute.knn()]
#'
#' Tunes the parameters of the [impute::impute.knn()] method for imputing missing values.
#' This function evaluates the performance of different parameter combinations through repeated
#' NA injection and imputation, across various feature and sample groupings.
#'
#' @inheritParams tune_prep
#' @param parallel Logical. If `TRUE`, enables parallel processing using the `furrr` package. Default is `FALSE`.
#' @param progress Logical. If `TRUE`, displays a progress bar for the tuning process. Default is `FALSE`.
#' @param ... Arguments passed to [impute::impute.knn()]
#' @return A data frame containing the tuning parameters and a column called `tune`
#' that contains the tuning results for each parameter combination.`tune` can be used
#' for performance measurement (i.e. with [yardstick::rmse()]).
#'
#' @examples
#' X <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' row.names(X) <- 1:100
#' colnames(X) <- 1:10
#' group_sample <- data.frame(sample_id = colnames(X), group = rep(1:2, length.out = 10))
#' group_feature <- data.frame(feature_id = rownames(X), group = rep(1:2, length.out = 100))
#' hp <- data.frame(k = c(3, 5), maxp_prop = c(0.1, 0.2))
#' tune_knn(
#'   X,
#'   group_sample,
#'   group_feature,
#'   hp,
#'   rep = 2,
#'   parallel = FALSE,
#'   progress = TRUE,
#'   c_thresh = 1,
#'   r_thresh = 1
#' )
#' @export
tune_knn <- function(
    X,
    group_sample,
    group_feature,
    hp,
    rep = 1,
    prop = 0.05,
    c_thresh = 0.9,
    r_thresh = 0.9,
    max_iter = 1000,
    parallel = FALSE,
    progress = FALSE,
    ...) {
  if (parallel) {
    lapply_fn <- furrr::future_pmap
  } else {
    lapply_fn <- purrr::pmap
  }

  stopifnot(
    "`hp` has to be a data.frame with 2 columns `k` and `maxp_prop`" = is.data.frame(hp) && all(c("k", "maxp_prop") %in% names(hp)),
    "`k` has to be positive integers" = all(as.integer(hp[["k"]]) == hp[["k"]]) &&
      min(hp[["k"]]) > 0,
    "`maxp_prop` has to be real numbers between 0 and 1 (inclusive)" = is.numeric(hp[["maxp_prop"]]) &&
      min(hp[["maxp_prop"]]) >= 0 && max(hp[["maxp_prop"]]) <= 1,
    "`hp` cannot have duplicated rows" = nrow(hp) == nrow(dplyr::distinct(hp))
  )

  extra_args <- list(...)
  # prep parameters
  params <- tune_prep(
    X = X,
    group_sample = group_sample,
    group_feature = group_feature,
    rep = rep,
    hp = hp,
    prop = prop,
    c_thresh = c_thresh,
    r_thresh = r_thresh,
    max_iter = max_iter,
    transpose = FALSE
  )

  # impute
  params$tune <- lapply_fn(
    .l = list(
      X = list(X),
      sample_id = params$sample_id,
      feature_id = params$feature_id,
      k = params$k,
      maxp_prop = params$maxp_prop,
      seed = params$seed,
      na_loc = params$na_loc
    ),
    .f = function(X,
                  feature_id,
                  sample_id,
                  k,
                  maxp_prop,
                  seed,
                  na_loc,
                  ...) {
      do.call(
        impute.knn1,
        c(
          list(
            X = X,
            feature_id = feature_id,
            sample_id = sample_id,
            k = k,
            maxp_prop = maxp_prop,
            seed = seed,
            na_loc = na_loc
          ),
          extra_args
        )
      )
    },
    .progress = progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  params$method <- "knn"
  return(params)
}

#' Wrapper for [impute::impute.knn()]
#'
#' @inheritParams inject_na
#' @param k Integer specifying the number of nearest neighbors to use.
#' @param maxp_prop A numeric value specifying the proportion of features used per block.
#' @param na_loc A vector of location for amputation. Generated by [inject_na()].
#' @param seed An integer seed for random number generation.
#' @param ... arguments passed to [impute::impute.knn()]
impute.knn1 <- function(
    X,
    feature_id,
    sample_id,
    k,
    maxp_prop,
    seed,
    na_loc,
    ...) {
  set.seed(seed)
  stopifnot(is.matrix(X))
  M <- X[feature_id, sample_id]
  truth <- M[na_loc]
  M[na_loc] <- NA
  maxp <- ceiling(nrow(M) * maxp_prop)
  obj <- impute::impute.knn(data = M, k = k, maxp = maxp, rng.seed = seed, ...)
  estimate <- obj$data[na_loc]
  return(dplyr::tibble(truth = truth, estimate = estimate))
}

# missMDA::methyl_imputePCA ----
#' Tune [missMDA::methyl_imputePCA()]
#'
#' Tunes the parameters of the [missMDA::methyl_imputePCA()] method for imputing missing values.
#' This function evaluates the performance of different parameter combinations through repeated
#' NA injection and imputation, across various feature and sample groupings.
#'
#' @inheritParams tune_prep
#' @param parallel Logical. If `TRUE`, enables parallel processing using the `furrr` package. Default is `FALSE`.
#' @param progress Logical. If `TRUE`, displays a progress bar for the tuning process. Default is `FALSE`.
#' @param ... Arguments passed to [missMDA::methyl_imputePCA()]
#' @return A data frame containing the tuning parameters and a column called `tune`
#' that contains the tuning results for each parameter combination.`tune` can be used
#' for performance measurement (i.e. with [yardstick::rmse()]).
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' row.names(X) <- 1:10
#' colnames(X) <- 1:10
#' group_sample <- data.frame(sample_id = colnames(X), group = rep(1:2, each = 5))
#' group_feature <- data.frame(feature_id = rownames(X), group = rep(1:2, each = 5))
#' hp <- data.frame(ncp = 1)
#' tune_imputePCA(
#'   X,
#'   group_sample,
#'   group_feature,
#'   hp,
#'   rep = 2,
#'   parallel = FALSE,
#'   progress = TRUE,
#'   c_thresh = 1,
#'   r_thresh = 1
#' )
#' @export
tune_imputePCA <- function(
    X,
    group_sample,
    group_feature,
    hp,
    rep = 1,
    prop = 0.05,
    c_thresh = 0.9,
    r_thresh = 0.9,
    max_iter = 1000,
    parallel = FALSE,
    progress = FALSE,
    ...) {
  if (parallel) {
    lapply_fn <- furrr::future_pmap
  } else {
    lapply_fn <- purrr::pmap
  }

  stopifnot(
    "`hp` has to be a data.frame with 1 column `ncp`" = is.data.frame(hp) && all(c("ncp") %in% names(hp)),
    "`ncp` has to be positive integers" = all(as.integer(hp[["ncp"]]) == hp[["ncp"]]) &&
      min(hp[["ncp"]]) > 0,
    "`hp` cannot have duplicated rows" = nrow(hp) == nrow(dplyr::distinct(hp))
  )

  extra_args <- list(...)
  # prep parameters
  params <- tune_prep(
    X = X,
    group_sample = group_sample,
    group_feature = group_feature,
    rep = rep,
    hp = hp,
    prop = prop,
    c_thresh = c_thresh,
    r_thresh = r_thresh,
    max_iter = max_iter,
    transpose = TRUE
  )

  # impute
  params$tune <- lapply_fn(
    .l = list(
      X = list(X),
      sample_id = params$sample_id,
      feature_id = params$feature_id,
      ncp = params$ncp,
      seed = params$seed,
      na_loc = params$na_loc
    ),
    .f = function(X,
                  feature_id,
                  sample_id,
                  ncp,
                  seed,
                  na_loc,
                  ...) {
      do.call(
        "imputePCA1",
        c(
          list(
            X = X,
            feature_id = feature_id,
            sample_id = sample_id,
            ncp = ncp,
            seed = seed,
            na_loc = na_loc
          ),
          extra_args
        )
      )
    },
    .progress = progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  params$method <- "imputePCA"
  return(params)
}

#' Wrapper for [missMDA::methyl_imputePCA()]
#'
#' @inheritParams inject_na
#' @param ncp number of PCs.
#' @param na_loc A vector of location for amputation. Generated by [inject_na()]
#' @param seed An integer seed for random number generation.
#' @param ... arguments passed to [missMDA::methyl_imputePCA()]
imputePCA1 <- function(
    X,
    feature_id,
    sample_id,
    ncp,
    seed,
    na_loc,
    ...) {
  set.seed(seed)
  stopifnot(is.matrix(X))
  # imputePCA have samples in the rows so you have to T
  M <- t(X[feature_id, sample_id])
  truth <- M[na_loc]
  M[na_loc] <- NA
  obj <- missMDA::methyl_imputePCA(M, ncp = ncp, ...)
  estimate <- obj$completeObs[na_loc]
  return(dplyr::tibble(truth = truth, estimate = estimate))
}

# other ----
utils::globalVariables(c("sample_id", "feature_id"))
