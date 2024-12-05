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
#' @param rep A vector of unique elements. The length of the vector is the times the tuning is repeated.
#' @param hp A data frame. Optional hyperparameter grid to be used for tuning. Default is `NULL`.
#' @param prop A numeric value between 0 and 1. Proportion of values to inject as missing in the subset of the matrix. Default is `0.05`.
#' @param min_na An integer. Minimum number of `NA` injected. Default is `50`. Used to calculate prediction performance.
#' @param max_na An integer. Maximum number of `NA` injected. Default is `5000`. Used to calculate prediction performance.
#' @param c_thresh A numeric value between 0 and 1. Maximum allowed proportion of missing values per column. Default is `0.9`.
#' @param r_thresh A numeric value between 0 and 1. Maximum allowed proportion of missing values per row. Default is `0.9`.
#' @param max_iter A positive integer. Maximum number of iterations for injecting missing values. Default is 1000.
#' @param transpose Transpose before NA injection or not. Default to FALSE
#' @param parallel Logical. If `TRUE`, enables parallel processing using the `furrr` package. Default is `FALSE`.
#' @param progress Logical. If `TRUE`, displays a progress bar for the tuning process. Default is `FALSE`.
#' @param temp_path Storage location for matrices. Options:
#'   * `NULL`: Save in temporary directory (cleared at session end)
#'   * Character string: Save in specified directory path
#'   * `NA`: Store in memory (may significantly increase memory usage)
#' @param ampute Logical. inject `NA` or not (for tuning). Default to `TRUE`.
#'
#' @return A data frame
#' @export
tune_prep <- function(
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
    ampute = TRUE) {
  # input validation
  stopifnot("'X' has to be a matrix" = is.matrix(X))
  if (ncol(X) > nrow(X)) {
    warning("'X' is expected to be long")
  }
  if (anyDuplicated(rep)) {
    stop("'rep' has to contain unique values")
  }
  if (!is.null(hp)) {
    stopifnot(
      "'hp' has to be a data frame" = is.data.frame(hp),
      "'hp' needs at least 1 row" = nrow(hp) > 0
    )
  }
  if (is.null(group_sample)) {
    group_sample <- data.frame(sample_id = colnames(X), group = "all_samples")
  }
  if (is.null(group_feature)) {
    group_feature <- data.frame(feature_id = rownames(X), group = "all_features")
  }
  stopifnot(
    "'group_sample' has to be a data frame with 2 columns: 'sample_id' and 'group'" = is.data.frame(group_sample),
    "'group_sample' needs at least 1 row" = nrow(group_sample) > 0,
    "'group_sample' must contain unique 'sample_id'" = !anyDuplicated(group_sample$sample_id)
  )
  stopifnot(
    "'group_feature' has to be a data frame with 2 columns: 'feature_id' and 'group'" = is.data.frame(group_feature),
    "'group_feature' needs at least 1 row" = nrow(group_feature) > 0,
    "'group_feature' must contain unique 'feature_id'" = !anyDuplicated(group_feature$feature_id)
  )
  stopifnot(
    "'X' has to have row names (feature_id)" = !is.null(row.names(X)),
    "'X' has to have col names (sample_id)" = !is.null(colnames(X)),
    "'X''s row names must match 'group_feature$feature_id'" = all(row.names(X) %in% group_feature$feature_id),
    "'X''s column names must match 'group_sample$sample_id'" = all(colnames(X) %in% group_sample$sample_id)
  )
  stopifnot(
    "'ampute' has to be a logical" = is.logical(ampute) && length(ampute) == 1
  )
  sample_group_counts <- table(group_sample$group)
  if (any(sample_group_counts == 1)) {
    stop("Group sample with only 1 sample detected.")
  }
  feature_group_counts <- table(group_feature$group)
  if (any(feature_group_counts == 1)) {
    stop("Group feature with only 1 feature detected.")
  }
  for (i in c(parallel, progress)) {
    stopifnot(
      "'parallel' and 'progress' must be logical values of length one" =
        is.logical(i) && length(i) == 1
    )
  }
  if (parallel == TRUE && future::nbrOfWorkers() == 1) {
    stop("'parallel' is TRUE, but only one worker is available.")
  }
  # end validation
  ## create combination of the group and features
  df <- tidyr::expand_grid(
    sample_group = unique(group_sample$group),
    feature_group = unique(group_feature$group)
  )
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
  df[["sample_id"]] <- lapply(df[["sample_id"]], \(x){
    x[["sample_id"]]
  })
  df[["feature_id"]] <- lapply(df[["feature_id"]], \(x){
    x[["feature_id"]]
  })
  # save the sub matrices. This is crucial for saving memories to avoid passing the full
  # beta matrix to the workers in parallelization
  df <- save_sub_matrix(df = df, X = X, temp_path = temp_path, progress = progress)
  df[, c("sample_id", "feature_id")] <- NULL

  ## add rep
  df$rep <- list(rep)
  df <- tidyr::unnest(df, cols = dplyr::all_of("rep"))
  df$seed <- sample.int(nrow(df) * 10, size = nrow(df))

  if (ampute) {
    # ampute. Each unique rep gets a unique NA pattern. But same NA pattern is used on
    # different hyper parameters.
    ampute_args <- list(
      prop = list(prop),
      c_thresh = list(c_thresh),
      r_thresh = list(r_thresh),
      min_na = list(min_na),
      max_na = list(max_na),
      max_iter = list(max_iter),
      transpose = list(transpose)
    )
    if (is.null(temp_path) || is.character(temp_path)) {
      ampute_args$X <- df$path
      ampute_args$hash <- df$hash
    } else if (is.na(temp_path)) {
      ampute_args$X <- df$sub_matrix
      ampute_args$hash <- list(NULL)
    } else {
      stop("Unexpected Error")
    }
    lapply_args <- list(.progress = progress)
    if (parallel) {
      lapply_args <- c(lapply_args, list(.options = furrr::furrr_options(seed = TRUE)))
      lapply_fn <- furrr::future_pmap
    } else {
      lapply_fn <- purrr::pmap
    }
    df$na_loc <- do.call(
      lapply_fn,
      c(
        list(
          .l = ampute_args,
          .f = inject_na
        ),
        lapply_args
      )
    )
    df <- tidyr::unnest(df, cols = dplyr::all_of("na_loc"))
  }
  df$ampute <- ampute

  ## add hp
  if (!is.null(hp)) {
    df$params <- list(hp)
    df <- tidyr::unnest(df, cols = dplyr::all_of("params"))
  }

  return(df)
}

#' Save Submatrices to Files or Memory
#'
#' Takes a data frame containing feature and sample IDs, extracts corresponding
#' submatrices from a source matrix, and either saves them to disk or keeps them
#' in memory.
#'
#' @param df data.frame
#' @inheritParams tune_prep
save_sub_matrix <- function(df, X, temp_path, progress) {
  if (!is.data.frame(df) || !all(c("feature_id", "sample_id") %in% names(df))) {
    stop("'df' must be a data frame with 'feature_id' and 'sample_id' columns")
  }
  if (!is.null(temp_path) && !is.na(temp_path)) {
    if (length(temp_path) != 1 || !is.character(temp_path)) {
      stop("'temp_path' must be a single string if not 'NULL' or 'NA'.")
    }
  }
  if (is.null(temp_path)) {
    temp_path <- tempdir()
  }
  if (is.na(temp_path)) {
    df$sub_matrix <- purrr::map2(
      df$feature_id,
      df$sample_id,
      function(feature_id, sample_id) {
        return(X[feature_id, sample_id])
      },
      .progress = progress
    )
    return(df)
  }
  if (!fs::dir_exists(temp_path)) {
    warning("'temp_path' doesn't exists. Creating folder")
    fs::dir_create(temp_path)
  }
  df$path <- fs::fs_path(
    paste(
      temp_path,
      paste(
        df$sample_group,
        df$feature_group,
        paste(rlang::hash(Sys.time()), sample.int(nrow(df) * 10, size = nrow(df)), sep = "_"),
        "qs2",
        sep = "."
      ),
      sep = "/"
    )
  )
  df$hash <- purrr::pmap_chr(
    .l = list(
      feature_id = df$feature_id,
      sample_id = df$sample_id,
      path = df$path
    ),
    .f = function(feature_id, sample_id, path) {
      M <- X[feature_id, sample_id]
      qs2::qs_save(M, path)
      return(rlang::hash(M))
    },
    .progress = progress
  )
  return(df)
}


#' Inject Missing Values
#'
#' Injects missing values `NA` into a subset of a matrix while adhering to column-wise and row-wise
#' missingness thresholds.
#'
#' @inheritParams tune_prep
#' @param hash hash of the saved object used to check hash.
#' @param debug return extra object for debugging.
#'
#' @return A vector of indices indicating the locations of injected `NA` in the subset matrix.
#' @export
inject_na <- function(
    X,
    hash = NULL,
    prop,
    min_na,
    max_na,
    c_thresh = 0.9,
    r_thresh = 0.9,
    max_iter = 1000,
    transpose,
    debug = FALSE) {
  if (is.character(X)) {
    X <- safe_read(X, hash = hash)
  }
  stopifnot(is.logical(transpose), length(transpose) == 1)
  for (i in list(c_thresh, r_thresh)) {
    if (!is.numeric(i)) stop("Thresholds must be numeric values")
    if (length(i) != 1) stop("Thresholds must be single values, not vectors")
    if (i < 0) stop("Thresholds must be non-negative")
    if (i > 1) stop("Thresholds must be less than or equal to 1")
  }

  stopifnot(
    "'max_iter' must be a positive integer" =
      is.numeric(max_iter) && length(max_iter) == 1 && max_iter > 0 && as.integer(max_iter) == max_iter,
    "'prop' must be a single numeric value between 0 and 1 (inclusive)" =
      is.numeric(prop) && length(prop) == 1 && 0 < prop && prop <= 1,
    "'min_na' must be a positive integer" =
      is.numeric(min_na) && length(min_na) == 1 && min_na > 0 && as.integer(min_na) == min_na,
    "'max_na' must be a positive integer" =
      is.numeric(max_na) && length(max_na) == 1 && max_na > 0 && as.integer(max_na) == max_na,
    "'max_na' must be greater than or equal to 'min_na'" =
      max_na >= min_na
  )

  # subset the matrix to the specified features and samples
  if (transpose) {
    X <- t(X)
  }
  na_mat <- !is.na(X)
  not_na <- which(na_mat)
  # calculate `na_size` and make sure falls within [min_na, max_na]
  if (min_na > length(not_na)) {
    stop(
      sprintf(
        "'min_na' (%d) exceeds the number of available non-NA elements (%d).
        Adjust 'min_na' or increase feature/sample group size.",
        min_na, length(not_na)
      )
    )
  }
  na_size <- ceiling(length(not_na) * prop)
  na_size <- max(min(na_size, max_na), min_na)
  if (na_size == 0 || na_size > length(not_na)) {
    stop(
      sprintf(
        "Invalid number of NAs to inject: calculated na_size = %d, available non-NA elements = %d.
        Adjust 'prop', 'min_na', 'max_na', or increase feature/sample group size.",
        na_size, length(not_na)
      )
    )
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
  results <- dplyr::tibble(na_loc = list(na_loc))
  if (debug) {
    results$X <- list(X)
    results$r_names <- list(row.names(X))
    results$c_names <- list(colnames(X))
  }
  return(results)
}

# Tune Functions ----
# model specific validation rules
validation_rules <- list(
  impute_knn = function(hp) {
    stopifnot(
      "'hp' has to be a data.frame with 2 columns 'k' and 'maxp_prop'" =
        is.data.frame(hp) && all(c("k", "maxp_prop") %in% names(hp)),
      "'k' has to be positive integers" =
        all(as.integer(hp[["k"]]) == hp[["k"]]) && min(hp[["k"]]) > 0,
      "'maxp_prop' has to be real numbers between 0 and 1 (inclusive)" =
        is.numeric(hp[["maxp_prop"]]) && min(hp[["maxp_prop"]]) >= 0 &&
          max(hp[["maxp_prop"]]) <= 1,
      "'hp' cannot have duplicated rows" =
        nrow(hp) == nrow(dplyr::distinct(hp))
    )
  },
  impute_PCA = function(hp) {
    stopifnot(
      "'hp' has to be a data.frame with 1 column 'ncp'" =
        is.data.frame(hp) && all(c("ncp") %in% names(hp)),
      "'ncp' has to be positive integers" =
        all(as.integer(hp[["ncp"]]) == hp[["ncp"]]) && min(hp[["ncp"]]) > 0,
      "'hp' cannot have duplicated rows" =
        nrow(hp) == nrow(dplyr::distinct(hp))
    )
  },
  impute_methyLImp2 = function(hp) {}
)

# helper function for common parameters
prepare_params <- function(params, extra_args = list()) {
  main_args <- list(
    seed = params$seed,
    ampute = params$ampute
  )
  if ("na_loc" %in% names(params)) {
    main_args$na_loc <- params$na_loc
  }
  if ("sub_matrix" %in% names(params)) {
    main_args$X <- params$sub_matrix
  } else {
    main_args$X <- params$path
    main_args$hash <- params$hash
  }
  return(c(main_args, extra_args))
}

# model specific parameters
argv_list <- list(
  impute_knn = function(params) {
    prepare_params(params, extra_args = list(k = params$k, maxp_prop = params$maxp_prop))
  },
  impute_PCA = function(params) {
    prepare_params(params, extra_args = list(ncp = params$ncp))
  },
  impute_methyLImp2 = function(params) {
    prepare_params(params)
  }
)

# model configuration
model_config <- list(
  impute_knn = list(
    fun_name = "impute_knn",
    transpose = FALSE
  ),
  impute_PCA = list(
    fun_name = "impute_PCA",
    transpose = TRUE
  ),
  impute_methyLImp2 = list(
    fun_name = "impute_methyLImp2",
    transpose = TRUE
  )
)

# factory function for creating tuning functions
create_impute_function <- function(method, ampute) {
  function(X,
           group_sample = NULL,
           group_feature = NULL,
           hp = NULL,
           rep = 1,
           prop = 0.05,
           min_na,
           max_na,
           c_thresh = 0.9,
           r_thresh = 0.9,
           max_iter = 1000,
           parallel = FALSE,
           progress = FALSE,
           temp_path = NULL,
           ...) {
    ampute <- force(ampute)
    # set up parallel processing function
    lapply_fn <- if (parallel) {
      furrr::future_pmap
    } else {
      purrr::pmap
    }
    # use withr
    if (is.null(temp_path)) {
      temp_path <- withr::local_tempdir()
    }
    # validate hyperparameters data.frame
    validation_rules[[method]](hp)
    if (!ampute && (!is.null(hp) && nrow(hp) > 1)) {
      stop("For imputation, 'hp' can only have 1 row")
    }
    # capture extra args passed to imputation functions
    extra_args <- list(...)
    # prep parameters
    prep_args <- list(
      X = X,
      group_sample = group_sample,
      group_feature = group_feature,
      hp = hp,
      transpose = model_config[[method]]$transpose,
      parallel = parallel,
      progress = progress,
      temp_path = temp_path
    )
    if (ampute) {
      prep_args$rep <- rep
      prep_args$prop <- prop
      prep_args$c_thresh <- c_thresh
      prep_args$r_thresh <- r_thresh
      prep_args$max_iter <- max_iter
      prep_args$min_na <- min_na
      prep_args$max_na <- max_na
      prep_args$ampute <- TRUE
    } else {
      prep_args$ampute <- FALSE
      prep_args$rep <- 1
    }
    params <- do.call(tune_prep, prep_args)
    # delete temp files on exit
    if ("path" %in% names(params)) {
      initial_paths <- params$path
      on.exit(
        {
          existing_files <- initial_paths[file.exists(initial_paths)]
          unlink(existing_files)
        },
        add = TRUE
      )
    }
    # prepare model-specific parameters
    impute_args <- argv_list[[method]](params)
    lapply_args <- list(.progress = progress)
    if (parallel) {
      lapply_args <- c(lapply_args, list(.options = furrr::furrr_options(seed = TRUE)))
    }
    # impute
    params$tune <- do.call(
      "lapply_fn",
      c(
        list(
          .l = impute_args,
          .f = function(...) {
            do.call(model_config[[method]]$fun_name, c(list(...), extra_args))
          }
        ),
        lapply_args
      )
    )
    params$method <- method
    if (ampute) {
      return(clean_tune(params))
    } else {
      return(clean_impute(params))
    }
  }
}

clean_tune <- function(params) {
  if ("na_loc" %in% names(params)) {
    params$na_injected <- sapply(params$na_loc, length)
  }
  exclude <- c("sample_id", "feature_id", "na_loc", "r_names", "c_names", "path", "hash")
  return(params[, setdiff(names(params), exclude)])
}

clean_impute <- function(params) {
  if (nrow(params) == 1) {
    return(params$tune[[1]])
  }
  results <- dplyr::summarize(
    params,
    tune = list(purrr::reduce(tune, rbind)),
    .by = "sample_group"
  )
  return(purrr::reduce(results$tune, cbind))
}

# tune_knn ----
# Create specific tuning functions
#' Tune [impute::impute.knn()]
#'
#' Tunes the parameters of the [impute::impute.knn()] method for imputing missing values.
#' This function evaluates the performance of different parameter combinations through repeated
#' `NA` injection and imputation, across various feature and sample groupings.
#'
#' @inheritParams tune_prep
#' @param ... Arguments passed to [impute::impute.knn()]
#' @return A data frame containing the tuning parameters and a column called `tune`
#' that contains the tuning results for each parameter combination.`tune` can be used
#' for performance measurement (i.e. with [yardstick::rmse()]).
#'
#' @export
tune_knn <- create_impute_function("impute_knn", TRUE)

# tune_imputePCA ----
#' Tune [missMDA::imputePCA()]
#'
#' Tunes the parameters of the [missMDA::imputePCA()] method for imputing missing values.
#' This function evaluates the performance of different parameter combinations through repeated
#' `NA` injection and imputation, across various feature and sample groupings.
#'
#' @inheritParams tune_prep
#' @param ... Arguments passed to [missMDA::imputePCA()]
#'
#' @inherit tune_knn return
#'
#' @export
tune_imputePCA <- create_impute_function("impute_PCA", TRUE)

# tune_methyLImp2 ----
#' Tune [methyLImp2::methyLImp2()]
#'
#' Tunes the parameters of the [methyLImp2::methyLImp2()] method for imputing missing values.
#' This function evaluates the performance of different parameter combinations through repeated
#' `NA` injection and imputation, across various feature and sample groupings.
#'
#' @inheritParams tune_prep
#' @param ... Arguments passed to [methyLImp2::methyLImp2()]
#'
#' @inherit tune_knn return
#'
#' @export
tune_methyLImp2 <- create_impute_function("impute_methyLImp2", TRUE)

#' Wrapper for [impute::impute.knn()]
#'
#' @inheritParams inject_na
#' @inheritParams tune_prep
#' @param k Integer specifying the number of nearest neighbors to use.
#' @param maxp_prop A numeric value specifying the proportion of features used per block.
#' @param na_loc A vector of location for amputation. Generated by [inject_na()].
#' @param seed An integer seed for random number generation.
#' @param ... arguments passed to [impute::impute.knn()]
impute_knn <- function(
    X,
    hash = NULL,
    k,
    maxp_prop,
    seed,
    na_loc,
    ampute,
    ...) {
  set.seed(seed)
  if (is.character(X)) {
    X <- safe_read(X, hash = hash)
  }
  stopifnot(is.matrix(X))
  if (ampute) {
    truth <- X[na_loc]
    X[na_loc] <- NA
  }
  maxp <- ceiling(nrow(X) * maxp_prop)
  # 'tryCatch' object. Failure to fit will result in NA
  obj <- tryCatch(
    impute::impute.knn(
      data = X,
      k = k,
      maxp = maxp,
      ...
    ),
    error = error_fn
  )
  if (ampute) {
    if (is.null(obj)) {
      return(dplyr::tibble(truth = NA, estimate = NA))
    }
    estimate <- obj$data[na_loc]
    return(dplyr::tibble(truth = truth, estimate = estimate))
  } else {
    return(obj$data)
  }
}

#' Wrapper for [missMDA::methyl_imputePCA()]
#'
#' @inheritParams inject_na
#' @inheritParams impute_knn
#' @param ncp number of PCs
#' @param ... arguments passed to [missMDA::methyl_imputePCA()]
impute_PCA <- function(
    X,
    hash = NULL,
    ncp,
    seed,
    na_loc,
    ampute,
    ...) {
  set.seed(seed)
  if (is.character(X)) {
    X <- safe_read(X, hash = hash)
  }
  stopifnot(is.matrix(X))
  # imputePCA have samples in the rows so you have to T
  X <- t(X)
  if (ampute) {
    truth <- X[na_loc]
    X[na_loc] <- NA
  }
  obj <- tryCatch(
    missMDA::methyl_imputePCA(X, ncp = ncp, ...),
    error = error_fn
  )
  if (ampute) {
    if (is.null(obj)) {
      return(dplyr::tibble(truth = NA, estimate = NA))
    }
    estimate <- obj$completeObs[na_loc]
    return(dplyr::tibble(truth = truth, estimate = estimate))
  } else {
    return(t(obj$completeObs))
  }
}

#' Wrapper for [methyLImp2::methyLImp2()]
#'
#' @inheritParams inject_na
#' @inheritParams impute_knn
impute_methyLImp2 <- function(
    X,
    hash = NULL,
    seed,
    na_loc,
    ampute) {
  set.seed(seed)
  if (is.character(X)) {
    X <- safe_read(X, hash = hash)
  }
  stopifnot(is.matrix(X))
  # `methyLImp2::mod_methyLImp2_internal` have samples in the rows so you have to t()
  X <- t(X)
  if (ampute) {
    truth <- X[na_loc]
    X[na_loc] <- NA
  }
  obj <- tryCatch(
    methyLImp2::mod_methyLImp2_internal(
      dat = X,
      min = 0,
      max = 1,
      skip_imputation_ids = NULL,
      minibatch_frac = 1,
      minibatch_reps = 1,
    ),
    error = error_fn
  )
  if (ampute) {
    if (is.null(obj) || is.character(obj)) {
      return(dplyr::tibble(truth = NA, estimate = NA))
    }
    estimate <- obj[na_loc]
    return(dplyr::tibble(truth = truth, estimate = estimate))
  } else {
    return(t(obj))
  }
}

# final imputation ----
# fit_knn ----
#' Fit Tuned KNN Imputation
#'
#' Fits a K-Nearest Neighbors (KNN) imputation model after tuning with [tune_knn()].
#' This function imputes missing values in the beta matrix.
#'
#' @inheritParams tune_prep
#' @param ... Additional arguments passed to [impute::impute.knn()]
#' @return A beta matrix with imputed values.
#'
#' @export
fit_knn <- create_impute_function("impute_knn", FALSE)

# fit_imputePCA ----
#' Fit Tuned imputePCA Imputation
#'
#' Fits a [missMDA::imputePCA()] imputation model after tuning with [tune_imputePCA()].
#' This function imputes missing values in the beta matrix.
#'
#' @inheritParams tune_prep
#' @param ... Additional arguments passed to [missMDA::imputePCA()]
#' @return A beta matrix with imputed values.
#'
#' @export
fit_imputePCA <- create_impute_function("impute_PCA", FALSE)

# fit_methyLImp2 ----
#' Fit methyLImp2 Imputation
#'
#' Fits a [methyLImp2::methyLImp2()] imputation model.
#' This function imputes missing values in the beta matrix.
#'
#' @inheritParams tune_prep
#' @param ... Additional arguments passed to [methyLImp2::methyLImp2()]
#' @return A beta matrix with imputed values.
#'
#' @export
fit_methyLImp2 <- create_impute_function("impute_methyLImp2", FALSE)

# other ----
utils::globalVariables(c("sample_id", "feature_id", "tune"))
error_fn <- function(e) {
  message("Error: ", e$message)
  return(NULL)
}

safe_read <- function(obj, hash = NULL) {
  obj <- qs2::qs_read(obj, validate_checksum = TRUE)
  if (!is.null(hash)) {
    stopifnot("'hash' mismatch" = rlang::hash(obj) == hash)
  }
  return(obj)
}
