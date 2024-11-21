inject_na <- function(X, feature_id, sample_id, prop) {
  stopifnot(is.matrix(X))
  # subset the matrix to the group and feature
  M <- X[feature_id, sample_id]
  na_mat <- !is.na(M)
  not_na <- which(na_mat)
  na_size <- floor(length(not_na) * prop)
  if (na_size == 0) {
    stop("Increase the proportion of missing values (prop).")
  }
  # init while loop
  c_miss <- TRUE
  r_miss <- TRUE
  na_loc <- NULL
  # inject NA and ensure no row or column is entirely missing
  while (c_miss || r_miss) {
    na_mat_test <- na_mat
    na_loc <- sample(not_na, size = na_size)
    na_mat_test[na_loc] <- FALSE
    c_miss <- any(colSums(na_mat_test) == 0) # cols with all miss will have sum zero
    r_miss <- any(rowSums(na_mat_test) == 0) # rows with all miss will have sum zero
  }

  # return(list(na_loc = na_loc, na_mat = na_mat)) # for debug if need
  return(na_loc)
}

impute.knn1 <- function(X, feature_id, sample_id, k, maxp_prop, seed, na_loc) {
  stopifnot(is.matrix(X))
  M <- X[feature_id, sample_id]
  truth <- M[na_loc]
  M[na_loc] <- NA
  maxp <- floor(nrow(M) * maxp_prop)
  obj <- impute::impute.knn(data = M, k = k, maxp = maxp, rng.seed = seed)
  estimate <- obj$data[na_loc]
  return(dplyr::tibble(truth = truth, estimate = estimate))
}

tune_prep <- function(X, group_sample, group_feature, rep = 1, hp = NULL, prop) {
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

  stopifnot(
    "prop must be a single numeric value between 0 and 1 (inclusive)" =
      is.numeric(prop) && length(prop) == 1 && 0 < prop && prop <= 1
  )

  df <- tidyr::expand_grid(
    sample_group = unique(group_sample$group),
    feature_group = unique(group_feature$group),
    rep = seq_len(rep)
  )

  if (!is.null(hp)) {
    df$params <- list(hp)
    df <- tidyr::unnest(df, cols = params)
  }

  group_sample_collapsed <- tidyr::nest(
    group_sample,
    sample_id = sample_id,
    .by = "group"
  )
  group_feature_collapsed <- tidyr::nest(
    group_feature,
    feature_id = feature_id,
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
    \(x, y) {
      inject_na(X = X, feature_id = x, sample_id = y, prop = prop)
    }
  )

  return(df)
}

tune_knn <- function(X, group_sample, group_feature, hp, rep = 1, prop = 0.05, parallel = FALSE, progress = FALSE) {
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

  # prep parameters
  params <- tune_prep(X = X, group_sample = group_sample, group_feature = group_feature, rep = rep, hp = hp, prop = prop)

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
      impute.knn1(
        X = X,
        feature_id = feature_id,
        sample_id = sample_id,
        k = k,
        maxp_prop = maxp_prop,
        seed = seed,
        na_loc = na_loc
      )
    },
    .progress = progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  params$method <- "knn"
  return(params)
}
