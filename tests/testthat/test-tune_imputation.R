test_that("imputation works", {
  s <- sim_mat(500, 50, perc_col_NA = 0.5)
  n_feat <- 0.1 * nrow(s$input)

  set.seed(1234)
  output_sink <- withr::local_tempfile()
  sink(output_sink)
  m_knn <- fit_knn(
    X = s$input,
    hp = data.frame(k = 5, maxp_prop = 0.1),
    rng.seed = 1234
  )

  set.seed(1234)
  m1_knn <- impute::impute.knn(s$input, k = 5, maxp = n_feat, rng.seed = 1234)$data
  sink()
  withr::deferred_run()

  m_pca <- fit_imputePCA(X = s$input, hp = data.frame(ncp = 1))
  m1_pca <- t(missMDA::methyl_imputePCA(X = t(s$input), ncp = 1)$completeObs)

  m_methyLImp2 <- fit_methyLImp2(X = s$input)
  m1_methyLImp2 <- t(methyLImp2::mod_methyLImp2_internal(
    t(s$input),
    min = 0,
    max = 1,
    skip_imputation_ids = NULL,
    minibatch_frac = 1,
    minibatch_reps = 1
  ))

  expect_identical(m_knn, m1_knn)
  expect_identical(m_pca, m1_pca)
  expect_identical(m_methyLImp2, m1_methyLImp2)
})

test_that("imputation by group works", {
  s <- sim_mat(
    n = 500,
    50,
    perc_col_NA = 0.1,
    ngrp = 2,
    perc_NA = 0.1
  )
  n_feat <- 0.1 * nrow(s$input)

  set.seed(1234)
  output_sink <- withr::local_tempfile()
  sink(output_sink)
  m_knn <- fit_knn(
    X = s$input,
    group_sample = s$group_sample,
    hp = data.frame(k = 5, maxp_prop = 0.1),
    rng.seed = 1234
  )
  m_knn <- m_knn[row.names(s$input), colnames(s$input)]

  m1_knn_list <- lapply(c("grp1", "grp2"), \(x){
    set.seed(1234)
    impute::impute.knn(
      s$input[, subset(s$group_sample, group == x)$sample_id],
      k = 5,
      maxp = n_feat,
      rng.seed = 1234
    )$data
  })
  m1_knn <- cbind(m1_knn_list[[1]], m1_knn_list[[2]])[row.names(s$input), colnames(s$input)]

  sink()
  withr::deferred_run()

  expect_identical(m_knn, m1_knn)
})
