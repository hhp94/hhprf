test_that("is_unique works on vectors", {
  rep_vec <- c(rep(1, 2), 3, 4)
  expect_false(is_unique(rep_vec))
  u_vec <- 1:10
  expect_true(is_unique(u_vec))
})

test_that("is_unique works on data.frame", {
  rep_df <- rbind(mtcars, mtcars)
  expect_equal(nrow(rep_df), nrow(mtcars) * 2)
  expect_false(is_unique(rep_df))
  expect_true(is_unique(mtcars))
})

test_that("is_unique works on matrices", {
  mat_mtcars <- as.matrix(mtcars)
  expect_true(is.matrix(mat_mtcars))
  rep_mat <- rbind(mat_mtcars, mat_mtcars)
  expect_equal(nrow(rep_mat), nrow(mat_mtcars) * 2)
  expect_false(is_unique(rep_mat))
  expect_true(is_unique(mat_mtcars))
})
