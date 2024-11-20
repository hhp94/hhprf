test_that("head1 works", {
  m <- matrix(rep(NA, times = 10^2), ncol = 10)
  h <- head1(m)
  expect_equal(nrow(h), 6)
  expect_equal(ncol(h), 6)
})
