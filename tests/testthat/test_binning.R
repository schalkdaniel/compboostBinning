context("Binning works")

test_that("Binning works", {

  # Simulate Data:
  nsim = 1e6L
  nunique = trunc(sqrt(nsim))

  x = rnorm(nsim)

  x_bin = binVectorCustom(x, nunique)
  expect_length(x_bin, nunique)

  x_bin2 = binVector(x)
  expect_equal(x_bin, x_bin2)

  idx = calculateIndexVector(x, x_bin)
  expect_true(all(idx <= nunique))
  expect_length(idx, nsim)
})
