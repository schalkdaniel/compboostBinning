context("Demmler-Reinsch-Orthogonalization works")

test_that("Demmler-Reinsch-Orthogonalization can be computed", {

  # Simulate Data:
  nsim = 1e6L
  nunique = trunc(sqrt(nsim))

  xunique = runif(n = nunique, min = 0, max = 10)
  k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)

  X = poly(x = xunique, degree = 20L)
  X_full = X[k,]

  xtx_full = t(X_full) %*% X_full
  xtx = binnedMatMult(X = X, k = k-1, w = w)

  expect_true(all.equal(xtx_full, xtx, check.attributes = FALSE))

  w = sample(x = seq_len(20L), size = nsim, replace = TRUE)

  xtwx_full = t(X_full * w)  %*% X_full
  xtwx = binnedMatMult(X = X, k = k-1, w = w)

  expect_equal(all.equal(xtwx_full, xtwx, check.attributes = FALSE))

  y = runif(nsim)

  out = binnedMatMultResponse(X = X, y = y, k = k-1, w = 1)
  out_full = t(X_full) %*% y

  expect_true(all.equal(out, out_full, check.attributes = FALSE))

  out = binnedMatMultResponse(X = X, y = y, k = k-1, w = w)
  out_full = t(X_full * w) %*% y

  expect_true(all.equal(out, out_full, check.attributes = FALSE))



})
