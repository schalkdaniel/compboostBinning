context("Matrix multiplication on binned matrices works")

test_that("same as above", {

  # Simulate Data:
  nsim = 1e6L
  nunique = trunc(sqrt(nsim))

  xunique = runif(n = nunique, min = 0, max = 10)
  k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)

  X = poly(x = xunique, degree = 20L)
  X_full = X[k,]

  xtx_full = t(X_full) %*% X_full

  expect_silent({ xtx = binnedMatMult(X = X, k = k-1, w = 1, FALSE) })
  expect_true(all.equal(xtx_full, xtx, check.attributes = FALSE))
  expect_silent({ xtx_fast = binnedMatMult(X = X, k = k-1, w = 1, TRUE) })
  expect_equal(xtx, xtx_fast)

  w = sample(x = seq_len(20L), size = nsim, replace = TRUE)
  xtwx_full = t(X_full * w)  %*% X_full

  expect_silent({ xtwx = binnedMatMult(X = X, k = k-1, w = w, FALSE) })
  expect_true(all.equal(xtwx_full, xtwx, check.attributes = FALSE))
  expect_silent({ xtwx_fast = binnedMatMult(X = X, k = k-1, w = w, TRUE) })
  expect_equal(xtwx, xtwx_fast)

  y = runif(nsim)

  out_full = t(X_full) %*% y

  expect_silent({ out = binnedMatMultResponse(X = X, y = y, k = k-1, w = 1) })
  expect_true(all.equal(out, out_full, check.attributes = FALSE))

  out_full = t(X_full * w) %*% y

  expect_silent({ out = binnedMatMultResponse(X = X, y = y, k = k-1, w = w) })
  expect_true(all.equal(out, out_full, check.attributes = FALSE))




  x_sp = rep(0, 20 * nunique)
  x_sp[sample(seq_len(20 * nunique), size = 10 * nunique)] = rnorm(10 * nunique)
  X_sp = Matrix::Matrix(x_sp, ncol = 20L, sparse = TRUE)
  X_full = X_sp[k,]

  xtx_full = as.matrix(Matrix::t(X_full) %*% X_full)

  X_sp = Matrix::t(X_sp)

  expect_silent({ xtx = binnedSparseMatMult(X = X_sp, k = k-1, w = 1) })
  expect_true(all.equal(xtx_full, xtx, check.attributes = FALSE))

  w = sample(x = seq_len(20L), size = nsim, replace = TRUE)
  xtwx_full = as.matrix(Matrix::t(X_full * w)  %*% X_full)

  expect_silent({ xtwx = binnedSparseMatMult(X = X_sp, k = k-1, w = w) })
  expect_true(all.equal(xtwx_full, xtwx, check.attributes = FALSE))



  y = runif(nsim)

  out_full = as.matrix(Matrix::t(X_full) %*% y)

  expect_silent({ out = binnedSparseMatMultResponse(X = X_sp, y = y, k = k-1, w = 1) })
  expect_true(all.equal(out, out_full, check.attributes = FALSE))

  out_full = as.matrix(Matrix::t(X_full * w) %*% y)

  expect_silent({ out = binnedSparseMatMultResponse(X = X_sp, y = y, k = k-1, w = w) })
  expect_true(all.equal(out, out_full, check.attributes = FALSE))
})
