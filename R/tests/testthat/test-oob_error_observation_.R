test_that("oob_error_observation_ returns per-observation OOB squared errors", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  set.seed(123)
  n <- 60
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y  <- 1 + 2 * X1 - 1.5 * X2 + rnorm(n, sd = 0.5)

  df <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y)

  model <- ranger::ranger(
    Y ~ .,
    data = df,
    num.trees = 30,
    keep.inbag = TRUE
  )

  out <- oob_error_observation_(model, df[, c("X1", "X2", "X3")], response = df$Y)

  expect_true(is.list(out))
  expect_true(all(c("oob_error_obs", "list_inbag") %in% names(out)))

  expect_length(out$oob_error_obs, n)
  expect_true(is.numeric(out$oob_error_obs))
  expect_true(all(out$oob_error_obs[is.finite(out$oob_error_obs)] >= 0))

  expect_true(is.list(out$list_inbag))
  expect_length(out$list_inbag, model$num.trees)
})

test_that("oob_error_observation_ errors if keep.inbag was not used", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  set.seed(1)
  n <- 40
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  df$Y <- df$X1 + rnorm(n)

  model_no_inbag <- ranger::ranger(
    Y ~ .,
    data = df,
    num.trees = 10,
    keep.inbag = FALSE
  )

  expect_error(
    oob_error_observation_(model_no_inbag, df[, c("X1", "X2")], response = df$Y),
    "keep.inbag"
  )
})