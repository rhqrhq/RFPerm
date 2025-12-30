test_that("permrf returns a valid p-value in [0,1]", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  set.seed(101)
  n_tr <- 90
  n_te <- 70
  p <- 3

  df_train <- as.data.frame(replicate(p, rnorm(n_tr)))
  names(df_train) <- paste0("X", 1:p)
  df_train$Y <- df_train$X1 - df_train$X2 + rnorm(n_tr, sd = 0.5)

  df_test <- as.data.frame(replicate(p, rnorm(n_te)))
  names(df_test) <- paste0("X", 1:p)
  df_test$Y <- df_test$X1 - df_test$X2 + rnorm(n_te, sd = 0.5)

  pval <- permrf(df_train, df_test, B = 60, tr_prop = 0.7)

  expect_true(is.numeric(pval))
  expect_length(pval, 1)
  expect_true(is.finite(pval))
  expect_true(pval >= 0 && pval <= 1)
})

test_that("permrf validates tr_prop", {
  df_train <- data.frame(X1 = rnorm(10), Y = rnorm(10))
  df_test  <- data.frame(X1 = rnorm(10), Y = rnorm(10))

  expect_error(permrf(df_train, df_test, tr_prop = 0))
  expect_error(permrf(df_train, df_test, tr_prop = 1))
  expect_error(permrf(df_train, df_test, tr_prop = -0.2))
})
