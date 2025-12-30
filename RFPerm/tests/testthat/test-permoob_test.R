test_that("perm_oobtest returns a valid p-value in [0,1]", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  set.seed(42)
  n_tr <- 80
  n_te <- 70
  p <- 4

  Xtr <- replicate(p, rnorm(n_tr))
  Xte <- replicate(p, rnorm(n_te))

  df_train <- as.data.frame(Xtr)
  names(df_train) <- paste0("X", 1:p)
  df_train$Y <- rowSums(df_train[, 1:p]) + rnorm(n_tr, sd = 0.5)

  df_test <- as.data.frame(Xte)
  names(df_test) <- paste0("X", 1:p)
  df_test$Y <- rowSums(df_test[, 1:p]) + rnorm(n_te, sd = 0.5)

  pval <- perm_oobtest(df_train, df_test, B = 80, scaling = FALSE)

  expect_true(is.numeric(pval))
  expect_length(pval, 1)
  expect_true(is.finite(pval))
  expect_gte(pval, 0)
  expect_lte(pval, 1)
})

test_that("perm_oobtest works with scaling=TRUE", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  set.seed(43)
  n_tr <- 60
  n_te <- 60

  df_train <- data.frame(X1 = rnorm(n_tr), X2 = rnorm(n_tr))
  df_train$Y <- 2 * df_train$X1 + rnorm(n_tr)

  df_test <- data.frame(X1 = rnorm(n_te), X2 = rnorm(n_te))
  df_test$Y <- 2 * df_test$X1 + rnorm(n_te)

  pval <- perm_oobtest(df_train, df_test, B = 50, scaling = TRUE)

  expect_true(is.numeric(pval))
  expect_true(pval >= 0 && pval <= 1)
})

test_that("perm_oobtest validates inputs", {
  df1 <- data.frame(X1 = 1:5, Y = 1:5)
  df2 <- data.frame(X1 = 1:5, Y = 1:5)

  expect_error(perm_oobtest(list(), df2), "data.frames")
  expect_error(perm_oobtest(df1, list()), "data.frames")

  df_bad <- data.frame(X1 = 1:5, Z = 1:5)
  expect_error(perm_oobtest(df_bad, df2), "response column")
  expect_error(perm_oobtest(df1, df_bad), "response column")

  df3 <- data.frame(X1 = 1:5, X2 = 1:5, Y = 1:5)
  expect_error(perm_oobtest(df1, df3), "same number of columns")
})
