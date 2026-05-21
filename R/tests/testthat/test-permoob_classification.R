test_that("permoob_classification returns a valid p-value in [0,1]", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  set.seed(202)
  n1 <- 120
  n2 <- 100

  make_df <- function(n) {
    X1 <- rnorm(n)
    X2 <- rnorm(n)
    p  <- 1 / (1 + exp(-(1.2 * X1 - 0.8 * X2)))
    Y  <- rbinom(n, size = 1, prob = p)
    data.frame(X1 = X1, X2 = X2, Y = Y)
  }

  df_exist <- make_df(n1)
  df_new  <- make_df(n2)

  pval <- permoob_classification(df_exist, df_new, B = 80, scaling = TRUE, num.trees = 50)

  expect_true(is.numeric(pval))
  expect_length(pval, 1)
  expect_true(is.finite(pval))
  expect_true(pval >= 0 && pval <= 1)
})

test_that("permoob_classification errors if Y is not coded 0/1", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")

  df_exist <- data.frame(X1 = rnorm(30), X2 = rnorm(30), Y = sample(0:2, 30, TRUE))
  df_new   <- data.frame(X1 = rnorm(25), X2 = rnorm(25), Y = sample(0:1, 25, TRUE))

  expect_error(permoob_classification(df_exist, df_new, B = 20), "coded as 0/1")
})
