test_that("mars_classification generator returns expected dimensions/columns", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("LaplacesDemon") # if you rely on rbern()
  set.seed(1)

  n <- 50
  n_nuis <- 3
  df <- mars_classification(n = n, n_nuisance = n_nuis, eps = 0.5)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), n)
  expect_true("Y" %in% names(df))
  expect_true(all(df$Y %in% c(0, 1)))

  # X_1..X_5 + nuisance + Y
  expect_equal(ncol(df), 5 + n_nuis + 1)
})

test_that("lm_generation_classification generator returns Y in {0,1} and correct columns", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("LaplacesDemon") # if you rely on rbern()
  set.seed(2)

  n <- 40
  beta <- c(1, -1, 0.5)
  p <- length(beta)
  n_nuis <- 2

  df <- lm_generation_classification(
    n = n,
    beta_hat = beta,
    mean_shift = 0.2,
    variance_shift = 1.0,
    cor = 0.3,
    n_nuisance = n_nuis,
    eps = 0.1
  )

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), n)
  expect_true(all(df$Y %in% c(0, 1)))

  # X1..Xp + X_nuis1..X_nuisK + Y
  expect_true(all(paste0("X", 1:p) %in% names(df)))
  if (n_nuis > 0) expect_true(all(paste0("X_nuis", 1:n_nuis) %in% names(df)))
})
