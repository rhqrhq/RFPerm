test_that("permTest_ returns correct structure and lengths", {
  set.seed(1)
  e1 <- c(1, 2, 3)
  e2 <- c(2, 4, 6, 8)

  out <- permTest_(e1, e2, B = 20)

  expect_true(is.list(out))
  expect_true(all(c("d0", "d_list") %in% names(out)))
  expect_type(out$d0, "double")
  expect_type(out$d_list, "double")
  expect_length(out$d_list, 20)

  # d0 definition: mean(error_list2) - mean(error_list1)
  expect_equal(out$d0, mean(e2) - mean(e1))
})

test_that("permTest_ validates B", {
  e1 <- c(1, 2, 3)
  e2 <- c(2, 4, 6)

  expect_error(permTest_(e1, e2, B = 0))
  expect_error(permTest_(e1, e2, B = -5))
  expect_error(permTest_(e1, e2, B = "10"))
})