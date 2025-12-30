#' Permutation test for difference in mean errors (internal)
#'
#' @details
#' INPUT CONTRACT:
#' \itemize{
#'   \item \code{error_list1}, \code{error_list2} are numeric vectors (finite values recommended).
#'   \item \code{B} is an integer scalar >= 1.
#' }
#'
#' Observed statistic: \code{d0 = mean(error_list2) - mean(error_list1)}.
#' Under the null (exchangeability), we permute group membership on pooled errors.
#'
#' @param error_list1 Numeric vector (e.g., training/OOB errors).
#' @param error_list2 Numeric vector (e.g., test errors).
#' @param B Integer scalar; number of permutations.
#'
#' @return A list with \code{d0} (numeric scalar) and \code{d_list} (numeric vector length B).
#'
#' @keywords internal
permTest_ <- function(error_list1, error_list2, B = 1000) {
  if (!is.numeric(B) || length(B) != 1 || B < 1) stop("`B` must be an integer scalar >= 1.")
  error_list1 <- as.numeric(error_list1)
  error_list2 <- as.numeric(error_list2)

  d0 <- mean(error_list2) - mean(error_list1)
  n1 <- length(error_list1)
  n2 <- length(error_list2)

  whole_list <- c(error_list1, error_list2)
  d_list <- numeric(B)

  for (i in 1:B) {
    sample_idx <- sample(seq_len(n1 + n2), n1, replace = FALSE)
    list1 <- whole_list[sample_idx]
    list2 <- whole_list[-sample_idx]
    d_list[i] <- mean(list1) - mean(list2)
  }

  list(d0 = d0, d_list = d_list)
}