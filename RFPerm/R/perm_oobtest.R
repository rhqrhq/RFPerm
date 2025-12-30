# ======================================================================================
# Exported: RFPerm (Comparing the mean of OOB error vs test error) via permutation tests
# ======================================================================================

#' RFPerm OOB-vs-test permutation test using random forests (ranger)
#'
#' @description
#' Fits a random forest on a reference/training batch, compares the training OOB squared
#' errors with the test-batch squared prediction errors, and returns a one-sided
#' permutation p-value for increased test error (drift).
#'
#' @details
#' INPUT CONTRACT:
#' \itemize{
#'   \item \code{df_train}, \code{df_test} must be \code{data.frame}s.
#'   \item Response column must be named \code{Y}.
#'   \item Predictors are assumed to be in columns \code{1:(ncol(df)-1)} and response in the last column.
#'   \item \code{df_train} and \code{df_test} must have the same predictor columns (names/order recommended).
#'   \item Predictors should be numeric (or compatible with \code{ranger}).
#' }
#'
#' The p-value is computed as the upper-tail probability under the permutation null:
#' \deqn{p = 1 - \hat{F}_{\mathrm{perm}}(d_0)}
#' where \eqn{d_0 = \bar{e}_{test} - \bar{e}_{train}}.
#'
#' @param df_train Data frame containing predictors and response \code{Y}.
#' @param df_test Data frame containing the same predictors and response \code{Y}.
#' @param classification Character flag retained for compatibility (currently unused).
#' @param B Integer scalar; number of permutations, by default 5000.
#' @param scaling Logical; if TRUE, scales predictors within each dataset using \code{scale()}.
#'
#' @return Numeric scalar p-value (one-sided; small values indicate increased test error).
#'
#' @examples
#' # df_train and df_test must contain response column named Y
#' # pval <- perm_oobtest(df_train, df_test, B = 2000)
#'
#' @export
perm_oobtest <- function(df_train, df_test,
                         classification = "FALSE",
                         B = 5000,
                         scaling = FALSE) {
  if (!is.data.frame(df_train) || !is.data.frame(df_test)) stop("`df_train` and `df_test` must be data.frames.")
  if (!("Y" %in% names(df_train)) || !("Y" %in% names(df_test))) stop("Both data frames must contain a response column named `Y`.")
  if (ncol(df_train) != ncol(df_test)) stop("`df_train` and `df_test` must have the same number of columns.")
  if (!all(names(df_train)[-ncol(df_train)] == names(df_test)[-ncol(df_test)])) {
    warning("Predictor column names/order differ between train and test; ensure alignment.")
  }

  n_col <- ncol(df_train)

  if (isTRUE(scaling)) {
    df_train[, 1:(n_col - 1)] <- as.data.frame(scale(df_train[, 1:(n_col - 1)], scale = TRUE))
    df_test[, 1:(n_col - 1)]  <- as.data.frame(scale(df_test[,  1:(n_col - 1)], scale = TRUE))
  }

  # Heuristic hyperparameters
  min_node_size_train <- round(sqrt(nrow(df_train)) / 2)
  mtry_fit <- round(n_col / 3)

  # Fit RF; keep.inbag is required for OOB bookkeeping
  model1 <- ranger::ranger(
    Y ~ .,
    data = df_train,
    num.trees = 150,
    min.node.size = min_node_size_train,
    mtry = mtry_fit,
    sample.fraction = 0.8,
    keep.inbag = TRUE
  )

  oob_error <- oob_error_observation_(
    model1,
    df_train[, 1:(n_col - 1)],
    response = df_train[, n_col]
  )$oob_error_obs

  prediction_error <- (df_test$Y - predict(model1, df_test[, 1:(n_col - 1)])$predictions)^2

  result <- permTest_(oob_error, prediction_error, B = B)
  p_val <- 1 - stats::ecdf(result$d_list)(result$d0)
  p_val
}