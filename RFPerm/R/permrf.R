# ==============================================================================
# Optional exported: validation-split RF permutation test
# ==============================================================================

#' Random-forest permutation test using a train/validation split within df_train
#'
#' @details
#' INPUT CONTRACT:
#' \itemize{
#'   \item \code{df_train}, \code{df_test} are data.frames with response column \code{Y}.
#'   \item Predictors are columns \code{1:(ncol(df)-1)} and \code{Y} is the response.
#'   \item \code{tr_prop} is in (0,1) and controls the split within \code{df_train}.
#' }
#'
#' Fits RF on a subset of df_train, compares validation squared errors vs test squared errors,
#' and returns a one-sided permutation p-value.
#'
#' @param df_train Data frame with predictors and response \code{Y}.
#' @param df_test Data frame with predictors and response \code{Y}.
#' @param B Integer scalar; number of permutations.
#' @param tr_prop Numeric scalar in (0,1); training proportion within df_train.
#'
#' @return Numeric scalar p-value.
#'
#' @export
permrf <- function(df_train, df_test, B = 5000, tr_prop = 0.7) {
  if (!is.data.frame(df_train) || !is.data.frame(df_test)) stop("`df_train` and `df_test` must be data.frames.")
  if (!("Y" %in% names(df_train)) || !("Y" %in% names(df_test))) stop("Both data frames must contain `Y`.")
  if (!is.numeric(tr_prop) || length(tr_prop) != 1 || tr_prop <= 0 || tr_prop >= 1) stop("`tr_prop` must be in (0,1).")

  n_col <- ncol(df_train)

  min_node_size_train <- round(sqrt(nrow(df_train)) / 2)
  mtry_fit <- round(n_col / 3)

  tr_ind <- sample(seq_len(nrow(df_train)), round(tr_prop * nrow(df_train)), replace = FALSE)
  df_tr <- df_train[tr_ind, ]
  df_val <- df_train[-tr_ind, ]

  model1 <- ranger::ranger(
    Y ~ .,
    data = df_tr,
    num.trees = 150,
    min.node.size = min_node_size_train,
    mtry = mtry_fit,
    sample.fraction = 0.8,
    keep.inbag = TRUE
  )

  validation_error <- (df_val$Y - predict(model1, df_val[, 1:(n_col - 1)])$predictions)^2
  prediction_error <- (df_test$Y - predict(model1, df_test[, 1:(n_col - 1)])$predictions)^2

  result <- permTest_(validation_error, prediction_error, B = B)
  p_val <- 1 - stats::ecdf(result$d_list)(result$d0)
  p_val
}
