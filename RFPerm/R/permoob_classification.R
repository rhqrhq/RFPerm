library(MASS)
library(ranger)

# ================================================================================
# Internal helper: OOB + test Brier errors for Random forests on Binary Response
# ================================================================================

#' Compute per-observation Brier errors for OOB (train) and test (new) (internal)
#'
#' @details
#' Input form (contract):
#' \itemize{
#'   \item \code{model} is a \code{ranger::ranger} classification model fit with:
#'     \code{probability = TRUE} and \code{keep.inbag = TRUE}.
#'   \item \code{data} is the predictor data used for training (n rows).
#'   \item \code{response} is integer 0/1 or logical; length n.
#'   \item \code{data_new} is predictor data for the new/test sample (m rows).
#'   \item \code{response_new} is integer 0/1 or logical; length m.
#' }
#'
#' This function extracts per-tree predicted probability for class "1" and computes
#' the Brier score \eqn{(p - y)^2}:
#' - OOB: average Brier score across trees where observation is OOB
#' - Test: average Brier score across all trees
#'
#' @param model A ranger probability forest (classification) with inbag counts.
#' @param data Predictor data.frame/matrix for the training sample.
#' @param response Vector length n, values 0/1.
#' @param data_new Predictor data.frame/matrix for the test sample.
#' @param response_new Vector length m, values 0/1.
#'
#' @return list(oob_error = numeric(n), pred_error = numeric(m))
#'
#' @keywords internal
.oob_brier_errors <- function(model, data, response, data_new, response_new) {
  if (is.null(model$inbag.counts)) stop("`model` must be fit with keep.inbag = TRUE.")
  response <- as.integer(response)
  response_new <- as.integer(response_new)

  n <- model$num.samples
  n_tree <- model$num.trees
  if (length(response) != n) stop("`response` must have length model$num.samples.")
  if (nrow(data) != n) stop("`data` must have nrow(data) == model$num.samples.")

  # Per-tree class-prob predictions on training data: [n, n_class, n_tree]
  pred_train <- predict(model, data, predict.all = TRUE)$predictions

  # Identify which class index corresponds to label "1"
  class_levels <- dimnames(pred_train)[[2]]
  if (is.null(class_levels)) stop("Could not read class labels from ranger predictions.")
  idx1 <- match("1", class_levels)
  if (is.na(idx1)) {
    # Fallback: if Y was factor with levels c("0","1"), this should work.
    stop("Class label '1' not found in ranger prediction dimnames; ensure Y is coded 0/1.")
  }

  # Extract p_hat(class=1) for each tree
  p_train <- matrix(NA_real_, nrow = n, ncol = n_tree)
  for (k in 1:n_tree) p_train[, k] <- as.numeric(pred_train[, idx1, k])

  # Mark in-bag predictions as NA for OOB computation
  list_inbag <- lapply(model$inbag.counts, function(x) which(x > 0))
  for (k in 1:n_tree) p_train[list_inbag[[k]], k] <- NA_real_

  # OOB Brier error per observation: mean over OOB trees
  oob_error <- numeric(n)
  for (i in 1:n) oob_error[i] <- mean((p_train[i, ] - response[i])^2, na.rm = TRUE)

  # Test/new sample: per-tree p_hat(class=1), use all trees (no OOB masking)
  m <- nrow(data_new)
  pred_new <- predict(model, data_new, predict.all = TRUE)$predictions

  p_new <- matrix(NA_real_, nrow = m, ncol = n_tree)
  for (k in 1:n_tree) p_new[, k] <- as.numeric(pred_new[, idx1, k])

  pred_error <- numeric(m)
  for (i in 1:m) pred_error[i] <- mean((p_new[i, ] - response_new[i])^2, na.rm = TRUE)

  return(list(oob_error = oob_error, pred_error = pred_error))
}


# ==============================================================================
# Exported: classification RFPerm (OOB vs test) using Brier score
# ==============================================================================

#' RFPerm for classification using Brier score (OOB vs test) with ranger
#'
#' @description
#' Fits a probability random forest on an "existing/reference" sample and compares:
#' \itemize{
#'   \item OOB Brier errors on the reference sample; vs
#'   \item prediction Brier errors on the new/test sample.
#' }
#' A permutation test on the pooled error vectors yields a one-sided p-value for
#' increased test error (drift).
#'
#' @details
#' Input form (contract):
#' \itemize{
#'   \item \code{df_exist}, \code{df_new} are data.frames.
#'   \item Must contain response column \code{Y} coded as 0/1 (integer or numeric).
#'   \item Predictors are all non-\code{Y} columns (used via \code{as.factor(Y) ~ .}).
#'   \item Train and test should have the same predictor set.
#' }
#'
#' @param df_exist Existing/reference data.frame with response \code{Y} in {0,1}.
#' @param df_new New/test data.frame with response \code{Y} in {0,1}.
#' @param B Integer scalar; number of permutations.
#' @param scaling Logical; if TRUE, scales predictors within each dataset using \code{scale()}.
#' @param num.trees Integer; number of trees.
#'
#' @return Numeric scalar p-value (upper-tail; small values indicate drift).
#'
#' @export
permoob_classification <- function(df_exist, df_new, B = 5000, scaling = TRUE, num.trees = 150) {
  if (!is.data.frame(df_exist) || !is.data.frame(df_new)) stop("`df_exist` and `df_new` must be data.frames.")
  if (!("Y" %in% names(df_exist)) || !("Y" %in% names(df_new))) stop("Both data frames must contain a response column named `Y`.")
  if (!is.numeric(B) || length(B) != 1 || B < 1) stop("`B` must be an integer scalar >= 1.")
  if (!is.numeric(num.trees) || length(num.trees) != 1 || num.trees < 1) stop("`num.trees` must be >= 1.")

  # Validate Y coding
  y1 <- as.integer(df_exist$Y)
  y2 <- as.integer(df_new$Y)
  if (any(!is.finite(y1)) || any(!is.finite(y2))) stop("`Y` must be finite and coded 0/1.")
  if (!all(y1 %in% c(0L, 1L)) || !all(y2 %in% c(0L, 1L))) stop("`Y` must be coded as 0/1.")

  pred_cols_exist <- setdiff(names(df_exist), "Y")
  pred_cols_new   <- setdiff(names(df_new), "Y")
  if (!setequal(pred_cols_exist, pred_cols_new)) warning("Train/test predictor sets differ; ensure consistent predictors.")

  # Optional scaling (predictors only)
  if (isTRUE(scaling)) {
    df_exist[, pred_cols_exist] <- as.data.frame(scale(df_exist[, pred_cols_exist], scale = TRUE))
    df_new[,   pred_cols_new]   <- as.data.frame(scale(df_new[,   pred_cols_new],   scale = TRUE))
  }

  n_exist <- nrow(df_exist)
  n_col <- ncol(df_exist)

  min_node_size <- round(sqrt(n_exist) / 2)
  mtry_fit <- max(1L, round(sqrt(n_col)))  # heuristic retained from your original code

  # Fit probability forest. Note: We explicitly factorize Y with levels 0/1.
  model <- ranger::ranger(
    stats::as.formula("as.factor(Y) ~ ."),
    data = df_exist,
    num.trees = as.integer(num.trees),
    min.node.size = min_node_size,
    mtry = mtry_fit,
    sample.fraction = 0.5,
    keep.inbag = TRUE,
    probability = TRUE
  )

  errs <- .oob_brier_errors(
    model,
    data = df_exist[, pred_cols_exist, drop = FALSE],
    response = df_exist$Y,
    data_new = df_new[, pred_cols_new, drop = FALSE],
    response_new = df_new$Y
  )

  # Permutation test:
  res <- permTest_(errs$oob_error, errs$pred_error, B = as.integer(B))
  p_val <- 1 - stats::ecdf(res$d_list)(res$d0)
  p_val
}
