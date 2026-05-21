
#' Observation-level OOB squared error from a ranger model (internal)
#'
#' @details
#' INPUT CONTRACT:
#' \itemize{
#'   \item \code{model} must be a \code{ranger::ranger} model fit with \code{keep.inbag = TRUE}.
#'   \item \code{data} must be a data.frame/matrix of predictors with \code{nrow(data) == model$num.samples}.
#'   \item \code{response} must be a numeric vector of length \code{nrow(data)}.
#' }
#'
#' For each observation, the function averages predictions over trees where the observation
#' was OOB (in-bag predictions are set to NA), then returns squared error
#' \eqn{(y_i - \bar{\hat{y}}_{i,\mathrm{OOB}})^2}.
#'
#' @param model A fitted \code{ranger::ranger} regression model (must have \code{keep.inbag = TRUE}).
#' @param data Predictor data used for OOB predictions (data.frame or matrix).
#' @param response Numeric vector of true responses aligned to \code{data}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{oob_error_obs}: numeric vector length \code{nrow(data)} (may contain NA if no OOB trees)
#'   \item \code{list_inbag}: list length \code{model$num.trees} giving in-bag indices per tree
#' }
#'
#' @keywords internal
oob_error_observation_ <- function(model, data, response) {
  if (is.null(model$inbag.counts)) stop("`model` must be fit with keep.inbag = TRUE.")
  response <- as.numeric(response)
  if (length(response) != nrow(data)) stop("`response` must have length nrow(data).")

  n_tree <- model$num.trees
  pred_matrix <- predict(model, data, predict.all = TRUE)$predictions
  list_inbag <- lapply(model$inbag.counts, function(x) which(x > 0))

  # Mark in-bag predictions as NA; remaining values are OOB-only predictions.
  for (k in 1:n_tree) pred_matrix[list_inbag[[k]], k] <- NA_real_

  mean_prediction <- apply(pred_matrix, 1, function(x) mean(x[!is.na(x)]))
  oob_error_obs <- (response - mean_prediction)^2

  list(oob_error_obs = oob_error_obs, list_inbag = list_inbag)
}