library(ranger)
library(MASS)
library(Rlab)


#' Generate MARS-like regression data (Uniform covariates) with optional nuisances
#'
#' @description
#' Generates a regression dataset with a nonlinear mean function and Gaussian noise.
#' The returned data frame contains the noise-free signal `Y1`, covariates `X*`, and
#' noisy response `Y`.
#'
#' @param n Integer scalar. Number of observations (must be >= 1).
#' @param p_nuisance Integer scalar. Number of nuisance covariates (must be >= 0).
#' @param sigma Numeric scalar. Standard deviation of additive Gaussian noise (must be >= 0).
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item \code{Y1}: numeric, noise-free signal
#'   \item \code{X1..X_{5+p_nuisance}}: numeric covariates
#'   \item \code{Y}: numeric response
#' }
#'
#' @examples
#' df <- mars_model_df(n = 200, p_nuisance = 10, sigma = 1)
#' str(df)
#'
#' @export
mars_model_df <- function(n, p_nuisance, sigma) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("`n` must be an integer scalar >= 1.")
  if (!is.numeric(p_nuisance) || length(p_nuisance) != 1 || p_nuisance < 0) stop("`p_nuisance` must be >= 0.")
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma < 0) stop("`sigma` must be >= 0.")

  X1 <- stats::runif(n, 0, 1)
  X2 <- stats::runif(n, 0, 1)
  X3 <- stats::runif(n, 0, 1)
  X4 <- stats::runif(n, 0, 1)
  X5 <- stats::runif(n, 0, 1)

  Y1 <- 0.1 * exp(4 * X1) +
    4 / (1 + exp(-20 * (X2 - 0.5))) + 3 * X3 +
    2 * X4 + X5

  Y <- Y1 + stats::rnorm(n, 0, sigma)

  if (p_nuisance > 0) {
    X_nuiss <- matrix(NA_real_, n, p_nuisance)
    for (i in 1:p_nuisance) X_nuiss[, i] <- stats::rnorm(n, 0, 1)
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <- c("Y1", paste0("X", 1:(5 + as.numeric(p_nuisance))), "Y")
  } else {
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <- c("Y1", paste0("X", 1:5), "Y")
  }

  data.frame(data_mat)
}

#' Generate MARS-like regression data (Normal covariates) with optional nuisances
#'
#' @param n Integer scalar. Number of observations (must be >= 1).
#' @param p_nuisance Integer scalar. Number of nuisance covariates (must be >= 0).
#' @param sigma Numeric scalar. Standard deviation of additive Gaussian noise (must be >= 0).
#'
#' @return A data.frame with columns \code{Y1}, \code{X*}, \code{Y}. See \code{\link{mars_model_df}}.
#'
#' @export
mars_model_df_normal <- function(n, p_nuisance, sigma) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("`n` must be an integer scalar >= 1.")
  if (!is.numeric(p_nuisance) || length(p_nuisance) != 1 || p_nuisance < 0) stop("`p_nuisance` must be >= 0.")
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma < 0) stop("`sigma` must be >= 0.")

  X1 <- stats::rnorm(n, 0, 1)
  X2 <- stats::rnorm(n, 0, 1)
  X3 <- stats::rnorm(n, 0, 1)
  X4 <- stats::rnorm(n, 0, 1)
  X5 <- stats::rnorm(n, 0, 1)

  Y1 <- 50 * sin(pi * X1 * X2) + 5 * (X3 - 0.05)^2 + 5 * X4 + 5 * X5
  Y <- Y1 + stats::rnorm(n, 0, sigma)

  if (p_nuisance > 0) {
    X_nuiss <- matrix(NA_real_, n, p_nuisance)
    for (i in 1:p_nuisance) X_nuiss[, i] <- stats::rnorm(n, 0, 1)
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <- c("Y1", paste0("X", 1:(5 + as.numeric(p_nuisance))), "Y")
  } else {
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <- c("Y1", paste0("X", 1:5), "Y")
  }

  data.frame(data_mat)
}

#' Generate correlated Gaussian design with linear signal and target SNR
#'
#' @param n Integer scalar. Number of observations (>= 1).
#' @param p Integer scalar. Number of covariates (>= 1).
#' @param snr Numeric scalar. Signal-to-noise ratio (must be > 0).
#' @param rho Numeric scalar. Correlation parameter (typically in [0,1)).
#'
#' @return A data.frame with columns \code{Y1}, \code{X1..Xp}, \code{Y}.
#'
#' @export
correlation_model_df <- function(n, p, snr, rho = 0.35) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("`n` must be an integer scalar >= 1.")
  if (!is.numeric(p) || length(p) != 1 || p < 1) stop("`p` must be an integer scalar >= 1.")
  if (!is.numeric(snr) || length(snr) != 1 || snr <= 0) stop("`snr` must be > 0.")
  if (!is.numeric(rho) || length(rho) != 1) stop("`rho` must be a numeric scalar.")

  corr_matrix <- matrix(NA_real_, p, p)
  for (i in 1:p) for (j in 1:p) corr_matrix[i, j] <- rho^(abs(i - j))

  X_design <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = corr_matrix)
  beta <- matrix(1, nrow = p, ncol = 1)

  Y1 <- X_design %*% beta
  sigma_noise_sq <- as.numeric(t(beta) %*% corr_matrix %*% beta) / snr
  Y <- Y1 + stats::rnorm(n, 0, sqrt(sigma_noise_sq))

  data_mat <- cbind(Y1, X_design, Y)
  colnames(data_mat) <- c("Y1", paste0("X", 1:p), "Y")
  data.frame(data_mat)
}

# ==============================================================================
# Synthetic data generator 1: MARS-like classification
# ==============================================================================

#' Generate MARS-like binary classification data with nuisance features
#'
#' @description
#' Generates a binary outcome \code{Y} from a nonlinear logistic model in five core
#' covariates plus optional nuisance covariates. This is useful for benchmarking
#' drift tests in classification settings.
#'
#' @details
#' Input form (contract):
#' \itemize{
#'   \item \code{n}: integer scalar >= 1
#'   \item \code{n_nuisance}: integer scalar >= 0
#'   \item \code{eps}: numeric scalar >= 0; noise injected into the logit
#'   \item \code{mean_seq}, \code{sd_seq}: numeric vectors of length 5 (means and SDs of X1..X5)
#' }
#'
#' The class probability is:
#' \deqn{p(x) = \mathrm{logit}^{-1}\left(\frac{10\sin(\pi X_1 X_2) + 20(X_3-0.05)^2 + 10X_4 + 5X_5 - 20}{3} + \varepsilon\right)}
#' where \eqn{\varepsilon \sim N(0, \mathrm{eps}^2)}.
#'
#' @param n Integer scalar; number of observations.
#' @param n_nuisance Integer scalar; number of nuisance features (>= 0).
#' @param eps Numeric scalar; logit noise SD (>= 0).
#' @param mean_seq Numeric vector length 5; means for X1..X5.
#' @param sd_seq Numeric vector length 5; standard deviations for X1..X5 (must be > 0).
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item \code{X_1..X_5}: numeric core features
#'   \item \code{X_nuis_1..X_nuis_k}: numeric nuisance features (Normal(0,1))
#'   \item \code{Y}: integer in \{0,1\}
#' }
#'
#' @export
mars_classification <- function(
  n,
  n_nuisance,
  eps = 1,
  mean_seq = c(0, 0, 0, 0, 0),
  sd_seq   = c(1, 1, 1, 1, 1)
) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("`n` must be an integer scalar >= 1.")
  if (!is.numeric(n_nuisance) || length(n_nuisance) != 1 || n_nuisance < 0) stop("`n_nuisance` must be >= 0.")
  if (!is.numeric(eps) || length(eps) != 1 || eps < 0) stop("`eps` must be >= 0.")
  if (!is.numeric(mean_seq) || length(mean_seq) != 5) stop("`mean_seq` must be numeric length 5.")
  if (!is.numeric(sd_seq) || length(sd_seq) != 5 || any(sd_seq <= 0)) stop("`sd_seq` must be numeric length 5 with all entries > 0.")

  X1 <- stats::rnorm(n, mean_seq[1], sd_seq[1])
  X2 <- stats::rnorm(n, mean_seq[2], sd_seq[2])
  X3 <- stats::rnorm(n, mean_seq[3], sd_seq[3])
  X4 <- stats::rnorm(n, mean_seq[4], sd_seq[4])
  X5 <- stats::rnorm(n, mean_seq[5], sd_seq[5])

  # Latent logit (plus optional additive noise)
  eta <- (10 * sin(pi * X1 * X2) +
            20 * (X3 - 0.05)^2 +
            10 * X4 + 5 * X5 - 20) / 3 +
    stats::rnorm(n, 0, eps)

  p <- 1 / (1 + exp(-eta))
  Y <- rbern(n, p)

  # Nuisance features: standard Normal
  if (n_nuisance > 0) {
    X_nuis <- matrix(stats::rnorm(n * n_nuisance, 0, 1), nrow = n, ncol = n_nuisance)
    colnames(X_nuis) <- paste0("X_nuis_", seq_len(n_nuisance))
    df <- data.frame(X_1 = X1, X_2 = X2, X_3 = X3, X_4 = X4, X_5 = X5, X_nuis, Y = Y)
  } else {
    df <- data.frame(X_1 = X1, X_2 = X2, X_3 = X3, X_4 = X4, X_5 = X5, Y = Y)
  }

  df
}

# ==============================================================================
# Synthetic data generator 2: Correlated Gaussian + logistic transform
# ==============================================================================

#' Generate correlated Gaussian features and a logistic binary outcome (with shift)
#'
#' @description
#' Draws a correlated Gaussian design \eqn{X \in R^p} and generates binary outcome
#' \code{Y} via a logistic transform of a linear predictor plus Gaussian noise.
#' Includes optional mean/variance shifts to create covariate shift scenarios.
#'
#' @details
#' Input form (contract):
#' \itemize{
#'   \item \code{n}: integer scalar >= 1
#'   \item \code{beta_hat}: numeric vector length p (linear coefficients)
#'   \item \code{mean_shift}: numeric scalar (applied to all feature means)
#'   \item \code{variance_shift}: numeric scalar > 0 (scales covariance)
#'   \item \code{cor}: numeric scalar (typically in [0,1)) correlation parameter
#'   \item \code{n_nuisance}: integer scalar >= 0
#'   \item \code{eps}: numeric scalar >= 0, noise SD on the latent score
#' }
#'
#' @param n Integer scalar; number of observations.
#' @param beta_hat Numeric vector length p; linear coefficients.
#' @param mean_shift Numeric scalar; shift in feature means.
#' @param variance_shift Numeric scalar > 0; scales feature covariance (variance inflation).
#' @param cor Numeric scalar; AR(1)-like correlation parameter.
#' @param n_nuisance Integer scalar >= 0; number of nuisance features (Uniform(0,1)).
#' @param eps Numeric scalar >= 0; latent noise SD.
#'
#' @return A data.frame with columns \code{X1..Xp}, nuisance \code{X_nuis*}, and \code{Y} in {0,1}.
#'
#' @export
lm_generation_classification <- function(
  n,
  beta_hat,
  mean_shift = 0,
  variance_shift = 1,
  cor,
  n_nuisance,
  eps
) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("`n` must be an integer scalar >= 1.")
  if (!is.numeric(beta_hat) || length(beta_hat) < 1) stop("`beta_hat` must be a numeric vector of length >= 1.")
  if (!is.numeric(mean_shift) || length(mean_shift) != 1) stop("`mean_shift` must be a numeric scalar.")
  if (!is.numeric(variance_shift) || length(variance_shift) != 1 || variance_shift <= 0) stop("`variance_shift` must be a numeric scalar > 0.")
  if (!is.numeric(cor) || length(cor) != 1) stop("`cor` must be a numeric scalar.")
  if (!is.numeric(n_nuisance) || length(n_nuisance) != 1 || n_nuisance < 0) stop("`n_nuisance` must be >= 0.")
  if (!is.numeric(eps) || length(eps) != 1 || eps < 0) stop("`eps` must be >= 0.")

  p <- length(beta_hat)

  # Toeplitz covariance: Sigma_{ij} = cor^{|i-j|}
  corr_matrix <- matrix(NA_real_, p, p)
  for (i in 1:p){
    for (j in 1:p){
      corr_matrix[i, j] <- cor^(abs(i - j))
    } 
  }

  # MASS::mvrnorm expects a covariance matrix; we scale it by variance_shift.
  X_design <- MASS::mvrnorm(
    n,
    mu = rep(mean_shift, p),
    Sigma = corr_matrix * variance_shift
  )

  score <- as.matrix(X_design) %*% matrix(beta_hat, nrow = p) + stats::rnorm(n, 0, eps)
  pY <- 1 / (1 + exp(-score))
  Y <- rbern(n, pY)

  # Nuisance features: Uniform(0,1) (as in your original code)
  if (n_nuisance > 0) {
    X_nuis <- matrix(stats::runif(n * n_nuisance, 0, 1), nrow = n, ncol = n_nuisance)
    colnames(X_nuis) <- paste0("X_nuis", seq_len(n_nuisance))
  } else {
    X_nuis <- NULL
  }

  df_return <- data.frame(X_design)
  colnames(df_return) <- paste0("X", seq_len(p))

  if (!is.null(X_nuis)) df_return <- cbind(df_return, X_nuis)
  df_return$Y <- Y

  df_return
}

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
