

########################################################################################
#' RFPerm for the binary response:
#'@param df_exist, df_new:         Numeric matrix/data.frame of dimensions
#                                  (n_exist, p) and (n_new, p) respectively.
#'@param B:                        Number of Permutations B, by default = 5000
#'@param scaling:                  Scaling the covariates or not.
#'@return p_val:                   Numeric p-value for the test
########################################################################################

permoob_classification <- function(df_exist, df_new, B = 5000,
  scaling = TRUE){
  n_col <- ncol(df_exist)
  if(scaling == TRUE){
    df_exist[, 1:(n_col - 1)] <- as.data.frame(scale(df_exist[, 1:(n_col - 1)], scale = TRUE))
    df_new[, 1:(n_col - 1)] <- as.data.frame(scale(df_new[, 1:(n_col - 1)], scale = TRUE))
  }
  n_exist <- nrow(df_exist)
  min_node_size <- round(sqrt(n_exist)/2)
  mtry_fit <- round(sqrt(n_col))
  model1 <- ranger(as.factor(Y)~., data = df_exist,
    num.trees = 150, min.node.size = min_node_size,
    mtry = mtry_fit, sample.fraction = 0.5, keep.inbag = TRUE,
    probability = TRUE)
  error <- oob_error_rate_brier(model1, 
    data = df_exist[, 1:(n_col - 1)], response = df_exist$Y,
    data_new = df_new[, 1:(n_col - 1)], response_new = df_new$Y)
  oob_error <- error$oob_error
  pred_error <- error$pred_error
  result <- permTest_(oob_error, pred_error, B = B)
  d_null <- result$d0 
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(p_val)
}
