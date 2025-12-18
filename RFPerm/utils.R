#permutation test for the covariate shift:
library(glmnet)
library(ranger)
library(densratio)
library(knockoff)
library(MASS)
library(ranger)
library(randomForest)
library(pROC)
library(LaplacesDemon)
library(xgboost)
library(caret)
library(lightgbm)
mars_model_df <- function(n, p_nuisance, sigma){
  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)
  X3 <- runif(n, 0, 1)
  X4 <- runif(n, 0, 1)
  X5 <- runif(n, 0, 1)
  Y1 <- 0.1 * exp(4 * X1) + 
    4/(1+exp(-20 * (X2 - 0.5))) + 3 * X3 +
    2 * X4 + X5
  random_noise <- rnorm(n, 0, sigma)
  Y <- Y1 + random_noise
  if(p_nuisance > 0){
    #nuisance parameters with standard normal random noise:
    X_nuiss <- matrix(NA, n, p_nuisance)
    for(i in 1:p_nuisance){
      X_nuiss[,i] <- rnorm(n, 0, 1)
    }
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <-
      c("Y1", paste("X", c(1:(5+as.numeric(p_nuisance))), sep = ""), "Y")
  }
  else{
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <- 
      c("Y1", paste("X", c(1:5), sep = ""), "Y")
  }
  return(data.frame(data_mat))
}

mars_model_df_normal <- function(n, p_nuisance, sigma){
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  X3 <- rnorm(n, 0, 1)
  X4 <- rnorm(n, 0, 1)
  X5 <- rnorm(n, 0, 1)
  #Y1 <- 0.1 * exp(4 * X1) +
  #  4/(1+exp(-20 * (X2 - 0.5))) + 3 * X3 +
  #  2 * X4 + X5
  Y1 <- 50 * sin(pi * X1 * X2) + 5 * (X3 - 0.05)^2 + 5 * X4 + 5 * X5
  random_noise <- rnorm(n, 0, sigma)
  Y <- Y1 + random_noise
  if(p_nuisance > 0){
    #nuisance parameters with standard normal random noise:
    X_nuiss <- matrix(NA, n, p_nuisance)
    for(i in 1:p_nuisance){
      X_nuiss[,i] <- rnorm(n, 0, 1)
    }
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <-
      c("Y1", paste("X", c(1:(5+as.numeric(p_nuisance))), sep = ""), "Y")
  }
  else{
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <-
      c("Y1", paste("X", c(1:5), sep = ""), "Y")
  }
  return(data.frame(data_mat))
}
correlation_model_df <- function(n, p, snr, rho = 0.35){
  corr_matrix <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in 1:p){
      corr_matrix[i, j] <- rho ^ (abs(i - j))
    }
  }
  X_design <- mvrnorm(n, mu = rep(0, p), corr_matrix)
  beta <- as.matrix(rep(1, p))
  Y1 <- X_design %*% beta
  sigma_noise_square <- t(beta) %*% corr_matrix %*% beta / snr
  Y <- Y1 + rnorm(n, 0, sqrt(sigma_noise_square))
  data_mat <- cbind(Y1, X_design, Y)
  colnames(data_mat) <- c("Y1", paste("X", c(1:as.numeric(p)), sep = ""), "Y")
  return(data.frame(data_mat))
}




#Gettign the tree-specific error in the out-of-bag error
oob_error_result_ <- function(model, data, response){
  n <- model$num.samples
  n_tree <- model$num.trees
  oob_pred <- matrix(0, nrow = n, ncol = n_tree)
  list_oob <- lapply(model$inbag.counts, function(x){
    which(x > 0)
  })
  pred_matrix <- predict(model, data, predict.all = TRUE)$predictions
  oob_error <- numeric(n_tree)
  for(i in 1:n_tree){
    oob_error[i] <- mean((response[list_oob[[i]]] - pred_matrix[list_oob[[i]], i])^2)
  }
  return(list(oob_error = oob_error, list_oob = list_oob))
}
#getting the oob error for each of the observation:
oob_error_observation_ <- function(model, data, response){
  n <- model$num.samples
  n_tree <- model$num.trees
  oob_pred <- matrix(0, nrow = n, ncol = n_tree)
  list_inbag <- lapply(model$inbag.counts, function(x){
    which(x > 0)
  })
  pred_matrix <- predict(model, data, predict.all = TRUE)$predictions
  #getting the observations-specific OOB error:
  for(k in 1:n_tree){
    #pred_matrix[, k] <- (response - pred_matrix[, k])^2
    pred_matrix[list_inbag[[k]], k] <- NA
  }
  mean_prediction <- apply(pred_matrix, 1, function(x){
    mean(x[!is.na(x)])
  })
  oob_error_obs <- (response - mean_prediction)^2
  #oob_error_obs <- apply(pred_matrix, 1, function(x){
  #  mean(x[!is.na(x)])
  #})
  return(list(oob_error_obs = oob_error_obs, list_inbag = list_inbag))
}

#Conduct the permutation test
permTest_ <- function(error_list1, error_list2, B = 1000){
  d0 <- mean(error_list2) - mean(error_list1)
  n1 <- length(error_list1)
  n2 <- length(error_list2)
  whole_list <- c(error_list1, error_list2)
  d_list <- numeric(B)
  for(i in 1:B){
    sample_idx <- sample(c(1:(n1+n2)), n1, replace = FALSE)
    list1 <- whole_list[sample_idx]
    list2 <- whole_list[-sample_idx]
    d_list[i] <- mean(list1) - mean(list2)
  }
  return(list(d0 = d0, d_list = d_list))
}


#################
#
perm_oobtest <- function(df_train, df_test, 
  classification = "FALSE", B = 5000, scaling = FALSE){
  n_col <- ncol(df_train)
  if(scaling == TRUE){
    df_train[, 1:(n_col - 1)] <- as.data.frame(scale(df_train[,1:(n_col - 1)], scale = TRUE))
    df_test[, 1:(n_col - 1)] <- as.data.frame(scale(df_test[,1:(n_col - 1)], scale = TRUE))
  }
  #setting of the hyperparameters by default
  min_node_size_train <- round(sqrt(nrow(df_train))/2)
  mtry_fit = round(n_col/3)
  #Fit the random forest by default parameter setting:
  model1 <- ranger(Y~., data = df_train, 
    num.trees = 150, min.node.size = min_node_size_train,
    mtry = mtry_fit, sample.fraction = 0.8, keep.inbag = TRUE)
  oob_error <- oob_error_observation_(model1, df_train[,1:(n_col - 1)], response = df_train[, n_col])$oob_error_obs
  prediction_error <- (df_test$Y - predict(model1, df_test[, 1:(n_col - 1)])$predictions)^2
  result <- permTest_(oob_error, prediction_error, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(p_val)
}

#Output by recording the mean squared error results:
perm_oobtest_withMSE <- function(df_train, df_test,
  classification = "FALSE", B = 5000, scaling = FALSE){
  n_col <- ncol(df_train)
  if(scaling == TRUE){
    df_train[, 1:(n_col - 1)] <- as.data.frame(scale(df_train[,1:(n_col - 1)], scale = TRUE))
    df_test[, 1:(n_col - 1)] <- as.data.frame(scale(df_test[,1:(n_col - 1)], scale = TRUE))
  }
  #setting of the hyperparameters by default
  min_node_size_train <- round(sqrt(nrow(df_train))/2)
  mtry_fit = round(n_col/3)
  #Fit the random forest by default parameter setting:
  model1 <- ranger(Y~., data = df_train, 
    num.trees = 150, min.node.size = min_node_size_train,
    mtry = mtry_fit, sample.fraction = 0.8, keep.inbag = TRUE)
  oob_error <- oob_error_observation_(model1, df_train[,1:(n_col - 1)], response = df_train[, n_col])$oob_error_obs
  prediction_error <- (df_test$Y - predict(model1, df_test[, 1:(n_col - 1)])$predictions)^2
  #Return the list of oob_error and prediction_error:
  mean_oob <- mean(oob_error)
  sd_oob <- sd(oob_error)
  mean_pred <- mean(prediction_error)
  sd_pred <- sd(prediction_error)
  result <- permTest_(oob_error, prediction_error, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(list(p_val = p_val,
    MSE_val = oob_error,
    MSE_test = prediction_error))
}

#adapting XGBoost as the component model to conduct the 
permXGB_nontuning <- function(df1, df2, train_prop = 0.5, tune_prop = 0.2, B = 5000){
  n_1 <- nrow(df1)
  n_2 <- nrow(df2)
  n_tr <- round(n_1 * train_prop)
  n_tu <- round(n_1 * tune_prop)
  n_col <- ncol(df1)
  n_val <- n_1 - n_tr - n_tu
  samp_ind1 <- sample(c(1:n_1), n_tr, replace = FALSE)
  samp_ind2 <- sample(c(1:n_1)[-samp_ind1], n_tu, replace = FALSE)
  df_tr <- df1[samp_ind1, ]
  df_tu <- df1[samp_ind2, ]
  df_val <- df1[-c(samp_ind1, samp_ind2), ]
  #comparing the validation error with the test error for the xgboost:
  dtrain <- xgb.DMatrix(data = as.matrix(df_tr[,1:(n_col-1)]), label = df_tr$Y)
  dtune <- xgb.DMatrix(data = as.matrix(df_tu[,1:(n_col-1)]), label = df_tu$Y)
  ####
  watchlist <- list(train=dtrain, test=dtune)
  xgb_model <- xgb.train(param <- list(max_depth = 4, eta = 0.2, subsample = 0.8, colsample_bytree = 0.8,
    lambda = 0.2, 
    eval_metric = "rmse"
    ), dtrain, watchlist = watchlist, nrounds = 150, early_stopping_rounds = 5)

  #doing the prediction on the df_val, conduct the prediction here:
  pred_val <- predict(xgb_model, xgb.DMatrix(data = as.matrix(df_val[,1:(n_col-1)])))
  pred_test <- predict(xgb_model, xgb.DMatrix(data = as.matrix(df2[,1:(n_col-1)])))
  MSE_val <- (df_val$Y - pred_val)^2
  MSE_test <- (df2$Y - pred_test)^2
  #record the full MSE:
  result <- permTest_(MSE_val, MSE_test, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(list(p_val = p_val, MSE_val = MSE_val, MSE_test = MSE_test))
}
#Let's start without tuning and see the robustness of the power:

permXGB_tuning <- function(df1, df2, train_prop = 0.5, tune_prop = 0.2, B = 5000,
  depth_grid = c(3,4,5,6,7), eta_grid = c(0.1,0.2,0.5,1),
  lambda_grid = c(0.1,0.2,0.3,0.5)){
  n_1 <- nrow(df1)
  n_2 <- nrow(df2)
  n_tr <- round(n_1 * train_prop)
  n_tu <- round(n_1 * tune_prop)
  n_col <- ncol(df1)
  n_val <- n_1 - n_tr - n_tu
  samp_ind1 <- sample(c(1:n_1), n_tr, replace = FALSE)
  samp_ind2 <- sample(c(1:n_1)[-samp_ind1], n_tu, replace = FALSE)
  df_tr <- df1[samp_ind1, ]
  df_tu <- df1[samp_ind2, ]
  df_val <- df1[-c(samp_ind1, samp_ind2), ]
  param_grid <- expand.grid(depth = depth_grid, eta = eta_grid, lambda = lambda_grid)
  param_grid$best_score <- 0
  #comparing the validation error with the test error for the xgboost:
  dtrain <- xgb.DMatrix(data = as.matrix(df_tr[,1:(n_col-1)]), label = df_tr$Y)
  dtune <- xgb.DMatrix(data = as.matrix(df_tu[,1:(n_col-1)]), label = df_tu$Y)
  ####
  watchlist <- list(train=dtrain, test=dtune)  
  for(i in 1:nrow(param_grid)){
    max_depth <- param_grid$depth[i]
    eta <- param_grid$eta[i]
    lambda <- param_grid$lambda[i]
    xgb_model <- xgb.train(param <- list(max_depth = max_depth, eta = eta, subsample = 0.8, colsample_bytree = 0.8,
    lambda = lambda, 
    eval_metric = "rmse"
    ), dtrain, watchlist = watchlist, nrounds = 150, early_stopping_rounds = 5)
    param_grid$best_score[i] <- xgb_model$best_score
    param_grid$best_iter[i] <- xgb_model$best_iteration
  }
  #The best model: 
  mi <- which.min(param_grid$best_score)
  best_depth <- param_grid$depth[mi]
  best_eta <- param_grid$eta[mi]
  best_lambda <- param_grid$lambda[mi]
  best_iter <- param_grid$best_iter[mi]
  xgb_model <- xgb.train(param <- list(max_depth = best_depth, eta = best_eta, subsample = 0.8, colsample_bytree = 0.8,
    lambda = best_lambda, eval_metric = "rmse"), dtrain, watchlist = watchlist, nrounds = best_iter)
  ####
  pred_val <- predict(xgb_model, xgb.DMatrix(data = as.matrix(df_val[,1:(n_col-1)])))
  pred_test <- predict(xgb_model, xgb.DMatrix(data = as.matrix(df2[,1:(n_col-1)])))
  MSE_val <- (df_val$Y - pred_val)^2
  MSE_test <- (df2$Y - pred_test)^2
  result <- permTest_(MSE_val, MSE_test, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(p_val)
}

#permutation test for random forest on separate validation set:
permrf <- function(df_train, df_test, B = 5000, tr_prop = 0.7){
  n_col <- ncol(df_train)
  #setting of the hyperparameters by default
  min_node_size_train <- round(sqrt(nrow(df_train))/2)
  mtry_fit = round(n_col/3)
  tr_ind = sample(1:nrow(df_train), round(tr_prop * nrow(df_train)), replace = FALSE)
  df_tr <- df_train[tr_ind, ]
  df_val <- df_train[-tr_ind, ]
  model1 <- ranger(Y~., data = df_tr, 
    num.trees = 150, min.node.size = min_node_size_train,
    mtry = mtry_fit, sample.fraction = 0.8, keep.inbag = TRUE)
  validation_error <- (df_val$Y - predict(model1, df_val[, 1:(n_col - 1)])$predictions)^2
  prediction_error <- (df_test$Y - predict(model1, df_test[, 1:(n_col - 1)])$predictions)^2
  result <- permTest_(validation_error, prediction_error, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(p_val)  
}


#adapting XGB for classification:
permXGBC_nontuning <- function(df1, df2, train_prop = 0.5, tune_prop = 0.2, B = 5000){
  n_1 <- nrow(df1)
  n_2 <- nrow(df2)
  n_tr <- round(n_1 * train_prop)
  n_tu <- round(n_1 * tune_prop)
  n_col <- ncol(df1)
  n_val <- n_1 - n_tr - n_tu
  samp_ind1 <- sample(c(1:n_1), n_tr, replace = FALSE)
  samp_ind2 <- sample(c(1:n_1)[-samp_ind1], n_tu, replace = FALSE)
  df_tr <- df1[samp_ind1, ]
  df_tu <- df1[samp_ind2, ]
  df_val <- df1[-c(samp_ind1, samp_ind2), ]
  #comparing the validation error with the test error for the xgboost:
  dtrain <- xgb.DMatrix(data = as.matrix(df_tr[,1:(n_col-1)]), label = df_tr$Y)
  dtune <- xgb.DMatrix(data = as.matrix(df_tu[,1:(n_col-1)]), label = df_tu$Y)
  ####
  watchlist <- list(train=dtrain, test=dtune)
  xgb_model <- xgb.train(param <- list(max_depth = round(sqrt(n_col-1)), eta = 0.075, subsample = 0.8, colsample_bytree = 0.8,
    lambda = 0.2, 
    objective = "binary:logistic", eval_metric = 'auc'
    ), dtrain, watchlist = watchlist, nrounds = 150, early_stopping_rounds = 5)
  #doing the prediction on the df_val, conduct the prediction here:
  pred_val <- predict(xgb_model, xgb.DMatrix(data = as.matrix(df_val[,1:(n_col-1)])))
  pred_test <- predict(xgb_model, xgb.DMatrix(data = as.matrix(df2[,1:(n_col-1)])))
  MSE_val <- (df_val$Y - pred_val)^2
  MSE_test <- (df2$Y - pred_test)^2
  #record the full MSE:
  result <- permTest_(MSE_val, MSE_test, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(list(p_val = p_val, MSE_val = MSE_val, MSE_test = MSE_test))
}

#adapting the lightgbm for regression:
permLGB_nontuning <- function(df1, df2, train_prop = 0.5, tune_prop = 0.2, B = 5000){
  n_1 <- nrow(df1)
  n_2 <- nrow(df2)
  n_tr <- round(n_1 * train_prop)
  n_tu <- round(n_1 * tune_prop)
  n_col <- ncol(df1)
  n_val <- n_1 - n_tr - n_tu
  samp_ind1 <- sample(c(1:n_1), n_tr, replace = FALSE)
  samp_ind2 <- sample(c(1:n_1)[-samp_ind1], n_tu, replace = FALSE)
  df_tr <- df1[samp_ind1, ]
  df_tu <- df1[samp_ind2, ]
  df_val <- df1[-c(samp_ind1, samp_ind2), ]
  dtrain <- lgb.Dataset(data = as.matrix(df_tr[,1:(n_col-1)]), label = df_tr$Y)
  dtest <- lgb.Dataset.create.valid(dtrain, as.matrix(df_tu[,1:(n_col-1)]), label = df_tu$Y)
  param_list <- list(objective = "regression", metric = "l2",
    learning_rate = 0.15, max_depth = 4
    )
  valids <- list(test = dtest)
  lgb_model <- lgb.train(params = param_list, data = dtrain, nrounds = 150,
    valids = valids, early_stopping_rounds = 5)
  pred_val <- predict(lgb_model, newdata = as.matrix(df_val[,1:(n_col-1)]))
  pred_test <- predict(lgb_model, newdata = as.matrix(df2[,1:(n_col-1)]))
  MSE_val <- (df_val$Y - pred_val)^2
  MSE_test <- (df2$Y - pred_test)^2
  #record the full MSE:
  result <- permTest_(MSE_val, MSE_test, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(list(p_val = p_val, MSE_val = MSE_val, MSE_test = MSE_test))
}


permLGBC_nontuning <- function(df1, df2, train_prop = 0.5, tune_prop = 0.2, B = 5000){
  n_1 <- nrow(df1)
  n_2 <- nrow(df2)
  n_tr <- round(n_1 * train_prop)
  n_tu <- round(n_1 * tune_prop)
  n_col <- ncol(df1)
  n_val <- n_1 - n_tr - n_tu
  samp_ind1 <- sample(c(1:n_1), n_tr, replace = FALSE)
  samp_ind2 <- sample(c(1:n_1)[-samp_ind1], n_tu, replace = FALSE)
  df_tr <- df1[samp_ind1, ]
  df_tu <- df1[samp_ind2, ]
  df_val <- df1[-c(samp_ind1, samp_ind2), ]
  dtrain <- lgb.Dataset(data = as.matrix(df_tr[,1:(n_col-1)]), label = df_tr$Y)
  dtest <- lgb.Dataset.create.valid(dtrain, as.matrix(df_tu[,1:(n_col-1)]), label = df_tu$Y)
  param_list <- list(objective = "binary", metric = "auc",
    learning_rate = 0.15, max_depth = 4
    )
  valids <- list(test = dtest)
  lgb_model <- lgb.train(params = param_list, data = dtrain, nrounds = 150,
    valids = valids, early_stopping_rounds = 5)
  pred_val <- predict(lgb_model, newdata = as.matrix(df_val[,1:(n_col-1)]))
  pred_test <- predict(lgb_model, newdata = as.matrix(df2[,1:(n_col-1)]))
  MSE_val <- (df_val$Y - pred_val)^2
  MSE_test <- (df2$Y - pred_test)^2
  #record the full MSE:
  result <- permTest_(MSE_val, MSE_test, B = B)
  #the null value:
  d_null <- result$d0
  h2 <- ecdf(unlist(result$d_list))
  d2 <- result$d0
  h_val <- h2(d2)
  p_val <- 1 - h_val
  return(list(p_val = p_val, MSE_val = MSE_val, MSE_test = MSE_test))
}

LM_generation <- function(n, beta_hat, cor, n_nuisance, eps){
  p <- length(beta_hat)
  corr_matrix <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in 1:p){
      corr_matrix[i, j] <- cor ^ (abs(i - j))
    }
  }
  X_design <- mvrnorm(n, mu = rep(0, p), corr_matrix)
  Y1 <- as.matrix(X_design) %*% as.matrix(beta_hat, nrow = p)
  random_error <- rnorm(n, 0, eps)
  #adding the nuisance random noise features:
  X_nuiss <- matrix(0, n, n_nuisance)
  for(i in 1:n_nuisance){
    X_nuiss[,i] <- rnorm(n, 0, 1)
  }
  Y <- Y1 + random_error
  df_return <- data.frame(cbind(Y1, X_design, X_nuiss, Y))
  ncol_df <- ncol(df_return)
  colnames(df_return) <- c("Y1", paste("X", c(1:p), sep = ""), paste("X_nuis", c(1:n_nuisance), sep = ""), "Y")
  X_return <- df_return[,2:(ncol_df - 1)]
  return(list(df_return = df_return, X_return = X_return))
}












