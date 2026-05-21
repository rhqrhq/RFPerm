#GBDT long term sequence:
#The input dataframe is the result of that part of the, we are aiming to illustrate the issue here:
import xgboost 
import numpy as np 
import pandas as pd 
import statistics
import random
import os  
import re 
#Using the Toeplitz correlation structure:
def LMGeneration(n = 1000, p = 20, beta = [1,1,1,1,1,1,1], cor = 0.3, eps_noi = 1, mean = 0):
    n_signal = len(beta)
    mean_vector = [mean] * n_signal
    cov_mat = np.zeros((n_signal, n_signal))
    for i in range(n_signal):
        for j in range(n_signal):
            cov_mat[i, j] = cor ** abs(i - j)
    cov = np.dot(cov_mat, cov_mat.T)
    X_design = np.random.multivariate_normal(mean_vector, cov, size=n)
    #generate the response:
    Y = np.dot(X_design, beta) + np.random.normal(loc = 0, scale = eps_noi, size = n)
    X_noise = np.zeros((n, p - n_signal))
    for k in range(p - n_signal):
        X_noise[:, k] = np.random.normal(loc = 0, scale = 1, size = n)
    X_design = np.concatenate((X_design, X_noise), axis = 1)
    df_design = np.concatenate((X_design, Y.reshape(-1, 1)), axis = 1)
    return df_design
#XGBoost dependence:
'''
Input parameters:
df_tr:         Training DataFrame X, Y
df_te:         Prediction DataFrame X, Y
n_boost_total: Total Number of Boosting Rounds(by default can be 150-200)
start:         The starting iteration for the model prediction.
'''
def sequence_dependence_XGB(df_tr, df_te, n_boost_total, start, B = 200):
    df_tr_X = df_tr[:, :-1]
    df_tr_Y = df_tr[:, -1]
    df_te_X = df_te[:, :-1]
    df_te_Y = df_te[:, -1]
    n_tr = df_tr.shape[0]
    n_te = df_te.shape[0]
    df_pred_mat_tr = np.zeros((n_tr, start))
    df_pred_mat_te = np.zeros((n_te, start))
    df_MSE_mat_tr = np.zeros((n_tr, start))
    df_MSE_mat_te = np.zeros((n_te, start))
    #Prediction Matrix for the original dataset:
    df_pred_mat_tr_B = np.zeros((B, n_tr, n_boost_total - start))
    df_pred_mat_te_B = np.zeros((B, n_te, n_boost_total - start))
    df_MSE_mat_tr_B = np.zeros((B, n_tr, n_boost_total - start))
    df_MSE_mat_te_B = np.zeros((B, n_te, n_boost_total - start))
    params = {'objective': 'reg:squarederror', 'verbose': False,
          "random_state": 42, "colsample_bytree": 0.85,
          "colsample_bylevel": 0.85, "gamma": 0.0,
          "max_depth": 6, "reg_alpha": 0.0, 
          "learning_rate": 0.15, "reg_lambda": 0.25,
          "tol": 0.001}
    xg_train = xgboost.DMatrix(df_tr_X, label=df_tr_Y)  
    model_fit_initial = xgboost.train(params, xg_train, start)
    for k in range(start):
        df_pred_mat_tr[:, k] = model_fit_initial[k].predict(xgboost.DMatrix(df_tr_X))
        df_pred_mat_te[:, k] = model_fit_initial[k].predict(xgboost.DMatrix(df_te_X))
        df_MSE_mat_te[:, k] = (df_te_Y - df_pred_mat_te[:, k])**2
    for i in range(B):
        B_ind = np.random.choice(n_tr, n_tr, replace = True)
        B_X = df_tr_X[B_ind, :]
        B_Y = df_tr_Y[B_ind]
        xgb_new = model_fit_initial
        new_xg_train = xgboost.DMatrix(B_X, label = B_Y)
        model_new = xgboost.train(params, new_xg_train, n_boost_total - start, xgb_model = model_fit_initial)
        for j in range(0, n_boost_total - start):
            df_pred_mat_tr_B[i, :, j] = model_new[j+start].predict(xgboost.DMatrix(B_X))
            df_pred_mat_te_B[i, :, j] = model_new[j+start].predict(xgboost.DMatrix(df_te_X))
        for j in range(0, n_boost_total - start):
            df_MSE_mat_te_B[i, :, j] = (B_Y - df_pred_mat_te_B[i, :, j]) ** 2
        print(i)
    #return all of the predictions:
    return df_pred_mat_tr, df_pred_mat_te, df_pred_mat_tr_B, df_pred_mat_te_B, df_MSE_mat_te, df_MSE_mat_te_B
#LightGBM dependence:
'''
Input parameters:
df_tr:         Training DataFrame X, Y
df_te:         Prediction DataFrame X, Y
n_boost_total: Total Number of Boosting Rounds(by default can be 150-200)
start:         The starting iteration for the model prediction.
'''
def sequence_dependence_LGB(df_tr, df_te, n_boost_total, start, B = 200):
    df_tr_X = df_tr[:, :-1]
    df_tr_Y = df_tr[:, -1]
    df_te_X = df_te[:, :-1]
    df_te_Y = df_te[:, -1]
    n_tr = df_tr.shape[0]
    n_te = df_te.shape[0]
    df_pred_mat_tr = np.zeros((n_tr, start))
    df_pred_mat_te = np.zeros((n_te, start))
    #Prediction Matrix for the original dataset:
    df_pred_mat_tr_B = np.zeros((B, n_tr, n_boost_total - start))
    df_pred_mat_te_B = np.zeros((B, n_te, n_boost_total - start))
    params = {'objective': 'regression',
          "random_state": 42, "colsample_bytree": 0.85,
          "colsample_bylevel": 0.85, "gamma": 0.0,
          "max_depth": 6, "reg_alpha": 0.0, 
          "learning_rate": 0.15, "reg_lambda": 0.25,
          "tol": 0.001}
    xg_train = lightgbm.Dataset(df_tr_X, label=df_tr_Y)  
    model_fit_initial = lightgbm.train(params, xg_train, start)
    for k in range(start):
        df_pred_mat_tr[:, k] = model_fit_initial.predict(df_tr_X, start_iteration = k, num_iterations = 1)
        df_pred_mat_te[:, k] = model_fit_initial.predict(df_te_X, start_iteration = k, num_iterations = 1)
    for i in range(B):
        B_ind = np.random.choice(n_tr, n_tr, replace = True)
        B_X = df_tr_X[B_ind, :]
        B_Y = df_tr_Y[B_ind]
        xgb_new = model_fit_initial
        new_xg_train = lightgbm.Dataset(B_X, label = B_Y)
        model_new = lightgbm.train(params, new_xg_train, n_boost_total - start, init_model = model_fit_initial)
        for j in range(n_boost_total - start):
            df_pred_mat_tr_B[i, :, j] = model_new.predict(B_X, start_iteration = j+start, num_iterations = 1)
            df_pred_mat_te_B[i, :, j] = model_new.predict(df_te_X, start_iteration = j+start, num_iterations = 1)
        print(i)
    #return all of the predictions:
    return df_pred_mat_tr, df_pred_mat_te, df_pred_mat_tr_B, df_pred_mat_te_B



    