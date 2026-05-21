from Utils import *
import xgboost 
import numpy as np 
import pandas as pd 
import statistics
import random
import os  
import re 

#50 covariates with 8 signal features, most obvious one:
df_tr = LMGeneration(n = 5000, p = 50, beta = np.ones(8), cor = 0.3, eps_noi = 1, mean = 0)
df_te = LMGeneration(n = 5000, p = 50, beta = [1,1,1,1,0,0,0,0], cor = 0.3, eps_noi = 1, mean = 0)
for i in range(10, 200, 10):
    df_pred_mat_tr, df_pred_mat_te, df_pred_mat_tr_B, df_pred_mat_te_B, df_MSE_mat_te, df_MSE_mat_te_B = sequence_dependence_XGB(df_tr, df_te, n_boost_total = 250, start = i,  B = 200)
    df_pred_mat_mean_tr = np.zeros((i, 250 - i))
    df_pred_mat_sd_tr = np.zeros((i, 250 - i))
    df_pred_mat_mean_te = np.zeros((i, 250 - i))
    df_pred_mat_sd_te = np.zeros((i, 250 - i))
    r1 = df_pred_mat_sd_te.shape[0]
    c1 = df_pred_mat_sd_te.shape[1]
    for j in range(r1):
        for k in range(c1):
            list1 = []
            list2 = []
            for l in range(200):
                list1.append(np.corrcoef(df_pred_mat_tr[:, j].tolist(), df_pred_mat_tr_B[l, :, k].tolist())[0,1])
                list2.append(np.corrcoef(df_pred_mat_te[:, j].tolist(), df_pred_mat_te_B[l, :, k].tolist())[0,1])
            df_pred_mat_mean_tr[j, k] = np.mean(list1)
            df_pred_mat_sd_tr[j, k] = np.std(list1)
            df_pred_mat_mean_te[j, k] = np.mean(list2)
            df_pred_mat_sd_te[j, k] = np.std(list2)  
            print(k)
    df_MSE_mean_te_B = df_MSE_mat_te_B.mean(axis = 0)
    df_MSE_sd_te_B = df_MSE_mat_te_B.std(axis = 0)
    np.savetxt("df_pred_mean_tr_start_" + str(i) + "_mat_xgb1.csv", df_pred_mat_mean_tr, delimiter = ",")
    np.savetxt("df_pred_sd_tr_start_" + str(i) + "_mat_xgb1.csv", df_pred_mat_sd_tr, delimiter = ",")
    np.savetxt("df_pred_mean_te_start_" + str(i) + "_mat_xgb1.csv", df_pred_mat_mean_te, delimiter = ",")
    np.savetxt("df_pred_sd_te_start_" + str(i) + "_mat_xgb1.csv", df_pred_mat_sd_te, delimiter = ",")   
    np.savetxt("df_MSE_mean_te_start" + str(i) + "_mat_xgb1.csv", df_MSE_mean_te_B, delimiter = ",")
    np.savetxt("df_MSE_sd_te_start" + str(i) + "_mat_xgb1.csv", df_MSE_sd_te_B, delimiter = ",")


#50 Covariates with 8 covariates, most obvious pattern:
df_tr = LMGeneration(n = 5000, p = 50, beta = [1,1,1,1,1,1,1,1], cor = 0.3, eps_noi = 1, mean = 0)
df_te = LMGeneration(n = 5000, p = 50, beta = [1,1,1,1,0,0,0,0], cor = 0.3, eps_noi = 1, mean = 0)
for i in range(20, 125, 5):
    df_pred_mat_tr, df_pred_mat_te, df_pred_mat_tr_B, df_pred_mat_te_B = sequence_dependence_LGB(df_tr, df_te, n_boost_total = 150, start = i, 
        B = 500)
    df_pred_mat_mean_tr = np.zeros((i, 150 - i))
    df_pred_mat_sd_tr = np.zeros((i, 150 - i))
    df_pred_mat_mean_te = np.zeros((i, 150 - i))
    df_pred_mat_sd_te = np.zeros((i, 150 - i))
    for j in range(i):
        for k in range(150 - i):
            list1 = []
            list2 = []
            for l in range(500):
                list1.append(np.corrcoef(df_pred_mat_tr[:, j].tolist(), df_pred_mat_tr_B[l, :, k].tolist())[0,1])
                list2.append(np.corrcoef(df_pred_mat_te[:, j].tolist(), df_pred_mat_te_B[l, :, k].tolist())[0,1])
            df_pred_mat_mean_tr[j, k] = np.mean(list1)
            df_pred_mat_sd_tr[j, k] = np.std(list1)
            df_pred_mat_mean_te[j, k] = np.mean(list2)
            df_pred_mat_sd_te[j, k] = np.std(list2)  
        print(j)          
    np.savetxt("df_pred_mean_tr_start_" + str(i) + "_mat_lgb1.csv", df_pred_mat_mean_tr, delimiter = ",")
    np.savetxt("df_pred_sd_tr_start_" + str(i) + "_mat_lgb1.csv", df_pred_mat_sd_tr, delimiter = ",")
    np.savetxt("df_pred_mean_te_start_" + str(i) + "_mat_lgb1.csv", df_pred_mat_mean_te, delimiter = ",")
    np.savetxt("df_pred_sd_te_start_" + str(i) + "_mat_lgb1.csv", df_pred_mat_sd_te, delimiter = ",")   


#50covariates, 30/15 signal features, less obvious:
#LightGBM version:
beta_train = np.ones(30)
beta_test = np.ones(15)
beta_test += np.array([1] * 15)
df_tr = LMGeneration(n = 5000, p = 50, beta = beta_train, cor = 0.3, eps_noi = 1, mean = 0)
df_te = LMGeneration(n = 5000, p = 50, beta = beta_test, cor = 0.3, eps_noi = 1, mean = 0)
for i in range(20, 125, 5):
    df_pred_mat_tr, df_pred_mat_te, df_pred_mat_tr_B, df_pred_mat_te_B = sequence_dependence_LGB(df_tr, df_te, n_boost_total = 150, start = i, 
        B = 500)
    df_pred_mat_mean_tr = np.zeros((i, 150 - i))
    df_pred_mat_sd_tr = np.zeros((i, 150 - i))
    df_pred_mat_mean_te = np.zeros((i, 150 - i))
    df_pred_mat_sd_te = np.zeros((i, 150 - i))
    for j in range(i):
        for k in range(150 - i):
            list1 = []
            list2 = []
            for l in range(500):
                list1.append(np.corrcoef(df_pred_mat_tr[:, j].tolist(), df_pred_mat_tr_B[l, :, k].tolist())[0,1])
                list2.append(np.corrcoef(df_pred_mat_te[:, j].tolist(), df_pred_mat_te_B[l, :, k].tolist())[0,1])
            df_pred_mat_mean_tr[j, k] = np.mean(list1)
            df_pred_mat_sd_tr[j, k] = np.std(list1)
            df_pred_mat_mean_te[j, k] = np.mean(list2)
            df_pred_mat_sd_te[j, k] = np.std(list2)  
        print(j)          
    np.savetxt("df_pred_mean_tr_start_" + str(i) + "_mat_lgb2.csv", df_pred_mat_mean_tr, delimiter = ",")
    np.savetxt("df_pred_sd_tr_start_" + str(i) + "_mat_lgb2.csv", df_pred_mat_sd_tr, delimiter = ",")
    np.savetxt("df_pred_mean_te_start_" + str(i) + "_mat_lgb2.csv", df_pred_mat_mean_te, delimiter = ",")
    np.savetxt("df_pred_sd_te_start_" + str(i) + "_mat_lgb2.csv", df_pred_mat_sd_te, delimiter = ",") 

#XGBoost Version:
df_tr = LMGeneration(n = 5000, p = 50, beta = beta_train, cor = 0.3, eps_noi = 1, mean = 0)
df_te = LMGeneration(n = 5000, p = 50, beta = beta_test, cor = 0.3, eps_noi = 1, mean = 0)

for i in range(5, 100):
    df_pred_mat_tr, df_pred_mat_te, df_pred_mat_tr_B, df_pred_mat_te_B = sequence_dependence_XGB(df_tr, df_te, n_boost_total = 150, start = i, 
        B = 500)
    df_pred_mat_mean_tr = np.zeros((i, 150 - i))
    df_pred_mat_sd_tr = np.zeros((i, 150 - i))
    df_pred_mat_mean_te = np.zeros((i, 150 - i))
    df_pred_mat_sd_te = np.zeros((i, 150 - i))
    for j in range(i):
        for k in range(150 - i):
            list1 = []
            list2 = []
            for l in range(100):
                list1.append(np.corrcoef(df_pred_mat_tr[:, j].tolist(), df_pred_mat_tr_B[l, :, k].tolist())[0,1])
                list2.append(np.corrcoef(df_pred_mat_te[:, j].tolist(), df_pred_mat_te_B[l, :, k].tolist())[0,1])
            df_pred_mat_mean_tr[j, k] = np.mean(list1)
            df_pred_mat_sd_tr[j, k] = np.std(list1)
            df_pred_mat_mean_te[j, k] = np.mean(list2)
            df_pred_mat_sd_te[j, k] = np.std(list2)  
        print(j)          
    np.savetxt("df_pred_mean_tr_start_" + str(i) + "_mat_xgb2.csv", df_pred_mat_mean_tr, delimiter = ",")
    np.savetxt("df_pred_sd_tr_start_" + str(i) + "_mat_xgb2.csv", df_pred_mat_sd_tr, delimiter = ",")
    np.savetxt("df_pred_mean_te_start_" + str(i) + "_mat_xgb2.csv", df_pred_mat_mean_te, delimiter = ",")
    np.savetxt("df_pred_sd_te_start_" + str(i) + "_mat_xgb2.csv", df_pred_mat_sd_te, delimiter = ",")   


'''
Save the dataframe and reformulate into R:
power_df <- data.frame(
  start55_mean = numeric(58),
  start60_mean = numeric(58),
  start65_mean = numeric(58), 
  start70_mean = numeric(58),   
  start75_mean = numeric(58),
  start80_mean = numeric(58),
  start85_mean = numeric(58),
  start90_mean = numeric(58),
  start55_sd = numeric(58),
  start60_sd = numeric(58),
  start65_sd = numeric(58),
  start70_sd = numeric(58),  
  start75_sd = numeric(58),    
  start80_sd = numeric(58),
  start85_sd = numeric(58),  
  start90_sd = numeric(58)
  )
power_df$start55_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_55_mat_xgb1.csv"))[1, 1:58])
power_df$start60_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_60_mat_xgb1.csv"))[1, 1:58])
power_df$start65_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_65_mat_xgb1.csv"))[1, 1:58])
power_df$start70_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_70_mat_xgb1.csv"))[1, 1:58])
power_df$start75_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_75_mat_xgb1.csv"))[1, 1:58])
power_df$start80_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_80_mat_xgb1.csv"))[1, 1:58])
power_df$start85_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_85_mat_xgb1.csv"))[1, 1:58])
power_df$start90_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_90_mat_xgb1.csv"))[1, 1:58])
power_df$start55_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_55_mat_xgb1.csv"))[1, 1:58])
power_df$start60_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_60_mat_xgb1.csv"))[1, 1:58])
power_df$start65_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_65_mat_xgb1.csv"))[1, 1:58])
power_df$start70_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_70_mat_xgb1.csv"))[1, 1:58])
power_df$start75_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_75_mat_xgb1.csv"))[1, 1:58])
power_df$start80_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_80_mat_xgb1.csv"))[1, 1:58])
power_df$start85_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_85_mat_xgb1.csv"))[1, 1:58])
power_df$start90_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_90_mat_xgb1.csv"))[1, 1:58])
power_df <- data.frame(as.matrix(power_df))
power_df.to_csv("xgb_longterm_correlation_results_startfrom55.csv")

power_df <- data.frame(
  start5_mean = numeric(100),
  start10_mean = numeric(100),
  start15_mean = numeric(100),
  start20_mean = numeric(100),
  start25_mean = numeric(100),  
  start30_mean = numeric(100),  
  start35_mean = numeric(100),  
  start40_mean = numeric(100),  
  start5_sd = numeric(100),
  start10_sd = numeric(100),
  start15_sd = numeric(100),
  start20_sd = numeric(100),
  start25_sd = numeric(100), 
  start30_sd = numeric(100), 
  start35_sd = numeric(100), 
  start40_sd = numeric(100)
  )
power_df$start5_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_5_mat_xgb2.csv"))[4, 1:100])
power_df$start10_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_10_mat_xgb2.csv"))[9, 1:100])
power_df$start15_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_15_mat_xgb2.csv"))[14, 1:100])
power_df$start20_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_20_mat_xgb2.csv"))[19, 1:100])
power_df$start25_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_25_mat_xgb2.csv"))[24, 1:100])
power_df$start30_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_30_mat_xgb2.csv"))[29, 1:100])
power_df$start35_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_35_mat_xgb2.csv"))[34, 1:100])
power_df$start40_mean <- c(as.matrix(read.csv("df_pred_mean_te_start_40_mat_xgb2.csv"))[39, 1:100])
power_df$start5_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_5_mat_xbg2.csv"))[4, 1:100])
power_df$start10_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_10_mat_xbg2.csv"))[9, 1:100])
power_df$start15_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_15_mat_xbg2.csv"))[14, 1:100])
power_df$start20_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_20_mat_xbg2.csv"))[19, 1:100])
power_df$start25_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_25_mat_xbg2.csv"))[24, 1:100])
power_df$start30_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_30_mat_xbg2.csv"))[29, 1:100])
power_df$start35_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_35_mat_xbg2.csv"))[34, 1:100])
power_df$start40_sd <- c(as.matrix(read.csv("df_pred_sd_te_start_40_mat_xbg2.csv"))[39, 1:100])
power_df <- data.frame(power_df)

power_df.to_csv("XGB_correlation_start_from_5th_tree.csv")
'''










