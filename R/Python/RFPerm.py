import os
import statistics
import sklearn
import numpy as np
import pandas as pd
import xgboost
from scipy.stats import combine_pvalues
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from typing import Iterable, Optional, Tuple, Union
from scipy.stats import multivariate_normal, chi2, t
from model_registry_class import ModelRegistry

#Generalize this to a general purposed permutation test where the mean of the validation error and the
#prediction error are compared.
def LM_generation(n, beta_hat, cor, n_nuisance, eps, mean_X = 0, var_X=1):
    beta_hat = np.asarray(beta_hat)
    p = len(beta_hat)
    corr_matrix = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            corr_matrix[i, j] = (cor ** abs(i - j)) * var_X
    X_design = multivariate_normal.rvs(
        mean=[mean_X] * p,
        cov=corr_matrix,
        size=n
    )
    Y1 = X_design @ beta_hat
    random_error = np.random.normal(0, eps, n)
    X_nuiss = np.random.normal(0, 1, size=(n, n_nuisance))
    Y = Y1 + random_error
    df_return = pd.DataFrame(
        np.column_stack([Y1, X_design, X_nuiss, Y])
    )
    colnames = (
        ["Y1"] +
        [f"X{i+1}" for i in range(p)] +
        [f"X_nuis{i+1}" for i in range(n_nuisance)] +
        ["Y"]
    )
    df_return.columns = colnames
    X_return = df_return.iloc[:, 1:-1]
    return {
        "df_return": df_return,
        "X_return": X_return
    }




def permTest(mse_list_val, mse_list_pred, n_perm = 5000, seed = 2026):
    '''
    Permutation test for one-sided or two-sided p-values to compare the mean of the 
    two sets of MSE
    '''
    rng = np.random.default_rng(seed)
    val = np.asarray(mse_list_val, dtype = float)
    pred = np.asarray(mse_list_pred, dtype = float)
    mean_pred = float(np.mean(val))
    mean_exist = float(np.mean(pred))
    mean_dif_original = mean_pred - mean_exist
    n_val = int(val.shape[0])
    n_pred = int(pred.shape[0])
    n = int(n_val + n_pred)
    gen = np.random.default_rng(seed = seed) if seed else np.random.default_rng()
    mean_dif_seq = np.zeros(n_perm)
    pooled_mse = np.empty(n, dtype = float)
    pooled_mse[:n_val] = val
    pooled_mse[n_val:] = pred
    pooled = pooled_mse.copy()
    for i in range(n_perm):
        pooled_ = rng.permutation(pooled)
        mean_dif_seq = np.mean(pooled[n_val:]) - np.mean(pooled[:n_val])
    #One-sided p-value:
    p_val_onesided = (1 + np.sum(mean_dif_original > mean_dif_seq))/(n_perm + 1)
    #Two-sided p-value:
    p_val_twosided = (1 + np.sum(np.abs(mean_dif_original) > np.abs(mean_dif_seq)))/(n_perm + 1)
    return p_val_onesided, p_val_twosided


def infer_response_type(Y):
    """
    Infer the types of response in Y:
    Y should either be Numeric or Object
    Return one of the "binary", "continuous" for the downstream modeling procedure.
    """
    if Y is None:
        raise ValueError("Y must not be None")
    Y_ = np.array(Y).reshape(-1)
    if Y_.size == 0:
        raise ValueError("Y must contain at least one observation")
    if pd.isna(Y_).any():
        raise ValueError("Y contains missing values. Please remove or impute them first.")
    dtype = Y_.dtype
    #Checking whether it's object or numeric, otherwise, ValueError
    is_numeric = np.issubdtype(dtype, np.number)
    is_object_like = dtype == object or np.issubdtype(dtype, np.str_)
    if not (is_numeric or is_object_like):
        raise TypeError(
            'Y must be numeric or categorical/object-like.'
            f"Got dtype={dtype}."
        )
    len_class = np.unique(Y_).size
    if len_class == 2:
        return 'binary'
    elif len_class > 6:
        return 'continuous'
    else:
        return 'multinomial'
    #checking the forms and conditions of this:




def PermValTest(df_exist, df_new, model_component, loss, test_size = 0.3, alpha = 0.05, B = 200, n_perm = 5000, seed = 2026):
    """
    Permutation Test for Distribution Shift via comparing the 
    mean of validation error and test error:
    The inner permutation test supplies one-sided p-value here.
    Parameters:
    - df_exist:         Existing Batch of data, (n_exist, p), dataframe.
    - df_new:           New Batch of data, (n_new, p), dataframe.
    - B:                Number of repeatances for the 
    - loss:             Notion of the loss function
    - alpha:            Significance Level
    - n_perm:           Number of the permutation
    - seed:             Random Seed
    """
    model = model_class[model_component]
    n_exist, p = df_exist.shape
    n_new = df_new.shape[0]
    X_exist = df_exist.iloc[:, :-1]
    Y_exist = df_exist.iloc[:, -1]
    X_new = df_new.iloc[:, :-1]
    Y_new = df_new.iloc[:, -1]
    response_type = infer_response_type(Y_exist)
    #Checking whether it aligns, or it would stop:
    pval_list = np.zeros(B)
    for i in range(B):
        X_train, Y_train, X_val, Y_val = train_test_split(
            X_exist, Y_exist, test_size = test_size,
            random_state = seed
        )
        fitted_model = model['fit'](
            X_train, Y_train, 
        )
        pred_val = model['predict'](
            fit = fitted_model,
            X_new = X_val
        )
        pred_test = model['predict'](
            fit = fitted_model, 
            X_new = X_new
        )
        p_one, _ = permTest(
            loss(Y_new, pred_test),
            loss(Y_val, pred_val)
        )
        pval_list[i] = p_one
    power = float(np.mean([1 if t < alpha else 0 for t in pval_list]))
    return power

def MSE(y_pred, y_target):
    return ((np.array(y_pred) - np.array(y_target))**2).tolist()

#The multinomial l2 metric is that argmin_{t}(Y_{t} - Y_target) ** 2 + \sum_{t}Y_{k!=t}-Y\
#y_pred: (n, n_c)
#y_target: ranges from 1 to n_c
def l2_multinomial(y_pred, y_target):
    n = y_pred.shape[0]
    c = y_pred.shape[1]
    y_true = np.zeros((n, c))
    y_target_index = (y_target - 1).astype(int)
    y_true[np.arange(n), y_target_index] = 1
    l2_loss = np.sum((y_true - y_pred) ** 2, axis = 1)
    return l2_loss

#The p-value aggregation procedure:
'''
Lancaster's procedure for the aggregation of the p-values:
(1) p-value lists: chi(sum(chi(1-p_value, weights)))
(2) weights for each of the p-value
'''
def lancaster_agg(p_values, weights):
    p_values = np.array(p_values)
    weights = np.array(weights)
    chi_stats = stats.chi2.ppf(1 - p_values, df = weights)
    total_stats = np.sum(chi_stats)
    total_df = np.sum(weights)
    p_value_combined = 1 - stats.chi2.cdf(total_stats, df = total_df)
    return p_value_combined


def bidir_adj(df_exist, df_new, )    
    


#Testing, for continuous 
df_exist = LM_generation(n = 1000, beta_hat = [1,1,1,1,1,1,1,1],
    cor = 0.3, n_nuisance = 12, eps = 3)['df_return']
df_new = LM_generation(n = 1000, beta_hat = [1,1,1,1,0,0,0,0],
    cor = 0.3, n_nuisance = 12, eps = 3)['df_return']
RFPerm(df_exist, df_new, loss = MSE, B = 100)







#Binary Outcome:
df_exist = LM_generation(n = 1000, beta_hat = [1,1,1,1,1,1,1,1],
    cor = 0.3, n_nuisance = 12, eps = 3)['df_return']
df_new = LM_generation(n = 1000, beta_hat = [1,1,1,1,0,0,0,0],
    cor = 0.3, n_nuisance = 12, eps = 3)['df_return']
df_exist.iloc[:, -1] = np.concatenate([np.zeros(500), np.ones(500)])
df_new.iloc[:, -1] = np.concatenate([np.ones(500), np.zeros(500)])
RFPerm(df_exist, df_new, loss = MSE, B = 100)


#Multinomial Outcome:
df_exist = LM_generation(n = 1000, beta_hat = [1,1,1,1,1,1,1,1],
    cor = 0.3, n_nuisance = 12, eps = 3)['df_return'].iloc[:, 1:]
df_new = LM_generation(n = 1000, beta_hat = [1,1,1,1,0,0,0,0],
    cor = 0.3, n_nuisance = 12, eps = 3)['df_return'].iloc[:, 1:]
df_exist.iloc[:, -1] = np.random.choice(np.arange(1, 5), size=1000)
df_new.iloc[:, -1] = np.random.choice(np.arange(1, 5), size=1000)
RFPerm(df_exist, df_new, loss = l2_multinomial, B = 100)








