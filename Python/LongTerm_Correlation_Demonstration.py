#
import os
import copy
import numpy as np
from scipy import stats
import random
from xgboost import XGBClassifier, XGBRegressor
import lightgbm
import statistics
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer, KNNImputer, MissingIndicator
from sklearn.preprocessing import StandardScaler, OneHotEncoder, OrdinalEncoder

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



'''
- Pre-train Rounds: 
- n_boost:          Total number of boosting rounds.
- df_te
'''
def sequence_prediction_XGB(df_tr, df_te, pre_train_round = 10, n_boost = 150, B = 20, tuning = True):
    df_tr_X = df_tr[:, :-1]
    df_tr_Y = df_tr[:, -1]
    df_te_X = df_te[:, :-1]
    df_te_Y = df_te[:, -1]
    n_tr = df_tr.shape[0]
    n_te = df_te.shape[0]
    df_pred_mat_tr = np.zeros((n_tr, n_boost))
    df_pred_mat_te = np.zeros((n_te, n_boost))
    if tuning:
        param_grid = {
          'XGB__n_estimators': [50, 100, 150, 200],
          'XGB__max_depth': [3, 6, 9],
          'XGB__learning_rate': [0.01, 0.05, 0.1, 0.2],
          'XGB__min_child_weight': [1, 5, 9]
        }
        pipeline = Pipeline(steps = [
            ('preprocessing': ColumnTransformer(
                ('column': OneHotEncoder(),
                 ''
                    )))
            ('xgb', XGBRegressor(
                n_estimators = 200,
                max_depth = 5,
                learning_rate = 0.05,
                subsample = 0.75,
                colsample_bytree = 0.75,
                random_state = 20,
                n_jobs = -1
                ))
            ])
        grid_search = GridSearchCV(
            estimator = pipeline, param_grid = param_grid,
            cv = 5, scoring = 'neg_mean_squared_error', n_jobs = -1)
        grid_search.fit(df_tr_X, df_tr_Y)
        best_pipeline = grid_search.best_estimator_
    else:
        params = {'objective': 'reg:squarederror', 
            "random_state": 42, "colsample_bytree": 0.85,
            "colsample_bylevel": 0.85, "gamma": 0.0,
            'n_estimators': pre_train_round,
            "max_depth": 6, "reg_alpha": 0.05, 
            "learning_rate": 0.05, "reg_lambda": 0.25}
        best_pipeline = Pipeline(steps = [
            ('preprocessing', StandardScaler()),
            ('xgb', XGBRegressor(**params))
            ])
        best_pipeline.fit(df_tr_X, df_tr_Y)
    #initial prediction:
    pred_initial = best_pipeline.predict(df_te_X)
    #Extracting the fitted XGBoost model:
    old_XGB_model = best_pipeline.named_steps['xgb']
    #For each of the model:
    pred_matrix = np.zeros((n_te, n_boost - pre_train_round, B))
    corr_mat = np.zeros((n_boost - pre_train_round, B))
    for i in range(B):
        cur_pipeline = copy.deepcopy(best_pipeline)
        np.random.seed(i)
        inB = np.random.choice(np.arange(n_tr), n_tr, replace = True)
        B_X = df_tr_X[inB, :]
        B_Y = df_tr_Y[inB]
        #The correlation for the later trees here:
        B_X_processed = cur_pipeline.named_steps['preprocessing'].transform(B_X)
        #imposing the learning scheduler:
        for j in range(n_boost - pre_train_round):
            cur_pipeline.named_steps['xgb'].learing_rate = 0.02
            cur_pipeline.named_steps['xgb'].n_estimators = 1
            cur_pipeline.named_steps['xgb'].fit(B_X, B_Y, xgb_model = cur_pipeline.named_steps['xgb'].get_booster())
            pred_matrix[:, j, i] = cur_pipeline.predict(df_te_X)
        print(i)
    for j in range(n_boost - pre_train_round):
        for i in range(B):
            corr_mat[j, i], _ = stats.pearsonr(pred_initial, pred_matrix[:, j, i])
    return corr_mat.mean(axis = 1), corr_mat.std(axis = 1)


df_tr = LMGeneration(n = 2500, p = 50, beta = [1,1,1,1,1,1,1,1], cor = 0.3, eps_noi = 1, mean = 0)
df_te = LMGeneration(n = 2500, p = 50, beta = [1,1,1,1,0,0,0,0], cor = 0.3, eps_noi = 1, mean = 0)


sequence_prediction_XGB(df_tr, df_te, pre_train_round = 10, n_boost = 30, B = 20, tuning = False)


def loss_function(linear_pred, y, task = 'classification'):
    if task == 'classification':
        return - y @ linear_pred + cp.sum(cp.logistic(linear_predictor))
    else:
        return cp.sum_squares(linear_pred - y)
    else:
        raise ValueError('task must be classification or regression')

def estimate_delta(X_T, Y_T, lambda_reg, fitted_model = None,
    f_x = None, task = 'classification'):
    n_T, n_p = X_T.shape
    delta = cp.Variable(p)
    if f_x is None:
        f_x = get_offset_f_x(X_T, fitted_model, task = task)
    linear_predictor = f_x + X_T @ delta
    raw_loss = loss_function(linear_predictor, Y_T, task = task)
    loss_avg = raw_loss/n_T
    #Impose the L1 regularization:
    l1_penalty = lambda_reg * cp.norm1(deltas)
    total_loss = loss_avg + l1_penalty
    problem = cp.Problem(cp.Minimize(total_loss))
    problem.solve(solver = cp.SCS, verbose = False)
    if delta.value is not None:
        return delta.value
    else:
        raise ValueError('optimization did not converge')





        